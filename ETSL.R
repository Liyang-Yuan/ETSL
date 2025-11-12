# ======================== 依赖管理 ======================== #
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  BiocManager, xcms, MSnbase, S4Vectors, dplyr, tidyr, data.table,
  fuzzyjoin, furrr, future, readr, rlang, purrr, tools, checkmate
)

if (!require("MSnbase", quietly = TRUE)) {
  BiocManager::install("MSnbase", update = FALSE)
}

# ======================== XCMS处理模块 ======================== #
process_ms1 <- function(file_path, 
                        peakwidth = c(5, 30), 
                        noise = 1000, 
                        snthresh = 10) {
  assertFileExists(file_path, extension = "mzML")
  assertNumeric(peakwidth, len = 2, lower = 0)
  
  raw_data <- readMSData(file_path, mode = "onDisk")
  cwp <- CentWaveParam(
    peakwidth = peakwidth,
    noise = noise,
    snthresh = snthresh
  )
  findChromPeaks(raw_data, param = cwp)
}

extract_ms2 <- function(xdata) {
  assertClass(xdata, "XCMSnExp")
  
  ms2_idx <- which(msLevel(xdata) == 2)
  if (length(ms2_idx) == 0) {
    warning("No MS2 scans detected")
    return(NULL)
  }
  
  list(
    precursor_mz = precursorMz(xdata)[ms2_idx],
    precursor_rt = rtime(xdata)[ms2_idx],
    spectra = spectra(xdata[ms2_idx])
  )
}

match_ms1_ms2 <- function(xdata, ms2_data, rt_tol = 5, ppm = 10) {
  if (is.null(ms2_data)) stop("MS2 data required")
  
  ms1_peaks <- as.data.table(chromPeaks(xdata))
  total <- length(ms2_data$precursor_rt)
  message("Matching ", total, " MS2 scans...")
  
  plan(multisession)
  
  results <- future_map(seq_along(ms2_data$precursor_rt), function(i) {
    rt <- ms2_data$precursor_rt[i]
    mz <- ms2_data$precursor_mz[i]
    mz_window <- mz * c(1 - ppm*1e-6, 1 + ppm*1e-6)
    
    matched <- ms1_peaks[
      rtmin <= rt + rt_tol & rtmax >= rt - rt_tol &
        mzmin <= mz_window[2] & mzmax >= mz_window[1]
    ]
    
    if (nrow(matched) > 0) {
      best <- matched[which.max(into)]
      ms2_spec <- ms2_data$spectra[[i]]
      if (length(ms2_spec@mz) == 0) return(NULL)
      
      data.table(
        RT = best$rt,
        mz_MS1 = best$mz,
        Intensity_MS1 = best$into,
        mz_MS2 = ms2_spec@mz,
        Intensity_MS2 = ms2_spec@intensity
      )
    }
  }, .progress = TRUE)
  
  final <- rbindlist(Filter(Negate(is.null), results))
  message("\nMatched ", nrow(final), " features")
  final
}

format_output <- function(matched_data) {
  dt <- as.data.table(matched_data)
  dt[, ms2_num := rowid(RT, mz_MS1, Intensity_MS1)]
  melted <- melt(dt, 
                 id.vars = c("RT", "mz_MS1", "Intensity_MS1", "ms2_num"),
                 measure.vars = c("mz_MS2", "Intensity_MS2"))
  dcasted <- dcast(melted, ... ~ variable + ms2_num, value.var = "value")
  dcasted[, ID := paste0("ID", .I)]
  setcolorder(dcasted, c("ID", "RT", "mz_MS1", "Intensity_MS1"))
  dcasted[]
}

process_single_file <- function(file_path, output_dir) {
  message("Processing: ", basename(file_path))
  tryCatch({
    xdata <- process_ms1(file_path)
    ms2_info <- extract_ms2(xdata)
    
    if (!is.null(ms2_info)) {
      matched_data <- match_ms1_ms2(xdata, ms2_info)
      if (nrow(matched_data) > 0) {
        matched_data_wide <- format_output(matched_data)
        output_name <- paste0(tools::file_path_sans_ext(basename(file_path)), "_Results.csv")
        fwrite(matched_data_wide, file.path(output_dir, output_name))
      }
    }
  }, error = function(e) {
    handle_errors(e, paste("Processing failed:", basename(file_path)))
  })
}

process_mzML_files <- function(input_dir, output_dir) {
  assertDirectoryExists(input_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  mzml_files <- list.files(input_dir, pattern = "\\.mzML$", full.names = TRUE)
  walk(mzml_files, ~ process_single_file(.x, output_dir))
}

# ======================== 匹配模块 ======================== #
read_and_preprocess <- function(path, is_library, ms1_col, compound_col, rt_col, intensity_col, ms2_prefix) {
  assertFileExists(path, extension = "csv")
  
  # 使用纯 data.table 语法
  df <- fread(path)
  df[, row_id := .I]  # 直接添加行号列
  
  # 列存在性验证
  required_cols <- if (is_library) {
    c(ms1_col, compound_col, grep(paste0("^", ms2_prefix), names(df), value = TRUE))
  } else {
    c(ms1_col, rt_col, intensity_col)
  }
  
  assertNames(names(df), must.include = required_cols)
  
  # 使用 data.table 方式过滤 NA
  if (anyNA(df[[ms1_col]])) {
    warning("Removing ", sum(is.na(df[[ms1_col]])), " NA values in ", ms1_col)
    df <- df[!is.na(get(ms1_col))]
  }
  
  # 确保返回纯 data.table
  return(df)
}

perform_ms1_matching <- function(lib, exp_data, ms1_tol, ms1_col, compound_col) {
  # 使用 set 系列函数避免复制
  set(lib, j = "min_mz", value = lib[[ms1_col]] - ms1_tol)
  set(lib, j = "max_mz", value = lib[[ms1_col]] + ms1_tol)
  
  set(exp_data, j = "exp_mz_low", value = exp_data[[ms1_col]])
  set(exp_data, j = "exp_mz_high", value = exp_data[[ms1_col]])
  
  # 设置键前确保数据类型正确
  setkey(lib, min_mz, max_mz)
  setkey(exp_data, exp_mz_low, exp_mz_high)
  
  # 执行区间连接
  matches <- foverlaps(exp_data, lib, nomatch = NULL)
  
  # 计算质量差
  matches[, MS1_Delta := abs(get(paste0("i.", ms1_col)) - get(ms1_col))]
  setorder(matches, MS1_Delta)
  
  # 分组筛选
  matches[, .SD[1:min(.N, 1000)], by = c(compound_col)] %>% 
    .[, .(row_id.lib = row_id, row_id.exp = i.row_id, Compound = get(compound_col), MS1_Delta)]
}

check_ms2 <- function(lib_ms2, exp_ms2, tol) {
  if (length(lib_ms2) == 0 || all(is.na(lib_ms2))) return(list(match = NA, delta = NA_real_))
  if (length(exp_ms2) == 0 || all(is.na(exp_ms2))) return(list(match = FALSE, delta = NA_real_))
  
  diffs <- abs(outer(lib_ms2, exp_ms2, "-"))
  min_diff <- min(diffs)
  list(
    match = min_diff <= tol,
    delta = if (min_diff <= tol) min_diff else NA_real_
  )
}

extract_ms2_values <- function(row, ms2_prefix) {
  cols <- grep(paste0("^", ms2_prefix), names(row), value = TRUE)
  unlist(row[, ..cols]) %>% as.numeric() %>% na.omit()
}

perform_ms2_validation <- function(ms1_matches, lib, exp_data, ms2_tol,
                                   ms1_col, compound_col, rt_col, intensity_col,
                                   lib_ms2_list, exp_ms2_list) {
  future_map2_dfr(
    ms1_matches$row_id.lib,
    ms1_matches$row_id.exp,
    ~ {
      lib_idx <- .x
      exp_idx <- .y
      
      data.table(
        Compound = lib[row_id == lib_idx][[compound_col]],
        MS1_Delta = ms1_matches[row_id.lib == lib_idx & row_id.exp == exp_idx, MS1_Delta],
        MS2_Match = check_ms2(lib_ms2_list[[lib_idx]], exp_ms2_list[[exp_idx]], ms2_tol)$match,
        MS2_Delta = check_ms2(lib_ms2_list[[lib_idx]], exp_ms2_list[[exp_idx]], ms2_tol)$delta,
        RT = exp_data[row_id == exp_idx][[rt_col]],
        Intensity_MS1 = exp_data[row_id == exp_idx][[intensity_col]]
      )
    },
    .progress = TRUE
  )
}

generate_report <- function(results, ms1_tol, intensity_col) {
  results[, Confidence := fcase(
    is.na(MS2_Match) & MS1_Delta < 0.5 * ms1_tol, "Medium (MS1 only)",
    MS1_Delta < 0.5 * ms1_tol & MS2_Match == TRUE, "High",
    MS1_Delta < ms1_tol & MS2_Match == TRUE, "Medium",
    default = "Low"
  )]
  
  results[order(-MS2_Match, MS1_Delta), .SD[1], by = Compound][
    , .(Compound, MS1_Delta, MS2_Delta, Confidence, RT, Intensity_MS1)]
}

# ======================== 主流程 ======================== #
main_process <- function(
    input_library = "library.csv",
    input_data_folder = "data_folder/",
    output_dir = "results/",
    ms1_tol = 0.05,
    ms2_tol = 0.05,
    workers = min(parallel::detectCores() - 1, 8),
    ms1_col = "mz_MS1",
    compound_col = "Compound",
    rt_col = "RT",
    intensity_col = "Intensity_MS1",
    ms2_prefix = "mz_MS2_") {
  
  start_time <- Sys.time()
  plan(multisession, workers = workers)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    # [1/5] 加载标准库
    message("[1/5] Loading library...")
    lib <- read_and_preprocess(input_library, TRUE, ms1_col, compound_col, rt_col, intensity_col, ms2_prefix)
    lib_ms2_list <- lapply(1:nrow(lib), function(i) extract_ms2_values(lib[i], ms2_prefix))
    
    # [2/5] 处理实验数据文件
    data_files <- list.files(input_data_folder, pattern = "\\.csv$", full.names = TRUE)
    if (length(data_files) == 0) stop("No CSV files found")
    
    walk(data_files, function(data_file) {
      file_start <- Sys.time()
      message("\n>> Processing: ", basename(data_file))
      
      tryCatch({
        # [2/5] 加载实验数据
        exp_data <- read_and_preprocess(data_file, FALSE, ms1_col, compound_col, rt_col, intensity_col, ms2_prefix)
        exp_ms2_list <- lapply(1:nrow(exp_data), function(i) extract_ms2_values(exp_data[i], ms2_prefix))
        
        # [3/5] MS1匹配
        ms1_matches <- perform_ms1_matching(lib, exp_data, ms1_tol, ms1_col, compound_col)
        
        # [4/5] MS2验证
        results <- perform_ms2_validation(
          ms1_matches, lib, exp_data, ms2_tol,
          ms1_col, compound_col, rt_col, intensity_col,
          lib_ms2_list, exp_ms2_list
        )
        
        # [5/5] 生成报告
        report <- generate_report(results, ms1_tol, intensity_col)
        output_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(data_file)), "_result.csv"))  # 修复括号
        fwrite(report, output_path)
        message("Saved to: ", output_path, " | Time: ", round(difftime(Sys.time(), file_start, units = "secs")), "s")  # 补充逗号
      }, error = function(e) handle_errors(e, basename(data_file)))
    })
    
    message("\n[Global] Processed ", length(data_files), " files | Total time: ", 
            round(difftime(Sys.time(), start_time, units = "mins")), " mins")
  }, error = function(e) handle_errors(e), finally = plan(sequential))  # 修复括号
}

# ======================== 错误处理 ======================== #
handle_errors <- function(e, current_step = NULL) {
  error_msg <- sprintf(
    "%s [ERROR] %s - %s\nStack trace:\n%s",
    Sys.time(), 
    current_step %||% "Global",
    conditionMessage(e),
    paste(sys.calls(), collapse = "\n")
  )
  writeLines(error_msg, "error.log")
  message(error_msg)
  NULL
}

# ======================== 执行入口 ======================== #
execute_full_workflow <- function(
    mzml_input_dir = "data/",
    library_path = "library.csv",
    xcms_output_dir = "MSresults/",
    final_output_dir = "final_results/",
    workers = min(parallel::detectCores() - 1, 8)) {
  
  tryCatch({
    message("=== XCMS Processing ===")
    process_mzML_files(mzml_input_dir, xcms_output_dir)
    
    message("\n=== Library Matching ===")
    main_process(
      input_library = library_path,
      input_data_folder = xcms_output_dir,
      output_dir = final_output_dir,
      workers = workers
    )
    
    message("\n=== Workflow Completed ===")
  }, error = function(e) handle_errors(e, "execute_full_workflow"))
}

# ======================== 示例执行 ======================== #
execute_full_workflow(
  mzml_input_dir = "D:/博士文件/文章/煎炸文章/Pseudotargeting/匹配算法文件/data/",
  library_path = "D:/博士文件/文章/煎炸文章/Pseudotargeting/匹配算法文件/library/library.csv",
  xcms_output_dir = "D:/博士文件/文章/煎炸文章/Pseudotargeting/匹配算法文件/results/V5/",
  final_output_dir = "D:/博士文件/文章/煎炸文章/Pseudotargeting/匹配算法文件/results/V5/"
)

