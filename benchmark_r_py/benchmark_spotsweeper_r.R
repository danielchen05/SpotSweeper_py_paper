library(SpotSweeper)
library(SpatialExperiment)
library(S4Vectors)


# =====================================================================
# Configuration
# =====================================================================

OUTPUT_DIR <- paste0(
    "/users/xchen5/revisions_SpotSweeper-py/SpotSweeper_py_paper/benchmark_r_py"
)

INPUT_PATH <- file.path(
    OUTPUT_DIR,
    "visium_benchmark_input.csv"
)

K <- 36
CUTOFF <- 3
WORKERS <- 1
N_REPEATS <- 10


# =====================================================================
# Read the input previously exported from Python
# =====================================================================

message("Reading shared benchmark input: ", INPUT_PATH)

benchmark_input <- read.csv(
    INPUT_PATH,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

message(
    "Number of spots: ",
    nrow(benchmark_input)
)


required_columns <- c(
    "spot_id",
    "sample_id",
    "x",
    "y",
    "log_total_counts",
    "n_genes_by_counts",
    "pct_counts_mt"
)

missing_columns <- setdiff(
    required_columns,
    colnames(benchmark_input)
)

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Missing required columns:",
            paste(missing_columns, collapse = ", ")
        )
    )
}

if (anyDuplicated(benchmark_input$spot_id)) {
    stop("spot_id values are not unique.")
}

if (anyNA(benchmark_input[, required_columns])) {
    stop("Shared benchmark input contains missing values.")
}


# =====================================================================
# Construct minimal SpatialExperiment
# =====================================================================

#
# Expression values themselves are not needed for localOutliers().
# We therefore construct a one-row dummy assay and place the exact
# benchmark QC values in colData.
#

n_spots <- nrow(benchmark_input)

dummy_counts <- matrix(
    0L,
    nrow = 1,
    ncol = n_spots,
    dimnames = list(
        "dummy_gene",
        benchmark_input$spot_id
    )
)


col_data <- S4Vectors::DataFrame(
    sample_id = benchmark_input$sample_id,
    log_total_counts = benchmark_input$log_total_counts,
    n_genes_by_counts = benchmark_input$n_genes_by_counts,
    pct_counts_mt = benchmark_input$pct_counts_mt,
    row.names = benchmark_input$spot_id
)


spatial_coords <- as.matrix(
    benchmark_input[, c("x", "y")]
)

storage.mode(spatial_coords) <- "double"

rownames(spatial_coords) <- benchmark_input$spot_id


spe_base <- SpatialExperiment::SpatialExperiment(
    assays = list(
        counts = dummy_counts
    ),
    colData = col_data,
    spatialCoords = spatial_coords
)


message(
    "SpatialExperiment constructed: ",
    ncol(spe_base),
    " spots"
)


# =====================================================================
# Benchmark specification
# =====================================================================

metric_names <- c(
    "log_total_counts",
    "n_genes_by_counts",
    "pct_counts_mt"
)

directions <- c(
    log_total_counts = "lower",
    n_genes_by_counts = "lower",
    pct_counts_mt = "higher"
)


# =====================================================================
# Run SpotSweeper R once and save spot-level results
# =====================================================================

spe_results <- spe_base


for (metric in metric_names) {

    message(
        "Running SpotSweeper R: ",
        metric,
        ", direction=",
        directions[[metric]]
    )

    spe_results <- SpotSweeper::localOutliers(
        spe_results,
        metric = metric,
        direction = directions[[metric]],
        n_neighbors = K,
        samples = "sample_id",
        log = FALSE,
        cutoff = CUTOFF,
        workers = WORKERS
    )
}


results_metadata <- as.data.frame(
    SummarizedExperiment::colData(
        spe_results
    )
)


r_results <- data.frame(
    spot_id = rownames(results_metadata),

    log_total_counts_z =
        results_metadata$log_total_counts_z,

    log_total_counts_outlier =
        results_metadata$log_total_counts_outliers,

    n_genes_by_counts_z =
        results_metadata$n_genes_by_counts_z,

    n_genes_by_counts_outlier =
        results_metadata$n_genes_by_counts_outliers,

    pct_counts_mt_z =
        results_metadata$pct_counts_mt_z,

    pct_counts_mt_outlier =
        results_metadata$pct_counts_mt_outliers,

    stringsAsFactors = FALSE
)


R_RESULTS_PATH <- file.path(
    OUTPUT_DIR,
    "r_results.csv"
)

write.csv(
    r_results,
    R_RESULTS_PATH,
    row.names = FALSE
)

message(
    "Wrote R results: ",
    R_RESULTS_PATH
)


# =====================================================================
# Runtime benchmark
# =====================================================================

run_one_metric <- function(
    spe,
    metric,
    direction
) {

    elapsed <- system.time(
        invisible(
            SpotSweeper::localOutliers(
                spe,
                metric = metric,
                direction = direction,
                n_neighbors = K,
                samples = "sample_id",
                log = FALSE,
                cutoff = CUTOFF,
                workers = WORKERS
            )
        )
    )[["elapsed"]]

    as.numeric(elapsed)
}


timing_list <- list()
record_index <- 1


for (metric in metric_names) {

    direction <- directions[[metric]]

    # -------------------------------------------------------------
    # Warm-up run: not recorded
    # -------------------------------------------------------------

    message(
        "Warm-up: ",
        metric
    )

    invisible(
        run_one_metric(
            spe_base,
            metric,
            direction
        )
    )


    # -------------------------------------------------------------
    # Timed repetitions
    # -------------------------------------------------------------

    for (repeat_id in seq_len(N_REPEATS)) {

        elapsed <- run_one_metric(
            spe_base,
            metric,
            direction
        )

        timing_list[[record_index]] <- data.frame(
            implementation = "SpotSweeper R",
            metric = metric,
            repeat_id = repeat_id,
            runtime_seconds = elapsed,
            n_spots = n_spots,
            seconds_per_1000_spots =
                elapsed / n_spots * 1000,
            k = K,
            cutoff = CUTOFF,
            workers = WORKERS,
            stringsAsFactors = FALSE
        )

        message(
            metric,
            ": repeat ",
            repeat_id,
            "/",
            N_REPEATS,
            ": ",
            sprintf("%.6f", elapsed),
            " s"
        )

        record_index <- record_index + 1
    }
}


r_timing <- do.call(
    rbind,
    timing_list
)


R_TIMING_PATH <- file.path(
    OUTPUT_DIR,
    "r_timing.csv"
)

write.csv(
    r_timing,
    R_TIMING_PATH,
    row.names = FALSE
)

message(
    "Wrote R timing: ",
    R_TIMING_PATH
)


# =====================================================================
# Save environment/version information
# =====================================================================

r_environment <- data.frame(
    implementation = "SpotSweeper R",
    r_version = R.version.string,
    spotsweeper_version =
        as.character(
            packageVersion("SpotSweeper")
        ),
    spatialexperiment_version =
        as.character(
            packageVersion("SpatialExperiment")
        ),
    platform = R.version$platform,
    stringsAsFactors = FALSE
)


write.csv(
    r_environment,
    file.path(
        OUTPUT_DIR,
        "r_environment.csv"
    ),
    row.names = FALSE
)


# =====================================================================
# Console summary
# =====================================================================

message("\nSpotSweeper R outlier counts:")

for (metric in metric_names) {

    outlier_column <- paste0(
        metric,
        "_outlier"
    )

    n_outliers <- sum(
        r_results[[outlier_column]]
    )

    message(
        "  ",
        metric,
        ": ",
        n_outliers,
        "/",
        n_spots,
        " (",
        sprintf(
            "%.3f",
            100 * n_outliers / n_spots
        ),
        "%)"
    )
}


message("\nSpotSweeper R median runtimes:")

print(
    aggregate(
        cbind(
            runtime_seconds,
            seconds_per_1000_spots
        ) ~ metric,
        data = r_timing,
        FUN = median
    )
)


message("\nR benchmark complete.")