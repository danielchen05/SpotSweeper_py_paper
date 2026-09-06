from pathlib import Path
import platform
import sys
import time
from importlib.metadata import version, PackageNotFoundError

import numpy as np
import pandas as pd
import scanpy as sc

import spotsweeper.local_outliers as lo


# =====================================================================
# Configuration
# =====================================================================

VISIUM_PATH = Path("/users/xchen5/SpotSweeper-py-paper/visium.h5ad")

OUTPUT_DIR = Path(
    "/users/xchen5/revisions_SpotSweeper-py/SpotSweeper_py_paper/benchmark_r_py"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

K = 36
CUTOFF = 3.0
WORKERS = 1
N_REPEATS = 10


# =====================================================================
# Load data
# =====================================================================

print(f"Reading: {VISIUM_PATH}")
adata = sc.read_h5ad(VISIUM_PATH)

print(f"Visium shape: {adata.shape}")
print(f"Available obs columns: {adata.obs.columns.tolist()}")


# =====================================================================
# Prepare shared input used by both Python and R
# =====================================================================

# We benchmark the three primary manuscript metrics:
#
#   1. log total counts       -> lower-tail outliers
#   2. genes detected         -> lower-tail outliers
#   3. % mitochondrial       -> higher-tail outliers
#
# All transformations are completed BEFORE calling either implementation.
# Therefore log=False is used in both R and Python.

if "log_total_counts" not in adata.obs.columns:
    if "total_counts" not in adata.obs.columns:
        raise KeyError(
            "Neither 'log_total_counts' nor 'total_counts' is present."
        )

    adata.obs["log_total_counts"] = np.log1p(
        adata.obs["total_counts"].to_numpy()
    )

required_obs = [
    "log_total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
]

missing = [
    column
    for column in required_obs
    if column not in adata.obs.columns
]

if missing:
    raise KeyError(
        f"Missing required QC columns: {missing}\n"
        f"Available columns: {adata.obs.columns.tolist()}"
    )

if "spatial" not in adata.obsm:
    raise KeyError("'spatial' not found in adata.obsm")


# Single Visium tissue section.
# Preserve an existing sample_id if available; otherwise create one.
if "sample_id" not in adata.obs.columns:
    adata.obs["sample_id"] = "visium"


coords = np.asarray(adata.obsm["spatial"])

if coords.ndim != 2 or coords.shape[1] < 2:
    raise ValueError(
        "adata.obsm['spatial'] must have at least two coordinate columns."
    )


shared_input = pd.DataFrame(
    {
        "spot_id": adata.obs_names.astype(str),
        "sample_id": adata.obs["sample_id"].astype(str).to_numpy(),
        "x": coords[:, 0],
        "y": coords[:, 1],
        "log_total_counts": adata.obs[
            "log_total_counts"
        ].to_numpy(dtype=float),
        "n_genes_by_counts": adata.obs[
            "n_genes_by_counts"
        ].to_numpy(dtype=float),
        "pct_counts_mt": adata.obs[
            "pct_counts_mt"
        ].to_numpy(dtype=float),
    }
)

if shared_input["spot_id"].duplicated().any():
    raise ValueError("spot_id values are not unique.")

if shared_input.isna().any().any():
    raise ValueError(
        "Shared benchmark input contains missing values."
    )


shared_input_path = OUTPUT_DIR / "visium_benchmark_input.csv"
shared_input.to_csv(shared_input_path, index=False)

print(f"Wrote shared input: {shared_input_path}")
print(f"Number of spots: {len(shared_input):,}")


# =====================================================================
# Benchmark specification
# =====================================================================

METRICS = {
    "log_total_counts": {
        "direction": "lower",
    },
    "n_genes_by_counts": {
        "direction": "lower",
    },
    "pct_counts_mt": {
        "direction": "higher",
    },
}


# =====================================================================
# Run SpotSweeper-py once and save spot-level results
# =====================================================================

adata_results = adata.copy()

for metric, spec in METRICS.items():

    print(
        f"Running SpotSweeper-py: "
        f"{metric}, direction={spec['direction']}"
    )

    lo.local_outliers(
        adata_results,
        metric=metric,
        direction=spec["direction"],
        n_neighbors=K,
        sample_key="sample_id",
        log=False,
        cutoff=CUTOFF,
        workers=WORKERS,
        coord_key="spatial",
    )


python_results = pd.DataFrame(
    {
        "spot_id": adata_results.obs_names.astype(str),
    }
)

for metric in METRICS:

    python_results[f"{metric}_z"] = (
        adata_results.obs[f"{metric}_z"]
        .to_numpy(dtype=float)
    )

    python_results[f"{metric}_outlier"] = (
        adata_results.obs[f"{metric}_outliers"]
        .astype(bool)
        .to_numpy()
    )


python_results_path = OUTPUT_DIR / "python_results.csv"
python_results.to_csv(
    python_results_path,
    index=False,
)

print(f"Wrote Python results: {python_results_path}")


# =====================================================================
# Runtime benchmark
# =====================================================================

def run_one_metric(base_adata, metric, direction):
    """
    Run one SpotSweeper-py local-outlier analysis.

    The AnnData copy is performed before timing so that copying the input
    object is not part of the benchmark.
    """

    tmp = base_adata.copy()

    start = time.perf_counter()

    lo.local_outliers(
        tmp,
        metric=metric,
        direction=direction,
        n_neighbors=K,
        sample_key="sample_id",
        log=False,
        cutoff=CUTOFF,
        workers=WORKERS,
        coord_key="spatial",
    )

    elapsed = time.perf_counter() - start

    return elapsed


timing_records = []

for metric, spec in METRICS.items():

    direction = spec["direction"]

    # -------------------------------------------------------------
    # Warm-up run: not recorded
    # -------------------------------------------------------------

    print(f"Warm-up: {metric}")

    _ = run_one_metric(
        adata,
        metric,
        direction,
    )

    # -------------------------------------------------------------
    # Timed repetitions
    # -------------------------------------------------------------

    for repeat in range(1, N_REPEATS + 1):

        elapsed = run_one_metric(
            adata,
            metric,
            direction,
        )

        timing_records.append(
            {
                "implementation": "SpotSweeper-py",
                "metric": metric,
                "repeat": repeat,
                "runtime_seconds": elapsed,
                "n_spots": adata.n_obs,
                "seconds_per_1000_spots": (
                    elapsed / adata.n_obs * 1000
                ),
                "k": K,
                "cutoff": CUTOFF,
                "workers": WORKERS,
            }
        )

        print(
            f"{metric}: repeat {repeat}/{N_REPEATS}: "
            f"{elapsed:.6f} s"
        )


python_timing = pd.DataFrame(timing_records)

python_timing_path = OUTPUT_DIR / "python_timing.csv"
python_timing.to_csv(
    python_timing_path,
    index=False,
)

print(f"Wrote Python timing: {python_timing_path}")


# =====================================================================
# Save minimal environment/version information
# =====================================================================

def package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


environment = pd.DataFrame(
    [
        {
            "implementation": "SpotSweeper-py",
            "python_version": sys.version.split()[0],
            "spotsweeper_version": package_version(
                "spotsweeper"
            ),
            "scanpy_version": package_version("scanpy"),
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
            "platform": platform.platform(),
        }
    ]
)

environment.to_csv(
    OUTPUT_DIR / "python_environment.csv",
    index=False,
)


# =====================================================================
# Console summary
# =====================================================================

print("\nSpotSweeper-py outlier counts:")

for metric in METRICS:
    n = python_results[
        f"{metric}_outlier"
    ].sum()

    print(
        f"  {metric}: "
        f"{n}/{adata.n_obs} "
        f"({100 * n / adata.n_obs:.3f}%)"
    )


print("\nSpotSweeper-py median runtimes:")

print(
    python_timing
    .groupby("metric")[
        [
            "runtime_seconds",
            "seconds_per_1000_spots",
        ]
    ]
    .median()
)


print("\nPython benchmark complete.")