from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
LIFETIME_ROOT = REPO_ROOT / "yield_tables__2024" / "rearranged___"
OUTPUT_DIR = REPO_ROOT / "figs" / "yield_plots" / "stellar_lifetimes_by_dataset"

DATASETS = {
    "K10, D14 and N13 without HNe": {
        "directory": LIFETIME_ROOT / "setllar_lifetime_from_Nomoto_sup_agb",
        "base_color": "tab:green",
    },
    "K10, D14 and LC18 R300": {
        "directory": LIFETIME_ROOT / "setllar_lifetime_from_Limongi_R300_sup_agb",
        "base_color": "tab:orange",
    },
    "K10 and LC18 R300": {
        "directory": LIFETIME_ROOT / "setllar_lifetime_from_Limongi_R300",
        "base_color": "tab:blue",
    },
        "K10, C22 and LC18 M300": {
        "directory": LIFETIME_ROOT / "setllar_lifetime_from_Limongi_M300",
        "base_color": "tab:brown",
    },
}

INVALID_LIFETIME_YEARS = 1e20
YEARS_PER_MYR = 1e6


def load_data_with_names(file_path: Path) -> dict[str, list[float]]:
    data_dict: dict[str, list[float]] = {}
    with file_path.open("r") as handle:
        lines = handle.readlines()

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("#") and idx + 1 < len(lines):
            row_name = stripped[1:].strip()
            data_dict[row_name] = list(map(float, lines[idx + 1].strip().split()))
    return data_dict


def valid_mass_and_lifetime(data: dict[str, list[float]]) -> tuple[np.ndarray, np.ndarray]:
    mass = np.asarray(data["mass"], dtype=float)
    lifetime = np.asarray(data["lifetime"], dtype=float)
    valid = (
        np.isfinite(mass)
        & np.isfinite(lifetime)
        & (mass > 0)
        & (lifetime > 0)
        & (lifetime < INVALID_LIFETIME_YEARS)
    )
    return mass[valid], lifetime[valid]


def metallicity_sort_key(file_path: Path) -> float:
    data = load_data_with_names(file_path)
    return float(data["metallicity"][0])


def plot_stellar_lifetimes_by_dataset() -> list[Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rc("font", family="serif")

    saved_paths: list[Path] = []

    for dataset_name, config in DATASETS.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        files = sorted(config["directory"].glob("*.txt"), key=metallicity_sort_key)
        colors = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

        for color, file_path in zip(colors, files):
            data = load_data_with_names(file_path)
            mass, lifetime = valid_mass_and_lifetime(data)
            if mass.size == 0:
                continue

            metallicity_label = data["metallicity"][0]
            ax.plot(
                np.log10(mass[:]),
                np.log10(lifetime[:] / YEARS_PER_MYR),
                color=color,
                linestyle="-",
                marker="o",
                markersize=3.5,
                linewidth=2,
                alpha=0.9,
                label=f"Z={metallicity_label:g}",
            )
        
        ax.axhspan(np.log10(30.0), np.log10(110.0), color="magenta", alpha=0.25, zorder=0)
        ax.axhline(np.log10(110.0), color="magenta", linestyle="--", linewidth=1.5, alpha=0.8)
        ax.axhline(np.log10(70.0), color="black", linestyle="--", linewidth=2, alpha=1)
        ax.axhline(np.log10(30.0), color="magenta", linestyle="--", linewidth=1.5, alpha=0.8)
        ax.text(0.95, 2.1, 'GN-z11 age from Jiang21', fontsize = 14, color ="purple")
        ax.set_xlabel(r"log stellar mass [M$_\odot$]", fontsize=14, fontfamily="serif")
        ax.set_ylabel("log stellar lifetime [Myr]", fontsize=14, fontfamily="serif")
        ax.set_title(
            f"Stellar lifetimes for {dataset_name}",
            fontsize=16,
            fontfamily="serif",
        )
        ax.legend(prop={"family": "serif", "size": 12})

        fig.tight_layout()
        output_path = OUTPUT_DIR / f"stellar_lifetimes_{dataset_name}.pdf"
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.show()
        saved_paths.append(output_path)

    return saved_paths


if __name__ == "__main__":
    saved_paths = plot_stellar_lifetimes_by_dataset()
    for saved_path in saved_paths:
        print(f"Saved plot to {saved_path}")
