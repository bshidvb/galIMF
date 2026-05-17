import argparse
import glob
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

import IMFs.Kroupa_IMF as kroupa_imf


def load_data_with_names(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    data_dict = {}
    for i, raw_line in enumerate(lines):
        line = raw_line.strip()
        if line.startswith("#"):
            row_name = line[1:].strip()
            if i + 1 < len(lines):
                data = list(map(float, lines[i + 1].strip().split()))
                data_dict[row_name] = data
    return data_dict


def get_kroupa_imf(alpha3=None):
    """Return a Kroupa IMF callable configured with the requested high-mass slope."""
    if alpha3 is not None:
        kroupa_imf.alpha3 = float(alpha3)
        kroupa_imf.integrated_mass = quad(kroupa_imf.mass_function, 0.08, 150, limit=50)[0]
    return kroupa_imf.custom_imf


def _cumulative_trapezoid(y, x):
    """Return cumulative trapezoidal integrals with the same length as the inputs."""
    if y.size == 0:
        return np.array([], dtype=float)
    if y.size == 1:
        return np.array([0.0], dtype=float)

    delta_x = np.diff(x)
    traps = 0.5 * (y[:-1] + y[1:]) * delta_x
    return np.concatenate(([0.0], np.cumsum(traps)))


def _prepare_mass_and_yield(mass, yield_m, mmin=None, mmax=None):
    mass = np.asarray(mass, dtype=float)
    yield_m = np.asarray(yield_m, dtype=float)

    mask = np.isfinite(mass) & np.isfinite(yield_m)
    mass = mass[mask]
    yield_m = yield_m[mask]
    if mass.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    lower = mmin if mmin is not None else mass.min()
    upper = mmax if mmax is not None else mass.max()
    select = (mass >= lower) & (mass <= upper)
    mass = mass[select]
    yield_m = yield_m[select]
    if mass.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    order = np.argsort(mass)
    return mass[order], yield_m[order]


def imf_weighted_yield_from_table(mass, yield_m, imf_func=None, mmin=None, mmax=None):
    """Integrate a tabulated stellar yield against an IMF and normalize by formed stellar mass."""
    if imf_func is None:
        imf_func = kroupa_imf.custom_imf

    mass, yield_m = _prepare_mass_and_yield(mass, yield_m, mmin=mmin, mmax=mmax)
    if mass.size == 0:
        return np.nan

    xi = np.array([imf_func(m) for m in mass], dtype=float)
    num = np.trapz(yield_m * xi, mass)
    den = np.trapz(mass * xi, mass)
    return np.nan if den == 0 else num / den


def imf_weighted_yield_curve(mass, yield_m, imf_func=None, mmin=None, mmax=None, weighting_mode="cumulative"):
    """Return a single-yield curve under the requested IMF-weighting mode."""
    if weighting_mode not in {"cumulative", "pointwise"}:
        raise ValueError(f"Unknown weighting_mode: {weighting_mode}")

    if imf_func is None:
        imf_func = kroupa_imf.custom_imf

    mass, yield_m = _prepare_mass_and_yield(mass, yield_m, mmin=mmin, mmax=mmax)
    if mass.size == 0:
        empty = np.array([], dtype=float)
        return {
            "mass": empty,
            "yield_curve": empty,
            "xi": empty,
            "weighted_yield": empty,
            "cum_weighted_yield": empty,
            "formed_mass_norm": np.nan,
            "weighting_mode": weighting_mode,
        }

    xi = np.array([imf_func(m) for m in mass], dtype=float)
    weighted_yield = yield_m * xi
    formed_mass_norm = np.trapz(mass * xi, mass)

    if weighting_mode == "pointwise":
        return {
            "mass": mass,
            "yield_curve": weighted_yield,
            "xi": xi,
            "weighted_yield": weighted_yield,
            "cum_weighted_yield": np.array([], dtype=float),
            "formed_mass_norm": formed_mass_norm,
            "weighting_mode": weighting_mode,
        }

    cum_weighted_yield = _cumulative_trapezoid(weighted_yield, mass)
    if formed_mass_norm > 0:
        yield_curve = cum_weighted_yield / formed_mass_norm
    else:
        yield_curve = np.full_like(cum_weighted_yield, np.nan, dtype=float)

    return {
        "mass": mass,
        "yield_curve": yield_curve,
        "xi": xi,
        "weighted_yield": weighted_yield,
        "cum_weighted_yield": cum_weighted_yield,
        "formed_mass_norm": formed_mass_norm,
        "weighting_mode": weighting_mode,
    }


def _label_from_yield_key(yield_key):
    return yield_key.replace("_eject_mass", "")


def _default_output_path(title, yield_key, alpha3, weighting_mode):
    safe_title = title.lower().replace(" ", "_").replace("/", "_")
    element = _label_from_yield_key(yield_key)
    return Path("./figs/yield_plots") / f"{safe_title}_{element}_{weighting_mode}_alpha3_{alpha3:.2f}.pdf"


def plot_imf_weighted_yields_vs_mass(
    file_paths,
    yield_key,
    alpha3,
    title,
    output_path=None,
    start_index=5,
    mmin=None,
    mmax=None,
    imf_func=None,
    apply_imf_weights=True,
    weighting_mode="cumulative",
):
    """Plot a single stellar yield against stellar mass with optional IMF weighting."""
    if not file_paths:
        print("No matching yield files found.")
        return []

    if not apply_imf_weights:
        imf_func = lambda mass_value: 1.0
    else:
        imf_func = imf_func or get_kroupa_imf(alpha3)

    paths = sorted(file_paths)
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.rc("font", family="serif")
    results = []

    for path in paths:
        data = load_data_with_names(path)
        if "mass" not in data:
            raise KeyError(f"'mass' row not found in {path}")
        if yield_key not in data:
            raise KeyError(f"'{yield_key}' row not found in {path}")

        curve_data = imf_weighted_yield_curve(
            data["mass"],
            data[yield_key],
            imf_func=imf_func,
            mmin=mmin,
            mmax=mmax,
            weighting_mode=weighting_mode,
        )
        if curve_data["mass"].size == 0:
            continue

        metallicity_vals = data.get("metallicity")
        metallicity = float(metallicity_vals[0]) if metallicity_vals else np.nan
        label_text = metallicity_vals[0] if metallicity_vals else os.path.basename(path)

        log_mass = np.log10(curve_data["mass"])
        with np.errstate(divide="ignore", invalid="ignore"):
            log_yield = np.log10(curve_data["yield_curve"])

        indices = np.arange(log_yield.size)
        valid = np.isfinite(log_mass) & np.isfinite(log_yield)
        if start_index is not None:
            valid &= indices >= start_index

        if not np.any(valid):
            continue

        ax.plot(log_mass[valid], log_yield[valid], label=f"Z={label_text}", lw=2)

        integrated_yield = imf_weighted_yield_from_table(
            data["mass"],
            data[yield_key],
            imf_func=imf_func,
            mmin=mmin,
            mmax=mmax,
        )
        results.append(
            {
                "metallicity": metallicity,
                "label": label_text,
                "log_mass": log_mass[valid],
                "log_yield": log_yield[valid],
                "integrated_yield": integrated_yield,
                "weighting_mode": weighting_mode,
            }
        )

    if not results:
        print("No valid data found for plotting.")
        return []

    element = _label_from_yield_key(yield_key)
    ax.set_xlabel("log stellar mass [M$_\\odot$]", fontsize=14)
    if weighting_mode == "cumulative":
        ax.set_ylabel(f"log cumulative IMF-weighted {element} yield", fontsize=14)
    else:
        ax.set_ylabel(f"log IMF-weighted {element} yield", fontsize=14)

    if apply_imf_weights:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f"{title} IMF weighting disabled", fontsize=16)
    ax.legend(fontsize=8.5)
    plt.tight_layout()

    if output_path is None:
        output_path = _default_output_path(title, yield_key, alpha3, weighting_mode)
    else:
        output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.show()
    return results


def summarize_imf_weighted_yields(file_paths, yield_key, alpha3, mmin=None, mmax=None, imf_func=None, apply_imf_weights=True):
    """Return integrated IMF-weighted single yields for each metallicity grid."""
    if not apply_imf_weights:
        imf_func = lambda mass_value: 1.0
    else:
        imf_func = imf_func or get_kroupa_imf(alpha3)

    summaries = []
    for path in sorted(file_paths):
        data = load_data_with_names(path)
        if "mass" not in data:
            raise KeyError(f"'mass' row not found in {path}")
        if yield_key not in data:
            raise KeyError(f"'{yield_key}' row not found in {path}")

        integrated_yield = imf_weighted_yield_from_table(
            data["mass"],
            data[yield_key],
            imf_func=imf_func,
            mmin=mmin,
            mmax=mmax,
        )
        metallicity_vals = data.get("metallicity")
        metallicity = float(metallicity_vals[0]) if metallicity_vals else np.nan

        summaries.append(
            {
                "metallicity": metallicity,
                "label": metallicity_vals[0] if metallicity_vals else os.path.basename(path),
                "integrated_yield": integrated_yield,
            }
        )

    return summaries


def parse_args():
    parser = argparse.ArgumentParser(description="Plot IMF-weighted single-element stellar yields from rearranged yield tables.")
    parser.add_argument(
        "--pattern",
        default="yield_tables__2024/rearranged___/setllar_N_eject_mass_from_Limongi_R000_sup_agb/*.txt",
        help="Glob pattern for the yield files to plot.",
    )
    parser.add_argument(
        "--yield-key",
        default="N_eject_mass",
        help="Row name to read from each yield file, for example N_eject_mass or O_eject_mass.",
    )
    parser.add_argument(
        "--title",
        default="LC18 R000 with super-AGB",
        help="Plot title prefix.",
    )
    parser.add_argument(
        "--alpha3",
        type=float,
        default=1.0,
        help="High-mass slope for the Kroupa IMF.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional output filename. Defaults to figs/yield_plots/<generated>.pdf",
    )
    parser.add_argument(
        "--start-index",
        type=int,
        default=5,
        help="Skip the lowest-mass points to match the notebook behavior.",
    )
    parser.add_argument("--mmin", type=float, default=None, help="Optional lower stellar-mass limit.")
    parser.add_argument("--mmax", type=float, default=None, help="Optional upper stellar-mass limit.")
    parser.add_argument(
        "--weighting-mode",
        choices=("cumulative", "pointwise"),
        default="cumulative",
        help="How to build the plotted yield curve.",
    )
    parser.add_argument(
        "--disable-imf-weights",
        action="store_true",
        help="Disable IMF weighting and plot the raw single-yield curves instead.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    file_paths = glob.glob(args.pattern)

    results = plot_imf_weighted_yields_vs_mass(
        file_paths=file_paths,
        yield_key=args.yield_key,
        alpha3=args.alpha3,
        title=args.title,
        output_path=args.output,
        start_index=args.start_index,
        mmin=args.mmin,
        mmax=args.mmax,
        apply_imf_weights=not args.disable_imf_weights,
        weighting_mode=args.weighting_mode,
    )

    if not results:
        return

    summaries = summarize_imf_weighted_yields(
        file_paths=file_paths,
        yield_key=args.yield_key,
        alpha3=args.alpha3,
        mmin=args.mmin,
        mmax=args.mmax,
        apply_imf_weights=not args.disable_imf_weights,
    )
    for row in summaries:
        print(f"Z={row['label']} | integrated {args.yield_key}={row['integrated_yield']:.3e}")


if __name__ == "__main__":
    main()
