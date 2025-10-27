#!/usr/bin/env python3
"""Toy driver for `galevo_nitrogen.galaxy_evol` with a fixed, linear SFR.

The helper script wires up a flat star-formation history, redirects all
simulation products into `playing-with-galevo/outputs` and
`playing-with-galevo/plots`, then visualises:
- gas-phase log10(N/O) versus 12+log10(O/H)
- the star-formation history applied during the run
"""
import argparse
import builtins
import os
import sys
from contextlib import contextmanager
from os import PathLike
from pathlib import Path
from typing import Iterator, Union

import matplotlib.pyplot as plt
import numpy as np

# Make sure the repo root is importable when the script lives in a subdirectory.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import galevo_nitrogen  # noqa: E402  (import after sys.path tweak)
from galevo_nitrogen import galaxy_evol  # noqa: E402
from element_abundances_solar import function_solar_element_abundances  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "outputs"
PLOT_DIR = SCRIPT_DIR / "plots"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
SFH_PATH = REPO_ROOT / "SFH.txt"
_SIM_RESULTS_TOKEN = "simulation_results_from_galaxy_evol"


def _redirect_path(path: Union[PathLike, str], output_root: Path, plot_root: Path) -> Path:
    """Translate legacy output paths into the toy script folder structure."""
    if not isinstance(path, (str, PathLike)):
        return Path(path)

    original = Path(os.fspath(path))
    path_str = str(original)
    if _SIM_RESULTS_TOKEN not in path_str:
        return original

    remainder = Path(path_str.split(_SIM_RESULTS_TOKEN, 1)[1].lstrip("/\\"))
    parts = remainder.parts
    if "plots" in parts:
        idx = parts.index("plots")
        before = Path(*parts[:idx]) if idx > 0 else Path()
        after = Path(*parts[idx + 1:]) if idx + 1 < len(parts) else Path()
        return (plot_root / before / after).resolve()

    return (output_root / remainder).resolve()


@contextmanager
def redirect_simulation_outputs(output_root: Path, plot_root: Path) -> Iterator[None]:
    """Patch file-system calls so galevo writes into the toy output folders."""

    original_open = builtins.open
    original_makedirs = os.makedirs
    original_exists = os.path.exists
    original_isdir = os.path.isdir
    original_isfile = os.path.isfile

    def redirected_open(file, mode='r', buffering=-1, encoding=None, errors=None,
                        newline=None, closefd=True, opener=None):
        if isinstance(file, (str, PathLike)):
            redirected = _redirect_path(file, output_root, plot_root)
            if redirected != Path(file):
                if any(flag in mode for flag in ('w', 'a', 'x', '+')):
                    redirected.parent.mkdir(parents=True, exist_ok=True)
                file = str(redirected)
        return original_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

    def redirected_makedirs(path, mode=0o777, exist_ok=False):
        redirected = _redirect_path(path, output_root, plot_root)
        return original_makedirs(redirected, mode=mode, exist_ok=exist_ok)

    def redirected_exists(path):
        redirected = _redirect_path(path, output_root, plot_root)
        return original_exists(redirected)

    def redirected_isdir(path):
        redirected = _redirect_path(path, output_root, plot_root)
        return original_isdir(redirected)

    def redirected_isfile(path):
        redirected = _redirect_path(path, output_root, plot_root)
        return original_isfile(redirected)

    builtins.open = redirected_open
    galevo_nitrogen.os.makedirs = redirected_makedirs  # type: ignore[attr-defined]
    galevo_nitrogen.os.path.exists = redirected_exists  # type: ignore[attr-defined]
    galevo_nitrogen.os.path.isdir = redirected_isdir  # type: ignore[attr-defined]
    galevo_nitrogen.os.path.isfile = redirected_isfile  # type: ignore[attr-defined]
    try:
        yield
    finally:
        builtins.open = original_open
        galevo_nitrogen.os.makedirs = original_makedirs  # type: ignore[attr-defined]
        galevo_nitrogen.os.path.exists = original_exists  # type: ignore[attr-defined]
        galevo_nitrogen.os.path.isdir = original_isdir  # type: ignore[attr-defined]
        galevo_nitrogen.os.path.isfile = original_isfile  # type: ignore[attr-defined]


@contextmanager
def temporary_sfh(sfr: float, epochs: int, destination: Path = SFH_PATH) -> Iterator[None]:
    """Create a flat SFH file (linear SFR in Msun/yr) for the duration of the context."""
    if epochs < 1:
        raise ValueError("Number of epochs must be at least 1")
    if sfr <= 0:
        raise ValueError("SFR must be positive")

    previous_contents = destination.read_text() if destination.exists() else None
    try:
        values = [0.0] + [sfr] * epochs
        target_length = 1302  # legacy code expects at least this many bins
        if len(values) < target_length:
            values.extend([0.0] * (target_length - len(values)))

        with destination.open("w") as handle:
            for val in values:
                handle.write(f"{val:.8e}\n")
            handle.write("# Flat SFH generated by playing-with-galevo.py\n")
        yield
    finally:
        if previous_contents is None:
            try:
                destination.unlink()
            except FileNotFoundError:
                pass
        else:
            destination.write_text(previous_contents)


@contextmanager
def pushd(path: Path) -> Iterator[None]:
    """Temporarily change the working directory."""
    original_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original_cwd)


def parse_epoch_count(value: str) -> int:
    """Parse `--epochs` while enforcing an integer number of 10 Myr bins."""
    try:
        numeric = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"--epochs must be an integer count of 10 Myr bins (received {value!r})"
        ) from exc

    if not numeric.is_integer():
        raise argparse.ArgumentTypeError(
            f"--epochs must be an integer count of 10 Myr bins (received {value!r})"
        )

    count = int(numeric)
    if count < 1:
        raise argparse.ArgumentTypeError("--epochs must be >= 1")
    return count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run a toy galIMF evolution with a constant SFR.")
    # --sfr : absolute star-formation rate in Msun/yr, applied to each 10 Myr epoch (not logarithmic).
    parser.add_argument("--sfr", type=float, required=True,
                        help="Constant SFR in Msun/yr applied to every 10 Myr epoch (linear value).")
    # --yield-table : massive-star yield grid recognised by galevo (e.g. portinari98, WW95, Kobayashi06,
    #                 Limongi_M000, Limongi_R000_sup_agb).
    parser.add_argument("--yield-table", type=str, required=True,
                        help="Massive-star yield table key (e.g. portinari98, WW95, Kobayashi06, Limongi_M000).")
    # --epochs : number of 10 Myr epochs using the supplied SFR (epochs=30 -> 300 Myr of activity).
    parser.add_argument("--epochs", type=parse_epoch_count, default=30,
                        help="Integer count of 10 Myr epochs with the chosen SFR (default: 30, i.e. 300 Myr).")
    parser.add_argument("--initial-metallicity", type=float, default=1.34e-7,
                        help="Initial gas metallicity Z_0 as a mass fraction (default matches galevo_nitrogen).")
    parser.add_argument("--solar-table", type=str, default="Anders1989",
                        help="Solar abundance scale for abundance conversions (default: Anders1989).")
    parser.add_argument("--output", type=Path, default=None,
                        help="Optional path to save the abundance/SFH figure (inside or outside the toy folder).")
    parser.add_argument("--no-show", action="store_true",
                        help="Generate the figure without opening a window.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.sfr <= 0:
        parser.error("--sfr must be strictly positive")
    if args.epochs < 1:
        parser.error("--epochs must be at least 1")

    with pushd(REPO_ROOT):
        with temporary_sfh(args.sfr, args.epochs), redirect_simulation_outputs(OUTPUT_DIR, PLOT_DIR):
            galaxy_evol(
                imf='Kroupa',
                STF=0.5,
                SFEN=args.epochs,
                Z_0=args.initial_metallicity,
                solar_mass_component='Anders1989_mass',
                str_yield_table=args.yield_table,
                IMF_name='Kroupa',
                steller_mass_upper_bound=150,
                time_resolution_in_Myr=1,
                mass_boundary_observe_low=1.5,
                mass_boundary_observe_up=8,
                SFH_model='provided',
                SFE=0.05,
                SNIa_ON=False,
                SNIa_yield_table='Thielemann1993',
                solar_abu_table=args.solar_table,
                high_time_resolution=None,
                plot_show=False,
                plot_save=False,
                outflow=None,
                check_igimf=False,
                tau_infalle9=0,
            )

    o_over_h = np.asarray(galevo_nitrogen.O_over_H_list)
    n_over_o = np.asarray(galevo_nitrogen.N_over_O_list)
    if o_over_h.size == 0 or n_over_o.size == 0:
        raise RuntimeError("No abundance data returned by galaxy_evol; check the simulation setup")

    solar_o = function_solar_element_abundances(args.solar_table, "O")
    solar_n = function_solar_element_abundances(args.solar_table, "N")
    metallicity_track = o_over_h + solar_o
    log_no_track = n_over_o + (solar_n - solar_o)

    sf_records = galevo_nitrogen.all_sfr[:args.epochs]
    sf_times = np.array([entry[1] for entry in sf_records]) if sf_records else np.array([])
    sf_values = np.array([entry[0] for entry in sf_records]) if sf_records else np.array([])

    fig, (ax_track, ax_sfh) = plt.subplots(1, 2, figsize=(10, 4.5))

    ax_track.plot(metallicity_track, log_no_track, marker='o', markersize=2.5, linewidth=0.75)
    ax_track.set_xlabel("12 + log10(O/H)")
    ax_track.set_ylabel("log10(N/O)")
    ax_track.set_title("Gas abundance track")

    if sf_times.size > 0:
        ax_sfh.plot(sf_times, sf_values, drawstyle='steps-post')
    ax_sfh.set_xlabel("Time [Gyr]")
    ax_sfh.set_ylabel("SFR [Msun/yr]")
    ax_sfh.set_title("Star-formation history")
    if sf_values.size > 0:
        ax_sfh.set_ylim(bottom=0, top=sf_values.max() * 1.1)

    fig.suptitle(f"Toy galIMF run: SFR={args.sfr:g} Msun/yr, yields={args.yield_table}", fontsize=11)
    fig.tight_layout()

    if args.output is not None:
        fig.savefig(args.output, dpi=200)
    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
