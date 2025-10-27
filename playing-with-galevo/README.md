# playing-with-galevo

A compact driver around `galevo_nitrogen.galaxy_evol` for exploring closed-box galaxy chemical evolution with a constant star-formation rate (SFR). The script focuses on the nitrogen–oxygen diagnostic plane, while still capturing the full ejecta bookkeeping done by the parent galIMF code.

---

## 1. Conceptual Overview

### 1.1 What the Script Does

`playing-with-galevo.py` automates four tasks:

1. **Star-formation history (SFH) setup** – Builds a temporary `SFH.txt` containing a flat SFR over a user-specified number of 10 Myr epochs.
2. **Model execution** – Invokes `galevo_nitrogen.galaxy_evol` with a simplified Kroupa-IMF configuration, no SNIa, and user-selected massive-star yields.
3. **Results redirection** – Forces every output the core code would normally drop inside `simulation_results_from_galaxy_evol/` to instead land in `playing-with-galevo/outputs/` (text dumps) or `playing-with-galevo/plots/` (raw plot data).
4. **Toy visualisation** – Produces a two-panel figure showing gas-phase `log10(N/O)` vs. `12 + log10(O/H)` and the implied SFR history.

### 1.2 Physics Wrapped by the Script

Under the hood, the galIMF / galevo framework solves a **closed-box, single-zone** chemical evolution problem with the following ingredients:

- **Gas reservoir** – Starts with mass set by the ratio between the desired stellar mass (`STF` parameter; fixed at 0.5) and the total SFR integral implied by the SFH.
- **Stellar IMF** – Here fixed to a classic Kroupa IMF (piecewise power law). The more elaborate IGIMF machinery in `galevo_nitrogen` is bypassed for this toy run, but the infrastructure remains.
- **Massive-star yields** – The `--yield-table` argument selects which tabulated yields to use for Type II SN and massive-star winds. Examples include
  - `portinari98` – Portinari et al. (1998) yields
  - `WW95` – Woosley & Weaver (1995)
  - `Kobayashi06` – Kobayashi et al. (2006)
  - `Limongi_M000`, `Limongi_R300`, `Limongi_R300_sup_agb`, etc.
- **Element tracking** – The code tracks total gas masses and element-specific masses (H, He, C, N, O, Mg, ...), deriving number ratios such as `O_over_H_list` (log10 relative to solar).
- **N/O diagnostic** – The script converts the stored `[O/H]` and `[N/O]` outputs into the commonly plotted `12 + log(O/H)` vs. `log(N/O)` plane.

Assumptions baked into this toy driver:

- **Closed-box**: no inflows or outflows. (The main code supports infall/outflow, but the toy script does not expose those toggles.)
- **Instantaneous mixing**: ejecta mix uniformly with the ISM each timestep.
- **No SNIa**: `SNIa_ON=False` to keep focus on massive-star yields.
- **Fixed time resolution**: 10 Myr epochs set by the `galevo_nitrogen` timestep structure.

---

## 2. Folder Layout

```
playing-with-galevo/
├── playing-with-galevo.py   # the toy runner
├── README.md                # this documentation
├── outputs/                 # redirected text outputs
└── plots/                   # redirected plot-friendly text dumps
```

All paths in this README are relative to the repository root. The script assumes the galIMF repository structure is intact (yield tables, `SFH.txt`, etc.) and will temporarily change into the repo root when executing `galaxy_evol`.

---

## 3. Running the Script

### 3.1 Basic Invocation

```bash
python3 playing-with-galevo/playing-with-galevo.py \
    --sfr 1.0 \
    --yield-table portinari98 \
    --epochs 30 \
    --output playing-with-galevo/plots/portinari98_sfr1_epochs30.png
```

### 3.2 Command-line Options

| Option | Meaning | Notes |
| --- | --- | --- |
| `--sfr` | Constant star-formation rate in linear Msun yr⁻¹ applied to each 10 Myr epoch. | Must be positive. This is **not** log10(SFR). |
| `--yield-table` | Key of the massive-star yield grid to use. | Examples: `portinari98`, `WW95`, `Kobayashi06`, `Limongi_M000`, `Limongi_R000_sup_agb`, `Limongi_R300_sup_agb`. |
| `--epochs` | Integer count of 10 Myr bins at the specified SFR. | Default: 30 (300 Myr). Non-integer values are rejected by `parse_epoch_count`. |
| `--initial-metallicity` | Initial gas metallicity mass fraction `Z_0`. | Default `1.34e-7` (very metal-poor). |
| `--solar-table` | Solar abundance scale for translating log ratios into absolute abundances. | Default `Anders1989`. `Asplund2009` is also supported. |
| `--output` | Optional PNG path for the diagnostic figure. | Can live anywhere; parent folders are not auto-created. |
| `--no-show` | Skip opening a GUI window. | Useful on headless systems; figure is still saved if `--output` is provided. |

### 3.3 Example Workflows

1. **High SFR, brief burst**
   ```bash
   python3 playing-with-galevo/playing-with-galevo.py \
       --sfr 60 \
       --yield-table Limongi_R300_sup_agb \
       --epochs 2 \
       --output playing-with-galevo/plots/L300_sfr60_epochs2.png
   ```

2. **Moderate SFR, 1 Gyr evolution**
   ```bash
   python3 playing-with-galevo/playing-with-galevo.py \
       --sfr 5 \
       --yield-table portinari98 \
       --epochs 100 \
       --output playing-with-galevo/plots/portinari98_sfr5_epochs100.png \
       --no-show
   ```

After each run, inspect:

- `playing-with-galevo/outputs/` – contains all text logs that `galaxy_evol` would normally place in `simulation_results_from_galaxy_evol/...`.
- `playing-with-galevo/plots/` – contains the text files (e.g. `N_over_O_time.txt`) useful for post-processing plus any PNGs you requested via `--output`.

---

## 4. Detailed Code Walkthrough

### 4.1 Imports and Path Bootstrapping
- The script lives in a subdirectory, so it inserts the repository root (`REPO_ROOT`) into `sys.path` before importing `galevo_nitrogen` and other galIMF modules.
- `OUTPUT_DIR` and `PLOT_DIR` are created up front; all redirected files end up there.

### 4.2 Output Redirection Mechanism
- `_redirect_path` intercepts any filename containing the literal `simulation_results_from_galaxy_evol` and remaps it under `outputs/` or `plots/`.
- `redirect_simulation_outputs(...)` monkey-patches `builtins.open`, `galevo_nitrogen.os.makedirs`, and `galevo_nitrogen.os.path` helpers (`exists`, `isdir`, `isfile`) to pass through `_redirect_path`.
- Because the patches happen **inside a context manager**, they are restored as soon as the galaxy run finishes, preventing side-effects on future Python code.

### 4.3 Temporary SFH Writer
- `temporary_sfh(...)` writes `SFH.txt` with a leading zero row (historical quirk of the original code) followed by the chosen constant SFR repeated for `epochs` lines.
- The legacy code expects at least 1302 lines; the script pads with zeros to preserve this requirement.
- On exit, the context manager restores the previous `SFH.txt` (or deletes it if it did not exist).

### 4.4 Working Directory Control
- Some internal routines in `galevo_nitrogen` rely on `os.getcwd()` to resolve relative paths (e.g. `yield_tables/...`).
- The helper `pushd(REPO_ROOT)` temporarily switches to the repo root so all file lookups succeed no matter where you invoke the script from.

### 4.5 CLI Parsing (`build_parser`)
- Inline comments describe the physical meaning of each CLI parameter.
- `parse_epoch_count` ensures the number of epochs is an integer and produces user-friendly errors when it is not.

### 4.6 Model Execution (`main`)
1. Parse CLI arguments and validate `--sfr` & `--epochs`.
2. Enter `pushd(REPO_ROOT)`.
3. Enter both context managers: write temporary SFH and redirect outputs.
4. Call `galaxy_evol(...)` with fixed parameters:
   - `imf='Kroupa'`, `SNIa_ON=False`, `SFE=0.05`, `tau_infalle9=0`, etc.
   - `str_yield_table` comes directly from `--yield-table`.
5. After the run, convert the module-level lists `O_over_H_list` and `N_over_O_list` into plotting coordinates:
   - `metallicity_track = O_over_H_list + log10_solar(O) = 12 + log(O/H)`
   - `log_no_track = N_over_O_list + (log10_solar(N) - log10_solar(O))`
6. Gather `all_sfr` records (SFR and age at each epoch) for the SFH plot.
7. Generate the two-panel matplotlib figure; honour `--no-show` and `--output`.

---

## 5. Scientific Interpretation

### 5.1 Interpreting the N/O vs O/H Track
- **X-axis (`12 + log10(O/H)`)** – Gas-phase oxygen abundance expressed in the standard astronomical scale. Tracks the global metallicity buildup.
- **Y-axis (`log10(N/O)`)** – Nitrogen-to-oxygen ratio relative to solar. Distinguishes between primary nitrogen production (flat trend) and secondary contributions (upturn at higher metallicities).
- The track is time-ordered: markers follow the simulation chronology, though the figure does not label individual ages by default.

### 5.2 Star-Formation History Panel
- Shows the constant SFR applied during the active epochs. Because the SFH is constructed in 10 Myr bins, the plot uses a step function.
- The `all_sfr` capture also includes zero-SFR padding beyond the last active epoch; the plot is truncated to the user-specified number of epochs so it reflects the driver intent.

### 5.3 Caveats and Limitations
- **No gas flows**: inflows/outflows are disabled. Their inclusion would alter metallicity evolution and could be exposed with extra CLI options if needed.
- **No SNIa channel**: iron enrichment is muted; the figure focuses on N/O so this has limited impact, but the text logs still record Fe quantities.
- **Global state**: `galaxy_evol` uses module-level globals. Each run overrides previous results; do not expect persistent storage across runs without saving logs.
- **Time resolution**: the smallest timestep is 10 Myr (epoch length). Shorter-term phenomena are not resolved.

---

## 6. Extending the Script

- **Different IMFs**: expose `imf` and `IMF_name` parameters to compare Kroupa vs. Salpeter vs. IGIMF.
- **Enable SNIa**: set `SNIa_ON=True` and surface the corresponding delay-time distribution options.
- **Variable SFH**: replace `temporary_sfh` with a routine that writes arbitrary SFH arrays (e.g. exponential decline, bursty histories).
- **Infall/Outflow**: surface `tau_infalle9` and `outflow` to explore open-box scenarios.
- **Additional plots**: the redirected logs in `outputs/` include time-series for many species; they can be parsed to make, for example, [α/Fe] vs. [Fe/H] diagrams.

---

## 7. Troubleshooting Checklist

| Symptom | Likely Cause | Fix |
| --- | --- | --- |
| `TypeError` mentioning `|` in type annotations | Running on Python 3.9 or lower with older script version. | Update to current script (uses `Union[...]`). |
| `FileNotFoundError` for `yield_tables/...` | Running the script from inside `playing-with-galevo/` with older version lacking `pushd`. | Update to current script; it now CD's into repo root. |
| No plots produced | Forgot `--output` or figure window closed immediately. | Add `--output` and/or remove `--no-show`. |
| Outputs still appear in `simulation_results_from_galaxy_evol` | Redirection context may have been bypassed (likely due to custom modifications). | Ensure `redirect_simulation_outputs` is active around the `galaxy_evol` call. |

---

## 8. References

- Yan, Z. (2019), *galIMF: Galaxy-wide Initial Mass Function* (original code base).
- Kroupa, P. (2001), *On the variation of the initial mass function*. MNRAS.
- Portinari, L. et al. (1998), *Modelling stellar yields*. A&A.
- Limongi, M., & Chieffi, A. (2018), *Presupernova evolution and explosion of rotating massive stars*. ApJS.
- Kobayashi, C. et al. (2006), *Galactic chemical evolution simulations with nucleosynthesis yields*. ApJ.

For deeper dives into the IGIMF theory and the full parameter space, consult `The_IGIMF_theory.md` and `galevo_nitrogen.py` in the repository.

---

Happy experimenting!
