# Cat states in one- and two-mode \(\mathbb{Z}_3\) Rabi models

This code is used to generate the figures in the paper **“Cat states in one- and two-mode \(\mathbb{Z}_3\) Rabi models”**: https://arxiv.org/abs/2509.08603

## Repository layout

- `src/Z3Rabi.m` – Mathematica package that defines operators, Hamiltonians, eigensolvers, displacement/parity operators, and entropy utilities for one- and two-mode \(\mathbb{Z}_3\) Rabi models.
- `src/WignerFunctions.m` – Mathematica package for computing and plotting Wigner functions (boson cat, two-boson \(\mathbb{Z}_3\) cat, and qutrit-boson Wigner data).
- `scripts/` – Figure-generation scripts used in the paper:
  - `Figure1-One-mode-Z3-Rabi-model-spectrum.wl`
  - `Figure2-Two-mode-Z3-Rabi-model-spectrum.wl`
  - `Figure3-Z2-and-Z3-cat-states.m`
  - `Figure4-Q2B-Wigner-function-of-QB-Z3-cat-state.m`
  - `Figure5-Q2B-Wigner-functions-comparison.wl`
  - `Figure6-Numerical-calculation-of-Q2B-Wigner-functions.wl`
- `output_images/` – Generated figure artifacts (PDF/PNG).

## Requirements

- Wolfram Mathematica / Wolfram Language.
- `MaTeX` package (used by several scripts for TeX-quality labels).

## Usage

Run scripts from the repository root so the package paths resolve correctly.

Example:

```bash
wolframscript -file scripts/Figure1-One-mode-Z3-Rabi-model-spectrum.wl
```

Each script writes its output to `output_images/` using `GetOutputPath[...]` from `WignerFunctions.m`.

## Reproducing figures

- **Figure 1**: one-mode \(\mathbb{Z}_3\) Rabi-model low-energy spectrum (numerical + analytic comparison).
- **Figure 2**: two-mode \(\mathbb{Z}_3\) Rabi-model low-energy spectrum.
- **Figure 3**: comparison of \(\mathbb{Z}_2\) bosonic cat and \(\mathbb{Z}_3\) two-boson cat Wigner functions.
- **Figure 4**: qutrit-two-boson (Q2B) Wigner-function component comparisons.
- **Figure 5**: Q2B cat vs mixed state vs two-boson cat Wigner-function rows.
- **Figure 6**: numerically computed Q2B Wigner functions from eigenstates of the two-mode model at different coupling values.

## Notes

- Photonic Hilbert-space truncation is controlled by `SetNphotTrunc[...]` in `Z3Rabi.m`; most scripts use `SetNphotTrunc[50]`.
- Heavy computations use parallel kernels (`LaunchKernels[]`) in several scripts.
