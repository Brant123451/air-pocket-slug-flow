# LASSI — Lagrangian Adaptive Slug Simulator (Fortran reproduction)

Reproduction of the **LASSI** scheme described in the PhD thesis
*"A Lagrangian Slug Capturing Scheme for Gas–Liquid Flows in Pipes"* (Chapter 3).

The full algorithm specification is at `../LASSI_ALGORITHM_SPEC.md`.

## Build

Requirements: Intel Fortran (`ifx`/`ifort`) **or** gfortran (>= 9). Free-form Fortran 2008.

```powershell
# From E:\PHD Thesis\lassi
.\build.ps1                 # builds lassi.exe with gfortran by default
.\build.ps1 -Compiler ifx   # use Intel ifx
.\build.ps1 -Test           # build + run unit tests
```

## Project layout

```
lassi/
  src/
    lassi_kinds.f90        ! KIND, constants (g, π, ε)
    lassi_geom.f90         ! γ(β), h_l, S_l, S_g, S_i, dα/dh, IKH κ
    lassi_eos.f90          ! ρ_g(p), ∂ρ_g/∂p
    lassi_friction.f90     ! Darcy-Weisbach friction, body force F, Bendiksen
    lassi_objects.f90      ! TYPE Object — linked-list node (Bubble/Section/Slug/Inlet/Outlet)
    lassi_grid.f90         ! Linked-list operations (insert/remove/split/merge)
    lassi_riemann.f90      ! Modified-SW Riemann solver (Tables A.1-A.10) + Bendiksen interface
    lassi_press_mom.f90    ! Implicit Pressure-Momentum solve + Thomas + gas-mass correction
    lassi_voidwave.f90     ! Void wave step (Riemann-driven cell-conservative remap)
    lassi_listmgmt.f90     ! Split / Merge / SlugInit / SlugShed / BC
    lassi_io.f90           ! Input parsing, monitor, snapshot writers
    lassi_main.f90         ! Driver: time loop and orchestration
  tests/
    test_geom.f90          ! γ(β) inversion, IKH κ sign
    test_riemann.f90       ! Tables A.1-A.10 vs analytic dam-break
  cases/
    INPUT_KL2018_case4.txt ! K-L slug benchmark
    INPUT_dambreak.txt     ! Modified-SW dam break (verification)
  build.ps1                ! build script
```

## Implementation status

| Module | Status |
|---|---|
| `lassi_kinds`, `lassi_geom`, `lassi_eos`, `lassi_friction` | done |
| `lassi_objects`, `lassi_grid` | done |
| `lassi_riemann` (Tables A.1–A.8, slug-side A.9–A.10) | done |
| `lassi_voidwave` (border move + cell remap) | done |
| `lassi_listmgmt` (split/merge/slug-init) | done |
| `lassi_press_mom` (implicit P-M + gas-mass correction) | **stub** |
| `lassi_io`, `lassi_main` | done |
| Test: `test_riemann` (dam break) | done |

The pressure-momentum implicit step is currently a stub that holds gas pressure and gas velocity steady; the void-wave step alone already reproduces the modified-SW dam break with high accuracy and is sufficient to verify the Riemann/topology core.

## Notation cheat sheet

| Code | Math |
|---|---|
| `beta` | β (liquid holdup) |
| `Ul`, `Ug`, `Ugs`, `Um` | U_l, U_g, U_g^S, U_m |
| `kappa` | κ (IKH coefficient) |
| `pres` | p |
| `phi` | φ (pipe inclination) |
| `D`, `Apipe` | pipe diameter, total area |

