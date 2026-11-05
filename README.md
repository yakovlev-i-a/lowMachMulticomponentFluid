# lowMachMulticomponentFluid

A low-Mach-number, pressure-split solver for transient, reacting, multicomponent
flows. It is implemented as a solver module for OpenFOAM (the OpenFOAM Foundation
"org" branch, version 14).

The solver extends the standard `multicomponentFluid` module. It keeps the species
transport, energy equation, reaction model, turbulence and thermophysical transport
exactly as they are in `multicomponentFluid`, and changes only the
pressure–velocity coupling. The goal is to remove acoustic waves from the solution
so that the time step is limited by the convective Courant number instead of the
speed of sound.

---

## 1. What problem it solves

In a low-Mach (weakly compressible) reacting flow the acoustic and thermodynamic
pressure scales are very different. A fully compressible solver must resolve fast
acoustic waves even when they are physically irrelevant, which forces a small time
step. This solver uses the classical low-Mach splitting

```
p = pRef + pd
```

where

- `pRef` — thermodynamic reference pressure, constant in time and space;
- `pd`  — dynamic (mechanical) pressure, the only pressure field actually solved.

The density is obtained from the equation of state at the constant reference
pressure:

```
rho = rho(pRef, T, Y_i)
```

so the density depends on temperature and composition only, not on the dynamic
pressure. This is analogous to the "incompressible ideal gas" option in Fluent and
to the `p0/pd` splitting of `lowMachBuoyantPimpleFoam`.

---

## 2. What is inherited unchanged

`lowMachMulticomponentFluid` derives from `multicomponentFluid` (which itself
derives from `isothermalFluid` and `basicFluidSolver`). The following are **not**
re-implemented and behave exactly as in `multicomponentFluid`:

- species transport equations for `Y_i`,
- the energy (enthalpy) equation,
- the reaction model (`reactionModel`),
- momentum and thermophysical transport models (turbulence, heat/mass diffusion),
- the PIMPLE outer loop and mesh-motion handling.

Only the momentum predictor and the pressure correction are overridden.

---

## 3. Governing equations actually changed

### 3.1 Momentum predictor

The momentum equation is advanced with the gradient of the **dynamic** pressure:

```
ddt(rho, U) + div(phi, U) + MRF.DDt(rho, U) + divDevTau(U)
    = - grad(pd)                       (without buoyancy)
    = - grad(pd) + rho*g               (with buoyancy, through ghf)
```

For buoyant cases the buoyancy force is assembled consistently with the
`lowMachBuoyantPimpleFoam` convention.

### 3.2 Dynamic pressure equation (strict Poisson)

The acoustic term `psi*ddt(p)` is removed. The remaining low-Mach pressure
equation is

```
ddt(rho) + div(phiHbyA) - laplacian(rhorAUf, pd) = fvModels().sourceProxy(rho, pd)
```

with the mass-flux correction

```
phi = phiHbyA + pdEqn.flux()
```

### 3.3 Density

There is no `thermo.correctRho(p - p0)` step, because the density must not respond
to changes of the dynamic pressure. The density is taken from the EOS,

```
rho = thermo.rho()          // evaluated at p = pRef
```

and then a continuity equation is solved so that `rho` stays discretely consistent
with the face mass flux `phi`:

```
ddt(rho) + div(phi) = fvModels().source(rho)
```

### 3.4 Velocity correction

The cell velocity is reconstructed from the face mass-flux correction rather than
from a cell-centred pressure gradient:

```
U = HbyA + rAU * reconstruct((phi - phiHbyA)/rhorAUf)
```

(plus the buoyancy flux term when buoyancy is active). This keeps the reconstructed
velocity consistent with the corrected mass flux on a collocated grid under
variable density.

---

## 4. Differences from `multicomponentFluid` in one table

| # | Change | Why |
|---|--------|-----|
| 1 | `p = pRef + pd`; `pRef` is constant | filters acoustics, matches low-Mach limit |
| 2 | acoustic term `psi*ddt(p)` removed from the pressure equation | removes sound waves |
| 3 | momentum source uses `-grad(pd)` | consistent with the pressure split |
| 4 | velocity corrected via `reconstruct((phi - phiHbyA)/rhorAUf)` | mass-flux-consistent velocity update |
| 5 | no `thermo.correctRho(p - p0)`; density from EOS at `pRef` | density independent of `pd` |

The standard OpenFOAM interpolation `interp(rho) * interp(HbyA)` is used for
`phiHbyA` (i.e. `rhof*fvc::flux(HbyA)`), the same as in `isothermalFluid`.

---

## 5. Build

Prerequisites:

- OpenFOAM 14, OpenFOAM Foundation ("org") branch, installed and sourced.

Build the solver library with:

```bash
wmake libso
```

This produces `liblowMachMulticomponentFluidSolver.so` in the user library path.

---

## 6. Run

Run a case with:

```bash
foamRun -solver lowMachMulticomponentFluid
```

---

## 7. Case requirements

The case is a standard `multicomponentFluid` case with two additional fields:

| File | Type | Required | Meaning |
|------|------|----------|---------|
| `0/pd` | `volScalarField` | yes (`MUST_READ`) | dynamic pressure, the solved field |
| `constant/pRef` | `uniformDimensionedScalarField` | optional (`READ_IF_PRESENT`) | thermodynamic reference pressure, default 0 |

All other fields (`U`, `p`, `rho`, `T`, `Y_i`, etc.) are the usual
`multicomponentFluid` fields.

---

## 8. Boundary conditions (recommended working setup)

The dynamic pressure `pd` needs a proper reference level. A robust choice is:

- `pd` at the outlet: `fixedValue uniform 0`;
- `pd` at walls: `fixedFluxPressure` (or `zeroGradient` for a simple case);
- `pd` at the inlet: `zeroGradient`.

The velocity at the outlet is typically `pressureInletOutletVelocity`.

A nearly-Neumann-like outlet condition (for example `waveTransmissive`) on `pd`
was found to make the strict Poisson equation poorly conditioned; with a Dirichlet
`pd = 0` at the outlet the strict Poisson form is well posed. This is why
`fixedValue 0` is recommended.

---

## 9. Recommended solver settings

The coupled thermo–flow iteration converges cleanly with tight inner tolerances and
implicit under-relaxation of the pressure equation. A validated configuration is:

```
PIMPLE
{
    momentumPredictor  no;
    nOuterCorrectors   6;
    nCorrectors        2;
    nNonOrthogonalCorrectors 1;
}

solvers
{
    pd
    {
        solver          GAMG;
        smoother        DICGaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }

    "(U|h|k|epsilon)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-10;
        relTol          0.0;
    }
}

relaxationFactors
{
    equations
    {
        U   0.7;     // momentum equation under-relaxation
        pd  0.3;     // dynamic-pressure equation under-relaxation
    }
}
```

### Why the `pd` equation needs under-relaxation

The strict low-Mach Poisson equation for `pd` has no diagonal acoustic damping.
On non-orthogonal meshes or with sharp density gradients (a hot wall adjacent to
cold gas) the pressure matrix is poorly conditioned and the velocity corrector
diverges on the very first step. The fix is **implicit under-relaxation of the
`pd` equation** (not of the field), applied in the solver by `pdEqnN.relax()`.
It adds the diagonal term (Fluent Theory Guide eq. 18.4-23)

    a_P * p / alpha = sum_nb(a_nb * p_nb) + b + (1-alpha)/alpha * a_P * p_old

with alpha = 0.3 for `pd`. The added diagonal entry restores diagonal dominance
without changing the converged solution. These are exactly the Fluent
pressure-based-solver defaults (`Pressure 0.3`, `Momentum 0.7`).

### Important distinction

- `equations { pd 0.3; }` — **implicit** relaxation of the pressure *equation*
  (matrix diagonal); this is the stabiliser. ✅
- `fields { pd 0.3; }` — **explicit** relaxation of the pressure *field* (applied
  after solving); this does **not** help the pressure-velocity coupling. ❌

Do not put `pd` under `fields` — put it under `equations`.

These are case-specific recommendations, not hard-coded into the solver. The number
of outer correctors and the relaxation factor depend on the flow, mesh and chemistry.

## 10. Validation case

The solver was validated on a quasi-one-dimensional laminar methane/air deflagration
in a 50 mm channel with 2000 cells:

- chemistry: two-step global mechanism `twoSCH4Bfer`
  (`CH4 + 1.5 O2 -> CO + 2 H2O`, `CO + 0.5 O2 -> CO2`),
- left boundary: no-slip wall, `U = 0`,
- right boundary: open outlet,
- `pRef = 1e5 Pa`,
- `pd` at outlet: `fixedValue uniform 0`.

Representative results on the mature flame front:

| Quantity | Value |
|----------|-------|
| burnt-gas velocity plateau | ~1.64 m/s |
| laminar flame speed | ~0.286 m/s |
| adiabatic flame temperature | ~2009 K |
| maximum CO2 mass fraction | ~0.121 |
| cold / hot density | 1.117 / 0.167 kg/m3 |

The same case was also run with the standard compressible `multicomponentFluid`
solver. At `maxCo = 0.1` both solvers give essentially the same flame (plateau
~1.64 m/s, no velocity defect). The advantage of the low-Mach solver appears at
higher Courant numbers, where the compressible solver is limited by acoustics while
the low-Mach solver remains stable until the chemical ODE integration becomes the
limiting factor.

---

## 11. License

GPLv3 (see `LICENSE`). The solver is a derivative of OpenFOAM modules
(`multicomponentFluid`, `isothermalFluid`, `basicFluidSolver`), which are
Copyright the OpenFOAM Foundation and distributed under the same license.
