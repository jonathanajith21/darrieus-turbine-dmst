# Darrieus VAWT: Design and Aerodynamic Analysis

A two-blade Darrieus vertical-axis wind turbine (VAWT) of the curved-blade
"egg-beater" type, designed in FreeCAD and analyzed with Double Multiple
Streamtube (DMST) theory. The central finding of the project is that this
rotor's solidity is too high for DMST to produce valid results, established by
having two independent solvers (a from-scratch MATLAB implementation and the
validated QBlade tool) diverge in the same way.

## What this project is

This was a project I decided to revisit. I designed the rotor in FreeCAD back in my senior year of high school, then set out to predict its performance with
DMST theory. I first wrote a DMST solver from scratch in MATLAB to learn the
method; it diverged to non-physical power coefficients. To test whether that was
my coding error or something inherent to the problem, I rebuilt the analysis in
QBlade, an open-source, validated VAWT tool. QBlade diverged the same way. That
reproduction is the key result: the cause is the rotor's high solidity, not an
implementation bug, and a higher-fidelity method is required for this geometry.

The full story is in the writeup (`darrieus_writeup.docx`).

## Turbine geometry

| Parameter | Value | Source |
|-----------|-------|--------|
| Blade count | 2 | Confirmed from CAD model |
| Chord | 0.067 m | Measured from airfoil part |
| Blade height | 1.0 m | Inferred from part extents |
| Max radius | 0.30 m | Inferred — equatorial (widest) radius |
| Airfoil | NACA 0018 | Symmetric, standard for VAWT blades |
| Solidity (Nc/R) | ~0.45 | High — at/above DMST's reliable limit |

The rotor is a curved-blade troposkein design: the blades bow outward and the
radius is largest at the mid-height equator, tapering toward the top and bottom
mounts.

## The finding

DMST models a VAWT as flow through two actuator disks in series (upwind and
downwind half). For this rotor, the upwind half extracts so much momentum that
the downwind streamtube balance has no valid solution: the induction factor
saturates at 1.0 (local wind velocity driven to zero) across the downwind pass,
and the predicted power coefficient goes strongly negative. Lowering the
relaxation factor and raising the iteration count did not resolve it. Both the
MATLAB solver and QBlade show this. Solidity (Nc/R) here is ~0.45; low-solidity
VAWTs where DMST works well are typically ~0.1-0.2, and ~0.3-0.4 is the
documented onset of streamtube breakdown.

## Repository contents

```
darrieus-turbine-dmst/
├── Darrieus_Motor_Design.FCStd   # FreeCAD rotor design
├── darrieus_writeup.docx         # Full writeup (design, both analyses, finding)
├── first_principles_solver/      # Hand-written MATLAB DMST (learning exercise)
│   └── darrieus_simulation.m
├── results/                      # QBlade plots: Cp vs TSR, induction factor
├── LICENSE
└── README.md
```

## The two solvers

### First-principles MATLAB solver (first_principles_solver/)

A DMST solver written from scratch in base MATLAB: upwind/downwind streamtube
coupling, induction-factor iteration, NACA 0018 polar interpolation over the
full angle-of-attack range, and torque integration. It diverged to non-physical
power coefficients. Kept here as a documented learning exercise — on its own its
divergence could have been a coding error, which is exactly why the QBlade
cross-check matters.

### QBlade cross-check

QBlade (TU Berlin) implements Paraschivoiu's DMST formulation with per-streamtube
induction factors and integrates XFOIL for polar generation plus 360-degree
extrapolation. Run on the same rotor, it diverged in the same way — saturating
induction factor on the downwind pass, negative Cp. Because an independent
validated tool reproduced the failure, the cause is the rotor and the method,
not the implementation.

## Next steps

- Confirm the true equatorial radius from the CAD model. The radius is inferred;
  a larger true radius lowers the solidity and could bring the rotor back into
  DMST's valid range. This is the first thing to check.
- Analyze with a higher-fidelity method (vortex model or CFD) that does not rely
  on the streamtube momentum balance that fails here.
- Optionally explore a lower-solidity design (smaller chord or larger radius) and
  re-run DMST.

## License

MIT — see LICENSE.
