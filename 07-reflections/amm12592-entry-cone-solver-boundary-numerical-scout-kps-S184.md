# AMM 12592 entry cone: a solver-boundary numerical scout

Status: **NUMERICAL-SCOUT / NO VERDICT.**  The golden-floor profile, feed,
caps, rule-A recurrence, and the `R=32` F1--F3 control are checked with exact
integer arithmetic.  Every LP solve and every residual obtained from it is
floating point.  This note claims no rational point, integer witness,
infeasibility certificate, AMM closure, or deadline improvement.

## Inheritance and salvage decision

The nearest proved mechanism is THM-3332's post-feed F1--F3 cone: an integer
causal prefix entering that cone supplies a certified tail.  The canonical
hostile example is THM-3374's exact `R=512,D0=0` rule-A death at row `107`,
whose 277-bit fatal residual cannot be repaired by changing only the previous
57 rows.  The underused sidecar was an untracked kps-S179 script that built a
smaller real LP only through first feed-free entry.

That old script was worth salvaging because its model is distinct from the
full-epoch sparse LP in kps-S169:

- first feed-free entry at `R=512,D0=0` is row `130`, degree `383`;
- the ordinary F1--F3 entry model has `44,750` variables, `89,528`
  inequalities, and `320,095` nonzeros;
- the stronger zero-entry model has the same variables, `89,146`
  inequalities, and `319,332` nonzeros.

The solver interpretation was not salvageable as evidence.  The replacement
companion
`04-computation/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.py`
therefore treats solver output only as a conditioning diagnostic.  It adds an
exact rule-A control, exact-input and sparse-pattern hashes, dependency
versions, and both scaled and reconstructed-native residuals.  The recorded
transcript is
`05-knowledge/results/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.out`.

## Frozen environment and exact controls

The transcript used Windows 11, CPython `3.13.14`, NumPy `2.3.2`, SciPy
`1.16.1` with embedded HiGHS `1.8.0`, CVXPY `1.9.2`, and Clarabel `0.11.1`.
These versions are printed by every run because numerical statuses can change
with solver and presolver versions.

At `R=32,D0=0`, exact rule A survives to entry and satisfies

```text
F1=True, F2=True, F3=True,
support_max=2, F2_margin=5, F3_min_margin=68.
```

This exact integer check is the positive control.  Merely embedding the same
point into the floating matrix already produces a reconstructed-native row
residual about `2.32e-8`; that residual is rounding noise, while the direct
integer check is authoritative.

At `R=512,D0=0`, the exact control dies at row `107` with forced top residual

```text
-199181334599768561751691151979867147451295075845970943582846950031770442839710071820,
```

independently reproducing the input constant of THM-3374.  It is a hostile
control, not an infeasibility statement about alternative prefixes.

## What the solvers actually say

The following table records reports, not mathematical verdicts.

| model | solver mode | solver report | max scaled violation | reconstructed-native violation |
|---|---|---|---:|---:|
| `R=32` F1--F3 | HiGHS, no presolve | optimal | `8.49e-8` | row `3.04e3` |
| `R=32` F1--F3 | HiGHS, presolve | infeasible | no point | no point |
| `R=32` F1--F3 | Clarabel | optimal | `5.64e-12` | bound `1.11e-4` |
| `R=512` F1--F3 | HiGHS, presolve | infeasible | no point | no point |
| `R=512` F1--F3 | HiGHS, no presolve | status `4` / not set | no point | no point |
| `R=512` F1--F3 | Clarabel | optimal | `1.23e-13` | row `4.35e192`, bound `1.91e195` |
| `R=512` zero entry | Clarabel | optimal | `1.60e-13` | row/bound about `5.91e191` |

The presolved HiGHS declaration at `R=32` is false because the exact rule-A
control above is feasible.  It follows that the analogous `R=512` presolve
status carries no evidentiary weight.  Unpresolved HiGHS also labels the
`R=32` run optimal while missing the requested scaled tolerance and violating
the reconstructed native rows by thousands.

Clarabel's `R=512` points are boundary artifacts in the precise numerical
sense relevant here.  Their normalized residuals are around `1e-13`, but the
largest diagonal scales have more than 200 decimal digits.  Undoing the
scaling turns tiny positive errors in forced sign/zero constraints into
values around `1e191`--`1e195`.  Consequently these arrays are not approximate
rational witnesses in the original coordinates.  This observation also does
**not** prove that an exact real point is absent: an exact feasible point may
exist on the same thin face while the floating solver fails to resolve it.

## Honest use and next exact target

The instrument has two legitimate uses.

1. It is a compact, independently hashed constructor for the first-feed-free
   F1--F3 relaxation, reducing the full-epoch dimension by roughly a factor of
   five.
2. It is a hostile test for any proposed numerical-to-exact upgrade.  Such an
   upgrade must pass the exact `R=32` point and remain valid after pulling
   residuals back through the full diagonal scale ledger.

The sharp next move is exact active-face extraction, not another feasibility
solve.  One may use repeated floating runs only to nominate a support.  Then
one must solve the nominated equalities over `Q`, verify every inequality in
exact arithmetic, and either obtain an exact F1--F3 entry point or an exact
Farkas certificate for that *particular* face.  Raw solver bases and statuses
must not be hashed as canonical data.  THM-3374 adds a causal constraint on
this search: a successful prefix cannot agree with rule A through row `49`,
so an exact reconstruction confined to the late boundary is structurally
incapable of finding the needed repair.

## Reproduction

```text
python 04-computation/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.py 32 0 --exact-rule-a-control
python 04-computation/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.py 32 0 --solver highs --presolve
python 04-computation/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.py 512 0 --solver highs
python 04-computation/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.py 512 0 --solver clarabel
python 04-computation/amm12592_entry_cone_solver_boundary_numerical_scout_kps_s184.py 512 0 --solver clarabel --zero-entry
```

Timings and candidate arrays are environment-dependent.  The direct integer
control, model dimensions, exact-input hashes, and sparse-pattern hashes are
the stable diagnostics.
