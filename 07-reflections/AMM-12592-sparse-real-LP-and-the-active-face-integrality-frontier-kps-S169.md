# AMM 12592: sparse real LP and the active-face integrality frontier

Status: **DERIVATION PROVED; NUMERICAL EVIDENCE ONLY at R=512.**  This note
does not claim a rational point, an integer witness, an AMM upper bound, or a
closure of the golden gap.

## Inheritance and question

THM-3329 and the later exact sidecars close the exact golden floor through
`R=256`; the current existence frontier is `(R,D0)=(512,0)`.  The
nonnegative transportation formulation identifies the real relaxation but
records its `R=512` status as open.  The clamped-Pascal conjugacy T2/T6 gives
a sparse causal representation of exactly the same relaxation.

The nearest proved mechanism is the T2 generalized junk flow.  The hostile
example is the plain-clamp death at `R=512,D0=0`, which cannot imply
infeasibility because the same rule dies at the known-feasible `R=256,D0=0`.
The least-used sidecar is the terminal Bernstein identity: if
`j_{R-2,d}=0`, decoding `bern(j,d)/x` at degree `d_last` is again the Pascal
kernel `(1,1)` or `(1,2,1)`.

## Sparse equivalence

For every nonterminal row let `j_{i,t}` be the real junk after selecting an
in-box correction `u_{i,t}`.  With `K_delta=(1,1)` for no degree increase and
`K_delta=(1,2,1)` for a degree increase,

```text
w_i = K_delta j_(i-1) + feed_i,
j_i = w_i-u_i,
-2 C(d_i-1,t) <= u_(i,t) <= 2 C(d_i-1,t-1).             (1)
```

Row zero uses the exact T4 closed form.  Divisibility by `x` is exactly
`j_(i,d_i)=0`, so those variables can be eliminated.  At the terminal row,
the exact Bernstein identity gives

```text
0 <= (K_delta j_(R-2))_t <= 2 C(d_last,t).               (2)
```

Equations (1)--(2), after positive diagonal normalization, are a sparse LP:
each row has at most four nonzero entries.  The maps are invertible affine
changes of the transportation variables, so this is the full real
relaxation, not an early-capture subpolytope.

The exact identity behind (2) was independently tested on pseudorandom
integer Bernstein vectors: whenever `j_d=0`, `bern(j,d)/x`, decoded at
degree `d+delta`, equals convolution by `K_delta` coefficientwise.

## Numerical frontier probe

The companion
`04-computation/amm12592_sparse_real_lp_frontier_kps_s169.py` builds the
sparse matrix using exact integer floor/feed data and floating values only
after positive diagonal scaling.  Clarabel reports feasible points on all
four controls/frontiers:

| R | variables | inequalities | nonzeros | checked max violation |
|---:|---:|---:|---:|---:|
| 32 | 856 | 1,810 | 6,088 | `1.34e-13` |
| 128 | 14,441 | 29,288 | 103,746 | `8.24e-12` |
| 256 | 58,276 | 117,366 | 418,814 | `7.60e-13` |
| 512 | 234,117 | 469,866 | 1,683,470 | `1.96e-13` |

The first three are positive controls with known exact integer witnesses.
The `R=512` row is therefore strong **numerical evidence that the real
transportation relaxation is feasible at the exact floor profile**.  It
localizes the observed construction gap to integer/parity structure rather
than real capacity, but it does not prove that localization.

## Hostile numerical audit and stopping boundary

Three apparently natural upgrades fail and are part of the result:

1. The installed HiGHS presolver falsely declares even the known `R=32`
   control infeasible after scaling.  Unpresolved HiGHS accepts small
   controls but becomes slow; no HiGHS infeasibility verdict is evidence.
2. Nearest even-integer rounding of the Clarabel point dies by rows `8--20`,
   including every known-feasible control.  Relative errors in central
   binomial coordinates are astronomically large in absolute units.
3. Fixed-denominator dyadic projection, including a causal interval clamp,
   also dies on known controls.  The feasible set lies on many forced faces;
   its uniform normalized interior radius is numerically zero.  Generic
   rational rounding therefore has no certified margin.

The exact next target is not another floating solve.  It is one of:

- recover an exact active-face basis and solve that sparse rational system;
- prove a face-preserving rounding theorem in centered junk coordinates; or
- add parity cuts to the sparse system and solve a structured lattice-flow
  problem.

This is the same compression lesson seen in the LRC phase carriers and the
U-spine projective fold: a positive scalar projection preserves feasibility
geometry but loses the discrete lift.  Here the missing sidecar is the
parity/active-face section.
