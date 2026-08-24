---
id: THM-3991
title: "Periodic unimodular toric cusps have factorial Euler characteristic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a rank-n
  lattice with a lattice-periodic unimodular
  triangulation, quotienting the associated toric central fibre by an
  index-d translation sublattice gives Euler characteristic d*n!: only
  top-simplex orbits contribute, and their number is forced by volume.
  The central fibre has d irreducible components. Consequently a compact
  complex (n+1)-fold fibred by n-tori, with this as its only nonzero-Euler
  fibre, can have the homology of S^(2n+2) only when d*n!=2. In the
  irreducible case this forces n=2. The conclusion is scoped to the smooth
  unimodular translation-quotient grammar; extra Euler fibres, finite
  quotients, singular cones, or a different fan can escape it.
source: root + s6_scout / Hopf-S6 higher-rank fan audit, 2026-08-24
audit: >
  PASS (root plus independent hostile audit, 2026-08-24). The audit checked
  the orbit
  dimensions, absence of translation stabilizers for bounded cells,
  compact-stratification additivity, normalized-volume count, component
  count, all integer solutions of d*n!=2, and the finite-quotient and
  nonunimodular failure boundaries. A separate 1,801-gate exact path verifies
  dimensions `1..5`, all `d<=7`, the `A_2` cell count, the nonunimodular
  five-tetrahedron hostile, and the factorial equation. It repaired ordinary
  versus compactly supported Euler wording, made `n>=1` and preservation of
  `Y_0` explicit, and separated compactness from finite orbit count. The
  theorem is independent of the claimed analytic existence of the S6
  manuscript's quotient.
depends_on: []
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
independent_audit_script: 04-computation/hopf_toric_cusp_factorial_euler_thm3991_independent_audit.py
independent_audit_output: 05-knowledge/results/hopf_toric_cusp_factorial_euler_thm3991_independent_audit.out
independent_audit_script_sha256: 2fcf3788c902dbc3a1a33917ffa7d946af4500fdb3a8afc1ab290d39489dd385
independent_audit_output_sha256: 21a65442a25925e7df74d26a87170483c797c1d35133f3268d84c809ecc004fa
independent_audit_semantic_sha256: ea1b877bb921a870651d13e32b12a51e0fd29c0546f2f4642f098803faf34208
hash_basis: raw LF bytes
---

# THM-3991 -- periodic unimodular toric cusps have factorial Euler characteristic

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** The theorem
isolates the six-dimensional Hopf/S6 fan grammar by an elementary invariant
which is insensitive to the period functions and clutch parameters. It
neither constructs a complex manifold nor verifies the 2026 S6 manuscript.

## 1. Periodic toric central fibres

Let `L` be a rank-`n` lattice with `n>=1`, let `L_R=L tensor R`, and let `T`
be a locally finite `L`-periodic triangulation of `L_R` such that

1. every vertex belongs to `L`;
2. every top-dimensional simplex is unimodular, hence has Euclidean volume
   `covol(L)/n!`.

Form the fan `Sigma_T` in `(L direct_sum Z)_R` whose cones are the cones over
the cells `P times {1}` of `T`, together with their faces, and let `Y` be the
associated smooth toric `(n+1)`-fold. The last lattice coordinate defines a
toric function `t`; write

```text
Y_0={t=0}.                                                  (1)
```

Let `Gamma subset L` be a translation sublattice of finite index

```text
d=[L:Gamma].                                                (2)
```

Allow the translation action to be multiplied by holomorphic torus factors
depending on `t`, but assume:

1. `Gamma` preserves `Y_0` and sends the toric orbit `O_P` to
   `O_(P+lambda)`;
2. its action on a small invariant neighborhood `Y_epsilon` is free and
   proper;
3. the quotient

   ```text
   W=Y_0/Gamma                                               (3)
   ```

   is compact.

Then

```text
chi(W)=d*n!,                                                (4)
```

and `W` has exactly `d` irreducible components. In particular, `W` is
irreducible if and only if `d=1`.

## 2. Orbit-stratification proof

A cell `P` of dimension `r` determines a toric orbit

```text
O_P isomorphic to (C*)^(n-r).                               (5)
```

No nonzero lattice translation stabilizes a bounded cell. Hence passage to
the `Gamma` quotient merely identifies the distinct strata in each
translation orbit; local finiteness, `L`-periodicity, and finite index make
the orbit set finite:

```text
W=disjoint_union_[P in T/Gamma] O_P.                        (6)
```

Use compactly supported Euler characteristic for this finite locally closed
stratification. It is additive; it agrees with ordinary Euler characteristic
on the compact space `W`, and the two values below vanish or equal one in
either convention. Since

```text
chi((C*)^m)=0 for m>0,             chi(point)=1,           (7)
```

only the `n`-simplex orbits contribute. Thus

```text
chi(W)=number of n-simplices in T/Gamma.                   (8)
```

The interiors of those simplices tile the real torus `L_R/Gamma` up to their
measure-zero faces. Its volume is `d*covol(L)`, while every simplex has volume
`covol(L)/n!`. Therefore the number in `(8)` is exactly `d*n!`, proving `(4)`.

The irreducible components of a toric central fibre correspond to the rays
through its vertices. Because every vertex lies in `L` and `T` is
`L`-periodic, the vertex set is one `L`-orbit; modulo `Gamma` it has exactly
`[L:Gamma]=d` orbits. This proves the component assertion.

For the `A_2` triangulation used in the S6 manuscript, `n=2,d=1`: the quotient
has one vertex orbit, three edge orbits, and two triangle orbits, so `(4)`
recovers `chi(W)=2` without using its attaching maps or homology claims.

## 3. One-Euler-fibre corollary

Let

```text
f:X -> C                                                     (9)
```

be a proper holomorphic map from a compact complex `(n+1)`-fold to a compact
curve. Assume it is a smooth `n`-torus bundle off finitely many points, one
special fibre is the `W` above, and every other special fibre has ordinary
topological Euler characteristic zero. Then

```text
chi(X)=d*n!.                                               (10)
```

Indeed, use compactly supported Euler characteristic on the regular locus and
special-point strata of the base. The regular part contributes zero because
the torus fibre has Euler characteristic zero; all special contributions
except `W` vanish by assumption, and `(4)` gives `(10)`. Since `X` is compact,
this compactly supported value is its ordinary Euler characteristic.

If `X` has the integral homology of `S^(2n+2)`, then `chi(X)=2`. Hence

```text
d*n!=2.                                                    (11)
```

The positive integer solutions are

```text
(n,d)=(1,2) or (2,1).                                      (12)
```

In particular, if the central fibre is irreducible, then `d=1` and the sphere
homology condition forces `n=2`. Thus the complex-threefold/S6 dimension is
isolated inside this exact grammar; no same-form construction in complex
dimension at least four can be a homology sphere.

## 4. Exact escape routes

The hypotheses are load-bearing.

1. **Extra singular fibres.** Other nonzero Euler contributions can change
   `(10)`.
2. **Free finite quotient.** If a finite group `H` acts freely, then
   `chi(W/H)=d*n!/|H|`; sphere Euler would require `|H|=d*n!/2`. This is only
   an arithmetic condition; no such free action is asserted.
3. **Fixed-point quotient.** With fixed points, ordinary Euler characteristic
   obeys the Burnside fixed-locus formula, not simple division. Orbifold Euler
   weights do not replace the coarse topology or smoothness audit.
4. **Nonunimodular fan.** Normalized simplex volumes still sum to `d*n!`, but
   ordinary Euler characteristic counts simplex orbits, not their volumes.
   Changing the count requires singular toric cones and then invoices their
   resolution corrections. A five-tetrahedron cube triangulation with
   normalized volumes `(1,1,1,1,2)` gives, after an index-two translation
   quotient, top-cell count `10` while the volume sum is `12`.
5. **Nonreduced structure.** Multiplicity alone does not alter the ordinary
   topological Euler characteristic of the reduced fibre.
6. **Different degeneration grammar.** A non-translation quotient, a
   non-simplicial fan, or a central fibre not stratified as in `(6)` lies
   outside the theorem.
7. **Noncompact/noncocompact quotient.** Compactness is used to identify
   compactly supported and ordinary Euler characteristic. If finite-index
   cocompactness is also dropped, `T/Gamma` may be infinite and `(8)` no
   longer gives a finite ordinary Euler count.
8. **Period/vertex mismatch.** The component count uses both `L`-periodicity
   and vertices in `L`, which make the nonempty vertex set one `L`-orbit. If
   the vertex lattice and translation-period lattice differ, several vertex
   orbits may survive and the assertion "exactly `d` components" must be
   recomputed.

The obstruction is therefore not a universal no-go for complex spheres. It
is a sharp stopping rule for the most immediate higher-rank extension of the
periodic unimodular toric cusp used in the S6 manuscript.

## 5. Reproduction

The independent companion checks the orbit and volume formulas across exact
dimension/index tables, the `A_2` control, the factorial equation, the
five-tetrahedron nonunimodular hostile, and the arithmetic finite-quotient
boundary. These gates audit the displayed combinatorics; they do not assert
that any finite escape action or global complex manifold exists.

```bash
python3 04-computation/hopf_toric_cusp_factorial_euler_thm3991_independent_audit.py
python3 -O 04-computation/hopf_toric_cusp_factorial_euler_thm3991_independent_audit.py
sha256sum 04-computation/hopf_toric_cusp_factorial_euler_thm3991_independent_audit.py \
  05-knowledge/results/hopf_toric_cusp_factorial_euler_thm3991_independent_audit.out
```

The Hopf/S6 global claim remains a preprint claim under its separate audit.
