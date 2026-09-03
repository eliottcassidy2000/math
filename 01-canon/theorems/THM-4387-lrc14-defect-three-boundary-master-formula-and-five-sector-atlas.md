---
id: THM-4387
title: "LRC14 coefficient-one defect-three boundary formula and sharp five-sector atlas"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; LRC(14) OPEN. On the coefficient-one l1=16 boundary, the complete
  lift address is the chart-labelled affine state (delta,k), with delta in
  {-3,0,3} and 3 not dividing k. An exact six-term roof formula and analytic
  tail give sharp maxima for all five coefficient-one patterns. The two
  remaining l1=16 coefficient shapes, arbitrary nonresonance, entry, and
  LRC(14) are outside this theorem.
source: root + lrc_defect_three + lrc_defect_three_cleanroom / continuation session, 2026-09-03
depends_on:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
related:
  - THM-4384-lrc14-small-defect-short-relation-master-formula-and-sharp-sector-atlas
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
primary_script: 04-computation/lrc14_defect_three_coefficient_one_boundary_thm4387.py
primary_output: 05-knowledge/results/lrc14_defect_three_coefficient_one_boundary_thm4387.out
primary_script_sha256: d4fcd040ee514a6d3fbe755c47fd41b6625b827b45356a4e2be32ed56420bc46
primary_output_sha256: 1a1286e86454dee247fef1f6d152ab1e93b5bb28472221f215f980b61ba6fdbc
independent_referee_script: 04-computation/lrc14_defect_three_coefficient_one_boundary_independent_referee_thm4387.py
independent_referee_output: 05-knowledge/results/lrc14_defect_three_coefficient_one_boundary_independent_referee_thm4387.out
independent_referee_script_sha256: b257ea62151cdc375e28d154a5a22731b7d4031052b42af9406fda74b14a38fd
independent_referee_output_sha256: 261294b0efbb9b5e2c09e7f7f4a54870ad70cc4eb60dc4da6c5768215f4bddb6
hash_basis: raw LF bytes
audit: >
  PASS. A dependency-free clean-room implementation independently enumerates
  all seven primitive full-support ternary-unit l1=16 shapes, isolates the
  five coefficient-one rows, rederives the affine torsor and six-term roof,
  proves the unimodal tail, checks the five sharp finite universes, intersects
  4,608 physical triples directly, audits endpoint torsion and chart changes,
  and runs 19,757,759 explicit checks. Normal, optimized, and hash-seeded
  streams match. It explicitly excludes shapes (2,7,7) and (4,5,7).
---

# THM-4387 -- The first coefficient-one defect-three boundary

**PROVED ELEMENTARY RELATIVE TO THM-4373 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED; LRC(14) OPEN.** This is exactly the coefficient-one slice of the
primitive full-support `l1=16` boundary. It proves no universal nonresonance,
seam entry, all-tail synchronization, or LRC(14) conclusion.

## Outcome

Put `r=3/14`.  Let `h+m=15`, `1<=h<m`, and `3` divide neither `h` nor
`m`.  Let `p,b,q` be primitive, distinct, positive odd three-units, and
suppose

```text
m b = s h p + t q,                 s,t in {+1,-1}.       (1)
```

This is the **coefficient-one slice** of the first `l1=16` coefficient
boundary.  It is not the full primitive full-support boundary: the additional
absolute coefficient shapes `(2,7,7)` and `(4,5,7)` are genuine and remain
outside the theorem and verifier below.

For a physical scale-three failure point use the nearest-integer data

```text
n_w=nint(wy),  e_w=wy-n_w,  |e_w|<r,
o_w=-w^(-1)n_w mod 3,  {o_p,o_b,o_q}=F_3.                (2)
```

Define

```text
delta=m n_b-s h n_p-t n_q,          k=b n_p-p n_b.       (3)
```

Then

```text
delta in {-3,0,3},                  3 does not divide k. (4)
```

Conversely, every nonempty lift interval with a pair `(delta,k)` satisfying
(4) is one physical failure component.  The quotient address is the pair
`(delta,k)`, not an endpoint determinant alone.

For `delta in {-3,0,3}` and `k in Z`, put

```text
L(delta,k) = positive_part(min(
    2r/p,
    2r/b,
    2r/q,
    r/p+r/b-|k|/(pb),
    r/p+r/q-|m k+p delta|/(pq),
    r/b+r/q-|s h k+b delta|/(bq))).                      (5)
```

The exact Haar formula is

```text
mu(F_(p,b,q))
   = sum_(delta in {-3,0,3}) sum_(k in Z, 3 does not divide k)
       L(delta,k).                                      (6)
```

It is a finite sum, and `L(-delta,-k)=L(delta,k)`.  Formula (6) is the
correct defect-three replacement for THM-4384's defect-zero trapezoid.

There are **five**, not four, parity-compatible coefficient-one patterns on
this first boundary slice:

```text
(h,m)=(1,14),(2,13),(4,11),(5,10),(7,8).                (7)
```

Their exact coefficient-one sharp atlas is:

| `(h,m,1)` | sharp maximum | unique maximizing speed set | one maximizing presentation | `delta=-3,0,+3` masses | width cutoff | presentations / triples |
|---|---:|---|---|---|---:|---:|
| `(1,14,1)` | `4/133` | `{1,5,19}` | `14*1=19-5` | `2/133, 0, 2/133` | 205 | 298 / 149 |
| `(2,13,1)` | `228/11165` | `{1,11,145}` | `13*11=-2*1+145` | `48/11165, 12/1015, 48/11165` | 409 | 1295 / 1295 |
| `(4,11,1)` | `218/10465` | `{7,13,115}` | `11*13=4*7+115` | `31/10465, 12/805, 31/10465` | 372 | 1270 / 1270 |
| `(5,10,1)` | `16/715` | `{1,11,65}` | `10*1=-5*11+65` | `23/5005, 6/455, 23/5005` | 335 | 1105 / 1104 |
| `(7,8,1)` | `36/1309` | `{1,11,85}` | `8*1=-7*11+85` | `5/1309, 26/1309, 5/1309` | 257 | 792 / 792 |

The first row is structurally important: its extremizer has **zero** mass in
the old defect-zero state.  Treating `delta=+/-3` as a perturbative error
would lose the sharp example.

## Inheritance pass and concept board

The closest proved mechanism is THM-4384's small-defect quotient and exact
period-three quadrature.  Its canonical hostile examples already show that
`delta=+/-3` is possible at `m+h=15`.  The corrected near miss is the
four-row coefficient-one boundary list: `(1,14)` also satisfies every
coefficient parity and three-unit condition and is a genuine, in fact
strongest, fifth sector.  This does not classify the coefficient shapes
`(2,7,7)` and `(4,5,7)`.  The least-used sidecar is the `Z/g` torsion in the
noncoprime endpoint quotient.

The live board is: bounded defect state; affine triple quotient; ternary
owner gate; three-interval roof; endpoint gcd torsion; unimodal sampling
error; and changes of relation chart on multiply resonant speed rays.

## 1. Why precisely three defect states occur

Equation (1) and `n_w=wy-e_w` give

```text
delta=s h e_p+t e_q-m e_b,
|delta|<(m+h+1)r=24/7<4.                               (8)
```

Modulo three, put `A=s h p` and `C=tq`.  The relation says `mb=A+C`.
All three terms are nonzero in `F_3`, so `A=C`.  Since `n_w=-w o_w`
modulo three,

```text
delta=A(o_p+o_b+o_q) mod 3.                             (9)
```

Distinct owners sum to zero, hence `3|delta`.  Combining this with (8)
proves the first assertion in (4).  This is exactly the point where the
THM-4384 mechanism changes: its strict bound was `|delta|<3`, whereas the
new strict bound permits the two adjacent congruence states.

## 2. Each defect sheet is a `Z`-torsor, and the state must be retained

Primitivity and (1) imply `gcd(p,b)=1`: a common divisor would also divide
`q=t(mb-shp)`. On a fixed affine defect plane, `k=b n_p-p n_b` is therefore
onto `Z`; two lifts have the same `k` exactly when their difference is a
physical translation in `Z(p,b,q)`. Thus

```text
{integer lifts with fixed delta}/Z(p,b,q)  --k-->  Z       (10)
```

is a bijection. For `delta!=0` the source is an affine `Z`-torsor, not
canonically a quotient group: adding two lifts does not preserve the defect
sheet. The integer `k` is a complete coordinate only after `delta` and the
relation chart are retained.

The three pairwise determinants are

```text
b n_p-p n_b = k,
q n_p-p n_q = t(m k+p delta),
q n_b-b n_q = t(s h k+b delta).                         (11)
```

Also `o_b-o_p=k/(pb)` modulo three.  Equation (9) continues to say that
the owner sum is zero for every allowed defect.  Hence two owners differ,
and therefore all three differ, exactly when `3` does not divide `k`.
This proves the second assertion in (4).  The owner gate stays simple; it is
the interval geometry that becomes three-dimensional in its bookkeeping.

For a section of (10), write `z=y-n_p/p`.  The three interval centers are

```text
0,  -k/(pb),  -t(mk+p delta)/(pq)                       (12)
```

with radii `r/p,r/b,r/q`.  The length of the intersection of three real
intervals is the least of their three widths and their three pair-overlap
lengths, clipped at zero.  Equations (11)--(12) give (5), and nearest-integer
uniqueness makes different `(delta,k)` classes disjoint.  This proves (6).

When `delta=0`, endpoint eligibility makes the middle interval redundant:
`e_b=(s h e_p+t e_q)/m` and `m>=h+1`.  Thus the zero state alone recovers
THM-4384's old one-determinant formula.  It is false for either nonzero state.

## 3. Bulk integral and a rigorous finite tail

For fixed `delta`, the real extension `k -> L(delta,k)` is a continuous,
compactly supported unimodal roof. Indeed, before clipping it is the minimum of affine
functions; its positive support is an interval.  Its height is at most

```text
2r/W,                    W=max(p,b,q).                   (13)
```

The change of variables

```text
e_p=pz,                 e_b=bz+k/p                       (14)
```

has Jacobian one.  Therefore the integral of the roof is the area in the
`(e_p,e_b)` square cut out by

```text
|e_p|<r, |e_b|<r, |delta-s h e_p+m e_b|<r.              (15)
```

For `delta=0`, the strip is fully contained across the square because
`m>=h+1`, and its area is `4r^2/m=9/(49m)`.  For either sign of `delta=3`,
(15) is a corner triangle with legs `2r/m` and `2r/h`, so its area is
`2r^2/(mh)=9/(98mh)`.

For any nonnegative unimodal `f` and either residue `a=1,2`, monotone
sum-integral comparison on the lattice `a+3Z` gives

```text
sum_j f(a+3j) <= (1/3) integral_R f + sup f.             (16)
```

Apply (16) to both allowed owner residues and all three defect roofs.  Using
(13)--(15) yields the uniform bound

```text
mu(F_(p,b,q))
 <= 6(h+1)/(49mh) + 18/(7W).                            (17)
```

The bulk constants for the five rows are respectively

```text
6/343, 9/637, 15/1078, 18/1225, 6/343.                  (18)
```

The displayed width cutoffs are the smallest integers `N` for which the
right side of (17), with `W=N`, is strictly below that row's proposed
maximum.  Hence an exhaustive scan of `p,b,q<N` proves the global sharp
claim; no unproved discovery bound remains.

The verifier enumerates every signed presentation in those finite boxes,
computes (6) term by term in exact rational arithmetic, and independently
recomputes the physical Haar measure by cutting the `y` circle at every
nearest-integer wall.  The boxes contain 4,760 relation presentations and
4,610 sector-labeled triples.  Two speed sets lie in two new sectors, so the
actual physical union has 4,608 triples; every one agrees.
There is one duplicate presentation within the asymmetric rows:

```text
{5,17,35}: (p,b,q)=(17,5,35) and (35,17,5) in (5,10),
```

and both presentations give `18/4165`.

Their total measure agrees, but their chart-labelled defect masses do not:
they are respectively `(1,16,1)/4165` and `(9,0,9)/4165` on the
`(-3,0,+3)` sheets. Thus even the defect split is not an intrinsic invariant
of the physical speed set.

## 4. The endpoint gcd seam

Let

```text
g=gcd(p,q).
```

Because the triple is primitive, `gcd(g,b)=1`; (1) then gives `g|m`.
Conversely `gcd(p,m)` divides `q`, so in fact

```text
g=gcd(p,m),       P=p/g, Q=q/g, M=m/g,       gcd(P,M)=1. (19)
```

The endpoint quotient is not cyclic when `g>1`:

```text
Z^2/Z(p,q) isomorphic to Z direct_sum Z/g.               (20)
```

If `alpha P+beta Q=1`, complete endpoint coordinates are

```text
j=(q n_p-p n_q)/g=t(Mk+P delta),
a=alpha n_p+beta n_q mod g.                              (21)
```

The owner gate can still be read from the determinant:
`3 does not divide j` iff `3 does not divide k`.  But `j` alone is not a
complete address.  Two defect states can have the same `j` precisely when

```text
M divides delta-delta'.                                  (22)
```

Since the only nonzero differences are `3` and `6` and `3` does not divide
`M`, determinant collisions occur exactly as follows:

| sector | possible nontrivial `g` | `M` | collision |
|---|---:|---:|---|
| `(1,14)` | 7 | 2 | `delta=-3` with `+3` |
| `(2,13)` | 13 | 1 | all three states |
| `(4,11)` | 11 | 1 | all three states |
| `(5,10)` | 5 | 2 | `delta=-3` with `+3` |
| `(7,8)` | none | 8 | none |

In every collision the torsion coordinate `a` separates the endpoint
classes.  Exact witnesses in the verifier include determinant coordinate
`j=1` with three distinct torsion values for the `M=1` seams and two distinct
torsion values for the `M=2` seams.  This is why `(delta,k)` is the safer
natural-number address: normalizing the endpoint determinant and discarding
torsion loses load-bearing information.

## 5. A relation-chart transition explains the strongest new row

The `(1,14)` winner `{1,5,19}` is also the unique positive three-unit ray at
the intersection with the older `(1,4)` coefficient sector:

```text
14*1=19-5,                  4*5=1+19.                    (23)
```

Use `(p,b,q)=(19,1,5)` in the first chart and `(1,5,19)` in the second.
If `(delta,k)` is the first chart's address and `k'` the defect-zero address
in the second, direct elimination gives

```text
k'=-k-delta.                                              (24)
```

The two nonempty new-chart components are

```text
(delta,k)=(-3,4),(3,-4), each of length 2/133,
```

and (24) sends them to `k'=-1,+1`, the two old-chart components.  Thus the
new sharp example is not a new physical component: it is an old component
whose ordinal address splits across two nonzero defect sheets after changing
relation charts.  This is a concrete instance of the repo-wide warning that
an integer rank is meaningful only together with the quotient and sidecars
that make it complete.

The exact coefficient cross-product audit finds thirteen nonempty old/new
sector pairs and two nonempty new/new pairs.  In particular `(7,8)` is
disjoint from all ten THM-4384 sectors and from the other four new sectors at
the level of primitive distinct odd three-unit rays.  Sector measures must
therefore not be added; the coefficient atlas is a cover with sparse, exactly
classifiable chart transitions.

## 6. Exact hostile controls

The verifier checks owner-distinct interior cells with nonzero defect for all
five rows.  One witness per row is:

```text
(1,14): (p,b,q)=(19,1,5), y=211/266, delta=+3;
(2,13): (5,1,23),         y=34/161,  delta=-3;
(4,11): (19,7,1),         y=213/1862,delta=+3;
(5,10): (13,7,5),         y=509/1274,delta=+3;
(7,8):  (23,17,25),       y=2189/5474,delta=+3.
```

These are hostile to any attempt to extend THM-4384 by silently setting
`delta=0`.

## Scope firewall and next tasks

This closes five signed coefficient-one sectors only.  It does not close the
full `l1=16` boundary: the primitive full-support coefficient shapes
`(2,7,7)` and `(4,5,7)` need a separate quotient and exact-measure analysis.
It also does not prove that an arbitrary nonresonant triple enters one of the
five sectors, give a universal triple-comb bound, synchronize several tail
triples, retain a body-component address, or prove a seam-entry/all-tail
lemma.  LRC(14) remains **OPEN**.

The highest-value generated tasks are:

1. Analyze the omitted `(2,7,7)` and `(4,5,7)` shapes; the coefficient-one
   Bezout argument does not automatically supply the same cyclic quotient.
   An independent incoming scout reports provisional targets
   `304/12397` at `{1,23,77}` and `2178/91945` at `{5,37,71}` respectively;
   these are signal for a separate audit, not results of this report.
2. Package (6) as a reusable finite-state roof lemma for general odd
   `m+h`, with the state count controlled by the archimedean defect width.
3. Replace the coarse `18/(7W)` envelope by an exact period-three roof
   discrepancy; this should shrink future coefficient boxes substantially.
4. Build the complete coefficient-chart incidence complex, labeling every
   overlap ray by its integral address-change matrix and gcd torsion.
5. Test whether a seam-entry argument can use the full `(delta,k,a)` packet
   without scalarizing away the body phase.

## Reproduction and frozen evidence

```powershell
python -B 04-computation/lrc14_defect_three_coefficient_one_boundary_thm4387.py
python -O -B 04-computation/lrc14_defect_three_coefficient_one_boundary_thm4387.py
python -B 04-computation/lrc14_defect_three_coefficient_one_boundary_independent_referee_thm4387.py
python -O -B 04-computation/lrc14_defect_three_coefficient_one_boundary_independent_referee_thm4387.py
```

Both verifiers use exact `Fraction` arithmetic, exhaustive finite universes,
definition-level circle cells, and explicit runtime checks that remain active
under optimized Python. Normal, optimized, and fixed-hash-seed streams match
the checked-in outputs. Canonical raw-LF hashes are recorded in the front
matter; no discovery script belongs to the proof surface.
