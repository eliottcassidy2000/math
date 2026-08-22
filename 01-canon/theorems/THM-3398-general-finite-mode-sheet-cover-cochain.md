---
id: THM-3398
title: "General finite-mode sheet covers and complete affine cochains"
status: >
  PROVED analytic all-q exact selected-block mode/cochain criterion plus
  FINITE-EXACT local q=2,...,29 and pair/global q=2,...,11 companion audit.
  INDEPENDENTLY EVENT-, CRT-, AND PROOF-AUDITED.
  Every owner is expanded into
  finitely many consecutive phase-block modes.  A family has a full cyclic
  sheet cover iff some mode blocks cover and their affine centre lattices
  carry one complete closed cochain.  This resolves wraparound and duplicate
  unwrapping, includes the m=1 core mode, recovers the q<=7 one-coset theorem
  without an owner-count restriction, gives an exact certificate-preserving
  dilation law, and gives an exact q=9 triple-domino zero-cochain control.  No
  refined-ledger decrement or LRC(14) follows.
source: root-2608-q8-multiblocker-2026-08-14
audit: self-contained grid-arc and generalized-CRT/Helly proof; 42328 modes and 47200 exact phase cells through q=29, including block lengths three through five; 40187 independent pair-lattice comparisons through q=11; 1669 bounded global subsets and 569 exact event sweeps through q=11; q7/q8 strict boundary, q8 zero-tie and nonclosure hostiles, q9 zero-tie, m=1 control; independent audit re-derived every direction and checked 3520 fixed-mode event/cochain cases through q=12, 892914 grid memberships through first length-11 regimes, 35748 sheet modes through q=25, and 50625 generalized-CRT instances
depends_on:
  - THM-3395-small-sheet-typed-cover-star-cochain
related:
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
  - THM-3388-three-sheet-phase-triangle-cover-clutter
  - THM-3389-four-sheet-typed-cover-clutter
artifacts:
  - 04-computation/lrc14_q8_q14_finite_mode_clutter_probe_20260814.py
  - 05-knowledge/results/lrc14_q8_q14_finite_mode_clutter_probe_20260814.out
  - 04-computation/lrc15_first_effective_triphase_mode_probe_20260814.py
  - 05-knowledge/results/lrc15_first_effective_triphase_mode_probe_20260814.out
  - 04-computation/lrc_q8_q15_full_physical_clutter_audit_20260815.py
  - 05-knowledge/results/lrc_q8_q15_full_physical_clutter_audit_20260815.out
script: 04-computation/lrc_general_finite_mode_sheet_cover_cochain_thm3398.py
output: 05-knowledge/results/lrc_general_finite_mode_sheet_cover_cochain_thm3398.out
script_sha256: 82929cbf6903701533c1b1f6ebed143e5c8f9edc570dfe2895cf8db70e478da9
output_sha256: ab25331039813f8c83626a66d0d0d8157e8b3826a76fccc690452a2cdad3169b
semantic_sha256: b6d3663a25b45b732d4acc597dc2772a5adfd9860934dab7a26c4174f268f6ff
hash_basis: LF-normalized bytes
---

# THM-3398 -- general finite-mode sheet covers and complete affine cochains

**PROVED analytic all-`q` exact selected-block mode/cochain criterion +
FINITE-EXACT local `q=2,...,29` and pair/global `q=2,...,11` companion
audit + INDEPENDENTLY EVENT-, CRT-, AND PROOF-AUDITED.**

## 1. Inheritance and connection contract

THM-3395 proves the one-kernel-coset criterion through seven sheets.  Its
sharp hostile is already the right instruction: at `q=8`, one odd speed can
block two adjacent phase classes.  The repair is to retain a finite
consecutive-block **mode** for each owner, not to discard the owner or force a
tournament.

| field | exact connection |
|---|---|
| source | the danger sets of positive owner speeds on a cyclic `q`-sheet fibre |
| target | a finite selected-block cover clutter decorated by an affine integral complete `1`-cochain |
| map | record `(g,m,a;r,s)`, the sheet block, centre residue `h`, width `w`, and lifted pair gaps `p_ij` |
| preserved | owner identity, simultaneous multi-sheet firing, cyclic order, strict endpoints, tooth lattices, and common source phase |
| destroyed by the bare block cover | centre lift, interval radius, affine congruence, and triangle closure |
| required sidecar | `(h,w)` for every mode and one compatible star, equivalently one complete closed cochain |
| cheapest decisive tests | the `q=8` zero-tie partition and the `q=8` `4+2+2` pairwise-compatible nonclosure hostile |

A mode below means **at least this block fires**.  It does not assert that
the block is the owner's exact firing set.  Extra neighbouring hits are
allowed and harmless for a cover.  This monotone typing is load-bearing at
mode boundaries.

## 2. The finite mode bank

Fix an integer `q>=2` and a positive speed `u`.  For source time `t` put

```text
D_(q,u)(t)={ell in Z/qZ: ||u(t+ell/q)||<1/14}.        (1)
```

Write

```text
g=gcd(u,q),                 m=q/g,
a=u/g mod m.                                           (2)
```

For `m>1`, `a` is a unit modulo `m`; for `m=1` use its unique residue zero.
The `q` sheets map to the `m` equally spaced phase classes by

```text
ell |-> a ell mod m,                                  (3)
```

with every phase class having `g` sheet preimages.

A **selected-block mode** is a pair

```text
r in Z/mZ,              1<=s<=ceil(m/7).              (4)
```

Its phase and sheet blocks are

```text
C(r,s)={r,r+1,...,r+s-1} mod m,
B_(q,u)(r,s)={ell:a ell mod m in C(r,s)}.              (5)
```

Thus `|B_(q,u)(r,s)|=gs`.  Define the centre and width residues

```text
h=-g(2r+s-1) mod 2q,
w=q-7g(s-1).                                          (6)
```

Condition `(4)` is equivalent to `w>0`.  The exact source-time set on which
every sheet in `(5)` is dangerous is

```text
I_(q,u)(r,s)
 = union_(n in Z) (
     (n+h/(2q))/u-w/(14qu),
     (n+h/(2q))/u+w/(14qu)).                          (7)
```

In particular each owner has exactly

```text
m ceil(m/7)                                           (8)
```

selected-block modes.  When `m=1`, there is one mode: it blocks all `q`
sheets, has `h=0,w=q`, centre lattice `(1/u)Z`, and source radius `1/(14u)`.
This is the algebraically complete core mode.  Transverse-cover applications
may omit it by requiring `q` not to divide `u`.

## 3. Grid-arc lemma and wrap uniqueness

Let `y=ut mod 1`.  The danger arc `(-1/14,1/14)` has length `1/7<1/2`.
Therefore its intersection with the `m`-grid is either empty or one
consecutive cyclic block.  If it contains `s` grid points, their unique short
unwrapping has span `(s-1)/m`, so

```text
(s-1)/m<1/7,
7(s-1)<m,
s<=ceil(m/7).                                         (9)
```

Conversely, unwrap `C(r,s)` as the real list

```text
r/m,(r+1)/m,...,(r+s-1)/m.                            (10)
```

The common translate into the danger arc is an open interval centred at

```text
-(2r+s-1)/(2m) mod 1                                  (11)
```

with phase radius

```text
R(m,s)=1/14-(s-1)/(2m)
      =(q-7g(s-1))/(14q)=w/(14q).                     (12)
```

Equations `(11)`--`(12)`, followed by division by `u`, give `(6)`--`(7)`.

This also closes the wraparound gap.  Since the span in `(9)` is less than
`1/2`, the short lift `(10)` is unique up to adding one common integer.  A
wrapped block merely lets `(10)` pass through `m`; reducing its sheets modulo
`m` changes `(11)` by an integer and leaves `h mod 2q` unchanged.  Since
`s<m` for `m>1`, the cyclic set has a unique start `r`; at `m=1` there is
only the single mode already described.  There is no duplicate unwrapping.

At an endpoint of `(7)`, one extreme phase is exactly at distance `1/14`, so
the open inequality fails.  In particular, `m=7,s=2` would have zero width
and is impossible, whereas `m=8,s=2` has phase radius `1/112`.  This is the
strict `q=7/q=8` boundary.

Finally, an actual nonempty set `D_(q,u)(t)` is its unique maximal block
`B_(q,u)(r,s)` from `(5)`.  Conversely, membership in `(7)` guarantees every
sheet in the selected block, though it may also fire a neighbouring sheet.

## 4. Complete affine mode cochains

Choose finitely many owners `u_0,...,u_(d-1)` and one mode
`(r_i,s_i)` for each.  Write `B_i,h_i,w_i` for its data.  A centre lift is

```text
x_i=(n_i+h_i/(2q))/u_i,              n_i in Z.        (13)
```

For an oriented pair define

```text
p_ij=2q u_i u_j(x_i-x_j).                            (14)
```

The exact pair fibre is

```text
p_ij == h_i u_j-h_j u_i  (mod 2q gcd(u_i,u_j)),
7|p_ij|<w_i u_j+w_j u_i.                              (15)
```

The congruence is exactly the condition that `(14)` comes from two centre
lattices `(13)`.  The inequality is exactly strict overlap of the two mode
intervals, whose radii are `w_i/(14qu_i)` and `w_j/(14qu_j)`.

Call antisymmetric integers `p_ij=-p_ji` a **complete mode cochain** when
they satisfy `(15)` for every pair and

```text
u_h p_ij+u_i p_jh+u_j p_hi=0                         (16)
```

for every triple of distinct owners.  Equivalently,

```text
delta_ij=p_ij/(2q u_i u_j)                            (17)
```

is a rational coboundary `delta_ij=z_i-z_j`.

As in THM-3395, only a star needs to be searched.  Choose `p_0i` from its
pair fibre.  Every remaining value is forced by

```text
p_ij=(u_i p_0j-u_j p_0i)/u_0.                        (18)
```

The star is compatible exactly when every forced value is integral and lies
in its fibre `(15)`.

## 5. Exact all-`q` cover theorem

Let `U` be any finite set of positive speeds.  Then

```text
there exists t with union_(u in U) D_(q,u)(t)=Z/qZ

iff

there are selected owners I subset U and modes (r_i,s_i), i in I,
whose blocks B_i cover Z/qZ and which admit a complete mode cochain.
                                                               (19)
```

More precisely, the full-cover locus itself is the union, over the block
covers in `(19)`, of the intersections of their mode sets `(7)`.

### Necessity

Take a full-cover time `t` and retain the owners that fire.  Section 3 turns
each actual nonempty firing set into one maximal consecutive mode block; those
blocks cover all sheets.  Choose the centre tooth in `(13)` whose interval
contains a real lift of `t`.  Formula `(14)` then gives `(15)`.  Actual centre
differences telescope, giving `(16)`.

### Sufficiency: closure, CRT, then Helly

Suppose the right side of `(19)` holds.  Closure `(16)` gives rational
potentials `z_i` with `(17)`.  Put

```text
L_i=(Z+h_i/(2q))/u_i.                                 (20)
```

The congruence in `(15)` is equivalent to

```text
delta_ij in L_i-L_j.                                  (21)
```

Consequently the rational lattice cosets `L_i-z_i` meet pairwise.  Multiply
all cosets by an integer divisible by every `u_i` and every offset
denominator.  They become ordinary integer congruence classes.  Pairwise
generalized-CRT compatibility is sufficient for a simultaneous solution, so
there is one shift `c` with

```text
x_i=z_i+c in L_i                                      (22)
```

for every owner.

The strict inequalities in `(15)` make the real open intervals centred at
the `x_i` pairwise intersect.  One-dimensional Helly gives a common real time
`t`.  By `(7)` that time fires every selected block, and their union covers
all sheets.  Reducing `t` modulo one proves the left side of `(19)`.

The lift is already coherent before passing to the circle.  Therefore there
is no circular-wrap exception and no owner-count restriction.  In
particular, when `q<=7`, every `m<=7` has only `s=1`; `(19)` recovers the
one-kernel-coset theorem of THM-3395 and extends its analytic statement from
six owners to an arbitrary finite owner set.

## 6. Boundary, positive, and hostile controls

### First domino and a zero-tie cover

At `q=8`, the four odd speeds `(1,3,5,7)` have domino modes

```text
1 -> {0,1},      3 -> {3,6},
5 -> {2,7},      7 -> {4,5}.                          (23)
```

They partition the sheets.  Their centre lattices all contain `t=-1/16`, so

```text
p_ij=0                                             for every pair. (24)
```

This is a complete tie, not a tournament.  A sign orientation would add
gauge while discarding the zero gap that makes the cover work.

### Pairwise overlap is still not closure

At `q=8`, use modes

```text
12 -> {0,2,4,6},      10 -> {1,5},      14 -> {3,7}.  (25)
```

They partition the sheets as `4+2+2`, and all three pair fibres are

```text
P_01=P_02=P_12={-16,16}.                              (26)
```

But closure would require

```text
14p_01+12p_12-10p_02=0,                               (27)
```

which no choice of signs in `(26)` satisfies.  Thus block cover plus every
pair overlap is not sufficient; the integral cycle class is load-bearing.

### A new `q=9` mixed-mode zero tie

At `q=9`, the speeds `(6,1,5,7)` have modes

```text
6 -> {0,3,6},       1 -> {1,2},
5 -> {5,7},         7 -> {4,8}.                       (28)
```

At `t=5/6` these blocks partition all nine sheets and all six gaps are zero.
This is the first exact mode beyond the `q=8` domino census: a triple kernel
coset and three dominoes glue through the same centre.  It illustrates why
the general carrier needs both gcd-type and block-length coordinates.

## 7. Pure dilation is an exact certificate-preserving action

Let `d>=1` be an integer and transform one sheet problem by

```text
(q,u,t) |-> (dq,du,t/d).                              (29)
```

For a sheet `k in Z/dqZ`,

```text
(du)(t/d+k/(dq))=u(t+k/q).                            (30)
```

Thus a degree-`dq` block is exactly the inverse image of the old block under
`Z/dqZ -> Z/qZ`.  The phase count `m` and phase step `a` are unchanged, while

```text
g' = dg,       h'=dh mod 2dq,       w'=dw,
centre lattice and radius are divided by d.           (31)
```

For corresponding centre lifts `x_i'=x_i/d`, the cochain scales as

```text
p_ij'=2dq(du_i)(du_j)(x_i'-x_j')=d^2 p_ij.            (32)
```

The congruence modulus and both sides of the overlap inequality `(15)` also
scale by `d^2`; closure `(16)` scales by `d^3`.  Hence block cover, complete
cochain, and physical realization correspond in both directions on the pure
`d`-multiple stratum.  This is an exact certificate-preserving dilation law
with its mode sidecar retained, not a similarity inferred from matching
counts.

The independent finite-mode artifact in the frontmatter verifies the first
three nontrivial controls:

```text
q=4 -> 8:    (1,3,5,7)   -> (2,6,10,14),
q=5 -> 10:   (1,2,3,4,7) -> (2,4,6,8,14),
q=6 -> 12:   (2,3,5,7)   -> (4,6,10,14).              (33)
```

Mixed gcd strata need not descend: forgetting the mode destroys the
sidecar that distinguishes them.

## 8. Verification and scope

Three independently developed finite artifacts now supply convergence
controls without entering the proof dependencies.  The `q=8,...,14` compiler
reconstructs the literal THM-3387 slices by independent event and mode routes,
with **rank-at-most-five** minimal-edge counts `(32,22,18,0,8,0,0)`, exact row
vector `(1152,1205,1269,1287,1271,1287,1287)`, and no new core rescue.  The
`q=15` boundary model finds no cover through rank five, exactly `157` minimal
rank-six edges, and `155` that require at least one three-phase mode.  For its
fixed canonical edge `{1,2,3,4,5,7}`, independent all-witness minimization
confirms that all four unit owners need trimodes; two other rank-six edges are
domino-sufficient.  A mode-independent full-subset audit finds higher-rank
physical edges at every `q=8,...,15`, correcting the former global-looking
profile labels while leaving all low-rank counts and row consequences intact.
All three artifacts are **FINITE-EXACT** and outside any new LRC(14) ledger
decrement.

The standard-library theorem companion performs the following exact checks.

- For `q=2,...,29` and `1<=u<=2q`, it constructs all `42,328` modes of
  `868` owner types, checks all `2,444` wrapped blocks and strict endpoints,
  and recovers the unique maximal mode on all `47,200` boundary/mid-cell
  phase samples.  This reaches the first modes of lengths three, four, and
  five at `q=15,22,29`.
- It compares `(15)` with an independent rational centre-lattice distance on
  `40,187` mode pairs.
- It compares mode-cochain search with an independent exact event sweep on
  `1,669` bounded subsets.  The global universes use speeds through
  `min(q+3,10)` and ranks through five for `q<=8`, then speeds through eight
  and ranks through four for `q=9,10,11`; `569` capacity-relevant subsets are
  swept cell by cell.
- It reconstructs a source time from each positive cochain by the generalized
  CRT and real interval Helly, and freezes `(23)`--`(28)` plus the `m=1`
  all-sheet control.

An independent audit re-derived the mode interval, wrap uniqueness, `m=1`
case, exact-state necessity, pair congruence, weighted closure, star formula,
generalized CRT, and open-Helly sufficiency.  It found no mismatch in `3,520`
independent fixed-mode event/cochain cases through `q=12`; checked `227,430`
grid-event memberships through `m=35` and another `665,484` at the first
length-`6,...,11` regimes; checked `35,748` sheet modes through `q=25`; and
exhausted `50,625` four-coset CRT instances.  The audit also reproduced all
hashes, normal/optimized output, strict controls, and the security boundary.

There are no floating literals or optimization-dependent assertions.
Reproduce with

```text
python 04-computation/lrc_general_finite_mode_sheet_cover_cochain_thm3398.py
python -O 04-computation/lrc_general_finite_mode_sheet_cover_cochain_thm3398.py
```

Both runs LF-normalized-byte-match the stored output.

The theorem exactly characterizes the cyclic sheet-cover carrier.  It does
not classify every finite mode clutter, decide whether a transverse cover is
rescued by a descended core clock, transport owner current into the refined
LRC ledger, give a ledger decrement, or prove LRC(14).
