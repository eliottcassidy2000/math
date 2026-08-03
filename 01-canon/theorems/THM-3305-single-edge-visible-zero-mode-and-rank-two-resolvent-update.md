---
id: THM-3305
title: "Single-edge visible zero mode and rank-two resolvent update"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In THM-3287's static
  maximizing-witness graph over THM-3277's selected eleven-edge tree, add
  each of the eleven remaining full-core row edges separately.  The complete
  signature census is the table below.  Exactly the edge `(7,18)`, carrying
  `Q-7 <-> Q+1`, makes the pre-existing adjacency kernel visible to the
  all-ones observer: its mass is exactly `3/5`, its observer/Hankel order is
  `15`, and its even/odd half-orders are `8/7`.  The mechanism is the local
  kernel relation `3 e_(18,Q+1) - sum_(six leaves) e_leaf`, and the new
  adjacency is a literal rank-two update.  Its complete walk series follows
  from an exact `2 x 2` Woodbury correction.  This is a theorem about a
  finite static compatibility graph, not chronological response composition,
  a tournament, or an FC/GMC/LRC mechanism.
audit: >
  The primary companion transitively replays the frozen THM-3287 relation
  source, exhausts all eleven single-edge additions, derives every full and
  half-transfer Krylov order, computes exact orthogonal kernel projections,
  proves the localized head and symbolic Woodbury identity, checks Hankel
  minimality, and carries active-node birth controls.  The independent audit
  imports or executes neither the primary companion nor the intermediate
  half-transfer scout: it pins and replays only THM-3287's theorem, relation
  source and frozen output, then separately rebuilds all graphs, recurrences,
  projections, the three-versus-six incidence certificate, the rational
  series, and the reduced rank-two correction.  Normal and optimized modes
  byte-match both stored outputs; both sources have zero assertion nodes and
  zero floating-point literals.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut
  - THM-3288-maximizing-witness-lifted-walk-rational-series
related:
  - THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas
  - THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary
script: 04-computation/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py
output: 05-knowledge/results/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.out
script_sha256: 5704ee8d6641a12e50d6a6179ccf58f911ba782cfe369ac70a32db51887565b2
output_sha256: b52db8c18d3182fb79ed5a230d9f75775155f438c7f32e974ab1ff3fcb15e556
independent_audit_script: 04-computation/gmc_single_edge_zero_mode_resolvent_update_independent_audit_20260803.py
independent_audit_output: 05-knowledge/results/gmc_single_edge_zero_mode_resolvent_update_independent_audit_20260803.out
independent_audit_script_sha256: f7f7fdfb32e26eff7d64d6e654916111d3b227122eaf3d69ba7b7501dd4cdb43
independent_audit_output_sha256: 2fca320ccbf92fcf47a4a5eda9ed7c22c07a02843599c50eef4a29bdffd63acf
hash_basis: LF-normalized bytes
---

# THM-3305 -- single-edge visible zero mode and rank-two resolvent update

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem isolates the smallest edge update in the static relation-walk
system of THM-3287/3288 that creates an observer-visible zero mode.  The
result is exhaustive in its declared eleven-edge universe and gives an exact
closed-form update for the resulting infinite sequence.

## 1. Graph, observer, and exhaustive universe

Let `T_*` be THM-3277's selected eleven-edge row tree and let `C` be
THM-3287's twenty-two-edge response core.  For a row-edge family `E`, form the
static decorated graph `Gamma(E)` as in THM-3288: for every maximizing-witness
pair

```text
(s,t) in R_(u->v),
```

retain the active vertices `(u,s)` and `(v,t)` and insert both directed
adjacency entries between them.  Write `A_E` for the resulting symmetric
zero-one adjacency matrix and `1_E` for its all-ones column.

For a matrix `A`, define

```text
o(A,w) = dim_Q span{w,Aw,A^2w,...},
P_0    = orthogonal projection onto ker(A),
mu_0   = 1^T P_0 1.
```

The exact observer signature is

```text
Sigma(A)=(active vertices, directed arrows,
          o(A,1), o(A^2,1), o(A^2,A1), dim ker(A), mu_0).          (1)
```

The selected-tree control is

```text
Sigma(A_T*)=(23,40,14,7,7,9,0).                                  (2)
```

There are exactly eleven row edges in `C \ T_*`.  Adding each one separately
to `T_*` gives the complete census below.  `Birth` records an active decorated
vertex absent from the selected-tree graph.

| added edge | maximizing relation | birth | exact signature `(1)` |
|---|---|---|---|
| `(2,10)` | `Q+4 <-> Q-1` | none | `(23,42,14,7,7,9,0)` |
| `(3,22)` | `Q-1 <-> Q+4` | none | `(23,42,14,7,7,9,0)` |
| `(7,18)` | `Q-7 <-> Q+1` | none | `(23,42,15,8,7,9,3/5)` |
| `(10,16)` | `Q-8 <-> Q+1` | `(10,Q-8)` | `(24,42,14,7,7,10,0)` |
| `(10,22)` | `Q-1 <-> Q+4` | none | `(23,42,14,7,7,9,0)` |
| `(11,13)` | `Q+1 <-> Q-8` | `(13,Q-8)` | `(24,42,14,7,7,8,0)` |
| `(13,18)` | `Q-1 <-> Q+1,Q+2,Q+4,Q+5` | none | `(23,48,14,7,7,9,0)` |
| `(13,22)` | `Q-1 <-> Q+4` | none | `(23,42,14,7,7,9,0)` |
| `(16,21)` | `Q+1 <-> Q-8` | `(21,Q-8)` | `(24,42,14,7,7,10,0)` |
| `(17,22)` | `Q-1 <-> Q+4` | none | `(23,42,14,7,7,9,0)` |
| `(19,22)` | `Q-1 <-> Q+4` | none | `(23,42,12,6,6,9,0)` |

Thus `(7,18)` is the **unique** one-edge addition in `C \ T_*` with
`mu_0 != 0`.  The table also proves two boundary facts.  Three additions
change the active-start vector by birthing one vertex, while `(19,22)` changes
no vertex and lowers the observer order from `14` to `12`.  Scalar recurrence
order is therefore not monotone under these edge additions and is not
determined by vertex/arrow census alone.

## 2. The three-versus-six zero-mode mechanism

For the exceptional edge, put

```text
u=(7,Q-7),                 c=(18,Q+1),
r=(19,Q-1).
```

Both `u` and `c` are already active in `Gamma(T_*)`.  The new undirected edge
adds exactly the two symmetric entries

```text
e_u e_c^T + e_c e_u^T,                                      (3)
```

which has matrix rank two.  In the updated graph the local neighbour profiles
are

```text
N(c)={u,r},
N(18,Q+2)=N(18,Q+4)=N(18,Q+5)={r},
N(22,Q+1)=N(22,Q+2)=N(22,Q+5)={u}.                         (4)
```

Consequently, if `L` is the six-element set of sibling vertices in `(4)`,

```text
v = 3 e_c - sum_(ell in L) e_ell                             (5)
```

satisfies

```text
Av=0,          1^T v=-3,          v^T v=15.                  (6)
```

Exact kernel row reduction shows that every other kernel direction is
orthogonal to `1`; hence

```text
P_0 1 = -v/5,                  mu_0=1^T P_0 1=3/5.           (7)
```

This is not an observer-invisible twin difference: the primitive kernel
coefficients have nonzero sum.  The recurrence head is the integral multiple

```text
h=36v.                                                               (8)
```

More precisely, with

```text
p(t)=(t-2)(t^2-5t+3)(t^4-14t^3+62t^2-94t+30),               (9)
```

one has

```text
h=p(A^2)1,       Ah=0,       1^T h=-108,
h^T h=19440,     p(0)=-180,  h/p(0)=P_0 1.                  (10)
```

The support of `h` is exactly `108` at `c` and `-36` at the six siblings in
`L`.  Thus the visible zero mode is completely localized by `(4)`.

## 3. Exact sequence and rational closed form

For the updated adjacency `A`, define

```text
a_n=1^T A^n 1,                       n>=0.                   (11)
```

The sequence begins

```text
23, 42, 116, 252, 692, 1602, 4404, 10518, 28966,
70344, 194026, 475722, 1313706, 3239490, 8953208,
22155444, 61266026, 151938198.                                (12)
```

Its minimal prefix-valid characteristic polynomial is

```text
q(z)=z p(z^2)
    =z(z^2-2)(z^4-5z^2+3)
       (z^8-14z^6+62z^4-94z^2+30).                           (13)
```

The exact Hankel rank and observer order are both `15`.  Nevertheless the
rational generating function has denominator degree `14`, because the
zero-eigenvalue component contributes the initial atom `(10)` but no positive
power of `A`:

```text
G_A(x)=sum_(n>=0) a_n x^n = N_A(x)/D_A(x),                    (14)

D_A(x)=(1-2x^2)(1-5x^2+3x^4)
       (1-14x^2+62x^4-94x^6+30x^8),                          (15)

N_A(x)=23+42x-367x^2-630x^3+2235x^4+3576x^5
       -6528x^6-9600x^7+9436x^8+12456x^9
       -6284x^10-7032x^11+1668x^12+1368x^13-108x^14.         (16)
```

The degree-seven recurrence from `p` fails on the even subsequence only at
its first possible index, by `-108`, and holds thereafter.  This scalar
boundary debt is exactly the head `(10)`.

## 4. Exact rank-two resolvent update

Let `B=A_T*`, let `U=[e_u,e_c]`, and let

```text
C=[[0,1],[1,0]],                  A=B+U C U^T,
R=(I-xB)^(-1),                    s=U^T R 1,
K=U^T R U.                                                   (17)
```

The Woodbury identity gives the complete infinite-sequence update through a
`2 x 2` inverse:

```text
G_A(x)-G_B(x)=x s^T (C-xK)^(-1) s.                            (18)
```

This is verified symbolically in `Q(x)`, independently of any truncated
sequence.  Put

```text
d_old=1-13x^2+54x^4-77x^6+24x^8,
d_new=1-14x^2+62x^4-94x^6+30x^8.                             (19)
```

Every entry of `s` and `K` has reduced denominator `d_old`, while `(18)` is
the reduced fraction

```text
G_A-G_B = 2x W(x)/(d_old d_new),                              (20)

W(x)=216x^15+144x^14-1101x^13-750x^12+1549x^11+1368x^10
     -120x^9-1086x^8-725x^7+358x^6+409x^5-28x^4
     -84x^3-7x^2+6x+1.                                      (21)
```

The numerator in `(20)` is coprime to its denominator.  The two remaining
observer factors pass through unchanged:

```text
D_B=(1-2x^2)(1-5x^2+3x^4)d_old,
D_A=(1-2x^2)(1-5x^2+3x^4)d_new.                              (22)
```

Thus one node-preserving compatibility edge changes exactly one local secular
factor, creates the visible head `(10)`, and updates the full closed form by a
two-dimensional calculation.

## 5. Scope and proof audit

The quantifier in this theorem is exactly:

```text
for each one-edge family T_* union {e}, with e in C \ T_*.
```

It does not cover arbitrary row edges, simultaneous or sequential additions,
or a graph whose active vertex convention differs from THM-3288.  The arrows
are symmetric static maximizing-witness compatibilities.  They are not the
one-pole chronological transitions excluded by THM-3287, and they have no
intrinsic orientation gauge, so they do not form a tournament.  Nothing here
proves an FC, GMC, LRC, or positivity statement.

The primary companion reconstructs all eleven cases and verifies the
relations, active-node changes, Krylov closures, exact kernel projections,
Hankel ranks, rational functions, local incidence identity, and Woodbury
formula.  The independent audit starts only from THM-3287's pinned relation
artifact and repeats those calculations with separate graph, Krylov,
projection, rational-observer, and resolvent routines.  It imports or executes
neither the primary companion nor its half-transfer predecessor.  Both paths
give the table and equations `(3)--(22)` exactly.

Reproduce with

```text
python3 04-computation/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py
python3 -O 04-computation/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py
python3 04-computation/gmc_single_edge_zero_mode_resolvent_update_independent_audit_20260803.py
python3 -O 04-computation/gmc_single_edge_zero_mode_resolvent_update_independent_audit_20260803.py
```

Each mode byte-matches its declared frozen output.  All arithmetic is exact;
both sources have zero assertion nodes and zero floating-point literals.

**QED.**
