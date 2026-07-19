---
id: THM-1161
title: The local three-tooth spacing factor has sharp universal infimum one
status: PROVED — exact infinite legal counterfamily; the local factor proposed in THM-1160 is refuted, uniform r=5 still open
source: codex-2026-07-18-S75 (exact referee of THM-1160's sparse-gap three-tooth target)
depends_on:
  - THM-1160  # proposed local factor corrected here
  - THM-1140  # legal-killer and four-comb setting
related: [THM-1141, THM-1147, THM-1148]
script: 04-computation/lrc14_three_tooth_sharp_one_codex_20260718.py
output: 05-knowledge/results/lrc14_three_tooth_sharp_one_codex_20260718.out
---

# THM-1161 — the local three-tooth factor is sharply one

THM-1160 proposed that three foreign teeth in a gap of the smallest
killer have a uniform spacing advantage greater than `1.295` over its
equal-split bound.  That statement is false.  In fact the best universal
multiplicative factor available from this one-gap datum (when the crude
equal-split bound is positive) is exactly **one**.

The correction is an exact infinite family, not a limiting numerical search.

## 1. Statement

At threshold `1/14`, write the closed tooth of comb `k` with index `q` as

```text
T(k,q) = [q/k - 1/(14k), q/k + 1/(14k)].                 (1)
```

Let

```text
P = {1,2,3,5,6,7,9,10},
K = 0 (mod 4),                 K >= 132,
(k1,k2,k3,k4) = (K,K+1,K+2,K+3).                        (2)
```

Put `m=K/4` and take the closed hull of the `K`-gap following tooth `m`:

```text
G_K = [A,B]
    = [1/4 + 1/(14K), 1/4 + 13/(14K)].                  (3)
```

The strict safe gap is `(A,B)`; adjoining its two wall endpoints changes no
length.  We use the closed hull only to state endpoint and intersection
identities without repeated limiting notation.

Then all four killers are legal above `P`, the whole gap lies in a
core-safe component, and each of the three foreign combs has exactly one
tooth meeting `G_K`.  Those teeth are disjoint and divide `G_K` into four
pieces

```text
p0 = 3(K-2)  / [28K(K+3)],
p1 = 3(K+6)  / [28(K+2)(K+3)],
p2 = (3K+22) / [28(K+1)(K+2)],
p3 = (3K+26) / [28K(K+1)],                              (4)
```

with

```text
0 < p0 < p1 < p2 < p3.                                  (5)
```

For THM-1160's crude equal-split lower bound

```text
E_K = (6/(7K) - 3/(7(K+1)))/4
    = 3(K+2)/[28K(K+1)],                                (6)
```

the exact factor is

```text
p3/E_K = (3K+26)/(3K+6)
        = 1 + 20/(3K+6)  ->  1.                         (7)
```

Even if “equal split” means the average of the **actual** remaining
length rather than the conservative quantity (6), the factor tends to one.
Thus no uniform factor `c>1` follows from the local three-tooth hypotheses.
Since the elementary longest-piece pigeonhole gives `c=1`, the sharp
universal infimum is exactly one.

## 2. The core and legality are exact

The interval

```text
J = [29/126,27/98]                                      (8)
```

is a core-safe component for `P`.  The following table gives a direct
certificate: the last column lies inside the displayed safe unit sector.

| `p` | sector `n` | `pJ` | containing sector |
|---:|---:|---:|---:|
| 1 | 0 | `[29/126,27/98]` | `[1/14,13/14]` |
| 2 | 0 | `[29/63,27/49]` | `[1/14,13/14]` |
| 3 | 0 | `[29/42,81/98]` | `[1/14,13/14]` |
| 5 | 1 | `[145/126,135/98]` | `[1+1/14,1+13/14]` |
| 6 | 1 | `[29/21,81/49]` | `[1+1/14,1+13/14]` |
| 7 | 1 | `[29/18,27/14]` | `[1+1/14,1+13/14]` |
| 9 | 2 | `[29/14,243/98]` | `[2+1/14,2+13/14]` |
| 10 | 2 | `[145/63,135/49]` | `[2+1/14,2+13/14]` |

The left endpoint is a `p=9` wall and the right endpoint is a `p=7`
wall, so (8) is maximal as well as safe.  Moreover

```text
A > 1/4 > 29/126,
27/98 - (1/4 + 13/(14*132)) = 239/12936 > 0.             (9)
```

The right side of (3) decreases with `K`; hence `G_K` is strictly inside
`J` for every member of (2).  It is therefore strictly core-safe away from
its own killer walls.

The legal-killer condition used in the four-comb reduction is
`k>13 max(P)`.  Here `max(P)=10`, while `K>=132>130`; consequently each of
`K,K+1,K+2,K+3` is legal and the four speeds are distinct from the core.

## 3. Exactly one foreign tooth, with complete wall chronology

Fix `d` in `{1,2,3}` and put `k=K+d`.  A tooth `T(k,q)` meets `G_K` only
if its center lies in the radius-expanded gap.  After multiplying by `k`,
the possible integer `q` must lie strictly between the two bounds

```text
m + d/4 + d/(14K),
m + 1 + d/4 + 13d/(14K).                                (10)
```

For `d=1,2,3` and `K>=132`, the first bound lies strictly between `m` and
`m+1`, while the second lies strictly between `m+1` and `m+2`.  Thus the
only possible index is

```text
q=m+1.                                                    (11)
```

Its endpoints are

```text
L_d = 1/4 + (26-7d)/[28(K+d)],
R_d = 1/4 + (30-7d)/[28(K+d)].                           (12)
```

The complete local wall order is

```text
A < L_3 < R_3 < L_2 < R_2 < L_1 < R_1 < B.              (13)
```

There is no hidden phase or endpoint assumption in (13): its seven
successive differences are exactly

```text
3(K-2)/[28K(K+3)],          1/[7(K+3)],
3(K+6)/[28(K+2)(K+3)],      1/[7(K+2)],
(3K+22)/[28(K+1)(K+2)],     1/[7(K+1)],
(3K+26)/[28K(K+1)].                                  (14)
```

They are all positive already for `K>2`.  Alternating survivor and tooth
lengths in (14) proves containment, disjointness, the exact formulas (4),
and that there is precisely one contributing tooth from each foreign comb.

Finally, direct subtraction gives

```text
p1-p0 = 3(3K+2)/[14K(K+2)(K+3)],
p2-p1 = (5K+24)/[14(K+1)(K+2)(K+3)],
p3-p2 = (5K+26)/[14K(K+1)(K+2)].                         (15)
```

This proves the strict ordering (5), so `p3` really is the longest piece.

## 4. Both meanings of “equal split” tend to the same obstruction

The three foreign teeth have widths `1/[7(K+d)]`.  Their actual surviving
mean is therefore

```text
M_K = (p0+p1+p2+p3)/4
    = (3K^3+24K^2+55K+36)
      / [28K(K+1)(K+2)(K+3)].                            (16)
```

Equations (4) and (16) give

```text
p3/M_K
 = 1 + (17K^2+93K+120)/(3K^3+24K^2+55K+36)
 -> 1.                                                    (17)
```

For an arbitrary local three-tooth configuration, removing at most three
intervals creates at most four surviving components.  Therefore the longest
piece is at least one quarter of the actual surviving length (and at least
the mean if fewer than four pieces occur).  Each foreign tooth has width at
most `1/(7k2)`, so this quarter-length is at least THM-1160's conservative
equal-split quantity.  This proves the universal factor `>=1`.  Equations
(7) and (17) show that the family approaches one under either normalization,
proving sharpness.

At the first legal member `K=132`, already

```text
p3/E_K = 211/201 = 1.049751... < 1.295 < 4/3.             (18)
```

This also corrects a literal conflation in the sparse-gap THM-1160.  The constants `4/3`
and `1.295` are not equivalent: `4/3` is the stronger proposed
nonuniformity factor, whereas `1.295` was a separate smaller target inferred
from the standing sampled row.  The exact family refutes **both** as
universal local factors.  THM-1148 independently refutes the ordinary
`4/3` heuristic in another four-residue configuration.

For perspective, the selected gap itself does not meet the desired
four-comb threshold:

```text
7(K+3)p3 = (K+3)(3K+26)/[4K(K+1)] -> 3/4,                (19)
```

and it equals `9495/11704<1` at `K=132`.  This does **not** refute the
global four-comb theorem: another `K`-gap or another core-safe component may
still supply the required survivor.  It refutes only the claim that every
three-tooth gap, or one chosen without a global selection argument, has a
fixed spacing bonus.

## 5. Exact replay and tournament / carrier audit

The deterministic companion script uses only `fractions.Fraction`.  It
checks the safe-sector table, legality, unique tooth indices, all wall and
piece identities, the first 992 legal family members through `K=4096`, and
the two limiting ratios.  Normal and optimized Python runs are byte-identical.

Frozen hashes:

```text
source  d78cf903e17b8af9f17ad33d05936c53ecf2332b748748cb65797a9d17a21eb6
output  87f8cd278a4eea06699c3db79d7df9468531b538a8279abf243768cf9bfca489
```

The useful tournament is not naturally on the runners.  Take the four
surviving pieces as vertices and use the antisymmetric observable

```text
Delta(i,j) = length(p_j)-length(p_i),
i -> j iff Delta(i,j)>0.                                 (20)
```

Sign reversal is the switch/gauge.  The label order `p0,p1,p2,p3` is the
tie Hamiltonian path (there are no ties in this family).  The tournament is
transitive: score histogram `(0,1,2,3)`, no directed triangles, four
singleton SCCs, one Hamiltonian path, and zero edge flips along the family.
It preserves the local max/mean question but forgets which tooth owns a
boundary.

The assumption challenge identifies the faithful carrier:

- runner vertices retain the four speeds but lose the phase of the chosen
  gap, so they cannot decide the local factor;
- the three tooth centers retain their order but lose both boundary pieces
  and the tooth widths;
- the four-piece quotient retains exactly the factor but loses ownership;
- the eight ordered wall events
  `A,L3,R3,L2,R2,L1,R1,B` retain the full local predicate.

Thus the obstruction is a wall-chronology object.  The next valid route to
uniform `r=5` must compare or select among many such stalks globally; a
single transitive piece tournament contains no further dispersion debt.

## Consequence for the frontier

The conditional `m=2` reduction in THM-1160 remains valid.  Its 800-row
absence-of-`m=2` observation remains a finite experiment.  The proposed
universal `m=3` spacing rescue is now exactly dead.  Uniform `r=5` and the
global four-comb theorem remain open.
