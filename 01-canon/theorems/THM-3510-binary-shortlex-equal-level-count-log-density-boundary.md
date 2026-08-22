---
id: THM-3510
title: "Binary shortlex equal level counts do not determine logarithmic density"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In binary shortlex,
  choosing the first or last half of every positive level gives the same
  count 2^(n-1) but logarithmic densities log(3/2)/log(2) and
  log(4/3)/log(2).  Alternating the two choices in stages of length
  r_k=2^(2^k) preserves every level count while forcing those two values as
  distinct complete-stage subsequential limits, so the resulting arbitrary
  language has no logarithmic density.  At depth two the four words, not the
  six pair edges, are the vertices of the lexicographic tournament; the
  three nonzero F_2-linear forms label its three perfect-matching edge
  partitions.  No ancestry or cross-domain transport is asserted.
source: codex/thm3499-boundary-independent-audit/2026-08-16
depends_on:
  - THM-3499-regular-shortlex-languages-have-logarithmic-density
related:
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
script: 04-computation/shortlex_equal_count_log_density_boundary_audit_20260816.py
output: 05-knowledge/results/shortlex_equal_count_log_density_boundary_audit_20260816.out
script_sha256: 54cd8c333d5fa2bf03a9f0742b9a42f19c19bc6de9ff8e690c724d07a6f4ff6e
output_sha256: 84d33c86d6478ea4a1506ef75d001358c49a5510110633ee9ea5d89375cbda21
exact_ledger_sha256: e0f67a9bdda9e6e158969c54c023b390c17fed8c4ce68b0ddeffee4936824513
hash_basis: LF-normalized UTF-8 bytes for files; exact rational semantic ledger
---

# THM-3510 -- the equal-count boundary for binary shortlex density

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and exact indexing

Use the shortlex convention of THM-3499: the empty binary word has index
one.  The words of length `n` therefore occupy exactly

```text
I_n=[2^n,2^(n+1)-1].                                    (1)
```

For every `n>=1`, split this interval into

```text
I_n^L=[2^n,3*2^(n-1)-1],
I_n^R=[3*2^(n-1),2^(n+1)-1].                            (2)
```

Both intervals have exactly `2^(n-1)` integers.  Let `A_L` choose `I_n^L`
on every positive level and let `A_R` choose `I_n^R` on every positive
level.  Then both logarithmic densities exist and

```text
delta_log(A_L)=log(3/2)/log 2,
delta_log(A_R)=log(4/3)/log 2.                           (3)
```

Thus the complete level-count sequence does not determine the logarithmic
density, even among regular binary prefix languages.

More sharply, put

```text
r_k=2^(2^k)                         (k>=0),
R_(-1)=0,       R_k=sum_(j=0)^k r_j.                    (4)
```

Stage `k` consists of the positive word levels

```text
R_(k-1)+1,...,R_k.                                      (5)
```

Define `A_osc` by choosing `I_n^L` throughout every even stage and `I_n^R`
throughout every odd stage.  It still has

```text
a_n=|A_osc intersect I_n|=2^(n-1)       for every n>=1. (6)
```

At the complete-stage cutoffs

```text
N_k=2^(R_k+1)-1,                                         (7)
```

the normalized harmonic sums have the two limits

```text
lim_(k->infinity, k even)  1/log(N_k) sum_(m<=N_k,m in A_osc) 1/m
  =log(3/2)/log 2,

lim_(k->infinity, k odd)   1/log(N_k) sum_(m<=N_k,m in A_osc) 1/m
  =log(4/3)/log 2.                                      (8)
```

The two values are distinct, so `A_osc` has no logarithmic density.  This is
the sharp arbitrary-language boundary of THM-3499's regular-language
existence theorem.

## 2. One-level harmonic masses

Write

```text
F_n=sum_(m in I_n^L) 1/m,
G_n=sum_(m in I_n^R) 1/m.                               (9)
```

For `h=2^-n` and `f(x)=1/(1+x)`, the endpoints in (2) give the exact left
Riemann sums

```text
F_n=h sum_(j=0)^(2^(n-1)-1) f(jh),
G_n=h sum_(j=2^(n-1))^(2^n-1) f(jh).                   (10)
```

Consequently

```text
F_n -> integral_0^(1/2) dx/(1+x)=log(3/2),
G_n -> integral_(1/2)^1 dx/(1+x)=log(4/3).             (11)
```

The endpoint convention in (2) is load-bearing: the first block stops at
`3*2^(n-1)-1`, while the second starts at `3*2^(n-1)`.  There is neither a
gap nor an overlap, and each count is exactly `2^(n-1)`.

The decreasing function `f` also gives summable one-sided errors.  On each
mesh interval, bound the integral between the left and right rectangles and
sum the resulting telescoping differences.  This yields

```text
0 <= F_n-log(3/2) <= 2^-n(f(0)-f(1/2)) = 2^-n/3,
0 <= G_n-log(4/3) <= 2^-n(f(1/2)-f(1)) = 2^-n/6.        (12)
```

In particular, for any choice of a left or right half on every positive
level, replacing the exact level masses by their corresponding limits incurs
a total error bounded independently of the number of completed levels; the
coarse uniform bound `sum_(n>=1)2^-n/3=1/3` suffices.

## 3. Fixed choices and arbitrary cutoffs

The final integer through level `D` is

```text
M_D=2^(D+1)-1.                                          (13)
```

Equations (11)--(12) imply

```text
sum_(m<=M_D,m in A_L)1/m=D log(3/2)+O(1),
sum_(m<=M_D,m in A_R)1/m=D log(4/3)+O(1).               (14)
```

Also

```text
log M_D=(D+1)log 2+log(1-2^(-(D+1)))=(D+1)log 2+o(1).  (15)
```

Dividing (14) by (15) proves (3) at complete-level endpoints.  A partial
final level has total harmonic mass at most

```text
sum_(m=2^n)^(2^(n+1)-1)1/m <= 2^n/2^n=1.               (16)
```

Thus it changes the numerator by `O(1)`, while its cutoff has logarithm
`n log 2+O(1)`.  The same limits therefore hold at arbitrary integer
cutoffs.  Inclusion or exclusion of the empty word likewise changes only a
bounded initial term.

## 4. Alternating stages and the normalization audit

Let

```text
mu_L=log(3/2),              mu_R=log(4/3),
epsilon_k=L for k even,     epsilon_k=R for k odd.       (17)
```

The summable error (12) gives the uniform complete-stage formula

```text
sum_(m<=N_k,m in A_osc)1/m
  =sum_(j=0)^k r_j mu_(epsilon_j)+O(1).                 (18)
```

The newest stage dominates every earlier stage.  Indeed,

```text
R_(k-1)/r_k
 <= k r_(k-1)/r_k
 =k*2^(-2^(k-1)) ->0.                                  (19)
```

Subtracting `mu_(epsilon_k) R_k` from (18), the contribution of all earlier
stages is at most

```text
|mu_L-mu_R| R_(k-1)+O(1)=o(R_k).                        (20)
```

Hence the numerator in (18), divided by `R_k`, tends to `mu_L` along even
`k` and to `mu_R` along odd `k`.

The exact endpoint (7) supplies the final off-by-one check:

```text
log N_k
 =(R_k+1)log 2+log(1-2^(-(R_k+1))),
log N_k/R_k -> log 2.                                  (21)
```

Combining (18)--(21) proves (8).  Normalizing by `R_k log 2` would give the
same limit, but (21) verifies the actual logarithmic-density denominator.

## 5. Depth-two K4 typing

Identify the four depth-two binary words with

```text
V=F_2^2={00,01,10,11}.                                  (22)
```

There are `binom(4,2)=6` unordered comparisons.  Orienting each comparison
from the lexicographically earlier word to the later word gives the
transitive tournament on the four word vertices, with exactly six arcs.
The orientation here is the declared lexicographic order; it is not supplied
by harmonic mass.

The three nonzero linear forms on `F_2^2` have two fibres of size two.  The
two within-fibre edges form the three perfect matchings

```text
x_1:       00--01 | 10--11,
x_2:       00--10 | 01--11,
x_1+x_2:   00--11 | 01--10.                             (23)
```

They partition the six edges of `K4`.  Equivalently, for every distinct
`u,v`, the nonzero difference `u-v` has a unique nonzero annihilating
covector, so the edge `{u,v}` occurs in exactly one row of (23).

The six edges are therefore comparisons or edge-objects, not six vertices
of this tournament.  If they were promoted to a separate six-vertex carrier,
a tournament would require `binom(6,2)=15` pairwise comparisons and an extra
orientation.  THM-2606 and THM-2753 describe the canonical association and
matching structure on those six edge-objects; no automatic six-vertex
tournament is inherited here.

## 6. Independent exact hostile

The companion is independent of both THM-3499 programs.  It uses exact
`Fraction` arithmetic to verify:

- all binary level endpoints, half seams, equal counts, and harmonic masses
  through level nine;
- exact rational enclosures of `log(3/2)`, `log(4/3)`, and `log 2` from the
  positive `2*atanh(z)` series, including a rational tail bound;
- complete-stage normalized rational intervals through `k=8`, with the
  late even and odd intervals disjoint and certified near the two targets;
- all four depth-two words, six lexicographic arcs, three nonzero covectors,
  their three `2+2` perfect matchings, and exact coverage of `E(K4)`; and
- the hostile count `binom(6,2)=15` for a genuinely six-vertex tournament.

Run

```bash
python3 04-computation/shortlex_equal_count_log_density_boundary_audit_20260816.py
python3 -O 04-computation/shortlex_equal_count_log_density_boundary_audit_20260816.py
```

Both runs line-match the stored output and end in `FAILED CHECKS: NONE`.
The exact semantic ledger is

```text
e0f67a9bdda9e6e158969c54c023b390c17fed8c4ce68b0ddeffee4936824513.
```

## 7. Scope and loss ledger

```text
source object:       binary words with their positive shortlex levels;
target object:       a subset of the positive integers with harmonic weights;
map:                 word -> its shortlex index;
preserved:           level and within-level left/right address;
forgotten by a_n:    every within-level address;
needed sidecar:      the selected rank interval on each level;
decisive hostile:    superdominant alternating left/right stages.           (24)
```

The result concerns only arbitrary harmonic subsets at the boundary of
THM-3499.  It does not recover word ancestry, attach a physical current,
transport the finite `K4` bookkeeping to another problem, or imply any LRC,
Jacobian, or automatic tournament statement.

QED.
