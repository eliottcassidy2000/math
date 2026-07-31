---
id: THM-2966
title: "Spine normal form for critical-run fair extractors"
status: >
  PROVED + VERIFIED-EXACT. Let T(n) >= n+1 be a deadline envelope. An
  exactly fair biased-coin extractor whose pathwise stop time is at most
  T(n) on critical value n exists if and only if there are integer vectors
  0 <= w_{m,k} <= binom(d_m,k), 0 <= v_{m,k} <= binom(d_m,k), with
  d_m = T(m)-m-1, whose spine polynomials satisfy
  sum_{m>=1} p^m q W_m(p) + sum_{m>=1} q^m p V_m(p) = 1/2 on (0,1).
  Consequently C* = inf{C: some T(n) <= Cn+D admits an extractor} equals
  1 + gamma*, where gamma* is the minimal linear degree-growth rate in the
  identity. The parity corollaries: every non-dyadic budget row forces
  half-integer deficits; a half-deficit at (z,o) dyadically split-jumps to
  (z+2^a,o) and (z,o+2^a) mod integer absorptions; ratio-2 shells are
  exactly the case where the two forced corner deficits collide and cancel.
source: death-star-2026-07-30-coinC
depends_on: []
related:
  - THM-2160
  - THM-2225
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_spine_normal_form_referee_thm2966.py
output: 05-knowledge/results/amm12592_spine_normal_form_referee_thm2966.out
---

# THM-2966 -- spine normal form for critical-run fair extractors

Bits are independent with `P(0)=p`, `P(1)=q=1-p`, `0<p<1` unknown. For a
nonconstant stream, the critical value `n` is the length of the maximal
constant initial run. An extractor is a deterministic stopping rule with
labels in `{H,T}`, pathwise total on nonconstant streams, with
`P_p(H)=1/2` for every `p`. Fix an envelope `T(n)>=n+1` and demand the
stop time be `<= T(n)` on every stream of critical value `n`.

## 1. Statement

Define relative depths

```text
d_m = T(m) - m - 1 >= 0.                                     (1)
```

**Theorem.** A `T`-deadline fair extractor exists iff there exist integers

```text
0 <= w_{m,k} <= binom(d_m,k),   0 <= v_{m,k} <= binom(d_m,k)  (2)
```

such that, with

```text
W_m(p) = sum_{k=0}^{d_m} w_{m,k} p^(d_m-k) q^k,
V_m(p) = sum_{k=0}^{d_m} v_{m,k} q^(d_m-k) p^k,               (3)
```

the spine identity holds:

```text
sum_{m>=1} p^m q W_m(p) + sum_{m>=1} q^m p V_m(p) = 1/2       (S)
```

for every `p in (0,1)`.

## 2. Proof

**Extractor implies (S).** For a finite word `x` let `G_x(p)` be the
probability that the eventual output is `H` given that the stream extends
`x`. The chain rule gives `G_x = p G_{x0} + q G_{x1}`, and `G_x` is the
constant `0` or `1` whenever a prefix of `x` has stopped. Every extension
of `0^m 1` has critical value exactly `m`, so all its continuations stop
by absolute length `T(m)`; hence `G_{0^m 1}` is a decided-tree H-share of
relative depth `<= d_m`. Refining every stopped node to exact relative
depth `d_m` (labels constant on subtrees) writes `G_{0^m 1}` in the form
(3) with `w_{m,k}` = the number of relative depth-`d_m` words of weight
`k` lying in `H`-labelled subtrees, which satisfies (2). Iterating the
chain rule along both constant spines,

```text
G_root = sum_{m=1}^{N} [p^m q G_{0^m 1} + q^m p G_{1^m 0}]
         + p^(N+1) G_{0^(N+1)} + q^(N+1) G_{1^(N+1)},          (4)
```

and the boundary terms vanish as `N -> infinity` since `G` is bounded.
`G_root = 1/2` is exactly fairness. This yields (S).

**(S) implies extractor.** Given (2)-(3), stop as follows: upon the first
alternation, which reveals the critical value `m` and the first symbol,
read through absolute flip `T(m)`; order the `binom(d_m,k)` relative
continuation words of weight `k` lexicographically and label the first
`w_{m,k}` (resp. `v_{m,k}`) of them `H`, the rest `T`. The rule is
deterministic, prefix-free, stops at exactly `T(m) = m+1+d_m <= T(n)`
pathwise, and its heads probability is the left side of (S), hence `1/2`
for every `p`. Nonconstant streams are decided almost surely; the two
constant rays are null. QED.

**Budget box is exact.** Any integer vector obeying (2) is realizable
(label any `w_{m,k}` of the weight-`k` leaves `H`), and any decided tree
refines to the full-depth box: the achievable set of `W_m` is exactly the
box. Referee check A verifies both inclusions exhaustively for depths
`<= 3`.

## 3. Consequences

1. **The minimal-C problem is a degree-growth problem.** With
   `T(n) = Cn + D`, `d_m = (C-1)m + D - 1`: an envelope with slope `C` is
   feasible iff (S) is solvable with degree growth `gamma = C-1`. Hence
   `C* = 1 + gamma*` (HYP-9061's question (Q)).
2. **Deficit flow.** Writing `W_m = 1/2 + Delta_m`, (S) becomes
   `sum p^m q Delta_m + sum q^m p Delta'_m = 0`. Cells `(m,k)` with
   `binom(d_m,k)` odd have half-integer deficit coordinates; all others
   integer.
3. **Dyadic split-jump.** For `Delta = 2^a`, over `F_2`,
   `(1+u)^N = (1+u)^(N-Delta) (1+u^Delta)`, i.e.
   `binom(N,K) = binom(N-Delta,K) + binom(N-Delta,K-Delta) mod 2`
   (referee check C, with non-dyadic hostile failures frozen): pushing a
   half-deficit at `(z,o)` down `Delta` levels leaves half-deficits at
   `(z+Delta,o)` and `(z,o+Delta)` only, all interior binomial mass being
   integer (absorbable by in-cone budget cells, which is where the
   quantitative rate problem lives).
4. **Why ratio two is the classical answer.** A sharp shell
   `[n, n+l)` with dyadic tail `l=2^a` has all interior tail classes even
   and exactly one forced corner deficit `+-(1/2) p^n q^l` (0-side; mirror
   on the 1-side). At ratio two (`l=n`) the two corners land on the same
   monomial `(n,n)` and cancel -- this is THM-2160's middle pair
   `0^h 1^h / 1^h 0^h` seen structurally. For any ratio below two the
   corner monomials are pairwise distinct (referee check D) and must be
   routed by split-jumps within budget: the open rate gate of HYP-9061.
5. **Checksum boundary shares.** For the THM-2225 extractor the spine
   partial sums close in an exact polynomial identity
   `sum_{m<M} [p^m q W_m + q^m p V_m] = 1/2 - (p^M+q^M)/2` for dyadic `M`
   (referee check B): the undecided spine states at dyadic boundaries
   carry conditional H-share exactly `1/2`.

## 4. Scope

This is a normal form, not a resolution of (Q): it neither proves nor
refutes `C* < 2`. The identity (S) is over the closed infinite sum; all
truncation arguments must carry the boundary terms of (4). The referee's
checks are finite instantiations of each load-bearing step. QED.
