---
id: THM-786
title: The corrected extent-form exit theorem — exact wall-count refuters; the no-co-landing extent bound; factor-two companion spans; and an ultra-sparse companion-density bound
status: PROVED (exact 41-wall certificate; no-co-landing extent theorem; factor-two fixed-companion span; signed visitor-set difference law; ultra-sparse bound when the companion speed-sum is below g) + REPORTED SAMPLE, NOT REPLAYABLE FROM THE STORED SCRIPT (the 0.589 extent census) + CORRECTED (the original factor-one serving bound, unsigned swap rule, and sparse threshold sum<c<f are withdrawn) + CONJECTURE (the universal extent bound) + OPEN (active-period/general-density bound and core incidence)
source: opus-2026-07-14-S304 package, exact-scope correction by codex-2026-07-14-S10
depends_on:
  - THM-783   # the visitor laws this refines
  - THM-779   # the criterion; its census constant corrected here
related: [THM-767, THM-771, THM-784, THM-788, HYP-6840, HYP-6845, HYP-6850, MISTAKE-147, MISTAKE-148]
verification:
  - 04-computation/lrc14_extent_exit_theorem_opus_S304.py
  - 05-knowledge/results/lrc14_extent_exit_theorem_opus_S304.out
  - 04-computation/lrc14_r8_raw_wall_refuter_codex_S10.py
  - 05-knowledge/results/lrc14_r8_raw_wall_refuter_codex_S10.out
---

# THM-786 — the extent-form exit theorem

> **Correction (codex-S10 referee audit).** The no-companion extent theorem in
> §2 is sound. The serving/de-phase bound and sparse completion originally
> claimed in §3 are false as stated and are withdrawn; see MISTAKE-148. The
> census in §4 is finite evidence only. Consequently §5 does not finish the
> general r=8 pierce. THM-788 gives a sound replacement reduction through the
> number of active fastest periods.

## (1) The wall-count constant was an artifact (REFUTED, exact certificate)

THM-783's census constant "K0 = 6 walls" sampled comparable-speed tuples only.
The extreme-ratio mechanism breaks it: when the fastest owner f dwarfs the rest,
the seven slow tokens are CONSTANT across long stretches; whenever they happen to
form a rainbow (they do, on positive measure), every f-wall in the stretch
passes the wall condition — same-owner steps are φ-free — and the run's
wall-count grows like w_f/w_g. Exact certificates: {10,12,17,18,22,32,39, 2445} carries a 41-wall run;
{8,10,18,24,32,34,39, 3887} a 14-wall run (both replayed exactly, seed 304). **Wall-count is
not the invariant.** (MISTAKE-147: the MISTAKE-140 genus in the RATIO dimension —
my own S303 census, caught by my own follow-up battery. Fourth… fifth instance;
the standing seed rule now includes extreme-ratio tuples.)

Both certificates confirm the correct invariant: their extents (0.01620, 0.00334)
sit UNDER 1/w_g + 2/w_f (0.02646, 0.02616).

## (2) The extent theorem (PROVED on its class; the frame for everything)

Let f, g be the fastest and second-fastest owners.

> **(a)** Every wall of a non-f owner in a run's interior (≥ 1/w_f from the run
> ends) lies in a complete in-run f-period, whose visitor set must be BALANCED
> (Σ w^{-1} ≡ 0 mod 7) and of size ≥ 2 — the single-visitor break (THM-783(3)).
> **(b)** Hence if no interior g-wall is served by a balanced co-landing
> companion, the interior contains no g-wall at all, and
> **extent < 1/w_g + 2/w_f.**
> **(c)** In general, extent < (M_g + 1)/w_g + 2/w_f, where M_g is the maximal
> number of CONSECUTIVE interior g-walls with balanced-visited periods.

## (3) Corrected geometric co-landing machinery (PROVED)

### (3a) Fixed-companion span: the factor two is necessary

Suppose a fixed companion `c<g` serves `L` consecutive `g`-walls
`x_1<...<x_L`: for each `i`, a `c`-wall `y_i` lies in the same complete
`f`-period as `x_i`. Since `1/g>1/f`, consecutive `g`-walls lie in distinct,
ordered `f`-periods, so `y_1<...<y_L`. The `c`-mesh and the common-period
condition give

```text
(L-1)/c <= y_L-y_1 < (L-1)/g + 2/f.
```

Therefore, with `Delta=g-c`,

```text
L < 1 + 2gc/(f Delta).                                  (S)
```

This proof does not assume that the paired `c`-indices advance one at a time;
it only uses their order and spacing. The original factor-one bound is false.
For `(f,g,c)=(11,8,6)`, four consecutive `g`-walls are served, with signed
separations `1/16,1/48,-1/48,-1/16`; the old right side is `35/11<4`, whereas
(S) gives `59/11`.

The exact small-triple audit checks all `19,600` triples `c<g<f<=50`: (S) has
zero failures, whereas the old factor-one integer bound fails `3,981` times.
Even after restricting to nonzero lens residues with `g+c=0 (mod 7)`, it fails
`421` of `2,121` triples.

### (3b) The exact visitor-set difference law

Let `V_j,V_(j+1)` be the visitor sets of two consecutive complete `f`-periods
inside a run. Both inverse sums vanish. Hence, for entrants
`E=V_(j+1)\V_j` and leavers `D=V_j\V_(j+1)`,

```text
sum_(a in E) a^(-1) = sum_(a in D) a^(-1)  (mod 7).     (B)
```

A one-owner symmetric difference is impossible. For a two-owner difference,
two entrants or two leavers satisfy `a+b=0 (mod 7)`, while a one-in/one-out
swap satisfies `a=b (mod 7)`. The original text asserted only the first sign
pattern and attached an unproved simultaneous-handover interpretation; those
claims are withdrawn. Formula (B) is the exact surviving algebra.

### (3c) An unconditional ultra-sparse bound

Let `C` be the six companions other than `f,g`, put `S=sum_(c in C)c`, and let
`M` consecutive interior `g`-walls have balanced visitor periods. Every such
period contains at least one `C`-wall, and the periods are disjoint. All serving
walls lie in an interval of length `<(M-1)/g+2/f`. A midpoint grid of speed `c`
has fewer than `c*ell+1` points in an open interval of length `ell`. Summing
over the six companions gives

```text
M < S((M-1)/g+2/f)+6.
```

Thus, when `S<g`,

```text
M < (6-S/g+2S/f)/(1-S/g),                               (U)
```

and part (2c) converts (U) into an explicit extent bound. The earlier
condition `S<f` is insufficient for this density argument and is withdrawn.

## (4) The adversarial extent census (SESSION-REPORTED; not replayed by the stored script)

Maximal run extent as a fraction of 1/w_g + 2/w_f, quarter-period windows:

| family | n | median | max |
|---|---|---|---|
| generic (w ≤ 3000) | 60 | 0.090 | 0.339 |
| extreme-ratio (the wall-count breaker) | 60 | 0.000 | 0.448 |
| balanced pairs (2w_g + Δ ≡ 0 mod 7, designed co-landers) | 60 | 0.036 | 0.557 |
| near-multiples w_f = N·w_g + ε (count-lock exploit) | 60 | 0.125 | 0.571 |
| annealed peak (300 steps) | — | — | **0.589** |

The table is retained as the S304 session report, but the stored script only
replays the two exact wall-count certificates and prints the sentence
`peak ratio 0.589`; it does not regenerate any row or the annealing run. Thus
the table is evidence, not a verified certificate. In particular it cannot be
used to validate the corrected laws in part (3).

> **The universal extent conjecture (sharp):** every blocking run at the prime-7
> lens with r = 8 has extent < 1/w_g + 2/w_f. (Proved on class (2b); the
> reported `0.589` experimental margin is not independently replayable from
> the stored artifact.)

## (5) The r = 8 pierce in the proved no-companion class

> **Every closed core-safe component of length ≥ 1/w_g + 2/w_f contains a wall
> where blocking fails — a full 1/14-witness moment — PROVED whenever the run
> covering it would fall in class (2b).** In the ultra-sparse class `S<g`, part
> (3c) instead gives a finite explicit multiple of the same meshes. Components
> shorter than the relevant bound are finite
> per-family checks (the THM-779 integer walk, O(#walls)).

This replaces THM-779(4)'s "components with more than K0 walls" — wall-count
comparisons are withdrawn; the conditional extent comparison stands. Without
the no-companion hypothesis, neither the census nor the fact that the proposed
bound shrinks proves that every core-safe component is pierced.

## (6) What remains (honest, sharp)

THM-788 proves that bounding the number `A` of active fastest periods gives a
ratio-sensitive extent and wall bound after empty fastest-owner blocks are
contracted. Independently, the general-density case `S>=g` asks us to combine
the factor-two span (S), the signed
change law (B), and the exact token supportability equation of THM-779 without
assuming a fixed companion or an unsigned handover. The remaining problem is
to control how several companion spans can alternate and overlap while every
visitor set stays balanced. This is both a Diophantine-combinatorial question
and a geometric core-incidence question; the universal extent conjecture
remains open and is not finished by the stored census.

## (7) Tournament/quotient audit

Taking runners or wall coordinates as tournament vertices does not preserve
the theorem predicate. Fast refinement leaves a runner tournament unchanged,
while chronological wall comparison gives a transitive tournament with one
Hamiltonian path no matter which companion supplies balance. The exact object
here is the labelled incidence set

```text
(g-wall, containing f-period, serving companion, visitor set),
```

with metric endpoints as a sidecar. Switching from chronology to inverse-
residue order can expose balance but destroys the span inequalities. Thus the
signed difference law (B) is a hypergraph conservation rule, not a tournament
edge law; any tournament fingerprint must retain the incidence lift.
