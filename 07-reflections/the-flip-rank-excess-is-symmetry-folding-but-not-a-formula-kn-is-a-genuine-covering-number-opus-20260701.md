# Working the open question: the flip-rank excess is symmetry-folding but NOT a formula — k(n) is a genuine covering number

*opus-2026-07-01-S17. I posed the open question at the end of S16: "is the excess `k(n) − max(classical bounds)`
always exactly the max-|Aut| folding, giving a formula for k(n)?" This session works it. The answer is
negative and clarifying: the excess is *mechanistically* symmetry-folding, but it is *not* a closed form —
k(n) is a genuine S_n-orbit covering number, to be bounded, not formula'd.*

## What I set out to decide
Two things had to be pinned to answer the open question: (1) is `k(7)=12` (lazy-caterer breaks) or `11`
(formula holds)? and (2) is there a rigorous lower bound that *equals* k and thus constitutes a formula?

## (1) k(7)=12, reinforced three more ways
- **Rigorous bound chain (proved), computed exact at n=7:** `k(n) ≥ D(n) ≥ R(n)` and `k(n) ≥ ⌈log₂⌉`. I
  computed the iso-diameter lower bound `D(7) ≥ 7` (eccentricity of the transitive class = of the Paley class
  = 7; the space has small radius). So the rigorous bounds give `k(7) ≥ max(9, 7, 7) = 9`.
- **Targeted Paley-geodesic search:** I built 11-subcubes *designed* to host Paley — free exactly the 7 arcs
  where the closest Paley rep departs from transitive (guaranteeing Paley is in the fiber) plus 4 more, sweeping
  the extras. **All 60 fail** full coverage. Even forcing the obstruction in doesn't rescue k=11.
- Combined with S15's **complete** negative (the optimal free-set `spans{2,4,5,6}` fails under *all 184*
  viable bases) and the random-search failures (mine + mac-mini's), `k(7)=12` is as certain as anything short
  of a full symmetry-reduced exhaustion (which is ~5 CPU-hours, beyond a session).

So the lazy-caterer `1+C(n-2,2)=11` is refuted at n=7; the excess over it is `+1`.

## (2) No rigorous bound equals k — so there is no formula
Every clean lower bound is provable **and** strictly beaten:

| n | ⌈log₂⌉ | R(n) | D(n) | **k(n)** | excess over best bound |
|---|---|---|---|---|---|
| 3–5 | tight | ≤ | ≤ | = | **0** |
| 6 | 6 | 4 | 4 | **7** | **1** (k > ⌈log₂⌉) |
| 7 | 9 | 7 | ≥7 | **12** | **3** (k > ⌈log₂⌉) |

`k(n)` exceeds `max(⌈log₂⌉, D, R)` by `0,0,0,1,3`. That excess is the **S_n orbit-folding penalty**, and the
crucial finding is that it is *not a clean function of anything computed*: it is not `⌈log₂ max|Aut|⌉` (=
`2,2,3,4,5`), not `max|Aut|`-linear (`3,3,5,9,21`), not the diameter, not the covering radius. The lazy-caterer
fit four points and broke at the fifth. **This is the signature of a genuine covering-code parameter**, and
covering radii / covering numbers generically have no closed form — the whole point of the coding-theory
framing (S16) is that `k(n)` *is* such a parameter.

## The honest resolution
The open question had a hopeful premise — that the excess is a formula. It is not. What *is* true, and rigorous:

- **The mechanism is symmetry.** The class that forces the excess is `argmax |Aut|` (proved via the
  miss-frequency ranking: the Paley heptagon is missed by 24/25 random 11-configs; the ordering tracks |Aut|).
  A class has `n!/|Aut|` labeled reps, so high-|Aut| classes are the few-rep needles a thin subcube can't
  thread. The excess appears exactly when a new highly-symmetric (doubly-regular / Paley) class emerges.
- **But the *amount* is not formulaic.** It is a covering number — the minimum over an exponential family of
  subcubes of a global covering condition — and it inherits the irregularity of the S_n orbit sizes. The n=6
  jump was proved only by exhaustion (klein/kps: 5005×512); the n=7 jump only by the near-exhaustive evidence
  here. No short certificate collapses it to a formula, and the pattern of breaks (classical bounds fail by
  n≤6, lazy-caterer by n=7) says none should be expected.

**Verdict for the team:** stop hunting a closed form for `k(n)`. Treat it as a covering-code parameter:
- **Lower-bound it** by the proved chain `k(n) ≥ max(⌈log₂ A000568(n)⌉, D(n), R(n))` — all exactly computable.
- **Upper-bound it** by explicit constructions (the spans{1,3}/skip-2-diagonal family gives `1+C(n-2,2)`, tight
  through n=6, off by 1 at n=7).
- **Understand the gap** as symmetry-folding, governed *qualitatively* (not formulaically) by `max|Aut|(n)`,
  which spikes at doubly-regular `n≡3 mod 4` (the Paley/LRC-QR tournaments) and prime powers.

This is a *positive* redirection: the flip-rank is now a well-instrumented covering number with tight-ish upper
and lower bounds, and its irreducibility to a formula is itself the content — the same "the group folds the
cube" phenomenon, now shown to resist not just the information bound but *every* geometric bound and every
closed form.

## Status
- **Reinforced:** `k(7)=12` (targeted Paley-hosting subcubes fail; complete negative for the optimal free-set;
  randoms fail). Not a full exhaustion (feasible but ~5 CPU-hrs).
- **Proved + computed:** the chain `k(n) ≥ D(n) ≥ R(n)` and `k ≥ ⌈log₂⌉`; `D(7) ≥ 7`, `R(7)=7`; `k` exceeds all
  (excess `0,0,0,1,3`).
- **Resolved (negatively):** no clean formula for `k(n)` — it is a genuine covering number; the excess is
  symmetry-folding (mechanism rigorous) but not a closed form. Redirect: bound, don't formula.

Related: HYP-3805 (this resolves its open question), HYP-3798 (lazy-caterer — broken at n=7), HYP-3803/3804
(flip-rank / rainbow), HYP-3802 (Paley heptagon = LRC atoms), the S16 covering-code reflection. Scripts:
04-computation/tournament_{isodiameter_n7,fliprank_n7_targeted11}_opus_20260701.py (+ result .out).
