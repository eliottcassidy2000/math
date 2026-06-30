# The metagraph spectral toolkit, the cyclicity slow mode, and why the LRC obstruction is level 2

*klein-2026-06-29-S3. The owner asked to extend the tournament-metagraph work, define new small useful objects, read arXiv:2507.05905 (Siegel-transform congruence second moments), and find how metagraph ideas help LRC proofs. Builds on THM-584/587/588, opus-S268 (merged-metagraph eigenvalues), oracle-S519 (walk-on-metagraph), and mac-mini's HYP-3543/3551.*

## Two metagraphs, now cleanly separated

"The tournament metagraph" has been two different objects, and conflating their spectra caused a
conjecture where a theorem was available:

- **opus-S268's merged-wiggly metagraph** `G_n/Z_2`: a sparse (unweighted-ish) graph on the `V_merged`
  complement-classes. Messy/semicircular spectrum; `H ≈ 2nd eigenvector` (a degree-`(n-1)` observable);
  Markov gap only *conjectured* `~ c/n`.
- **The arc-flip (dominance-reversal) metagraph** (THM-584/587): the canonical `d`-regular *weighted*
  graph on all `A000568` classes, one edge per single dominance reversal. EXACT integer spectrum `d-2k`
  (THM-584), multiplicities = the signed cycle index (THM-587).

THM-588 fixes the canonical one and replaces the conjecture: its Laplacian spectrum is exactly `{2k}`, and
its **algebraic connectivity is 4 for every `n>=3`** with **walk gap exactly `4/d = 4/C(n,2)`**. The two
metagraphs even have different slow modes — S268's is the `H`-gradient, the arc-flip's is the **cyclicity**
(3-cycle count), and the arc-flip one is the exactly-solvable sibling.

## A small new toolkit (all from `n!` permutations, no `2^{C(n,2)}` enumeration)

Everything below is a coefficient identity in `mult(k)` (THM-587), so it is computable at any `n`:

| object | value / formula | note |
|---|---|---|
| Laplacian spectrum | `{2k : mult(k)>0}` | integer, exact |
| `mult(1)` | `0` (proved) | no first-order invariant |
| `mult(2)` | `1` (= cyclicity) | the single quadratic invariant |
| algebraic connectivity | `4` (all `n>=3`) | Fiedler = 3-cycle count |
| arc-flip walk gap | `4/d` exactly | mixing `O(n^2 log n)` |
| neutral-flip trace `tr(A)` | `2,6,16,60,328,3160,54928` (n=3..9) | total silent-mutation (self-loop) weight; `= sum mult(k)(d-2k)` |
| complexity `kappa(n)` | `(1/N) prod_{k>=1}(2k)^{mult(k)}` | weighted spanning trees; closed form |
| heat trace `Theta(b)` | `sum_k mult(k) e^{-2kb}` | `Theta(1)->1.023` |
| spectral zeta `zeta(s)` | `sum_{k>=1} mult(k)/(2k)^s` | metagraph zeta |

The two cleanest are the **algebraic connectivity = 4** (a universal constant of the family — surprising
for a graph whose size `A000568(n)` explodes) and the **cyclicity Fiedler mode**: the slowest observable
to equilibrate under random dominance reversals is the count of 3-cycles. Local cyclic structure is the
last thing the dynamics forgets.

## The structural punchline: no linear term, one quadratic term

`mult(1)=0` and `mult(2)=1` say something sharp: **the lowest nontrivial content of the metagraph is a
single pairwise (quadratic) mode, and there is no first-order mode at all.** The transposition `(i j)`
negates `chi_{arc ij}`, so no single-arc invariant survives — the antipodal/signed structure kills the
linear level. The first surviving invariant is the quadratic one, the cyclicity. This is not a curiosity;
it is a *template*.

## How this helps LRC (the transfer, HYP-3552)

mac-mini's "one R, three spectra" (HYP-3543) says the LRC cap matrix `M` shares the metagraph's
`R`-even/`R`-odd structure. Read THM-588 through that bridge:

1. **The LRC obstruction is level 2 — a single pairwise second moment.** The metagraph has *no* linear
   invariant and exactly *one* quadratic one. Transferred: the cap's first moment `S_1` vanishes
   (matching "the floor is purely `R`-even"), and the binding content is the single **pairwise**
   `S_2 = -tr(M_odd)/2` — the cyclicity-analogue. The proof effort belongs entirely at the second-moment
   level; there is no linear term to chase and no need for a higher-moment hierarchy.

2. **arXiv:2507.05905 is the analytic engine for that one second moment.** Han–Lee give first/second
   moment formulas for Siegel transforms **with congruence conditions** in dimension 2 — second moments of
   lattice-point counts under `mod q` constraints. The LRC `S_2` is exactly such an object (speed lattice,
   danger region, `mod 14`). So: metagraph says *what* to bound (one mod-14 second moment); the Siegel
   second-moment formula says *how* to bound it; HYP-3551's anti-Littlewood positive floor is the
   positivity it must yield.

3. **The walk picture (S519) gets a rate.** oracle-S519 showed the runner system is a closed `t`-ordered
   walk on the menu metagraph (each arc flipping `2|v_i-v_j|` times). THM-588's exact gap `4/d` and
   "slowest mode = cyclicity" give the *random*-walk decorrelation rate; S519 correctly warns the
   realizable walk is arithmetic, not random — but the metagraph still says the binding *coordinate* is the
   quadratic one, which is where the realizability constraint bites hardest.

The honest scope: the metagraph side is proved (THM-588); the LRC transfers (1)–(3) are a testable program
(HYP-3552), not theorems. But the diagnosis is concrete and falsifiable: **if the LRC cap has no linear
obstruction and a single quadratic one, the metagraph explains why, and a congruence-conditioned Siegel
second moment should close it.**

## The one-line version

The dominance-reversal metagraph forgets everything except 3-cycles last (algebraic connectivity 4,
Fiedler = cyclicity), it has no first-order invariant and exactly one quadratic one, and that
"no-linear / single-quadratic" law — proved here on the clean tournament side — is the structural reason
the LRC obstruction is a single pairwise second moment, the kind arXiv:2507.05905 evaluates with `mod 14`
congruences.

See THM-588, THM-587, [[eigenvalues-of-the-merged-metagraph]] (the other metagraph),
[[the-per-level-signed-cycle-index-borsuk-ulam-ky-fan]] (THM-587), [[the-one-involution-three-spectra]]
(HYP-3543), [[lrc-walk-on-metagraph-proof-attempt-s519]]. Hypotheses: HYP-3552 (LRC transfer),
HYP-3544 (the homology). Reference: arXiv:2507.05905 (Han–Lee).
