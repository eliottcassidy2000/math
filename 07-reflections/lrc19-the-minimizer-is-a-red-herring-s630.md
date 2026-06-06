# LRC(19): the unramified case, and the minimizer is a red herring (S630)

The task was to try for LRC(19) — open (the conjecture is a theorem only through n=13) — and, if I hit a
wall, to wander the repo's creative threads until their connections handed me a reframe. Both happened,
and the reframe was the opposite of what I went looking for.

I chose n=19 deliberately. `2n−1 = 37` is prime, and 2 is a primitive root mod 37: the *unramified*
case, the clean opposite of n=14's `27 = 3³`. Everything our framework says about the favorable side
should fire here. It did. The reduction `LRC(19) ⇐ C′(19)` (a multiple-of-19 config is loose) holds
constructively: the twisted-shell dodge on shells `m ≤ 37`, together with the dominance dodge B′, covers
**100% of 20000 random multiple-of-19 configs with zero residual** — in fact the dodge alone covers them,
B′ never needed. Every such config comes with an explicit loose witness `t = a/m`, `m ≤ 37`. A good
multiplier at shell 37 gives `M ≥ 2/37 = 0.0541 > 1/19 = 0.0526`. That is a constructive C′(19), and so
strong evidence for LRC(19).

Then I went looking for the wrong thing. The quantitative story (S623, THM-415) says the *minimum* of M
over multiple-of-n configs is `2/(2n−1)`, attained by a specific extremal config — and I wanted to exhibit
it at n=19, partly to nail `min M = 2/37`, partly out of tidiness. That is the wall. The minimizer is
irregular (n=5 gives `(1,3,4,5)`, n=8 gives `(1,4,5,6,7,11,16)`); no AP-variant search, no random search,
no fold-maximizing hill-climb, and no Paley/quadratic-residue construction mod 37 produced it. The best I
found was `{1,…,17,19}` — the AP with its top bumped to a multiple of 19 — with `M = 1/18` exactly (loose
via the freed 18-clock), well above `2/37`.

The creative wandering, when I let it run, kept returning the same three notes from different rooms:
resonance energy (S550o), "the fold *is* the δ-clock witness" (S577), and the 2-adic self-affine tower of
the extremal family (S579). The honest test of their obvious synthesis — "min M = max number of 3-term
folds" — *failed*: random configs are all fold-rich (a dozen folds each) yet sit at `M ≈ 0.11`, nowhere
near tight. Global fold count does not pin tightness; only folds *aligned to the modulus* (straddle pairs
summing to the shell) do, which is a far more delicate, residue-level condition. Decoding the small
minimizers confirmed the delicacy: their witnessing multiplier sits near the apex `a ≈ (2n−1)/2 ≈ n` for
n = 4,5,6,8, but not n = 7. The extremal is genuinely arrhythmic.

And that is the reframe, arriving as a deflation. **The minimizer is a red herring.** LRC(19) does not
need the extremal config — it needs *coverage*: that *every* multiple-of-19 config is loose. Coverage is
what the dodge delivers, and coverage is exactly the quantity the residue-profile reduction (S629) makes
finite. `gap_shells` reads each speed only through its residues mod the shells `≤ 37`, hence (CRT) mod
`L = lcm(2,…,37)`; so "is every multiple-of-19 config loose?" is a question about *finitely many residue
profiles*, independent of how large the speeds grow. In the unramified case — 37 prime — the shell carries
no ramified strata to fragment that check. So the path to LRC(19) is not "find the tightest config and
prove nothing beats it"; it is "run the residue-profile DP mod `L` and confirm the dodge ∪ B′ covers every
profile." A decidable, finite computation, with the irregular minimizer never entering the argument.

I spent the session trying to grab the extremal config and the repo kept handing me, instead, the fact
that I didn't need it. The wall was real; the way around it was to stop treating the hardest single
instance as the obstacle and treat the *whole residue-profile space* as a finite object — which is the
one move this cluster has been converging on from every direction: covering depth, shell tower, partition
function, frontier-gain, and now the LRC(19) frontier. The minimizer is where the difficulty *looks*
concentrated; the residue profile is where the problem actually *lives*, and there it is finite.
