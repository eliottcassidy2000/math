# Characterization compounds into proof: the value and the obstruction are one object

*klein-2026-07-01-S68. A synthesis reflection (HYP-3800) on how S65–S67 + the concurrent mac-mini /
kind-pasteur threads assemble, and what turns accumulated understanding into a proof.*

For a dozen sessions the covering-min problem looked like an existence search: is there some covering set
— possibly using enormous, cleverly tuned speeds — that beats `n/Phi6`? Existence searches over unbounded,
disordered spaces are where proofs go to die. This session the pieces finally arranged themselves so that
the search is gone, replaced by one analytic inequality on a finite object. The move that did it was not a
new computation but the *characterization* of properties we already half-knew — and I want to record how
characterization, done patiently, compounds into a proof scaffold.

**Three characterizations, three reductions.** Each session added an adjective to the far element, and each
adjective removed a dimension of the search.

- S66/S67 established the far element is *narrow*: a huge speed's danger is a needle, coupling to the
  lonely set only as `O(1/w)`. That killed the single-far case — a lone huge speed cannot cover `L_C`.
- S67 established the multi-far interaction is *self-similar*: pairwise correlations peak exactly when the
  comb difference lands in `(n-1)Z`, the same lattice at every scale. That collapsed the tower of orders
  `r=2,3,4,…` into one recurring cell.
- S68 (today) established that cell is *translation-invariant*: `corr2(a, a+δ)` depends only on the
  difference `δ`, not on where `a` sits — verified constant (`≈0.10`) from `a=300` to `a=2000`. That
  collapsed the *scale* itself. The infinite axis of "how huge is the speed" simply drops out; what
  remains is the difference-set and its phases in `Z/Phi6`, a finite arithmetic object.

Narrow, self-similar, translation-invariant. Each is a property, not a theorem; but stacked, they turn
"search all speed-configurations" into "evaluate a phase-histogram on `Z/Phi6`." That is what
characterization buys — not the proof, but the *shape* of the proof, and a domain small enough to prove
things about.

**Fixing the sign is worth more than bounding the size.** The subtlest gift was S67's observation, sharpened
today: the strong correlations are *redundant* (positive) — resonant far combs lock together and
double-cover the same part of `L_C`. Redundancy is the favorable sign: it makes survival *larger*, not
smaller. So a proof of "a lonely time always survives" does not have to defeat the resonances at all; the
resonances are on its side. It has only to bound the *spreading* (negative) correlations — the combs that
would cover different parts of `L_C` — and those are exactly the weak, off-lattice, equidistributed ones.
Knowing the *sign* of the enemy told us the enemy is an ally, and left a much smaller adversary. This is a
general lesson: before bounding a correction, learn its sign; a signed quantity you understand is worth
more than an absolute quantity you can merely majorize (the absolute sum here diverges — MISTAKE-078 — and
every session that tried to bound `Σ|·|` hit that wall; the sign is what walks around it).

**The value and the obstruction are the same object.** The covering-min *value* is `n/Phi6` — the
hexagonal point, the Eisenstein denominator, born from the killer `n(n-1) ≡ -1 mod Phi6`. This session the
*obstruction* — the invariant that decides whether a far speed couples to the lonely set — turned out to be
`p(w) = nw mod Phi6`, a residue in the very same `Z/Phi6`, with the very same killer identity making `n-1`
its unit step, and the very same continued fraction `[0; n-1, n]` (true for all `n`) organizing its
scales. The answer to the problem and the reason for the answer are literally the same arithmetic. This is
the recurring signature of this project — everything is the triangle, the hexagon, `Phi6` — but here it
has teeth: because the obstruction lives in `Z/Phi6`, it factors by CRT over the Eisenstein primes of
`Phi6(n)`, and the prime instruction from last session lands exactly: `2` is the antipode `t ↔ 1-t` (the
sign), `n/2` gates (the `7`-vanishing), and the primes of `Phi6` carry the phases. The proof, when it
comes, will be a statement about covering `Z/Phi6` modulo `±1`, and it will factor over those primes.

**What is left, stated honestly.** The scaffold has three steps: bounded speeds (done, ILP), single-and-few
far (done, Fourier + total-variation), and multi-far (open). The open step is now a single inequality —
that the spreading correlations sum to less than the independence term `(1-2r)^r` — over a finite
phase-histogram, with three tools pointed straight at it (kind-pasteur's signed Cauchy–Schwarz and moment
relaxation, mac-mini's per-speed decay, and the additive-energy control that says redundancy is harmless).
I have not closed it. But the difference between "search the huge-speed jungle" and "bound one signed sum
on `Z/Phi6` mod `±1`, prime by prime" is the difference between a problem and a proof-in-progress.
Characterization did that. The job the owner set — new definitions, key angles, structure — is exactly the
job that converts the one into the other: you cannot bound what you have not named, and once the far
element had enough adjectives, the thing to bound named itself.
