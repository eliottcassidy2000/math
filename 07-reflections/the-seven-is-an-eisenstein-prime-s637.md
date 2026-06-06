# The 7 is an Eisenstein prime — and the density face is already exhausted at 5

*S637 reflection. On the chromatic number of the plane: what is actually provable, and where the
difficulty genuinely lives.*

The user asked me to spend a long session on `χ(ℝ²) ∈ {5, 6, 7}` and to attempt any novel statement,
even something small like eliminating one of the three. I want to be honest first: I did not eliminate
one. Nobody has. Doing so would be a famous result, and a long session of clever reframing does not
manufacture a 6-chromatic unit-distance graph or a 6-coloring out of nothing. So I went looking
instead for the two things I *could* nail down precisely — what the numbers **are**, and what an
elimination would actually **cost** — and both turned out to land squarely on the arc.

**The 7 was never a free parameter.** The hexagonal upper bound is not "seven happens to work." It is
reduction modulo the Eisenstein prime of norm 7. `7 = N(3+ω)`, `ℤ[ω]/(3+ω) ≅ 𝔽₇`, and the seven
hexagon colors are the seven residues. The six sides of the hexagon — the six nearest-neighbour
directions, the six Eisenstein units, `6 = 2·3` — map exactly onto `𝔽₇*`, and the center plus its six
neighbours fill out all of `𝔽₇`. A hexagon and its ring is an `𝔽₇`-torsor. *That* is why seven is
both necessary and sufficient for this tiling. And it is the **same 7** the arc has been circling for
twenty sessions: `7 = Φ₃(2)` (the forbidden tournament value), `7 = M₃` (the Mersenne), `7 = N(3+ω)`
(the chromatic bound) — three readings of one prime, because `Φ₃` is the minimal polynomial of `ω` and
a prime splits in `ℤ[ω]` exactly when `Φ₃` has a root mod it. The witnessing root mod 7 is `2`. The
chromatic upper bound and the tournament's forbidden gap are the same arithmetic.

Then the small shock that made the session worth it: **the evaluation point is a cube root of unity
mod 7.** The tournament invariant is `H = I(Ω, 2)` — the independence polynomial at `x = 2`. The
forbidden value is `Φ₃(2) = 7`. Reduce `Φ₃` at 2 *modulo its own value*: `2² + 2 + 1 = 7 ≡ 0`. So `2`
is a primitive cube root of unity in `𝔽₇`, and so is `ω = 4 = 2²`. The two cube roots of unity mod 7
are exactly `{2, 4}`, inverse to each other. The number you evaluate the partition function at, the
cube root that builds the hexagon, and the root of the cyclotomic polynomial are *the same two
elements* of the one field where they all live. I did not expect the evaluation point to be the
resonance itself. It is.

**The other half is a no-go, and that is also progress.** I kept asking: could the density method — the
LRC `p₀`, the 1-avoiding density `m₁`, the face of the unified object I have leaned on all arc —
eliminate one of {5,6,7}? The answer is a clean, provable *no*, and stating it sharply is itself the
result. A measurable `k`-coloring forces some class to density `≥ 1/k`, each class is 1-avoiding so
density `≤ m₁`, hence `χ_m ≥ 1/m₁`. With `m₁ ≤ 0.247` this proves `χ ≥ 5`. But a 1-avoiding set of
density `0.2293` *exists* (Croft), so `m₁ ≥ 0.2293`, so `1/m₁ ≤ 4.36 < 6` — *forever*. To reach 6 you'd
need a 1-avoiding set to be impossible above density `1/6`, and one exists at `0.229`. The density bound
is capped at 5 by a construction, not by our cleverness. So the face of the object I trust most is
**provably exhausted** at exactly the bound it already gives.

That is the real lesson, and it is the same lesson as LRC(14). There, the first-moment / union bound is
vacuous and all the content lives in the correlations between arcs (HYP-2195). Here, the single-class
density bound is vacuous past 5 and all the content lives in the correlations between color classes —
how they must interlock — which is precisely the off-diagonal, multi-class Krawtchouk/Delsarte dual I
flagged as the open lever for LRC(14). Plane and circle, same wall, same door. The difficulty in
Hadwiger–Nelson is not the number theory of the construction (that is settled, it is an Eisenstein
prime) and not the density of one class (settled, capped at 5). It is the pair-correlation of classes —
the same off-diagonal certificate the whole arc keeps reducing to.

So: no elimination. But the 7 is named, the evaluation point is the resonance, and the density face is
closed with a proof. The thing I can honestly hand forward is a sharper map of where the hard part is —
and it is exactly where the arc already said it would be.
