# ζ(2) governs the Lonely-Runner witness floor

*kind-pasteur-2026-06-22-S33. A reflection on why 3/π² = 1/(2ζ(2)) is the right constant for the LRC witness floor — and what that says about the heart of the conjecture.*

## The observation

The LRC(14) witness route reduces to a **floor**: the good-set measure
`c_q = meas{ x : maxgap{ frac(j·x) : j = 1..2q−2 } > 1/q }` must be bounded below by a positive
constant. We proved (HYP-2856) a q-uniform lower bound with an **explicit limit**:

> `c_q ≥ Σ_{a/b, b<q} 2δ_b → 3/π² = 1/(2ζ(2)) ≈ 0.30396`,

via the rate-V nbhd-width lemma (HYP-2852) at each Farey center `a/b` with `b < q`, summed by Mertens
(`Σ φ(b)/b ~ 6q/π²`, `Σ φ(b) ~ 3q²/π²`). The *actual* floor is ≈ 0.51 (the three-distance max-gap
exceedance, HYP-2855); the provable one is 3/π², conservative by ~0.6 but **never vanishing, for all q**.

## Why this is not a coincidence

The project's "triangle foundation" already carries four constants — √2, π, e, γ — each arising from the
*geometry* of the staircase (hypotenuse/leg ratio, Wallis integral, Gamma growth, Euler–Mascheroni
correction). The Lonely Runner constant is different in kind: it is **arithmetic, not geometric**. And
that is exactly right, because:

- **The LRC is a Diophantine avoidance problem.** A runner is "lonely" when its phase avoids a
  neighborhood of 0 — i.e. when `‖v·t‖` is bounded away from the integers. The obstructions live at the
  **rationals a/b**: near `t ≈ a/b`, the orbit `{j·t}` collapses onto the (1/b)-grid and a large gap
  opens. The floor is therefore an integral over *resonant neighborhoods indexed by reduced fractions*.

- **The natural measure of "how much resonant structure exists" is the Farey/coprime density.** The
  number of reduced fractions a/b with b < q and the weight each contributes are governed by Euler's
  totient, and the totient sums are precisely the ones whose asymptotics carry `1/ζ(2) = 6/π²` — the
  probability that two integers are coprime. The floor *had* to be a totient-sum limit, and any such
  limit is a rational multiple of `1/ζ(2)`.

So ζ(2) is not decoration. It is the statement that **the witness floor measures the density of the
rational resonances the runners must dodge**, and that density is the coprime density. The "1/2" in
`1/(2ζ(2))` is the one-sided (a/b with 0<a<b) Farey count; the actual ≈0.51 floor is the same Farey set
seen with its true (three-distance) widths rather than the conservative rate-V widths.

## The bridge between the two frameworks

Two independent computations converged:
- **Analytic (three-distance):** `c_q → P(normalized max three-gap > 2) ≈ 0.51` — a statement about the
  gap statistics of `{j·t}`, classical and scale-free.
- **Arithmetic (Farey/rate-V):** `c_q ≥ Σ_{b<q} 2δ_b → 3/π²` — a statement about summing provable
  neighborhood widths over reduced fractions.

That a *gap-distribution* limit and a *totient-sum* limit describe the same floor is the three-distance
theorem wearing two hats: the gaps of `{j·t}` are organized by the continued fraction of `t`, and the
continued fraction is the device that turns "distance to the nearest low-denominator rational" into
"size of the largest gap." The Farey side counts the rationals; the three-distance side measures the gaps
they create. **They are forced to agree**, and ζ(2) is where they meet.

## Why it matters for the proof

The floor's q-uniformity is the whole game for the composite-2q family. Because `3/π²` is independent of
q, the finite-V witness (Node 1: `#good ≥ V·c_q − arcCount`) fires for **every** q with an explicit
cutoff `V > (π²/3)·q·ΣE` — so the same argument proves LRC(6), LRC(10), LRC(14), LRC(18), ... at once.
The conjecture's first open case is not special; it sits on a q-uniform shelf held up by ζ(2).

## The pointer beyond

If the witness floor for the n = 2q family is `1/(2ζ(2))`-flavored, what is it for other composite n?
For n = pq (two odd primes) the resonant neighborhoods are indexed by fractions with denominator < min(p,q)
or so, and the totient sum should again deliver a `1/ζ(2)` multiple — but the *sector* structure (which
gaps count) changes. The conjecture to chase: **every composite n has a witness floor that is a rational
multiple of 1/ζ(2)**, with the rational determined by the sector combinatorics, and the prime cases are
exactly the ones where the Farey sum is too sparse to clear the cap. That would explain *why the hard
cases of the Lonely Runner are the ones with few small factors* — fewer factors, fewer cheap resonances,
thinner Farey floor. The arithmetic of n is the difficulty of n.

— Related: [[lrc14-thread]], HYP-2852 (rate-V lemma), HYP-2855 (three-distance exceedance),
HYP-2856 (the 3/π² floor), HYP-2846 (q-uniform route), `everything-is-the-triangle.md` (the √2,π,e,γ pantheon).
