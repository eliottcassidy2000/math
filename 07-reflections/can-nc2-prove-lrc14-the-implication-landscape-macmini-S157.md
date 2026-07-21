# Can NC2 prove LRC(14)? The implication landscape, honestly mapped

*mac-mini-2026-07-20-S157. The owner: "if we assume the two-dimensional nullcone statement
conjectured by Derksen–van den Essen–Zhao is true, what can we prove? Can we make LRC(14) be
proven by that conjecture ⟹ GMC(2) ⟹ LRC(14)? Or another path — get creative." This maps the
implication landscape rigorously. The short answer is a precise NO with one real arrow, and the
creative content is what the shared structure does and does not buy.*

---

## The objects, precisely

- **NC2** (klein THM-1510 C): *in two real Gaussians, the only `P` with all moments `E[P^m]=0`
  are the purely holomorphic or antiholomorphic ones* — i.e. the moment nullcone is exactly the
  charge-one-sided locus. This is the "2-D nullcone conjecture." It equals my atom-covering
  statement `V(all atom forms) = one-sided` (THM-1780).
- **GMC(2)** (Derksen–van den Essen–Zhao): `ker(E)` is a Mathieu–Zhao subspace of `ℂ[Z,W]`.
- **LRC(14)**: for any 13 distinct positive integer speeds, `max_t min_j ‖v_j t‖ ≥ 1/14`.

## The one real arrow: NC2 ⟹ GMC(2)

This is **proved** (klein THM-1510 C, three lines by weight counting). So *assuming NC2 gives
GMC(2)* — the last open Gaussian moment case — for free. That much of the owner's chain is real.

## GMC(2) ⟹ LRC(14): no — a fourfold obstruction

Both problems live on the **same relation lattice** `{k ∈ ℤⁿ : k·v = 0}` = the charge lattice
(boxeph THM-1820; my S156, where the minimal relations are pairs). But the implication fails,
and not for lack of work:

1. **Different measure / different weights.** GMC(2) is the Gaussian on `ℂ` — the relation-
   lattice weights are **factorials** (`E[V^k]=k!`, Laplace). LRC is Lebesgue on the torus —
   the weights are **sincs** (`sin(2πkδ)/(πk)`, from the indicator's Fourier series). Verified:
   the factorial tower is **monotone** — it never cancels, so the moment nullcone is clean
   (one-sided only). The sinc tower **oscillates and changes sign** — so the good-set sum *can*
   vanish (a covering), which is the entire difficulty of LRC. A nullcone theorem for one weight
   system says nothing about the other.
2. **Powers vs products.** GMC constrains `E[P^m]` — one polynomial, `m` identical copies. LRC
   constrains `∫₀¹ ∏_j f_j(v_j t) dt` — `n` *different* functions, one copy each. The LRC
   product `∏ f_j` is a single multinomial term of `(Σ f_j)ⁿ`; GMC/NC2 constrains the *sum* of
   all such terms, never that one term alone.
3. **Opposite classification.** The GMC nullcone is *nonempty* (one-sided members exist); LRC
   asserts the covering nullcone below `1/(n+1)` is *empty*. Same lattice shape, opposite
   answer (boxeph THM-1820).
4. **LRC(14) is number-theoretic.** It is about specific integer speed sets, and the repo's
   program (density floor THM-663 closed; finite `V_max` glue THM-527-A; Lean) routes through
   covering combinatorics, not Gaussian analysis.

> **Verdict.** `NC2 ⟹ GMC(2)` is real; `GMC(2) ⟹ LRC(14)` is **not** — a genuine obstruction,
> not missing work. The two are *parallel* relation-lattice pairing problems; the link is a
> shared lattice (technique transfer), not a logical implication.

## What assuming NC2 actually buys

Honestly, not much beyond GMC(2) itself:

- **JC is dead as a downstream.** `GMC(n) ⟹ JC(n)` needs GMC for *all* `n`; but `GMC(n≥3)` is
  now **false** (the Gaussian counterexamples, THM-1500), and JC(3) is itself false (THM-1300).
  So GMC(2) gives nothing toward the Jacobian Conjecture — that cascade is severed.
- **TNC is already a theorem** (= DvdK, THM-1630), independent of NC2.
- So NC2's payoff is **exactly** the 2-variable Gaussian Mathieu–Zhao subspace (GMC(2)). It is
  a beautiful capstone of the Gaussian moment story, but it is a *leaf*, not a hub — nothing of
  weight hangs below it.

## The creative content: what IS shared, and one honest speculative path

What genuinely transfers between NC2/GMC(2) and LRC is **machinery, not theorems**:

- The **pair-reduction** — the relation lattice is pair-generated on both sides (S156).
- The **transitivity discriminant** — the pair-in-radical closure is the Vandermonde =
  signed tournament sum (THM-1815/1805), and LRC tightness is "relation-richest = the deepest
  nullcone point" (boxeph). Both extremals are read through the same discriminant.
- The **first-return / renewal** structure (THM-1770) is literally the LRC covering/three-gap
  toolkit seen on the moment side.

So the productive direction is *not* "assume NC2, deduce LRC" but "the discriminant/first-return
machinery that closes the moment nullcone is the same machinery LRC's covering needs" — a
tool-sharing, which is where boxeph, klein, opus and I have converged.

**The one honest speculative path** (stated so it is not mistaken for a proof): a *weight-
agnostic* relation-lattice conjecture — "for every admissible weight system, the pairing has the
transitive/one-sided locus as its extremum" — would specialize to NC2 (factorial weights). But
it would **still not give LRC**, because LRC is a *positivity* (the good set is nonempty just
below threshold), not a nullcone, and the sinc weights' sign changes are exactly what such a
conjecture cannot control. So even the strongest unified statement stops at the measure barrier.

## The honest bottom line for the owner

- **Yes**, assuming NC2 proves GMC(2) (real, klein THM-1510).
- **No**, GMC(2) does not prove LRC(14), and neither does NC2 by any known path — the
  obstruction (factorial vs sinc, powers vs products, opposite classification) is genuine.
- **LRC(14) does not need NC2**: it is nearly closed by its own program (density floor done;
  finite glue + Lean remaining), and that program is combinatorial/number-theoretic.
- **The real payoff of the connection** is technique transfer via the shared relation lattice,
  pair-reduction, and transitivity discriminant — which the fleet is already exploiting.

## Honest scope

- The obstruction is argued and verified (factorial monotone vs sinc oscillating), not a formal
  impossibility theorem; "no known path" is the honest status, and I would not claim more.
- `NC2 ⟹ GMC(2)` and `NC2 = the atom-covering` are cited (THM-1510, THM-1780), not re-derived.
- The JC-cascade-is-dead point depends on `GMC(n≥3)` false (THM-1500) and JC(3) false
  (THM-1300), both in canon.
- I make **no** claim on LRC(14) itself here — only on what does and does not imply it.

---

*Cross-links: THM-1510 (NC2 ⟹ GMC(2)), THM-1780 (NC2 = atom-covering), THM-1815/1805 (the
transitivity discriminant), THM-1820 + my S156 reflection (LRC ⟷ moment dictionary, pair-
generated lattice), THM-1500/1300 (GMC(n≥3) and JC(3) false), THM-663/527-A (the LRC(14)
program). Artifacts: `04-computation/nc2_gmc2_lrc_implication_macmini_S157.py` (+out).*
