# The sixth coordinate deletes the realization leg

**klein-2026-07-09-S209** (following S208's "one node, five names"; owner directive:
attack the mid-band on pair-sum rulers)

## What happened when the mid-band met the modular frame

S208 localized the open content of LRC(14) to "realization for speeds in
(Vmax/14, 9Vmax/14)": runners too fast to ride the snap window, too slow to join the
confined cluster. Five prior coordinate systems named the same residual node. The
pair-sum frame (THM-668) was supposed to be a sixth name for it.

It is — but with a difference that matters more than the naming: **in the modular
coordinate the realization leg does not exist.** A live pair (q, p) — all thirteen
residues v_l·p mod q inside the middle band [⌈q/14⌉, q−⌈q/14⌉] — IS the witness:
t = p/q, done (kps's `LRCPairSumDispatch.mreach_ge_of_pairsum_band` already consumes
it, sorry-free). There is no fast phase to place, no drift to absorb, no snap window
for a mid-band runner to fail to ride. The mid-band dichotomy ("can't ride, can't
cluster") was an artifact of the τ-line frame: mod q, a mid-band speed is just a
residue like any other.

The node itself does not vanish — nothing makes it vanish. It reappears as the
lattice-theta bound on the exact relations n·v = 0, which is the density-floor
decorrelation object the fleet has already closed in measure form. What vanishes is
the OTHER half of the proof: the passage from measure to witness. In every previous
coordinate that passage was a separate open leg (ρ_K→ρ*, ET grid-hit, sqrt
cancellation, (H1) erosion, mid-band safety). Here the witness is the counted object.

## The measurements that forced this reading

- On the k≥8 mid-band-heavy family (the S208 residual), the live-multiplier supply
  is abundant: ~70 live rulers per instance, max LM scaling linearly in V with
  density ≈ (6/7)^13 — the iid product law, visibly.
- The self-hosting guess (witness pair contains the mid-band member) is REFUTED:
  entanglement 63% vs base 65%. The mechanism is not clever pairing; it is bulk
  supply.
- The union-type certificates (C1 ledger, my C4 Hunter upgrade) certify the whole
  census (76/76) but are adversarially breakable to zero while the exact supply
  stays fat. Low-order truncations cannot see a product-density truth — the same
  lesson as every truncated proxy this month (E2, R3, arc counts).
- Universal supra-Vmax liveness holds at V ≳ 91 (2882/2882 moduli, adversarial
  instances included) and FAILS at V ≤ 18 (6148/16947 dead moduli): the law is
  asymptotic, with a finite exact-check floor — precisely the [counting above V0] +
  [native_decide below V0] architecture the project already trusts.

## The aggregated route (HYP-5732)

Sum the supply over all moduli q ∈ (V, 2V]: main term ≈ (6/7)^13·(3/2)V². The
corrections split by mechanism:
- n·v = 0 (exact relations): hit every q — leading-order, the known node, closed in
  measure form by THM-661/663; needs a DISCRETE port (Dirichlet-kernel coefficients
  decay only 1/|j|, so the port must reuse the moment-LP/box machinery, not naive
  absolute convergence — the D3 lesson).
- n·v = mq, m ≠ 0 (ruler-specific): each small (n, m) pins at most ‖n‖₁ moduli —
  total O(V) against a V² main term. New, and cheap.

V²-versus-V slack is the first genuinely widening margin any coordinate has offered:
every τ-line instrument fought at constant-versus-constant (m_P vs 1/14, ∫W vs 1/7,
2.8 vs 1.41). Aggregation over moduli buys a polynomial gap.

## The meta-lesson

We spent a month building instruments to WALK from a measure statement to a witness
(embed, snap, erode, ride). The modular frame shows the walk was the artifact of
insisting the witness lives on the τ-line parametrized by one privileged ruler. Count
witnesses across ALL rulers at once and the walk disappears; what remains is the one
irreducible object (the relation-lattice theta bound) standing alone, with room to
spare. When a proof has two legs and one keeps breaking, sometimes the fix is not a
better leg — it is a coordinate system in which that leg was never needed.
