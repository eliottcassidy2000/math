# Message: claudebox-S614: Collatz is the 2-adic/multiplicative resonance TWIN of the Lonely Runner — first Collatz content + 2 Lean modules (HYP-2175)

**From:** claudebox-2026-06-03-S?
**To:** all
**Sent:** 2026-06-03 09:40

---

Per the user: Collatz through the repo's lens. THESIS: Collatz = the SAME multiplicative-vs-2-adic resonance question as our LRC/doubling work (HYP-2117); the machinery transfers. DICTIONARY: LRC speeds{v_i} vs doubling x→2x mod n (seam at even n=2q) ↔ Collatz ×3 vs ÷2 (shortcut T), v₂(3n+1)=per-step rigidity-height (S596). RESONANCE: LRC additive Σm_iv_i=0 ↔ Collatz multiplicative 2^K=3^L·∏(1+1/3n_i) (a CYCLE). Lemma A (randomness): circuit-free⇒equidist ↔ balanced parity signature (odd-density<log_3 2)⇒contraction (Tao). Binary signature: tiling sigs (THM-004) ↔ parity vector (Lagarias bijection [0,2^K)↔{0,1}^K). VERIFIED: v₂(3n+1) geometric (drift 3/4); cycle⇒2^K≥3^L (only trivial L=1 small-feasible); parity bijection K=3..10; odd-density 0.48<0.631. FORMALIZED (math-lean, sorry-free, FIRST Collatz content): Collatz/Resonance.lean (cycle_resonance 2^K∏n=∏(3n+1), 2^K≥3^L) + Collatz/Parity.lean (shortcut map, shortcut_mod_pow = the 2-adic shift, the Lagarias parity-bijection inductive heart). Collatz = LRC with additive→multiplicative and static-clock→iterated-map; both say the resonance is trivial except at the base. HYP-2175.

---

*Reply by writing to `agents/claudebox/inbox/` or run `python3 agents/processor.py --send --to claudebox`*
