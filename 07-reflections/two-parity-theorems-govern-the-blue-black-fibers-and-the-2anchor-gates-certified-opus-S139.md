---
source: opus-2026-07-07-S139
status: THM-644 (lemma PROVED; fiber law verified 72/72; ANTI-REDEI conjecture opened, verified n<=7;
  parity law proved) + the four 2-anchor handoffs delivered (finite dip exact, certified limit table,
  rigorous double-cover chain sufficing for k=13 on the limit class, observer-rank prefilter)
tags:
  - tournaments
  - metagraph
  - blue-black
  - redei
  - anti-redei
  - lonely-runner
  - two-anchor
---

# Two parity theorems govern the blue/black fibers — and the 2-anchor gates, certified

**opus-2026-07-07-S139.** Owner two-front directive. Fleet context: mac-mini-S46 (THM-643
reserve) and klein-S161 took the line-side parity/node-type formulas; my differentiated lane
was the **fiber allocation** — which particular tilings map to which class — plus the four
LRC 2-anchor handoffs. Both fronts landed; the metagraph front landed something unexpected.

## The metagraph front: the picture is TWO Rédei-type parities

The census (n = 4..7) + targeted verification produced a chain that assembles exactly the
"grand restricting picture" the owner asked for (all in THM-644):

1. **Grid symmetry is the reversal anti-automorphism** (2-line proof): `t` gridsym ⟺
   `ρ(i) = n+1−i ∈ AntiAut(T_t)`. The mysterious explorer predicate is a first-class
   tournament symmetry; "non-transpose-self ⟹ pure black" becomes a theorem.
2. **The fiber law** `g(C)·|Aut| = H_anti(C)` (anti-symmetric Hamiltonian paths) — verified
   on every class n ≤ 6; it is LEM-003 restricted to the ρ-equivariant locus. So the
   allocation is orbit-stabilizer twice: `N = H/|Aut|` counts all tilings; `g = H_anti/|Aut|`
   counts the blue layer.
3. **The parity law**: `N(C)` is ALWAYS odd — proof in one line: `H` odd (**Rédei**, the
   project's founding theorem) and `|Aut|` odd (tournament automorphism groups have odd
   order). And `g(C)` is odd-or-zero, zero exactly off the transpose-self classes — granted:
4. **THE ANTI-RÉDEI CONJECTURE** (new): `H_anti(T)` is odd for every `T` admitting an
   anti-automorphism. Verified on all self-converse classes n ≤ 6 and via the all-odd
   n = 7 g-spectrum `{1,3,5,7,9}`. This is the σ-equivariant refinement of Rédei — the exact
   completion THM-587's frame ("Rédei is σ-odd, lonely is σ-even") predicted. It implies
   "every self-converse class contains an anti-symmetric Hamiltonian path" (no
   transpose-self class is pure black) — existence as a parity corollary, pure Rédei style.
5. Line closed-forms (`blue = 2^{(m+f)/2−1}`, `black = 2^{m−1} − blue`) verified; node-type
   counts tabulated (pure-blue/mixed/pure-black = 4/84/368 at n = 7); blue class-self lines
   vanish at odd n (new alternation conjecture); the flip-partner class map is non-constant
   on most classes — the line-multigraph's distribution is the remaining descriptive layer.

The owner's "blues contribute odd amounts and blacks even" is now precise: per class, the
blue count is odd wherever nonzero (Anti-Rédei), and the black count is even exactly there
(odd N minus odd g). **Two parities — one classical, one new — restrict the whole
structure.** Natural next: prove Anti-Rédei by a ρ-twisted parity involution (the repo's
home move); then the fiber side of the metagraph is complete modulo the flip-multigraph
distribution.

## The LRC front: the four handoffs

- **Near-limit sweep (exact):** the finite dip below the T² limit is real and resonant —
  k=8 min PA₂ = **0.761046 at (a,d) = (3,4)** (not the limit's 0.8017; perturbations don't
  deepen); k=9: 0.6864 at (5,8); k=10: 0.5697 at (1,2). Margins over T_k stay ≥ +0.14. The
  rigidity lemma's honest constants are these finite-resonance values.
- **Certified 2-anchor limit table:** grid + Lipschitz error bound (each config point moves
  at speed ≤ K): limits 0.8013/0.6944/0.5976/0.4978/0.4292/0.3665 with certified margins
  **≥ +0.1825 … +0.3092** at every k. First rigorous-error version of the table.
- **The double-cover chain (rigorous, and it bites):** `gap∋0(E,2x) ≤ 2·min(gap∋0,gap∋½)(E,x)`
  (0/20000 pointwise violations; the ½-preimage argument is a proof), so
  `PA₂ ≥ P(min > 1/7) ≥ P(gap∋0 > 2/7)`. Exact 2/7-origin-mass on the limit class:
  2321/10290, …, **8087/64680 = 0.1250 at k=13 — which CLEARS T₁₃ = 0.0565**: the chain
  alone rigorously handles the k=13 leg on the binding class, reducing it to a
  **single-anchor tail at the fatter threshold 2/7** — a strictly simpler object than the
  two-anchor 1/7 statistic. (Lossy at k ≤ 12; the 2-anchor route stays live there.)
- **Observer-rank prefilter:** rank ≤ 1 predicts ledger-zone membership at ~79%, rank ≥ 2
  predicts spread at ~70% — a cheap prefilter, not a sharp classifier; reported honestly.

Also closed this session: the S137 a=2 pinch corner completed exactly — **corner minimum
12680829120251011/23152045800495150 ≈ 0.5477 = 9.7× m_P** at the far-element block (1¹¹,66);
all arrangements at the worst cell tie. kps-S62's corner is settled.

## The convergence note

The two fronts touched: the Anti-Rédei conjecture is the tournament-side σ-parity; the
½-anchor/double-cover is the runner-side ℤ/2 structure; THM-587 already said these are the
two isotypic halves of one equivariant object. Today each half acquired its concrete
theorem-shaped core: H_anti-oddness on one side, the 2/7-single-anchor reduction on the
other. Both are parity statements about a reversal.
