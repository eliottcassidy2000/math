---
id: THM-791
renumber_note: "Renumbered from THM-790 by kind-pasteur cont.6: opus-S304/305's leg-law file claimed THM-790 first (rewrite 19:03 on 07-14 vs my push ~23:xx; different filenames masked the collision at claim time). Content unchanged."
title: The H-COMPANION to the transitivity flow — the same blue/black complement-line flow measured on the project's Rédei invariant H (Hamiltonian-path count) instead of the degree-derived C3/E4 axes of THM-785/787. (a) VERIFIED LAW (even n): blue lines preserve H mod 4 (ΔH ≡ 0 mod 4; n=4,6 exhaustive), while at odd n ΔH ≡ 0 mod 2 only (Rédei-trivial; ΔH=2 occurs at n=5 — the mod-4 law is genuinely a parity-of-n phenomenon, the H-analogue of THM-785's ΔE4 ≡ 0/8 mod 16 dichotomy); (b) the blue |ΔH| distribution is near-flat vs the decaying black one; (c) score majorization-comparability: s(t) and s(t̄) are majorization-COMPARABLE for every tiling (zero incomparable pairs, exhaustive n ≤ 6) — finer than the 1-D C3 flux, following the exact affine law s(t̄) = C − s(t), C = (n−2, n−1, …, n−1, n) (= THM-785 §2's degree transformation); (d) the pure-blue census formula pure_blue(n) = ⌊(n+1)/2⌋ − [n even] (kps-S66 conjecture) is CONFIRMED through n=7 (2,1,3,2,4; n=7 count 4 = codex THM-785 §6's independent atlas)
status: (a) PROVED for all even n (cont.6, referee-verified every link on all 1096 tilings n=4,5,6 — see The Proof below), in the STRONGER STRATUM FORM — at even n every iso class containing a grid-symmetric tiling has H ≡ 1 (mod 4); the blue ΔH ≡ 0 (mod 4) law and the certificate "H ≡ 3 (mod 4) ⟹ pure black (even n)" are corollaries; n ≥ 8 needs no computation (the proof is general). (b) VERIFIED n=4..7. (c) CORRECTED at n=7: comparability is a BLUE law — blue lines have ZERO incomparable pairs through n=7 (176 toward-regular / 80 toward-transitive at n=7), while BLACK incomparability first appears at n=7 (268 of 16,128 lines ≈ 1.7%); blue-comparability proof target: grid-sym forces d_i + d_{n−1−i} = n−1 (anti-palindromic scores) + s̄ = C − s. (d) CONFIRMED n≤7 and REFINED at n=7: the four pure-blue classes are (H=1, transitive), (H=3, (0,1,3,3,3,5,6)), (H=9, (1,1,1,3,5,5,5)), (H=15, (0,3,3,3,3,3,6)) — the REGULAR class is NOT pure-blue at n=7 (n=5's regular-is-pure-blue was accidental; the census formula ⌊(n+1)/2⌋−[n even] still exactly matches 2,1,3,2,4)
source: kind-pasteur-2026-07-14-S128 (cont.5; owner directive: trace the flow of transitivity, quantify symmetry/imbalance, build the ordering)
depends_on: []
related:
  - THM-785 (codex-S9: C3 flux, closed blue binomial law, categorical interface (13)/(14), oriented quotient imbalance, the flow address (15)) — the frame this companions; my independent n≤6 run reproduced its interface law and node censuses (PB,M,PK) = (3,5,2), (2,10,22) exactly before reading it
  - THM-787 (opus-S304: E4 axis, blue-avoids-pure-black proved, parity/max conjectures proved by 785)
  - HYP-6855 (this session's log entry); kps-S66/HYP-4997 (blue/black line counts; the pure-blue census conjecture now confirmed)
  - MINIMAL-INVARIANT ORDERING remark: the 4-stage order (score-seq via (Σs² desc, lex), phase PB<MX<PBk, H, canonical word) resolves 3/10/28 of 3/10/34 classes at n=4/5/6 before the single arbitrary stage (canonical word breaks 0/1/6 residual ties) — a low-coordinate alternative to THM-785's 8-coordinate flow address when minimal arbitrariness is the goal; H does real separating work in it (stage gains +0/+2/+11)
---

# THM-791 — the H-companion laws

THM-785/787 put the transitivity flow on degree-derived axes (C3, E4). H — the Hamiltonian-path
count, the invariant this project is built on (Rédei; tilings·|Aut| = H) — is not determined by
the degree sequence, so the flow's H-behaviour is independent structure. Measured exhaustively
(all complement-tiling lines, strict explorer conventions):

**(a) The blue mod-4 law (even n).** Along every blue line at n=4, 6: H(cls(t̄)) ≡ H(cls(t)) mod 4.
At odd n only the Rédei-trivial mod 2 holds (n=5 has blue ΔH = 2, 10). The even/odd-n dichotomy
mirrors THM-785's ΔE4 ≡ 0 vs 8 (mod 16) — but for H, where no degree-sequence argument applies.
n=4: |ΔH| ∈ {0,4}; n=6: {0,4,8,12,16,20,24,28} (all ≡ 0 mod 4, distribution 5,6,5,6,6,1,2,1).

**(b) Blue is flat, black decays.** The blue |ΔH| distribution is near-uniform across its support
(n=6 above); black's decays (n=6: 48,38,48,66,50,42,36,44,12,30,24,18,10,10,4 at |ΔH|=0..28 even).
The H-axis sees the same "blue = long symmetric strides, black = short diffusive steps" picture
as the C3 axis, on an invariant the C3 laws cannot reach.

**(c) Majorization comparability.** s(t̄) = C − s(t) exactly (C above; = THM-785's §2 transform
summed). Over all lines at n ≤ 6: the score multisets of the two ends are ALWAYS majorization-
comparable (no incomparable pair exists; blue n=6: 16 toward-regular / 12 level / 4 toward-
transitive; black: 240/108/132). So the flow is not merely variance-monotone (C3) — it moves
along CHAINS of the majorization lattice, with the blue layer 4:1 regular-directed at n=6.

**(d) Pure-blue census.** pure_blue(n) = ⌊(n+1)/2⌋ − [n even] = 2,1,3,2,4 for n=3..7 (confirmed;
n=7 via THM-785 §6's (4,84,184)). The pure-blue nodes are the flow's sources; at odd n they
include BOTH endpoints of the spectrum (transitive AND regular — n=5: {transitive, (0,2,2,2,4),
regular}), at even n only the transitive end plus interior nodes (n=6: {transitive, (1,1,1,4,4,4)}).

## THE PROOF of (a) — the even-n grid-sym stratum law (cont.6)

**Theorem.** Let n be even and let K be an n-tournament iso class containing a grid-symmetric
tiling. Then the total number of directed odd cycles of K is even, and H(K) ≡ 1 (mod 4).
More generally, any tournament admitting a FIXED-POINT-FREE involutive anti-automorphism has
c_odd even and H ≡ 1 (mod 4).

*Proof.* Three steps, each referee-verified exhaustively (mod4_blue_law_proof_referee_kps_S128c6.py:
n=4,5,6, all 1096 tilings, zero exceptions).

1. **(OCF mod-4, canon.)** THM-466 (= Grinberg–Stanley Thm 7.1): H = Σ_k α_k·2^k with α_k the
   number of k-element sets of pairwise vertex-disjoint directed odd cycles; hence
   H ≡ 1 + 2·α₁ (mod 4), α₁ = c_odd = total directed odd cycles.
2. **(Grid symmetry = reversal anti-automorphism.)** A tiling t is grid-symmetric ⟺ ρ(i) = n+1−i
   is an anti-automorphism of T(t) (u→v ⟺ ρ(v)→ρ(u)): ρ maps the base path to itself reversed and
   tile (x,y) to tile (n−y+1, n−x+1) reversed, so tile-map invariance of t is exactly
   anti-invariance of T(t). (Referee checked the equivalence both ways on every tiling.)
3. **(Fixed-point-free parity kill.)** σ(C) := reverse(ρ(C)) is an involution on the set of
   directed odd cycles of T (ρ anti ⟹ ρ(C) is a cycle of T^op ⟹ its reverse is a cycle of T;
   ρ² = id). A σ-fixed cycle has ρ-invariant support; at even n, ρ is a fixed-point-free
   involution, so every ρ-invariant vertex set has even size — no odd-size support qualifies.
   Hence σ is fixed-point-free on odd cycles and c_odd is even. By (1), H ≡ 1 (mod 4). ∎

**Corollary 1 (the blue mod-4 law, all even n).** Both endpoints of a blue line contain
grid-symmetric tilings, so both have H ≡ 1 (mod 4): ΔH ≡ 0 (mod 4).
**Corollary 2 (a pure-black certificate).** At even n, H ≡ 3 (mod 4) forces the class PURE BLACK
(n=6: all 24 such classes are PBk; the 44 PBk split 20/24 between residues, so it is one-way).
**Odd-n contrast (why the law dies).** ρ fixes the middle vertex; σ-fixed odd cycles (ρ-invariant
support through the centre) exist and c_odd parity is free: at n=5 the 16 grid-sym tilings split
8 even / 8 odd exactly, and grid-sym classes realize both H mod 4 residues.

The score-side (opus THM-790 leg law: reversal + leg defect) and this H-side stratum law are the
two halves of the d=m line's invariant story: degrees move by the legs; H's 2-adic second digit
(THM-466's frame) is frozen by the fixed-point-free reversal symmetry.

## Evidence log

- [x] n=4,5,6 exhaustive (transitivity_flow_merged_metagraph_kps_S128c5.py + .out): all laws as stated;
      independent reproduction of THM-785's interface law + node censuses before reading it
- [x] n=7 exhaustive (same script, 1569.7s): censuses (4,84,184) = codex's atlas (triple-certified);
      blue ΔH shows both mod-4 residues (odd n, as the proof predicts); blue majorization
      incomparables = 0, black = 268 (first failures); the four pure-blue classes identified
- [x] (a) PROVED — referee mod4_blue_law_proof_referee_kps_S128c6.py (+.out): OCF corollary,
      the grid-sym ⟺ ρ-anti-automorphism equivalence, even-n c_odd parity (64/64 even at n=6,
      4/4 at n=4; 8/8 split at n=5), stratum law + certificate tables
- [ ] proof target remaining: (c) blue majorization-comparability via anti-palindromic scores + s̄ = C − s
