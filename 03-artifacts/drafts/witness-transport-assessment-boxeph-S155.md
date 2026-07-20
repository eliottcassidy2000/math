# Is Witness Transport Novel and Meaningful? — Assessment and the Proof/Lean Queue

**boxeph-2026-07-20-S155 (HYP-8240).** Owner question: JC was disproven
externally by a single counterexample. Is transporting F through the de
Bondt–van den Essen symmetric reduction to explicit VC/IC/Mathieu witnesses
novel and meaningful? Should it (and which other repo results) become rigorous
proofs + Lean writeups? "Think Hadwiger."

## 1. Verdict on witness transport: YES — novel, meaningful, publishable standalone

**Novelty (checked):** the equivalences JC ⟺ VC (Zhao 2007), IC ⟹ JC (Zhao
2010), Mathieu-subspace framework (Zhao; de Bondt–van den Essen 2005 symmetric
reduction) are truth-transfer theorems. With ¬JC external, the literature now
knows ¬VC and ¬IC *abstractly* — but NO explicit witness exists anywhere
(S149 sweep: untouched in-repo and unclaimed in the fleet; the PROBLEM-LEDGER
lists A4/A6/A7/A7b as targets). The first explicit polynomial P with
Δ^m(P^m) = 0 for all m but Δ^m(P^{m+1}) ≠ 0 for some computed m would be a
new mathematical OBJECT, not a corollary.

**Meaningfulness — the owner's Hadwiger point, made precise:** de Bruijn–Erdős
says the chromatic number of the plane is witnessed by a finite graph;
that abstract witness-existence sat inert for decades, and it was de Grey's
EXPLICIT 1581-vertex graph that moved the field (Polymath16, record chasing,
new structure theory). Equivalences are normally run truth-only; running them
witness-effectively (i) produces the first concrete failure objects of the
Mathieu program, (ii) measures the CONSTRUCTIVE COST of the equivalences —
the final dimension N and the degree of P are new constants quantifying "how
far" VC is from JC in witness distance, and (iii) tests the thesis the owner
named: equivalence proofs preserve more than truth values — they are (or fail
to be) witness functors. If some step of de Bondt–van den Essen turns out
witness-ineffective (e.g., needs a generic linear change with no canonical
choice), THAT is a finding about the equivalence itself.

**Publishability:** standalone paper: "Explicit counterexamples to the
vanishing conjecture, the image conjecture, and a Mathieu-subspace property."
Content: the transport chain (Keller F → Yagzhev cubic-linear form → gradient
map F = x − ∇P, P homogeneous quartic → VC witness → IC/Mathieu witness),
each step effective, each intermediate object printed, all identities
machine-verified. Risk (stated): dimension blow-up — Yagzhev linearization
grows variables with monomial count; expect N in the tens-to-hundreds; the
Δ-power check needs one explicit m — effective but possibly large. Mitigation:
do the transport mod p first (cheap dry run pinning N and degrees), then exact.

## 2. The same question for the rest of the repo — the proof/Lean queue

Criteria: (a) is the statement NEW (not a re-derivation)? (b) does it survive
the external ¬JC news (i.e., not merely "JC is false", which is no longer
ours to claim)? (c) is it Lean-feasible with today's repo patterns (citation
hypotheses for classical inputs, kernel-pure finite checks)?

RANKED QUEUE:
1. **THM-1370 ({7,21} h-gap, all n)** — new statement (web-checked), clean
   proof (multiplicativity + insertion lemma + Moon [cite] + finite floors).
   Lean: HIGH feasibility — the n ≤ 6 exhaustion is kernel-decidable, the
   insertion lemma is constructive, Moon enters as citation hypothesis (the
   LRC(≤13) pattern). PLUS the completeness half (§3.1) would make it "the
   h-spectrum is exactly odd ∖ {7,21}" — a complete characterization paper.
2. **Witness transport (this memo, §1)** — proof = the executed chain; Lean =
   identity checking (ring/decide), staged: det JK = −2 + collision (small,
   do FIRST — the first formal verification of the JC-disproof input), then
   the DC endomorphism identities (THM-1300/S141), then the VC witness.
3. **Mod-27 rung → LRCMod27Spread.lean** — direct mirror of Mod19/23/25
   files; the ladder-completion statement (rungs ⟺ pairs(q) ≤ 13) is a
   two-line arithmetic lemma on top. LOW effort, immediate.
4. **Anti-order laws (S150/S151)** — anti-orders ≡ 2 mod 4 (Moon [cite] +
   two lines); the 2μ(m) threshold (orbit-pairing + explicit construction);
   |Aut| divides h (free action). All elementary, all Leanable, together a
   tidy paper: "Dihedral rigidity of tournament symmetry."
5. **THM-1365 Galois reduction** — Lean with Kambayashi + Keller-birational
   as citation hypotheses; the Cartan-freeness lemma is formalizable. The
   FP(n) bridge statement is the paper's hook.
6. **Sliver = BU-index threshold (S153)** — LRC-community note; Lean of the
   rung-termination arithmetic is trivial; the framing is the value.
NOT recommended for standalone write-up: the JC counterexample verification
itself (external priority), the Moon-Busch floor (known theorem — our
contribution is only the exhaustive confirmation + the {7,21} use).

## 3. Attackable targets (community-inspiring), from this session's vantage

3.1 **Spectrum completeness** (finish THM-1370 into a characterization):
prove every odd ≥ 23 is attained. Route: strong(9) already covers
[75, 2881] ∖ {77,79,83,...tiny}; padding + products cover below; need ONE
constructive family filling [23, 75) gaps (25,45,75 strong floors + composites
27,33,35,...: remaining at n ≤ 9: NONE below 2883 except 7,21 — so the
missing piece is a family showing every odd ≥ 23 stays attained for ALL
larger... it already does: padding is monotone. The ONLY open part is odds
≥ 2883 at their first n. Conjecture: strong(n) ⊇ [f(n), c·3^n] with holes
only near the ceiling; prove strong-interval density via block chains +
single-arc toggles. Effort M; payoff: complete characterization.)
3.2 **Reflection-equivariant h-monoid** (S150 row 7): anti-Rédei ⟷
odd-degree-conjecture transfer — if the h-monoid dictionary extends to
reflections, the odd-degree conjecture for Keller counterexamples gains a
tournament-side proof template. Effort L, payoff high.
3.3 **ind_Z/2 of the ×7 gate tower** — the p=43 ±-section is now verified
(S155 (1)); next: compute the index of the band section in the 7-torsion
filtration — reprices S-T wall (b). Effort M.
3.4 **HS(B) proximity structure** — new object (S154), rib-concentration law
(S155 (2)): even-class disjointness concentrates on SC-NS pairs (n=5: 100%
of NS-NS pairs share an even class; n=6: enrichment 42% vs 34%). First
attack: prove the n=5 statement (NS-NS ⟹ shared even class) from the
cycle-space bijection. Effort S.
3.5 From the niche sweep (agent report appended when it lands — network
retried once).

### 3.5 The niche sweep (agent-mined; bodies marked UNVERIFIED were not re-read this session)
Tier A (verified computationally in-repo, unproved, standalone-paper-shaped):
cpA=>cpK spectral conjecture (T1546; 0/116 splits; rank-one resolvent attack, M);
Moser 32nd-region deficit law + A362193 Grassmannian bridge (T1533, M);
Sequence D not-in-OEIS + A060488 tricovering diagonal (T1532, S/M);
Faulhaber shear constant (Lambert-W?) + 4 unsubmitted OEIS sequences (T1554, S/M);
1001 = three sixties (ord_1001(2) = Pisano(10) = |A5| = 60; T1534, S — quotable).
Tier B (bridges to known problems): Erdos #592 axis R(n,2) = 2n+1 (T831, $1000, L);
Alcuin escapes minor-closure — finite forbidden-subgraph set without WQO (T828/829, M);
Heegner<->Rabinowitsch dictionary (T718/719, M, UNVERIFIED); toothpick A139250 corner
recursion (T1545, M); n=24 Niemeier rung open (T1534, M/L); GLMY-Paley court case
(HYP-8230, M); Farey-Fibonacci power-payload flip (T1083, S/M).
Tier C (one-offs, UNVERIFIED): Cayley-Delannoy near-paper (.tex — closest to
submittable), Grinberg-Stanley note, Napolitano gauge draft, hard-core/idea-gas
H = 3^alpha, 987-amplituhedron, label-drift coincidence, adelic geometry.
RECOMMENDED next niche session: (1) Alcuin forbidden-subgraph set k=3;
(2) Sequence D + shear-constant OEIS submissions (cheap, visible);
(3) cpA=>cpK proof attempt; (4) read + triage the Cayley-Delannoy tex.
