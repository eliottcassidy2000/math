# First occurrence: audited inheritance and constructive transfer

Read-only source snapshot: session222-stage commit 59ff4cd82c66. Parent coordinates checking later remote changes. No repository edits or pushes by this worker.

## 33: exact prior work, independently checked

THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md routes to the 2026-10-27 atlas. Its frozen n6 output gives five strong isomorphism classes, 2,640 labelled tournaments with H=33. Order5 maximum is15, so minimum order is6. The complete strong inventory is P33(z)=5z^6+3z^7; floor at8 is45, so no larger strong class exists.

Independent standard-library replay check33.py enumerates the full labelled n6 universe via S6 orbits and evaluates each representative by subset DP. Canonical masks in increasing pair order (0,1),(0,2),..., bit1 means i→j:
- 88 and377: converse pair, automorphism order3 each.
- 156 and282: converse pair, rigid each.
- 145: self-converse, rigid.
Therefore33 has three converse-orbits at birth: one fixed class and two pairs. Equal H and equal score sequences do not classify these objects: masks88 and156 have scores(2,2,2,2,3,4) but different automorphism orders.

Prior audited S50 source 04-computation/experiments/h_spectrum_multiplicities_20261021/h_multiplicities.out gives N(n,33)=n²−5n−1 for every n≥6. This follows from SCC uniqueness:
F33(z)=(5z^6+3z^7)/(1−z)^2+2z^8/(1−z)^3.
The second summand is the two orders of the unique strong H3 and H11 classes, of orders3 and5, with singleton gaps. Thus first-order census is distinct from total strong census (8 classes), and from arbitrary-order census (unbounded).

Sources: THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md; 05-knowledge/results/h_spectrum_multiplicities_20261021.md; 04-computation/experiments/h_spectrum_sweep_20261027_n6.out. Historical gn-blue-black-structure.md mentions the Aut3 converse pair but is provenance, not the proof.

## Spectrum status and corrections

THM-1370-h-spectrum-omits-7-21-all-n.md proves7,21 impossible at every order using SCC multiplicativity, insertion monotonicity, Moon's strong deletion theorem, and finite small spectra. Its priority correction credits the2016 MathOverflow discussion; do not claim a new exclusion theorem. The exact full range odd positive integers minus{7,21} remains OPEN.

Chan–Pak, Computational complexity of counting coincidences, TCS1015(2024)114776, arXiv:2308.10214v2, §6.5 Conjecture6.2, records the proposed range, crediting bof and Gordon Royle(2016). The source audit 05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md clarifies Royle's coverage is strictly below80,557, safely inclusive80,555, not80,557. URLs: https://arxiv.org/abs/2308.10214v2 and https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament . Strong realization is a stronger variant and should not be silently equated with the ordinary range conjecture.

THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md gives finite-exact audited coverage to2^32. Its addendum gives2^35 but explicitly says one run, not independently recomputed at2^35. Neither closes Chan–Pak.

## Constructive mechanisms (elementary proofs)

1. For ordinal sums A→B, every Hamiltonian path visits all A before all B; hence H(A→B)=H(A)H(B). Converse reverses component order and converses each component. Adding a source/sink preserves H but cannot create a first occurrence for h>1.

2. Let m(h) be least order and s(h) least strong order, with infinity for unattained values and excluding singleton h1. For h>1,
m(h)=min over ordered factorizations h=h1⋯hk, hi>1, of Σs(hi).
Proof: the SCC condensation is a total order and H is the product. A minimal witness cannot contain a singleton SCC; replacing any nonminimal strong component shortens it. Conversely ordinal sums attain every candidate. This gives a min-plus construction algorithm from the strong inventory; it does not imply uniqueness of minimizers. Prime h forces a strong minimal witness. For33 the disconnected3·11 candidate costs3+5=8, while strong witnesses cost6.

3. General substitution needs more than H. For each block B let c_k(B) count ordered k-tuples of nonempty vertex-disjoint directed paths covering B. For a quotient tournament Q, let W_Q(k1,...,kr) count words with ki appearances of letter i, no equal adjacent letters, and every consecutive distinct pair following an arc of Q. Then
H(Q[B1,...,Br])=Σ_k W_Q(k) ∏_i c_ki(Bi).
Proof: uniquely split a Hamiltonian path into maximal consecutive block-runs. Their block letters give the quotient word; each block's runs, in occurrence order, give the ordered path cover. This is a bijection. In a transitive quotient each letter can occur only once, recovering product law. In a cyclic quotient reentry is possible and the higher path-cover coordinates are indispensable. This is a transfer principle, not a claim of novelty.

4. On every birth fibre converse is an involution, so its class count is a+2b, with a self-converse classes and b nontrivial pairs. No involution theorem forces a+b=1. At33 the exact pair(a,b)=(1,2). At23 and29 the prior n6 atlas already has4 classes at first occurrence. Even at15 there are two self-converse classes at first occurrence (the n5 maximum), so 'exactly two classes' need not mean a distinct converse pair. The right atlas retains the full pair(a,b), not a binary label.

5. Moon interval chains provide a more structured family: for composition p with parts≥2, H(IC(p))=∏(2^(p_i−1)+1), distinct compositions give distinct classes (sum p≥3), converse reverses p, self-converse iff p is palindromic. This is THM-4744, with product law credited Moon1972 and finite extremal classification credited Moon–Yang2007. Commutative numerical product forgets the noncommutative composition word: the precise mechanism for many classes with the same H. A construction can certify an upper bound on m(h), but a first-occurrence claim additionally needs exclusion of all smaller orders.

## Recovered all-order converse refinement (already prior work)

The complete SCC/square-divisor/Burnside grammar is explicitly in 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md, equations(1)–(3); it credits Sessions22/23 for the architecture. The sibling session24-multiplicities-review/REVIEW.md states APPROVED mathematically and computationally, with independent full augmentation replay and all low-H controls. PROOF.md's initial 'review still required' banner is stale relative to this review.

For strong ordinary polynomials P_m, strong self-converse Q_m, ordinary all-order F_m, and self-converse all-order D_m:
D_m(z)=sum over d²e=m of F_d(z²)([e=1]+Q_e(z)), taking Q1=z and F1=1/(1-z).
The first half of the SCC palindrome is an ordinary tournament, not a converse-orbit. For33, the independent S50 audit inventory_375.json has n7 strong rows(code,H,self-converse,Aut): (1025,33,true,1), (133188,33,false,1), (1344001,33,false,1). This is the auditor's different code convention, not our lex-mask convention. Together with the n6 replay it gives Q33=z6+z7, hence D33=z6/(1-z), exactly one self-converse33 class at every order≥6. Nontrivial converse-pair count is(n²−5n−2)/2 and total merged count is(n²−5n)/2. These are recovered prior results, not new theorems from this session.
