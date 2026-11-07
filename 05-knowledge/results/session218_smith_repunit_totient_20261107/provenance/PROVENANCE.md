# Carmichael and the exponent fragment 14,26,34,38: provenance and lawful scales

Read-only source snapshot: `session216-stage`, Git HEAD `88ff4b0afd38315a4ad96db6164b80b55488e649`, clean at inspection. No source files changed. This is a synthesis of available repository evidence, not an independent replay of the full historical experiments or a literature-priority claim.

## 1. Resolve the four meanings first

1. **Carmichael's totient conjecture (OPEN):** no value of Euler's phi has exactly one positive-integer preimage. This is about fibers of an arithmetic function.
2. **Carmichael numbers:** composite squarefree N with p−1 dividing N−1 for every prime p|N. Equivalently every residue a satisfies a^N=a mod N. This is a universal modular-power predicate, not the totient conjecture.
3. **Carmichael lambda:** exponent of the unit group modulo N. For squarefree N, lambda(N)=lcm_{p|N}(p−1). Thus a Carmichael number has lambda(N)|(N−1), but lambda is neither phi nor the tournament function.
4. **Camion tournament functions:** s(n) counts strong labelled tournaments; c(n)=2^((n−1)/2) for odd n, c(n)=0 for even n≥2, c(1)=1. Set D=s−c. A Camion pseudoprime is a composite n satisfying n|D(n). Strong-singleton convention gives s(1)=1. This name comes from Camion's Hamiltonian-cycle theorem; do not spell it Carmichael.

## 2. Ambiguous prefix: nontotients first, a Camion overlap second

The fragment alone does NOT identify one sequence. In the context of Carmichael's totient conjecture, the natural candidate is the even nontotient sequence

14,26,34,38,50,62,68,74,76,...,

meaning even integers not in the image of phi. An independent exact sieve here verifies this prefix through 76: phi(n)>=sqrt(n/2), so any preimage of a target<=76 has n<=11552; sieving that complete universe suffices. This is zero-fiber data, whereas Carmichael's conjecture forbids singleton fibers: related but distinct predicates.

Repository provenance precedes the October Camion packet. The June29 reflection `07-reflections/tournaments-are-number-theoretic-the-paley-bridge-and-the-multiplicative-spine.md`, line31, explicitly says totients miss nontotients14,26,34; HYP-3541's historical index entry at `05-knowledge/hypotheses/INDEX-HISTORICAL-THROUGH-2026-07-21.md`, line2490, records the same multiplicative-gap analogy. This is provenance, not authority for its neighboring obsolete assertion that63 has no strong realizer: current THM-4884 records m(63)=2. The June6 canon `01-canon/theorems/THM-416-cm-density-quantum-totient-ladder.md`, remarks, already uses phi(m)=14 having no solution to derive max{m:phi(m)|14}=6. Import only this elementary claim, not every surrounding number-field assertion.

The October Camion exponent set below is a DIFFERENT object containing the four supplied numbers. Its full initial list differs before14 and at50. Direct recurrence gives D(123)=229 mod251, and 2^50−1=123 mod250 with251 dividing2^50−1. Thus50 is an even nontotient but not a Camion-passing exponent. Any overlap theorem requires proof beyond the four-term match.

### Camion exponent context and its correction

Primary canon:
`01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md`, notation, sections 4–7 and status.
Original research source:
`05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md`, section 2.
Independent historical audit:
`05-knowledge/results/tournament_pseudoprimes_audit_20261021.md`, Major 2, lines 123–133.
Correction authority:
`01-canon/MISTAKES.md`, MISTAKE-651, beginning at line 20991.

The mathematically defined exponent set is

E = {k≥2 : 2^k−1 is composite and (2^k−1) divides D(2^k−1)}.

The corrected FINITE-EXACT repository claim is

E intersect [2,64] = {4,6,10,12,14,26,34,38,62}.

Consequently 14,26,34,38 is a fragment of E, not its complete prefix. It is simultaneously the initial prefix of even nontotients, which is a more natural interpretation of the user's totient-conjecture context. It is not the initial sequence of Carmichael numbers or of attained totient values.

The earliest relevant located checkpoint is f002fc4253 (2026-10-21, identified by the mistake ledger). Its erroneous table asserted passes 4,10,14,26,34,38,62 and left k=59 undecided. The script had skipped all 14 squareful Mersenne numbers in this range, and capped reduced residues at 3000. Correction adds k=6 and12 (63 and4095 pass), rejects k=59 at prime179951 (reduced index29737, D residue46737 independently obtained by two methods), and removes the undecided status. The omitted squareful exponents were 6,12,18,20,21,24,30,36,40,42,48,54,60,63. 'Non-squarefree k' in some old prose really means non-squarefree 2^k−1, not non-squarefree exponent k.

No bounded prefix establishes an all-order classification. The known sufficient family does not classify E.

## 3. The exact half-scale and repunit bridge (PROVED algebra)

Let q≥5, M=2^q−1 and W=(2^q+1)/3 both prime. Then M=3W−2 and

N=2^(2q)−1=3MW.

The local residues are N=3 mod(M−1) and N=3 mod(W−1); the factor3 imposes no obstruction. Since D(3)=0 and the odd part is squarefree, THM-4746 proves N is Camion. Hence 2q belongs to E. The repository lists double-prime exponents q=5,7,13,17,19,31,61,127,... as examples, not as an independently verified infinite statement here. This explains the fragment through k/2=7,13,17,19, and adds k=10,62 within the audited range. q=3 is exceptional because W=3, N=63 is not squarefree; the mod9 continuation proves it passes. k=4 and12 are not supplied by the stated q≥5 double-prime theorem.

For repunits R_l(b)=(b^l−1)/(b−1), exact identities are

R_(2q)(2)=3R_q(4),
R_q(4)=MW,
W=Phi_(2q)(2) for odd prime q.

The scale k→k/2 changes a binary length2q repunit into a base4 lengthq repunit only after dividing its value by3. This distinction is essential: N=1023 passes Camion but R_5(4)=341 fails at31 because its reduced index is11 and31 does not divide D(11). Multiplication by3 moves the obstructing residue W to the free triangle residue3. Thus this is a residue-changing repair with a proof, not generic scale invariance. At q=3, R_3(4)=21 passes Camion while failing the base2 Fermat test; 341 passes the latter while failing Camion. These are hostile controls in both directions.

Moreover for q≥5 in this double-prime family, N is NOT a Carmichael number: N≡3 mod(M−1), so N−1≡2 mod(M−1), with M−1>2. The family therefore separates Camion from Carmichael rather than identifying them.

A direct totient bridge also exists at the exponent level: for every odd q>1,

phi(2q)=phi(q), with 2q≠q.

So 14↔7,26↔13,34↔17,38↔19 are explicit equal-totient pairs; their common values are6,12,16,18. This holds for every odd q, not because of Mersenne/Wagstaff primality. It resolves these particular fibers and gives no progress on the unresolved general conjecture. Keep three separate coordinates: exponent k; integer value2^k−1; and function value phi(k).

For arbitrary bases the safe starting identity is R_(uv)(b)=R_u(b)R_v(b^u). It proves divisibility/factorization, not Camion preservation or Smith digit-sum equality. A length/base change needs the full factor and residue sidecar. The binary factor3 arises specifically from b+1=3 and the fixed tournament free residue3; replacing2 by any b is not justified by the binary proof.

## 4. Current Camion mechanism and continuation

THM-4746 proves, for odd primes p and m≥1,

D(m+p−1) ≡ (2/p) D(m) mod p.

For a composite n with squarefree odd part, this yields the iff criterion p|D(r_p) for all odd p|n, with r_p the representative of n mod(p−1) in1..p−1. The 2-part is automatic by a free rotation action. D(1)=D(2)=D(3)=0 supplies every Carmichael number, every2p and every3p. None of these assertions concerns the totient conjecture.

The original false claim that this periodicity was unique to linear strong-component condensation was repaired by MISTAKE-651: connected labelled graphs obey the same law from m≥2. Its real mechanism is a finite exponential family closed under Hurwitz inversion and shift. Rooted counts m*s(m), not connected graphs, give the appropriate different-period control.

The same correction replaces a false mod-p² nonperiodicity claim with quasi-period p(p−1) from m≥3. Session25 continuation gives quasi-period phi(p^a) modulo p^a from m≥2a−1, allowing unrestricted prime-power CRT. Source route:
`05-knowledge/results/session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md` and linked review (recovered read-only with git show HEAD:path from the sparse checkout).
The current canonical summary records prime-power towers as initial segments with frozen p-adic valuation, and its finite p≤10^4 computation leaves only prime powers9 and169 beyond primes. This is not an all-prime assertion: later THM-4750 records Camion-Wieferich prime22409 through10^5.

Further current canon:
`01-canon/theorems/THM-4750-the-camion-function-half-speed-interpolation-and-the-camion-lift.md`.
Here 'half-speed' means p-adic Mahler divisibility rate ceil(i/2), arising from the quadratic edge exponent binom(n,2). It is NOT the exponent halving k→k/2 above. Its Lambda_p and Camion lift are not Carmichael lambda.

## 5. Earlier Carmichael-number work

`05-knowledge/results/carmichael_successor_blocks_20261014b.md` gives the verified interface (N,L,e), L=lambda(N), e=(N−1)/L, together with distinct prime factors. The exact one-new-prime classifier is

Nq Carmichael iff q=1+dL is prime, q not dividing N, d a positive divisor of e.

The update is (N,L,e)→(Nq,dL,N+e/d). This is a rigorous scale N→Nq with explicit preservation guards; arbitrary multiplication is not allowed. The same note gives a finite complete two-prime classifier and catalogs: root561 has no one-prime extension but exactly one two-prime extension561*401*3041=684106401; neither intermediate product is Carmichael. Root1105 has0/13 one-/two-prime children;1729 has5/21. These are reported historical results, not rerun here.

Its finite product-one compiler keeps coprime L,k, distinct primes p=dk+1 with d|L and p not dividing L, and a chosen product N≡1 mod L. Then CRT gives N≡1 mod kL, sufficient for Korselt. Dropping coprimality fails at L=6,k=2, primes7,13:91≡1 mod6 andmod2 but13−1 does not divide90. This is a strong reusable 'preserve conductor and factorization' bridge to all-base repunit experiments.

`05-knowledge/results/apollonius_enumerative_20261019_galois_parity_cube_ford_descartes_windows.md` maps a user-supplied October2026 Carmichael-counting paper into character-product methods and the (3/2)-pseudoprime sequence A073631. Its header says pre-audit. Do not promote the attached paper's claimed C(X)=X^(1−o(1)) or its transferred exponent as an independently established external theorem. The elementary implication itself is safe: a Carmichael N coprime to6 satisfies3^(N−1)≡2^(N−1)≡1 mod N, hence N divides their difference. The converse fails; this does not characterize Carmichael numbers.

## 6. Totient-inspired tournament fibers, corrected current state

Primary current source:
`01-canon/theorems/THM-4884-definability-lenses-korselt-congruence-of-spectrum-polynomial-strong-multiplicity-census-minimal-rank-sublanguages.md`.
Detailed provenance:
`05-knowledge/results/definability_lenses_spectrum_diophantine_carmichael_20261105.md`, especially sections4/4b.
Correction: `01-canon/MISTAKES.md`, MISTAKE-723, lines23089–23095.

There are TWO different fibers:
- m(h): strong tournament ISOMORPHISM CLASSES, over all orders, with Hamiltonian-path count h.
- pf(h): distinct numerical source PROFILES (D,E,U,d), reachable by O/Z/J from(1,0,2,1), with C1 readout D−E+U=h. This D is a profile coordinate, not the Camion difference s−c.

The source moves are O(D,E,U,d)=(D,0,2U,d), Z=(D−E+U,d,2U−d,d), J=(D−E+U,D−E,D+U,D). Profile, word, labelled-object and isomorphism-class multiplicity are not interchangeable.

Reported FINITE-EXACT census for odd3≤h<625:
- m(h)=0 exactly at7,21;
- m(h)=1 exactly at3,5,11,13;
- m(h)=2 exactly at9,17,19,35,39,47,63;
- m(h)=3 exactly at25;
- every attained h in(63,625) has m(h)≥6.
The earliest Session24 inventory and the S50 multiplicity note already provided the data below73; the November5 lens interprets it and extends via THM-4744's inventories. The pre-audit checkpoint841dc8e7f6 incorrectly gave m(41)=3; correction886693f53a restores m(41)=6. Never propagate41 as a singleton or triple fiber.

pf has singleton values3,11 through50001, whereas m has singleton5 too (two profiles, one class), and m(13)=1 while pf(13)=0. The source reports exactly29 C1 omissions through50001, largest2647; this is finite coverage only, not a cofinite theorem. Session200's primary source `session200/carmichael-fibers/PROOF.md` was recovered by git show. It proves explicit physical twins F(O(s))=F(J(s))=D+2U, with distinct endpoint coordinates, and dyadic continuation splitting by (2^(k+1)−2)(U−d). Its two-run words O^a Z^b have F=2+(2^(a+1)−1)(2^(b+1)−1), a direct Mersenne-product bridge. It explicitly imports the odd-input totient pairing and isolates the rigid nonbalanced chamber rather than claiming global multiplicity.

Important quantifier repair: the note loosely says a Ford-style 'every multiplicity≥2 occurs' analogue fails because multiplicities9,13,21,22,34,37 are absent below625. That inference is not proved by a bounded census. Only absence in this target range is established. In contrast, an exact singleton such as m(3)=1 genuinely refutes the universal no-singletons analogue. m(1) needs a separate singleton convention, so the summary above explicitly starts at3.

The orbit-counting Korselt lens of THM-4884 is a real bridge: the spectrum polynomial P_n(x)=sum_T x^H(T) reduces coefficientwise mod p to the sum over tournaments fixed by rotation n/p. It supplies congruences, not a positive-spectrum coverage theorem or a totient-fiber identity. Its mechanism predates the November5 packet (THM-4746, THM-4747, LEM-003).

## 7. Recommended decisive experiments and guardrails

1. Store a typed row (q,k=2q,binary N,base4 R_q(4),factorization,local Camion residues,phi(k),phi(q)) rather than aligning bare prefixes. Include q=3,5,7,11 and k=12 as exceptional/hostile controls.
2. Separate a double-prime-family flag from the full E predicate. Any extension of E must check every prime power in the factorization, report undecided residues as undecided, and preserve threshold data.
3. To test Smith numbers, add base, digit sums and prime factors with multiplicity. Neither phi-pairing nor Camion congruence determines the digit-sum predicate. Use carry accounting rather than a presumed scaling law.
4. For Carmichael scaling use the exact (N,lambda,e) classifier; for Camion scaling recompute local zero classes; for exponent rescaling use exact repunit identities. These are three different operations.
5. Keep LRC(14), the {7,21} spectrum problem, and the general totient conjecture OPEN. Shared integers or digit lengths do not transport their predicates.

## Sparse-source verification and limits

The checkout is sparse: absence of an on-disk path was NOT absence from the repository. `git ls-tree HEAD` confirms and `git show HEAD:path` recovered the Session200 source-fiber proof, Session25 additive-prime-power proof, and historical audit `04-computation/experiments/tournament_pseudoprimes_20261021/audit/a5_mersenne.out`. The latter records the corrected k<=64 decisions, including50 failing at251 and59 failing at179951. Reading those receipts is not an independent reproduction of the full historical census. This artifact independently verifies the small even-nontotient prefix and the k=50 local obstruction, and directly checks the displayed elementary algebra. The source worktree remains unchanged.

## Audited stronger bridge

The [scaling proof](../totient/PROOF.md) and [independent core audit](../audit-core/CORE_REVIEW.md) establish that simultaneous Mersenne/Wagstaff primes at q>5 force 2q to be an even nontotient. This supplies a genuine common subfamily, not a prefix coincidence. The exception q=5 remains: 10=phi(11), while 2^10−1 passes Camion. The [all-base classification](../smith/BASE_COMPLETENESS.md) separately supplies a Smith base for every even nontotient, but that compatibility is too broad to characterize this sequence.
