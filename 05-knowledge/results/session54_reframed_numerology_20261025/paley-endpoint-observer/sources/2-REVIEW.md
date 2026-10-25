# Session 25: independent congruence and symmetry audit

## Subsequent additive closure (Session25)

A separate argument, now fully written and independently audited in ADDITIVE_PRIME_POWER_PROOF.md, proves s(n+phi(p^r))=(2/p)s(n) modulo p^r for every odd p, r>=1, and n>=2r-1, with the same law for D. Thus the incoming additive conjecture is now proved by this new filtered-Hurwitz/Newton argument. The earlier caution that the multiplicative theorem alone does not settle it remains correct. The new proof also gives arbitrary-unit-base closure and the sufficient connected-graph bound n>=2r. No minimal-threshold or literature-priority claim is made.

## Current reconciliation: incoming audit takes precedence

Incoming commit `3596c5d0b1990b05ae960e980a863956b47cc6b5` and MISTAKE-651 supersede the initial pre-audit corrections below. Its revised THM-4746 already repairs the connected-graph control, proves the prime-square equivalence by additive quasi-periodicity modulo p^2, fixes the Mersenne table (including 63 and 4095), and qualifies Z_3. Revised THM-4747 already records 262 classes and classical rotation credits. Preserve those changes; the old proposed patch is withdrawn and retained only as HISTORICAL_SUPERSEDED_CANON_CORRECTIONS.txt.

The connected-graph correction is stronger than our initial finite-scope caution: connected labelled graphs DO satisfy the same quasi-periodicity for every m>=2. Thus ordered strong-component condensation is not the distinguishing cause of quasi-periodicity. See INCOMING_AUDIT_RECONCILIATION.md for a direct contraction verification and incoming source credit.

Our general multiplicative prime-power law supplements, rather than replaces, incoming Lemma B. It does not settle the incoming additive higher-power conjecture. Larger incoming numerical searches and external bibliography claims are source-attributed, not independently rerun or reverified here.

## Scope and verdict

Target: immutable incoming PRE-AUDIT commit f002fc425336; exact targets:
- THM-4746, `01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md`.
- THM-4747, `01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md`.

The main congruences and group-action proofs pass independent mathematical audit. The scope/census corrections below record the original pre-audit review; incoming MISTAKE-651 now supplies the repaired canon for overlapping items. The entire long frontmatter is NOT thereby independently verified. Repo remained read-only; reproducible independent artifacts are this directory's audit.py and audit.out. Original experiment files were read from the immutable git object because the working checkout does not materialize those paths. No literature novelty claim is made; these are elementary orbit arguments, and a bibliography search was not performed in this audit.

## 1. A direct replacement proof of the strong-count congruence (PROVED)

Let p be an odd prime and m >= 1. On m+p-1 fixed labels, let g cyclically permute a selected p-element block B and fix all remaining m-1 labels. Let s(n) count strongly connected labelled tournaments, with the conventional s(1)=1.

The order-p action on the finite set of strong tournaments gives s(m+p-1) congruent modulo p to its fixed-point count. A fixed tournament has:
1. a circulant tournament on B, with exactly 2^((p-1)/2) choices;
2. constant adjacency between each outside vertex and B;
3. arbitrary tournament arcs between the outside vertices.

Every circulant tournament is strong: its condensation is a linearly ordered set of strong components; automorphisms act trivially on that line, whereas translations act transitively on the vertices, so there can be only one component. Alternatively at prime order any nonzero step already produces a directed Hamiltonian cycle.

Contract B to one distinguished vertex. Strongness of the expanded tournament is equivalent to strongness of the resulting m-vertex tournament, because B is internally strong and all outside adjacencies are constant. This remains valid at m=1. Thus the exact fixed-point count is 2^((p-1)/2)s(m). Euler's criterion gives

s(m+p-1) = (2/p)s(m) mod p.

The same relation for c(m) and D(m)=s(m)-c(m) follows immediately. No formal power series division, finite-field extension, or group-ring injectivity assumption is needed. The original group-ring proof also appears valid; its homomorphism need not be injective.

Boundary discipline: m=0 is excluded and cannot simply be added; p=2 is excluded; s(1)=1 is the strong-singleton convention, not an ordinary loop-free Hamiltonian cycle count. Camion's Hamiltonicity identification needs the usual nontrivial-order qualification.

## 2. Korselt criterion and exact 341 comparison (PROVED)

For n composite with squarefree odd part, each odd prime p dividing n gives D(n)=epsilon^k D(r_p) modulo p, where 1 <= r_p <= p-1 and n=r_p+k(p-1). The scalar epsilon^k is invertible. For even n, the full 2-primary rotation subgroup acts freely on tournaments: any nontrivial subgroup contains the half-turn, which would reverse an edge. Thus the entire 2-part of n divides s(n), while c(n)=0. CRT proves the stated iff criterion.

Independent modular Wright-recurrence results:
- n=21: s(n)=16, c(n)=16 modulo21.
- n=341: s(n)=23, c(n)=1 modulo341, hence D(n)=22 modulo341.
- r_31=11 and D(11)=22 modulo31. The residue-class index is11; the discrepancy value is22. These must not be conflated.
- n=1023: s(n)=2, c(n)=2 modulo1023.
- n=9: s(n)=c(n)=7 modulo9, explicitly covering the exceptional p=3 endpoint of “every3p”.

The families 2p, 3p (p odd prime), Carmichael numbers, and 3b(3b-2) with b>=5 and both displayed factors prime are valid. For 3p with p>=5, use CRT; p=3 is the separately checked finite n=9 case. The Mersenne–Wagstaff conditional family is valid when its hypotheses and b>=5 hold; q=3 lies outside that family and gives63, which must be treated separately if included. Small-prime “Z_p contains {1,2,3}” language should instead restrict to the representatives that actually lie in1..p-1; at p=3 the representative3 does not exist.

D(r)>0 for r>=4 is justified because all circulants are strong and a nonregular strong tournament exists: substitute a transitive tournament of size r-2 into one vertex of a directed triangle. This is strong and nonregular.

## 3. Individual cycle and path actions (PROVED)

For an individual tournament T, Aut(T) acts freely on its ordered directed Hamiltonian paths, since a fixed ordered path fixes every vertex. Consequently |Aut(T)| divides H(T).

This is NOT true for unrooted directed Hamiltonian cycles. Their stabilizers inject into the cyclic rotation group C_n. A nonidentity stabilizer element is semiregular and has no fixed vertex; a prime-order automorphism with a fixed vertex therefore forces its prime to divide C(T).

For T=Cay(Z_p,S), p odd prime, and affine subgroup A of order p*mu:
- a nontranslation affine map fixes a vertex, so fixes no directed Hamiltonian cycle;
- a cycle fixed by a nonidentity translation is exactly an arithmetic step cycle;
- there are (p-1)/2 such cycles, each stabilized by the full translation subgroup and no nonidentity multiplier;
- all other A-orbits of cycles have size p*mu.

Therefore C(T)=(p-1)/2 modulo p*mu. The proof does not require A=Aut(T).

The signed affine action on paths also passes: every path-fixing anti-automorphism is an involution; the only such affine elements are the p reflections. Their fixed-path counts are equal. Burnside, free ordinary action, and odd H give rho_0=mu modulo2mu. This cannot be obtained by pretending the signed action is free.

The Wilson pair-count proof passes for odd prime p. Prime p=2 must not be inserted into the displayed tournament pair formula.

## 4. Exact bridge to metagraph presentations (PROVED, with quotient scope)

Let X_cl be the set of closed Hamiltonian paths, namely paths whose last vertex points to their first. Each directed Hamiltonian cycle supplies exactly p such paths. The free A-action therefore gives

|X_cl/A| = p*C(T)/|A| = C(T)/mu = (p-1)/(2mu) modulo p.

This is an exact closed-presentation congruence for the A-quotient. If and only if the quotient in use is by A=Aut(T), it is the genuine closed metagraph-presentation count. If A is proper in Aut(T), do not silently identify these quotients.

For a prime Paley tournament with the usual verified equality Aut(T)=A, the closed-presentation count is1 modulo p. For the supplied small cases:
- P3: total presentations1, closed1, open0.
- P7: total9, closed8, open1.
- P11: total1729, closed1101, open628.
These follow from C=1,24,5505 and mu=1,3,5. They refine the rotation/closure accounting of THM-4741 and THM-4729, not an aggregate s(n) statement.

There is also a canonical involution on H(T)/A induced by reflection. Its fixed A-orbits number rho_0/mu, an odd integer: count the p reflection-fixed path sets and divide by |A|, using uniqueness of any path-fixing anti-automorphism. Again this is a signed quotient statement, not a free signed action.

General obstruction: an automorphism cycle stabilizer is cyclic and semiregular, with order dividing n. Thus an automorphism subgroup G whose every nonidentity element fixes some vertex acts freely on Hamiltonian cycles, giving |G| divides C(T). This is a useful symmetry obstruction for any proposed metagraph cycle-coordinate compression.

## 5. Independent computation and concrete corrections

Reproduce from workspace root: python3 session25-congruence-review/audit.py

Fresh Wright-recurrence implementation checked 24,000 congruences: all odd primes p<=97 and every1<=m<=1000. A separate subset DP counts paths/cycles, and direct signed-pair permutation enumeration counts reflected paths, on all62 labelled connection sets at n=3,5,7,9,11. All required congruences passed.

The source checker reads264 data rows, but these represent only262 distinct multiplier classes. The separate Paley19 and Paley23 checks duplicate classes already in the census; the Paley23 connection set is a different multiplier-equivalent representative. Counts by odd order3..23 are1,1,2,4,4,6,16,16,30,88,94. Therefore “all264 circulant classes” should be “264 checked rows representing262 multiplier classes.” At composite orders, multiplier classes are not automatically full isomorphism classes; retain that qualification.

HISTORICAL PRE-AUDIT FINDING, now superseded by incoming Lemma C: the theorem headline's broad connected-graph claim was not proved by its evidence. The actual control output checks only p=5,7,11,13,17,19,23 and only periods/antiperiods L<2p over its tested range. Do not promote this to a general no-period theorem. The incoming audit now proves the positive tail law for m>=2; the initial obstruction is confined to m=1. The semicomplete extension's group-ring proof works for odd p!=3; p=3 needs a separate argument (the recurrence gives s_3(n)=(-1)^(n+1)n! mod3, hence the claimed zero multiplier for n>=3). Its arithmetic is a control, not a new graph-structural equivalence.

The prime-square search through p<=3000, full high-order Paley numbers, and all large exponent lists were not independently reproduced by this audit. They should remain source FINITE-EXACT or separately audited. UPDATE: the formerly unaudited equivalence of D(p)=0 modp^2 with p^2 passing is now PROVED for every odd prime by the separately supplied and independently audited PRIME_SQUARE_PROOF.md. It follows from D(p^2)=(2/p)D(p) modulo p^2, not from the prime-modulus argument alone.

## 6. Separation ledger and novelty

s(n) counts strong labelled tournament objects across an entire universe. H(T) counts paths inside one tournament; C(T) counts its unrooted directed Hamiltonian cycles; H(T)/Aut(T) counts presentations. These are different types. The contraction map preserves strongness and labels, but destroys internal path multiplicities, endpoint data, stabilizers, and the metagraph wiring. Restoring those requires a path/endpoint/stabilizer sidecar, not merely D(n).

The ordered strong-component decomposition remains a genuine tournament factorization mechanism, and K1/C3 are genuine boundary atoms. It is not a distinguishing explanation of quasi-periodicity, which incoming Lemma C also proves for unordered connected graphs. The integer21 passing the aggregate congruence does not explain H-spectrum gaps7/21, establish S2/S3 spectral identities, or identify any individual wiring. The exact usable new connection here is the closed-presentation orbit count and its stabilizer obstruction. Novelty remains unchecked, and no claim of new named pseudoprime theory or literature priority is justified by this audit.

### Supplementary prime-square hostile probe (FINITE-EXACT only)

Direct independent recurrence at every odd prime p<=31 confirms that D(p) and D(p^2) vanish modulo p^2 simultaneously, exactly at p=3,13 in this range. In fact the checked residues obey D(p^2)=(2/p)D(p) modulo p^2. These raw pairs are in audit.out. They supported the initial probe; the subsequent PRIME_SQUARE_PROOF.md now proves the displayed lifting relation and equivalence for every odd prime. The finite search bounds elsewhere remain separately qualified.

### General lifting closure (PROVED)

The subsequent PRIME_POWER_PROOF.md proves the stronger all-order statement
s(mp^r)=(2/p)^m s(mp^(r-1)) modulo p^r, and identically for D, for every odd prime p and all m,r>=1. Its coefficient-valuation argument uses the free action of odd semiregular stabilizers on unordered pairs. Independent checks include even m, p dividing m, r=1, and p=3. Normalizing by (2/p)^(-mr) yields p-adic Cauchy sequences. Do not infer a higher-prime-power passing criterion from D(p) alone. This general proof supersedes the original narrow prime-square gap qualification; the original census and literature-scope cautions remain in force.
