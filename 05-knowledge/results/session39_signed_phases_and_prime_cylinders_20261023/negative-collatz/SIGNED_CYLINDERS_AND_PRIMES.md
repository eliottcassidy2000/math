# Signed Collatz cylinders, primes, and the negative sheet

Session 39, 23 October 2026. Status: proved elementary identities and corollaries of established prime-distribution theorems; independently executable finite controls. No Collatz proof or historical novelty claim.

## 1. Maps and precise relation

For ε in {+1,−1}, define Cε(n)=n/2 for even n and 3n+ε for odd n. On positive odd integers define Uε(n)=(3n+ε)/2^v2(3n+ε). Both positive domains are invariant. Negation J(n)=−n satisfies C+(J(n))=J(C−(n)), and likewise for odd acceleration. Thus positive 3x−1 is exactly the negative-integer sheet of 3x+1, not positive 3x+1 with approximately the same dynamics.

The known positive U− cycles are:

- 1 → 1, valuation word (1).
- 5 → 7 → 5, valuation word (1,2).
- 17 → 25 → 37 → 55 → 41 → 61 → 91 → 17, valuation word (1,1,1,2,1,1,4).

These correspond to negative cycles under U+. Their existence disproves any claim that every positive 3x−1 orbit reaches 1. The conjecture that every positive orbit reaches one of these three cycles is unresolved; it is not implied by the ordinary positive Collatz conjecture. Angeltveit, February 2026, §7 explicitly states the integer extension conjecture and implements the positive 3x−1 shortcut map: https://arxiv.org/html/2602.10466v1#S7 . Distinguish the shortcut map (3n−1)/2 on odd n from U−, which removes all powers of two.

## 2. Signed affine monodromy and exact cycles (classical/inherited)

Let w=(a0,...,a(k−1)) with every ai≥1. Put Sj=sum_(i<j) ai, A=Sk,

B(w)=sum_(j=0)^(k−1) 3^(k−1−j) 2^Sj, and D(w)=2^A−3^k.

Composing the indicated rational affine branches gives

Fε,w(x)=(3^k x+εB)/2^A.

This follows by induction: appending a letter a updates (A,B) to (A+a,3B+2^A). D is nonzero and odd, and B is odd. The unique fixed point is x=εB/D. Consequently an actual positive integral cycle with the supplied word exists exactly when D divides B and εD>0. Repetitions of a shorter cycle are allowed in this statement.

There is no extra unverified admissibility assumption in this closed-word criterion. Section 3 proves that a positive odd n realizes w exactly when 3^k n+εB≡2^A mod 2^(A+1). If n=εB/D is integral, n is odd, and the fixed-point equation gives that congruence. Alternatively, cyclically indexed carries obey 3B_j+D=2^a_j B_(j+1); divisibility propagates because D is odd, and all resulting orbit points are odd. Thus the indicated valuations are exact.

For U− the three certificates (k,A,B,D) are (1,1,1,−1), (2,3,5,−1), and (7,11,2363,−139), with 2363=17·139. Positive U+ cycles require A/k>log_2(3); positive U− cycles require A/k<log_2(3). This is an affine fixed-point condition, NOT the assertion that the full affine monodromy is the identity. It must not be conflated with flat holonomy in a frame representation.

For an actual U− cycle with minimum m and maximum M, telescoping the product of branch ratios gives

δ=k log 3−A log 2=−sum_i log(1−1/(3xi)),

k/(3M) < δ < k/(3m−1).

The strict inequalities follow from t<−log(1−t)<t/(1−t) for 0<t<1. These are necessary cycle bounds, not exclusion of all further cycles.

## 3. Exact finite cylinders and mirror symmetry

**Cylinder theorem.** For each nonempty positive valuation word w and each sign ε, its positive odd starting values are exactly one residue class modulo M=2^(A+1):

rε(w) ≡ 3^(−k)(2^A−εB(w)) (mod M).

The representative rε in [1,M) is odd. Moreover r−≡−r+ (mod M).

Proof: if the actual first k valuations are w, the odd final endpoint gives 3^k n+εB=2^A·odd, which is precisely the displayed congruence. For sufficiency, let Nj=3^j n+εB(prefix_j). The total congruence and the decomposition B(w)=3^(k−j)B(prefix_j)+2^Sj B(suffix_j) show that 2^Sj divides Nj, so all candidate intermediate xj=Nj/2^Sj are integers. Work backward from xk odd: 3xj+ε=2^aj x(j+1) with aj≥1 forces xj odd, hence the exact valuation is aj. Starting positive, each branch stays positive, so these really are Uε iterates. Invertibility of 3 modulo M proves uniqueness. B odd proves rε odd. Finally −2^A≡2^A modulo M, proving reflection. ∎

The last bit is essential: modulus 2^A ensures divisibility but does not ensure an odd endpoint. Replacing the residue by rε+2^A changes the last oddness bit and cannot preserve the full valuation word.

Examples: w=(1,2) has r−=5 and r+=11 modulo 16. The 17-cycle word has r−=17 and r+=4079 modulo 4096. These are finite shadows, not a positive-sign cycle through 4079.

## 4. A valid prime connection, and an equally exact limitation

**Prime universality corollary.** Every finite positive valuation word occurs for infinitely many prime starting values under each sign. Its prime starting values have asymptotic proportion 2^(−A) among all primes.

Proof: rε is coprime to 2^(A+1), so Dirichlet's arithmetic progression theorem gives infinitude. The prime number theorem in arithmetic progressions gives

#{p≤X: p prime and Uε has prefix w at p} ~ 2^(−A) X/log X.

Source for both external theorems: Kedlaya, *Notes on analytic number theory*, Chapter 4, Theorems 4.11 and 4.12, https://kskedlaya.org/ant/chap-primes-in-ap.html . This is a fixed-word/fixed-modulus limit; no uniform error bound as A grows with X is being claimed. The elementary cylinder proof is ours here; the prime-distribution theorems are established inputs, not newly proved.

**Finite-valuation primality obstruction.** Every such cylinder also contains infinitely many composites. If its representative r>1, the integers r(1+Mt), t≥1, suffice. If r=1, use (1+Mt)^2. Hence no classifier depending only on a fixed finite valuation prefix can distinguish prime from composite inputs. More strongly, no always-terminating decision procedure that sees only successive valuations and stops after a finite, possibly input-dependent, number can correctly decide primality for all positive odd inputs: any terminal decision applies to a cylinder containing both kinds. This statement does NOT prohibit primality algorithms that also inspect the input value, orbit magnitudes, modular factors, or additional data.

For example 4079 is prime, 16711663=4079·4097 is composite, and both have the U+ prefix (1,1,1,2,1,1,4). Likewise 17 and 69649=17·4097 share that U− prefix. The words capture neither prime factorization nor primality.

Therefore “Collatz represents every prime” can be made true as a finite-symbolic universality statement, but that is weaker than an algebraic characterization of primes, and different from representing every prime in a single convergent basin. Denominator-prime divisibility D|B is another valid arithmetic bridge, with a different quantifier: all prime-power tests must hold for one common word.

## 5. Repository inheritance and boundaries

Read-only inspection used snapshot e5d18e55 in session35-publication/clean-repo. The closest predecessors inspected were:

- `05-knowledge/results/collatz_local_prime_power_gates_20261020/LOCAL_CARRY_AND_PRIME_POWER_GATES.md`: exact binary-word integral-cycle recognition, chronological rotation content, two-one exclusion, prime-power and common-witness cautions. Its denominator 1805=5·19² and carry 475 illustrate why prime support alone is insufficient.
- `05-knowledge/results/collatz_one_swap_negative_exclusion_20261020/NEGATIVE_ONE_SWAP_EXCLUSION.md`: claimed audited exclusion/classification within a specific one-swap mechanical family, not arbitrary negative cycles. Not re-proved or invoked as a premise here.
- `05-knowledge/results/collatz_minus17_return_20261003.md`: established negative 17 cycle, carry 2363, and precise positive U+ return cylinders with a last-oddness-bit guard. The guard and affine algebra are inherited, not new.
- `agents/broadcast/MSG-3231-from-mac-mini-2026-10-20-mac-mini-2026-10-20-debtcones.md`: debt-cone/153N annulus work explicitly leaves first-exit payment and all-scale annulus control open; threshold-two cone universality is explicitly false. No premise about those unproved bounds is used here.

Targeted title and content searches were broad enough to find these close predecessors, but not a proof of historical or repository priority. In particular this report is a clean signed/prime synthesis, not a priority claim for the well-known affine cycle formula or residue-cylinder method.

Open obligations: prove global convergence on either positive sheet; classify all additional integer cycles; obtain effective useful least-prime bounds for growing prescribed words; and construct any proposed correspondence to tournament/symmetry objects with explicit maps and invariants. None follows from finite cylinders, local valuation tests, common formal vocabulary, or finite simulation.

## 6. Reproduction

Run `python3 verify_signed_cylinders.py` and `python3 -O verify_signed_cylinders.py`. Both outputs match the supplied receipts. No mathematical repository code is imported; checks use explicit exceptions, so optimization does not disable them.

Universe: all words of lengths 1 through 5 over {1,2,3,4}, both signs, four explicit representatives per cylinder: 10,912 endpoint/valuation checks, 2,728 hostile last-bit checks, and 14 positive integral-cycle word checks (including repetitions). Independently enumerate every odd residue for all words of lengths 1 through 3 over {1,2,3}, both signs: 78 complete residue checks. The three known cycle words and six prime/composite example pairs are also replayed exactly. These are finite validation controls, not evidence of universal convergence.
