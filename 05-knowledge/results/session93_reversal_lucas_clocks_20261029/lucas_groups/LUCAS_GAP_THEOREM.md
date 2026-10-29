# The strict reflection gap: Lucas closure beyond the eleven-point exception

**2026-10-29; Session93.** Status: **PROVED relative to the cited classification** for the universal group theorem; **PROVED elementary** for projective exclusion, Fibonacci gap and parity; **FINITE-EXACT** for explicit Jordan certificates. No novelty priority claim, no assertion of infinitely many Lucas primes.

Source pin: `59cb22d0e93b6c364c4ec5d2954d9fe14ed468d8`. Inheritance: `05-knowledge/results/lucas_mirror_gap_spaces_20261028/LUCAS_MIRRORS.md`, especially §§4–6. That packet proves the degree11 closure PSL(2,11), its twelve matching frames and a degree29 Jordan certificate. The all-odd-period reader is inherited from `01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`, with the all-period repair in `05-knowledge/results/session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md`.

Concept board: fixed-point gaps; projective descent; signed reflection; guarded-versus-unguarded closure; short group certificates. Anchor: universal Lucas closure. Niche: projective fixed-subgeometry obstruction. Wildcard: very short Jordan witnesses in degrees199 and521. Hostile controls: dihedral fixed-point endpoint, PSL(2,11), outer field involutions.

## 1. Precise group theorem

Let p be an odd prime and G≤S_p a transitive permutation group. Suppose G contains an involution r with exactly f fixed points and

> **1 < f < sqrt(p).**

Then either G contains A_p, or p=11 and G is PSL(2,11) or M11 in its eleven-point action. In particular, for p≠11 the only possibilities are A_p and S_p.

Transitivity suffices: it implies primitivity at prime degree; Cauchy's theorem also supplies a p-cycle. Alternatively a specified full cycle supplies both facts. The statement with **1≤f** is false for every odd prime: shift and negation generate D_p, with negation fixing exactly one point. The strict upper endpoint is automatic whenever f²<p and is the useful exact-arithmetic test.

### Classification input [CITED, CFSG-dependent]

Jones, *Primitive permutation groups containing a cycle*, [Theorem1.2, case k=0](https://arxiv.org/pdf/1209.5169), lists the non-alternating-containing primitive full-cycle groups as affine groups between C_p and AGL1(p), projective groups between PGL_d(q) and PΓL_d(q) of degree (q^d−1)/(q−1), and the exceptional actions PSL2(11), M11, M23 of degrees11,11,23. This is a classification input, not Jordan's elementary small-cycle criterion. Jones explicitly notes its dependence on CFSG.

Affine involutions are x↦−x+b, hence fix one point. The projective alternatives are excluded below. The [ATLAS eleven-point M11 action](https://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/M11G1-p11B0) has one involution class, fixing3 points. The [ATLAS twenty-three-point M23 action](https://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/M23G1-p23B0) has one involution class, fixing7 points; 49>23 excludes it. PSL2(11) remains a genuine exception, realized by the inherited Lucas-five reflection fixing3 points.

### Projective involution exclusion [PROVED]

Write S_d(q)=1+q+...+q^(d−1). If S_d(q) is prime then d is prime: a proper divisor e of d gives a proper factor S_e(q). Thus d=2 or d is odd.

**Linear case, d odd.** Represent a projective involution by A with A²=λI. In odd characteristic, λ must be a square: otherwise the irreducible polynomial X²−λ makes the vector-space dimension even. After scaling, A²=I; one of its two eigenspaces has dimension at least (d+1)/2. In characteristic two λ has a square root; after scaling A=I+N with N²=0, so ker N has dimension at least (d+1)/2. In either case there are at least S_((d+1)/2)(q) fixed projective points. Its square exceeds S_d(q), because the polynomial square has every coefficient at least1, and some strictly greater than1. Thus f²>p.

**Linear case, d=2.** Primality of q+1 forces q even. A nontrivial projective involution is unipotent and fixes exactly one point. This fails the lower strict inequality.

**Outer semilinear case.** The image of an involution in the field automorphism group has order2, so q=s². Lift it to a σ-semilinear J, σ:x↦x^s, with J²=λI. Associativity gives σ(λ)=λ, so λ∈F_s*. Surjectivity of the norm F_q*→F_s* permits rescaling J until J²=I.

For completeness, descent needs no unverified fixed-subgeometry assertion. Choose θ∈F_q\F_s. For every v set

    w=(Jv−v)/(σθ−θ),    u=(σθ v−θ Jv)/(σθ−θ).

Then Ju=u, Jw=w and v=u+θw. Hence fixed vectors span V over F_q; select a fixed F_q-basis. Relative to this basis J acts coordinatewise by σ. A projective point fixed by J has Jv=cv with cσ(c)=1; finite-field multiplicative Hilbert90 (the norm-one subgroup is the image of t↦t/σt) rescales v to a fixed vector. Therefore exactly S_d(s) projective points are fixed. Since S_d(s)²>S_d(s²), this again lies above sqrt(p). The dual hyperplane action has the same argument.

This closes the projective and semilinear loopholes; no projective exceptional prime degrees remain in the strict gap.

## 2. All odd Lucas periods with prime state count

Let n=2m+1≥5. Let W_n be the named cyclic binary words with no cyclic adjacent ones, including zero. Let r reflect positions i↦−i. The count |W_n|=L_n is classical; the inherited odd-period golden reader identifies it bijectively with O/(φ^n−1), O=Z[φ]. If p=L_n is prime, that quotient ring is F_p. Define t as arithmetic translation x↦x+1 through this reader and set G_n=⟨t,r⟩.

The reflection has exactly

    f_n=F_(m+2)=F_((n+3)/2)

fixed words. Indeed a fixed word is specified by m+1 reflection-orbit bits; the last orbit comprises the adjacent positions m,m+1 and must be zero. The remaining m bits form an ordinary path independent-set problem, giving F_(m+2).

The useful exact identity is

    L_(2m+1) − F_(m+2)^2 = 2 F_m F_(m−1).

For m≥2 the right side is positive and f_n>1. It follows directly from F_(2m)=F_m(2F_(m+1)−F_m) and F_(2m+2)=F_(m+1)(2F_m+F_(m+1)). Thus every such Lucas reflection falls strictly in the gap, not merely asymptotically.

**Conclusion.** For odd n≥5 with L_n prime:

- n=5: G_n=PSL(2,11), inherited exact exceptional case.
- n≠5: G_n contains A_(L_n).
- More precisely G_n=S_(L_n) if n≡11 mod12, and G_n=A_(L_n) otherwise.

The last statement is a sign computation, not another classification assumption. Translation is an odd-length cycle, so is even. Reflection has (L_n−F_((n+3)/2))/2 transpositions. Lucas numbers and Fibonacci numbers modulo4 both have period6. For n≡1,5,7,11 mod12 the resulting sign is respectively +,+,+,−. Multiples of3 have even L_n and cannot occur here with L_n prime, since n≥5. The exceptional n5 reflection is even but its closure is properly smaller than A11.

No assertion is made about how often L_n is prime, and no extrapolation is made to composite ring size. Primality is a hypothesis, not a prime-distribution conclusion.

## 3. Short independent finite certificates

`check_lucas_gap.py` uses only the standard library and explicit exceptions, hence checks remain active under `python -O`. It performs **2340 checks**, including direct word enumeration at odd periods5..19, the universal identities checked through1001 as regression tests (not proofs), projective integer bounds, and three finite Jordan certificates. Normal and optimized output coincide.

Use labels ell(w)=sum w_i α^i modp, T:x↦x+1, R the reflected word. For an exponent list (a1,...,ak), define g0=id and gj=R T^aj g_(j−1). The certificates are:

| n | p | α | exponents | cycle lengths of g | power | isolated cycle |
|---|---:|---:|---|---|---:|---:|
|7|29|24|7,1,9,26,25,16,27,7,18,5|2,13,14|14|13|
|11|199|62|161|4,14,23,41,51,66|644028|23|
|13|521|422|168,442,105|2,2,3,7,15,209,283|1774410|7|

Each isolated prime cycle fixes at least3 points. Jordan's theorem therefore forces alternating containment, and generator signs give **A29, S199, A521**. These finite conclusions do not require the CFSG classification, only Jordan's classical theorem. In degree199 a **single translated reflection** supplies the certificate after powering. The search uses a fixed seed, but the witness is explicit and checked independently of any group-order package.

## 4. What this says about retained tournament information

Source: guarded Lucas words; target: their arithmetic translation/reflection closure; map: odd golden reader. Preserved: individual state identity, named reflection, exact transition permutations. Lost on passing to the unrestricted group: the original guarded toggle relation and any chosen Paley matching frame. Minimal restored data is observer-dependent, not determined by permutation-group order.

The inherited period5 example restores a Paley decoder by twelve matching frames, with only two surviving all fifteen native arcs. The present theorem shows that this exceptional compact incidence closure does not persist at other prime Lucas state counts. Indeed A_p (p≥4) is transitive on ordered distinct point pairs; an invariant directed binary relation is therefore only empty or complete. Consequently **no tournament orientation on those p states can be invariant under G_n** once n≠5. This elementary obstruction is stronger and cleaner than treating the large group as tournament automorphisms.

This is an exact bridge to reversal marking: adding observational operations can enlarge the symmetry of the unmarked carrier while destroying the marked directed relation. It is not an identification of Lucas reflection with tournament converse, and it does not classify black/blue metagraph edges. A tournament decoder must carry symmetry-breaking side data. The precise smallest sidecar at general prime Lucas periods is **OPEN** here.

Reproduce:

    python3 check_lucas_gap.py
    python3 -O check_lucas_gap.py

Source provenance and exact finite receipts remain separate from the universal proof. The theorem is a deduction from classical classification plus elementary bounds, not a proposed solution to a major open problem.

## 5. New decoder lower bound: a large group forces many marking states

**PROVED, elementary.** Suppose a finite group G acts faithfully on a named point set and C is a nonempty finite G-set of decoder states, with a G-equivariant map D:C→{tournaments on that point set}. For each c, its stabilizer G_c lies in Aut(D(c)). Every tournament automorphism group has odd order: otherwise Cauchy's theorem supplies an involution, which swaps some pair of vertices and reverses its arc, impossible. Hence G_c is odd and

    2^v2(|G|) divides |G.c|.

In particular any such closed equivariant carrier has at least 2^v2(|G|) states. This lower bound concerns a **finite carrier closed under the entire unguarded group**. It is not a bound on a one-off certificate, a single marked decoder, or a carrier closed only under guarded native operations.

For G=A_p this gives at least 2^(p−s2(p)−1) states; for G=S_p at least 2^(p−s2(p)), where s2 is binary digit sum and Legendre's factorial valuation gives v2(p!)=p−s2(p). Thus the explicit closures force:

| State count p | Closure | Decoder-state lower bound | Fixed-length binary storage lower bound |
|---:|---|---|---:|
|29|A29|2^24|24 bits|
|199|S199|2^194|194 bits|
|521|A521|2^517|517 bits|

At the exceptional degree11 closure PSL2(11), |G|=660 gives only a factor4 lower bound; the inherited twelve-frame carrier is consistent with it. This supplies a precise reason the compact eleven-point matching repair cannot remain comparably small after moving to the higher prime-Lucas closures. It does **not** prove the displayed bounds are attained.

The mechanism is useful beyond Lucas words: an even symmetry group can act on undecoded information, but every stabilizer of a fully decoded tournament must have odd order. This is a 2-primary cost of exact equivariant orientation, not a numerological connection between particular vertex counts.

### Sharp exceptional repair: twelve states are necessary and sufficient

At n=5 the minimum nonempty finite equivariant tournament-decoder carrier has **exactly twelve states**. Indeed G=PSL2(11) is simple of order660. Every decoder-state stabilizer is odd. The only odd divisor of660 larger than55 is165. A subgroup of order165 would have index4; its coset action would map the simple group G injectively into S4 (the action is nontrivial), impossible since660>24. Thus every decoder orbit has at least660/55=12 states.

This is attained by the **actual orbit of the labelled Paley11 tournament**, not merely a putative matching-frame orbit. The checker enumerates ⟨T,R⟩ as660 point permutations, applies all of them to the55 directed Paley arcs, and obtains exactly12 distinct tournaments with stabilizer55. This independently verifies the carrier represented by the twelve inherited matching frames. The lower bound uses the classical simplicity of PSL2(11), whereas the orbit and stabilizer counts are finite-exact. No sharper general minimum for A29, S199 or A521 is asserted.

## Concurrent provenance at publication

The incoming packet `lucas_prime_reflection_towers_20261028` at commit
`4368e3847072d9a5d88ddb7d80470f9861c97b93` supplies overlapping reflection/holonomy
results and a general port-moment recurrence. See the parent
[incoming comparison](../incoming/CONCURRENT_COMPARISON.md) for credited overlap,
distinct refinements and independently replayed finite controls. This note
retains its initial source pin and direct proof; no priority claim is implied.
