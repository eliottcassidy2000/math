# Session 41: additive supports, two norm algebras, and both coloring branches

Status: **PROVED** for elementary identities and obstructions below; **FINITE-EXACT** only for the enumerated arithmetic controls. No new Goldbach or sum-of-two-squares classification is claimed. This package is independent of edits to the frozen `session41-fibonacci` package.

## 1. Inheritance and the distinction that must survive

The inherited sources are `duck_zeckendorf_20260925.md` (Z3–Z14), `zeckendorf_guard_automaton_20261003.md` (E1–E4 and arithmetic registers), and `additive_basis_goldbach_zeckendorf_s501.out`, all under `05-knowledge/results/`. The first already proves the additive charge and its binary carry. The second explicitly corrects a conflation of that charge with the owner's stripe. The third is a finite arithmetic audit, not an infinite coverage proof.

There are three different objects:

* The historical 35-letter owner stripe, translated R→C and K→A, whose letter 35 is B.
* The canonical recursive candidate from Session 41, whose letter 35 is C. Its letter at n is C when the smallest standard Zeckendorf index is 2, A when that index is odd, and B when it is even and at least 4.
* The four-state charge Q(n) mod 2. This is neither three-letter stripe and has a genuine zero state.

Both stripes agree through 34. Nothing here chooses which infinite continuation the owner intends. In particular, the earlier row-35 choice remains legal as a typed representation. With a fixed canonical block dictionary the compressed-zero D grammar fails at row36; the retained-zero R grammar rescues row36 but fails at row37. These are scoped grammar obstructions, not rejection of the supplied row35. The natural encompassing space is the set of typed representations and row addresses, before any single row-selection rule is imposed.

A concrete additive consequence keeps both ends visible. The pair 1+35=36 is in color bin CC for the canonical candidate and CB for the supplied historical stripe. Its sum, prime status, norm status and canonical Fibonacci carry are unchanged. The older stripe must not be silently extended beyond the supplied data. Our finite prime-pair tables explicitly use only the canonical candidate.

### Incoming complete branch classification and the selection-rule guard

The audited sibling `session41-fibonacci/BRANCH_ADDENDUM.md`, Theorems4–5, now identifies the larger occurrence-dependent universe exactly. Its antecedent C→B repair is credited there to incoming `odd_cell_gluing_20261023/ZERO-SLOT.md`, commit0429c7d8fdf. With independent recursive choices per occurrence, every retained-zero R stripe is obtained by the coordinatewise substitutions

    canonical A -> A; canonical B -> A or B; canonical C -> B or C.

The compressed D grammar permits only the B→A freedom. Thus the supplied old35B and canonical35C are both globally realizable in R; the fixed dictionary is a stricter synchronization requirement. No infinite D stripe in this universe can have B at35.

A direct consequence is important for additive number theory: for ANY subset E of the canonical C positions there is an R stripe whose C positions are exactly E. Keep C precisely on E and change the other canonical C positions to B, leaving canonical B positions unchanged. This follows immediately from the full Cartesian language classification.

In particular, choosing E to be the prime canonical-C positions makes “every C is prime” true by construction, not by discovering a theorem about primes. Choosing E by Gaussian-norm membership is equally possible. Therefore arithmetic correlations have mathematical content only after the selector rule is fixed independently of the property being tested, or when a theorem quantifies over the entire selector family. The shared-support and quadratic-algebra results below do not rely on choosing either stripe as the owner's sole intention.

## 2. Shared digits are the actual additive attachment ports

Write f_i=F_(i+2), i≥0, and let Z(n) be the binary no-adjacent-one support on f_i. Put α=(3−√5)/2 and

    b(n)=floor((n+1)α),   Q(n)=(n−2b(n), b(n)).

The inherited theorem says that Q(n) is the sum of v_i over ANY distinct Fibonacci representation, where v_0=(1,0), v_1=(0,1), v_(i+2)=v_i+v_(i+1). For ordinary addition,

    Q(a)+Q(b)−Q(a+b)=(-2δ,δ),
    δ=b(a)+b(b)−b(a+b) ∈ {-1,0,1}.

Disjoint supports imply δ=0, because their union is a distinct representation. The converse fails: 3+3 has completely overlapping support {f_2}, but b(3)+b(3)−b(6)=1+1−2=0. Thus even the exact zero carry does not certify absence of shared atoms. The raw support or an overlap sidecar is genuinely stronger.

### PROVED: k-summand carry range and tree independence

For k≥2 positive integer summands n_1,...,n_k, put N=Σn_i and D=Σb(n_i)−b(N). If u_i={(n_i+1)α}, then

    D=−floor(Σu_i−(k−1)α),
    ceil((k−1)α)−k ≤ D ≤ ceil((k−1)α).

Proof: expand each floor, noting that N+1=Σ(n_i+1)−(k−1); then use 0<Σu_i<k and the irrationality of (k−1)α. Also

    ΣQ(n_i)−Q(N)=(-2D,D).

Every binary addition tree produces exactly this same total correction: its internal b-values cancel. The range is sharp when the summands vary independently, since fractional parts of positive multiples of an irrational number are dense and sums of k independent such phases are dense in (0,k). Only the bound and tree independence, not sharpness, are used by the checker. For k=1, D=0 separately.

This is a useful all-arity repair, but it is still a charge statement rather than a primality statement.

### PROVED: an exact multivariate shared-support transfer matrix

For k summands, retain at digit position i the subset S_i⊆{1,...,k} of summands using f_i. The exact legal transition is

    S_i ∩ S_(i+1) = empty.

Assign the state weight ∏_(j∈S_i) x_j^(f_i). Start with the empty previous state; multiply these weighted transition matrices over the desired digit positions and sum all terminal states. There are 2^k states. The coefficient of x_1^n_1...x_k^n_k is exactly one if every n_j has a canonical representation in the retained positions, and zero otherwise. This follows by the bijection between a state path and k binary canonical digit strings, followed by uniqueness of each Zeckendorf representation.

Restricting the state alphabet to empty and singleton sets gives k+1 states and counts precisely the tuples with pairwise disjoint supports. For k=2 its unweighted adjacency matrix in order empty,{1},{2} is

    [1 1 1]
    [1 0 1]
    [1 1 0].

This is an exact generating-function compiler for the shared-digit constraint. Multiplying separate marginal counts without the shared subset state would lose the restriction. It is the same logical reason that separately legal tournament-cell attachments do not ensure compatibility at a shared port. It is not an identification between these two mathematical object classes.

The unrestricted matrix's extra state {1,2} retains overlap. A scalar carry alone cannot replace it, as the 3+3 example proves. Owner-color conditions can be imposed by retaining the smallest occupied index class in each summand, or simply filtering the uniquely represented coefficient indices. Prime and Gaussian-norm filters are additional arithmetic predicates, not consequences of this digit transition.

## 3. PROVED: a larger universe of quadratic charge algebras

For any integers s,u, form R_(s,u)=Z[t]/(t²−s t−u). Represent A+B t by (A,B). Then

    (A,B)*(C,D)=(AC+uBD, AD+BC+sBD),
    conjugate(A+B t)=A+B(s−t),
    N_(s,u)(A,B)=A²+sAB−uB².

Multiplying the two conjugates, or multiplying the corresponding 2×2 multiplication matrices and taking determinants, proves

    N_(s,u)(xy)=N_(s,u)(x)N_(s,u)(y).

Two meaningful ends of the analogy now occupy different parameters:

* Fibonacci/golden charge: s=1,u=1, discriminant 5, norm A²+AB−B². The atom v_i represents t^i and has norm (-1)^i.
* Gaussian integers: s=0,u=−1, discriminant −4, norm A²+B². Here multiplication is exactly the classical sums-of-two-squares composition law.

The quadratic family contains both, but preserves their distinct multiplication and norm forms. There is no need to discard either, and no license to replace one norm by the other.

### PROVED: why the numerical charge is not multiplicative

The integer represented by Fibonacci charge (A,B) is e(A,B)=A+2B. For golden multiplication,

    e(x*y)=e(x)e(y)−BD.

Thus for x=Q(a),y=Q(b), the exact numerical defect is b(a)b(b), not zero. The smallest positive hostile is a=b=2: Q(2)=(0,1), its golden square is (1,1), whereas Q(4)=(2,1). Their evaluated integers are 3 and 4 respectively.

More strongly, evaluation t→2 cannot define a unital ring homomorphism from the golden ring to ANY nonzero quotient Z/mZ: the defining relation evaluates to 4−2−1=1. This obstructs a multiplicative prime/square classification based on pretending Q(n) is a ring embedding of ordinary integers.

The additive charge DOES embed Fibonacci powers into a genuine ring. The wrong step is identifying the separate numerical evaluation A+2B with ring evaluation at a root of t²−t−1.

### PROVED: the mod-two obstruction is structural, not a numerical mismatch

Modulo 2, the golden algebra is F_2[t]/(t²+t+1)=F_4. Its three nonzero elements form a multiplicative cycle of order 3. There is no nonzero nilpotent.

Modulo 2, the Gaussian algebra is F_2[t]/((t+1)²). The nonzero element 1+t squares to zero. Hence the two four-element algebras are not isomorphic, and no relabeling of their colors can preserve both addition and multiplication. In particular, the golden order-three multiplication cannot be reinterpreted as Gaussian multiplication by i, whose mod-two square is 1.

The companion phase dictionary `session41-root/BINARY_PHASE_DICTIONARY.md` uses the positive side of the golden structure: multiplication by a nonzero F_4 element together with Frobenius gives the six permutations of the three nonzero states. In basis (1,t), Frobenius is (A,B)→(A+B,B), while multiplication by t and by t² are the two order-three companion generators. This is a transport/charge alphabet, not an identification with the original owner stripe; the Fibonacci incidence quotient loses an original letter distinction.

Their norm shadows also differ: the golden norm is 1 on all three nonzero F_4 states; the Gaussian norm is A+B modulo 2 and vanishes on the nonzero diagonal state (1,1). This makes the loss precise. An arbitrary bijection of three labels remains possible, but it does not preserve this arithmetic.

### PROVED: the positive other end, a prime-local norm bridge

The mod-two obstruction is NOT an obstruction at every prime. Modulo 3,

    t -> 2+i,       (A,B) -> (A+2B,B)

is an isomorphism from the golden algebra to the Gaussian algebra. Indeed
(2+i)²−(2+i)−1=0 modulo3, the coordinate map has determinant1, and

    (A+2B)²+B² = A²+AB−B² (mod3).

Consequently the map preserves multiplication, conjugation, and norm. Both sides are the nine-element field. The checker verifies all 81 products and all nine norms independently.

More generally, for an odd prime p≠5, each algebra is split F_p×F_p when its discriminant is a nonzero square, and is the quadratic field otherwise. Completing the square proves this classification: two distinct roots give evaluation into the product, while irreducibility gives the degree-two field. Hence the golden and Gaussian algebras are isomorphic over F_p exactly when 5 and −4 have the same square/nonsquare type. Field or product isomorphisms preserve conjugation and norm. At p=5 the golden algebra is nonreduced and the Gaussian algebra split, so they differ again.

This produces a local comparison program that genuinely considers both ends: identify the coefficient prime before transporting a quadratic invariant. The mod3 isomorphism does not fix the numerical evaluation defect e(x*y)−e(x)e(y)=−BD, which can still be nonzero modulo3. Thus it does not convert Q(n) into a multiplicative coloring of ordinary integers.

## 4. PROVED: the actual Fibonacci-to-two-squares bridge

For k≥1,

    F_(2k+1)=F_k²+F_(k+1)².

One elementary proof uses the Fibonacci addition identity
F_(m+n)=F_(m−1)F_n+F_m F_(n+1), itself proved by induction in n, with m=k+1,n=k. Thus every odd-index Fibonacci number has an explicit Gaussian norm representation. Consecutive Fibonacci numbers are coprime by repeated Euclidean subtraction.

Their Euclid triple is

    (F_(k+1)²−F_k², 2F_k F_(k+1), F_(2k+1)).

At k=1 this is the degenerate triple (0,2,2); nondegenerate triangle language applies only for k≥2. Its common content is 2 if k≡1 (mod3), and 1 otherwise. Indeed Fibonacci parity repeats 0,1,1; the two parameters are both odd exactly at k≡1 (mod3). Coprimality excludes any odd common divisor; in the odd-odd case the hypotenuse is 2 mod4, so the content is exactly 2. This is a legitimate parity/color connection to primitive Pythagorean normalization.

For the canonical owner stripe, each F_(2k+1), k≥1, has color A: its only Zeckendorf atom has an odd index. This supplies an infinite explicit sum-of-two-squares subsequence in one owner class. It does NOT characterize that class. For instance A-colored 54 is not a sum of two squares, while A-colored 2 is, despite both being 2 mod4.

CITED, not reproved here: the full Gaussian norm criterion is that every prime congruent to 3 modulo 4 has even valuation. The retained sidecar is therefore the appropriate factorization/valuation information, or an explicit Gaussian representation, not merely a three-letter color or even the value modulo4.

## 5. Arithmetic controls and the limits they certify

`check.py` uses exact integer arithmetic and explicit exceptions, with identical normal and optimized receipts. It checks:

* All 90,601 ordered pairs 0≤a,b≤300 for charge carries, disjoint-support implication, exact product defect and golden norm multiplication.
* All 31³ ordered triples through 30 for the ternary bound and tree-independent correction.
* All s,u∈{-2,-1,0,1,2}, and all pairs of coordinates in [-2,2]², for universal quadratic norm multiplication: 15,625 tests.
* Both four-element mod-two rings for square-zero elements, plus all 81 products under the explicit modulo-three isomorphism.
* Fibonacci norm and primitive-content identities for k=1,...,59.
* The full and disjoint two-summand transfer polynomials for nine digit positions, comparing every coefficient against independent greedy supports.
* Every even n from 4 through 2000 for a complete ordered prime-pair partition by canonical owner-color pair and exact carry δ. Selected full histograms are in `receipt.json`. No target is missing in this finite interval, which is not an infinite theorem.

The following pairs have the SAME canonical owner color AND same value mod4, but different arithmetic membership:

| Predicate | A | B | C |
|---|---|---|---|
| primality | 2 prime / 10 composite | 3 prime / 55 composite | 9 composite / 17 prime |
| two squares | 2 yes / 54 no | 8 yes / 24 no | 4 yes / 12 no |

The two-squares hostiles are also elementary: squares mod3 are 0 or1, so if a²+b² is divisible by3 both a,b are divisible by3 and the sum is divisible by9. This excludes 24 and12. For54, both squares would be divisible by9, reducing to a²+b²=6, impossible mod3 by the same argument. The positive cases are explicit.

For odd prime summands p>2, Q(p) mod2 has first coordinate 1. Only two of the four charge states occur there. This is an immediate parity consequence, not extra distributional evidence. An owner-stripe partition is another coordinate entirely.

## 6. Literature and next obligations

Primary sources checked 2026-10-23:

* [Oxford, Rings & Arithmetic 9: The Gaussian integers](https://people.maths.ox.ac.uk/flynn/genus2/alg0506/RandALect9.pdf), especially the norm definition and the sum-of-two-squares criterion. We use the criterion only as cited classical input.
* [Yohan Hosten, On the variation of the sum of digits in the Zeckendorf representation, v2](https://arxiv.org/abs/2309.14285v2). The abstract describes a Zeckendorf-adic model, an algorithm for digit-sum variation distributions, and mixing for block additions. This is a relevant route for digit-carry statistics. Its statistic is digit-sum variation, not the owner's stripe, the charge b(n), or prime-pair coverage. No theorem is transferred between these without a map.

OPEN next obligations:

1. The occurrence-dependent row-selection family is now classified in the cited branch addendum. The remaining stricter problem is a globally shared, possibly noncanonical block dictionary: one expansion per Fibonacci value reused everywhere. Do not cite the occurrence-dependent theorem as solving that synchronization problem.
2. Use the shared-support matrix with explicit residue/valuation sidecars to study restricted representation counts. A finite automaton for fixed moduli cannot by itself replace unbounded prime valuations or primality.
3. For a Goldbach target, proving every color bin nonempty is unnecessary and often false; proving at least one arithmetic-admissible bin nonempty uniformly is the real obligation. The tables do not address that uniform lower bound.

The useful common structure is not repeated small numbers: it is a faithful lifted state, an explicit quotient, and the compatibility or arithmetic data the quotient loses. Here the missing data are overlap, carry, boundary selector, and quadratic multiplication law. Their exact retention gives generating functions and obstructions; forgetting them creates attractive but false unifications.
