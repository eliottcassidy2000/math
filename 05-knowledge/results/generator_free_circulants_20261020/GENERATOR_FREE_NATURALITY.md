# A generator-free prime-circulant map and sharp naturality obstructions

2026-10-20. **PROVED elementary statements; finite controls supplied separately.** The first-order nondefinability statement for the particular parity-projection map additionally uses the cited classical circuit lower bound. No historical novelty claim is made.

## 1. Inheritance, definitions, and exact questions

The inherited all-prime construction is `circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md`, with bounded Fitting routing in `bounded_fitting_routing_20261020/BOUNDED_FITTING_ROUTING.md`. The basis-free Artin–Schreier analogue is `artin_schreier_routing_20261020/ARTIN_SCHREIER_AND_CRITICAL_PERIODS.md`. HYP-9340 is used only in its repaired form: affine-equivariant labelled **surjection**, type **bijection**, and no blanket odd-symmetry preservation claim.

Let C be an abstract cyclic multiplicative group of order 2m, with unique involution z. Let X(C) be the skew subsets S: exactly one of x,zx belongs to S. Let Y(C) consist of the z-invariant subsets A whose indicator, read on C/<z> in any cyclic coordinate, has even-weight primitive root. The predicate defining Y is independent of the generator: the stabilizer of A determines its primitive quotient, and its primitive weight is the number of selected stabilizer orbits in C/<z>.

Write Xbar=X(C)/C and Ybar=Y(C)/C for multiplication-orbit sets. Aut(C) acts on both. At 2m=p−1, p prime, the inherited prime Sylow argument identifies these with the actual isomorphism types of circulant tournaments and sign-even circulant graphs on F_p. Sign-even means that every automorphism reverses an even number of graph edges relative to a reference orientation; it does not mean even vertex degree.

There are four different demands:

1. A uniform terminating algorithm on arithmetically labelled prime fields.
2. A construction natural under all automorphisms of the **abstract multiplicative cyclic group**.
3. A construction natural only under **field** isomorphisms, with the tournament predicate supplied as input.
4. One fixed first-order formula in the ring language expanded by a unary input predicate S, defining output membership A(x).

These demands are not interchangeable. The inherited least-primitive-root algorithm already supplies (1). Over F_p, field automorphisms are trivial, whereas Aut(C_(p−1)) is generally not. Naturality is necessary for first-order definability, but not sufficient. Indeed the arithmetic enumeration 0,1,1+1,... is canonical in a prime field, so using its least primitive root also defines a field-isomorphism-natural algorithm. The phrase “independent of an externally supplied generator” by itself does not impose (2): that stronger demand forgets addition and retains only the multiplicative cyclic group.

The active concept board is: antipodal-complement quotient; augmentation projection; dyadic obstruction; semilinear Frobenius; and bounded-depth parity. The hostile examples are m=4 abstractly, p=13 in the prime family, and q=9 over extension fields.

## 2. Positive theorem: canonical parity projection when m is odd

**Theorem A.** If m is odd, there is an Aut(C)-natural, C-equivariant surjection F:X(C)→Y(C), with every labelled fiber exactly {S,zS}. It induces an Aut(C)-equivariant bijection Xbar→Ybar. It needs no generator, ordering, or local predecessor choice.

**Construction.** Let H be the unique subgroup of order m. Since m is odd, C=H×<z>. Set s(h)=1_S(h) and

    t(S) = sum_(h in H) s(h) in F_2,
    w(h) = s(h)+t(S).

Define A=F(S) by 1_A(h)=1_A(zh)=w(h).

**Proof.** Because |H| is odd, sum_H w=t+mt=0. If w has primitive period d, its repetition count m/d is odd, so its total weight and primitive-root weight have the same parity. Thus A∈Y(C).

Conversely, for any even-weight w on H, the two words s=w and s=w+1 give this output. They correspond precisely to S and zS. Translation by h∈H reindexes the sum t and w. Multiplication by z complements s and changes t to t+1, leaving w unchanged. Hence F is C-equivariant. Its fibers are a single two-element subset of one C-orbit, so it induces a bijection of orbit sets. Every automorphism of C preserves H and z and reindexes the parity sum; therefore F is Aut(C)-natural. ∎

**Prime-field formula.** For p≡3 mod4, H is the set of nonzero squares. If r=|S∩H|, then

    F(S) = {±h : h∈S∩H},             if r is even;
    F(S) = {±h : h∈H\S},             if r is odd.

The image graphs are exactly the sign-even circulants. Equivalently they are the symmetric connection sets with |A| divisible by four. The map is affine-equivariant and exactly two-to-one on translation-invariant labelled objects, bijective on actual graph/tournament types. Forming squares, scanning the input, and taking parity are finite arithmetic operations; there is no primitive-root search.

**Mechanism and connection.** On F_2^H let J(s)=(sum_H s)·1. Oddness gives J²=J, so P=I+J is the canonical projection onto the augmentation subspace ker(sum_H), with kernel <1>. Here J is the averaging/Reynolds projector onto the constants for the odd group H, and P=I+J is its complementary augmentation projector. The latter replaces the derivative I+R by a generator-free projection: quotienting out complementation and imposing even parity are canonically equivalent when |H| is odd. The resulting pairing need not equal the earlier derivative pairing.

## 3. Negative theorem: complete abstract-cyclic naturality dichotomy

**Theorem B.** If m is even and m≥4, no Aut(C)-equivariant injection Xbar→Ybar exists. In particular no type bijection and no labelled natural map inducing a type bijection exist. Combined with Theorem A and the one-type case m=2, an Aut(C)-natural type bijection exists exactly when m is odd or m=2.

**Proof.** Temporarily use an exponent coordinate C=Z/(2m). The unit a=m+1 defines an automorphism α(i)=ai because gcd(m+1,2m)=1. It is the identity modulo m. Every symmetric graph word is m-periodic, so α fixes **every labelled A∈Y(C)**, and hence all of Ybar.

Consider the antiperiodic tournament word f=0^m1^m. Since m is even,

    f_(ai) = f_i + (i mod 2).

Its decimated word has alternating bits within each of its two halves and matching bits at both half-boundaries. The cyclic transition count of f is 2; that of αf is 2m−2. The transition count is invariant under rotations. As m≥4, the two counts differ, so [αf]≠[f]. Equivariance would give F([αf])=αF([f])=F([f]), contradicting injectivity. ∎

For m=2, both orbit sets have one member; the target is the zero word. This is the exception p=5. Thus, among odd primes, the abstract-multiplicative-group demand has the exact answer:

- p≡3 mod4: the canonical parity projection works;
- p=5: the unique type maps to the empty graph;
- p≡1 mod4, p>5: no such natural type bijection exists.

The first prime obstruction is p=13, m=6, with transition counts 2 and 10. For p=17 the counts are 2 and 14. The obstructing automorphism is x↦x^(m+1), not generally an automorphism of the field. This theorem therefore **does not** rule out using field addition or a parameter-free ring formula over prime fields.

A stronger check, not needed for Theorem B, is that when 4 divides m the automorphism α has no fixed tournament type. To see this, a fixed type would have f_(i+k)=f_i+(i mod2). If k is odd, iterating m times gives f_(i+m)=f_i+(m/2 mod2)=f_i, a contradiction. If k is even, write r=v_2(m), s=v_2(k), with k=0 handled by the same odd-orbit argument as s>r. For s≤r the subgroup generated by k contains m; iterating the relation from any even i then contradicts antiperiodicity. For s>r, the order of k in Z/(2m) is odd; iterating from any odd i contradicts periodicity. Thus no k works.

## 4. Genuine field-automorphism obstruction over even-degree extensions

For K=F_(p^r), q=p^r odd, form Xbar and Ybar using C=K^* and m=(q−1)/2. These are **multiplier types**, and Y is the even-primitive-root necklace target. They are not automatically full graph-isomorphism types or fully sign-even graphs.

**Theorem C.** If r is even, Frobenius σ:x↦x^p has no fixed point on Xbar, whereas Ybar has the fixed empty type. Consequently there is no Galois-equivariant bijection Xbar→Ybar. In particular there is no one parameter-free ring-formula output construction on all odd finite fields that induces these multiplier-type bijections.

**Proof.** If the multiplier type of S were fixed, there would exist a∈K^* such that a·σ(S)=S. The semilinear map L(x)=a x^p would then be an automorphism of the Cayley tournament on the additive group of K with connection set S. Its image under the quotient homomorphism ΓL_1(K)→Gal(K/F_p) is σ of order r. Therefore the finite order of L is divisible by r, hence even. A nontrivial involutive power of L is impossible as a tournament automorphism: any transposed pair of vertices has its oriented edge reversed. Thus the supposed fixed type cannot exist. The empty graph type is fixed by every field automorphism. A Galois-equivariant bijection must preserve fixed-point counts, a contradiction. ∎

A formula φ(x;S) in the ring language with unary predicate S respects isomorphisms of the expanded structures: F(σS)=σF(S). Thus, if it also induces the asserted multiplier-type bijection, it would contradict Theorem C. The same obstruction applies to formulas with parameters in F_p; arbitrary extension-field parameters may break Frobenius and are outside the statement. In fact q=9 alone already rules out the parameter-free construction on that one field in this precise model.

Exact fixed-type controls are:

| q | Total source/target types | Frobenius-fixed source/target types |
|---|---:|---:|
| 9 | 2 / 2 | 0 / 2 |
| 25 | 172 / 172 | 0 / 46 |
| 27 | 316 / 316 | 16 / 16 |
| 49 | 349536 / 349536 | 0 / 3372 |
| 81 | 13743895360 / 13743895360 | 0 / 2522 |

The q=27 control shows that the even-degree hypothesis is material. The theorem makes no claim of nonexistence for every odd extension degree. The positive map from Theorem A is field-isomorphism natural whenever q≡3 mod4, still within the multiplier/necklace model.

### Why the model boundary is essential

For q=p^r, r>1, the translation group on K has order p^r and is not a Sylow p-subgroup of S_q. The prime proof identifying multiplier types with full graph types therefore does not apply. Additive K is (C_p)^r, not C_(p^r), so these are not circulants on Z/p^r either.

At q=9 the nonzero even-root graph type is the 3×3 rook graph K_3 square K_3. In F_9=F_3[t]/(t²+1), g=1+t has order eight; the graph with exponent word 00110011 has connection set {±t,±(1+2t)} and primitive half-word 0011 of even weight. In the basis e=t, f=1+2t its edges change one coordinate. The involution (u,v)↦(−u,v) swaps the endpoints of exactly three edges, so it has odd edge sign. Hence this graph is affine-multiplier-even but **not fully sign-even**. Passing to full graph types also changes the equivalence relation. Neither Theorem C nor these counts claim a nonexistence result for a different full-isomorphism target.

For rings Z/p^r, every ring automorphism fixes 1 and hence every element. There is no Frobenius field-automorphism argument there. Also nonzero elements occupy several unit-group orbits. Genuine prime-power circulants require a separate formulation and proof.

## 5. Root-choice torsors: useful obstructions, exact scope

**Primitive-element choice.** A parameter-free uniquely defined element of F_(p^r) is fixed by Frobenius and therefore lies in F_p. For r>1 it cannot generate K^*. Thus no parameter-free ring formula can uniformly choose a primitive element in nonprime finite fields.

**Square-root choice.** Fix any nonsquare b∈F_p and K=F_(p²). The two roots of x²=b lie outside F_p and are exchanged by Frobenius. There is no Frobenius-natural function choosing one of them from b. Thus no parameter-free ring-definable square-root selector can work on all nonzero squares in all odd finite fields, even allowing parameters in the fixed prime field.

**Skew-selector choice.** In every even-degree extension K/F_p, a nonzero element x of the trace-zero line in F_(p²)⊂K satisfies x^p=−x. A Frobenius-invariant subset cannot choose exactly one of x,−x. Therefore a parameter-free ring-definable skew subset of K^* does not exist in this case.

These are obstructions to selecting points of specific torsors, not to every construction that could avoid the selection. A supplied tournament predicate S may itself break a symmetry. Theorem C separately deals with the source-dependent **type bijection**, rather than incorrectly inferring it from the primitive-element obstruction.

## 6. The new positive map is not itself one first-order formula

**Theorem D (uses the classical parity circuit lower bound).** No single first-order formula φ(x;S) in the ring language with a unary input predicate S defines the parity-projection map of Theorem A over all primes p≡3 mod4. This remains true after adding any fixed finite collection of input-independent relations, constants, or field-dependent parameter values.

**Proof.** Fix a field F_p and list its m nonzero squares H. An arbitrary assignment of m bits s(h) defines a valid skew selector by S(h)=s(h), S(−h)=1−s(h), and S(0)=false. The output bit at x=1 is

    w(1)=s(1)+sum_(h∈H)s(h)=sum_(h∈H\{1})s(h).

It computes parity of m−1 independent input bits. Set s(1)=0 if desired.

Expand a fixed first-order formula over the finite universe: each existential/universal quantifier becomes an unbounded OR/AND over p elements; ring operations and all input-independent predicates are hardwired. Its size is p^O(1), with depth bounded by the formula, independent of p. Every atomic S-term becomes a constant, an input bit, or its negation after the above substitution. The resulting constant-depth, polynomial-size Boolean circuits would compute parity on m−1 inputs for arbitrarily large p≡3 mod4. The classical nonuniform AC^0 lower bound for parity forbids this. ∎

The circuit lower bound is attributed to the primary source Johan Håstad, *Almost Optimal Lower Bounds for Small Depth Circuits*, STOC 1986, DOI https://doi.org/10.1145/12130.12132; author-hosted full version: https://www.csc.kth.se/~johanh/largesmalldepth.pdf. The proof above imports only the theorem that fixed-depth polynomial-size AND/OR/NOT circuits cannot compute parity for unbounded input size.

This is not a blanket ring-definability impossibility for the prime-circulant problem. Different bijections can have different complexity. With a supplied generator g, the inherited p≡3 mod4 derivative map has the constant-local formula A(x) iff S(x) XOR S(g²x). The positive map here is generator-free and natural but globally parity-dependent. Naturality, algorithmic uniformity, and first-order definability demonstrably separate.

## 7. Reproducible controls and limits

Run `python3 verify_generator_free_naturality.py --output verification.json` and repeat with `python3 -O`. The checker exhaustively constructs the two orbit sets for m≤13, verifies the positive map on every odd-m input in that range, checks multiplier and automorphism equivariance, exact two-element labelled fibers, the obstruction and its exception, and compares every unit-decimation fixed count to an independently derived Burnside computation. It also computes the prime-power table without enumerating 2^24 or 2^40 words and checks the q=9 rook-graph hostile directly.

For the count computation, the coset Burnside identity gives |(X/C)^α|=(1/|C|)Σ_c|Fix(cα)|. Tournament fixed words are counted by signed component constraints. For graph words, ordinary affine-permutation cycle counts give both total and even-weight fixed words at every divisor period. Subtraction of smaller periods separates exact primitive even and odd weights, and summing primitive-even periods gives the target fixed count. Agreement with exhaustive enumeration validates the formulas on the stated small range; the elementary proofs above establish the general theorems.

No claim here settles: parameter-free first-order definability of some other bijection over all prime fields; a fully sign-even/type bijection for prime-power circulants; the general unlabelled tournament/even-graph natural-bijection problem; or historical priority.
