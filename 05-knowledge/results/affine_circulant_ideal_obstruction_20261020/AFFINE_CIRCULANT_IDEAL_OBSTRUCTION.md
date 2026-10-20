# The Boolean-affine boundary for prime-circulant type bijections

2026-10-20. **PROVED, with exact finite controls and completed independent audit.** No historical priority claim. This packet distinguishes an affine function of the input Boolean coordinates from a possibly nonlinear function equivariant under the affine group of a prime field.

## 1. Result and inheritance

Let p=2m+1 be an odd prime. There is a Boolean-affine, multiplier-equivariant map from labelled circulant tournaments to labelled sign-even circulant graphs inducing a bijection of isomorphism types **if and only if m is odd or m=2**. Thus it exists for p≡3 mod4 and p=5, and does not exist for p≡1 mod4, p>5.

The obstruction applies even when the rule may choose a primitive root or use field addition in its fixed coefficients. It does not impose naturality under all automorphisms of the abstract multiplicative group. It does not obstruct nonlinear AGL(1,p)-equivariant rules: the inherited derivative-path construction already supplies one for every odd prime.

Dependencies and inheritance:

- Closest proved mechanism: the necklace dictionary and all-prime nonlinear construction in [CIRCULANT_NECKLACE_BIJECTION.md](../circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md), especially §§2–3 and 5–7. That packet's completed audits are separate files; its initial pending-audit header is not a retraction.
- Least-used sidecar: the stable ideal and repeated x+1 primary factor in [BOUNDED_FITTING_ROUTING.md](../bounded_fitting_routing_20261020/BOUNDED_FITTING_ROUTING.md).
- Distinct earlier obstruction: [GENERATOR_FREE_NATURALITY.md](../generator_free_circulants_20261020/GENERATOR_FREE_NATURALITY.md) rules out abstract-multiplicative-automorphism naturality when m is even and m>2. That obstruction permits arbitrary nonlinear maps but requires more symmetry. The present obstruction permits fewer symmetries but requires Boolean-affineness.
- Corrected near miss: “affine-equivariant bijection” is not a labelled bijection. The relevant forward map is a labelled surjection inducing a type bijection. Already p=3 has two source objects and one target object.
- Hostiles: the all-one word, the adjacent-pair word 1+x, and the exceptional one-type cases m=1,2.

Anchor: classify affine rules. Niche: identify the maximal cyclic linear code inside the target. Wildcard: connect primitive-root parity to nilpotent height and cyclic order. The concept board is the antiperiodic affine torsor; augmentation hyperplane; cyclic ideals and their socles; nilpotent height; and labelled/type quotient losses.

## 2. Exact objects and affine convention

Let X_m be the affine space of binary words f on Z/(2m) satisfying

    f_(i+m)=f_i+1.

Let T rotate this word by one place. Choosing any antipodal transversal identifies X_m with F₂^m; changing the transversal permutes and possibly complements coordinates, so the meaning of a Boolean-affine map is independent of this choice.

Let V_m=F₂^m, with cyclic rotation R. Put

    B_m={w∈V_m : the primitive root of w has even Hamming weight}.

The source action is C_(2m) acting by T; its action on V_m factors through C_m acting by R. A Boolean-affine map means F(f)=b+L(f) in a transversal coordinate, with b constant and L linear over F₂. Multiplier-equivariance means F(Tf)=R F(f).

At p=2m+1 prime, X_m describes circulant tournaments and B_m describes sign-even circulant graphs. The inherited Sylow/normalizer argument identifies multiplier orbits with actual isomorphism types. Translations act trivially on these connection-set descriptions, so multiplier-equivariance is exactly AGL(1,p)-equivariance here. “Sign-even” concerns the sign of the induced action on oriented graph edges; it does not mean even vertex degree.

Write

    A_m=F₂[x]/(M),   M=x^m−1,
    m=e q,   e=2^v₂(m),   q odd,
    t=x+1,   M=t^e Q^e,   Q=(x^q−1)/(x+1),
    U=t^e A_m.

Here x is the formal polynomial variable; h(R) means evaluation of that polynomial at the rotation operator R. Identify a word w with its polynomial Σ_i w_i x^i, so multiplication in A_m corresponds to cyclic convolution. The ideal U is the invertible/stable summand for D=I+R. It has dimension m−e. For m=e, Q=1 and U={0}.

## 3. Every affine equivariant rule is one convolution after differentiation

**Theorem 1 (complete normal form).** Every Boolean-affine multiplier-equivariant map X_m→V_m is uniquely

    F_h(f)=h(R) Δf,     h∈A_m,

where Δf=(I+T)f, regarded as an m-periodic word. Consequently there are exactly 2^m such rules.

**Proof.** The word Δf has odd total weight by telescoping. Conversely every odd-weight m-word has exactly two antiperiodic integrations, differing by the all-one word. These integrations are interchanged by T^m.

Equivariance gives F(f+1)=F(T^m f)=R^m F(f)=F(f). Thus F descends to an affine map on the odd augmentation hyperplane

    H={u∈V_m : ε(u)=1},     ε(u)=sum_i u_i.

Any affine map on H extends uniquely to a linear map on V_m: if H=u₀+E with E=ker ε and the affine map is u₀+v↦b+A(v), define its extension by v+c u₀↦A(v)+c b. Uniqueness follows because H spans V_m. Rotation preserves H, so equivariance on H implies that the extension commutes with R. Every linear endomorphism commuting with the regular cyclic shift is a unique polynomial h(R), because its value on one cyclic basis vector determines all its columns. Conversely every such polynomial gives an affine equivariant map after Δ. ∎

When writing Δ² or a higher Δ power below, the first application is the affine map Δ:X_m→V_m and every later application is the ordinary cyclic derivative D=I+R on V_m; thus Δ^s means D^(s−1)Δ. This classification includes constants and arbitrary global linear combinations. It is not merely a test of monomial derivatives or bounded-local formulas.

## 4. The ideal obstruction and sharp existence boundary

**Lemma 2 (maximal cyclic code in the target).** The ideal U is the unique largest rotation-invariant linear subspace contained in B_m.

**Proof.** First U⊆B_m. If w has primitive period d|m, set E=2^v₂(d), and write its primitive root polynomial as v. Repetition gives

    w=v · (x^m−1)/(x^d−1).

The second factor has t-valuation e−E. If w∈U, then v has t-valuation at least E, unless w=0; in particular v(1)=0. The primitive weight is even.

Conversely let I be a rotation-invariant linear subspace contained in B_m. It is an ideal (g) of A_m for a divisor g of M. The all-one word is

    J=1+x+⋯+x^(m−1)=M/t,

and is not in B_m, since its primitive root is 1. If t^e does not divide g, then g divides M/t, so J∈I, a contradiction. Hence I⊆(t^e)=U. ∎

**Theorem 3 (no affine type-surjection for even m>2).** If m is even and m>2, no Boolean-affine equivariant map X_m→B_m is surjective on multiplier types.

**Proof.** The image of an affine map is an affine subspace. Surjectivity on types and equivariance imply surjectivity onto every labelled member of B_m, because the image contains every rotation of each one of its elements. In particular the image contains 0 and is therefore a rotation-invariant linear subspace. Lemma 2 forces it into U.

But w=1+x belongs to B_m: for m>2 the adjacent pair has least period m and primitive weight two. Its t-valuation is exactly one, whereas e≥2. Thus w∉U, a contradiction. ∎

There is an even shorter independent obstruction to a labelled affine surjection, without equivariance. B_m contains 0 and each adjacent-pair word x^(2j)(1+x), for 0≤j<m/2. Their sum is J, which is excluded. Thus B_m is not an affine subspace for even m>2. Equivariance is needed to turn type-surjectivity into labelled surjectivity; the argument does not prohibit a nonequivariant affine map that merely meets every type.

**Theorem 4 (all affine type-bijections).**

1. For odd m, F_h induces a type bijection if and only if

       gcd(h,x^m−1)=x+1.

   Every such map has exactly two preimages per labelled output, and the two inputs are complementary antiperiodic words.
2. For m=2, the sole affine type-bijection is the constant zero map. Its labelled fiber has size four.
3. For even m>2, there are none.

**Proof.** If m is odd, the repetition count of every primitive root is odd, so B_m=ker ε=U. Also x^m−1 is squarefree. If gcd(h,M)=t, multiplication by h has image U and kernel the one-dimensional space spanned by J. Since ε(J)=1, its restriction H→U is bijective. Combining with the two-to-one Δ proves the claim, including the type statement, since Δ identifies exactly complementary words within a source orbit.

Conversely a type-bijection is a labelled surjection. Thus F_h(H)=U. The full linear image hA_m is the span of F_h(H), because H spans A_m; therefore hA_m=U, equivalently gcd(h,M)=t. For m=1 this condition says h=0, so p=3 is included. For m=2, B_m={0}; uniqueness follows from Theorem 1. The remaining case is Theorem 3. ∎

In particular h=t gives the familiar F_h=Δ² for odd m. It is not the only affine rule: their number is |(F₂[x]/Q)^×|, equivalently the product of 2^deg(P)−1 over the distinct irreducible factors P of Q. The empty product at m=1 is one.

## 5. Full classification of target-valued affine rules

The no-go does not mean every affine rule fails to land in the target. Many do, and their exact limitation is controlled by one primary height.

Define, with v_t(0)=∞,

    k(w)=e−min(v_t(w),e),       0≤k(w)≤e.

This is the nilpotent height of w's t-primary component: the least k for which D^k w∈U. It is rotation-invariant. It also does not depend on the chosen generator of the cyclic coordinate: replacing x by x^a, gcd(a,m)=1, preserves the t-primary filtration. When m is even, a is odd and (x^a+1)/(x+1) is a unit at x=1; when m is odd, the capped valuation is simply augmentation parity.

**Lemma 5 (parity as an order/height comparison).** If w has least period d and E=2^v₂(d), then

    w∈B_m  iff  k(w)<E.

**Proof.** In the repetition factorization used in Lemma 2,

    min(v_t(w),e)=e−E+min(v_t(v),E).

Therefore k(w)=E−min(v_t(v),E). This equals E exactly when v(1)=1, and is smaller exactly when the primitive root has even weight. The zero word has period one and k=0, so it is covered. ∎

**Theorem 6 (complete admissibility classification).** For h∈A_m put

    k=e−min(v_t(h),e).

Then F_h(X_m)⊆B_m if and only if k=0 or k is not a power of two. Moreover

    0∈F_h(X_m)  iff  k=0.

**Proof.** The Chinese remainder decomposition is

    A_m ≅ F₂[t]/(t^e) × F₂[x]/(Q^e).

An odd word u has a unit in its first component; its second component is unrestricted. Thus every element of hH has the same first-component height k. If k=0, all outputs lie in U. If k is positive and not a power of two, Lemma 5 makes every output target-valued: k≤E and E is a power of two, so k<E.

If k is a power of two, choose u with t-primary component 1 and other component 0. It has odd augmentation. Then hu has only its t-primary component and has height k. On that component x=1+t, and (1+t)^k=1+t^k. Its least period is exactly k: no smaller power of two kills a vector of height k, and the order of 1+t is a power of two. Lemma 5 says hu is not in B_m. This supplies an explicit primary-component hostile.

Finally a unit times a nonzero t-primary component cannot vanish. If k=0, the same choice of u gives hu=0. ∎

**Corollary 7 (counts and ranks).** The number of admissible affine rules is

    2^(m−e) · [ 2^e − Σ_(j=0,...,log₂e) 2^(2^j−1) ].

Exactly 2^(m−e) of them hit zero. The largest affine image dimension is m−1 for odd m and m−2 for even m.

**Proof.** There is one zero t-primary component and 2^(k−1) elements of exact height k>0. The other primary component has 2^(m−e) choices. Sum over the allowable heights from Theorem 6. If k=0, the image is an ideal in U and has dimension at most m−e, achieved by h=t^e. If k>0, the image is an affine shell of dimension at most k−1+(m−e). For even m, the maximal allowable k is e−1 when e≥4, giving m−2, achieved by h=t. When e=2, only k=0 is allowable, again giving m−2. Odd m has e=1 and the usual rank m−1. ∎

Examples and sharp hostiles:

- m=4, h(x)=x+1: e=4 and k=3. The image consists of the four rotations of 0011, all target-valued, but omits 0000. This smallest non-power-of-two-height example explains why admissible affine maps need not have linear images. It is a necklace example, since 2m+1=9 is not prime.

- m=6, p=13: only height zero is allowable. There are 16 admissible affine rules among 64 equivariant affine rules. The stable image has 16 labelled words, while B₆ has 28.
- m=8, p=17: allowable heights are 0,3,5,6,7. There are 117 admissible affine rules among 256. The only one hitting the empty graph is the constant zero map, because U={0}. The largest safe affine image has 64 labelled words; B₈ has 117.
- At e≥4, Δ² is target-valued but never reaches the empty graph. At e=2, Δ² is not even target-valued. Thus “the fixed derivative fails” has two different mechanisms.
- Admissibility is not monotone under adding derivatives. For e=8, Δ², Δ³, Δ⁴ are target-valued (heights 7,6,5), but Δ⁵ is not (height 4).

## 6. Exact best type coverage when the empty graph is retained

An affine admissible map hitting zero has image in U, and h=t^e attains all of U. Thus the largest number of target types such a rule can cover is exactly the number of rotation orbits in U:

    N_U(m) = (1/m) Σ_(d|m) φ(m/d) · 2^(d−2^v₂(d)).

Indeed, words fixed by a shift with gcd(m,j)=d are repetitions of length-d words. Their intersection with U has dimension d−2^v₂(d), by the same repetition valuation calculation. Burnside's averaging yields the formula. This counting formula is independent of the no-go proof.

More generally, the shell C_k={w:k(w)=k} contains

    N_k(m) = (2^(k−1)/m)
             Σ_(d|m, 2^v₂(d)≥k) φ(m/d) · 2^(d−2^v₂(d))

rotation orbits for k>0. Each allowable shell is attained by an admissible affine rule: h=t^(e−k) attains C_k for k>0, and h=t^e attains C_0=U. All other rules of the same height have images contained in that shell. Hence the exact maximum number of types covered by any admissible affine rule is

    max { N_U(m), N_k(m) : 1≤k≤e and k is not a power of two }.

This is a sharp finite divisor formula, not an asymptotic estimate. For the first obstructed primes, the maximum with the empty graph retained is 4 of 6 types at p=13, and 1 of 16 at p=17. At p=17, allowing the empty graph to be missed raises the maximum to 8 of 16 types.

### A closed sharp coverage law

Let A(m) denote the common number of source and target necklace types. Direct
Burnside counting on odd-total words gives

    A(m) = (1/m) Σ_(d|m, v₂(d)=v₂(m)) φ(m/d) 2^(d−1).

**Corollary 8 (optimal affine type coverage).** Among all target-valued affine
equivariant rules, the largest number of target types covered is

- A(m), if m is odd;
- [A(m)+A(m/2)]/2, if v₂(m)=1;
- A(m)/2, if 4 divides m.

For 4|m the optimum is attained by h=t, namely Δ². It misses the empty graph,
except that at the necklace-only case m=4 the constant-zero rule ties its
one-type coverage. For v₂(m)=1 only zero-height rules are admissible, and
h=t² attains the optimum. In particular, this formula includes m=2.

**Proof.** The odd case is Theorem 4. If e=2, only height zero is admissible.
Writing m=2q, q odd, the two possible dyadic parts d=r and d=2r in the formula
for N_U give

    N_U(2q) = (1/(2q)) Σ_(r|q) φ(q/r)
                          [2^(r−1)+2^(2r−2)]
            = [A(2q)+A(q)]/2.

Now let e≥4. The largest admissible height e−1 has only the terms with
2^v₂(d)=e in its orbit formula, and hence N_(e−1)(m)=A(m)/2.
To prove maximality, write every divisor as d=Er with E|e and r|q. For a
positive admissible k let E₀ be the smallest power of two at least k. At each
fixed r its contribution to mN_k is bounded above by

    φ(q/r) 2^(k−1) 2^(e(r−1)) Σ_(E₀|E|e) φ(e/E)
      = φ(q/r) 2^(e(r−1)) 2^(k−1) e/E₀.

Since k is not a power of two, k<E₀, and e/E₀≤2^(e−k−1). The last expression
is at most the corresponding contribution
φ(q/r) 2^(e(r−1)) 2^(e−2) to mN_(e−1). For k=0 the same estimate has
Σ_(E|e)φ(e/E)=e, and e≤2^(e−2) for e≥4. Thus neither the zero shell nor any
other admissible shell exceeds N_(e−1). For e=4, equality for the zero shell
requires every r=1, hence q=1; for e≥8 its bound is already strict. ∎

## 7. The connection to bounded nonlinear routing

The same distinguished primary block provides three precise facts:

1. Its length e is the nilpotence depth of D=I+R on the transient summand.
2. Its vanishing ideal U=im(D^e) is the unique maximal cyclic linear code in B_m.
3. The primitive-root sign is the comparison k(w)<2^v₂(period(w)).

The inherited nonlinear derivative-path map can leave one affine shell and choose different derivative depths for different source orbits. An affine rule F_h cannot do so: every output has the one fixed t-primary height determined by h. A type-surjection must hit both the zero-height empty graph and the positive-height adjacent-pair graph, so this fixed-height restriction already identifies the missing coordinate.

Connection contract:

- Source: antiperiodic selector torsor with multiplier action.
- Map: Δ, then a cyclic convolution h(R), then passage to multiplier orbits.
- Preserved predicates: multiplier-equivariance; odd augmentation at the Δ stage; one fixed nilpotent-height shell after h.
- Lost information: complementation at Δ and, for nonunit h on the stable factor, additional linear coordinates.
- Needed nonlinear sidecar: a source-dependent derivative stopping/branch choice, or some other operation capable of moving between primary heights.
- Decisive hostile: empty graph plus the adjacent pair. The former forces height zero for an affine rule; the latter requires positive height when m is even.

This is a proof of necessary nonlinearity in the Boolean coordinates for the stated prime family. It is not a Boolean-degree lower bound stronger than two, a circuit-size bound, a first-order nondefinability theorem, or an obstruction to using field addition. Those would require further arguments.

## 8. Exact verification and remaining scope

The accompanying standard-library verifier independently enumerates every h and every input word for m≤12; it checks normal-form injectivity, equivariance, actual primitive roots, the height criterion, target landing, zero occurrence, ranks, full labelled surjectivity, induced type bijections, and two-element labelled fibers in the positive family. For m≤3 it additionally enumerates **every** matrix and translation of an arbitrary affine map, without using the polynomial normal form, and verifies that the equivariant subset is exactly the stated family.

The adjacent-pair non-affineness witness is checked separately. Finite tables are labelled FINITE-EXACT; the theorems are elementary all-m proofs. No finite cutoff is used as evidence for an unproved universal claim.

Reproduce:

    python3 verify_affine_circulants.py --max-m 12 --output verification.json > verification.out
    python3 -O verify_affine_circulants.py --max-m 12 --output verification.optimized.json > verification.optimized.out
    cmp verification.json verification.optimized.json

The checker uses explicit exceptions rather than assert, so optimized Python does not disable checks. See VALIDATION.md for the completed-run receipt and proof audit. The previous staged receipt was stale only in its check count and script hash; it has been replaced by the matching completed rerun. No numerical table row changed.

The necklace and affine classification holds for all m. Its interpretation as a bijection of actual graph/tournament isomorphism types is asserted only when 2m+1 is prime. No composite-order or prime-power graph-isomorphism extension is inferred.
