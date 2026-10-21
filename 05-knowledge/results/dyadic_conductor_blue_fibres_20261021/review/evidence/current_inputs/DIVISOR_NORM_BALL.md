# The divisor reflection chart is an exact quadratic-order norm ball

21 October 2026. **PROVED:** the all-exponent norm-ball identification, the conductor-ideal rotor obstruction, and the full geometric symmetry group. **KNOWN CLASSICAL SPECIALIZATION, RE-PROVED HERE:** the even-conductor idoneal classification. **FINITE-EXACT:** the literal tables and integer controls. No historical-priority, new general class-number, native-source, tournament, or Collatz claim is made.

## 1. What is inherited, and what this identification adds

The faithful exponent chart and its reflection counts are inherited from [trijection arithmetic, sections 3–4](../trijection_reset_20261020/arithmetic.md). The subsequent [reflection-resolution packet](../divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md) proves an obstruction to transporting those two reflections into an ordinary graph-complement/path-reversal chart and constructs a minimal source-preserving resolution. Neither result is superseded here.

The present extra structure is not another count of ten. The divisor chart is **all nonzero points in a specified norm ball of a specified quadratic order**, with an explicit inverse. The largest part retaining Eisenstein rotation is precisely its intersection with the conductor ideal. This identifies an arithmetic integrality obstruction, while retaining the same prime-labelled divisors and the actual two involutions.

Fix distinct primes p,q,r, an integer k≥1, and

    D_k = {p^a q^b r^c : 0≤a≤k, b,c∈{0,1}}
          minus {1,p^k q r}.

Keep the prime labels throughout. Put ω=(-1+sqrt(-3))/2, so ω²+ω+1=0, and

    O_k = Z + k Z ω,          N(x+tω)=x²-xt+t².

This is the conductor-k order in the Eisenstein field; its discriminant is −3k². Its conductor ideal in O_1 is k O_1. The order and conductor terminology is standard; see Conrad, Example 2.1 and Theorem 1.3 in [The conductor ideal of an order](https://kconrad.math.uconn.edu/blurbs/gradnumthy/conductor.pdf).

## 2. Exact ball theorem, including surjectivity

**Theorem 1.** The map

    Φ_k(a,b,c) = a−kc + k(b−c)ω

is a bijection

    D_k  →  B_k := {z∈O_k : 0<N(z)≤k²}.

It intertwines divisor complement C, prime exchange Q, and their product with

    C: z↦−z,       Q: z↦conjugate(z),       CQ: z↦−conjugate(z).

**Proof.** Write z=x+kyω. Completing the square gives

    N(z)=(x−ky/2)² + 3k²y²/4.

If N(z)≤k², then |y|≤2/sqrt(3)<2. Thus y=0,1,−1. Direct substitution gives precisely

    y=0:    x∈{−k,…,−1,1,…,k};
    y=1:    x∈{0,…,k};
    y=−1:   x∈{−k,…,0}.

These are 2k+(k+1)+(k+1)=4k+2 distinct points. Their unique inverse triples are

    y=1:          (a,b,c)=(x,1,0);
    y=−1:         (a,b,c)=(x+k,0,1);
    y=0, x>0:     (a,b,c)=(x,0,0);
    y=0, x<0:     (a,b,c)=(x+k,1,1).

Every inverse lies in D_k. The two omitted exponent-box endpoints both map to zero, and no other pair in the full box collides. Finally

    C(a,b,c)=(k−a,1−b,1−c),       Q(a,b,c)=(a,c,b).

Substitution proves Φ_k(Cv)=−Φ_k(v). Since conjugate(ω)=−1−ω, it also proves Φ_k(Qv)=conjugate(Φ_k(v)); the product follows. QED.

The inherited chart coordinates are exactly the real and imaginary components, with fixed scales:

    f_k=2a−k(b+c)=2 Re Φ_k,       g=b−c=2 Im Φ_k/(k sqrt(3)).

Thus this is an identification of the previously decoded geometry, not a replacement labelling chosen to fit a norm.

## 3. Reflection axes, parity, and the literal ten-point table

The real axis consists of the 2k nonzero integers in B_k, hence is exactly Fix(Q). On the imaginary axis, x=ky/2. With y=±1 there are two such points iff k is even; y=0 would be the excluded origin. Consequently

    Fix(C)=∅,       |Fix(Q)|=2k,
    Fix(CQ)={±(k/2)sqrt(-3)} when k is even, and ∅ otherwise.

The two CQ-fixed points have norm 3k²/4. Their literal divisors are p^(k/2)q and p^(k/2)r. The two sign partitions are

    sign(Im z):    (k+1), 2k, (k+1);
    sign(Re z):    2k, 2, 2k                    for even k,
                   2k+1, 0, 2k+1              for odd k.

For k=2, O_2=Z[sqrt(-3)]. Write j=sqrt(-3); the complete table is:

| Divisor | (a,b,c) | Φ_2 in Z[j] | Norm | C | Q | CQ |
|---|---|---|---:|---|---|---|
| p | (1,0,0) | 1 | 1 | pqr | p | pqr |
| q | (0,1,0) | −1+j | 4 | p²r | r | p²q |
| r | (0,0,1) | −1−j | 4 | p²q | q | p²r |
| p² | (2,0,0) | 2 | 4 | qr | p² | qr |
| pq | (1,1,0) | j | 3 | pr | pr | pq |
| pr | (1,0,1) | −j | 3 | pq | pq | pr |
| qr | (0,1,1) | −2 | 4 | p² | qr | p² |
| pqr | (1,1,1) | −1 | 1 | p | pqr | p |
| p²q | (2,1,0) | 1+j | 4 | r | p²r | q |
| p²r | (2,0,1) | 1−j | 4 | q | p²q | r |

The 3–4–3 split by imaginary sign and the 4–2–4 split by real sign are therefore literal perpendicular-axis decompositions of one ten-point set. The norm inventory is 1 twice, 3 twice, and 4 six times. This does not identify the several different inherited 3–4–3 arithmetic partitions: prime/squarefree-composite/nonsquarefree, rank, repeated-prime exponent, and imaginary sign are not the same observable.

**Corollary 1 (no hidden global rotation).** For every k≥2, the entire Euclidean isometry group of B_k is the displayed Klein four-group. The same is true of the real-linear norm isometries preserving O_k.

Indeed B_k has centroid zero by negation symmetry, so every Euclidean set isometry fixes zero. The only norm-one points of O_k are ±1: a nonreal point has norm at least 3k²/4≥3. An orthogonal map preserving B_k must therefore preserve the real axis and is one of z,−z,conjugate(z),−conjugate(z). All four do preserve it. The lattice argument uses the same shortest-vector pair. At k=1 the six-point unit hexagon has the larger dihedral symmetry, so this boundary is necessary.

## 4. The exact rotation gate is the conductor ideal

**Theorem 2.** For z=x+kyω∈O_k,

    ωz=−ky+(x−ky)ω ∈ O_k   iff   k divides x.

Hence

    {z∈O_k : ωz∈O_k}=k O_1.

This is the largest additive sublattice of O_k stable under multiplication by ω, and it is exactly the conductor ideal. Moreover

    B_k ∩ k O_1 = k O_1^×
                = {±k, ±kω, ±kω²}.

Its six literal divisors are

    {p^k, q, r, qr, p^kq, p^kr}.

**Proof.** The displayed product makes the divisibility criterion immediate. If L⊆O_k is ω-stable, every element of L satisfies this criterion, so L⊆kO_1. Conversely kO_1 is stable. Since O_1=Z[ω], requiring z and ωz to lie in O_k is also exactly zO_1⊆O_k, the definition of the conductor ideal. If z=ku lies in B_k, then 0<N(u)≤1. The norm is an integer, so u is one of the six Eisenstein units. Encoding the six indicated triples gives precisely those six points. QED.

This is an exact maximality result: outside the six-point core, even one application of ω leaves the original divisor carrier. For k≥2 the smallest witness is p↦1, while ω∉O_k. The outer core itself forms two three-cycles under ω:

    p^k → q → r → p^k,
    qr → p^kr → p^kq → qr.

A 60-degree unit rotation has the same integrality gate, since it differs from ω or ω² by a sign. The rotation preserves the field norm but not this order. Allowing denominator k erases the order obstruction, since O_k[1/k]=O_1[1/k]; doing so changes the retained integrality condition and does not restore the original finite divisor carrier.

## 5. An unconditional even-conductor idoneal classification

For even k, put d=3(k/2)², so disc(O_k)=−4d. Here “idoneal” means the classical integer d, equivalently the proper invertible ideal class group Cl(O_k) has exponent dividing two. This equivalence is Grube's theorem, with the group formulation in Kani, Theorem 6 and Remark 7(b), [Idoneal numbers and some generalizations](https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf). This is a classical family, not a proposed new classification in the literature.

**Theorem 3 (classical specialization, with a direct proof).** For even k,

    Cl(O_k) has exponent dividing two  iff  k∈{2,4,8}.

Equivalently the idoneal integers in the family d=3m² are exactly 3,12,48.

**Proof, step 1: odd factors cannot occur.** The standard order class-number formula is

    h(O_f)=f/3 ∏_(ℓ|f) (1−χ_(−3)(ℓ)/ℓ),        f>1.

The factor 1/3 is essential: O_1 has six units but O_f has only ±1 for f>1, while h(O_1)=1 (the sole reduced form of discriminant −3 is (1,1,1)). One source for the exact sequence yielding this formula is Conrad, Theorem 5.3 in the conductor notes linked above. To see the specialization directly, use O_f/(fO_1)=Z/fZ and CRT. For a prime power ℓ^e, the quotient of unit-group sizes of O_1/(ℓ^e) and Z/(ℓ^e) is ℓ^(e−1)(ℓ−χ_(−3)(ℓ)); this follows from the split, inert, and ramified reductions of X²+X+1. Dividing their product by the unit index three gives the formula.

Write

    f=2^a 3^b ∏_(ℓ≥5) ℓ^e_ℓ,       a≥1.

Because χ_(−3)(2)=−1, it simplifies to

    h(O_f)=2^(a−1) 3^b ∏_(ℓ≥5) ℓ^(e_ℓ−1)(ℓ−χ_(−3)(ℓ)).

For every prime ℓ≥5, χ_(−3)(ℓ)=1 for ℓ≡1 mod3 and −1 for ℓ≡2 mod3; hence 3 divides ℓ−χ_(−3)(ℓ). If b>0 or any such ℓ occurs, three divides h. A finite abelian group of exponent dividing two has power-of-two order. Thus f must be a power of two.

**Step 2: the first three dyadic levels work.** The complete primitive reduced positive-definite forms (A,B,C), representing Ax²+Bxy+Cy², are

| f | discriminant | d | complete reduced forms | class number |
|---:|---:|---:|---|---:|
| 2 | −12 | 3 | (1,0,3) | 1 |
| 4 | −48 | 12 | (1,0,12), (3,0,4) | 2 |
| 8 | −192 | 48 | (1,0,48), (3,0,16), (4,4,13), (7,2,7) | 4 |

Use |B|≤A≤C and B≥0 if |B|=A or A=C. Every class has one representative with this convention; A≤sqrt(|D|/3) bounds the lists. Each listed form is primitive and has B=0, |B|=A, or A=C, so its class equals its inverse. These are the standard ambiguous boundaries. Their counts also equal h(O_f)=f/2, independently certifying completeness. Reduction and its uniqueness for orders are recalled, for example, in [MIT 18.785, 2018 Problem Set 9, Problem 5(b)–(c)](https://math.mit.edu/classes/18.785/2018fa/ProblemSet9.pdf).

**Step 3: explicit non-self-inverse witnesses at every later level.** At f=16, the primitive form (13,8,16) has discriminant −768 and 0<8<13<16. For f=2^a with a≥5, use

    (A,B,C)=(16,8,1+3·4^(a−3)).

Its last coefficient is odd, so the form is primitive; its discriminant is

    8²−4·16·(1+3·4^(a−3))=−3·4^a.

Moreover 0<8<16<C. It is strictly reduced and its inverse (A,−B,C) is a different reduced representative. Thus the class is not self-inverse. This proves the entire remaining tower without a conductor-surjection assumption or a finite-search extrapolation. QED.

For explicit first-failure control, the complete f=16 list is

    (1,0,192), (3,0,64), (4,4,49), (7,−4,28),
    (7,4,28), (12,12,19), (13,−8,16), (13,8,16).

It has eight classes but only four self-inverse classes. Power-of-two class number alone is therefore insufficient.

**Inherited qualification is unchanged.** [The earlier conductor packet, section 6](../escape_renewal_idoneal_20261017d.md), and [THM-4566](../../../01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md) retain the repository's adopted qualification: global completeness of the 65 idoneal numbers is accepted modulo openai/math #003. The restricted theorem above neither invokes that input nor reopens or removes its qualification. Its proof is unconditional within classical quadratic-order and reduction theory. Kani's broader classical results, including Theorem 12, already cover such nonsquarefree specializations.

## 6. What the bridge preserves, and the comparisons it rules out

The following predicates are different, even though each now lives beside a concrete order:

| Predicate | Exact condition in this family |
|---|---|
| CQ has two fixed divisors | k is even |
| Proper-divisor balance F=S+U, k≥2 | k=2 |
| Cl(O_k) has exponent dividing two, k even | k=2,4,8 |
| Entire B_k is stable under ω | k=1 |
| A point of B_k survives ω in the same carrier | It is one of the six conductor-core points |

For the second row, F=4k+2, S=7 and U=3 for k≥2, so F−S−U=4(k−2). At k=1 the number is squarefree and S=6 instead; no extrapolation of that defect formula to k=1 is intended. For example k=4 is idoneal but fails the divisor balance; k=6 has the same two central reflection-fixed points but class number three, so it is not idoneal.

The [native cancellation alphabet](../eisenstein_circuit_completion_20261018g.md), section 1, consists of six norm-one and six norm-three vectors in **O_1**. It is not the ten-point ball in **O_2**. In particular, six points of B_2 have norm four, which never occurs in that alphabet. There is no norm-preserving map from all of B_2 into that alphabet. Nor can a common scalar rescaling fix the three-versus-two distinct norm levels. Removing the conductor integrality condition merely changes the target; it does not give the original ten-point set rotation closure.

The source-to-ball map preserves:

- each literal divisor, through the inverse and the retained prime labels;
- divisor complement, prime exchange, their product, and their fixed sectors;
- the specified full chart coordinates, now an exact field norm and conductor.

It does **not** transport multiplication as ring multiplication: Φ_k(p)=1 while Φ_k(p²)=2 for k≥2. Divisibility can be recovered by decoding exponents, not by ordinary planar coordinatewise order. A ball point is not an ideal class, and no canonical group law on the divisor set is supplied. The ideal class group uses invertible ideal classes, not a norm-ball point count. Finally no native source ID, endpoint, integer braid clock, or cancellation certificate is created by this static identification. Those sidecars remain necessary before any dynamical comparison is meaningful.

## 7. Exact reproduction and hostile controls

Run from this directory:

    python3 -B verify_norm_ball.py > /tmp/norm-ball.json
    python3 -B -O verify_norm_ball.py > /tmp/norm-ball.optimized.json
    cmp /tmp/norm-ball.json /tmp/norm-ball.optimized.json
    cmp /tmp/norm-ball.json verification.json

The checker has no third-party dependencies and uses explicit exceptions rather than removable assertions. It checks 72,384 live conditions. Both modes give byte-identical receipts. It exhausts B_k for k=1,…,64, tests direct integer-divisor enumeration for three marked prime assignments and k=1,…,5, verifies the entire-lattice rotor criterion on declared boxes, checks all primitive reduced-form lists for conductors 1,…,128 against the independent class-number product, and checks the symbolic witness identities at dyadic exponents 4,…,128. The all-k and all-conductor conclusions rely on the written proofs, not these finite ranges.

Hostiles retained in the checker include the omitted-endpoint collision, p↦1 leaving O_2 under rotation, the norm-four obstruction to the native alphabet, an imprimitive false reduced class, and f=16 having power-of-two class number without exponent two. The literal ten-row reflection table and complete small-conductor form lists are emitted in the receipt. Source hashes, exact source sections, and a publication manifest accompany this packet.
