# Distinct magic arrays: exact torsion counts and the doubled-image criterion

**Status:** PROVED elementary finite-abelian-group enumeration, independently checked by complete finite Smith certificates; RECOVERED/REFORMULATED existence classification, explicitly attributed below. No arithmetic-square solution or external priority claim. This third packet leaves the earlier base and branch packets unchanged.

## 1. The exact distinctness statistic

Let A be a finite abelian group, N=|A|, and h_d=|A[d]|. Let I(A) count ordered pairs(x,z)∈A² for which the nine offsets

    {i x+j z :i,j∈{−1,0,1}}

are pairwise distinct. The universal group-magic normal form from `RESULT.md` has precisely these offsets, in a fixed magic placement, about its center E.

**Theorem D1.**

    I(A)=N²−4(1+h_2)N−2−4h_2+4h_4+2h_5+4h_6
                       +2h_2²+h_2 h_4.                (1)

Consequently the number of labelled injective group-magic arrays is N·I(A); the number additionally having no zero entry is(N−9)·I(A). In particular I(A)=0 whenever N<9, so the latter expression is never a negative count.

### Direct human proof: four corners and two5-torsion intersections

First require the four corners±x±z to be distinct. This is equivalent to2x,2z,2x+2z,2x−2z all being nonzero. Put B=2A, M=|B|=N/h_2. The doubling map A[4]→B[2] is surjective with kernel A[2], so |B[2]|=h_4/h_2. For u∈B nonzero, the number of v avoiding0,u,−u is M−2 if u has order2, and M−3 otherwise. Thus the count of such(u,v) is

    (|B[2]|−1)(M−2)+(M−|B[2]|)(M−3)
      =M²−4M+|B[2]|+2.

Every pair(u,v) has h_2² lifts( x,z), so the number of corner-good pairs is

    C=N²−4h_2 N+2h_2²+h_2 h_4.

Among corner-good pairs, the only remaining collision conditions are the four equations

    z=2x,  z=−2x,  x=2z,  x=−2z.

On any one of these lines, corner-goodness excludes exactly A[4]∪A[6] in its free parameter. Its size is therefore N−h_4−h_6+h_2, since A[4]∩A[6]=A[2]. Of the six pair intersections, four force order dividing3 or4 and are already corner-bad. The only surviving intersections are(z=2x,x=−2z) and(z=−2x,x=2z). Each is parametrized by a nonzero element of A[5], so has h_5−1 points. The two intersections are disjoint and no triple intersection remains. Inclusion–exclusion gives

    I(A)=C−4(N−h_4−h_6+h_2)+2(h_5−1),

which is exactly(1). This proves the all-group formula without any finite computation and explains geometrically why5-torsion first appears on imposing distinctness. ∎

### Why eight collision characters suffice

A collision between two offsets is a zero value of a nonzero integer linear form ax+bz with a,b∈{−2,−1,0,1,2}. Up to simultaneous sign there are twelve such forms. The kernels for(1,0),(0,1),(1,1),(1,−1) are contained in the kernels for their doubled forms. Removing those redundant four leaves exactly

    D={(2,0),(0,2),(2,2),(2,−2),(2,1),(2,−1),(1,2),(1,−2)}.   (2)

This reduction is valid for EVERY abelian group: it uses the implication v=0⇒2v=0, never its converse. Thus I(A) is the complement of these eight kernels in A².

Apply Smith inclusion–exclusion, exactly as in `branch-addendum.md`. The nonzero signed table is

| Smith factors | Signed coefficient |
|---|---:|
|empty|1|
|1|−4|
|2|−4|
|1,1|−2|
|1,2|−4|
|1,4|4|
|1,5|2|
|1,6|4|
|2,2|2|
|2,4|1|

A type(d_1,...,d_r) contributes N^(2−r)∏h_(d_i), independently recovering(1). All256 subsets of(2) are certified, independently cross-checked against all4096 subsets of the unreduced twelve-form list. The rank-two invariant factors are calculated from the gcd of entries and gcd of pair determinants. A separate parent implementation uses only the reduced eight-row list and direct group enumeration.

For a successful gap pair, all E∈A give distinct entries by translation. Exactly nine centers place zero at one of the nine positions, proving the two center factors. No quotient by translations, D4, or group automorphisms is taken.

**What distinctness adds.** The earlier nonzero-entry formula uses h_2,h_3,h_4. Distinctness also detects h_5 and new mixed torsion terms. Because2 and3 are coprime, h_6=h_2 h_3 for EVERY abelian group;6-torsion is not an independent extra coordinate. For example C25 and C5² have the same order and the same2-,3-,4-torsion data, hence both have10680 nonzero-entry magic arrays, but their injective-gap counts are440 and480. Their total injective magic counts are11000 and12000 respectively. Nonzero entry counts do not determine injectivity.

## 2. Exact finite-field and cyclic consequences

For the additive group of F_q,

    I(F_q,+)=
      0                         if char=2;
      (q−1)(q−3)                if char=3;
      (q−1)(q−5)                if char=5;
      (q−1)(q−7)                if char≥7.             (3)

This can also be proved directly: in odd characteristic the eight rational collision directions are the projective slopes0,∞,±1,±2,±1/2. They remain eight directions for characteristic≥7, merge to six at5 and to four at3; a union of m distinct lines through zero in F_q² has1+m(q−1) points. In characteristic2 opposite offsets already coincide.

Thus over F9 the number of labelled injective arrays is9·8·6=432, recovering Session57's count without identifying F9 with Z/9Z. For C9, formula(1) gives I=24 and only216 labelled injective arrays. In F9 every injective array uses all nine field elements, so none has all entries nonzero, consistent with the factor q−9.

For a cyclic group C_n, substitute h_d=gcd(n,d) in(1). This gives a quasipolynomial in n of period60. It is not obtained by treating Z/nZ as a field. For odd n it simplifies to

    I(C_n)=n²−8n+1+2gcd(n,5)+4gcd(n,3).

For even n it is

    I(C_n)=n²−12n−2+6gcd(n,4)+2gcd(n,5)+4gcd(n,6).

For the generic multiplicative branch complement with pairwise distinct entry coordinates, put n=q−1. Its finite-field point count is

    (n−9) I(C_n).                                  (4)

The additive special-fiber count for distinct nonzero entries is(q−9) times(3). These are base-parameter counts, not rational points of the512-sheet ordinary-square cover. At q=11 they are8 generically and80 additively, a concrete test of the distinction.

## 3. A compact existence criterion—known classification, fresh formulation

**Theorem D2.** A finite abelian group A admits a3×3 magic array of nine distinct elements iff

    |A|≥9  and  |2A|≥4.                            (5)

Here2A={2a:a∈A} is the image of doubling, NOT A[2]. In terms of the notation above, |2A|=N/h_2.

**Necessity.** Nine distinct elements require N≥9. The four corner offsets±x±z are distinct only if2x,2z,2x+2z,2x−2z are all nonzero. No group of order1,2 or3 contains u,v with u,v,u+v,u−v all nonzero. Thus |2A|≥4.

**Sufficiency, with explicit witnesses.** Let m be the exponent of A. A finite abelian group has an element g of order m. If m≥9, take x=g,z=3g: the offsets are the nine consecutive multiples−4g,...,4g and are distinct.

It remains to consider m≤8. Exponents1 and2 violate |2A|≥4. If m∈{3,5,7}, the group is an elementary abelian p-group; N≥9 forces at least two independent cyclic factors, and choose x,z in those factors. If m=4, |2A|≥4 forces at least two C4 factors; their independent generators work. If m=6, |2A|≥4 forces at least two C3 factors; their independent generators work. If m=8, N≥9 forces A to contain C8×C2. In that subgroup choose x=(1,0),z=(2,1). The three offsets with second coordinate0 are(−1,0),(0,0),(1,0); the six with second coordinate1 have first coordinates1,2,3,5,6,7. They are all distinct. These cases exhaust the finite exponents below9. ∎

An injective array avoiding zero exists iff N≥10 and |2A|≥4: after a successful gap pair, choose a center outside its nine forbidden values.

**Attribution.** Bowerman–Fleece–Insall, *Introducing n-Magic Groups and Characterizing3-Magic Finitely Generated Abelian Groups*, arXiv:2308.01858v1, Theorem3.1, already proves the finite-abelian existence classification in a Sylow-factor case form. Definition1.1 explicitly requires distinct entries. Primary HTML/PDF checked at https://arxiv.org/html/2308.01858v1 and https://arxiv.org/pdf/2308.01858 . Propositions2.6,2.7,2.15 and2.17–2.18 give the relevant cyclic, product and forbidden-family mechanisms. Statement(5) is a compact doubled-image reformulation with a direct proof, NOT a new existence theorem. The exact torsion enumeration(1) is a separate calculation; no assertion that the literature lacks it is made.

## 4. Optional topological consumer: adding collision divisors

This is a **CITED application** of the Poincaré specializations already sourced precisely in `branch-addendum.md`, not an additional general theorem. Delete both the nine zero-entry divisors and the collision divisors. Over C, their additive characteristic polynomial is

    χ_inj,+(u)=(u−9)(u−1)(u−7)
              =u³−17u²+79u−63.

For the multiplicative arrangement, the fully split cyclic constituent has n divisible60. Then I(C_n)=n²−12n+56, so its arithmetic characteristic polynomial is

    χ_inj,×(u)=(u−9)(u²−12u+56)
              =u³−21u²+164u−504.

The corresponding complex Poincaré polynomials are

    P_inj,+(s)=1+17s+79s²+63s³,
    P_inj,×(s)=1+24s+209s²+690s³.

The extra four geometric toric components reflected in the coefficient21 rather than17 come from the four doubled collision characters: each kernel has a +1 and a −1 component over C. Counting only their tangent hyperplanes forgets those torsion components. Euler characteristics are0 and−504. This describes the open base arrangement, not the square-root cover's topology.

## 5. Reproduction and scope

Run `python3 check_distinct_magic.py` and `python3 -O check_distinct_magic.py`; only the standard library is needed. All tests use explicit exceptions. The source certifies the complete12-row and reduced8-row subset tables, verifies that every grid collision occurs in the twelve-row set, compares direct offsets with both kernel descriptions for75 finite groups (cyclic orders1–60 and15 noncyclic products), explicitly counts admissible centers for groups of order at most25, and tests the formula/existence criterion on every invariant-factor model of order at most512. The latter is a bounded consistency check, not the proof of(5).

The normal and optimized certificates/stdout agree byte for byte. The existence theorem is independently proved in §3; the all-group enumeration has the direct corner-and-line proof in §1 and the independent complete finite Smith certificate. Squarehood, positivity, distinct entry values versus distinct roots, and generic rational sections remain separate conditions.
