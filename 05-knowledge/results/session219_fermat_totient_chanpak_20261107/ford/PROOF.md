# Guarded multiplicity, Mersenne normal forms, and exact (3,7)-smooth scaling

Status: PROVED relative to the cited Bang–Zsigmondy theorem and inherited physical O/Z/J grammar. No historical novelty claim. Chan–Pak completeness remains OPEN.

## Inheritance and live board

Sources read: session218-stage AGENTS, research protocol, active guards, THM-4884-definability-lenses (Session218 correction), THM-4749-repunit-degree-law, THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness §6, and Session200 carmichael-fibers/PROOF.md with PHYSICAL_ADJOINT.md. Anchor: constructive tournament spectrum scaling. Niche: exact multiplicative fibers. Wildcard: Ford–Konyagin matching guards. Board: primitive prime support; Mersenne products; strong component factorization; forbidden7 and21; physical source profiles. Closest proved mechanism: ordered strong-component products and the two-run C1 compiler. Hostile: 63 cannot split into smaller strong counts although 63=3²7. Corrected near miss: finite missing multiplicities do not refute a global Ford analogue. Least-used sidecar: primitive-divisor index versus the number of strong components.

## 1. Which 1999 theorem, and what actually transfers?

Ford–Konyagin, *On two conjectures of Sierpinski concerning the arithmetic functions sigma and phi*, Number Theory in Progress (1999), 795–803, Theorem2 proves: if A(m)=k, then for every positive r some l satisfies A(lm)=rk. Thus every positive even multiplicity occurs. Ford's separate Annals1999 paper proves all multiplicities≥2, using a +2 induction. These are not one theorem.

FK Lemma1 forces permitted prime suppliers into a bipartite star-plus-diagonal incidence pattern. The only perfect matchings are the identity and the r−1 transpositions swapping the hub with one leaf: exactly r. Conditions(ii),(iii) exclude repeated suppliers and every unintended prime support. A sieve constructs the admissible arithmetic inputs. The transferable idea is **construct intended decompositions and exclude all unintended decompositions**, not scalar multiplication alone. The star has r perfect matchings because any moved leaf must map to the hub, so at most one leaf moves, forcing its transposition.

Primary source: https://www.ford126.web.illinois.edu/wwwpapers/sigma.pdf (Theorem2, Lemma1). Separate paper: https://www.kurims.kyoto-u.ac.jp/EMIS/journals/Annals/150_1/ford.pdf .

## 2. Complete multiplicative normal form for binary Mersennes

Let M_j=2^j−1 for j≥2. We count finite UNORDERED multisets of indices, including repetitions, whose Mersenne product is N. M_1=1 is excluded, otherwise every fiber is infinite.

**Theorem A.** The only relation needed to generate all equalities between such products is

    M_6 = M_2² M_3 = 63.

Normalize a multiset by replacing each6 with two2s and one3. The normalized multiset is uniquely determined by the product.

**Proof.** Existence of this normalization is elementary. Bang–Zsigmondy says every M_j, j≥2 except6, has a prime divisor absent from every earlier M_i. (At j=2 use3 directly.) Suppose two normalized products are equal. Cancel common factors M_j and take the largest remaining index j. It is not6. Its primitive prime divides the side containing M_j and no factor on the other side, contradiction. Thus the normalized multisets coincide. QED.

Cited primitive-divisor theorem: Glasby, Lübeck, Niemeyer, Praeger, *Primitive prime divisors and the nth cyclotomic polynomial*, J. Aust. Math. Soc.102 (2017), 122–135, introduction, DOI https://doi.org/10.1017/S1446788715000269 (recalls Bang–Zsigmondy). Only its binary specialization is used; the finite verifier does not reprove this all-index theorem.

If a normalized product has c2 copies of2 and c3 copies of3, all original factorizations are parametrized by

    t=0,...,min(floor(c2/2),c3),
    count(6)=t, count(2)=c2−2t, count(3)=c3−t,

with every other multiplicity unchanged. Hence the fiber has exactly

    1+min(floor(c2/2),c3)

members. In particular 63^(r−1) has exactly r unordered Mersenne factorizations for every r≥1. This is an explicit all-multiplicity theorem for this factorization map, NOT Euler phi or tournament classes.

**Homogeneous rigidity.** Normalization changes factor count by2 per replaced6. Thus at any prescribed number of factors there is at most one unordered factorization. Ordered representations are precisely permutations of that multiset. The n≥2 guard is essential: admitting M_1=1 creates the length-three collision M_6 M_1²=M_2² M_3. None of the homogeneous rigidity assertions admits that unit. Excluding index3 altogether makes the Mersenne generators multiplicatively independent: normalized count3 fixes count6, then count2 is fixed.

For Session200's two-run compiler F(O^a Z^b)=2+M_(a+1)M_(b+1), a,b≥1, each target in this sublanguage therefore has precisely one unordered exponent pair. Its source-profile multiplicity within this two-run family is1 on the diagonal and2 off it. This says nothing about other physical words, which provide many additional profiles (e.g. target51).

The loss ledger is decisive: mapping an exponent multiset to a number forgets its length and primitive-index support; normalization restores the latter and the parameter t restores the former. Mapping that number to a Hamilton-path count is NOT automatic: the sole arithmetic relation uses M_3=7, a forbidden count. The relation does not supply a legal tournament factorization of63.

## 3. Exact strong-factorization degree on the entire (3,7)-smooth plane

Let w(h) be the largest number of nontrivial strong tournament components whose Hamilton-path counts multiply to h; set w(1)=0 and w(h)=0 when no representation exists. This is THM-4749's strong-factorization degree, not strong-class multiplicity m(h). The only inherited exclusions required are the proved gaps7 and21. No global Chan–Pak assumption is used.

The inherited physical grammar starts at t=(1,0,2,1), with

    O(D,E,U,d)=(D,0,2U,d)
    Z(D,E,U,d)=(D−E+U,d,2U−d,d)
    J(D,E,U,d)=(D−E+U,D−E,D+U,D).

For every source word its strong C1 closing has H=F=D−E+U. Explicit sources are:

- 3: empty word, profile(1,0,2,1).
- 49: OJJOO, profile(9,0,40,5), closing order8.
- 343: OOJOZJJO, profile(97,0,246,61), closing order11.
- 63: OZZOO, profile(11,0,52,1), closing order8.

A source word of length l closes to a graph of order l+3 by the inherited physical compiler. All orders here are constructive upper bounds, not verified minima.

**Theorem B.** For integers a,b≥0:

- b=0: w(3^a)=a.
- b=1,a=0 or1: the target is7 or21, unattained.
- b=1,a≥2: w(3^a7)=a−1.
- b≥2: w(3^a7^b)=a+floor(b/2).

**Upper bound.** Every component count divides3^a7^b and hence equals3^u7^v. If v=0, then u≥1. If v=1, forbidden7 and21 force u≥2. If v≥2, floor(v/2)≥1. Thus every component consumes at least one unit of u+floor(v/2). Summing gives at most a+floor(b/2). When b=1 exactly one component consumes the7 and at least two3s, leaving at most a−2 other components: at most a−1 total.

**Construction.** For b=0 use a triangles. For b=1,a≥2 use one63 component and a−2 triangles. For even b≥2 use a triangles and b/2 copies of the explicit49 component. For odd b≥3 use a triangles, one explicit343 component, and (b−3)/2 copies of49. Ordered joins multiply H and preserve the strong components. Their counts reach every upper bound. QED.

The packaging construction itself is KNOWN: THM-4094 §6 uses exactly the49/63/343 carry atoms and these parity cases in its two-prime-lane reduction. The extension here is the sharp all-parameter optimality bound and resulting exact degree formula. For nonunit targets, the displayed explicit source words give total realization order 3a for b=0; 3a+2 for b=1,a≥2; 3a+4b for even b≥2; and 3a+4b−1 for odd b≥3. The unit target1 is realized by a singleton tournament, with w(1)=0. Thus63^k has a maximum-component realization of order10k for even k and10k−1 for odd k≥3 (63 itself has order8).

In particular:

    w(63)=1;  w(63^k)=floor(5k/2) for all k≥2,
    w(7^k)=floor(k/2) for all k≥2,
    w(21^k)=k+floor(k/2) for all k≥2.

So squaring repairs two genuinely missing target values:49 and441 are attained although7 and21 are not. More sharply, 63 itself is factorization-blocked, but 63²=3⁴·49 has maximal component length5, not2. Powers of a missing value and powers of an indecomposable attained value behave differently from ordinary exponent scaling.

By the inherited SCC polynomial law, for fixed a,b in the attained cases, the number of n-vertex tournament classes with H=3^a7^b eventually has exactly this degree in n. This supplies an all-parameter extension of the earlier w(63)=1 observation. It is not a claim about all-order strong multiplicities.

## 4. Concrete lesson for a Chan–Pak construction

The wildcard arithmetic factorization problem is completely solvable because primitive primes exclude every unintended source. The tournament plane is solvable because the only unavailable7 factors can be packaged into legal blocks49,63,343. This is the inherited **support-packaging** mechanism with a newly explicit optimality invariant, not a numerical analogy. THM-4094 already reduces global completeness to prime targets other than7 and seven-times-prime targets other than21. The useful next target is therefore direct construction in those two lanes, with valuation cost tracking the unavoidable loss, rather than rediscovering the same packaging. Multiplicative packaging alone cannot reach prime targets: those still require direct strong constructions. Ford-type multiplicity amplification also cannot certify that an empty target fiber becomes nonempty without a construction.

## 5. Verification and scope

Run `python3 verify.py` or `python3 -O verify.py`. The standard-library verifier checks star matchings through r=8; all Mersenne multisets on indices2..18 with length≤5 (26,334 multisets including empty); exact all-index factorization fibers of63^k for k≤5 by divisibility recursion; the two-run compiler for a,b≤20; explicit physical witnesses; and the degree formula on0≤a,b≤14 against independent integer-resource dynamic programming. The DP uses only the four certified component values3,49,63,343, so lower constructions are unconditional. The written upper proof covers all strong values, not just that dictionary. No exhaustive tournament census is claimed.
