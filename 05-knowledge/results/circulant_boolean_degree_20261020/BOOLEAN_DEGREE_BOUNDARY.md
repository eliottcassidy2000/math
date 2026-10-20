# Unbounded Boolean degree for equivariant circulant type bijections

2026-10-20. **PROVED**, with independent exact controls. No historical priority claim. The lower-bound argument is elementary and does not depend on the exhaustive quadratic search.

## 1. Result and inheritance

Let `X_m` be the binary antiperiodic selector torsor, let `T` be its twisted cyclic shift, and let `B_m` be the length-`m` binary words whose primitive root has even weight, with ordinary cyclic shift `R`. Consider maps

    F : X_m -> B_m,     F(Tx)=R F(x),

which induce bijections between source and target cyclic orbit sets. The degree of `F` is its Boolean algebraic degree in any antipodal transversal coordinates.

Write `m=e q`, with `q` odd and `e=2^a`. Then:

1. If `e>=4`, every such map has

       degree(F) >= ceil((e-1)/a).

2. The necessary degree is unbounded even when `2m+1` is required to be prime. This does **not** assume infinitely many Fermat primes.
3. At `m=6`, corresponding to `p=13`, the minimum degree is exactly two.
4. At `m=8`, corresponding to `p=17`, the minimum degree is exactly three. A seven-monomial cubic rule is given below.
5. At the necklace-only case `m=4`, the minimum degree is exactly two.

The first examples of the general lower bound are:

| Prime p | m=(p-1)/2 | Dyadic part e | Necessary degree |
|---|---:|---:|---:|
| 17 | 8 | 8 | 3 |
| 97 | 48 | 16 | 4 |
| 193 | 96 | 32 | 7 |
| 641 | 320 | 64 | 11 |
| 257 | 128 | 128 | 19 |
| 65537 | 32768 | 32768 | 2185 |

These are lower bounds, not claimed exact minima, except at `p=17`.

Inheritance and hostile controls:

- Closest proved mechanism: [CIRCULANT_NECKLACE_BIJECTION.md](../circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md), which supplies nonlinear equivariant type bijections for every `m` and identifies cyclic types with actual prime circulant isomorphism types.
- The immediately preceding affine boundary: [AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md](../affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md). Its maximal cyclic ideal obstruction proves degree at least two for even `m>2` but does not by itself force unbounded degree.
- Related established mechanism: [THM-4644, Fermat-prime circulant rigidity](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md). Here the free action at a dyadic length is used to force a truth-table weight, rather than to classify tournament automorphism groups.
- Least-used sidecar: the zero-output type, whose missing complement creates an exact, nonzero Walsh defect.
- Corrected near miss: source and target **type** bijection is not a labelled bijection. At `p=17` the successful cubic has labelled fibers of sizes two, four, and sixteen.
- Hostile boundaries: `m=2` has a constant-zero bijection and must be excluded from the positive-degree bound; `m=4` attains degree two, so a blanket quadratic obstruction from all even `m>2` is false; `m=6` is a prime quadratic positive control.

Anchor: minimum Boolean degree. Niche: dyadic fixed-subtorsor inheritance. Wildcard: Boolean weight divisibility and multiplication depth. The concept board is the complement quotient; free dyadic source orbits; complement-paired target types; Walsh divisibility; and fixed-subtorsor restriction.

## 2. Exact coordinate conventions and complete scalar parametrization

Identify `X_m` with `F_2^m`, using coordinates `x_0,...,x_(m-1)`, and define

    T(x)_0 = 1+x_(m-1),       T(x)_i=x_(i-1) for i>0,
    R(y)_0 = y_(m-1),         R(y)_i=y_(i-1) for i>0.

Thus `T^m x=x+1` and `T^(2m)=1`. A full antiperiodic word is the first `m` coordinates followed by their complements. Define

    z_i=x_i+x_(m-1),     0<=i<m-1.

**Lemma 1 (all equivariant maps).** Every equivariant map `F:X_m -> F_2^m` is uniquely determined by a complement-invariant scalar Boolean function `f`, with

    F_i(x)=f(T^(-i)x),        f(x+1)=f(x).

The function `f` is an arbitrary Boolean function of `z_0,...,z_(m-2)`. Its Boolean degree equals the maximum degree of the coordinates of `F`, and equals its degree as a function of the `z_i`.

**Proof.** Equivariance determines every coordinate from the zeroth. The relation `T^m x=x+1`, together with `R^m=1`, forces complement invariance. Conversely that invariance makes the displayed formula equivariant, including the wrap from coordinate zero to coordinate `m-1`. The change of coordinates from `x` to `(z,x_(m-1))` is invertible linear. Complementation fixes `z` and toggles the last coordinate, so the unique squarefree Boolean polynomial is independent of the last coordinate. Invertible affine substitutions preserve algebraic degree, and every `T^i` is affine. ∎

Consequently, the complete universe of degree-at-most-`d` equivariant rules has size

    2^[ sum_(j=0,...,min(d,m-1)) binomial(m-1,j) ].

At `m=8,d=2`, this is exactly `2^29=536,870,912`, not a restricted local-rule sample.

## 3. A forced Walsh defect at dyadic length

For a Boolean function `f:F_2^r -> F_2`, write

    W(f)=sum_z (-1)^f(z)=2^r-2*weight(f).

**Lemma 2 (dyadic orbit count).** If `m=2^a`, every source orbit in `X_m` has length `2m`, and the number of target types is

    A(m)=2^(m-1)/m.

**Proof.** Every proper orbit length of the cyclic group of order `2m` divides `m`. Such an orbit would be fixed by `T^m`, but `T^m` complements every coordinate and has no fixed points. Hence the source count is `2^m/(2m)`.

For a self-contained target count, if `d=2^b>1`, the number of primitive length-`d` even-weight words is `2^(d-1)-2^(d/2)`: all shorter-period words are repetitions from length `d/2` and have even total weight. There is one good period-one word, zero. Therefore

    #target types = 1 + sum_(b=1,...,a) [2^(2^b-1)-2^(2^(b-1))]/2^b
                  = 2^(m-1)/m,

by telescoping. ∎

**Lemma 3 (forced scalar weight).** Suppose `m=2^a>=4` and `F` is an equivariant type bijection. For its scalar `f`, viewed on the `m-1` complement-quotient variables,

    weight(f)=2^(m-2)-m/2,       W(f)=m.

**Proof.** Every nonzero target type has primitive period `d>1`, hence even `d`. Complementation preserves its primitive period and even primitive-root weight. Thus complementation is an involution on the nonzero target types. Two distinct paired types have length-`m` weights summing to `m`; a self-complementary type has weight `m/2`. Zero is the sole target type whose complement is excluded. It follows that the sum of the weights of one representative of every target type is

    S=(m/2)*(A(m)-1).

An entire source orbit has `2m` members. Equivariance makes the zeroth coordinate read every output coordinate exactly twice over that source orbit, regardless of the output's smaller primitive period. Consequently the weight of the zeroth coordinate over all `2^m` source words is `2S`.

Each complement-quotient point has exactly two source representatives, so the scalar weight on `F_2^(m-1)` is `S=2^(m-2)-m/2`. The displayed Walsh value follows. ∎

The self-complementary target types cause no exception. Nor may one assume all output orbits have length `m`; the twice-per-coordinate calculation deliberately avoids that mistake.

## 4. Elementary degree divisibility

**Lemma 4 (Boolean Walsh divisibility).** If `f` is a Boolean polynomial of degree at most `d`, where `1<=d<=r`, then

    2^ceil(r/d) divides W(f).

**Proof.** Write the algebraic normal form as `c + sum_(j in J) M_j`, where every nonconstant squarefree monomial `M_j` has support of size at most `d`. Pointwise, as an identity of integers on Boolean inputs,

    (-1)^f = (-1)^c * product_(j in J) (1-2M_j).

Choose a term using `k` nonconstant factors. If their union of supports has size `u`, its sum over all Boolean inputs is a signed power

    2^(k+r-u).

Here `u<=min(r,dk)`, so

    k+r-u >= max(k,r-(d-1)k) >= ceil(r/d).

For the last inequality, either `k>=ceil(r/d)`, or `k<r/d`, in which case `r-(d-1)k>r/d`; the exponent is an integer. The empty-factor term is `2^r` and also has the required divisibility. Every term in the expansion is therefore divisible by `2^ceil(r/d)`, as is their sum. ∎

This proof is included in full. No external divisibility theorem or floating-point calculation is being assumed.

**Theorem 5 (dyadic degree bound).** For `m=2^a>=4`, every equivariant type bijection satisfies

    degree(F) >= ceil((m-1)/a).

**Proof.** There are at least two types, so `F` is nonconstant and its degree `d` is positive. Lemmas 1 and 3 give a degree-`d` scalar on `r=m-1` variables with `W(f)=m=2^a`. Lemma 4 forces `ceil((m-1)/d)<=a`, which is equivalent to the claimed lower bound. ∎

In particular degree two is impossible for every dyadic `m>=8`. The full all-degree result is stronger than a quadratic-rank argument and uses the same exact missing-complement defect.

## 5. Restriction to the dyadic fixed subtorsor

**Theorem 6 (odd-repeat restriction).** Let `m=e q`, with `q` odd. Every equivariant type bijection at length `m` restricts to an equivariant type bijection at length `e`; the restricted Boolean degree is no larger.

**Proof.** Work with full antiperiodic words of length `2m`. The source fixed subset

    Fix(T^(2e)) subset X_m

consists exactly of repetitions of antiperiodic length-`2e` words. Indeed, after imposing period `2e`, the antiperiodic shift by `m=e q` is the shift by `e`, because `q` is odd. Thus this fixed subset is an affine embedded copy of `X_e`, and the restriction of `T` is its usual twisted shift. The embedding repeats or complements coordinates, so it does not increase Boolean degree.

If `x` lies in this subset, equivariance gives `R^(2e)F(x)=F(x)`. Because `R^m=1` and `gcd(2e,m)=e`, the output has period dividing `e`. These output words are precisely repetitions of length-`e` words, and their membership in `B_m` is equivalent to membership of the repeated block in `B_e`: the primitive word itself is unchanged.

The source embedding preserves and reflects cyclic orbit equivalence. The same is true of the target embedding. Global injectivity on types therefore gives an injection from the source types at length `e` to the target types at length `e`. Those finite sets have equal cardinality, by the inherited all-length necklace bijection, or by Lemma 2 when `e` is dyadic. Thus the restricted injection is a bijection. Taking its first `e` output coordinates cannot raise degree. ∎

**Corollary 7 (general dyadic-part bound).** For `m=2^a q`, with `q` odd and `a>=2`, every equivariant type bijection has

    degree(F) >= ceil((2^a-1)/a).

The non-dyadic coordinates cannot conceal the obstruction: the dyadic fixed subtorsor retains it.

## 6. Unbounded necessary degree over actual primes

**Corollary 8 (no bounded-degree prime family).** There is no fixed integer `D` such that, for every odd prime `p`, one can choose a multiplier-equivariant type bijection for prime circulant tournaments and sign-even circulant graphs whose Boolean degree is at most `D`.

**Proof.** For each integer `n>=2`, choose any prime divisor `p` of the Fermat number

    2^(2^n)+1.

Then `2^(2^n)=-1 mod p`, while `2^(2^(n+1))=1 mod p`. Thus the multiplicative order of `2` modulo `p` is exactly `2^(n+1)`. It divides `p-1`, so for `m=(p-1)/2`, its dyadic part is `e=2^a` with `a>=n`.

The function `(2^a-1)/a` is increasing for positive integers `a`: the numerator of the difference at successive indices is `(a-1)2^a+1>0`. Corollary 7 consequently gives

    degree(F) >= ceil((2^n-1)/n),

which is unbounded as `n` grows. The inherited prime-circulant dictionary identifies the cyclic types with actual isomorphism types, completing the prime interpretation. ∎

The argument uses **prime divisors of Fermat numbers**, not a claim that the Fermat numbers themselves are prime. It neither assumes nor proves infinitely many Fermat primes. The degree bound permits a different primitive root, different coefficient choices, and a completely different rule at each prime.

**Scoped circuit consequence.** If a rule is evaluated by a Boolean arithmetic circuit with constants, XOR gates, and binary AND gates, and the maximum number of AND gates on an input-to-output path is `h`, then its Boolean degree is at most `2^h`. Thus its multiplicative depth is at least

    ceil(log_2 ceil((e-1)/log_2(e))).

This is a consequence only for this gate model. It is not a general first-order nondefinability theorem, a lower bound for arbitrary choice of primitives, or a claimed sharp general circuit-complexity theorem.

## 7. Sharp small witnesses

All polynomials below are over `F_2`, with the `z_i` of §2. Set `F_i(x)=f(T^(-i)x)`.

### m=6, p=13: exact minimum degree two

Take

    f = z_0 + z_1 + z_0*z_2 + z_0*z_3.

Represent words as integers `sum x_i 2^i`. Its source representatives and labelled images are

    sources: 0, 2, 4, 6, 10, 12
    images: 40,27,6,53,15,0.

The images meet each of the six good target types once. Full verification on all 64 source words proves equivariance, target landing, labelled surjectivity, type bijectivity, and degree two. The affine obstruction proves optimality.

### m=8, p=17: exact minimum degree three

Take the seven-monomial polynomial

    f = z_0 + z_0*z_4 + z_1*z_5 + z_3*z_4 + z_4*z_5
        + z_1*z_2*z_3 + z_3*z_5*z_6.

Its source representatives and labelled images are

    sources: 0,2,4,6,8,10,12,14,18,20,22,24,26,36,42,44
    images:  226,210,130,198,48,0,154,207,232,190,222,153,201,240,9,172.

The target type representatives are

    0,3,5,9,15,23,27,29,39,43,45,51,53,63,95,111.

Full verification on all 256 source words proves every required property. Its labelled output set has 117 members: 112 have fibers of size two, four have fibers of size four, and zero has a fiber of size sixteen. The scalar has weight 120 on the full input space, weight 60 on the complement quotient, and quotient Walsh value eight. Theorem 5 excludes degree two, so the minimum degree is exactly three.

### m=4: necklace-only boundary

The polynomial `f=z_1+z_0*z_1` maps source representatives `0,2` to `0,3` and gives a quadratic type bijection. Here `2m+1=9` is composite, so no prime graph-isomorphism interpretation is asserted. This is a sharp positive control for the dyadic lower bound.

## 8. Independent computation and audit

The original degree-at-most-two enumeration gives:

| m | Polynomial dimension | Complete map universe | Target-valued rules | Type bijections |
|---|---:|---:|---:|---:|
| 4 | 7 | 128 | 25 | 8 |
| 6 | 16 | 65,536 | 2,368 | 432 |
| 8 | 29 | 536,870,912 | 22,136,841 | 0 |

The `m=8` target-valued total is the sum of the published coverage histogram; its sharp maximum coverage is fifteen of sixteen types, attained by 24,576 quadratic rules. The exhaustive search is corroboration, not a dependency of the impossibility proof.

Independent audit paths:

1. `audit_basis.py` independently reconstructs twisted source orbits, target primitive periods, complement-quotient monomials, and every bit of the input table used by the complete C++ enumeration. It reproduces `input8.txt` exactly.
2. The C++ enumeration was recompiled and rerun over all `2^29` coefficient vectors. Gray-code traversal toggles the basis indexed by the least significant set bit of the counter; this visits each coefficient vector exactly once. Its output matched `exhaust8.json` byte-for-byte.
3. An independent direct-coefficient traversal, not Gray-code traversal, matches all small-map counts through `m=6`; see `small_enumeration.audit.json`.
4. `verify_degree_boundary.py` reconstructs the witnesses through direct repeated forward shifts rather than the search's inverse-bit-shift implementation. It checks all source words, computes algebraic normal forms independently by Möbius inversion in both full and quotient variables, and verifies the fiber statistics.
5. The same verifier checks the Walsh divisibility statement for every Boolean function in up to four variables; the exact dyadic target count and weight identity at `m=4,8,16`; and odd-repeat source/target embeddings at `(e,q)=(4,3),(4,5),(8,3),(8,5)`.
6. Normal and optimized Python runs gave identical receipts, with 268,536 explicit checks. The checker uses no removable `assert` statements.

The sparse cubic was found by deterministic seeded search in the degree-three subspace satisfying the linear even-output trace constraints. Discovery was a search, not an exhaustive cubic census. The displayed full witness is independently exact; its minimal **degree** is proved, but minimal **monomial count** is not claimed.

Reproduce from the directory containing this note (all Python data paths are anchored to the script directory):

    python3 audit_basis.py
    python3 audit_small_enumeration.py
    g++ -O3 -std=c++17 exhaust8.cpp -o /tmp/exhaust8
    /tmp/exhaust8 input8.txt > exhaust8.audit.json
    cmp exhaust8.json exhaust8.audit.json
    python3 verify_degree_boundary.py --output verification.json
    python3 -O verify_degree_boundary.py --output verification.optimized.json
    cmp verification.json verification.optimized.json

Optional deterministic witness-discovery reproduction:

    python3 build_search_inputs.py
    g++ -O3 -std=c++17 search_sparse_cubic8.cpp -o /tmp/search_sparse_cubic8
    /tmp/search_sparse_cubic8 1000000000 cubic8.txt > cubic8.sparse.json

The search is unnecessary for validating the explicit seven-term witness; the standalone Python verifier does that directly.

## 9. Connection contract and open next questions

- Source: a low-degree Boolean scalar on the complement quotient of an antiperiodic selector torsor.
- Map: form its equivariant vector rule; restrict to the dyadic fixed subtorsor; sum the scalar sign over that quotient.
- Preserved predicates: cyclic equivariance, injectivity on source types, and nonincrease of Boolean degree under affine restriction.
- Lost information: the scalar Walsh sum forgets almost the entire type permutation and all output phases.
- Restoring sidecar: exact target complement pairing and free source orbit length force the one Walsh value needed for the obstruction.
- Decisive test: the missing complement of the zero target creates `W=e`, whose dyadic valuation is too small for bounded degree as `e` grows.

The sharp remaining question is whether the lower bound `ceil((e-1)/log_2 e)` is attained at further dyadic lengths, beginning at the necklace-only case `e=16` (lower bound four), or whether additional orbit-matching constraints force higher degree. For prime lengths, the next low-bound cases with `v_2(m)=1`, beyond the exact `p=13` witness, also remain open here.

No composite-order or prime-power graph-isomorphism extension is inferred. The universal bound is a theorem about these cyclic selector/necklace objects; its actual circulant graph/tournament interpretation is asserted only at prime `p=2m+1`. It does not obstruct the already-proved unrestricted nonlinear bijections, and it does not settle arbitrary first-order or generator-free definability.
