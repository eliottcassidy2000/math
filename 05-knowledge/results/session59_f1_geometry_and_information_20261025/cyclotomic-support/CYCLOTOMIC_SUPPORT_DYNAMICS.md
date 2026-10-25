# Cyclotomic support frames: exact semilinear fixed points and arithmetic boundaries

Status: **PROVED** by the arguments below; finite tests are separately **FINITE-EXACT**. The cyclic pointed-monoid/projective frame and its elementary count are established literature, not new discoveries. The present reusable contribution is a self-contained fixed-point compiler with diagonal twists, coordinate permutations and power maps; its dynamics/orbit consequences; and explicit guards on what the Session-58 arithmetic index bridge transports. No global hourglass exclusion, F1 descent of the hourglass equations, or literature-priority claim is made.

## 1. Model and the essential base-category guardrail

Let A_m = C_m ∪ {0}, m≥1, where C_m is cyclic of order m, multiplication extends by an absorbing zero, and maps preserve 0 and 1. Write C_m multiplicatively, or choose a generator ζ and identify its exponents with Z/m. The identity unit ζ^0 must never be confused with the absorbing zero. Let X_k be the standard pointed-monoid projective space on k homogeneous generators over A_1. Put

    P_m(k) = Hom_{A_1}(Spec A_m, X_k)
           = ((A_m)^k \ {(0,...,0)}) / C_m,
    N_k(m) = ((m+1)^k−1)/m.

For k=0 both point sets are empty and N_0(m)=0. For k≥1 the second equality is a natural bijection, not just equality of cardinalities. A homogeneous vector with coordinate i nonzero gives chart coordinates x_j/x_i ∈ A_m. Conversely a chart map is determined by those ratios; maps on overlapping charts identify exactly proportional vectors. The one-point source lies in some chart, so these constructions exhaust the Hom set.

Equivalently one can use **base-preserving** Hom_{A_m}(Spec A_m, P^{k−1}_{A_m}). If instead the target is already extended to A_n and morphisms need not preserve a specified base map, there is an extra coefficient homomorphism C_n→C_m. There are gcd(n,m) such homomorphisms (the image of a generator is any element killed by n). The same chart argument therefore gives

    |Hom_{A_1}(Spec A_m, P^{k−1}_{A_n})| = gcd(n,m) N_k(m).

In particular an unrestricted absolute Hom into P_{A_m} carries an extra factor m. Silently switching between relative and absolute Hom would change the answer.

### Sources and known/new distinction

Thas, *Projective Spaces over F1ℓ*, §§2 and 4, explicitly defines the cyclic group with zero, projective frames by proportionality, and the functor-of-points count; §4.4 gives the monomial symmetry group and support orbits. Deitmar, *Schemes over F1* (2006), §1 instead begins with unpointed monoids and F1={1}. We use the pointed convention throughout, not an unnoticed mixture of the two. The chart proof above makes the base condition explicit. Sources were read directly on 2026-10-25:

- [Thas, primary PDF](https://backoffice.biblio.ugent.be/download/8589084/8603796), §§2,4, especially equations (4), (17), (20), (21).
- [Deitmar, primary full text](https://arxiv.org/html/math/0404185v7), §1 and Theorem 1.1; [record](https://arxiv.org/abs/math/0404185).

Subsequent formulas are proved here from the finite model. We have not searched sufficiently to claim they are new to the literature.

## 2. Support fibers and three incompatible counts

For a nonempty S⊆{1,...,k}, the exact-support fiber is C_m^S/C_m. Normalizing one chosen coordinate to 1 identifies it with C_m^{|S|−1}. Hence

    |P_m(k)_S| = m^{|S|−1},
    |P_m(k)| = Σ_{s=1}^k binom(k,s)m^{s−1} = N_k(m).

The last identity is the binomial theorem after multiplication by m. The support projection is the map induced by A_m→A_1, and its section labels every nonzero coordinate by the identity unit. Coordinate-zero slices/deletion commute with these maps.

The underlying projective monoid spectrum has one point for every nonempty support, thus 2^k−1 topological points, independently of m; its closed/rank-zero supports are the k singletons. Literal Hom points at m=1 therefore number 2^k−1. They do not number k.

Over a finite field, the usual projective counting polynomial is

    P_k(q) = Σ_s binom(k,s)(q−1)^{s−1} = 1+q+...+q^{k−1}.

Our Hom count is P_k(m+1). At m=1 it is P_k(2), not P_k(1). Formal q=1 corresponds to formally setting m=0 in the polynomial support expansion, retaining singleton supports only. There is no finite cyclic group C_0 in this construction. Thus:

- literal pointed Hom at A_1: 2^k−1;
- closed/rank-zero points: k;
- formal q=1 evaluation: k.

A choice C_{q−1}≅F_q^× identifies the frame point set with P^{k−1}(F_q). It does not by itself add field addition, incidence beyond coordinate support, or an interpolation of integer gap equations.

## 3. A canonical power action and its complete finite dynamics

For integer r≥0 define F_r by 0↦0 and u↦u^r on units, including the explicit convention F_0(0)=0. It induces a projective map because scaling by c changes the image by c^r. The rule F_r F_s=F_{rs} holds. For k≥2 the induced projective map is an automorphism if and only if gcd(r,m)=1. For k=0 the projective point set is empty, and for k=1 it is a singleton; in both cases every power induces the identity permutation. On the coefficient monoid itself the iff is exact. Thus Aut(C_m)=(Z/m)^× acts canonically, without choosing a generator.

**Power fixed-point theorem.** Put d=gcd(m,r−1). Then

    Fix(F_r) ≅ P_d(k),       #Fix(F_r)=N_k(d).

Here C_d is the unique order-d subgroup of C_m. Proof: in the normalized chart on a fixed support, every relative unit must satisfy u^{r−1}=1. These are precisely the d units of C_d. In homogeneous coordinates the common proportionality scalar need not be 1; normalization is what correctly removes it. On each support the answer is d^{|S|−1}.

**Full tail/period theorem.** Fix r≥1. Let B(t,n) count points whose preperiod is at most t and eventual period divides n. For t≥0,n≥1,

    B(t,n)=N_k(gcd(m, r^t(r^n−1))).                         (1)

Indeed F_r^{t+n}(x)=F_r^t(x), and on normalized phase coordinates this says r^t(r^n−1)a=0 in Z/m. For any map of a finite set the displayed equality is equivalent to the stated tail/period condition. The number E(t,n) of points with exact preperiod t and exact eventual period n is

    E(0,n)=Σ_{d|n} μ(n/d) B(0,d),
    E(t,n)=Σ_{d|n} μ(n/d)[B(t,d)−B(t−1,d)]  (t≥1).         (2)

These follow by Möbius inversion, after taking a difference in the tail variable. In particular the number of actual cycles of length n is E(0,n)/n. This remains valid when r is not invertible, where a Burnside-only analysis would miss tails.

For example m=12,k=3,r=2 gives 183 total points. The (tail,period):point counts are (0,1):7, (0,2):14, (1,1):6, (1,2):30, (2,1):18, (2,2):108. These add to 183. Only the first two classes are periodic.

### Unit-group orbit formula

Let U_m=(Z/m)^×, with U_1 the trivial group and φ(1)=1. Burnside gives

    |P_m(k)/U_m| = (1/φ(m)) Σ_{r∈U_m} N_k(gcd(m,r−1)).    (3)

There is a second divisor form, useful for independent control:

    |P_m(k)/U_m|
      = Σ_{s=1}^k binom(k,s) Σ_{d|m} J_{s−1}(d)/φ(d),     (4)

where J_h(d)=Σ_{e|d}μ(d/e)e^h, including J_0(1)=1 and J_0(d)=0 for d>1. To prove (4), normalize a support and sort vectors in (Z/m)^{s−1} by the order d of the subgroup they generate. There are J_{s−1}(d) such vectors. Their scalar stabilizer is exactly the kernel of U_m→U_d, because their coordinates generate the subgroup. The reduction U_m→U_d is surjective: lift each prime-power unit locally, and choose units on primes dividing m but not d; then use CRT. Every orbit has size φ(d). This proves the divisor formula, without assuming that the whole U_m action is free.

## 4. Main compiler: diagonal twist, coordinate permutation, and power

Choose ζ only to write exponents. Let b∈(Z/m)^k, σ∈S_k, and integer r≥0. Define

    T_{b,σ,r}([x]) = [y],    y_i=ζ^{b_i} x_{σ^{-1}i}^r,

with a zero input coordinate still zero. The map is well-defined on projective classes. In additive exponent notation on present coordinates it is a↦b+rσa, with (σa)_i=a_{σ^{-1}i}. Its exact composition law is

    T_{b,σ,r} T_{c,τ,s}=T_{b+rσc, στ, rs}.                (5)

Diagonal b differing by a constant vector have the same projective action. For gcd(r,m)=1 this supplies the action of

    ((C_m)^k/C_m) ⋊ (S_k × U_m).

We claim this natural monomial-semilinear action, not the full arbitrary permutation group of the point set and not universal preservation of additive equations.

For a σ-cycle C=(i_0,i_1,...,i_{ℓ−1}), indexed so σi_j=i_{j+1}, put i_ℓ=i_0 and define

    g_C = gcd(m,r^ℓ−1),
    A_C = 1+r+...+r^{ℓ−1},
    B_C = Σ_{j=1}^ℓ r^{ℓ−j} b_{i_j}.

All congruences below are integer congruences, so representatives of b do not matter.

**Semilinear fixed-point theorem.** For every k≥1,m≥1,r≥0,

    #Fix(T_{b,σ,r})
      = (1/m) Σ_{t∈Z/m} [ Π_{C cycle of σ}
              (1+g_C·1_{g_C | A_C t−B_C}) − 1 ].          (6)

A support-restricted answer is obtained by retaining the product g_C indicators for any specified nonempty union of cycles. Supports not preserved by σ have no fixed points.

**Proof.** A fixed projective point has a representative with exponents a on a nonempty support S and a uniquely determined scalar exponent t satisfying

    b_i+r a_{σ^{-1}i}=a_i+t.                             (7)

Its support must be a union of σ-cycles. On one selected cycle, iterating the recurrence a_{i_j}=b_{i_j}+r a_{i_{j−1}}−t gives

    (r^ℓ−1)a_{i_0}=A_C t−B_C  (mod m).                  (8)

This scalar congruence has g_C solutions if the displayed divisibility holds and none otherwise. Every solution determines all remaining cycle coordinates by recurrence. This works even when r is not a unit; no inverse of r was used. Coordinates on distinct selected cycles are independent after t is fixed. A cycle may instead be absent, contributing 1. Multiplying the alternatives and subtracting 1 excludes the empty support. Finally each fixed projective point is represented exactly m times in the counted pairs (a,t): shifting a by c shifts t by (r−1)c, and the shift action on nonempty a is free. Divide by m. This proves (6), including integrality and independence of the cycle starting point. ∎

Setting b=0 and σ=id recovers the normalized power theorem. Setting m=1 gives 2^{c(σ)}−1, the nonempty σ-invariant supports. Setting r=1 gives the particularly transparent cycle constraints m | ℓt−Σ_{i∈C}b_i.

For any finite subgroup H of the invertible action, insert (6) in |P_m(k)/H|=(1/|H|)Σ_{h∈H}#Fix(h). Formula (5) computes all powers of an element exactly, and hence its cyclic orbit distribution by Burnside/Möbius inversion. It is not enough in general to know only r and the cycle lengths: the weighted diagonal cycle twists B_C are necessary data.

**Nontrivial verified example.** For m=4,k=3,b=(0,1,2), σ=(1 2 3), r=3, the induced map has order 6 on 31 points. It has one 1-cycle, four 3-cycles, and three 6-cycles, hence eight orbits. Dropping the twist is not part of the theorem; all iterate twists are composed using (5).

## 5. A finite-field descent corollary: totals conceal support twists

Let q≥2,d≥1,m=q^d−1 and suppose σ^d=id. Set T=T_{0,σ,q}. Then

    #Fix(T)=(q^k−1)/(q−1).                               (9)

This statement is purely cyclic and holds for any integer q≥2. When q is a prime power, a choice C_m≅F_{q^d}^× interprets it as a genuine finite-field Frobenius-permutation action.

Proof: T^d is the identity. If T(v)=λv for a representative vector, iteration gives λ^{1+q+...+q^{d−1}}=1. In C_m this is precisely the condition λ∈im(c↦c^{q−1}). Rescale v by c with c^{q−1}=λ^{-1}; then T(cv)=cv. Thus every fixed projective point has a genuinely fixed vector representative. On a permutation cycle of length ℓ|d, a fixed vector is determined by one coordinate satisfying x^{q^ℓ}=x; there are 1+gcd(q^d−1,q^ℓ−1)=q^ℓ choices, including zero. Across cycles there are q^k fixed vectors. Remove the zero vector and divide by q−1 fixed nonzero scalar choices. ∎

The total in (9) forgets the splitting of support strata. For k=2,q=3,d=2,m=8:

- σ=id: four fixed projective points, comprising two singleton-support points and two full-support points;
- σ=(12): also four fixed points, but zero singleton-support points and four full-support points.

A fixed support must be a union of Frobenius-permutation cycles. This is an exact split/nonsplit-support distinction invisible in the ordinary projective point count. It illustrates why an arithmetic quotient or toric model may need permutation/phase/Frobenius data as a sidecar. It does not prove any particular external quotient model; that requires specifying its action and descent.

## 6. Exact Session-58 arithmetic bridge at m=2

The inherited squarefree result fixes an ordered center e=∏_{i=1}^k p_i, distinct p_i≡1 (mod 4), and Gaussian primes π_i over p_i. For j∈{−1,0,1}^k\{0}, put

    z(j)=∏_i π_i^{1+j_i} conjugate(π_i)^{1−j_i},
    d(j)=|2 Re(z(j)) Im(z(j))|.

The Session-58 audited theorem gives a bijection j~−j ↔ D_e, where D_e consists of positive d such that e²±d are both squares. Identifying the units of A_2 with {±1} identifies this index quotient exactly with P_2(k), and reduction modulo 3 also identifies it with P^{k−1}(F_3). On support S it has 2^{|S|−1} points and

    gcd(e²,d(j))=(e/∏_{i∈S}p_i)².

Deleting inactive coordinates gives exact integer square division. This support/division diagram is the arithmetic content; the cyclic construction does not establish a new interpolation in m.

The power action is especially limited here: odd r acts identically, while even r sends [j] to [|j|], retaining only support and selecting the all-positive unit labels. It is a retraction onto 2^k−1 support representatives. Under the fixed Gaussian orientations this becomes an explicit map D_e→D_e, but the selected gap for a support depends on those orientations. Aut(C_2) itself is trivial. Thus nontrivial arithmetic permutations at m=2 come from coordinate signs/permutations, not a nontrivial cyclotomic Galois group.

### A rigorous nonadditivity control, without inventing a Schur triple

At e=65, choosing π_5=2+i and π_13=3+2i gives

    j=(1,0) ↦4056,   (0,1) ↦3000,
      (1,1) ↦3696,   (1,−1)↦2016.

The coordinate sign flip in the second coordinate fixes 4056 and 3000 and interchanges 3696 and 2016. It cannot extend to an additive endomorphism of the subgroup of Z generated by D_65: that subgroup is 24Z (the gcd of 4056, 3000, 3696, and 2016 is 24), whose additive endomorphisms are integer multiplications; fixing a nonzero element forces multiplier 1, contradicting the interchange. The even-power retraction also fixes 3000 but changes 2016 to 3696, so it likewise has no additive extension. Therefore these symmetries are not integer-linear/arithmetic-additive symmetries.

This does **not** exhibit an actual Schur triple in D_65 that is destroyed: D_65 contains no Schur triple at all. Claiming such a witness would be false. The justified statement is that the index action and fixed-point formula do not transport the additive predicate automatically. In particular a hypothetical distinct Schur triple entirely within one support would be collapsed by the even-power map to a single positive gap and would no longer satisfy x+y=z. Existence of such a triple is not asserted.

## 7. Why higher exponents are not larger projective scalar groups

For center e=∏p_i^{a_i}, inherited Gaussian allocation indices run through ∏[-a_i,a_i], modulo **simultaneous sign only**, and the count is

    |D_e|=(∏_i(2a_i+1)−1)/2.

When all a_i=a, replacing the sign quotient by the full cyclic scalar quotient C_{2a} would give ((2a+1)^k−1)/(2a), smaller by a. Already a=2,k=1,e=25 has two gaps: j=1 gives 600 and j=2 gives 336, while P_4(1) has one point. Their p=5 valuations are respectively 2 and 0; the extra scalar quotient destroys this proven arithmetic information.

If q=2a+1 is an odd prime, coordinate reduction modulo q identifies the allocation index set with (F_q^k\{0})/{±1}, not P^{k−1}(F_q). The natural projection to projective space is an a-to-one cover: every nonzero F_q-scalar orbit has 2a vectors and hence a sign classes. At a=1 only, the sign subgroup is already the full scalar group. For general a with q not prime, no integer-residue field identification is claimed. Arbitrary alphabet relabeling cannot repair the lost valuations or supply canonical arithmetic operations.

## 8. Verification, provenance, and limits

Run `python3 verify_cyclotomic.py` and `python3 -O verify_cyclotomic.py`. All gates explicitly raise errors; none relies on `assert`. Only integer arithmetic is used. Raw projective points are constructed by normalizing every nonzero tuple. Formula (6) is implemented independently by scalar congruences and permutation-cycle products. The tests cover:

- m=1,...,8, k=1,...,5: 40 point/support counts and 440 ordinary-power fixed counts;
- all twists b, all permutations, and r=0,...,m for m≤4,k≤3: 3,024 exhaustive semilinear fixed-count cases;
- m≤8,k=4, every coordinate permutation and three deterministic twist choices: 3,168 more cases;
- 2,009 pointwise composition checks and 21 cyclic Burnside orbit checks;
- 468 complete functional graphs at m≤12,k≤3, with 3,200 exact tail-period cells checked using (1)–(2);
- 32 unit-group orbit checks using raw orbits, Burnside, and the Jordan-divisor formula;
- coefficient-homomorphism and twisted-descent controls, reported in `verification.json`;
- exact Gaussian values at e=65 and e=25, with integer-square endpoint tests.

These finite checks validate implementations and hostile examples. The universal statements rest on the proofs, not extrapolation. No general F1 formalism is asserted to share this Hom count. A power map is only a “Frobenius-style” endomorphism until a finite-field realization and arithmetic Frobenius are explicitly supplied as in §5.

Inherited arithmetic source files (read-only):

- `05-knowledge/results/session58_hourglass_complement_geometry_20261025/projective-support-bridge/PROJECTIVE_SUPPORT_BRIDGE.md`;
- `05-knowledge/results/session58_hourglass_complement_geometry_20261025/hourglass-prime-support/RESULT.md`.

This worker performed no git operation. `normal.log`, `optimized.log`, `verification.json`, and source hashes accompany the proof/test packet; see `README.md` for environment-independent reproduction commands.
