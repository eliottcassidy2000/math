# Labelled square-magic densities: a field asymptotic and ring-content bridge

Status: PROVED from the cited Weil bound and the pinned session55 ring-count theorem; finite controls are supporting computations, not an asymptotic proof. No external novelty or first-asymptotic claim is made. This packet does not modify its source packets.

## Scope, inheritance and the connection contract

Count labelled 3-by-3 magic arrays whose nine entries are pairwise distinct square VALUES, with zero permitted. Root choices and dihedral orbit quotients are not counted. The common parameterization is

    c+u      c-u-v     c+v
    c-u+v    c         c+u-v
    c-v      c+u+v     c-u.

Every magic array has this unique form over every commutative ring, by the division-free argument in the pinned source proofs. The parameter space has dimension three.

Inheritance: the audited six-distinct-linear-factor slice (field source, Sections 2–4) and square-content descent (ring source, Sections 4–5). Anchor: replace an existence estimate with a uniform field asymptotic. Niche: expose the density effect of the two-step descent. Wildcard/hostile: distinguish the numerical 63 factor from the characteristic-two field F64, where distinct magic arrays are impossible. The corrected near miss is an independent-residue heuristic: the nine square constraints are correlated and must be bounded, not assumed independent. The retained sidecar is zero-center and collision data.

Concept board: labelled magic plane; quadratic-character slices; exceptional hyperplanes; square-content valuation; normalized densities; field/ring object type. Each participates below: slices count the plane, exceptional hyperplanes have lower dimension, valuation yields a recurrence, and normalization compares constants without identifying fields and rings.

Connection contract: square scaling transports nonzero-center FIELD arrays to center one and preserves all predicates, while forgetting the scaling factor unless retained. Square-content division transports all-nonunit RING arrays at level k bijectively to arrays at level k−2, preserving squares, magic equalities and distinctness, with the modulus as required sidecar. There is NO map here from F64 objects to these ring arrays, nor any MRB/modular-CM assertion, integer-square conclusion, or priority assertion.

## 1. Uniform finite-field asymptotic

THEOREM. For every odd prime power q, let F(q) be the number of arrays in the stated universe over the field F_q. Uniformly over odd characteristics and extension degrees,

    F(q) = q^3/512 + O(q^(5/2)).

In particular F(q)/q^3 tends to 1/512 as q tends to infinity through odd prime powers. An explicit conservative bound valid for every odd q is

    |F(q) − q^3/512| <= q^(5/2) + 14q^2.

Proof. Let chi be the quadratic character, extended by chi(0)=0. Let A be the set of v for which both 1+v and 1−v are nonzero squares. For

    W(v)=(1+chi(1+v))(1+chi(1−v))/4,

linear character sums vanish and the exact two-root quadratic identity gives

    sum_v W(v) = (q−chi(−1))/4.

Away from v=±1, W is the indicator of A. At those two points its total is (1+chi(2))/2. Thus, exactly,

    |A| = (q−chi(−1))/4 − (1+chi(2))/2,
    ||A|−q/4| <= 5/4.

Set E={0,±1,±2} and G=A\E. Then ||G|−q/4|<=25/4, uniformly even when E has repeated values in small characteristic.

Fix v in G. The six remaining entries at c=1 are

    1+u, 1−u, 1+u+v, 1+u−v, 1−u+v, 1−u−v.

Their roots ±1, ±(1+v), ±(1−v) are distinct. Expand 2^(-6) product_L(1+chi(L(u))). The empty subset contributes q; singleton sums vanish; the fifteen two-factor sums each have absolute value one. The Weil bound for squarefree products of d distinct linear factors gives error at most (d−1)sqrt(q). The aggregate for d>=3 is

    [20*2 + 15*3 + 6*4 + 1*5] sqrt(q) = 114 sqrt(q).

This is the cited standard bound, not an independence assumption. It applies in characteristics 3 and 5 as well; squarefreeness, not degree less than characteristic, is the hypothesis. Deleting the six roots costs at most 3 in weighted count. Distinctness deletes at most seven further u values, namely 0, ±v, ±2v, ±v/2. Let N_v count the remaining labelled arrays. Therefore

    |N_v−q/64| <= (15+114sqrt(q))/64 + 10.

Summing over G and using |G|<=q gives

    |sum_G N_v − q^2/256|
    <= (114/64)q^(3/2) + (655/64+25/256)q.

Let C(q) count ALL valid center-one arrays, allowing zero entries. Arrays with no zero entry outside the sum above must have v in A\G, costing at most 5q. Arrays with some zero entry lie on one of the eight noncentral affine lines c+su+tv=0 at c=1; each contains q pairs (u,v), including in small odd characteristic. Their union costs at most 8q. Consequently

    |C(q)−q^2/256| <= 2q^(3/2)+24q.

Every nonzero square center c has exactly C(q) arrays: multiply by c^(-1), which is a square. There are (q−1)/2 such centers. The zero-center branch has at most q^2 arrays simply because there are only q^2 pairs (u,v); this bound does not assume distinctness or square membership. Thus, writing 0<=Z(q)<=q^2,

    F(q) = ((q−1)/2) C(q) + Z(q).

Subtract q^3/512 and apply the preceding bounds. The error is at most q^(5/2)+(12+1+1/512)q^2, which is at most the displayed q^(5/2)+14q^2. This completes the proof.

The primary character-sum reference inherited from the source is Theorem 4.3 of L. Babai's researcher-authored notes: https://people.cs.uchicago.edu/~laci/HANDOUTS/weil.pdf . The local source gives the elementary exact quadratic identity used above. No fresh literature-priority survey is claimed.

## 2. General square-content recurrence, including the boundary

For ANY prime p and k>=1, write N_p(k) for the same labelled distinct square-value count over Z/p^k Z, and P_p(k) for those counted arrays with AT LEAST ONE unit entry. Put N_p(0)=0 (the one-element ring has no nine distinct entries), and for recurrence notation N_p(−1)=0. Then

    N_p(k)=P_p(k)+N_p(k−2), k>=1.

For k>=2, every nonunit square is divisible by p^2: if x=r^2 and p divides x, then p divides r. Multiplication by p^2 is a bijection from Z/p^(k−2) Z onto the ideal p^2(Z/p^k Z). A value p^2y is a square modulo p^k iff y is a square modulo p^(k−2): one direction squares a lifted root multiplied by p; the other divides a root divisible by p. Applied entrywise, this map preserves distinctness because it is injective onto the ideal, and preserves and reflects magic equalities because p^2(S−T)=0 modulo p^k iff S−T=0 modulo p^(k−2). Hence the all-nonunit part is in bijection with level k−2.

At k=1 the only nonunit is zero, so an all-nonunit array cannot have distinct entries. At k=2 the only nonunit square is again zero, with the same conclusion. Thus N_p(1)=P_p(1) and N_p(2)=P_p(2); these are the correct boundary cases, not invented negative-level rings.

For p>=5, P_p(k) need NOT mean all-unit. Hostile witness modulo 29, in the displayed parameterization with c=0,u=1,v=6, is

    1 22 6
    5  0 24
    23 7 28.

These are nine distinct squares modulo 29, and each line sums to zero. It is primitive in the at-least-one-unit sense but not all-unit. This distinction is essential in extending the recurrence beyond 2 and 3.

## 3. Why 63 and 1728 occur in this comparison

For p=2,3 ONLY, the local square-residue constraints force each grid into the all-unit or all-nonunit branch (pinned ring proof, Section 4). Its primitive part therefore equals its all-unit part. With s_2=3 and s_3=1, the affine maps x -> 1+p^s_p x give

    P_2(k)=A_2(k−3), k>=3;
    P_3(k)=A_3(k−1), k>=1,

where A_p(h) counts ordinary distinct magic arrays at modulus p^h. Low levels are zero as detailed in the source. The exact source formulas imply

    A_p(h)=p^(3h)+O(p^(2h)),
    lim P_2(k)/2^(3k)=2^(−9)=1/512,
    lim P_3(k)/3^(3k)=3^(−3)=1/27.

Iterating the content recurrence and normalizing by p^(3k), each two-exponent drop contributes p^(−6), since the magic parameter space has three coordinates. Summing the geometric series gives

    lim N_p(k)/p^(3k) = p^(−3s_p)/(1−p^(−6)), p in {2,3}.

For rigor, the source's error after summing is bounded by a constant times p^(2k)/p^(3k), and tends to zero; the leading truncated geometric sum converges to the displayed expression. Thus

    lim N_2(k)/2^(3k) = (1/512)*(64/63) = 1/504,
    lim N_3(k)/3^(3k) = (1/27)*(729/728) = 27/728.

Relative to the independently proved odd-field density 1/512, the ratios are

    (1/504)/(1/512) = 64/63,
    (27/728)/(1/512) = 1728/91.

The first factor 63=2^6−1 comes directly from square-content renewal in dimension three. The second arithmetic simplification uses 728=8*91 and 512*27/8=64*27=1728. Equivalently the ternary primitive density differs already by 512/27, and then acquires its separate renewal factor 729/728. The integer 1728 is NOT itself the ternary renewal denominator or an independently proved modular invariant here. These are precise density comparisons, not an identification of the underlying objects.

In particular F64 is a characteristic-two FIELD and has zero distinct magic arrays because opposite entries coincide. The local dyadic densities above are for RINGS Z/2^k Z as k tends to infinity; the factor 64/63 does not turn F64 into their carrier. No MRB, modular-CM, supersingularity, or integer-height connection follows from these constants.

## 4. Reproduction and scope of tests

Run `python3 -B verify.py --write` to regenerate certificate.json, then `python3 -B verify.py` and `python3 -O -B verify.py` to compare exactly. All load-bearing checks use explicit exceptions, not removable asserts. The verifier checks pinned source hashes; the exact admissible-slice identity and two-sided slice bound over all odd prime powers through 125; center counts against the source census; general content recurrence and primitive/all-unit distinction in bounded quotient rings; exact rational limiting constants and ratios; and normalized errors at finite large k. A portable copy of source code is included, with provenance and byte hashes in source-pins.json. The imported field implementation is reused for arithmetic, not advertised as an independent field-engine audit. The all-q and all-k statements are established by the proofs above.
