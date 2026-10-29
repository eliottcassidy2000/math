# Audit and integration of concurrent arithmetic extensions

2026-10-29. **PROVED elementary deductions / inherited group-order input / FINITE-EXACT controls.** This note credits and audits the additional mechanisms in [concurrent PROOFS.md](https://github.com/eliottcassidy2000/math/blob/6ce7838d11e510708b1ab9ed0de418fc4e281c49/05-knowledge/results/quasirandom_arithmetic_structure_20261029/PROOFS.md), read at immutable commit `6ce7838d11e510708b1ab9ed0de418fc4e281c49`, blob `3605deff872fdb5b87cad6f171f0815d253d0796`. The source credits Session107 checkpoint `966ab75d90` for the independently shared colored law and its stronger arbitrary-fiber bound. Conversely, the concrete additions below are credited to that concurrent packet, not claimed as discoveries of this audit. Independent structural audit PASS; see `../audit_structure/AUDIT.md`.

## 1. Two observers on the six-label bank

Direct arguments (4,2,5,3,7,13) have squareclass rank5 and32 colors, as our first checkpoint states. Applying sigma FIRST gives arguments (7,3,6,4,8,14), whose squareclasses use only2,3,7 and have rank3. Thus the sigma observer has8 colors. These are different maps, not contradictory rank computations.

Index bits c_b by the original labels b. The complete relations are

    c_3=0, c_2+c_5+c_7=0, c_4+c_7+c_13=0.

They follow by multiplying the corresponding arguments:4, 3*6*8=144, and7*8*14=784 are squares. The three relations are independent, and columns for b=2,4,7 already span the three prime coordinates, so there are no additional relations. The24 reduced classes q=3 mod4 modulo168 yield each color exactly three times. Fixed-progression PNT therefore gives limiting density1/8. sigma(11)=12 has the same squareclass as sigma(2)=3, so c_11=c_2. The incoming word “free” here means **zero additional storage**, not an independent free bit.

Three stored bits are necessary and sufficient to answer all static named queries: eight different complete answers occur. This says nothing about sufficiency for future arithmetic updates. Noncoprime atoms remain nonmultiplicative under sigma; in particular this character calculation does not legalize treating sigma(2*4) as sigma(2)sigma(4).

## 2. Hamming kernel and simplex character image

For the seven-prime bank (2,5,7,19,29,59,89), the squareclasses of p+1 are (3,6,2,5,30,15,10). As columns in the basis(2,3,5), these are precisely the seven nonzero vectors in F2^3. Let M be this3-by7 matrix.

A selection of columns has square product exactly when its incidence vector is in ker M. Rank3 gives dimension4; no column is zero and no two columns coincide, so no kernel vector has weight1 or2. There are weight3 vectors because each pair of distinct columns has a third column equal to its sum. This proves the [7,4,3] Hamming parameters without appealing to the code's name.

The character response space is image(M transpose), dimension3. For each nonzero z, the nonzero linear functional z dot(-) is1 on exactly four of the seven nonzero vectors. Thus every nonzero response has weight4, giving the [7,3,4] simplex parameters. Inclusion in (ker M)^perp follows from z dot Mw=0; equal dimensions give equality. This is the exact duality between square-product selections and character evaluation.

Distinct coordinate columns are linearly independent, so each pair of response bits takes all four values equally often. A triple whose columns sum to zero has parity zero identically. This realizes pairwise independence with exact triple dependence; it does not contradict quasirandom moving-vertex edges. The fixed-color theorem proved independently in both packets transports this constrained color law together with uniform finite edge patterns. No Fano-to-ordinary-graph-cycle identification is made.

## 3. All odd Sylvester terms, not only the prime subsequence

Define s_1=2, s_(j+1)=s_j^2-s_j+1. Inductively s_(j+1)-1=product_(i<=j)s_i. Hence terms are pairwise coprime and s_j=1 modulo s_i when i<j. Starting from s_2=3, the recurrence preserves residue3 modulo4 for every subsequent term.

For2<=i<j, the Jacobi symbol (s_j/s_i)=1 because its numerator is1 modulo the positive odd denominator. General Jacobi reciprocity applies to these coprime odd integers without primality or squarefreeness and gives (s_i/s_j)=-1. Thus the Jacobi-sign orientation is a legitimate tournament on **all** odd terms and is the reverse chronological linear order. It is an unconditional infinite transitive family. Unlike our earlier prime-subsequence statement, it needs no assumption about infinitely many prime terms. This strengthens that statement rather than refuting it.

For n selected consecutive odd terms split into earlier and later halves A,B, the signed cut has magnitude |A||B|. Testing the matrix on their indicator vectors gives operator norm at least sqrt(|A||B|), linear in n. Cyclic triangles have density0, rather than the quasirandom1/4. The recurrence enforces correlations by permanently setting future residues to1; it is not an averaging mechanism. A positive Jacobi symbol modulo a composite is **not** an iff criterion for being a square. That stronger interpretation is nowhere needed here.

## 4. A genuine composite-base cycle

Use the inherited audited Monster/Baby Monster order ratio R, with prime exponents

    v_2=5, v_3=7, v_5=3, v_7=4, v_11=1, v_13=2,
    v_29=v_41=v_59=v_71=1,

and all others zero. This audit checks deductions from that vector; source group-order verification is delegated to the group's existing audit.

For an integer b>=2 with prime factorization product p^a_p, define

    e_R(b)=min_(p|b) floor(v_p(R)/a_p).

This equals max{j:b^j divides R}, by comparing each prime exponent. In particular e_R(4)=floor(5/2)=2. The exact functional graph therefore contains the genuine cycle

    4 -> 2 -> 5 -> 3 -> 7 -> 4,

with13->2 and11->1. This uses composite-base divisibility, not prime valuation at “prime4,” and does not identify4 with13. At1 the maximum is not finite (1^j divides R for all j), so1 is explicitly terminal rather than part of the function domain.

For R/2 only v_2 changes, from5 to4. The result is4<->2, with5->3->7->4 and13->2 feeding it. This is a genuine new map resolving the earlier raw-six-cycle mismatch, not a claim that the previous prime-valuation map had such a cycle. Neither map is a group action. The sigma observer above is yet another separate operation on the same labels.

## 5. Fermat-prime rigidity and the two operators at3511

If p=2^(2^m)+1 is prime, then 2^(2^m)=-1 modulo p, and d=2^(m+1) is the exact order of2: the order divides this power of2 but does not divide its half. Moreover

    2^d-1=(p-1)^2-1=p(p-2).

Thus v_p(2^d-1)=1. Since d divides p-1 and (p-1)/d is not divisible by p, expanding (2^d)^((p-1)/d) modulo p^2 shows v_p(2^(p-1)-1)=1 too. Hence **no Fermat prime is Wieferich**, including any not currently known. The exact algebra proves the universal implication; testing the five familiar primes is only a guard. The concurrent packet credits this classical exclusion to Crandall–Dilcher–Pomerance; we make no novelty claim or independent historical-priority claim.

In contrast, the finite certificates ord_3511(2)=1755 and v_3511(2^1755-1)=2 show a delayed prime-power order lift. As1755=(3511-1)/2 is odd, the subgroup generated by2 is exactly the nonzero squares. Multiplication by2 therefore preserves the Paley tournament on F_3511.

The Paley signed matrix K_xy=chi(y-x) satisfies K1=0 and KK transpose=qI-J for every prime q=3 mod4. To check the only nontrivial identity, for x!=z translate and scale the correlation to sum_t chi(t(t-1)). The equation y^2=t(t-1) has q plus that sum solutions. The invertible change of variables (u,v)=(2t-1-2y,2t-1+2y) identifies its solutions with uv=1, of which there are q-1. Hence the correlation is-1. Diagonal entries are q-1, yielding the stated Gram matrix and norm sqrt(q).

This exact additive-difference cancellation operator is not the permutation operator multiplying units modulo q^k by2. The former is optimally mixing in its spectral scale while the latter has delayed cycle growth at3511; no contradiction results. At1093, the exact congruence2^182=-1 modulo1093^2 precludes a doubling-invariant tournament connection set: invariance would imply S=-S, incompatible with choosing exactly one of each nonzero pair{x,-x}. Thus the two Wieferich examples also differ in their allowable orientation symmetries.

## 6. Verification scope

Run `python3 incoming/check_concurrent_extensions.py` from the packet root. It writes `concurrent_extensions_receipt.json`. It uses explicit `require` guards (also active under Python optimization), exhaustive finite binary codes and CRT residue classes, nine Sylvester terms including composite1807, both composite-base exponent maps, exact Wieferich valuations, five Fermat-prime examples, and Paley Gram matrices at3,7,11,19,43. The proofs, not these bounded checks, establish all-size assertions.

This scoped audit does not re-audit the incoming nonlinear19 all-level clock, the full automorphism classification of arbitrary circulant tournaments, or representation-degree tables. Those independent claims are not prerequisites of the additions proved here.

## 7. Compact colored A000568 bridge (parent-session deduction)

Let q=2^r be the number of feasible fixed character colors, treated as **named** colors: we quotient by vertex relabelling, not by any permutation or automorphism of the color code. For a partition lambda=(l_1,...,l_t) of n into odd parts, define

    e(lambda)=sum_i (l_i-1)/2 + sum_(i<j) gcd(l_i,l_j),
    z_lambda=product_j j^(m_j) m_j!,

where m_j is the multiplicity of part j. The number of colored tournament isomorphism classes is

    A_n(q)=sum_(lambda odd partition of n) 2^e(lambda) q^t/z_lambda.   (E)

**Proof.** Burnside averages fixed objects over vertex permutations. If a vertex cycle has even length l, its l/2 power swaps opposite vertices, so no tournament is fixed. For an all-odd-cycle permutation, edge orbits cannot contain an edge and its reverse: the permutation has odd order. Inside a cycle of length l there are (l-1)/2 unordered-pair orbits. Between cycles of lengths a,b there are gcd(a,b) orbits. Each pair orbit permits two consistent orientations, giving2^e. A fixed vertex coloring must be constant on each vertex cycle, giving q^t possibilities independently. There are n!/z_lambda permutations of each cycle type, proving(E).

In particular

    A_3(q)=(4q^3+2q)/3,    A_4(q)=(8q^4+4q^2)/3.

Exact enumeration and an independent permutation-average calculation agree for n=1,...,4 and q=1,2. Their counts are respectively(1,1,2,4) and(2,4,12,48).

For any fixed colored isomorphism class H on n vertices, the colored prime-prefix theorem gives limiting induced-subset frequency

    n!/|Aut_color(H)| * 2^(-binom(n,2)) * q^(-n) > 0.

Here Aut_color preserves every named color. All named feasible color assignments occur, so every class counted by(E) is realized with positive asymptotic frequency. This is a finite-pattern statement, not uniform in growing n or color rank. Formula(E) transfers the ordinary tournament Burnside sum by adding a cycle weight q; it supplies **no faster enumeration algorithm** and does not quotient by the simplex/Hamming code's own symmetry group. This compact transfer was proposed by the parent session after the concurrent-source audit and is not attributed to the incoming source.
