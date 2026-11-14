# Growing exact (H,Q2) sources for the full 99-digit bank

Session 326. **Proved; separate literal-profile replay passed.** This supplies an infinite physical source alphabet and a complete target-entry algorithm for its ordinal closure. It does not supply arbitrary profiles or close Chan–Pak.

## 1. A local two-slot equality survives unbounded strong growth

For a binary word w of length n, let T(w) be the transitive chain 0→1→…→n−1 together with apex x, where x→i exactly when w_i=1. Write

    P_w(t)=Σ_{i<j, w_i=1, w_j=0} t^(j−i−1),
    B_w(t)=Σ_{w_i=0}t^i + Σ_{w_i=1}t^(n−1−i).

The exact S322 one-apex correlation formula is

    Q_T(w)(z)=z^(n+1)+z^(n−1)[P_w(1+1/z)+P_w(1−1/z)].       (1)

It is a polynomial identity; evaluation at z=1 includes the usual 0^0=1 polynomial convention. To recall its mechanism, every directed cycle uses x. A cycle with first/last chain vertices i,j chooses an even number of internal vertices precisely when it is odd. The odd-cycle packing identity then contributes

    z^(n−d−1)[(z+1)^d+(z−1)^d], d=j−i−1,

and there cannot be two nontrivial cycles in one packing. Thus (1) is an exact correlation formula, not a fitted recurrence.

Choose the two eleven-bit words

    w_- = 10101010100,
    w_+ = 10111011100.

The changed positions are 3 and 7 (zero-based), a symmetric pair. Hence B_+=B_-. Direct subtraction of the finite inversion polynomials gives

    P_+(t)−P_-(t)=t^5−t^4−t^2+2t−2,
    ΔP(1+x)+ΔP(1−x)=2(4x²−1)(x²+1).                       (2)

For k≥0 wrap either word by k leading ones and k trailing zeros:

    w_{±,k}=1^k w_± 0^k,   T_{±,k}=T(w_{±,k}).

Every such tournament is strong: the apex reaches the first chain vertex, the whole chain reaches the last, and the last returns to the apex. Its order is 12+2k.

**Wrapping lemma.** Put S_k(t)=1+t+…+t^(k−1), with S_0=0. Splitting inversions into old-old, new-old, old-new and new-new gives

    P_{1^k w 0^k}(t)=P_w(t)+B_w(t)S_k(t)+t^n S_k(t)^2.     (3)

Because B_+=B_-, the inversion-polynomial difference (2) is unchanged by wrapping. Therefore the entire optional-slot difference is

    Q_{+,k}(z)−Q_{-,k}(z)
       =2 z^(6+2k)(4−z²)(z²+1).                          (4)

In particular the sources have **exactly equal Q(2)** for every k, and their Hamilton counts differ by exactly 12. Their Q(3) values differ by −72900·9^k, so this is not equality of their full slot states.

## 2. Closed profiles, with unit Hamilton counts

The physical profiles are

    h_-(k)=2048·4^k−1026·2^k+119,
    h_+(k)=2048·4^k−1026·2^k+131,
    q(k) =354294·9^k−242340·6^k+50336·4^k−2052·2^k+2.     (5)

At k=0 these give (H,Q2)=(1141,160240),(1153,160240); at k=1 they give (6259,1931848),(6271,1931848).

These closed forms follow by substituting (3) into (1), not by interpolation. For clarity, the shared boundary data are

    B(2)=3070,
    2^11 B(3/2)=466248,
    2^11 B(1/2)=2048.

At z=1, the boundary under one further wrap obeys B_{k+1}(2)=2B_k(2)+2^(13+2k), while H_{k+1}−H_k=2^(11+2k)+B_k(2). Summing gives the two h formulas. At z=2, use S_k(3/2)=2((3/2)^k−1) and S_k(1/2)=2(1−2^(−k)) in (3); collecting 9^k,6^k,4^k,2^k,1 gives q in (5).

Every h_+(k) is coprime to 99. Indeed h_+≡1 modulo 3. Modulo 11 put x=2^k: h_+≡2x²−3x+10, whose discriminant is 6, a nonsquare modulo 11. For h_-, the corresponding polynomial 2x²−3x+9 has roots x=2,5. Since 2 has order ten modulo 11,

    gcd(h_-(k),99)=1  iff  k is not 1 or 4 modulo 10.       (6)

The two h values are coprime to each other for every k: their difference is 12, both are odd, and neither is divisible by 3.

## 3. An infinite physical alphabet and full-bank output

Take as letters all T_{+,k}, k≥0, and all T_{-,k} with k outside the two forbidden classes in (6). Form any finite ordinal multiset of these letters in a fixed canonical order, then add t≥5 singleton SCCs. The multiset may be empty. If its letters have parameters (h_i,q_i), the resulting literal source has

    H=∏h_i,    Q=2^t ∏q_i,    gcd(H,99)=1.                 (7)

Multiplicativity follows by concatenating each path slot through the totally ordered SCCs. Every letter is an explicit strong graph, not a number asserted to be attainable.

The S325 full bank now strongly realizes, for every odd c in 1,…,197,

    N=99·2^t ∏q_i + c∏h_i.                                (8)

The phase budget is paid: the c-template has shift r_c≤4, so remove r_c of the t singleton SCCs before substitution. Its coefficient 99·2^r_c then restores exactly (8). At least one singleton remains; the stipulated nonempty source condition holds. Strip the template's fixed singleton SCCs to retain its strong core, exactly as certified by S325.

This is an unbounded-H profile family. For any allowed fixed k and ℓ≥1 copies, varying a of the copies to the minus source gives ℓ+1 different H values

    h_-(k)^a h_+(k)^(ℓ−a),  0≤a≤ℓ,

at one exact Q=2^t q(k)^ℓ. The 99 digits give 99(ℓ+1) distinct targets at that Q. To see the latter, equality of two offsets with unequal a would force h_- or h_+ to divide a digit ≤197, impossible because the two bases are coprime and exceed 197. This statement uses k allowed for both letters; the plus family alone is allowed at every k.

Concrete target entries from a single plus letter and t=5 include:

| k | digit c | exact strong target | certificate |
|---:|---:|---:|---|
| 0 | 17 | 507659921 | prime |
| 0 | 41 | 507687593 | 7·72526799, prime cofactor |
| 1 | 19 | 6120213613 | prime |
| 1 | 127 | 6120890881 | 7·874412983, prime cofactor |

The local replay verifies the finite primality statements and the corresponding source profile and phase data. The original S325 replay independently validates the literal template profiles and their strong-core construction. It does not infer infinitely many primes on any one subfamily.

## 4. A complete target-first inverse, without factoring

Given an odd target N:

1. Enumerate canonical multisets of the above alphabet with core product q=∏q_i≤N/(99·32). Include the empty multiset (h,q)=(1,1). Only finitely many k occur below this bound, so generate their explicit profiles from (5).
2. At each multiset put h=∏h_i. There is a unique odd c in {1,3,…,197} with ch≡N modulo 198, because gcd(h,198)=1.
3. Test whether

       R=(N−ch)/(99q)

   is an integral power of two at least 32. If so its exponent is the required t. Assemble the multiset, add t singleton SCCs, remove the selected bank phase r_c, and use the literal c-template. Otherwise reject this multiset.
4. Reject the target in this grammar only if every bounded multiset fails.

This is necessary and sufficient for the entire named infinite-alphabet grammar, including arbitrary repetitions, not a finite seed subbank. Every accepted profile is physically supplied. Conversely an actual presentation (8) has q below the bound and appears as a multiset; its unique c and exact t are forced by the target. The inverse does not replace exact Q by a residue, and does not use integer factorization.

## 5. Subpolynomial work and the exact remaining limitation

Equation (3) has nonnegative added inversion coefficients. Thus q(k+1)≥4q(k), and q(0)=160240≥4, so q(k)≥4^(k+1). A multiset below target therefore obeys

    Σ_letters (k+1) ≤ (log₂N)/2.

There are at most two colors of each positive integer weight k+1. For weighted cost at most L, its generating function is bounded by

    F(e^(−u))=∏_{j≥1}(1−e^(−uj))^(−2),
    log F(e^(−u))
      =2 Σ_{r≥1} 1/[r(e^(ur)−1)] ≤ 4/u.

The number of multisets of cost ≤L is at most exp(uL+4/u); choosing u=2/√L gives exp(4√L). Missing minus letters only reduce this bound.

The catalogue has O(log N) letters. Canonical depth-first enumeration sorted by increasing q_i visits only admissible multisets, with at most one failed product bound per node. The remaining arithmetic and witness assembly are polynomial in log N. Thus the exact entry compiler uses

    exp(O(√log N)) arithmetic operations,

and every constructed source has order O(log N). All intermediate counts have O(log N) bits; polynomial bit costs do not change the displayed asymptotic form.

Allowing all t≥5 and all 99 digits contributes only O(log X) times this multiset bound to the number of distinct outputs through X. Hence the whole grammar supplies only exp(O(√log X)) targets through X. It cannot supply a positive-density odd tail, all ordinary primes, or all 7p atoms. For the last two statements, this bound implies convergence of the reciprocal sum of its outputs, contrary to prime reciprocal divergence if either entire prime lane were present.

The advance is therefore a verified growing physical profile supply and a complete, fast entry test for it. The full 99-digit bank is not being advertised as an absolute universal decoder or as a repeatable state-closed induction.

## Dependencies and certificate scope

The only inherited mathematical ingredients are the exact one-apex formula in S322 `ONE_APEX_CAPACITY.md`, ordinal-slot multiplicativity, and the S325 `FULL_99_DIGIT.md` bank and literal `corridor_witnesses.json`. The source wrapping identity, closed profiles, coprimality classes, infinite alphabet, inverse and complexity bounds are proved above. The source certificate is not a random-search-completeness assertion. The replay uses explicit checks that remain active under Python optimization.
