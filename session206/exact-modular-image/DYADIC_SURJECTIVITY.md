# Congruence completeness of O/Z tail images

Status: unconditional all-order theorem ruling out a specific obstruction strategy, not a covering theorem.

Let `s=(h,R,e,d)` be any physical unit-start seed covered by the session205 arithmetic formula, with h and d odd and R>d>0. For every k≥1 and every odd residue a modulo 2^k, its O/Z scalar image contains arbitrarily large numbers congruent to a modulo 2^k.

Thus no fixed dyadic residue obstruction can disprove eventual coverage by a proposed finite bank. Indeed this obstruction is unavailable even for one seed. The final theorem strengthens this to every modulus for the root seed. These statements do not imply any interval coverage or positive image density.

## Proof

For a positive integer t, put

    m_t = 3 sum_(j=0)^(t−1) 2^[k+j(k+2)].

Its binary expansion consists of t blocks `11`, separated by k zero bits, with at least k low zero bits. Therefore `m_t≡0 (mod 2^k)` and `m_t` is even. The pair formula for Q gives `Q(m_t)≡2t (mod 2^k)`: each internal adjacent pair contributes 2, while a pair in different blocks has distance at least k+1 and therefore contributes a multiple of 2^k. No additional adjacency term crosses blocks.

The exact tail formula consequently gives

    H_s(m_t)≡h−2dt (mod 2^k).

Since d is odd, t running through any complete residue system modulo `2^(k−1)` makes the right side run through every odd residue modulo 2^k. Choose a positive t realizing a. For every L≥0, Q(2^L m_t)=Q(m_t) and 2^L m_t remains divisible by 2^k. Hence

    H_s(2^L m_t)≡a (mod 2^k).

These scalar values tend to infinity, since R>0 and

    H_s(2^L m_t)=h+R 2^L m_t−dQ(m_t).

All are actual O/Z descendants by the arithmetic reduction. This completes the proof.

## Scope

The theorem concerns a fixed modulus at a time. It does not rule out an obstruction using moduli that grow with the target, or congruence constraints combined with the finite inverse bound. The strengthening below handles arbitrary moduli under a coprimality hypothesis.

## Stronger theorem: all fixed moduli for a seed coprime to their odd part

Let `M=2^k u`, where `k≥1` and u is odd. If `gcd(d,u)=1`, the same seed has arbitrarily large outputs in **every odd residue modulo M**. In particular, the root seed (`d=1`) has this property for every even modulus M, and meets every residue modulo every odd modulus. Therefore no obstruction based solely on omission of a fixed congruence class can exclude the root image, or a bank containing it.

Choose `b≥max(k,2)` such that `2^b≡1 (mod u)`; when u=1 any such b works. Put `x=2^b−1` and, for each positive t, define

    m_t=x sum_(j=0)^(t−1) 2^[k+j(b+k)].

These are t blocks of b consecutive ones, separated by k zeros. We have `m_t≡0 (mod M)`. A block's internal Q contribution is

    c=Q(2^b−1)=2^b−2.

This identity follows immediately by induction from the recurrence for Q. For two blocks whose starting positions differ by D≥b+k, their complete cross contribution is

    x^2 2^(D−b).

There is no cross-block adjacency. This cross contribution is divisible by u (because u divides x) and by `2^k` (because D−b≥k). Thus

    Q(m_t)≡tc (mod M),
    H_s(m_t)≡h−dtc (mod M).

As `b≥2`, c has exactly one factor of 2; as `2^b≡1 (mod u)`, we have `c≡−1 (mod u)`. Hence `gcd(dc,M)=2`. The residues `h−dtc` therefore run through all odd residues as t runs through a complete residue system modulo M/2. Multiplying any selected m_t by arbitrarily large powers of 2 preserves its divisibility by M and its Q value, and makes H_s tend to infinity. This proves the strengthened theorem.

The hypothesis on d matters: for example the physically reachable seed `(h,R,e,d)=(5,6,2,3)` satisfies `H(m)≡2−2(m mod 2) (mod 3)`, so its image omits the residue 1 modulo 3. The root has no such issue.

This does not exclude congruence arguments whose modulus grows with the target and which also exploit the finite inverse bound.

## Exact modular image for every seed

There is a complete version without the coprimality hypothesis. Fix any even modulus M and write

    g=gcd(R,d,M/2).

The image of H modulo M is **exactly** the residues satisfying either

    N≡h (mod 2g), or N≡h+R−e (mod 2g).

Every permitted residue has arbitrarily large actual O/Z outputs. The two classes may coincide. Since `R−e` is even, both are odd classes. Odd moduli are handled by projecting this statement for 2M onto M.

### Necessity

The fact that Q is even gives `dQ(m)≡0 (mod 2g)`. For even m, `Rm≡0 (mod 2g)`; for odd m, `Rm≡R (mod 2g)`. These prove the two necessary classes.

### A block-pumping lemma

For any fixed nonnegative integer a and even M, the image contains arbitrarily large outputs in each residue of the form

    H(a)+2d z (mod M), z in Z.

Use `M=2^k u`, b, x, and c as in the previous construction. Choose A at least k larger than the binary length of a, and set

    m_(a,t)=a+x sum_(j=0)^(t−1) 2^[A+j(b+k)], t≥1.

No high block overlaps a. The block-to-block cross contributions vanish modulo M as before. A cross term from a bit of a at position i to a bit of a high block beginning at B is part of the sum

    x sum_(i: bit_i(a)=1) 2^(B−i−1).

This is divisible by u because u divides x, and by `2^k` because `B−i−1≥k`. No cross adjacency term occurs. Hence

    m_(a,t)≡a (mod M),
    Q(m_(a,t))≡Q(a)+tc (mod M),
    H(m_(a,t))≡H(a)−dtc (mod M).

Since `gcd(c,M)=2`, varying positive t gives exactly the coset generated by 2d. Taking t in any fixed residue class modulo M preserves the output residue and makes m_(a,t) tend to infinity; the inverse lower bound makes H tend to infinity. This proves the lemma.

### Sufficiency

Put `q=gcd(2d,M)`. By the lemma, it suffices to identify the H(a) residues modulo q. Since q divides 2d and Q(a) is even, for a=2j and a=2j+1 respectively these residues are

    h+2Rj, and h+R−e+2Rj (mod q).

As j ranges over nonnegative integers, each progression fills its entire coset modulo

    gcd(2R,q)=gcd(2R,2d,M)=2g.

Thus the two necessary classes are also sufficient, with arbitrarily large outputs in each. This completes the exact characterization.

This result concerns modular projections, not equality to individual integers. A set can have full modular projections and still omit infinitely many integers. In particular, congruence completeness of the root does not settle the finite-bank question.
