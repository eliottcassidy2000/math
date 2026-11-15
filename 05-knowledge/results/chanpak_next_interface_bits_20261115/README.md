# More arithmetic bits of the unmarked tournament interface

The [integer interface theorem](../chanpak_symmetry_interfaces_20261115/OBSERVABLES.md) writes the complete state for arbitrary unmarked substitution as d_r=(r−1)!c_r, where c_r counts ordered nonempty r-path covers. Here we give exact states through modulus 128 and a much stronger restriction when two modules of the same order agree in their first slot readings.

**Main new restriction.** Suppose U,V have the same order and P_U(j)=P_V(j) for j=1,…,s, where P counts ordered slots including empties. The difference of their Hamiltonian counts in **any** common unmarked context is divisible by
```
2 s! (2s+1)!                  at odd module order,
2 (s+1)! (2s+1)!              at even module order.
```
At s=2 the bounds are 480 and 1,440 respectively. A literal strong order-eight pair and a strong context attain the even-order 1,440 difference exactly. Thus equality of H and the second slot does not mean equality of the full interface, but the first possible difference is much larger than two.

## 1. Integer factorization behind the stronger modulus

The [odd-cycle packing formula](https://arxiv.org/abs/2307.05569) for a nonempty n-vertex tournament is
```
P_U(k) = sum_C 2^j k^(n−2z),
    j=number of nontrivial odd cycles in the packing C,
    z=sum_{gamma∈C} (|gamma|−1)/2.
```
The empty packing gives the monic term k^n. All other coefficients are even. The polynomial has exponents of the same parity as n and no constant term.

Let γ=1 for odd n and 2 for even n. For two modules of the same order the monic terms cancel. Their difference equals 2k^γ times an integer polynomial in k². Agreement at 1,…,s forces distinct integer roots 1²,…,s² in this latter polynomial. Division by the corresponding monic factors stays integral. Therefore
```
P_U(k)−P_V(k) = 2 k^γ product_{j=1}^s (k²−j²) R(k²),    R in Z[x].
```
The factor k product(k²−j²) is a product of 2s+1 consecutive integers, so every integer evaluation is divisible by 2(2s+1)!. The differences in c_r are integer finite differences of these evaluations, and they vanish for r≤s. Multiplying by (r−1)! for r>s yields the odd-order bound, then the context theorem applies.

In even order there is an extra factor k. At r=s+1, the finite difference for c_r reduces to the single term k=s+1, giving the extra factor s+1. At r≥s+2 it is supplied by (r−1)!. This proves the even-order bound. If s≥floor((n−1)/2), the factor degree and cancellation of the leading monomial force equality of the *entire* polynomial at that fixed order.

**Sharp even-order example for s=2.** Take strong modules with outgoing rows
```
U = (120,85,89,114,192,86,128,47),
V = (158,92,64,116,196,23,161,46).
```
They have H=345 and P(2)=5,640 both; P_U(3)=44,505, P_V(3)=45,225. Their polynomial difference is exactly 2k²(k²−1)(k²−4), from the factorization and the evaluation at 3. The [level-three probe](../chanpak_symmetry_interfaces_20261115/OBSERVABLES.md) has coefficient 2 at c_3 and no higher reads. Hence it gives strong outputs whose counts differ by 1,440. This is also a pair with identical mod-32 contextual states but different mod-64 states.

## 2. Exact contextual states at successive powers of two

Write H=P(1), W=P(2), ε=n mod 2 and C_l for directed l-cycles counted up to rotation. Up to mod 16 only the first three coordinates survive:
```
d_1=H,    d_2=W−2H,    d_3=2[P(3)−3W+3H].
```
For mod 32 and 64, d_4=6[P(4)−4P(3)+6W−4H] is also potentially present. All higher d_r are zero at these moduli because r!(r−1)! divides d_r. The following states are **necessary and sufficient**, not only upper bounds on needed information:

| Modulus | Module orders | Exact state |
|---:|---|---|
| 8 | n≥2 | H mod 8, ε, χ, where W≡4χ mod 8 |
| 16 | n≥3 | H and W mod 16, ε |
| 32 | n≥3 | H and W mod 32, n mod 4 |
| 64 | n≥4 | H and W mod 64, n mod 4, floor(n/4)+J mod 2 |
| 128 | n≥5 | H and W mod 128, n mod 4, the two-bit drift below, and (only at odd order) floor(C_n/4) mod 2 |

Here χ is C_n mod 2 at odd order, zero at even order ≥4, and one at order 2. J is the parity of the number of directed cycles of lengths 3 mod 4 (all such lengths, not just spanning). The two exceptional sizes at the bottom may simply be treated by their literal c vectors; the small-order exceptions for the last two rows concern the reduction to these convenient invariants, not failure of the d interface.

**Modulo 8 and 16.** All exponents of P have parity ε, so P(3)≡3^ε H mod 8. Since W is divisible by 4 for n≥2, d_3≡4ε mod 8 and still determines ε mod 16. H and W are recovered from d_1,d_2. The packing expansion gives W≡4C_n for odd n≥5, W≡8C_(n−1) for even n≥4, and W≡8+4C_3 at n=3, all modulo 16: a single cycle has weight 2^(n−l+2), and two packed cycles have weight at least 16.

**Modulo 32.** For a packing of defect z, 3^(n−2z)≡3^n(1−8z) mod 16. Terms with nonzero z have even coefficient, and the empty packing has z=0; consequently P(3)≡3^n H mod 16. In d_3 this distinguishes all four n classes modulo 4 since 3 has period four mod 16 and H is odd. To calculate d_4 mod 32 only c_4 mod 16 is needed. For n≥3, P(4) is 8C_n at odd n and zero at even n mod 16. In the odd case this contribution cancels 6W mod 16, and −4P(3)−4H is zero. In the even case 6W vanishes and −4P(3)−4H is 8. So d_4≡16(1−ε) mod 32 adds no extra bit.

**Modulo 64.** Now P(3) is needed mod 32 and 9^(−z)≡1−8z mod 32. The empty term still has no correction, terms with two selected cycles are multiples of 32 after correction, and a single selected cycle contributes 16 exactly if its length is 3 mod 4. Thus
```
P(3) ≡ 3^n (H+16J)  (mod 32).
```
Since 3^4≡1+16 mod 32, this reads precisely n mod 4 and the mixed bit floor(n/4)+J mod 2. For d_4 mod 64 use c_4 mod 32. At odd n≥5, P(4)≡8C_n mod 32 and W≡4C_n mod 16 recovers the needed C_n mod 4. At even n≥4, P(4)≡0 mod 32. The remaining terms of c_4 are determined by H,W and the indicated P(3) state. The strong order-eight example in §1 certifies the mixed bit is not redundant.

**Modulo 128.** Put Z = Σ_C 2^j z, with packing size j and defect z as in §1. It is even, and Z/2 modulo 4 reads only single cycles and disjoint pairs of odd cycles. Since 9^(−z)≡1−8z mod 64,
```
P(3) ≡ 3^n (H−8Z)
     ≡ 3^(n mod 4) H [1+16(floor(n/4)−H^(-1) Z/2)]  (mod 64).
```
The expression in parentheses modulo 4 is the two-bit drift, with the inverse of odd H taken modulo 4. This is exactly the extra information in d_3 once H,W are known. Now d_4 reads P(4) mod 64. At even n≥6, P(4)≡32C_(n−1), already determined by W. At odd n≥5, P(4)≡8C_n and W determines C_n only modulo 4, leaving exactly its bit floor(C_n/4) mod 2 as a readout. Finally d_5=24c_5 only needs c_5 mod 16. All P(j) for j≤4 are already known to this precision and P(5)≡5^n H mod 16 by the same even-coefficient argument as for P(3). Thus d_5 adds no state coordinate. All r≥6 vanish modulo 128 by the factorial divisibility. The table calls these states necessary and sufficient for contextual equivalence; it does not claim every formal combination is realized.

More generally, if a is odd, write a²=1+8u. For modulus 2^m the defect expansion is
```
P_U(a) ≡ a^n sum_{r=0}^{floor((m−1)/3)} (8u)^r
                  sum_C 2^j binom(−z,r)  (mod 2^m).
```
This is the binomial expansion of (1+8u)^(−z). For even a the term instead vanishes once j+v_2(a)(n−2z)≥m. Thus odd evaluations read a bounded set of defect moments, while even evaluations read packings with few cycles and small uncovered remainder. This supplies a finite sufficient description at higher powers, without asserting its coordinates are all independent.

## 3. Simultaneous modular swaps: the first distinguishing bit

For a quotient Q and modules U_i of orders at least two, put h_i=H(U_i). Let A_i count quotient schedules with i repeated twice and others once, B_i those with i three times and others once, and A_ij those with i,j each twice. Then
```
H(Q[U_i]) ≡ (H(Q)−2 Σ_i A_i) Π_i h_i
          + 4 Σ_i [A_i χ_i+(B_i/2) ε_i] + 4 Σ_(i<j) A_ij   (mod 8).
```
The B_i are even by the excursion swap. Run counts at least four, or three labels repeated, vanish modulo 8 by c_r divisibility by r!. With two repeated labels and one visiting three times, swapping its two internal excursions again gives an even schedule coefficient. The excursions cannot coincide: each needs both an out-neighbor and an in-neighbor label, so coinciding would repeat at least *two* other labels, not just one. The two even cover factors therefore annihilate this term modulo 8. With both labels visiting twice, c_2≡2 mod 4 gives 4. The single-label terms use c_2≡4χ−2h mod 8 and c_3≡2ε mod 4; these give the displayed formula.

For a physical example take equal-H order-five modules with rows (24,17,3,6,12) and (16,9,3,5,14). Both have H=15 and their directed spanning cycle counts are respectively 2 and 3. Swap them in slots 0 and 1 of the strong quotient (6,12,8,1), keeping TT₂ in the other two slots. Here A_0=3 and A_1=2, so the output changes by 4 modulo 8. The exact strong order-14 counts are 2,270,565 and 2,286,945. This is a bit-level change, not an exact shift by four.

There is also a simultaneous **parity shadow modulo 3**. For any module, c_2≡0 if its order is odd and c_2≡−h if even; all c_r for r≥3 vanish. This is from P(−1)=(−1)^n H. Hence the count modulo 3 is the product of the ordinary h_i times the quotient count obtained by replacing every odd-order module by a singleton and every even-order module by TT₂. A swap that preserves the parity mask therefore also preserves H modulo 3. Together with mod-4 permutation invariance, any swap of two nonsingleton modules of the same parity in a fixed quotient preserves H modulo 12. In particular the physical swap above cannot be promoted to an exact shift by four in this way.

**Replay:** `python verify_nextbits.py` checks the two literal swap counts, strong roots and the sharp 1,440 probe difference, with independent subset path and slot counts. The all-order divisibility and state results follow from the proofs rather than these controls.
