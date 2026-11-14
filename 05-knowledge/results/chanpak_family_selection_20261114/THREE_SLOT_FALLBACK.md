# A third-slot fallback for exact singleton-join misses

The cyclic join of nonempty blocks A,B,C has a literal run decomposition.
Every maximal consecutive block of same-color vertices in a Hamiltonian path
is a nonempty path, and the run colors advance cyclically through A,B,C.
For each block X let g_j(X) be the number of ordered covers of its vertex
set by j nonempty paths; let F_j(X) allow empty paths in its j slots.
For a fixed run-color word with visit vector (x,y,z), its count is exactly
g_x(A)g_y(B)g_z(C), with no correction terms. This is a specialization of
the visit-kernel formula THM-5133; the rank-three kernel below is also
recorded in the independent audit of THM-5135.

## Replace the singleton with a two-vertex terminal reader

Write `C3[A,B,C]` for the cyclic join and `TT_n` for an n-vertex
transitive tournament. Since `TT_2` has g=(1,2,0,...), its kernel
`C3[A,B,TT_2]` in the nonempty-cover coordinates is

    [5 3 0]
    [3 7 2]
    [0 2 2].

Indeed three starting colors give the (k,k,k) patterns, one starting
color each gives the three (k+1,k,k) and three (k+1,k+1,k) permutations,
and the third coordinate admits k≤2 only. Taking B=TT_2 reads the first
block by

    H(C3[A,TT_2,TT_2]) = 11 g_1(A)+17 g_2(A)+4 g_3(A)
                      = 4 F_3(A)+5 F_2(A)-11 H(A).

The conversion uses F_2=2g_1+g_2 and F_3=3g_1+3g_2+g_3. With A=TT_m
the output is 4·3^m+5·2^m−11, and m=2,3,4 give 45,137,393 exactly.
They are strong by their cyclic block quotient.

## 45 and 137 admit no singleton cyclic join at any order

Let h=H(A), f=F_2(A). The singleton join of A and B has
`N=(f_A−h_A)(f_B−h_B)+2h_A h_B`. We use f−h≥2^n−1 from
the nonempty subset choices, and f≥(n+1)h by cutting every full
Hamiltonian path in any of n+1 places. At n≥2, f is divisible by four:
the odd-cycle expansion at two colors contributes
2^(n−|V(pack)|+2|pack|) per packing, with exponent at least two.
In addition f≡(−1)^n h mod 3, because the two-color weighting
relative to the one-color packing weight is 2^(n−|V(pack)|+|pack|),
and the additional exponent differs from n by an even integer.
The same specialization gives the general congruence
`F_q(T)≡ q^n H(T) mod (q^2−1)` for integers q≥2.

If both inputs are nontransitive, each has h≥3 and f−h≥9
(the triangle has f−h=9 and order≥4 implies f−h≥15).
Their count is at least 99, and if either has order≥4 it is at
least 153. The pair of two directed triangles gives exactly 99.
If both inputs are transitive, N−2 is a product of two Mersenne
numbers. For 45 this requires factoring 43 with a Mersenne factor
43, impossible. For 137 the target product is 135, whose only
Mersenne divisors are 1,3,15, and none has a Mersenne complement.

It remains to take the second input transitive of order k and the
first nontransitive, in which case N=2h+(2^k−1)(f−h).

- For 45, k≥3 already gives N≥69. For k=2 the budget forces the
  first block to have order≤3, giving 33. For k=1 the output is
  f+h. A nontransitive block of order≥5 has f≥48: in the positively
  weighted odd-cycle expansion its empty packing contributes at
  least 32 and one triangle contributes at least 16. Lower-order
  nontransitive choices give f+h in {15,27,37}.
- For 137, k≥4 already gives N≥141. For k=3 the budget gives
  f−h≤18, hence order≤4, and only the triangle fits (N=69).
  For k=2, 3f−h=137 forces h≡7 mod12. At n≥5 the cut bound
  implies h≤137/(3(n+1)−1)≤8, so only the forbidden count 7
  could fit, and the orders≤4 cannot meet the congruence.
  For k=1, f+h=137 forces n even by the mod-3 identity and n≤7
  since f≥2^n. Order≤4 is too small, so n=6. Write t_3,t_5,t_33
  for the numbers of triangle, 5-cycle, and disjoint-pair-of-triangles
  packings. Then
  `f+h=65+34t_3+10t_5+20t_33` cannot be 137: the excess 72
  forces `17t_3+5t_5+10t_33=36`, hence both t_3≡3 mod5 and t_3≤2.

For the low-order inputs above, the nontransitive three-vertex
block has (h,f)=(3,12); at order four the nontransitive
signatures are (3,24),(5,32), by the single-triangle contribution
and the fact there are one or two directed triangles.

## Why carrying only H and F₂ cannot simulate this reader

Two literal eight-vertex sources, in outgoing mask notation, are

    (44,49,50,38,41,192,31,95)
    (92,201,42,112,6,19,52,125).

They have identical (H,F_2)=(345,5640), but F_3=44505 and
45225, respectively. Thus every singleton cyclic-join reader
gives the same ordinary output to these inputs. The two-vertex
fallback reads 202425 and 205305, respectively. Adding s ordinal
singletons within each source preserves h, multiplies F_2 by 2^s
and F_3 by 3^s. The exact output difference becomes 2880·3^s
even though the singleton readers still cannot distinguish them.

Run `python replay_three_slot.py` to reproduce the five source
profiles and all full Hamiltonian counts by literal subset DP.
The third-slot reader is an exact fallback for specific targets,
not a proof of coverage for every admissible odd target.
