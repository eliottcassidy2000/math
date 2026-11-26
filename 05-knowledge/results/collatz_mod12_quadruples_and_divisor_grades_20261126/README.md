# Squares, two odd steps, and four-factor tournament grades

This builds on the [three Collatz rows](../arithmetic_braids_20260917_collatz.md), the [four-tile divisor fold](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md), the [balanced three-component selector](../chanpak_ternary_buckets_20261122/DIVISOR-TYPES.md), and the [class series through 73](../chanpak_fibers_to_73_and_collatz_triple_codes_20261126/README.md).

## 1. The literal \(9/4\) in Collatz and the six four-step words

On odd \(x\), reduction modulo 24 gives
\[
x^2\equiv x^4\equiv\begin{cases}1&3\nmid x,\\9&3\mid x,\end{cases}
\qquad x^3\equiv x.
\]
Thus at this resolution a square is a projection (onto two idempotents), a cube retains the residue, and the fourth power repeats the square. Doubling is a different operation: modulo 12 it reads only \(x\bmod6\). In general the additive differences of squares recover the product via \((x+y)^2-(x-y)^2=4xy\). This distinction between forgetting under powers and recovering a product by additive differences is important.

With one halving per step, \(T(x)=x/2\) for even \(x\) and \((3x+1)/2\) for odd \(x\). Two consecutive odd steps are allowed exactly on \(x\equiv3\bmod4\), and there
\[
T^2(x)+1=\frac94(x+1).
\]
Repeated all-odd blocks just square this dilation; there is no arithmetic-complexity implication.

A **four**-step word of weight **two** instead acts as \((9x+b)/16\), the square of the slope \(3/4\). The additive carry depends on the order of the same two odd steps. Its exact mod-12 readout is:

| word (time order) | \(b\) | legal input mod 16 | output mod 12 |
|---|---:|---:|---|
| 0011 | 20 | 12 | 2, 5, 8, 11 |
| 0101 | 14 | 2 | 2, 5, 8, 11 |
| 0110 | 10 | 6 | 1, 4, 7, 10 |
| 1001 | 11 | 13 | 2, 5, 8, 11 |
| 1010 | 7 | 1 | 1, 4, 7, 10 |
| 1100 | 5 | 3 | 2, 5, 8, 11 |

The carry recurrence when appending bit \(d\) at position \(j\) is \(b\mapsto3^d b+d2^j\). Along the legal input ray \(r+16k\) the output changes by \(9k\), so refining the input from mod 16 to mod 64 resolves all four output residues. Its residue mod 3 depends on the time of the last odd step, not the number of odd steps alone.

Closing a four-step word gives \(7x=b\), with \(7=4^2-3^2\). Exactly the alternating words close at integer 1 and 2, the two phases of the ordinary halved cycle. The other four are rotations of a single extension to rationals with odd denominator, \(5/7\to11/7\to20/7\to10/7\to5/7\), using numerator parity. These are not integer Collatz cycles. The denominator 7 has a literal origin here; it is not a reason that 7 or 21 is absent in the tournament spectrum.

## 2. Four-factor divisor geometry and the exact largest component grade

In the old \(F=U+S\) convention, \(F\) counts all proper nontrivial divisors, \(U\) its prime elements and \(S\) its squarefree elements. The balanced shapes \(p,p^3,p^2qr\) have token counts 1, 3, 4. If all involved primes are Hamiltonian heights, the largest numbers of nonsingleton strong components are exactly 1, 3, 4. Write \(G_d(t)\) for the polynomial counting strong types of height \(d\) by vertex number. In the three cases the singleton-free top-grade series are \(G_p,\ G_p^3,\ 12G_p^2G_qG_r\). The coefficient 12, rather than 24, is the permutation count after folding the two identical prime positions; the same involution produces the \((3,2,2)\) divisor box from the four-bit cube. There is also an additive reading: the fourth mixed difference of \(\binom{x}{4}\) in directions \((p,p,q,r)\) is \(p^2qr\), with lower moments cancelling. The two equal directions fold the signed subset terms, preserving their multiplicities. This does not equate a Pollock representation or a tiling with a prime factorization.

Here is a more general effect of the forbidden count **7**. Let \(h=7^e m\), \(7\nmid m\); assume every prime dividing \(m\) is itself realized. Write \(M=\Omega(m)\), counting prime tokens with multiplicity. For \(e\ge2\) the largest possible number of nonsingleton strong components is
\[
g=M+\lfloor e/2\rfloor.
\]
Moreover the maximal words have exactly the following block multisets, alongside the unused non-7 primes as individual blocks:

* \(e=2k\): \(k\) blocks of height 49, and no other blocks carrying 7.
* \(e=2k+1\), \(k\ge1\): either \(k-1\) blocks of height 49 and one of 343, or \(k\) blocks of height 49 and one of \(7p\) for a prime token \(p\mid m\) with \(G_{7p}\ne0\). In the second case remove one \(p\) from the unused list. The option \(p=3\) is excluded because \(G_{21}=0\).

Proof: if \(r\) blocks carry 7, at least \(\max(0,2r-e)\) must contain just one 7. Each such block must consume at least one non-7 token, so the grade is at most \(M+r-\max(0,2r-e)\). Equality at its maximum forces the stated token assignments, including no merging among the remaining non-7 primes. The known strong heights 49 and 343 establish attainability. More generally if \(e=1\), the largest grade is \(M+1-\min\{\Omega(d):d\mid m,\ G_{7d}\ne0\}\), with an empty minimum meaning no witness. Indeed a single block must absorb all the 7 tokens and the remaining primes can be used individually. This extends the earlier three-component criterion without assuming Chan–Pak.

For any listed multiset with counts \(r_d\) of each height \(d\), its maximal-grade numerator is \(g!\prod G_d^{r_d}/\prod r_d!\); divide by \((1-t)^{g+1}\) to restore arbitrary singleton padding. For instance at \(441=3^2 7^2\) the maximal grade is three, with the exact grade-three class series \((21t^{13}+3t^{14})/(1-t)^4\), using \(G_{49}=7t^7+t^8\) and \(G_3=t^3\). At \(1715=5\cdot7^3\) there are two mechanisms at maximum grade two, namely \(\{5,343\}\) and \(\{35,49\}\), since 35 is attainable. Compositeness and a repeated prime both matter, but a forbidden prime changes which grade can read them.

## 3. The genuine matrix connection: triangles as a tensor and a fourth-degree correction

The [November 24 matrix multiplication preprint](https://mathcommunity.github.io/oai-math-revised/papers/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026/paper.pdf) reports the arithmetic exponent \(9/4\). Its triangle tensor \(\sum_{i,j,k}e_{ij}\otimes e_{jk}\otimes e_{ki}\) is literally the weighted triangle readout \(\operatorname{tr}(XYZ)\) of three cyclic parts; taking four tensor copies groups four coordinate positions into side \(n^4\). At exponent \(9/4\) that side has arithmetic complexity scale \(n^{9+o(1)}\), not a finite identity with exactly \(n^9\) operations. This is an exact tensor connection to triangle readouts and separate from the numerical equality of the Collatz slope.

For a tournament adjacency matrix \(A\), a directed arc \(i\to j\) belongs to \(e_{ij}=(A^2)_{ji}\) cyclic triples; \(t_i=\sum_j A_{ij}(A^2)_{ji}\) counts cyclic triples through \(i\). If there are \(c_3\) cyclic triples, the number of disjoint triangle pairs is
\[
d_{33}=\binom{c_3}{2}-\sum_i\binom{t_i}{2}+\sum_{i\to j}\binom{e_{ij}}{2}.
\]
A pair meeting in two vertices gets subtracted twice by the vertex sum and added back once by the arc sum. A pair meeting in one vertex is subtracted just once. This is a direct square of local triangle incidences with overlap removed, not an independence model. The 5-cycle count is \(c_5=\operatorname{tr}(A^5)/5\), since any nonsimple closed directed 5-walk would have to split into two shorter directed cycles, each of length at least three. Thus a constant number of matrix products separates the two contributions \(2c_5+4d_{33}\) to the quadratic coefficient of the normalized path-cover polynomial. This directly computes what was hidden in the 17, 27 and 33 trades, without enumerating Hamiltonian paths. It does not accelerate the unrestricted path count by the same exponent.

## 4. A sharper mod-12 readout: keep the longest odd-cycle bit

Let \(P_T(z)\) count path covers with \(z\) labelled colors, one path per color including empty colors. Write \(s=1\) for odd \(n\ge3\), and \(s=2\) for even \(n\ge4\). For any \(b\ge1\), the [Rédei–Berge cycle-pack formula](https://arxiv.org/abs/2307.05569) gives
\[
\begin{split}
P_T(2^b)&\equiv 2^{bs+1}\,c_{n-s+1}\pmod {2^{bs+2}},\\
P_T(2^b)&\equiv (-1)^{bn}H(T)\pmod3.
\end{split}
\]
In each term the exponent of \(z\) is \(n-\sum_C(|C|-1)\); it has the parity of \(n\) and is at least 1 (odd) or 2 (even). Two or more cycles already contribute a coefficient divisible by 4. A single longest odd cycle attains the stated exponent and coefficient 2, and shorter single cycles cannot contribute at that bit. The empty term is divisible by \(2^{bs+2}\) in the indicated ranges. Mod 3, the exponents all have the same parity, yielding the second line.

At \(b=1\) this upgrades reduction mod 12 to **mod 24 at odd order** and **mod 48 at even order**, using one additional cycle-parity bit. For the six-vertex witnesses of height 33, chord masks 12 and 22 both give \(P_T(2)=360\), with 9 and 7 five-cycles. Mask 25 has \(P_T(2)=336\) with 6 five-cycles. All three agree mod 12, but the 24-unit difference is detected mod 48. This is a concrete restriction on which equal-count repairs can change the next readout bit.

The replay checks the parity-word table, 1,396 independent matrix, cycle and subset-path controls, and 48 arithmetic top-grade models with available primes 3 and 5. Run `python3 replay.py` in this directory.
