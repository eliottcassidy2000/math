# Reading tournament counts with Volchkov's kernel and decoding a perfect-number suffix

**PROVED.** The exact kernel in Volchkov's RH criterion extracts the Mellin moments of a tournament's ordered-slot path polynomial from its counting zeta product. A single rationally shifted reading encodes every coefficient by the Chinese remainder theorem. On the arithmetic side, for odd exponents at least five the four decimal digit parities of \(E_p=2^{p-1}(2^p-1)\) have least exponent period **100**. The exact four digits have least period **500**. At exponent 61 they are 2176. The four digits alone recover the exponent modulo 500 in that example.

The inheritance is the [pointed set of colored cycle packings](../chanpak_f1_pointed_rails_20261118/README.md), the [slot-polynomial parity and return operator](../chanpak_arc_fiber_grammar_20261122/SQUARE-EARS.md), and the [quadratic norm at the perfect neighbors](../chanpak_fermat_solids_seams_20261122/README.md). The directed square supplies two completions with the same ordinary height. The delayed [32-vertex marked example](../chanpak_arc_fiber_grammar_20261122/HIDDEN-FLIP.md) keeps even its two-slot reading. These are concrete tests for any statistic claiming to see the whole slot profile.

## Where the integers in the integral come from

**CITED.** Volchkov's 1995 paper proves that the equality in the question is equivalent to RH. A second primary treatment recovers it by integrating the argument of zeta by parts after the contour pole approaches the simple zeta pole at 1. The Laurent expansion \(\zeta'/\zeta(1+u)=-1/u+\gamma+O(u)\) supplies Euler's constant. In their formula with lower line \(b\), the finite correction accompanying \(\gamma\) is \(-3/(2(1-b))\). Putting \(b=1/2\) gives \(-3\). [Volchkov, original criterion](https://umj.imath.kiev.ua/index.php/umj/article/view/5435) · [Sekatskii, Beltraminelli and Merlini, Theorem 5a and Remark 2](https://umj.imath.kiev.ua/index.php/umj/article/download/2568/1895/6147)

Set \(a=1/2\) and \(g_a(t)=a^4t/(a^2+t^2)^2\). Then
\[
 g'_a(t)=\frac{a^4(a^2-3t^2)}{(a^2+t^2)^3}
         =\frac{1-12t^2}{(1+4t^2)^3}.
\]
The zero endpoint and the limit at infinity make the integral of this derivative zero. The repeated power 2, the resulting cubic denominator, \(4=a^{-2}\), and \(12=3a^{-2}\) follow from differentiation at the symmetry line. The scale producing 32 is calculated below. This explains the appearances of 0, 1, 2, 3, 4, 12 and 32 within the formula. The four-element square in the tournament grammar and the 32-vertex collision above are separate combinatorial objects.

These numerals already have precise, different roles in the tournament archive. At 0 the colored-packing evaluation keeps the empty packing. A transitive tournament has 1 Hamilton path, and an edited arc here has gap 2. The 3-vertex triangle receiver is excluded modulo 3 in the decimal note, while the surviving 4-vertex receiver at maximal order forces five terminal ones, meaning its address is 31 modulo 32. Earlier triangle-return work gives an odd-order signed-gap law modulo 12. Independently, the square-sum graph on the integers 1 through 32 is the first one in that family to have a Hamilton cycle and has a unique unoriented one. The integral calculation needs none of these numerical parallels. The map to the finite zeta product below explains exactly which observable carries between the two settings. [Pointed basepoint](../chanpak_f1_pointed_rails_20261118/README.md) · [Receiver order table](../chanpak_e61_fibonacci_target_20261122/SELECTOR.md) · [Square-return residues](../chanpak_arc_fiber_grammar_20261122/SQUARE-EARS.md) · [THM-4851, square-sum prehistory](../../../01-canon/theorems/THM-4851-the-prehistory-of-the-square-sum-numbers-unique-cycles-at-32-and-33-terminal-spectrum-with-the-hole-at-39.md)

**Source status, 2026-11-22.** Alpöge's November announcement reports a proof and describes chaining theta lifts to retain multiplicity. Our finite calculation asks which multiplicities can be recovered after collapsing many objects to a scalar observation. All the identities below follow directly from finite cycle packings, rational algebra and integration of a one-pole kernel. [Alpöge's announcement, pp. 1 and 3–4](https://alpo.ge/remark.pdf)

## The same kernel reads a counting polynomial over the field with one element

Let \(T\) have \(n\ge1\) vertices. Its ordered \(q\)-slot path-cover count, empty slots allowed, is a polynomial \(P_T(q)=\sum_{m=1}^n a_m q^m\) with nonnegative integer coefficients. Each packing of vertex-disjoint directed odd cycles contributes two choices per cycle. Each cycle and each uncovered vertex independently chooses one slot. Its exponent \(m\) is the number of these choices of slot, so \(m\equiv n\pmod2\), \(a_n=1\) and \(P_T(1)=H(T)\). This gives an actual disjoint union of affine spaces, one for each two-colored packing, with its coordinates indexed by the cycles and uncovered vertices. Counting its points over any finite field gives the same polynomial. It extends the earlier pointed-set construction by retaining independent slot coordinates. [Repository packing model](../chanpak_f1_pointed_rails_20261118/README.md) · [Grinberg and Stanley, cycle expansion](https://arxiv.org/abs/2307.05569)

The polynomial-counting prescription for an absolute zeta function gives \(\prod_m (s-m)^{-a_m}\). We use its reflected factors
\[
 Z_T(s)=\prod_m(s+m)^{-a_m},\qquad s>0 .
\]
The model above is a union of standard affine spaces. Connes and Consani ask about an absolute curve whose zeta function would be the completed Riemann zeta function. Here the counting object and its finite product are explicit. [Connes and Consani, §1.2 and §2](https://arxiv.org/abs/0903.2024)

For every real shift \(b\ge0\), the following three quantities agree:
\[
 \mathcal V_T(b)=\frac{32}{\pi}\int_0^\infty g'_{1/2}(t)
 \int_{1/2}^\infty \log\frac{|Z_T(\sigma+b+it)|}{Z_T(\sigma+b)}\,d\sigma\,dt
 =\sum_m\frac{a_m}{m+b+1}
 =\int_0^1 x^b P_T(x)\,dx .
\]
Indeed, for \(a,u>0\), integration by parts shows
\[
 \int_0^\infty g_a'(t)\log\frac{|u+it|}{u}\,dt
 =-a^4\int_0^\infty\frac{t^2\,dt}{(t^2+a^2)^2(t^2+u^2)}
 =-\frac{\pi a^3}{4(a+u)^2}.
\]
The last equality is the derivative in \(a\) of the elementary integral \(\int_0^\infty t^2/((t^2+a^2)(t^2+u^2))\,dt=\pi/(2(a+u))\). Set \(a=1/2\) and \(u=\sigma+b+m\). Integrating in \(\sigma\) gives \(-\pi/[32(m+b+1)]\). The negative exponent in the zeta product gives the displayed positive sign. The subtraction at \(t=0\) makes Fubini legitimate. For each mode the inner integral of the positive logarithm is at most \(\int_0^\infty\log(|u+it|/u)du=\pi t/2\), integrable against \(|g'_{1/2}|\). Thus each normalization, including 32, is fixed before choosing a tournament.

**Lossless single reading.** Put \(L=n!\) and \(d_m=1+(m+1)L\). For distinct exponents a common divisor of \(d_m,d_k\) divides \((m-k)L\), is coprime to \(L\), and hence divides \(m-k\). Since \(|m-k|\mid L\), that divisor is one. With \(D=\prod_m d_m\) over the possible exponents, the single integer
\[
 W=D\,\mathcal V_T(1/L)/L=\sum_m a_m D/d_m
\]
gives \(a_m\equiv W(D/d_m)^{-1}\pmod {d_m}\). The bound \(a_m\le H(T)\le n!\) and \(d_m\ge1+2L\) recover the actual coefficients. For two marked completions on the same vertices their single signed reading recovers every coefficient difference by centered residues, including its sum, the exact arc gap. Exact rational arithmetic is essential. The reading is an invariant of the polynomial, and two nonisomorphic hosts can share that polynomial.

**Small shifts instead.** At the same order put \(r=\lfloor(n-1)/2\rfloor\). In a difference the top monomial cancels, leaving \(r\) possible degrees in a set \(S\). With \(D(b)=\prod_{m\in S}(b+m+1)\), the numerator \(D(b)\Delta\mathcal V(b)\) has degree at most \(r-1\) with leading coefficient \(\Delta H\). Knowing the readings at \(b=0,\ldots,r-1\) determines it by interpolation and hence the entire slot difference by residues. In particular \(\Delta H=\Delta_b^{r-1}[D(b)\Delta\mathcal V(b)]_{b=0}/(r-1)!\) for \(r\ge1\), where \(\Delta_b\) is forward difference. At order 123 there are exactly 61 possible lower degrees. This is a reconstruction from supplied observations. It gives no tournament before those observations can be realized. In the six-vertex neutral pair of the square note the slot difference in the 376-minus-352 direction is \(2q^2(q^2-1)\), whose unshifted moment is \(-4/15\), even though its Hamilton-count difference is zero. [Square-return example](../chanpak_arc_fiber_grammar_20261122/SQUARE-EARS.md)

There is also a wholly combinatorial version in which the final ordinary Hamilton count of one tournament encodes the polynomial. Its width constraint directly excludes large transitive returns at exponent 61. The same analysis gives the smallest marked arc whose count stays equal through one return and differs after two. Its source has exactly ten vertices. [One physical count and the order-ten minimum](PHYSICAL-READING.md)

Finite-field evaluations give an arithmetic test of the same polynomial. One tournament realizes any prescribed joint profile at the seven odd primes through 19. A fixed marked six-vertex seed recovers the desired endpoint residues at the first six primes and exposes an obstruction at the seventh. [Joint realizability and the physical seed](FINITE-FIELDS.md)

The decimal proof and its arithmetic restriction on the order-123 receivers are in [DECIMAL.md](DECIMAL.md). Run `python3 -B check.py` for the finite residue tables, the three receiver congruences and an independent reconstruction from one rational reading on all tournament types through order five.
