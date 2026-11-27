# The 23-square exception: 364,723 tournament classes, only two with best route balance

Continue the [odd-modulus tournament classification](../collatz_modular_multiplication_tournament_classes_20261127/README.md). Set \(g=9/4\), and let \(h=\operatorname{ord}_p(g)\) be odd at the prime \(p\ne2,3\). A connection set consists of whole \(g\)-orbits, selecting exactly one of each pair of opposite orbits. Crucially, the first failure of the multiplicative order to grow changes both the graph-isomorphism problem and the best possible routing.

## 1. Full lifting gives a sequence of prime-order factors

Write \(a=v_p(9^h-4^h)\). By the binomial order-lifting law \(\operatorname{ord}_{p^b}(g)=h p^{\max(0,b-a)}\). Thus for a difference of valuation \(j\) mod \(p^e\), the allowed sign depends on its normalized unit mod \(p^{\min(a,e-j)}\): the corresponding subgroup is the entire inverse image of the subgroup at that smaller precision. In particular when \(a=1\), the sign depends only on the first nonzero base-\(p\) digit of the difference. All resulting tournaments are exactly the iterated lexicographic products of \(e\) independently selected admissible prime-order factors, in order from least to most significant digit.

There are \(c_p^e\) isomorphism classes, if \(c_p\) is the number at \(p\). For completeness, the maximal proper strong modules of a vertex-transitive tournament on a prime number of vertices have equal sizes (they form an invariant partition); hence they are singletons. The [modular decomposition theorem for tournaments](https://combinatorialpress.com/article/ars/Volume%20122/volume-122-paper-27.pdf) then gives indecomposability, since the alternative transitive prime-size quotient cannot be vertex-transitive. In a product of such factors an arbitrary module meeting more than one top fiber projects to a module in the outer factor: an omitted fiber must see all projected fibers uniformly. If the projection is the whole outer factor, then the original module contains each fiber entirely. Otherwise a vertex missing from a partly occupied fiber would see both directions from included vertices of other fibers, by regularity of the outer factor. Thus top fibers are canonical maximal proper modules and induction reads the factor sequence.

For regular outer adjacency \(A\) and inner \(B\), with inner order \(s\), write \(R(A)=A+A^2\). The product adjacency is \(A\otimes J_s+I\otimes B\), so
\[
R(A[B])=sR(A)\otimes J_s+I\otimes R(B).
\]
This uses \(BJ_s=J_s B=(s-1)J_s/2\). Each off-diagonal route count for a regular prime-size factor is in \([1,(p-1)/2]\). The upper bound follows by counting its out-neighbors, and if a destination beats the source and has no incoming two-step route from it, that destination would beat all the source's out-neighbors and the source, contradicting regularity. Consequently the route ranges at successive factors are disjoint and recover their individual histograms.

In particular there are \(2^e\) classes at \(19^e\); the row histogram identifies the entire digit sequence of Paley and non-Paley factors. At 25 all sixteen fixed-label orientations give the single class \(C_5[C_5]\). Its ordered Hamiltonian path count is exactly **2,398,512,090,439,546,875**. To replay, collapse a path to its word of consecutive fiber runs. For each outer vertex with \(j\) runs, the inner contribution is the number of ordered covers of \(C_5\) by \(j\) nonempty paths, namely \((15,90,210,240,120)\). A DP on the five run counts and last outer vertex sums the products of these weights. This avoids a 25-vertex subset DP.

## 2. Exact classification at a stalled square

When \(a\ge2\), put \(n=(p-1)/h\), an even integer. The unit layer at \(p^2\) is encoded by a sign word of length \(np\) with opposite signs half a turn apart, and the nonzero multiples of \(p\) by such a word of length \(n\). Common multiplication shifts both words by the same amount. However by [Klin, Liskovets and Pöschel, Theorem 4.15](https://www.mat.univie.ac.at/~slc/wpapers/s36klp.pdf), precisely when the unit layer is invariant under multiplication by \(1+p\), the two multipliers may be chosen independently. This invariant locus is the same as pulling back the entire unit layer from mod \(p\). In that locus the two sign words have effective length \(n\). Outside it the common multiplier is required. This is a graph isomorphism classification, not a count up to arbitrary choice of coordinates only.

For even \(N\) define \(b_N(j)=0\) if \(\gcd(N,j)\mid N/2\), and \(b_N(j)=2^{\gcd(N,j)/2}\) otherwise. This counts antipodal words fixed by a shift, since antipodes within one shift-cycle are impossible and otherwise cycles pair freely. The exact class count at the stalled square is therefore
\[
\frac{1}{np}\sum_{j=0}^{np-1} b_{np}(j)b_n(j)
-\frac{1}{n}\sum_{j=0}^{n-1}b_n(j)^2
+\left(\frac{1}{n}\sum_{j=0}^{n-1} b_n(j)\right)^2.
\]
The terms count common shifts everywhere, subtract common shifts on the invariant locus, and put back independent shifts there.

At \(p=23\), \(h=11,a=2,n=2\). The first stalled admissible square is therefore 529, with **\(2^{24}\) labeled orientations yielding exactly \(364723\) isomorphism classes**. One simple count normalizes the inner layer to Paley squares. There are \(2^{23}\) unit words; the residual shift group has order 23. The two invariant words both give the single lexicographic product. All other words have free orbits of size 23, giving \(1+(2^{23}-2)/23=364723\). Lower admissible primes through 19 have full lift.

## 3. Exact optimal balance at 529

The Paley product has one-row route histogram \(\{6:22,138:506\}\), with value to multiplicity. Reversing just the \(g\)-orbit of 1 (and its negative) gives \(\{28:22,134:44,136:154,137:44,138:242,140:22\}\). That already improves the sum of squared deviations from the forced mean 132.5, from 367356 to 251680.

The exact optimum over the entire Collatz-compatible circulant family is much better: **5588** per row, attained by exactly **two isomorphism classes**. The next distinct energy is 5896. One minimizer has route values ranging from 125 to 139. The scope is fixed \(g\)-compatible circulants at 529, not every 529-vertex tournament.

Here is the exhaustive computation in short. Fix the inner sign by a multiplier, leaving 23 independent unit signs \(s_i\). With signed skew matrix \(S=A-A^\top\), one has \(R(A)=(529J-I+S^2)/4\). Thus 16 times the row energy is \(\sum_{t\ne0}((S^2)_{0t}-1)^2\). For one representative \(t\) per antipodal orbit pair, \((S^2)_{0t}\) is a quadratic sign polynomial given by cyclic convolution of the first row with itself. The 24 pairs have size 22 each. Square the polynomials in the algebra \(s_i^2=1\), add with that multiplicity, and use the Walsh–Hadamard transform to evaluate the resulting exact degree-four polynomial at every sign assignment. The integer bounds are safe in binary64: the sum of absolute coefficients before evaluation is at most \(528\cdot529^2\), since a quadratic convolution has at most 528 summands and the offset contributes one more. The replay counts 46 normalized minimizers, in two orbits of size 23. It independently checks a minimizer against actual routes, and recovers both comparison energies. Thus the failed lift makes real extra balance possible, not merely more labels for the same tournament.

Perfect flatness is nevertheless impossible at a composite cyclic order. In a regular circulant \(AA^\top=kJ-R(A)\). Thus flat routing is equivalent to the connection set being a cyclic skew-Hadamard difference set, and [the cyclic classification quoted in Ó Catháin, Theorem 45](https://doras.dcu.ie/29113/1/NonaffineHadamard.pdf) allows only prime-order Paley examples.

Replay with two Node scripts: one verifies the orbit counts, histograms and exact 25-path count; the other is the Walsh enumeration of the 529 optimum.
