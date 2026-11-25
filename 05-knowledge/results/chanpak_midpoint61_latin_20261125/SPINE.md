# How much freedom must an exact Farey inverse have?

Let \(E_p=2^{p-1}(2^p-1)\), including composite exponents. On the lower side of a persistent two-flank construction write
\[
H(x,y)=k+qxy+tx+cy+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y)).
\]
Here \(f(0)=0,\ f(2a+1)=f(a),\ f(2a)=f(a)+a+[a\text{ odd}]\), and \(0\le f(a)\le a\). In the Farey charts below, reversing the same physical arc gives \(H+2\). Hence the exponent-61 question is the single integer equation \(H=E_{61}-1\). [Farey chart and marked edit](../chanpak_decimal_kernel_20261122/E15-FAREY.md)

## Independent repeat counts still give only finitely many exponents

**PROVED.** Assume \(q,h,t,c,k\) are positive integers, \(r,e\) are nonnegative integers and \(q\) has an odd prime divisor. Fix an integer \(d\) and two slender regular languages of binary addresses. Here slender means that the number of words of any given length has a fixed finite upper bound. Only finitely many exponents \(p\ge2\) admit words from these two languages with \(H=E_p+d\). A finite cutoff is computable from the chart and automata. This includes every fixed positive Farey chart built around the exponent-15 middle, since its two-slot coefficient contains the factor \(320\). The two word lengths and the exponent can be chosen independently. [Middle factorization and the bounded-run obstruction](../chanpak_decimal_kernel_20261122/RUN-COMPLEXITY.md)

For example, each flank may choose its own number of copies of a mixed block such as 01 or 01101. Such inputs have unboundedly many binary runs. The earlier mixed-block result covers an affine repetition schedule in a common parameter. The obstruction here also permits an arbitrary target-dependent scheduler, including a Turing machine choosing both repeat counts. The finite languages of candidate templates can be arbitrarily large. [Strict drift and the earlier schedule theorem](../chanpak_absolute_liouville_leech_20261122/UNEQUAL-FLANK-INVERSE.md)

**Proof.** A slender regular language is a finite union of single words and languages \(uv^*w\) with nonempty \(v\). To see this, trim a deterministic finite automaton to states on accepting paths. A strongly connected component with two internal outgoing edges from one state produces two distinct first-return walks with different first letters. Repeating each to a common length gives exponentially many distinguishable words of equal lengths by concatenation. Every cyclic component must therefore be a directed cycle. If an accepting route traverses two cyclic components, add laps to both while preserving their combined length in arbitrarily many ways. The exit position from the first component distinguishes the accepted words by determinism. Thus every accepting route crosses at most one cyclic component. There are finitely many routes with complete laps removed, giving the asserted description.

Fix two such templates. For a word with \(m\ge1\) copies of a block of length \(a\), set \(U=am\). Its numeric value is \(x=\alpha2^U+\beta\). The inherited append law gives
\[
 f(x)=\alpha_f2^U+\gamma m+\beta_f,
\]
with rational coefficients, \(\gamma\ge0\), and \(\gamma>0\) exactly for a block containing both digits. The parity stabilizes after the first copy. For instance, appending a block of value \(v\) and length \(a\) sends the numeric address \(z\) to \(2^a z+v\). The fixed suffix modifies the exponential coefficient and constant of the read. The cited strict drift proof supplies the sign of the linear coefficient. Treat zero copies, constant numeric values and singleton words as fixed-address cases.

For two growing templates write \(V=bn\) on the other side. Clearing denominators gives the necessary equation
\[
 A2^{U+V}+B2^U+C2^V+Z = D(2^{2p-1}-2^{p-1}),\qquad Z=Pm+Qn+G. \tag{1}
\]
Here \(D,A>0\), \(A/D=q\alpha\alpha'\), \(P,Q\ge0\), and all coefficients are integers. The offset \(d\) is absorbed in \(G\). Since \(H\ge(1+x)(1+y)\) and each growing template satisfies \(1+x\gg 2^U\), we have \(U+V\le2p+O(1)\). In particular \(|Z|=O(p+1)\).

First exclude a solution with \(2^{\min(U,V)}=O(p+1)\) along any unbounded sequence of exponents. Suppose \(2^U\) is the smaller power. Collect the terms in \(2^V\) on the left of (1). Their coefficient \(A2^U+C\) has magnitude \(O(p+1)\), while the right side after moving \(B2^U+Z\) has magnitude comparable to \(2^{2p}\). Therefore \(V\ge2p-O(\log p)\). Reducing (1) modulo \(2^{p+v_2(D)}\) now forces
\[
 B2^U+Z\equiv-D2^{p-1}\pmod {2^{p+v_2(D)}}.
\]
Every representative of that residue has magnitude at least \(2^{p-1+v_2(D)}\), whereas the left side is \(O(p+1)\). Exchanging the two words proves the other case. If a template is fixed, (1) has just one exponential and a linear remainder. Its exponent likewise equals \(2p+O(1)\) at any sufficiently large solution, and the same reduction excludes it. Two fixed words give a constant height.

If \(Z\ne0\), equation (1) forces \(2^{\min(U,V,p-1)}\mid Z\), giving the excluded small power for all sufficiently large \(p\). Consequently \(Z=0\). If either block is mixed, its strictly positive drift and the nonnegativity of both drifts bound the corresponding repeat count by the fixed constant \(|G|\). This too gives an excluded small power.

Both blocks must therefore be pure. For these blocks the positive leading coefficients \(\alpha,\alpha'\) are integers. Indeed, if the prefix has numeric value \(j\) and the suffix has length \(s\), the leading coefficient is \(2^s j\) for a repeated zero and \(2^s(j+1)\) for a repeated one. On any remaining unbounded solution sequence both \(U,V\) tend to infinity and
\[
 2^{2p-1-U-V}
 =\frac{A+B2^{-V}+C2^{-U}}{D(1-2^{-p})}\longrightarrow A/D.
\]
Integer powers of two are discrete near this positive limit. Hence \(A/D=q\alpha\alpha'\) must be a power of two, contrary to its odd prime divisor.

For effectiveness, the affine word formulas are rational and explicit. A lower bound for \((1+x)/2^U\) is the minimum of \(\alpha\) and its value at the first copy. This gives a computable linear majorant for \(|Z|\). In the pure case choose a constant threshold on both \(U,V\) from the positive distance between \(A/D\) and the two nearest powers of two. Below that threshold use the small-power argument. For nonzero \(Z\) use its linear majorant, and for mixed zero-drift-sum solutions use the constant bound on the number of repeats. All remaining comparisons are between an exponential and a known linear polynomial, or its square after clearing the small coefficient. Increasing \(p\) until domination and monotonicity hold gives a cutoff per template pair. Taking their maximum completes the proof.

**Two-scale form.** The same proof gives a criterion for other prescribed odd counts. Fix nonzero integers \(a,b\), with \(a>0\), and a constant offset \(d\). Suppose the target is
\[
                    a2^i+b2^j+d,
\]
where \(i>j\) and \(\min(j,i-j)/\log(i+2)\longrightarrow\infty\). If an odd prime divides \(q\) but does not divide \(a\), a pair of slender regular languages can attain only finitely many such targets along this sequence. In (1) the drift is now \(O(i+1)\). A small power on one side forces the larger exponent to reach \(i-O(\log i)>j+v_2(Db)\). The remaining small terms must have absolute value at least \(2^{j+v_2(Db)}\), which exceeds their linear bound. For two pure blocks the limit forces \(q\alpha\alpha'/a\) to be a power of two, excluded by the stated odd prime. Thus the obstruction comes from two separated active binary scales. Neither primality of the exponent nor perfectness of the center is needed.

**Automaton consequence.** If a fixed pair of regular address languages supplies exact crossings for infinitely many exponents, at least one of their trimmed deterministic automata must have either two cyclic components on one accepting route or two internal letter choices at one state of a cyclic component. In the first case the number of accepted words is at least linear in the word length on an arithmetic progression of lengths. In the second case it is at least exponential on such a progression. The equal-length lap constructions in the proof give these counts. This is a necessary branching cost for a regular inverse. It does not assert that each accepted word participates in a crossing.

## Five altered low bits on the smaller flank at exponent 61

**FINITE-EXACT.** Take either of the two known order-118 spine-four address pairs below, and either of the canonical spine-four or spine-eight charts. In every nonnegative integer solution to \(H=E_{61}-1\), the **numerically smaller** physical address must differ in at least **six** of its low 40 bits from the corresponding side of **each** reference pair. Every higher bit and the entire opposite address are free. If the addresses tie, both tests apply. The reference pairs are literal near witnesses only on spine four. [Marked reference arcs](../chanpak_odd_square_spines_20261122/NEAR-WITNESSES.md) · [Earlier two-sided 80-bit radius-two result](../chanpak_odd_square_spines_20261122/RATIONAL-CARRIES.md)

| Reference | Left address | Right address | Reference midpoint minus \(E_{61}\) on spine four |
|---|---:|---:|---:|
| 0 | 735435556549430 | 1575929833923500 | −844710 |
| 1 | 4814468855604222 | 240731609083994 | +184814588 |

Put \(N=E_{61}-1\). For a candidate smaller left address \(x\), positivity bounds \(x\le\lfloor\sqrt{(N-k)/q}\rfloor\). Define \(R=N-k-tx-r[x\text{ odd}]-hf(x)\) and \(A=qx+c\). Then every opposite address \(y\) satisfies
\[
  x\le y\le\lfloor R/A\rfloor,\qquad (A+h)y+e\ge R,\qquad
  Ay+e[y\text{ odd}]\equiv R\pmod h.
\]
For each parity solve the displayed linear congruence, take every resulting candidate in the interval, and evaluate \(Ay+hf(y)+e[y\text{ odd}]\) exactly. To treat a smaller right address, swap \((t,r)\) with \((c,e)\). Each run visits all \(\sum_{j=0}^5\binom{40}{j}=760{,}099\) suffixes and all eligible high words. The concatenation formula for \(f\) lets the high word advance cheaply. A floating estimate in one quotient branch is corrected in both directions using exact integer multiplication. On the tested balls the smallest possible \(x\) exceeds 32 million, and all products fit unsigned 128-bit arithmetic. [Reversed suffix identity](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md)

| Spine | Reference | Smaller side | Eligible address rows | Final opposite words tested | Hits |
|---:|---:|---|---:|---:|---:|
| 4 | 0 | left | 744145395 | 1147271 | 0 |
| 4 | 0 | right | 744211473 | 51888629 | 0 |
| 4 | 1 | left | 744145397 | 1703937 | 0 |
| 4 | 1 | right | 744137720 | 1224700 | 0 |
| 8 | 0 | left | 93500724 | 300 | 0 |
| 8 | 0 | right | 93574539 | 133775 | 0 |
| 8 | 1 | left | 93515635 | 443 | 0 |
| 8 | 1 | right | 93500724 | 286 | 0 |

The eight runs give 3,350,731,607 address rows counting multiplicity, and 56,099,341 final digital reads. The result excludes a repair even if all forty low bits of the opposite flank may change. Unrestricted spine-four and spine-eight solutions, and an exponent-61 physical crossing in other families, remain open. The all-exponent language theorem gives a separate constraint on a rule that would select such solutions across exponents.

From this directory, reproduce with:

~~~sh
c++ -std=c++17 -O3 spine_repair.cpp -o /tmp/chanpak61-spine-repair
for spine in 4 8; do
  for ref in 0 1; do
    for side in 0 1; do
      /tmp/chanpak61-spine-repair "$spine" "$ref" "$side" 5
    done
  done
done
/tmp/chanpak61-spine-repair 4 0 0 2 planted
/tmp/chanpak61-spine-repair 4 1 1 2 planted
~~~

The last two controls replace the target with the actual lower height of each reference arc and recover its two addresses using the same search.
