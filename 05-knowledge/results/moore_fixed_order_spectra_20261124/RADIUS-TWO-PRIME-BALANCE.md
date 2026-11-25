# Radius-two balance, prime readouts, and what survives for Chan–Pak

Use the labelled arc-flip cube, with all \(D=\binom n2\) arcs available; do not collapse to unweighted metagraph neighbours after isomorphism or converse merging. For \(H=H(T)\), let \(S_r(T)\) be the sum of Hamilton counts at exactly \(r\) distinct arc flips from \(T\). A **balanced shell** has \(S_r=\binom Dr H(T)\). A **balanced two-edit square** uses two specified arcs and has zero mixed difference on its four corners. The latter is a stronger labelled statement about a chosen pair of edits, not implied by an average across the full shell. This addresses the [Chan–Pak problem](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf) without assuming its converse or an all-size inverse.

## 1. Two genuinely balanced balls, and one square above 21

Number vertices from zero; at pair \(i<j\), bit 1 means \(i\to j\), in lexicographic pair order. With \(m=2\) on four vertices and \(m=153661494\) on eight, **every** shell is exactly balanced, about **3** and **315** respectively. This is exact arithmetic at every radius, so it also holds modulo every modulus for the **same** tournament. It does not say every neighbouring height is equal, and it is not coverage of arbitrary centres.

Here is an exact sphere formula. Let \(W_b\) count ordered vertex words with exactly \(b\) backward steps in \(T\). For a given word, a flip set making it Hamiltonian must flip exactly the \(b\) wrong consecutive pairs, no good consecutive pair, and any choice among the other \(D-n+1\) pairs. Thus
\[
          \sum_r S_r z^r=(1+z)^{D-n+1}\sum_b W_b z^b.
\]
Word reversal gives \(W_b=W_{n-1-b}\). So the first \(\lfloor(n-1)/2\rfloor+1\) shells, including radius zero, recover **all** shell sums over the integers. Full balance is equivalent to \(W_b=H(T)\binom{n-1}{b}\). The companion confirms this vector independently by subset DP and enumeration of all vertex words for both fixtures. Full balance also implies \(H(T)=n!/2^{n-1}\) by summing the shells, so for \(n\ge2\) it is possible only at a power of two (Legendre's formula). We do not infer sufficiency at larger powers from two cases.

For a *physical* two-edit square above the exceptional range, take the seven-vertex mask 158750 and mark arcs \(\{3,5\},\{2,6\}\), in that order. The heights are **131, 133, 135, 137** at no edit, first, second and both edits respectively. There are two ways to the fourth corner, but counting the endpoints once gives a perfectly additive 2 × 2 square. This is not a general construction of a square around each target.

Radius one alone cannot certify a marked square. At five-vertex mask 10, the common edit \(\{2,3\}\) gives 11 from 9. Both alternatives \(\{1,4\}\) and \(\{1,3\}\) give 13 from 9. With the common edit already made, however, they land at 15 and 11 respectively. So both marked radius-one readings are **(9,11,13)**, while their mixed differences are **0 and −4**. This is equality of the three actual integers, not just at one prime.

One can also separate the Moore analogy cleanly. In the labelled cube the two distinct orders of any two edits always give the same board. A Moore radius-two neighbourhood has unique physical routes, so it cannot be copied literally; the object to balance here is an endpoint, a marked square, or a shell sum, with explicit multiplicities.

## 2. Why radius two is enough data at the three smallest odd primes

For a marked set of \(k\) arc-sign bits, its Walsh character is an eigenvector of the shell operator \(S_r\), with eigenvalue
\[
 K_r(k)=[z^r](1-z)^k(1+z)^{D-k}.
\]
Indeed a flip changes the sign on exactly its intersection with the marked set. This is the usual [Hamming-scheme Krawtchouk calculation](https://doi.org/10.1016/j.indag.2025.05.011). Use the repo's [even-path forest Walsh expansion](../../../01-canon/theorems/THM-5313-switching-derivatives-of-hamiltonian-paths-and-cycles.md). For clarity, expanding the path products of \((1\pm s_e)/2\) shows that for a linear forest with \(k\) edges and \(j\) nontrivial paths, the coefficient has magnitude
\[
                  2^{j-(n-1)}(n-k)!.
\]
If any path is odd in edge length then reversing it cancels the coefficient; a forest of even-length paths has the stated magnitude because each block can be read in two orientations with the same sign, and its \(n-k\) blocks may be ordered arbitrarily. Other edge sets do not occur.

Modulo an odd prime \(p\), only the **even Walsh levels \(k\) with \(1\le n-k<p\)** survive. Conversely every such level between zero and \(n-1\) has a nonzero coefficient, by choosing just one even path or the empty forest. There are at most \(q=(p-1)/2\) surviving levels, with distinct residues mod \(p\). For a particular \(T\), let \(X_k\) be the contribution of its \(k\)-th Walsh level evaluated at \(T\). The shell data are exactly
\[
             S_r(T)\equiv\sum_k K_r(k)X_k \pmod p.
\]
On these distinct \(k\)'s the first \(q'\) rows, \(r=0,\ldots,q'-1\) for the actual number \(q'\) of levels, form an invertible interpolation system: each \(K_r(k)\) is a polynomial in \(k\) of degree \(r\) with leading coefficient \((-2)^r/r!\), a unit at \(p\). Thus the required radius is at most \((p-3)/2\). When the constant level is among them its contribution \(n!/2^{n-1}\) is already known and can be subtracted.

In particular **\(H,S_1,S_2\) determine every shell sum modulo 3, 5 and 7** for any one tournament. This is a theorem about the path count function, not about general Boolean functions. Modulo 2 all shells are trivial from the oddness of \(H\). An analogous support cutoff modulo \(p^a\) uses \(v_p((n-k)!)<a\); for powers one cannot automatically reuse the field interpolation argument.

**Sharp failure above 7.** Two nine-vertex masks 31070013745 and 45517909057 have the same first three shell sums mod 11, **(5,6,8)**, but radius three is respectively **1 and 2 mod 11**. The exact rows through radius three are
\[
(1215,44094,779358,8933310),\quad
(1589,56458,975994,10935850).
\]
Thus no nonlinear rule based only on those first three residues works at 11 either.

**Prediction is not balance.** At \(n=7,\ p=3\), the sole surviving Walsh level is 6, so every tournament's first two shells agree with the flat prediction mod 3. But \(K_3(6)\equiv0\) and \(\binom{21}3\equiv1\pmod3\): for the transitive tournament, \(H=1\) and \(S_3=33336\equiv0\). A perfectly balanced radius two mod 3 can fail one radius later even though the third sum is computable from the first three.

## 3. Every prime has a natural hierarchy of radii, and addresses are not integers

For a general function on the cube modulo odd \(p\), a repeated-stroll test is weaker than an exact sphere test. Write \(\sigma_i\) for the bit flip and \(L=\sum_i\sigma_i\). Since the \(\sigma_i\)'s commute and are involutions, \(L^p\equiv L\pmod p\). For \(0\le r<p\), the shell operator \(S_r\) is a polynomial in \(L\) (recurrence \(L S_r=(r+1)S_{r+1}+(D-r+1)S_{r-1}\)). At \(r+1=p\) this recursion cannot divide by \(p\). Modulo \(p\), the actual shells at radii \(1,p,p^2,\ldots\le D\) together distinguish all Walsh levels. To see this, use \((1\pm z)^{p^j}\equiv1\pm z^{p^j}\). In \(K_{p^j}(k)\), the coefficient is \((D-k)_j-k_j\) plus a term determined by the lower base-\(p\) digits of \(k\) and \(D-k\). With lower digits known, \((D-k)_j\) is \(D_j-k_j\) minus the known borrow. The digit \(k_j\) is thus recovered with slope \(-2\), a unit. The companion tests this on a 36-contact cube for six primes. This explains why radius two works as a special count-specific checksum and why the large-prime or higher-precision problem needs more observables.

The [two literal +2 families covering all powers of 13](THIRTEEN-ODD-ALPHABET.md) show why individual prime coverage is not the missing exact inverse. There is a general obstruction to promoting any fixed finite bank of diagonal dyadic widths to **simultaneous** coverage. If each bank member's lower height is an integer polynomial in at most \(d\) quantities \(2^{m_i}\), then modulo \(2^N-1\) it supplies at most \(N^d\) residues, because each power depends only on \(m_i\bmod N\). A bank of \(J\) such families supplies at most \(JN^d\), strictly fewer than \(2^N-1\) for large \(N\). This holds even if every member has zero mixed difference or a perfectly balanced marked square. The modulus is odd, so the missed residue classes contain odd integers. It does not assert the missed residues are missed by all tournaments.

The constructive test to carry forward is consequently twofold: keep a **specified** pair of commuting arc moves with verified mixed term, and give a **joint** modular lift whose heights stay within an explicit interval shorter than the joint modulus. A congruence to the target in that interval then forces the integer, whereas independent choices of chart for each prime do not. Nothing here supplies that final lifting mechanism at every target.

Run the companion for all exact fixtures and the modular interpolation controls.
