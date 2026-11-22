# First appearance is usually far from unique

**FINITE-EXACT** census through order nine; **PROVED** asymptotic consequence of the inherited birth concentration. Here \(H(T)\) counts ordinary directed Hamiltonian paths, and birth multiplicity is the number of **tournament isomorphism classes at the least order**, not labelled tournaments, Hamiltonian frames, or signed tilings.

## 1. The exact first-birth census

| Order \(n\) | New values | First-birth classes | New values with one class | Self-converse first-birth classes | Greatest multiplicity, and where |
|---:|---:|---:|---:|---:|---|
| 1 | 1 | 1 | 1 | 1 | 1 at 1 |
| 2 | 0 | 0 | 0 | 0 | none |
| 3 | 1 | 1 | 1 | 1 | 1 at 3 |
| 4 | 1 | 1 | 1 | 1 | 1 at 5 |
| 5 | 4 | 6 | 2 | 6 | 2 at 9 and 15 |
| 6 | 12 | 33 | 1 | 9 | 5 at 33 and 37 |
| 7 | 58 | 328 | 11 | 72 | 15 at 91 |
| 8 | 243 | 5,457 | 2 | 145 | 49 at 345 |
| 9 | 1,200 | 163,676 | 25 | 2,408 | 257 at 1,467 and 1,503 |

Both maximal order-nine fibers consist of one self-converse class **and 128 distinct converse pairs**, and are entirely strong. All other first births through order nine of order at least three have strong witnesses only as well. The smaller maxima at orders 6, 7, and 8 likewise have exactly one self-converse class each.

The very first multiplicity is not merely a tournament and its converse. At the order-five first birth of \(H=9\), two outgoing-mask lists are
\[
(16,1,3,7,14),\qquad (4,1,18,7,11).
\]
Both are individually self-converse under the same label reversal \(i\mapsto4-i\), and both have sorted score vector \((1,1,2,3,3)\). They are nevertheless nonisomorphic: their multiset of deletion counts \(H(T-v)\) is respectively \((1,1,5,5,5)\) and \((1,3,3,5,5)\). Thus even quotienting by converse **and** knowing the scores does not make this first witness unique. At the other first multiple birth, \(H=15\), the two classes already split by their scores: \((1,2,2,2,3)\) versus \((2,2,2,2,2)\).

The inherited [order-nine birth atlas](../chanpak_birth_atlas_metagraph_20261117/metagraph/PROOF.md) proves the class generation by new-vertex induction and reports new-value and class totals, with one explicit representative per attained value. Here we emit every class in the same complete induction, then tally all 1,520 individual first fibers and their converse structure. The new-vertex induction takes each complete preceding list, tests every adjacency pattern of an appended vertex, and canonicalizes by outdegree followed by all degree-preserving relabelings. The retained orders are thereby exact. Subset endpoint DP computes \(H\), and the additional anti-involution test backtracks over degree-complementary matchings and verifies every arc. An anti-involution exists for every self-converse tournament, as proved below, so the test is complete.

## 2. Structural causes and asymptotic multiplicity

There are two exact sources of multiplicity even before searching a strongly connected fiber. Converse preserves \(H\) because reversing the Hamiltonian words is a bijection. An ordered SCC decomposition multiplies the component counts, and its ordered class word determines its isomorphism class. Thus even a fixed component multiset can contribute several birth witnesses when that multiset achieves the minimum total size. The [inherited factor classification](../chanpak_birth_atlas_metagraph_20261117/factor/PROOF.md) gives the precise all-order recurrence by choosing the first strong SCC, plus the palindrome formula for the self-converse classes. Counting arbitrary binary splits overcounts the same ordered component word. For squarefree \(H>1\), a self-converse first witness must be strong, because an off-center nontrivial SCC in a converse palindrome would force a square divisor. None of this replaces strong-fiber enumeration: all first-birth multiplicities in the table above order two already occur **within** strong classes.

More broadly, let \(A_n\) denote all classes of order \(n\), and let \(m(h)\) be the number of classes at the least order of count \(h\). For all sufficiently large \(n\), **all but an \(\exp(-\Omega(\sqrt n))\) fraction of \(A_n\)** belong to entirely strong first-birth fibers with
\[
m(h)\ \ge\ e^{-\sqrt n}\frac{2^{\binom n2}}{(n!)^2}
      \ =\ 2^{\binom n2-O(n\log n)} .
\]
Moreover, in each fiber selected by this assertion, the proportion of self-converse classes is at most \(2^{-n^2/4+O(n\log n)}\). So almost every tournament class at large order is not just a first witness: its count first appears in a huge family of witnesses, overwhelmingly made of nonfixed converse pairs.

**Proof of the order and strength premise.** The [inherited concentration result](../chanpak_birth_atlas_metagraph_20261117/families/PROOF.md#what-the-typical-count-reads-the-outdegree-score-imbalance) puts \(H\) in \((\mu_n/2,3\mu_n/2)\), where \(\mu_n=n!/2^{n-1}\), outside a fraction \(\exp(-\Omega(\sqrt n))\) of classes. It uses the even Fourier layers and hypercontractivity, then transfers through Burnside. That argument uses [Divoux–Wang's maximum asymptotic](https://arxiv.org/html/2610.23755v1), which implies \(M(n-1)=O(\mu_n/n)\) for the maximum \(M(n-1)\) at the preceding order. Thus the typical counts exceed \(M(n-1)\); by ordinal padding all earlier orders have no larger maximum either. The same inherited note §1 excludes any decomposable realization above \(M(n-1)\) for large \(n\), by the factorial bound on a nontrivial product of two smaller SCC parts. Hence the whole order-\(n\) fiber, not just the sampled class, is first and strong in this range.

**Proof of the multiplicity and converse bounds.** There are at least \(2^{\binom n2}/n!\) tournament classes at order \(n\), and at most \(n!\) possible values of \(H\). The number of classes in fibers smaller than the displayed threshold is at most that threshold times \(n!\). Dividing by the class lower bound gives the exceptional fraction \(e^{-\sqrt n}\).

Every automorphism group has odd order: an involutory nontrivial automorphism would interchange endpoints of an arc and reverse it, so Cauchy's theorem precludes even group order. If \(\alpha\) is any anti-isomorphism from the tournament to itself, \(\alpha^2\) has odd order \(r\), so \(\alpha^r\) is an anti-involution. At most one vertex can be fixed, since an arc joining two fixed vertices would have to reverse itself. All such involutions on \(n\) vertices are conjugate to the same prescribed matching, with one unmatched vertex when \(n\) is odd. Relabel any self-converse class to use that prescribed witness. On unordered arcs there are exactly \(\lfloor n/2\rfloor\) fixed pairs, the matched pairs themselves, and all other orbits have size two. The orientation of a single member of an orbit can be chosen freely and forces the other. Hence there are at most
\[
2^{(\binom n2+\lfloor n/2\rfloor)/2}=2^{\lfloor n^2/4\rfloor}
\]
self-converse classes of order \(n\) in total. Divide by the displayed lower bound on the size of any selected birth fiber. The proportion claim follows.

This asymptotic is for a uniformly chosen **class**, not a uniformly chosen attainable count. It does not prove that almost all attainable values, much less almost all odd integers, have huge birth multiplicity. It also gives no all-order assertion of strongness for exceptional, low-count first births.

## 3. Tournaments in the same birth fiber are numerous, but no fiber dominates

In fact the probability of choosing any one count from uniform order-\(n\) isomorphism classes is \(O(n^{-2/5})\), uniformly in the count. In particular, there are \(\Omega(n^{2/5})\) **distinct** counts per sufficiently large order satisfying the huge, entirely strong birth-fiber conclusion in §2.

Here is a quantitative reason. For a uniform labelled tournament let \(X_e\) be the independent Rademacher signs on unordered edges, choose each edge column \(B_e\) of \(B\) to be the corresponding signed incidence vector of \(K_n\), and let \(s=BX\). The [inherited score expansion and Fourier bounds](../chanpak_birth_atlas_metagraph_20261117/families/PROOF.md#what-the-typical-count-reads-the-outdegree-score-imbalance) give
\[
Y_n=\sqrt{n/2}\bigl(H/\mu_n-1\bigr)=Q(X)+R,\qquad
Q(X)=-\gamma\bigl(\|BX\|^2-n(n-1)\bigr),\qquad
\gamma=\frac{\sqrt{n/2}}{n(n-1)}.
\]
The \(L^4\) norm of \(R\) is \(O(n^{-1/2})\): the inherited homogeneous layers \(f_k\) of degree \(2k\), \(k\ge2\), have \(L^2\) norm at most \((2/\sqrt n)^k/\sqrt{k!}\) and \(L^4\) norm at most \(3^k\) times that by hypercontractivity. The triangle inequality with the \(\sqrt{n/2}\) normalization gives the claim.

For replacement use the **multilinear** version \(Q(x)=-2\gamma\sum_{e<f}(B_e\cdot B_f)x_ex_f\), equal to the expression above on signs. Replace one Rademacher input at a time by an independent standard real Gaussian. With all other inputs fixed, the polynomial is \(U+x_e V\), where \(V\) is a sum of at most \(2(n-2)\) independent centered variables with coefficients of magnitude \(2\gamma=O(n^{-3/2})\): only edges meeting \(e\) interact. A fourth-moment expansion, using fourth moments at most 3, gives \(\mathbb E|V|^4=O(n^{-4})\). Taylor expand a smooth test function through order 3 around \(U\). The terms through degree 3 agree after averaging over the replaced variable, since both inputs have mean zero, variance one, and zero third moment. Each replacement therefore costs \(O(n^{-4}\|f^{(4)}\|_\infty)\); all \(\binom n2\) cost \(O(n^{-2}\|f^{(4)}\|_\infty)\).

For the Gaussian input \(g\), \(BB^\top=nI-J\), so \(\|Bg\|^2\) has distribution \(n\chi^2_{n-1}\). The Gaussian score expression \(-\gamma(\|Bg\|^2-n(n-1))\) differs from the **multilinear** \(Q(g)\) by the diagonal random term \(-2\gamma\sum_e(g_e^2-1)\). Its \(L^4\) norm is \(O(\gamma n)=O(n^{-1/2})\) by independence and a fourth-moment expansion. The Gaussian score expression itself has Kolmogorov distance \(O(n^{-1/2})\) to the standard normal by Berry–Esseen on the \(n-1\) independent centered squares; its variance is \(n/(n-1)\), an asymptotically negligible rescaling.

Sandwich each threshold indicator between smooth steps of transition width \(\epsilon=n^{-2/5}\), with fourth derivative \(O(\epsilon^{-4})\). Replacement costs \(O(n^{-2}\epsilon^{-4})\). Dropping either the Rademacher high-layer residual or the Gaussian diagonal costs \(O(n^{-2}\epsilon^{-4})\) by the fourth-moment Markov inequality, plus widening the threshold interval by \(\epsilon\). The continuous normal law supplies a bound \(O(\epsilon)\) on all such windows. Together these bounds give
\[
\sup_x\left|\Pr[Y_n\le x]-\Phi(x)\right|=O(n^{-2/5}).
\]
Continuity of \(\Phi\) bounds the mass of every atom by twice this distance. The same atom bound holds for uniform isomorphism classes: an invariant event's labelled proportion and class proportion differ by at most the exponentially small nonidentity Burnside contribution used in the inherited concentration argument. Finally the huge, entirely strong first-birth fibers from §2 carry \(1-o(1)\) of the class probability, and each single fiber has probability \(O(n^{-2/5})\). Thus at least \(\Omega(n^{2/5})\) distinct values carry that conclusion. This is an anti-concentration estimate on the class distribution, not a uniformity theorem for the attained integers.

## Replay

Compile the class emitter and feed its output to the converse analyzer:
    g++ -O3 -std=c++17 birth_classes.cpp -o /tmp/birth-classes
    /tmp/birth-classes 9 > /tmp/birth-classes.tsv
    python3 check_birth_multiplicities.py /tmp/birth-classes.tsv --fibers /tmp/birth-fibers.tsv
