# Session31: the mirror is a tree function and an extinction duality

**Status:** PROVED below, using the classical Gelfond–Schneider theorem where stated; FINITE-EXACT and numerical checks are explicitly separated. This is a self-contained synthesis of classical mechanisms, not a novelty claim or an open-problem solution.

## 0. Inheritance, scope, and the question actually answered

The supplied coefficient prefix is

\[
1,\quad\frac32,\quad\frac83,\quad\frac{125}{24},\quad\frac{54}{5},\quad\frac{16807}{720}.
\]

For an integer \(n\ge2\), let \(r_n>0\), \(r_n\ne n\), satisfy
\(r_n^{1/r_n}=n^{1/n}\). The exact structure is:

1. The prefix is the sequence of coefficients of \(z,z^2,\ldots,z^6\) in the classical tree quotient \(R(z)=-W_0(-z)/z\).
2. \(r_n=R(\log n/n)\) for every integer \(n\ge3\), but **not** for \(n=2\).
3. The rational mirrors are exactly \(r_2=4\), \(r_4=2\). All the other integer-indexed mirrors are transcendental.
4. For \(n\ge3\), \(r_n/n\) is exactly the extinction probability of a Poisson(\(\log n\)) Galton–Watson process. Conditioning on extinction changes its offspring mean to \(\log r_n\).
5. There is an honest tree/forest interpretation, but no established map transporting this structure to the tournament Hamiltonian-path gaps 7 and 21. A concrete forest-versus-path counterexample is given below.

A finite prefix alone never forces an infinite continuation. The coefficient identification becomes canonical because the mirror equation forces \(R=e^{zR}\), not because six rational numbers uniquely determine a function.

**Inheritance pass.** Closest proved mechanism: Lagrange inversion for \(T=ze^T\), and the classical extinction-conditioned branching law. Canonical hostile: \(n=2\), where the analytic branch returns the original number, not its mirror. Corrected near miss: the existing tournament 7/21 exclusion is not a new result (repository MISTAKE-647); matching arithmetic constants do not supply a target-preserving map (MISTAKE-658). Least-used sidecar here: the branch label and the extinction normalization.

**Concept board:** (i) two real inverse branches; (ii) labelled trees/rooted forests; (iii) Galton–Watson extinction and finite total progeny; (iv) rational versus algebraic-irrational powers; (v) tournament Hamiltonian paths. The first four are connected by explicit maps below. The fifth presently has a useful hostile, not a transported theorem.

## 1. Existence, uniqueness, and branch selection

Put \(f(x)=\log x/x\), \(x>0\). Then

\[
f'(x)=\frac{1-\log x}{x^2}.
\]

Consequently \(f\) increases from 0 to \(1/e\) on \([1,e]\), then decreases to 0 on \([e,\infty)\). For each \(x>1\), \(x\ne e\), there is exactly one distinct positive mirror. It lies on the opposite side of \(e\). There is no distinct positive finite mirror at \(x=1\) or \(x=e\); for \(0<x<1\), the negative level has only one positive preimage.

Write \(z=\log n/n\). The equation is \(\log r=zr\). Setting \(t=zr\) gives \(te^{-t}=z\), so its two positive solutions are

\[
r=-\frac{W_0(-z)}{z},\qquad r=-\frac{W_{-1}(-z)}{z}.
\]

For every integer \(n\ge3\), the input \(n>e\) occupies the \(W_{-1}\) branch, and the mirror occupies \(W_0\). For \(n=2\), these roles reverse: the principal branch returns 2, and the distinct mirror is 4. In particular \(z_2=z_4=\log 2/2\): the scalar \(z\) loses which member of the pair was the input.

The real mirror map, extended by \(r(e)=e\), is a decreasing involution on \((1,\infty)\). This is a branch exchange, not ordinary reflection across a straight line.

## 2. Complete arithmetic classification for integer inputs

**Theorem.** For integers \(n\ge2\), \(r_n\) is rational if and only if \(n\in\{2,4\}\). Otherwise \(r_n\) is transcendental.

**Proof of rational classification.** Write a rational mirror as \(r=a/b\) in lowest positive terms. The equality of logarithmic levels implies \(n^r=r^n\), hence

\[
n^a b^{nb}=a^{nb}.
\]

A prime dividing \(b\) divides the left side but not the right, a contradiction. Thus \(b=1\): a rational mirror of an integer must itself be an integer. The smaller member of a distinct positive integer pair lies in \((1,e)\), and so is 2. On the decreasing branch, 4 is the unique solution sharing the level of 2. This gives only \((2,4)\) and \((4,2)\).

**Proof of transcendence.** If a remaining \(r_n\) were algebraic, the preceding paragraph makes it algebraic irrational. Gelfond–Schneider makes \(n^{r_n}\) transcendental, because \(n\) is algebraic and is neither 0 nor 1. But \(n^{r_n}=r_n^n\) would be algebraic. Contradiction. The positive real power is one of the powers covered by the theorem. This argument proves transcendence of each individual mirror, not of infinite sums involving the mirrors.

Primary reference for the imported theorem: A. Gelfond, *Sur le septième problème de Hilbert* (1934), [original article record](https://www.mathnet.ru/rus/im4924). The primary record and indexed original text were checked; direct full-PDF retrieval timed out, so this is an invocation of the standard theorem, not a new audit of its proof.

## 3. Exact coefficients, convergence, and remainder bounds

Let \(T(z)\) be the analytic solution of \(T=ze^T\), \(T(0)=0\). Lagrange inversion gives

\[
[z^k]T(z)=\frac1k[u^{k-1}]e^{ku}=\frac{k^{k-1}}{k!}.
\]

Thus \(R(z)=T(z)/z=e^{T(z)}\), with the removable value \(R(0)=1\), has

\[
R(z)=\sum_{m\ge0}a_mz^m,
\qquad a_0=1,\quad
 a_m=\frac{(m+1)^{m-1}}{m!}\quad(m\ge1).
\]

Explicitly,

\[
R(z)=1+z+\frac32z^2+\frac83z^3+\frac{125}{24}z^4
 +\frac{54}{5}z^5+\frac{16807}{720}z^6
 +\frac{16384}{315}z^7+\cdots.
\]

This corrects two possible indexing errors: the displayed input starts at the coefficient of \(z\), not the constant coefficient; its numerator is \((m+1)^{m-1}\), not \(m^{m-1}\) with the same factorial.

An independent exact recurrence follows from \(R'=R^2/(1-zR)\):

\[
a_{m+1}=\frac{m+2}{2(m+1)}\sum_{j=0}^m a_ja_{m-j}.
\]

Stirling's formula gives

\[
a_m=\frac{e^{m+1}}{\sqrt{2\pi}\,m^{3/2}}
 \left(1-\frac{19}{12m}+O(m^{-2})\right).
\]

The radius of convergence is exactly \(1/e\), and the series converges absolutely on the entire closed circle \(|z|=1/e\). On the positive boundary, \(R(1/e)=e\), while the derivative diverges. Hence evaluating a convergent value at the branch point does not create two distinct roots there.

For effective evaluation use

\[
\frac{a_{m+1}}{a_m}=\left(1+\frac1{m+1}\right)^m<e.
\]

For \(0<z<1/e\), \(M\ge0\),

\[
0<R(z)-\sum_{m=0}^Ma_mz^m
 \le\frac{a_{M+1}z^{M+1}}{1-ez}.
\]

This is a proved real-number bound. The verification script evaluates it with high-precision decimals, not directed-rounding interval arithmetic. At the positive boundary the remainder instead has the slower asymptotic
\(e\sqrt{2/\pi}\,M^{-1/2}\).

The hardest integer input is \(n=3\), for which \(ez=e\log3/3\approx0.99544594\). Slow convergence here is explained by proximity to the branching singularity, not a failure of the series.

For fixed \(K\), as \(n\to\infty\),

\[
r_n=\sum_{m=0}^K a_m\left(\frac{\log n}{n}\right)^m
 +O\left(\left(\frac{\log n}{n}\right)^{K+1}\right).
\]

In particular \(r_n\downarrow1\) and
\(r_n=1+\log n/n+\tfrac32(\log n/n)^2+O((\log n/n)^3)\).

The Lambert-W/tree identities and their classical provenance are documented in Corless–Gonnet–Hare–Jeffrey–Knuth, *On the Lambert W Function*, especially §§1–2, [author-hosted published text](https://www.uwo.ca/apmaths/faculty/jeffrey/pdfs/W-adv-cm.pdf). The coefficient and remainder arguments above are given explicitly rather than inferred from numerical matching.

## 4. The exact combinatorial carrier: a forest with one extra root

The integer \(m!a_m=(m+1)^{m-1}\) counts labelled trees on \(\{0,1,\ldots,m\}\), with the label 0 distinguished. Delete 0 and root each remaining component at its former neighbor of 0. This is a bijection to rooted forests on \(\{1,\ldots,m\}\), with one chosen root in each component. Its inverse joins all component roots to the new vertex 0. Therefore

\[
R(z)=e^{T(z)}
\]

is literally the exponential generating function for a set of rooted labelled trees, including the empty forest.

The operation retains the full forest and its component roots. Forgetting those roots is a quotient and changes the count. The exponential factorials are label normalizations, not probabilities by themselves.

More generally, for any parameter \(s\), Lagrange inversion gives the polynomial identity

\[
[z^m]R(z)^s=\frac{s(m+s)^{m-1}}{m!}\quad(m\ge1).
\]

For positive integer \(s\), the numerator counts forests on \(m+s\) labelled vertices in which each component contains exactly one of \(s\) specified roots. This is the appropriate multi-root extension, rather than a numerological substitution of a tournament count for \(s\).

## 5. The probabilistic bridge: extinction exchanges the two branches

### 5.1 Exact theorem

**Theorem.** Let \(n\ge3\) be an integer and set

\[
\lambda=\log n>1,\qquad q=\frac{r_n}{n},\qquad
\mu=\lambda q=\log r_n\in(0,1).
\]

Start a Galton–Watson process from one individual, each individual having independently a Poisson(\(\lambda\)) number of children. Then:

- its extinction probability is exactly \(q\);
- conditioned on extinction, the entire process is Galton–Watson with Poisson(\(\mu\)) offspring;
- its finite total-progeny probabilities are

\[
p_\lambda(k)=\mathbb P(N=k)
 =e^{-\lambda k}\frac{(\lambda k)^{k-1}}{k!},\quad k\ge1,
\]

and \(\sum_{k\ge1}p_\lambda(k)=q\);
- the conditional distribution satisfies \(p_\lambda(k)/q=p_\mu(k)\), with mean total progeny \((1-\mu)^{-1}\).

These are classical branching identities with an explicit mirror parameterization, not a new general branching theorem.

### 5.2 Extinction and the offspring law

The offspring generating function is \(\phi(s)=e^{\lambda(s-1)}\). Iterating \(q_0=0\), \(q_{j+1}=\phi(q_j)\) gives extinction by successive generations. Its increasing limit is the least fixed point in \([0,1]\).

The mirror identity gives \(\log r_n=(\log n)r_n/n\), so

\[
\log q=\log r_n-\log n=\lambda(q-1).
\]

Since \(1<r_n<e<n\), we have \(0<q<1\). Convexity of \(\phi(s)-s\), together with its two roots \(q,1\), shows this is the least fixed point. Thus it is the extinction probability, not merely an arbitrary algebraic solution of the fixed-point equation.

If the ancestor has \(j\) children, all \(j\) independent descendant processes must become extinct. Hence

\[
\mathbb P(K=j\mid\mathrm{extinction})
 =\frac{e^{-\lambda}\lambda^j q^j}{j!q}
 =e^{-\mu}\frac{\mu^j}{j!}.
\]

Conditional independence of the descendant subtrees follows from conditioning independent subtree events separately. Applying the same calculation recursively proves the full conditioned branching law. The relation \(\mu=\log r_n<1\) proves subcriticality directly.

### 5.3 Total progeny and the coefficients

Let \(F_\lambda(s)=\mathbb E[s^N;N<\infty]\). The branching decomposition gives

\[
F_\lambda(s)=s\exp\{\lambda(F_\lambda(s)-1)\}.
\]

Lagrange inversion at \(s=0\) yields the stated formula for \(p_\lambda(k)\). With
\(z=\lambda e^{-\lambda}\),

\[
p_\lambda(k)=e^{-\lambda}\,a_{k-1}z^{k-1},\qquad
F_\lambda(1)=e^{-\lambda}R(z).
\]

For \(\lambda=\log n\), this reads \(F_\lambda(1)=R(\log n/n)/n=r_n/n\), precisely the extinction probability.

**Normalization boundary.** For \(0<\lambda\le1\), \(R(\lambda e^{-\lambda})=e^\lambda\) and the finite mass is 1. For \(\lambda>1\),

\[
R(\lambda e^{-\lambda})=e^\mu=e^{\lambda q},
\quad \mu e^{-\mu}=\lambda e^{-\lambda},
\]

and the finite mass is \(q<1\). Using \(R=e^\lambda\) indiscriminately would choose the wrong branch and falsely normalize the supercritical process.

Substituting \(q=e^{-\lambda(1-q)}\) in the point-mass formula gives
\(p_\lambda(k)/q=p_\mu(k)\) term by term. Differentiating the subcritical generating-function equation at 1 gives its mean \((1-\mu)^{-1}\).

There is also an exact relation between the elementary iterations. If \(y_0=1\), \(y_{j+1}=e^{zy_j}\), then

\[
y_j=nq_{j+1}\longrightarrow r_n.
\]

Thus the stable fixed-point iteration for the mirror is the extinction-by-generation iteration with its units changed. Its limiting derivative is \(zr_n=\mu<1\).

**Arithmetic corollary.** Among integer \(n\ge3\), the extinction probability \(q\) is rational only at \(n=4\), where \(q=1/2\); for all other such integers it is transcendental. Multiplication by a nonzero integer preserves the rational/transcendental alternatives already proved for \(r_n\). No classification of the arithmetic nature of \(\mu=\log r_n\) is asserted.

**Hostile \(n=2\).** The distinct mirror is 4, but \(4/2=2\) cannot be a probability. Poisson(\(\log2\)) is subcritical and has extinction probability 1. The principal tree series gives 2, not 4, exactly as required for that probability.

Classical primary context: Jagers–Lagerås, [*Supercritical general branching processes conditioned on extinction are subcritical*](https://arxiv.org/abs/0712.1872), which explicitly identifies the simple Galton–Watson case as already known; and Dwass, [*The total progeny in a branching process and a related random walk*](https://doi.org/10.2307/3212112), *Journal of Applied Probability* 6 (1969), 682–686. The specialization used here is fully proved above, so neither a general theorem's hypotheses nor an inaccessible proof is silently imported.

## 6. A scoped MRB-series connection

Set \(z_n=\log n/n\). For \(n\ge3\),

\[
r_n-n^{1/n}=R(z_n)-e^{z_n}
 =z_n^2+\frac52z_n^3+\frac{31}{6}z_n^4+O(z_n^5)>0.
\]

All coefficients after the common constant and linear terms are positive. Furthermore

\[
0<R(z)-e^z\le R(z)-1-z
 \le\frac{3z^2}{2(1-ez)}\quad(0<z<1/e).
\]

For integer \(n\ge3\), \(ez_n\le e\log3/3<1\); therefore
\(\sum_{n\ge3}|r_n-n^{1/n}|<\infty\), by convergence of
\(\sum(\log n)^2/n^2\).

This is a real analytic bridge: replacing \(n^{1/n}-1\) by \(r_n-1\) in the alternating tail preserves convergence and changes its sum by an absolutely convergent correction. It does not preserve the value, prove a relation to a named constant, or establish irrationality/transcendence of either infinite sum. The isolated \(n=2\) term must be handled separately.

## 7. Why this does not explain the tournament gaps 7 and 21

The permanent Hamiltonian-path omissions are already present in repository
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, with the explicit priority correction of 2026-10-20. That file must be cited by slug because THM-1370 has an unrelated legacy identifier collision. Its proof uses multiplication over strong components, insertion, strong-order lower bounds, and small exact spectra. The classical lower-bound input is Busch, [*A Note on the Number of Hamiltonian Paths in Strong Tournaments*](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), not a Lambert-W identity. Completeness of the infinite H-spectrum apart from those two omissions is **not proved here**.

There is a precise forest comparison, but it defeats the proposed direct inference. For the directed triangle \(C_3\), the out-degree Laplacian is

\[
L=\begin{pmatrix}1&-1&0\\0&1&-1\\-1&0&1\end{pmatrix},\qquad
\det(I+L)=7.
\]

Its directed rooted forests are the empty arc set, the three single arcs, and the three two-arc directed paths: exactly 7. Its Hamiltonian paths are the three cyclic starts: exactly 3. Thus a natural forest carrier attains the very integer 7 forbidden to the Hamiltonian-path invariant.

The connection contract is therefore:

- **Source:** labelled rooted forests, or a tournament together with its directed rooted forests.
- **Target:** spanning directed Hamiltonian paths in that tournament.
- **Map:** a Hamiltonian path determines a rooted directed spanning tree by choosing its sink as root; this is an injection into rooted forests.
- **Preserved:** the actual arcs and their directions for each such path.
- **Lost in a scalar forest count:** spanning/one-component/linear-degree constraints; the paths cannot be recovered from the total count alone.
- **Required sidecar:** component count, all local in/out-degree restrictions, and endpoint/root data, rather than just the forest total or tree-series parameter.
- **Decisive hostile:** \(C_3\) gives forest total 7 but path total 3.

The map does not claim that *no conceivable future bridge exists*. It shows exactly why this particular scalar count does not transport the gap predicate. Nor does divisibility by 7 characterize gaps: the script explicitly constructs an order-8 tournament with H=63, independently confirming the existing canonical hostile.

Numerically, \(r_7\approx1.53014011947\) and \(r_{21}\approx1.18794512385\) are ordinary points on the same analytic decreasing branch; no special branch, rationality, or extinction transition occurs at either index. Both are transcendental for the same reason as every other integer input except 2 and 4. The actual critical transition is at the real input \(e\), not at 7 or 21.

## 8. Reproduction and audit boundaries

Run from this directory:

```
python3 verify_mirror_tree.py --output mirror_verification.json
```

The script uses only the Python standard library. It checks:

- the supplied six coefficients, plus the independent exact convolution recurrence through degree 40;
- exhaustive rooted-forest parent functions for sizes 1 through 6, obtaining `1, 3, 16, 125, 1296, 16807`;
- 90-digit bisection versus the convergent series and its analytic tail bound at inputs 2,3,4,5,7,21,63,100,1000;
- extinction fixed points, offspring tilting through 25 children, and the total-progeny identity through size 80;
- the rational integer pairs through 30 as a positive control, not as proof of classification;
- the n=1 and n=2 branch hostiles;
- every labelled tournament through order 6, with exact subset-DP path counts, confirming no H=7 or H=21 in that finite universe;
- the directed-triangle forest/path hostile;
- a concrete H=63 tournament on 8 vertices, with outmasks `[50,252,249,241,224,192,129,1]`.

The exact finite tournament computation is not promoted to an all-order proof; the all-order 7/21 result is an inherited theorem with its own cited mechanism. High-precision decimal residuals are diagnostics, not certificates of transcendence, and not formal interval enclosures. No claim about Collatz, LRC(14), MRB irrationality, or complete tournament H-spectrum realization follows.

**Worktree discipline.** This worker read the bounded startup packet, guardrails, protocol, relevant canonical gap theorem, and relevant correction entries. The coordination-status command reported that `.machine-id` is absent; no registration or repository mutation was attempted. All deliverables are outside the repository in this session packet; publication belongs to the coordinating task.
