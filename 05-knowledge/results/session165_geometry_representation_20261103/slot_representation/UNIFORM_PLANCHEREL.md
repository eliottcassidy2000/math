# Every growing tournament has an asymptotically Plancherel slot module

**Status.** General theorem proved by cutting and joining path covers. It strengthens the scope of [the slot-representation note](NOTE.md): arbitrary tournament sequences have a uniform quadratic rate, while repeated C-moves retain the stronger exponential rate proved there. No new tournament-count spectrum claim is made.

## 1. A path-cover cut inequality

Let X be any directed graph on n≥1 vertices, and let p_k count its unordered covers by k nonempty directed paths. Put p_0=0 and p_k=0 for k>n.

**Lemma.** For 1≤k≤n,

\[
 (n-k)p_k\le k(k+1)p_{k+1}.                         \tag{1}
\]

**Proof.** Every k-path cover contains exactly n−k path arcs. Choose one arc and delete it, obtaining a (k+1)-path cover together with the ordered pair of paths that should be rejoined. Conversely, a fixed (k+1)-path cover has at most k(k+1) ordered pairs of distinct paths; only pairs whose last-to-first arc exists can be rejoined. Counting this same relation from its two sides proves (1). ∎

No tournament property is used here. The possible joins need not all exist, which is why the statement is an inequality. More precisely, the nonnegative deficit

\[
 \delta_k=k(k+1)p_{k+1}-(n-k)p_k
\]

counts a (k+1)-cover together with an ordered pair of distinct path blocks A,B for which the arc end(A)→start(B) is absent. Thus the right side of (3) is itself a positive cover-defect invariant.

Iteration gives a concrete lower bound in terms of Hamilton paths H=p_1:

\[
 p_k\ge \frac{H}{k!}\binom{n-1}{k-1}.               \tag{2}
\]

This is meaningful whenever H>0, but remains true when H=0.

## 2. A stronger normalized-slot inequality

Define the ordered-slot polynomial

\[
 F(r)=\sum_kp_k(r)_k.
\]

The cut deficits in (1) give the exact polynomial identity

\[
 rF(r+1)-(n+r)F(r)
 =\sum_{k\ge1}
    \big[k(k+1)p_{k+1}-(n-k)p_k\big](r)_k.           \tag{3}
\]

To verify it, use

\[
 (r+1)_k-(r)_k=k(r)_{k-1},\qquad
 r(r)_{k-1}=(r)_k+(k-1)(r)_{k-1}.
\]

Thus the right side has **nonnegative falling-factorial coefficients**. In particular, at each positive integer r,

\[
 rF(r+1)\ge(n+r)F(r).                              \tag{4}
\]

Equivalently,

\[
 \frac{F(r)}{\binom{n+r-1}{r-1}}
\]

is nondecreasing in the positive integer r. Whenever 1≤s≤r and F(r)>0,

\[
 \frac{F(s)}{F(r)}
 \le
 \frac{\binom{n+s-1}{s-1}}
      {\binom{n+r-1}{r-1}}.                        \tag{5}
\]

Taking s=1 recovers the useful bound

\[
 \boxed{\quad F(r)\ge H\binom{n+r-1}{r-1}.\quad}    \tag{6}
\]

This also follows from (2) and Vandermonde's identity:

\[
 \sum_k\binom{n-1}{k-1}\binom rk
 =\binom{n+r-1}{r-1}.
\]

The complete bidirected graph attains equality: its path-cover numbers are the unsigned Lah numbers

\[
 p_k=\frac{n!}{k!}\binom{n-1}{k-1},\qquad
 F(r)=n!\binom{n+r-1}{r-1}.
\]

Every ordered pair of blocks can be joined, so every deficit is zero. Thus (5) is sharp over general digraphs. The statement concerns slot counts on one fixed graph as the number of slots changes. It is not an assertion about a physically constructed new graph.

## 3. Universal Plancherel convergence

Let X be an n-vertex tournament, let r≥3 be fixed, and let V_X(r) be its complex slot-cover permutation representation. Write

\[
 D=\dim V_X(r)=F_X(r),\quad
 m_\lambda=\text{multiplicity of }S^\lambda,\quad
 d_\lambda=\dim S^\lambda.
\]

Here D>0: every tournament has a Hamilton path. For a nonidentity permutation σ∈S_r, the fixed-slot count is at most r−2. Since the character is F_X(fix(σ)), (5) gives

\[
 0\le\frac{\chi_{V_X(r)}(\sigma)}{D}
 \le\frac{F_X(r-2)}{F_X(r)}
 \le \varepsilon_{n,r},                            \tag{7}
\]

where

\[
 \boxed{\quad
 \varepsilon_{n,r}
 =\frac{(r-1)(r-2)}{(n+r-1)(n+r-2)}.
 \quad}                                             \tag{8}
\]

Derangements have zero character. The character inner product and |χ_λ(σ)|≤d_λ now yield the explicit bound

\[
 \left|\frac{m_\lambda}{D}-\frac{d_\lambda}{r!}\right|
 \le d_\lambda\left(1-\frac1{r!}\right)
       \varepsilon_{n,r}.                         \tag{9}
\]

Hence the probability that a uniformly chosen representation dimension lies in the λ-isotypic component satisfies

\[
 \left|\frac{d_\lambda m_\lambda}{D}
             -\frac{d_\lambda^2}{r!}\right|
 \le d_\lambda^2\left(1-\frac1{r!}\right)
       \varepsilon_{n,r}.                         \tag{10}
\]

**Theorem.** For every sequence of tournaments X_n with |V(X_n)|→∞, the fixed-r isotypic dimension distribution of V_(X_n)(r) converges to the Plancherel distribution of S_r. The error is O_r(n^−2), uniformly over all tournaments of order n.

There is no strength, guard, recursion, randomness, or profile-boundedness assumption. The same conclusion holds for any sequence of digraphs whose r-slot modules are nonzero, because (1) and the fixed-slot character formula hold there too.

The r=2 module is exactly a sum of regular representations for every nonempty tournament: its only nonidentity permutation is a derangement. The r=1 case is tautological.

### A particularly transparent S3 consequence

The exact multiplicities are

\[
 (m_{(3)},m_{(2,1)},m_{(1,1,1)})
 =\left(\frac{F(3)+3H}{6},\frac{F(3)}3,
              \frac{F(3)-3H}{6}\right).
\]

The isotypic dimension distribution is therefore

\[
 \left(\frac16+\frac{H}{2F(3)},\frac23,
              \frac16-\frac{H}{2F(3)}\right).
\]

Its total variation distance from S3 Plancherel is exactly H/(2F(3)), and (6) gives

\[
 d_{\mathrm{TV}}\le\frac1{(n+1)(n+2)}.               \tag{11}
\]

Thus all tournaments have the standard S3 sector at exactly 2/3 of their 3-slot dimension, with the trivial/sign imbalance uniformly shrinking at least quadratically.

## 4. Quadratic order is optimal even for tournaments

Choose a labeled tournament uniformly at random. The complete bidirected graph has

\[
 \frac{n!}{k!}\binom{n-1}{k-1}
\]

unordered k-path covers: order all vertices, choose k−1 cuts, and divide by the k! orders of the blocks. A specified cover requires n−k distinct arcs, whose orientations occur independently with probability 2^−(n−k). Hence

\[
 \mathbb E[p_k]
 =\frac{n!}{k!2^{n-k}}\binom{n-1}{k-1}.             \tag{12}
\]

In particular,

\[
 \mathbb E[F(3)]=(2n^2+1)\mathbb E[H].              \tag{13}
\]

To check the coefficient, divide the k=1,2,3 contributions by E[H]: they are 3, 6(n−1), and 4 binom(n−1,2), summing to 2n²+1.

The ratio E[H]/E[F(3)] is the F(3)-weighted average of the positive ratios H/F(3). Therefore, for every n, some n-vertex tournament satisfies

\[
 \frac{H}{F(3)}\ge\frac1{2n^2+1}.
\]

Its S3 isotypic distribution has

\[
 d_{\mathrm{TV}}\ge\frac1{4n^2+2}.                 \tag{14}
\]

Combining (11) and (14) gives the all-order extremal estimate

\[
 \boxed{\quad
 \frac1{4n^2+2}
 \le \max_{|T|=n}d_{\mathrm{TV}}(V_T(3),\mathrm{Plancherel})
 \le \frac1{(n+1)(n+2)}.
 \quad}                                             \tag{15}
\]

Thus the uniform quadratic order cannot be improved over tournaments, already at rank three. The lower bound is an existence theorem from averaging, not an explicit tournament family, and neither extremal constant is claimed optimal.

### The same optimal order holds among strong tournaments

A nonstrong tournament admits a nonempty proper source cut A→B. If |A|=k, then on this event every Hamilton path consists of an A-path followed by a B-path, so H(T)=H(A)H(B). The cross-cut event has probability 2^−k(n−k), independently of the two induced tournaments. Summing over possible source cuts gives the weighted union bound

\[
 \frac{\mathbb E[H\,\mathbf1_{\mathrm{nonstrong}}]}{\mathbb E[H]}
 \le 2\sum_{k=1}^{n-1}2^{-k(n-k)}
 \le(n-1)2^{2-n}.                                    \tag{16}
\]

For the first inequality, the k-th contribution before division by E[H] is

\[
 \binom nk2^{-k(n-k)}
       \frac{k!}{2^{k-1}}\frac{(n-k)!}{2^{n-k-1}}.
\]

Set α_n=1−(n−1)2^(2−n), which is positive for n≥4 and tends to 1. Since E[F(3)1_strong]≤E[F(3)], weighted averaging restricted to strong tournaments gives

\[
 \boxed{\quad
 \frac{\alpha_n}{4n^2+2}
 \le\max_{\substack{|T|=n\\T\ \mathrm{strong}}}
       d_{\mathrm{TV}}(V_T(3),\mathrm{Plancherel})
 \le\frac1{(n+1)(n+2)},\qquad n\ge4.
 \quad}                                               \tag{17}
\]

Thus the quadratic rate is optimal in order even for strong tournament sources. This is a probabilistic existence proof, not an identification of extremizers.

## 5. Relation to literal C-growth

The main note's C-chain theorem gives

\[
 O\!\left((\rho(r-2)/\rho(r))^j\right)
\]

after j one-vertex moves. It remains a genuine strengthening: on that physical family the universal quadratic convergence accelerates to a computable exponential rate. The two theorems answer different questions, uniformity over all tournaments versus the sharper dynamics of a specific construction.

## 6. Finite controls

Run `python3 check_uniform.py` and `python3 -O check_uniform.py`. The checker independently enumerates path covers by permutations and cuts, rather than the companion checker’s subset dynamic programming. It verifies the cut deficits, all displayed slot inequalities at ranks 1 through 8, and the exact S3 error on all 1,099 labeled tournaments through order five. It also checks all 64 loopless directed graphs on three vertices, five bidirected equality cases, and the random-tournament expectation and extremal-existence inequalities at every tested order. The general cut proof does not require the tournament restriction. The finite counts are controls; the double-counting and character arguments prove the all-order theorem.
