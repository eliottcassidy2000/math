# Multiplicand births and a bounded square-deletion inverse

Status: **PROVED** algebraic statements below; **INHERITED** graph definitions
and basic fibre counts; **FINITE-EXACT** controls through 512. No literature
priority claim is made.

## Inheritance and exact object

The older sources use several different objects. **THM-362**,
[`natural-operation-graph-shadows`](../../../../01-canon/theorems/THM-362-natural-operation-graph-shadows.md),
defines the weak loopless shadow as proper divisibility. **THM-2422**,
[`operation-fibres-summand-closure-and-twin-center-ancestry`](../../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
retains the factor-pair involution and its square diagonal. The September
22 [strict multiplicand note](../../collatz_mod6_20260922_lead_reflection2.md)
uses packets $1<a<b$, $ab=n$, contributing both arcs $a\to n,b\to n$.

We use precisely this strict nonunit graph, on vertices $2,3,\ldots,N$.
Its weak version permits $a=b$; neither version uses the unit pair.
The recent [S259 multiplicand packet](../../session259_arithmetic_graph_growth_20261111/multiplicand/RESEARCH.md)
already determines weak and strict-reachability persistence. We inherit that
work rather than rerunning its topology. The present addition retains the
full labelled incidence operator and gives an explicit inverse with uniformly
bounded coefficients.

Live concepts: factor packets, numeric births, square diagonal, divisor
incidence, prime-exponent binary digits, and source-dependent path signs.
The closest mechanism is ordinary Möbius inversion. Hostiles are $4$ and
$9$ as strict composite sources, and $6\to36$ as the first inverse weight
two. The corrected near miss is to read strict reachability as two-parent
generation. The least-used sidecar is the current absolute factor: deleting
$a\to a^2$ is not a condition on the edge ratio alone.

## 1. Adding the next positive integer

At the birth of $n\ge2$, all new loopless edges point into $n$; no old edge
changes. Write $\tau(n)$ for the divisor count and $\delta(n)$ for the
indicator that $n$ is a square. The exact counts are

| Convention | New arcs into $n$ | Unordered nonunit factor packets |
|---|---:|---:|
| Weak, unit present, loops removed | $\tau(n)-1$ | — |
| Weak, unit removed | $\tau(n)-2$ | $(\tau(n)-2+\delta(n))/2$ |
| Strict, unit removed | $\tau(n)-2-\delta(n)$ | $(\tau(n)-2-\delta(n))/2$ |

Divisors pair under $d\leftrightarrow n/d$. The fixed point is exactly
$d=\sqrt n$. Deleting the diagonal packet therefore deletes one shadow arc,
not two. Every strict packet has two different parents and contributes two
different arcs, so the strict birth indegree is always even. The weak
nonunit birth indegree is odd exactly at squares.

A square diagonal is not a self-loop: $a\to a^2$ has different endpoints
for $a\ge2$. If unit multiplication and loops are both retained, the packet
$1\cdot n=n$ supplies the loop $n\to n$ at every integer, not just squares.
This explains why the loop convention must be stated separately.

The weak nonunit sources are the primes. The strict nonunit sources are
the primes and prime squares. For a composite that is not a prime square,
its least prime divisor $p$ satisfies $1<p<n/p$, producing a strict packet.

Let $D(N)=\sum_{n\le N}\tau(n)$ and $s=\lfloor\sqrt N\rfloor$. The total
strict edge count and packet count are respectively

\[
E_{<}(N)=D(N)-2N-s+2,
\]
\[
P_{<}(N)=\sum_{a=1}^{s}\left\lfloor\frac Na\right\rfloor
           -N+1-\frac{s(s+1)}2=\frac{E_{<}(N)}2.
\tag{1}
\]

The first identity sums the birth indegrees. The second counts the choices
$a<b\le N/a$ directly, or uses the exact divisor-hyperbola identity
$D(N)=2\sum_{a\le s}\lfloor N/a\rfloor-s^2$. These are exact all-$N$
growth laws; no average-order estimate is required.

## 2. The retained incidence state and its online inverse

Give each vertex an arbitrary integer weight $f(n)$, with $f(1)=0$.
Read the current vertex together with its strict incoming neighbours:

\[
g(n)=f(n)+\sum_{\substack{2\le d<n\\d\mid n, d^2\ne n}}f(d)
    =\sum_{d\mid n}f(d)-\delta(n)f(\sqrt n).
\tag{2}
\]

This transform preserves much more information than an unweighted graph
statistic. Its finite matrix $W_N=I+B_N$ has rows indexed by targets and
columns by sources; $B_N(n,d)=1$ for a strict arc $d\to n$. It is lower
triangular with diagonal one. At a new integer its exact update is

\[
B_N=\begin{pmatrix}B_{N-1}&0\\v_N^T&0\end{pmatrix},\qquad
W_N^{-1}=\begin{pmatrix}W_{N-1}^{-1}&0\\-v_N^TW_{N-1}^{-1}&1\end{pmatrix}.
\tag{3}
\]

Thus one new incoming row determines one new inverse row; no earlier row
changes. In particular the weights are recovered online by

\[
f(n)=g(n)-\sum_{\substack{2\le d<n\\d\mid n, d^2\ne n}}f(d).
\tag{4}
\]

Every adjacency eigenvalue is zero, and every $W_N$ eigenvalue is one.
These eigenvalues alone cannot distinguish prime births from composite
births. The informative state is the incoming row or the path/response
operator, not the characteristic polynomial.

## 3. A finite dyadic Dirichlet inverse

Work with formal Dirichlet series; no convergence assumption is needed.
Let $F(s)=\sum_{n\ge2}f(n)n^{-s}$ and likewise $G$. Then (2) is exactly

\[
\boxed{G(s)=\zeta(s)F(s)-F(2s)}.
\tag{5}
\]

Here $\zeta$ is the formal divisor-sum series, and its inverse has the
usual Möbius coefficients. Substituting $s\mapsto2s$ repeatedly gives

\[
\boxed{F(s)=\sum_{j\ge0}
 \frac{G(2^j s)}{\zeta(s)\zeta(2s)\cdots\zeta(2^j s)}}.
\tag{6}
\]

This is a coefficientwise finite identity. A contribution at scale $j$
has source index $d^{2^j}$ with $d\ge2$, so it cannot reach target $n$
unless $2^{2^j}\le n$. The remainder after iteration consequently vanishes
at every fixed coefficient. The maximum correction depth at target $n$ is
$\lfloor\log_2\log_2n\rfloor$.

For $a=\prod_p p^{e_p}$ define

\[
\eta_j(a)=
\begin{cases}
(-1)^{\sum_p\operatorname{popcount}(e_p)},&e_p<2^{j+1}\text{ for every }p,\\
0,&\text{otherwise}.
\end{cases}
\tag{7}
\]

Here popcount counts the ones in the binary expansion. Indeed,

\[
\prod_{r=0}^j\zeta(2^r s)^{-1}
=\prod_p\prod_{r=0}^j(1-p^{-2^r s}).
\]

For each prime, choosing a subset of the powers $1,2,\ldots,2^j$ is the
unique binary expansion of an exponent between zero and $2^{j+1}-1$.
Its sign is precisely (7). Thus the explicit inverse matrix is

\[
\boxed{K(d,n)=\sum_{\substack{j\ge0\\d^{2^j}\mid n}}
 \eta_j\!\left(\frac{n}{d^{2^j}}\right)},\qquad
\boxed{f(n)=\sum_{d\mid n, d\ge2}K(d,n)g(d)}.
\tag{8}
\]

The notation $K(d,n)$ puts source first, although it is entry $(n,d)$ of
$W_N^{-1}$. Unlike ordinary Möbius inversion, it is not a function of
$n/d$ alone: $K(2,4)=0$ whereas $K(3,6)=-1$. Both ratios are two, but
only the first arc is a deleted square diagonal.

## 4. Uniformly bounded inverse entries and the exact two-scale law

**PROVED.** For every $2\le d\le n$,

\[
K(d,n)\in\{-2,-1,0,1,2\}.
\tag{9}
\]

There are at most two nonzero terms in (8), and if there are two, their
scale indices are adjacent and $d$ is squarefree. The two terms cancel
when $\omega(d)$ is odd and reinforce when $\omega(d)$ is even.

To prove this, set $a_p=v_p(d)$ and $e_p=v_p(n)$. A scale contributes if
and only if, for every prime,

\[
0\le e_p-2^j a_p<2^{j+1}.
\tag{10}
\]

Equivalently the power $2^j$ lies in the single half-open interval

\[
\max_p\frac{e_p}{a_p+2}<2^j\le
\min_{p:a_p>0}\frac{e_p}{a_p}.
\tag{11}
\]

Any one prime in $d$ bounds the ratio of the interval endpoints by at most
$(a_p+2)/a_p\le3$, so it contains at most two powers of two. If it contains
two adjacent powers $2^j,2^{j+1}$, then
$2^{j+1}a_p<2^j(a_p+2)$ for every $p\mid d$, forcing $a_p=1$.
For these primes write $e_p=2^{j+1}+b_p$ with $0\le b_p<2^j$.
The lower-scale residual exponent is $2^j+b_p$, while the upper-scale
one is $b_p$. Their popcount parities differ by one for each prime in
$d$ and agree for all other primes. This proves the cancellation/sign law
and (9). If the interval is empty, the inverse entry is zero.

An especially simple diagonal identity is

\[
\boxed{K(d,d^2)=1+\mu(d)}.
\tag{12}
\]

Only scales zero and one can occur, contributing $\mu(d)$ and 1.
Consequently a prime square has inverse entry zero, the square of a
squarefree two-prime source has entry two, and a nonsquarefree source has
entry one at its square. For example

\[
K(2,4)=0,\qquad K(4,16)=1,\qquad
K(6,36)=2,\qquad K(6,180)=-2.
\]

This bound concerns individual inverse entries; it does not bound the
number of nonzero terms in a complete inverse row.

## 5. Actual path cancellation and the additive partner

Since $B_N$ is nilpotent,

\[
W_N^{-1}=I-B_N+B_N^2-B_N^3+\cdots.
\tag{13}
\]

Thus $K(d,n)$ is the signed number of strict directed paths from $d$ to
$n$, including the empty path when $d=n$. Formula (9) proves a uniform
cancellation theorem for those paths. At $6\to36$ the direct square arc
was removed, but the two paths through 12 and 18 survive and contribute
$+2$. It is an actual chain count, not a visual comparison.

The companion [additive calculation](../summand/INVERSE.md), developed in
the sibling lane, is

\[
G_+(z)=\frac{F_+(z)}{1-z}-F_+(z^2),\qquad
F_+(z)=\sum_{j\ge0}G_+(z^{2^j})\prod_{r=0}^j(1-z^{2^r}).
\tag{14}
\]

This is an exact operation dictionary. Ordinary cumulative sum is replaced
by divisor sum; doubling an index is replaced by squaring it; binary signs
on integer positions are replaced by binary signs on each prime exponent.
The additive correction has depth $O(\log n)$; the multiplicative correction
has depth $O(\log\log n)$. Equal-parent deletion is the common source of
both formulas. The multiplicative inverse permits coefficients of magnitude
two because several independent prime exponents can reinforce; its additive
partner has only one integer coordinate.

More precisely, define the additive scale kernel

\[
\kappa_j(a,e)=
\begin{cases}
(-1)^{\operatorname{popcount}(e-2^j a)},&0\le e-2^j a<2^{j+1},\\
0,&\text{otherwise}.
\end{cases}
\]

Then the additive inverse is $K_+(a,n)=\sum_j\kappa_j(a,n)$, while our
multiplicative inverse is the synchronized product

\[
\boxed{K(d,n)=\sum_{j\ge0}\prod_{p\mid dn}
 \kappa_j(v_p(d),v_p(n))}.
\tag{14a}
\]

The same scale $j$ must be retained across every prime coordinate.
Summing separately in each coordinate and then multiplying is incorrect:
$K_+(1,2)=(-1)+1=0$, but $K(6,36)=(-1)^2+1^2=2$. This gives an exact
source, map, and lost sidecar for the additive/multiplicative connection:
the exponent-vector map preserves the common dyadic scale; independent
coordinate summation destroys it.

A nonnegative-input example makes that loss observable. Put weight one at
vertex 12 and zero elsewhere. Among divisors of 36, only $g(12)=g(36)=1$
are nonzero. The true inverse has $K(12,36)=-1$, so it recovers
$f(36)=-g(12)+g(36)=0$. The coordinate-product guess instead uses
$K_+(2,2)K_+(1,2)=1\cdot0=0$ and incorrectly returns one. The arc
$12\to36$ is legal even though the exponent of 3 doubles: a deleted square
arc requires every prime coordinate to double at the same time.

## 6. What transfers to factor words and what needs a witness

A path in the **weak** proper-divisor graph from $d$ to $n$ is an ordered
factorization of $n/d$ into factors greater than one. Its length generating
series is therefore the coefficient of $(n/d)^{-s}$ in

\[
\frac1{1-t(\zeta(s)-1)}.
\tag{15}
\]

The strict graph removes a step precisely when its ratio equals the
current vertex. That current vertex is the additional state lost by a
ratio-only word reader.

The old tournament SCC compiler has the same ordered-word mechanism after
replacing each arithmetic factor by an actual strong tournament, retaining
its vertex size and automorphism weight. Its product law transfers through
the **weak** factor graph. Equal numerical factors are allowed in that
compiler: two directed triangles really give $H=9$. Strict arithmetic
diagonal deletion cannot be imported as a prohibition on that construction.
See the [previous multiplicity proof](../../spectrum_multiplicity_signed_tables_20261111/multiplicity/PROOF.md).

Likewise graph reachability is not two-parent generation. Starting with all
primes and requiring two distinct already present factors generates exactly
the primes and the integers with at least two distinct prime divisors; no
pure prime power $p^a$, $a\ge2$, is generated. Yet $p^3$ has incoming arcs
in the full strict graph, because the graph contains the unavailable
companion $p^2$. Adjoining all prime squares repairs this: every $p^a$ for
$a\ge3$ follows by multiplying $p$ by the already available $p^{a-1}$.
Hence primes together with prime squares are the unique least seed set
generating every integer at least two under strict proper multiplication.
Necessity follows from their empty incoming fibres; sufficiency also follows
by induction using a strict factor pair at every other composite.
The sibling [finite-support monoid closure theorem](../summand/CLOSURE.md)
gives the corresponding exact description for arbitrary initial seed sets.

The connection contract is therefore explicit: labelled packets map to an
incoming-incidence observer, which preserves the chosen weighted sum; the
inverse retains the entire target-indexed observation sequence and restores
all input weights. Passing instead to topology, unlabelled support, or a
scalar count discards information required by this reconstruction. None of
the inverse identities alone supplies an absent arithmetic factor or a
missing tournament constructor.

## Reproduction

Run `python3 -B 05-knowledge/results/summand_multiplicand_growth_20261111/multiplicand/check.py`.
The receipt records source hashes, all birth packets through 36, every
inverse/divisor pair through 512, and 130,816 independent lower-triangular
basis equations checking the literal graph matrix times the explicit inverse.
It also checks arbitrary signed-weight recovery and the interval/sign
classification. The proofs establish all-order scope; these bounded controls
corroborate the formulas and their conventions.
