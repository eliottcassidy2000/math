# Code-constrained boxes and the rank-to-corner law

## 1. The retained arithmetic source

Let H be a binary r-by-n matrix of row rank r; C_s={c: Hc=s}. In the
arithmetic examples coordinate i names an actual prime p_i and column i is
the valuation-parity vector of sigma(p_i)=p_i+1. Then selection c denotes
N(c)=product p_i^{c_i}; sigma(N(c)) is square iff Hc=0. Distinct primes give
an injective source map. Replacing the primes by abstract columns preserves
this one predicate, not N, its aliquot successor, its Collatz orbit or its
quadratic class group.

Fix positive widths a_i and independent U_i uniform on[0,1]. Choose c uniformly
from one coset C_s, independently of U, and set X_i=(-1)^{c_i}a_i U_i and
S=sum X_i. This defines a new continuous companion to the arithmetic source;
it does not assert that primes or integer dynamics are sampled this way.
The retained source is (c,U), with c still identifying the integer.

Write A=sum a_i, a(J)=sum_{i in J}a_i, and h(c)=sum_{i:c_i=0}a_i. All
probabilities in this note use this specified measure. For an arbitrary
nonempty subset D of sign words the same construction uses uniform c in D.

## 2. Exact distribution from a weighted enumerator (PROVED)

Define

\[
B_a(z)=\frac1{n!\prod_i a_i}
 \sum_{J\subseteq[n]}(-1)^{|J|}(z-a(J))_+^n.
\]

Then for any nonempty D,

\[
\Pr_D(S>t)=\frac1{|D|}\sum_{c\in D}B_a(h(c)-t).                 (1)
\]

**Proof.** For fixed c, the deficit h(c)-S is a sum of n independent uniforms
on[0,a_i]: replace U_i by1-U_i on positive coordinates and leave U_i on
negative ones. Its distribution is independent of c. Inclusion-exclusion on
the upper faces of the positive orthant simplex gives B_a, the CDF of that
sum. Averaging the exact conditional CDF gives(1). Endpoints have measure
zero. This uses classical box-spline/inclusion-exclusion machinery; the
arithmetic-coset application and source-preserving assembly are the contribution
of this packet, without a literature-priority claim.

Let h*=max_D h(c), M be the number of maximizers, and g the gap from h* to
the next distinct height (infinity if absent). For
0<delta<=min(min_i a_i,g), equation(1) reduces to

\[
\Pr_D(S>h_* -\delta)=
 \frac{M}{|D|}\frac{\delta^n}{n!\prod_i a_i}.                 (2)
\]

This is the first corner: no upper face and no lower-height orthant has yet
contributed. Its validity needs the actual heights and widths, not just a
support count or low moments. Equality at a boundary is harmless by continuity.

For D=C_0, the all-zero word uniquely maximizes h at A. Every other word
loses at least min a_i, so for0<delta<=min a_i,

\[
\Pr_{C_0}(S>A-\delta)
 =2^{r-n}\frac{\delta^n}{n!\prod_i a_i}.                    (3)
\]

For the full independent-sign cube the coefficient is2^{-n}. Thus the ratio
is exactly2^r. This is the **rank-to-corner law**. The extra mass comes from
conditioning, whose probability was2^{-r}; it is not additional physical mass
and is not a free Collatz payment. With nonuniform sign weights, the coefficient
is the conditional mass of the maximizing sign word instead; rank alone no
longer determines it.

For a nonzero syndrome, h*=A-min_{Hc=s}sum a_i c_i. Thus the endpoint displacement
is a minimum-weight syndrome-decoding problem. Multiplicity M is the number of
minimum-weight representatives. Weighted decoding, support onset and corner
mass are exact views of one retained object.

## 3. The arithmetic Hamming family (PROVED)

Take H with all nonzero columns of F2^r, n=2^r-1, r>=2, and widths1. Every
nonzero dual word has weight w=2^{r-1}. Hence all cosets have identical uniform
coordinate marginals of order<w. This is the inherited code theorem, now lifted
by independent magnitudes to continuous uniform[-1,1] coordinate marginals.

Writing u=t/2, the characteristic functions satisfy

\[
\phi_s(t)=\operatorname{sinc}(u)^n\left[
 (\cos u)^n+K_s(i\sin u)^w(\cos u)^{n-w}\right],
\quad K_0=2^r-1,\quad K_s=-1\ (s\ne0).                    (4)
\]

**Proof.** Integrating the independent magnitudes gives
sinc(u)^n times the sign-word average of exp(iu sum(-1)^{c_i}). Expand each
factor as cos u+i(-1)^{c_i}sin u. A binary character has coset mean zero unless
its exponent is in row(H); for exponent bH that mean is(-1)^{b dot s}.
All nonzero bH have weight w, and sum_{b!=0}(-1)^{b dot s} is2^r-1 or-1.
This proves(4), including its normalization at t=0.

Consequently all scalar moments below degree w agree, and their first difference
is exactly

\[
\mathbb E_{C_0}S^w-\mathbb E_{C_s}S^w=2^{r-w}w!\quad(s\ne0). (5)
\]

For r=3 this difference is12 at degree4; for r=4 it is2520 at degree8.
The finite prime banks give unconditional arithmetic realizations at r=3,4.
The abstract binary-code identity holds for every r; general prime realization
is a separate open prime-values question, as Session85 explicitly notes.

For r=3 the primes in column order1,...,7 are

    7, 2, 5, 19, 89, 59, 29.

Both continuous cosets have mean0 and variance7/3, but their fourth moments
are259/10 and139/10. At t=13/2 their upper-tail masses are respectively
1/10321920 and0. Indeed C0 has endpoint7, unique maximum; every nonzero coset
has a unique weight-one representative, hence endpoint6. More generally,
for0<delta<=1 the tails at their respective endpoints are both

    delta^7 / (16 * 7!).

So the same marginal data, normalization and local tail shape coexist with a
one-unit difference in support. The inherited discrete marginal obstruction has
become an exact continuous obstruction to inferring a source tail from low-order
statistics.

## 4. Relations on trees, not graphic cycle spaces (CITED + direct proof)

At a coordinate cut E=L disjoint-union R, put V_L=span(H_i:i in L) and V_R
similarly. Among partial words that extend to a zero-syndrome word, the interface
state is H_L c_L in V_L intersect V_R. Its minimal number of binary dimensions is

    lambda(L)=rank(H_L)+rank(H_R)-rank(H).                   (6)

Every intersection value occurs. Two distinct values have disjoint right
completion sets; equal values have identical completion sets. Thus there are
exactly2^lambda distinct continuation languages. This is the classical minimal
trellis/tree-state mechanism, not a new coding theorem: see Kashyap,
[On Minimal Tree Realizations of Linear Codes](https://arxiv.org/abs/0711.1383),
and the incoming Session85 LCA/scope-loss theorem.

For multiplicative coordinate weights f_i(c_i), transmit the actual sums of
products indexed by these states and combine compatible states once. This
computes the weighted enumerator and(4) exactly. Complex Fourier weights are
allowed because this is a finite algebraic identity; positivity must not be
inferred from them. The values remain functions of t, not constant-size bit
strings. Under arbitrary non-product source weights, extra correlations must
be retained; rank alone does not bound that cost.

A Tanner/tree factorization can therefore represent the Fano code despite its
nongraphic matroid. What fails is the much narrower request that its seven
columns be exactly an ordinary graph's edge-cycle dependencies. Local relations,
shared variables and weights supply the missing structure. For composite atoms,
the coprimality conflict factors must be included as well.

## 5. Borwein is the unconditional box control (CITED)

With independent X_j uniform[-1/(2j+1),1/(2j+1)], j=1,...,m, let S=sum X_j.
Then the published Fourier identity gives

    I_m(t)=integral_0^infty (sin(tx)/x) product_j sinc(x/(2j+1)) dx
          =pi/2 - pi P(S>t).

If epsilon=sum_j 1/(2j+1)-t is positive and no more than2/(2m+1), the corner
formula gives

    pi/2-I_m(t)=pi epsilon^m /(2^m m! product_j (1/(2j+1))).

At the first crossing this condition holds. For t=1 the first bad m is7
(denominator15); t=2, equivalently the2cos(x) variant, first fails at m=56
(denominator113). Thus seven integrals indexed m=0,...,6 and fifty-six indexed
m=0,...,55 have the plateau value. The last good index is not the first bad one.

The sinc-sum/integral identity has a different threshold: the first nonzero
Poisson image at2pi, with first bad m=40249 and denominator80499. It is not an
80499-term plateau. These are published results in BB2001/BBB2008/AG2014 and
BB2016, not discoveries here. THM-4837's corrected audit supplies the current
repo route. The shared mechanism with(1) is positive-part box volume; no theorem
identifies the harmonic crossing index56 with an E7 representation dimension.

## 6. An inverse transform: small asymmetry can recover the whole sign law (PROVED)

Let pi(c) be any probability law on{0,1}^n, retaining independent uniform
magnitudes. Formula(1) works with coefficients pi(c). With Z=sum a_i U_i,

    phi_S(t)=phi_Z(-t) * sum_c pi(c) exp(i t h(c)).           (7)

The known factor phi_Z is nonzero in a neighbourhood of0. Thus the full scalar
law of S determines exactly the discrete measure sum pi(c) delta_{h(c)}.
To see uniqueness without a deconvolution assumption at its zeros, divide near0:
finite exponential sums with distinct frequencies are linearly independent,
for example by their Vandermonde derivative matrix.

It follows that **the scalar law determines every sign-word probability iff
all subset sums of the widths are distinct**. If heights collide, transferring
mass between two colliding sign words changes no scalar law; if they do not,
all masses are recovered. This is a complete inverse statement for this
construction, not merely a sufficient genericity claim.

An arbitrarily small perturbation of equal widths suffices at every fixed n:
choose a_i=1+epsilon*2^i for i=0,...,n-1 with
0<epsilon<2^(-n-1). Different cardinalities have height difference greater
than1/2; equal cardinalities have different binary sums. Choosing epsilon even
smaller keeps every width within any prescribed neighbourhood of1.

This is **informationally exact, not numerically free**. Separations can be
of order epsilon; the distribution is an infinite-precision object. At equal
widths, only the Hamming-weight distribution survives, and all nonzero Hamming
syndromes give the same scalar law. The construction therefore exhibits both
loss and a precisely priced way of removing it. It suggests active probe
design rather than passive accumulation of more coarse invariants.
