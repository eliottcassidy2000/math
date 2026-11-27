# Sharp weighted switch samples and the Hadamard readout

2026-11-27. The [marked-pants switching theorem](../balanced_pants_switch_fourier_20261127/README.md)
gives universal averages on half the cut switches. Half is
optimal even if the sample may be nonlinear and each switch
may be given an arbitrary positive weight. A single
Gram-measuring profile is much cheaper, while universal mixed
second moments need every switch.

A universal sampling law here is fixed independently of the labelled
starting tournament. It may have an arbitrary nonempty support and
nonnegative weights summing to one.

## 1. Half is necessary

Write \(G=\{\pm1\}^n/\{\pm\mathbf1\}\), \(N=|G|=2^{n-1}\),
and identify its characters with even labelled sets. For \(n=2\)
and every \(n\ge4\), the minimum support for the universal
Hamiltonian-path mean is \(N/2\). At \(n=1,3\) the minimum is
\(N\). A law working for Hamilton paths automatically works for
every ordered path-cover profile.

To prove this let \(f=N\mu\) be the sampling density and
\(g=f-1\). The sharp Fourier theorem says that averaging the
Hamiltonian-path function must kill every character of nonzero
even weight at most \(2\lfloor n/3\rfloor\), and only those.
For necessity apply the putative law to every switched
translation of a tournament whose Fourier coefficient does
not vanish at that character; orthogonality on the translate
forces the corresponding coefficient of \(\mu\) to vanish.
Every other path-cover profile has support within the same
cutoff. Thus \(g\) is mean zero and has support only at
weights at least \(d=2\lfloor n/3\rfloor+2\).
For \(n=1,3\) there are no allowed nonzero characters.

For distinct high-weight sets \(S,T\), their symmetric
difference has weight at most
\(2(n-d)\le d-2\). Therefore \(g^2\) has no high Fourier
coefficient, and under uniform expectation \(E[g^3]=0\).
Put \(v=E[g^2]\). Under the probability measure \(\mu\),
\[
 E_\mu[g]=v,\qquad E_\mu[g^2]=v,
\]
so variance nonnegativity gives \(v\le1\).
By Cauchy–Schwarz
\[
 |\operatorname{supp}\mu|\ \ge\ \frac{N}{E[f^2]}\ \ge\ N/2.
\]
A coset of the kernel of any high character attains equality;
such a character exists at \(n=2\) and all \(n\ge4\).
At equality both inequalities above are equalities, so
\(f=2\) on the support and \(f=0\) elsewhere. Thus an optimal
law is necessarily uniform on its support, though that
support need not be a coset. In particular a universal
*subgroup* cannot have codimension larger than one:
its annihilator is a subspace and two distinct nonzero high
characters sum to a low character.

## 2. Classifying optimal samples and constructing nonlinear ones

Equality means \(g\) takes values \(\pm1\). On the full sign cube,
this is exactly
\[
 g(x)=\left(\prod_{i=1}^n x_i\right) h(x)
\]
where \(h\) is a \(\{\pm1\}\)-valued function of Fourier degree at
most \(n-d\), satisfying \(h(-x)=(-1)^n h(x)\). Indeed,
multiplication by full parity complements the support sets;
the parity condition ensures \(g\) descends to the switch
group. Conversely these conditions force zero constant
coefficient for \(g\) and hence an optimal law uniform on
\(g=1\).

If \(n-d\le1\), the Boolean \(h\) must be constant or a signed
single coordinate (expand \(h^2=1\)). All optima there are
affine hyperplanes in \(G\). The first nonlinear examples
are \(n=8\) for even order and \(n=11\) for odd order:

* For even \(n\ge8\), take
  \(h=(x_1x_2+x_1x_3+x_4x_2-x_4x_3)/2\).
  This equals \(x_1x_2\) when \(x_2=x_3\), otherwise \(x_4x_2\).
* For odd \(n\ge11\), take
  \(h=(x_1+x_2+x_3-x_1x_2x_3)/2\), the majority sign on
  the three coordinates.

The minimum degrees required are two and three respectively,
and remain allowed in each larger order of the given parity.
In both cases there are four nonzero Fourier coefficients of
\(g\), not the single coefficient of an affine sample.
The replay checks the 64- and 512-state examples exactly.

## 3. A Hadamard-sized sample for the Gram profile

For \(\ell=(3,1,\ldots,1)\), the sharp support criterion is only
\(E_\mu[x_ix_j]=0\) for every \(i\ne j\). Choose representatives
with \(x_n=1\). This is strength-two weighted binary
orthogonality on the other \(n-1\) columns, the analogue of
the [Hadamard–orthogonal-array correspondence](https://arxiv.org/abs/2204.00835).

The \(n\) sign columns are orthonormal in \(L^2(\mu)\), so
at least \(n\) support states are needed. Equality is possible
if and only if a real Hadamard matrix of order \(n\) exists:
if there are \(n\) states then the square matrix of row
vectors scaled by \(\sqrt{\mu(x)}\) has orthonormal rows as
well as columns. Each row norm gives \(\mu(x)=1/n\).
Conversely normalize one column of a Hadamard matrix to
\(+1\) and take the rows uniformly. For arbitrary \(n\ge3\),
Sylvester doubling gives a sample of size at most
\(2^{\lceil\log_2 n\rceil}\): retain the constant and sufficient
independent generating columns to keep the rows distinct.
At \(n=4\) the two minimum sample sizes are both four;
at \(n=8\) they are exactly 8 versus 64.

More generally, for a profile of cutoff \(2t\), using \(x_n=1\)
turns universality into weighted binary strength \(2t\) on
\(n-1\) signs: a character of odd size there lifts by adjoining
\(n\). Orthogonality of all monomials of degree at most \(t\)
gives the weighted Rao bound
\(\sum_{i=0}^t\binom{n-1}{i}\) on the support size.
This is another way to see the Hadamard lower bound at \(t=1\).

## 4. The cost of a mixed second moment is the full switch group

Suppose the law must reproduce the uniform-group mean of
\[
 H(T^\varepsilon)\,H(U^\varepsilon)
\]
for **every pair** of labelled starting tournaments \(T,U\).
For \(n\ge3\), the only such law is uniform on \(G\). This
asserts a mixed-product result, not that a single-tournament
variance alone forces the same cost.

More generally consider an ordered list of path-cover profiles
whose cutoffs are \(D_1,\ldots,D_k\), and evaluate their counts
on independently chosen starting tournaments (all using the
same random switch). A law reproduces the uniform-group mean
of every such product exactly when it kills every nonzero
even character of weight at most \(\sum_iD_i\). For sufficiency,
the Fourier support of the product is contained in that range.
For necessity, split a target even set into disjoint even sets
\(S_i\) of sizes at most \(D_i\); empty sets are allowed.
For each \(S_i\) choose a starting tournament with nonzero
coefficient there, or use the fixed positive mean for an
empty set. Translate the \(k\) starting tournaments independently
by \(y_1,\ldots,y_k\in G\). The coefficient of
\(\prod_i\chi_{S_i}(y_i)\) in the difference of the two tested
means is the product of those nonzero coefficients times
\(\widehat\mu(\triangle_i S_i)\). It must vanish.
Taking two Hamiltonian profiles covers every even set since
\(4\lfloor n/3\rfloor\ge 2\lfloor n/2\rfloor\) for \(n\ge3\).
This gives three concrete experimental costs: a Hadamard-sized
sample for the Gram profile, half for the universal path mean,
and all switches for mixed Hamiltonian second moments.

Run the companion Node.js replay for the finite checks.
