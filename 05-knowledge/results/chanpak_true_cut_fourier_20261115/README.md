# Fourier analysis on the actual (sharp-cut) Chan–Pak carrier domain

**Proved for the displayed carrier family.** These are bounds for
**presentations** in a specified region, not distinct-count asymptotics
or pointwise coverage. The [168-form conjecture](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/hypotheses/HYP-9722-finite-carrier-conjecture-every-odd-n-from-333-is-a-value-of-one-of-168-digit-quadratic-carrier-forms.md)
remains open.

Write \(e(u)=\exp(2\pi i u)\), and let
\[
p(0)=0,\quad p(2u+1)=p(u),\quad p(2u)=p(u)+u+(u\bmod2)
\]
(the even equation is used for positive input). The independent
background estimate in §1 retains the reversed-bit interaction;
the carrier estimate in §2 instead uses the existing \(xy\) term
and keeps the **true low-digit cutoff** by interval completion.
At fixed rational centers §3 computes the actual arithmetic factor
and extends it to a bounded neighborhood.

## 1. A one-tape uniform Fourier estimate

For an interval \(I\) of \(N\) nonnegative integers, suppose
\(\theta=2\alpha=a/q+O_K(q^{-2})\) with \((a,q)=1\) and the error in
absolute value at most \(K/q^2\) for fixed \(K\). Then uniformly in
\(\beta\in\mathbb R\),
\[
\frac{1}{N}\left|\sum_{n\in I}e(\alpha p(n)+\beta n)\right|
\ll_K
\left(\frac1q+\frac{\log(2q)}{\sqrt N}+
      \frac{q\log(2q)}{N}\right)^{1/2}.
\]
On \(N^\epsilon\le q\le N^{1-\epsilon}\) this saves
\(N^{-\epsilon/2}\) up to a square root of a logarithm.

**Bilinear bound.** For intervals of respective lengths \(U,V\),
and arbitrary weights \(|a_u|,|b_v|\le1\), let
\(S=\sum_{u,v}a_u b_v e(\theta uv)\). Cauchy in \(u\) and geometric
summation give
\[
|S/(UV)|^2 \ll \frac1V+
\frac1{UV}\sum_{1\le d<V}\min(U,\|\theta d\|^{-1}).
\]
When \(|\theta-a/q|\le K/q^2\), divide the \(d\)'s into blocks
of length comparable to \(q/(2K+2)\), rounded down but at least
one. Distinct values in each block have circular distance at
least \(1/(2q)\) by the reduced rational and the perturbation
bound. Sorting by distance from zero bounds each block sum
by \(O_K(U+q\log(2q))\); small \(q\) use the trivial bound.
Consequently
\[
|S/(UV)|^2\ll_K
\frac1V+\frac1q+\frac{\log(2q)}U+
\frac{q\log(2q)}{UV}.                                      \tag{1}
\]
This works in translated intervals, without a restriction on the
origin of either diagonal twist.

For the one-tape sum, take a power of two \(V=2^r\asymp\sqrt N\),
write \(n=Vu+b\), and discard at most two boundary rows at cost
\(O(V)\). Use the exact
[binary split](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md)
\[
p(Vu+b)=p(u)+p(b)+(V-1-\operatorname{rev}_r b)\,u+
 [b<V/2](u\bmod2).
\]
It also follows by appending the \(r\) binary suffix bits using
the recurrence. After splitting \(u=2v+d\) by parity and replacing
\(b\) by its bit reversal, the only nonseparable phase is
\(-2\alpha v\operatorname{rev}_r b\). The remaining factors, including
the linear twist, are unit weights. Apply (1) with both lengths
comparable to \(\sqrt N\).

The generic bilinear bound cannot promise better than \(q^{-1/2}\)
in the first term. At \(\theta=1/q\) for odd \(q\), take both lengths
multiples of \(q\) and weights
\(a_u=e(\overline2 u^2/q)\), \(b_v=e(\overline2 v^2/q)\), with
\(2\overline2\equiv1\pmod q\). The phase becomes
\(e(\overline2(u+v)^2/q)\), and the complete sum has modulus exactly
\(UV/\sqrt q\). The one-variable Gauss sum has squared modulus \(q\):
change variables to sum and difference (invertible for odd \(q\)),
then use additive orthogonality. This is a limitation of ignoring
the structure of the weights, not a lower bound for the checksum
weights themselves.

## 2. Minor arcs in a true carrier sublevel or shell

For a form in the [published bank](../chanpak_guard_depth_20261115/forms_orders3_5.txt)
write
\[
F(x,y)=Qxy+tx+sy+a(p(x)+p(y))+n+f[x\text{ odd}]+g[y\text{ odd}].
\]
All 168 forms have \(a>0,\ Q\ge4a,\ 0\le f,g\le a\), with nonnegative
coefficients everywhere. These are the only coefficient inequalities
used below. Fix \(0<\eta<\delta<1/2\) and \(D=\lceil T^\delta\rceil\).
If \(\alpha=b/q+O(q^{-2})\) with reduced \(b/q\), absolute error at
most \(q^{-2}\), and
\[
T^\eta\le q\le T^{1-\eta},
\]
then uniformly in \(0\le L\le T\),
\[
\left|\sum_{\substack{x,y\ge D\\L<F(x,y)\le T}}
 e(\alpha F(x,y))\right|
\ll_{F,\eta,\delta} T^{1-\eta/2}(\log(2T))^{5/2}.             \tag{2}
\]
Separable bounded tapewise weights may be inserted. The excluded
strips \(x<D\) or \(y<D\) are **not** estimated by (2).

**The cutoff is an interval in the shorter coordinate.**
The digit recurrence gives \(p(2u+1)-p(2u)=-u-(u\bmod2)\).
For an odd-to-even step the same recurrence reduces to an increment
at roughly half the argument plus \(u+1+((u+1)\bmod2)\).
Induction therefore gives \(p(y+1)-p(y)\ge-y/2-1\) for all \(y\ge0\).
If \(1\le y\le x\), then
\[
F(x,y+1)-F(x,y)\ge Qx+s-a(y/2+1)-g>0.
\]
Thus in \(x\ge y\), fixing \(x\) cuts a genuine interval out of the
short tape \(y\), even with both bounds of a shell imposed. For
\(y>x\) exchange the roles. No smoothness of the cut in the **long**
coordinate is assumed.

**Complete each such interval on its dyadic box.** Divide into
dyadic boxes, with side lengths \(U,V\) and the *long* tape as the
first coordinate. For each \(u\), expand the indicator of its cut
interval of the \(V\)-long tape on the cyclic grid of length \(V\).
The Fourier coefficients \(c_u(j)\), in the \(1/V\)-normalized
convention, satisfy
\[
|c_u(j)|\ll (1+\min(j,V-j))^{-1},\qquad 0\le j<V.
\]
This is just the geometric-sum formula (the zero coefficient is
at most one). Divide each coefficient by its displayed uniform
bound and include it in the first tape's weight in (1); the character
of frequency \(j\) is absorbed into the second tape's weight. The
completion costs \(O(\log(2V))\). The original phase was already
separable outside \(Q\alpha xy\), so there is no independent
smooth-checksum replacement. If \(\alpha\) has the approximant
above, \(Q\alpha\) has reduced denominator \(q'=q/(q,Q)\) and an
error at most \(Q/(q')^2\), so (1) applies.

Each nonempty box has \(U,V\gg D\) and \(UV\ll_F T\), since
\(F\ge Qxy\). For every fixed \(U\), sum the permissible \(V\)'s
geometrically; overall \(\sum UV\ll T\log(2T)\) and
\(\sum\sqrt{UV}\ll\sqrt T\log(2T)\).
In (1), the first three square-root terms times \(UV\) are bounded
after summation by \(O(T^{1-\eta/2}(\log(2T))^{3/2})\).
The fourth term is \(\sqrt{UVq'\log(2q')}\), which meets the same
total bound by \(q'\le T^{1-\eta}\). Add the completion cost,
and the two halves together give (2).

## 3. Fixed-modulus factors and bounded major arcs

Let \(M\) be fixed and set \(d=(a,M)\). Remove the two checksum
terms and call the remaining polynomial-with-parities \(B(x,y)\),
so \(F=B+2a(p(x)/2+p(y)/2)\). For an odd target residue \(b\bmod2M\)
define
\[
\sigma_M(b)=\frac{d}{M}\frac{
 \#\{u,v\bmod 2d:B(u,v)\equiv b\pmod{2d}\}}{(2d)^2}.        \tag{3}
\]
For fixed \(M,F,\delta\) there is \(\kappa>0\) such that, uniformly
in \(0\le W\le T\), the number of presentations with
\(x,y\ge D,\ F(x,y)\le W\), and \(F(x,y)\equiv b\pmod{2M}\) is
\[
\sigma_M(b)\,\#\{x,y\ge D:Qxy\le W\}
+O_{M,F,\delta}(T^{1-\delta\kappa}).                          \tag{4}
\]
In (4) the modulus is genuinely fixed. This is not a bound with
a constant uniform in \(M\).

Here is the inherited mixing input and its application. Appending
a fair bit acts on states \((x\bmod2M,p(x)/2\bmod M)\).
Each state has two incoming transitions, since the last bit is
determined by the parity of the new first coordinate, the two
possible previous first coordinates differ by \(M\), and each
determines the previous checksum uniquely. The invariant law is
uniform. The [fixed-length suffix programmer](../chanpak_absolute_steering_20261115/PREFIX_PROGRAMMING.md)
makes a fixed power of the transition matrix strictly positive,
so it converges exponentially by the elementary uniform
minorization. Consequently a bounded state-dependent test on
any aligned dyadic interval of length \(2^k\) has discrepancy
\(O_M(2^{k(1-\kappa)})\). An arbitrary interval is a disjoint
union of at most two aligned dyadic intervals per scale; sum this
geometrically (take \(0<\kappa<1\)).

First, on \(x,y\ge D\), one has \(0\le p(u)\le u\) by induction;
for odd \(u\), the stronger \(p(u)\le(u-1)/2\) pays for the parity
increment of \(p(2u)\). Thus \(Qxy\le F(x,y)\le Qxy(1+O_F(1/D))\).
Passing between the two cutoff domains costs
\(O_F(T\log(2T)/D)\) by summing the hyperbolic boundary strip over
\(x\le T/(QD)\), including floor errors.
On the hyperbolic domain, for each \(x\) the \(y\)-test is a function
of its modular state only. Mix over that interval and sum errors:
\(\sum_{D\le x\ll T/D}(1+T/x)^{1-\kappa}\ll T D^{-\kappa}\).
The conditional expected test is a bounded function of the
\(x\)-state. Mix it using the nonincreasing row length as a
weight: on a dyadic \(x\)-block of scale \(U\) partial summation
bounds discrepancy by \(O_M((T/U)U^{1-\kappa})\).
These errors again sum to \(O_M(TD^{-\kappa})\).
Absorb the domain error by decreasing \(\kappa\) if necessary.
For independent uniform states, conditional on \(x,y\bmod 2M\)
the sum of the two checksums modulo \(M\) is uniform, so \(2a\)
times that sum has image exactly \(2d\mathbb Z\bmod2M\) and point
probability \(d/M\). This gives (3).

For \(\alpha=r/(2M)+\xi/T\) with \(M\) fixed, (4) and partial
summation give the fixed-center major-arc expression
\[
\sum_{\substack{x,y\ge D\\F(x,y)\le T}}e(\alpha F(x,y))
 =
 \left(\sum_{b\text{ odd}\bmod2M}\sigma_M(b)e(rb/(2M))\right)
 \frac1Q\int_{QD^2}^T\log\!\frac{u}{QD^2}\,e(\xi u/T)\,du
 +O_{M,F,\delta}((1+|\xi|)T^{1-\delta\kappa}).
\]
The integral is interpreted as zero if its lower endpoint exceeds \(T\).
It comes from integrating \(W/(Qx)-D\) over \(D\le x\le W/(QD)\)
and differentiating the resulting area; the lattice-area
discrepancy is \(O_F(T/D)\). In particular the central domain at
\(W=T\) has leading presentation mass
\((1-2\delta)T\log T/Q\).

**Only finitely many rational centers have a nonzero arithmetic
factor for any one carrier.** Under independent uniform states,
the average of the checksum phase at \(r/(2M)\) is zero unless
\(M\mid ra\). Equivalently the reduced denominator of the center
must divide \(2a\). This follows by summing the additive
character \(e(ra c/M)\) over the checksum residue \(c\).
A vanishing factor means cancellation at that fixed rational center
and on its bounded \(1/T\)-scale neighborhood. It is not a
uniform bound over growing denominators.

**Single-carrier example with a completely explicit nonvanishing local law.**
The five-vertex class-one mask-zero row has \(a=11,Q=100,s=0,t=89,
g=0,f=11,n=11\). Modulo 11 it reads \(x(y+1)\).
If \(11\nmid M\) the fraction (3) is exactly \(1/M\); otherwise
it is \(21/(11M)\) at \(b\equiv0\pmod {11}\), and \(10/(11M)\)
at the other residues. This follows from the fiber sizes 21
and 10 of multiplication over \(\mathbb F_{11}\) and the fact
that the form is always odd. In particular no fixed congruence
blocks this one carrier, without summing over the larger bank. Its factor at a nonzero
multiple of \(1/11\) is \(1/11\). The corresponding centers shifted
by \(1/2\) have factor \(-1/11\) because the values are odd; the
factors at \(0,1/2\) are \(1,-1\). Other reduced denominators
have zero factor.

The major-arc constants here are controlled only for fixed denominators,
and (2) is a supremum saving, still much too large if inserted
directly in the Fourier integral for one prescribed integer.
Even the explicit everywhere-positive fixed-modulus law plus (2)
does not imply that a target has a witness. Quantitative growing-modulus
control, integrated cancellation, or a new constructive repair is still
needed.

## Reproduction

`python verify.py` checks the published table's coefficient and
parity hypotheses, the interval cut on exact integer examples,
the exact digit factorization, and the displayed single-carrier
rational law. The analytic estimates themselves are proved above.
