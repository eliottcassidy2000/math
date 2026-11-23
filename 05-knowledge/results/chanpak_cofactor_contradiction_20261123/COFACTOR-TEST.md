# Cofactor-tolerant, additive absence witnesses

**PROVED scoped; witness supply OPEN.** Let \(d_i\) denote the
ordinary Hamilton count of the \(i\)-th SCC of a tournament,
so \(x=H(T)=\prod_i d_i\). Counts for singleton SCCs are one.
For positive odd \(m\), \(\Omega(m)\) counts prime factors with
multiplicity, and \(\Omega(1)=0\).

## Prime lane

Fix an odd prime \(p\), put \(r=\nu_p(x)\), \(m=x/p^r\), and
\(K=\#\{i:p\mid d_i\}\). If no tournament has count \(p\), then
for each of those \(K\) components its local cost
\[
  (\nu_p(d_i)-1)+\Omega(d_i/p^{\nu_p(d_i)})
\]
is at least one. Summing gives \(2K\le r+\Omega(m)\).
Equivalently the additive score \(2K-r-\Omega(m)\) is nonpositive
on every tournament under absence. A positive score implies a
literal component counted by \(p\); no factoring operation on
graphs is presumed.

In multiplicative terms the sharper envelope is
\[
  p^K\le \sqrt{x}\sqrt{p^{\Omega(m)}/m}.
\]
For \(m\le x^\varepsilon\) one may use
\(p^K\le x^{1/2+(\log p)/(2\log3)\varepsilon}\); if all
factors of \(m\) are at least \(p\), no error exponent is
needed. Thus a strict excess exponent for fixed \(p\) is a
sufficient certificate even with a controlled cofactor.
For example positive integers with \(m=1\) recover the pure
power-ray grade test. This improves the permitted witness
domain, not the known coverage.

## Single-seven lane with many simultaneous 7-carriers

Fix \(p\ne7\), set \(r=\nu_p(x)\), \(e=\nu_7(x)\), and
\(m=x/(p^r7^e)\). Retain \(K\) and mark
\(Z=\#\{i:\nu_7(d_i)=1,\ \nu_p(d_i)\ge1\}\).
If no tournament has count \(7p\), then
\[
   Z\le (r-K)+\Omega(m).
\]
Indeed each marked component that does not include a new
cofactor prime must carry at least two \(p\)'s. Extra copies
of \(p\) in the \(p\)-bearing components contribute exactly
\(r-K\), and the other cofactors across all components cost
at most \(\Omega(m)\). This condition avoids assuming that
all 7-carriers have multiplicity one or even carry \(p\).
A violation forces some actual component counted by \(7p\).
Both tests are additive on ordinal sums, so they can also be
used with nonnegative averaging over genuinely constructed
hosts, never over a set of unverified scalar labels.

If a procedure supplies \(Z\ge\delta r\),
\(K\ge(1-\varepsilon)r\) and \(\Omega(m)\le\varepsilon r\)
at a positive \(r\) with \(2\varepsilon<\delta\), the
hypothetical missing \(7p\) is impossible. The examples
\(p=3\), ordinary strong heights 63 and 105, expose the
two separate charges: 63 has one surplus 3 and no extra
cofactor, while 105 has no surplus 3 and one extra prime
5. In each case \(Z=1\) and equality holds. Height 63
has a single strong core since 7 and 21 are absent; for 105
use the explicit strong \(A_9\) in the [signed endpoint family](../chanpak_guarded_ternary_packing_20261123/PRODUCT.md). These are controls at a known absent \(7p\),
not evidence that such an excess exists for each allowed \(p\).

Together these tests give a sufficient route to rejecting
the first halt: for every still unexcluded prime \(p\ne7\)
one needs a witness of positive prime score, and for every
still unexcluded \(p\ne3,7\) one of positive seven score.
One can concentrate on whichever lane the assumed first
missing integer lies in. The supply of these witnesses
remains the substantive open step.

## Why near power hits on the original affine ray do not suffice

There is even an effective archimedean obstruction to using a finite
bank of fixed ordinal multipliers on the **original** dyadic rail
as a near-power fitter. Fix positive integers \(A,b\) with \(b\) odd
and fix an odd prime \(p\), and an integer \(c\).
For \(Y=bp^r\) one has, for some effective \(B,C>0\)
depending on \(A,b,p,c\) and all sufficiently large \(r\),
\[
 \inf_{s\ge1}|Y-(A2^s+c)|\ \ge\ C Y (1+r)^{-B}.
\]
If \(A2^s/Y\) lies outside \([1/2,2]\), the bound is immediate.
Otherwise \(s=O(r)\), and the linear form
\[
 r\log p-s\log2+\log(b/A)
\]
does not vanish: its vanishing would equate an odd integer to an
even integer. The fixed-logarithm specialization of
[Matveev's lower bound](https://doi.org/10.1070/IM2000v064n06ABEH000314)
bounds its magnitude below by a constant times \((1+r)^{-B}\).
On the compact ratio interval the logarithm is Lipschitz, giving
the corresponding absolute distance from \(A2^s\). The fixed
intercept \(c\) is eventually smaller than half the bound.
This includes \(b=1,7\) and any finite bank of fixed rescalings
of the given parent and its three children. In particular it
excludes a square-root-scale placement error for these rays at
all sufficiently late fixed-prime power targets. It makes no
claim about varying multipliers, arbitrary fillers, or the full
set of tournament witnesses.
