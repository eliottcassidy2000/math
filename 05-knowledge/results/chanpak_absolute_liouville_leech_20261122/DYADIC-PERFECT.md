# Why the original unit chart has no fixed dyadic formula at perfect centers

**PROVED, for a specified strategy; not an obstruction to individual
landings.** The original six-vertex two-sided core in
[Gap star ports](../chanpak_gap_star_ports_20261122/README.md) cannot
hit \(E_p+j\) at unbounded exponents with any fixed bounded-offset pair
of dyadic polynomial addresses, for any fixed odd \(j\). Here a dyadic
polynomial address means \(x=X(2^{p-1})\), \(y=Y(2^{p-1})\) for fixed
polynomials \(X,Y\) over \(\mathbb Z[1/2]\), taking positive integral
values at all exponents used. This includes both flanks of order
\(\sqrt{E_p}\) and one fixed flank with one of order \(E_p\).
In fact the one-fixed-flank case is impossible for **every** fixed core
with the same form of linear binary read and odd \(h\). Other cores
with two growing flanks, and addresses rebuilt from the entire target
at each exponent, are not ruled out.

## The binary read at the leading scale

Let \(q=2^m\). Let \(Y(q)\) be an integer-coefficient polynomial with
degree \(d\ge 1\), leading coefficient \(u>0\), and constant coefficient
\(w\in\mathbb Z\). Its leading binary-read coefficient is
\[
    p(Y(q))=\kappa(w)u q^d+O(q^{d-1}),\qquad
    \kappa(w)=
    \begin{cases}
     1-\operatorname{rad}_2(w), & w\ge 0,\\
     \operatorname{rad}_2(-w-1), & w<0.
    \end{cases}
\]
The function \(\operatorname{rad}_2\) is the base-two radical inverse:
\(\operatorname{rad}_2(\sum_i b_i2^i)=\sum_i b_i2^{-i-1}\).
In particular \(\kappa(w)\) is dyadic in \([0,1]\), and equals
1 precisely at \(w=0\) and 0 precisely at \(w=-1\).

To prove the formula, write \(Y=qV+r\) with \(0\le r<q\). Then
\(V=u q^{d-1}+O(q^{d-2})\) for \(d\ge2\), and \(V\) is eventually a
fixed integer for \(d=1\). The unrolled binary recurrence gives
\[
 p(qV+r)=p(V)+(q-1-\operatorname{rev}_m(r))V+p(r)
             +[r<q/2][V\ \text{odd}].
\]
For \(w\ge0\), \(r=w\), so \(p(r)=O(1)\) and the reversed term is
\(q\operatorname{rad}_2(w)\).
For \(w<0\), \(r=(q-1)-(-w-1)\), so \(q-1-\operatorname{rev}_m(r)\)
is \(q\operatorname{rad}_2(-w-1)\). When \(d\ge2\) the term
\(p(r)\le q\) is already lower order. For \(d=1\) one must retain
it: if \(w<0\), write \(r=2^\ell(2^{m-\ell}-1)+v\) for fixed
\(\ell\) so that \(v=2^\ell+w\) is in \([0,2^\ell)\). Applying the
same recurrence at width \(\ell\) shows
\(p(r)=q\operatorname{rad}_2(-w-1)+O(1)\), since the high word is all
ones. Then \(V=u-1\), not \(u\); these two terms sum to the coefficient
claimed. For \(d=1,w\ge0\) the fixed high word is simply \(u\).
Dyadic polynomials reduce to this proof after rescaling \(q=2^{p-1-s}\)
for fixed \(s\) that clears the nonconstant denominators. If the
resulting polynomial has integer values at unbounded exponents, its
remaining constant coefficient is also an integer.

### An exact polynomial test, not just an asymptotic one

The same proof computes the *whole* polynomial in finitely many steps.
Write \(Y(q)=qU(q)+w\). When \(w\ge0\), for sufficiently large \(q\) the
unrolled recurrence becomes
\[
p(Y)=p(U)+\big(q\kappa(w)-1\big)U+p(w)+[U\ \text{odd}].
\]
When \(w<0\), choose \(B=2^\ell>-w\), put \(v=B+w\) and
\(R=\operatorname{rev}_\ell(-w-1)\). The high word is \(V=U-1\) and the
low word is \(q+w>q/2\). The two nested unrolls give
\[
p(Y)=p(U-1)+q\kappa(w)U+p(v)-R+[v<B/2].
\]
The parity indicator in the first case is eventually constant, namely
the parity of \(U(0)\). In both cases the recurse argument has smaller
degree and is eventually nonnegative; the base case is one ordinary
integer evaluation. The margins used here, such as \(q>2|w|\) and
\(q\ge2B\), give an effective width threshold. Thus a fixed template
can be tested exactly by comparing coefficients in the two count
polynomials; agreement at infinitely many exponents is equivalent to
identity, and the identity gives all exponents above that threshold.

## Two growing flanks: the two-slot mass must be a power of two

In any fixed marked core, the mixed term is \(Qxy\), with all other
terms \(O(x+y+1)\) because \(0\le p(v)\le v\). If both addresses
grow, \(E_p=2(2^{p-1})^2-2^{p-1}\) forces the sum of their degrees
to be 2. Its leading coefficient 2 is the product of \(Q\) and
two positive dyadic leading coefficients. Thus **\(Q\) must be
a power of two**. The unit core in question has \(Q=320=5\cdot64\).
This also excludes the unrescaled two-sided cube-root chart from an
all-exponent exact dyadic template: changing low bits cannot repair
the leading mixed coefficient.

## One growing flank: the only eligible slopes are the two extremes

With one fixed address, the count in the prescribed parity chamber
has the form \(K+cY+h p(Y)\), with \(h\) odd. The chamber is eventually
fixed, because the parity of a dyadic polynomial at powers of two is
eventually the parity of its constant term. Since the count is
asymptotic to a positive multiple of \(Y\), the growing address
has degree two. After clearing denominators as above, the
leading equality requires \(u(c+h\kappa)=2^{2s+1}\).
Thus \(c+h\kappa\) is a power of two in the rationals; since
\(c\ge1\), it is an **integer** power of two. Since \(h\) is odd and
\(\kappa\) dyadic, integrality forces \(\kappa\in\{0,1\}\).
Accordingly, any fixed-core rule of this kind must use either
\[
       w=-1,\quad c\ \text{a power of two};\qquad
       w=0,\quad c+h\ \text{a power of two}.
\]
These are leading-term eligibility tests, not sufficient conditions.
A core that passes can be tested next on its linear and constant
coefficients; the parity of its edge must be retained separately.

For the two original flip cores \(h=29\). If the left flank \(x\)
is fixed, \(c=129+320x\) is odd and greater than 1; \(c+h=158+320x\)
has 2-adic valuation exactly one and is greater than 2. If the right
flank \(y\) is fixed, \(c=162+320y\) has valuation one and is greater
than 2; \(c+h=191+320y\) is odd and greater than 1. Thus no chamber
of this chart passes, regardless of which endpoint is designated.
Both bounded addresses give bounded counts and are impossible as well.

### The linear coefficient rules out both one-flank cases in any core

Take \(q=2^{p-1-s}\), and write the only remaining
address as \(Y=u q^2+v q+w\), with \(u,v\) integers.
The perfect-number polynomial is \(2^{2s+1}q^2-2^sq\).
Write \(2^r\) for the eligible leading slope of the chart, so
\(u=2^{2s+1-r}\). In the case \(w=-1\) we have \(p(Y)=p(uq+v-1)\).
The coefficient of \(q\) must be \(-2^s\). It forces \(v<0\), because
the binary-read leading coefficient is nonnegative. Write
\(-v=2^t J\)
with \(J>0\) odd, let \(k=r+t-s\), and use
\(\operatorname{rad}_2(2^tJ)=2^{-t}\operatorname{rad}_2(J)\).
The equation becomes
\[
             2^k J-h2^{1-k}\operatorname{rad}_2(J)=1. \tag{A}
\]

In the case \(w=0\),
\(p(Y)=(q-1)(uq+v)+p(uq+v)+[(uq+v)\ \text{odd}]\).
Now \(v=0\) would give zero linear coefficient. If \(v=2^tJ>0\)
the same calculation gives (A) with right-hand side \(-1\).
If \(v=-2^tJ<0\), the identity
\(1-\operatorname{rad}_2(2^tJ-1)
=2^{-t}(3/2-\operatorname{rad}_2(J))\) gives
\[
             2^k J+h2^{1-k}(3/2-\operatorname{rad}_2(J))=1.
\]
This last left-hand side exceeds 1 by AM–GM: its positive terms
have product \(2hJ(3/2-\operatorname{rad}_2(J))>1\).

For either sign in (A), let \(L\) be the binary length of \(J\) and
\(R\) its \(L\)-bit reversal. Both \(J\) and \(R\) are odd and in
\([2^{L-1},2^L)\). Then \(2^kJ\mp1=h2^{1-k-L}R\).
For \(k>0\), the left side is odd while the right is not even
integral. For \(k=0\) the left side is even (or zero), but the
right side has valuation \(1-L\le0\). For \(k=-\ell<0\),
the valuation instead gives \(L=2\ell+1\); multiplying out gives
\(J\mp2^\ell=hR\).
Since \(\ell\ge1\), the left side is strictly smaller than
\(3\cdot2^{2\ell}\), while \(R\ge2^{2\ell}\). The only
possible positive odd \(h\) is 1. But a binary word and its reversal
cannot differ by \(2^\ell\) in absolute value: if they differ, the
least differing position is strictly below the middle \(\ell\) by
reflection, and gives the exact valuation of their difference.
This contradiction rules out both remaining one-flank cases.
It is independent of the intercept and works with any fixed offset.

The useful distinction is between a repeatable fixed polynomial recipe
and a target-dependent inverse. A finite exact hit does not imply
the former, and these obstructions do not negate the latter.

Replay the exact coefficient recursion on fixed signed polynomials via
`python3 check_dyadic_read.py`. The impossibility itself follows from
the valuation argument and does not rely on a computational cutoff.
