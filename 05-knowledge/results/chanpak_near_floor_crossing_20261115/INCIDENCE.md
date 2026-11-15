# Carry-separated physical carriers: \(\Omega(T^{2/3})\) exact counts

Use the order-five class-one mask-zero physical carrier
\[
F(x,y)=100xy+89x+11(p(x)+p(y))+11+11(x\bmod2),\qquad x,y\ge1.
\]
Here \(p(0)=0\), \(p(2w+1)=p(w)\), and \(p(2w)=p(w)+w+(w\bmod2)\).
The [sharp-cut packet](../chanpak_true_cut_fourier_20261115/README.md)
has \(0\le p(z)\le z\), and its binary split states that for \(V=2^r\),
\(0\le b<V\),
\[
 p(Vu+b)=p(u)+p(b)+(V-1-\operatorname{rev}_r(b))u+
 {\bf1}_{b<V/2}(u\bmod2).
\]

**Theorem.** For \(V=2^r\ge8\), the \(V^4/16\) points
\[
 V^2\le u<2V^2,\quad V/2\le s<V,\quad 0\le b<V/8,
 \qquad x=Vu+b,\quad y=Vu+s-b
\]
give distinct absolute Hamiltonian path counts, all at most \(500V^6\).
Thus this one physical carrier has \(\Omega(T^{2/3})\) distinct
counts up to \(T\).

**Proof.** Set \(c=s-b\), so that \(0\le c<V\), and \(S=x+y=2Vu+s\).
For ordinary binary addition \(b+c=s\), let \(k_i\) be the carry out
of bit \(i\), with \(k_{-1}=k_{r-1}=0\). Reversing the per-bit
identities \(b_i+c_i=s_i+2k_i-k_{i-1}\) gives
\[
 \operatorname{rev}_r(b)+\operatorname{rev}_r(c)
    =\operatorname{rev}_r(s)+3K,\qquad
 K=\sum_{i=0}^{r-2}k_i2^{r-2-i}.
\]
Thus on a fixed line \((u,s)\), the large part of the checksum that
varies with \(b\) is precisely \(-3uK\). The low part \(p(b)+p(c)\)
is between \(0\) and \(s<V\); the parity boundary in the split lies
between \(0\) and \(2\).

For distinct \(K\)'s, the large checksum contribution to \(F\)
differs by at least \(33u\ge33V^2\). All the other variation on the
line has absolute difference at most
\[
 25V^2+(89/8+11)V+33<33V^2.
\]
Here the varying product part is \(bc\in[0,V^2/4]\), the varying
linear part is \(89b\), and 33 bounds both parity differences.
Thus different carries cannot collide. On the same carry, if \(b'>b\)
then
\[
 b'(s-b')-b(s-b)=(b'-b)(s-b'-b)\ge V/4,
\]
and \(89x\) increases. The product contribution increases by at
least \(25V\), which exceeds the adverse change \(11V+33\) in low
checksum and parities. This proves injection within each line.

For different lines use \(S\), which determines \((u,s)\). Relative
to the center \(L(S)=25S^2+(89/2)S\) the whole line lies in the band
\[
 [L(S)-25V^2-(89/2)V,\ L(S)+11S+(89/2)V+22].
\]
Indeed \(100xy=25S^2-25(x-y)^2\), \(|x-y|<V\), and
\(p(x)+p(y)\le S\). Ensure separation of successive integer values
\(S,S+1\) in the range \(S\ge2V^3\):
the center increase \(50S+25+89/2\) exceeds
\(11S+25V^2+89V+22\) for \(V\ge8\). Then the bands are ordered and
disjoint for any two integer values, not just occupied ones.
Finally \(S<4V^3\), so the displayed upper band is less than
\(400V^6+222V^3+(89/2)V+22<500V^6\). Choosing the largest power
of two \(V\le(T/500)^{1/6}\) yields the claim.

This is separation for a chosen subset, not an approximation or
representation of prescribed targets. The factor 3 is the difference
between doubling a carry and halving its reversed place value.

## Stronger support by averaging collisions across successive lines

**Theorem.** The same single carrier has \(\Omega(T^{3/4})\) distinct
counts up to \(T\). Specifically, for every power of two \(V\ge8\)
the points
\[
 V\le u<2V,\quad V/2\le s<V,\quad 0\le b<V/8,
 \qquad (x,y)=(Vu+b,Vu+s-b)
\]
give at least \(7V^3/128\) distinct counts at most \(500V^4\).
The theorem concerns counts, not just presentations.

The bands indexed by \(S=2Vu+s\) above stay disjoint. In fact their
separation inequality only needs \(S\ge2V^2\); at that bound
\(11S+25V^2+89V+22<50S+25+89/2\). Since \(S<4V^2\), the upper
band is below \(500V^4\). Keep all \(b\)'s, not just one carry
class. The identity \(p(b)+p(c)\) and product-difference argument
above show that fixed \((u,s)\) and fixed \(K\) give no collisions.

For \(b\ne b'\) with distinct carries \(K\ne K'\), fix \(s\) and
subtract the two counts on the line indexed by \(u\). The common
\(100(V^2u^2+Vus)\), \(89Vu\), \(22p(u)\) and shared reversed
checksum term cancel. The exact difference is
\[
 -33(K-K')u + A(s,b,b')+B(s,b,b')(u\bmod 2)
\]
for integers \(A,B\) independent of \(u\). In particular, for each
ordered triple \((s,b,b')\) there is at most one collision on even
lines and at most one on odd lines. This uses the full exact
integer checksum, not an approximation or genericity assumption.

There are \(P=V^3/16\) presentations. The total ordered
off-diagonal collision count \(C\) is at most
\(2\,(V/2)(V/8)(V/8-1)\le V^3/64\). Each value with \(k\) physical
presentations contributes \(k-1\le k(k-1)/2\) to the loss of support,
so the number of distinct values is at least
\(P-C/2\ge7V^3/128\). Taking the largest admissible power of two
\(V\le(T/500)^{1/4}\) proves the result at all sufficiently large
\(T\). The simpler pointwise collision-free proof used \(u\) of
order \(V^2\); the improvement comes from charging each opposing
carry pair to at most two absolute baselines, which allows \(u\)
of order only \(V\).

The [replay](verify_incidence.py) verifies all these incidence identities
and rational roots independently of sampled spectral sums. It is
not needed for the exponent proof. Density one for the full
carrier is still open.
