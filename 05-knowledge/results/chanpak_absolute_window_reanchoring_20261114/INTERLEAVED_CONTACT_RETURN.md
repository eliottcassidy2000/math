# First carries, interleaved exact switches, and a transitive-core return

This is a physical three-switch return on infinite arithmetic lines of source
parameters, not a proof that such a line can be chosen for every absolute
target. It uses an apex over a wholly transitive old tournament, rather than
a neutral multi-core bank.

## Exact single-spacer readout

Let a fixed marked \(A\) precede a variable transitive \(m\)-vertex spacer, then
a fixed marked \(B\), with all old cross-arcs forward. A spacer vertex receives
contact bit one when the new apex beats it. Write the bits in transitive order,
most significant first, as integer \(z\), retaining leading zeros for the
fixed spacer length. Let \(p(z)\) be the low bit, \(t_m(z)\) the high bit and
\(r_m(z)\) the \(m\)-bit reversal. Let \(g(z)\) denote the one-apex count with
just these transitive old vertices. Leading zeros do not affect \(g\), as they
are universal sources among the old vertices and apex. Then, up to a fixed
constant independent of \(z\), the full path count is
\[
h_AQ_B z -P_A h_B r_m(z) +h_A h_B g(z)
+h_A T_B p(z)-S_A h_B t_m(z).
\]
To derive it from the ordinal composition formula, the variable block has
\(h=1,f=2^m,P=z,D=r_m(z),S=p(z),b=t_m(z)\), and
\(E=1+z-p(z)-g(z)\). Sum the direct terms of \(K=h+P-S-E\)
in the composed core. In particular the coefficient of \(z\) is
\(h_A(f_B-h_B-D_B)=h_AQ_B\).

## Why separated multibit digits interact

Unrolling the recurrence for \(g\) in little-endian bits \(z_i\) gives
\[
g(z)=1+\sum_{i<j}(1-z_i)z_j\left(2^{j-i-1}+\mathbf1_{j=i+1}\right).
\]
The second part of the weight is the parity carry. Thus all mixed
differences of the full readout are governed by the distance matrix with
weights \(w(d)=2^{d-1}+\mathbf1_{d=1}\). In two switchable disjoint ordered
intervals with signed differences \(a_i,b_j\in\{-1,0,1\}\), the mixed
difference divided by \(-h_Ah_B\) is
\[
\Big(\sum_i a_i2^{-i-1}\Big)\Big(\sum_j b_j2^j\Big)+a_{\ell}b_{\ell+1},
\]
where the last term is present only for adjacent supports.
Both sums are nonzero since the dominant nonzero digit outweighs
all smaller ones. The exceptional adjacent carry can cancel the
leading product when one support is a singleton, but **two ordered
patches that each change at least two bits cannot cancel**.
See [the exact boundary and two-site classification](SWITCH_CLASSIFICATION.md).
The three-switch construction below is unaffected by this exception.

This is a precise analogue of a sampled old repo obstruction:
unit-distance reduction can preserve every residue edge while losing the
first arithmetic carry. Here the planar-distance carry is inspiration,
not an imported theorem: the mixed bit carry is explicitly counted.

## An interleaved exception with three exact choices

The distance matrix commutes with reversal, so any symmetric signed
switch is orthogonal to any antisymmetric one on disjoint supports.
For two distinct antisymmetric elementary mirror pairs with positions
\(i<j<j'<i'\), the cross term is twice \(w(j-i)-w(j'-i)\) up to sign.
The only equal distinct positive distances are \(w(1)=w(2)=2\).
Thus in the four central positions two antisymmetric moves can commute,
and an outer symmetric switch can commute with both. The parity carry
is exactly what makes this zero; with the bare exponential weight the
two anti switches would interfere.

Use the six-vertex words \(z=24+33a-14b-4c\) for independent
\(a,b,c\in\{0,1\}\). The outside pair is flipped from 00 to 11; the middle
two mirror pairs are swapped independently. Direct evaluation, or the
zero mixed coefficients and four corner evaluations, gives
\(g(z)=23-14b-4c\).
Prepend a transitive left wing encoding \(x\), and append the right wing
encoding \(y\) in the reversed complemented convention. Then the
eight counts are the fixed baseline
\[
H_0=64xy+57x+24y+g(x)+g(y)+21+p(x),
\]
plus
\[
a[33(y-x)+p(y)-p(x)]-(14b+4c)(x+y+1).
\]
In particular change all three switches simultaneously. The shift is
\(15y-51x-18+p(y)-p(x)\), and is exactly
- \(-2\) for \(x=10t,y=34t+1\);
- \(+2\) for \(x=10t+9,y=34t+32\).

For \(t\ge1\), both wings are nonempty, and taking their ordinary binary
representations makes the full tournament strong: the apex reaches the
first old vertex and the last old vertex reaches the apex. The \(x=0,
y=1\) calculation also gives a literal 47/45 pair; use an empty left wing.
For positive-height lines, the \(H_0\) term is quadratic in \(t\) plus
\(g\) terms. These are not a prescribed-target covering family.

## Replay

Run "python interleaved_controls.py". The first part independently builds
literal outgoing rows and counts small complete tournaments by subset DP,
including both arithmetic lines. The larger symbolic controls check the
four corners and the exact contact formula against the standalone
single-apex recurrence throughout the two arithmetic lines.

## A further decoder restriction

In a fixed spacer, moving a marked one down one position across a vacant
zero strictly lowers \(g\). With little-endian positions \(i,i+1\), write
\(q=\lfloor z/2^{i+2}\rfloor\) for the part above them. The local zero-payment
drop is \(q+2-p(q)\ge2\). If position \(i-1\) is a zero, its floor payment
decreases by one and its parity addend increases by one, exactly cancelling.
Payments at all still lower zero positions decrease, and higher payments
are unchanged. Thus moving a one from \(j\) to vacant \(i<j\) while holding
all other bits fixed decreases \(g\) by at least \(2(j-i)\), by decomposing
the change into adjacent swaps.

For the whole marked tournament one down swap decreases \(z\), increases
\(r_m(z)\), and therefore neither of the two main linear terms can undo the
drop. At an extreme bit, one endpoint correction can undo at most
\(h_Ah_B\). For \(m\ge3\) an adjacent swap meets at most one extreme.
So the full readout strictly decreases on every such move, and the fiber
on any fixed popcount is an antichain in the corresponding move order.
This provides an exact branch-and-bound order for target decoding, but
does not guarantee that its intermediate positive odd values all occur.
