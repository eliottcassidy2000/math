# Wrapped terraces: a new sparsity theorem and an unbounded channel ladder

Use the [single-wrap readout](https://github.com/eliottcassidy2000/math/blob/8051b02275faab24bffba7fbed5cc8d3afc29bac/05-knowledge/results/chanpak_target_first_terrace_density_20261114/SINGLE_WRAP.md):
put (u\to B\to v) and reverse their outer arc to (v\to u).
For nonempty (B) the total is (H(B)+L_B(1)). Here (L_B(1))
sums the product of the Hamiltonian path counts on two color
classes over all two-colorings of the middle, including empty
classes. These statements concern one specific wrapped family, not
arbitrary switched tournaments.

## A single wrap around arbitrarily many terrace blocks is still sparse

**Proved:** let (B) be an arbitrary ordinal sum of minimum-triangle
strong tournaments and singletons. The distinct (H(W(B))\le X)
number (O(X^{93/100})). Neither the order of the blocks nor their
number is bounded.

First consider a single normalized terrace with (n\ge4).
Its (n-4) unrestricted anchor bits are contiguous. A three-site
window of type `1?1`, with first and third sites both anchors, is
a directed triangle. For `101` there is one backward arc from
the third to the first and forward arcs through the middle. For
`111` there are two consecutive backward arcs and one forward arc.
Greedily take the first such window and restart just after it;
the resulting triangles are pairwise disjoint. Write (j) for
the number selected.

From the odd-cycle packing expansion
[
 L_B(1)=2^n\sum_{\mathcal C}\prod_{C\in\mathcal C}2^{2-|C|}.
]
Choosing subsets of the (j) disjoint triangles gives the lower
bound (2^n(3/2)^j) on (L_B(1)). No selection of local bits is
being mistaken for exact global count.

For (s=93/100) set (t=(2/3)^s). Read arbitrary bits while
remembering the last at most two unmatched bits. A successful
`1?1` window sends the scanner to the empty state with weight (t);
other steps have weight one and retain the last at most two bits.
Let (S_t) be the resulting nonnegative weighted adjacency matrix
on states `empty,0,1,00,01,10,11`, using a row per current state.
Then ((S_t^m\mathbf1)_{\mathrm{empty}}) is the sum of (t^j)
over all words of length (m). The following positive test vector
gives a strictly subcritical bound:
[
 b=(1000,1000,901,1000,901,886,827)^t,
 \qquad S_tb\le(1904/1000)b,
 \qquad 1904/1000<2^{93/100}.
]
The integer replay checks these with (t<686/1000), which
follows by raising the positive quantities to power 100.
Thus the sum of (L_B(1)^{-s}) over all terrace words is finite:
for each order use (2^{-ns}(S_t^{n-4}\mathbf1)_{\rm empty})
and then sum the geometric series. Include the two finite small
cases (singleton and order-three cycle) separately.

Now (H(B)) and (L_B(1)) both multiply over the blocks in any
ordinal sum, so the wrapped scalar depends only on their multiset.
Sum over all multisets of the constituent terrace classes and
singletons. The Euler product for the inverse (s)-powers of
their (L(1)) values converges, since every value is at least two
and the single-block sum was finite. The actual (H(W(B))\) is
at least that product of (L) values, giving the bound on
distinct values. This argument does not claim that the wrapped
tournament itself is determined by the multiset.

Therefore arbitrary length and ordinary ordinal gluing do not
make a single outer wrap sufficient to cover all prime targets.
The obstruction is spatial: nonoverlapping local triangles force
the auxiliary readout to pay more than the word's raw (2^n)
description cost, in the precise weighted-series sense above.

## Iterate the wrap: all color channels genuinely matter to linear transfer

For a nonnegative integer (m) let
[
 F_m(B)=\sum_{c:V(B)\to[m]}\prod_{i=1}^m H(B[c^{-1}(i)]),
]
with (H(\varnothing)=1). In particular (F_1=H), (F_2=L(1));
(F_0(B)) is one exactly when (B) is empty. These functions are
multiplicative under ordinal sum. They obey the exact recurrence
[
 F_m(W(B))=
 m^2F_m(B)+m\bigl(F_{m+1}(B)-F_{m-1}(B)\bigr),\qquad m\ge1.
]

To prove it, color the two new distinguished vertices. If they are
different colors their appearances do not affect any color-class
count, since each is a single terminal vertex. This gives (m(m-1))
copies of (F_m(B)). If they are equal in a given color (i), the
count on that color uses the wrap formula on its part (S) of the
middle. The original term gives one more (F_m(B)) per (i).
Its correction for (S\ne\varnothing) splits color (i) into two
distinguished color classes and gives (F_{m+1}(B)-F_{m-1}(B)):
the subtraction removes the case when both distinguished classes
are empty. Multiply by the (m) choices of (i). This proof also
handles an empty middle: then the only change is that the two new
vertices form a transitive pair after the reversal.

In particular for an (n)-vertex transitive middle,
(F_m(TT_n)=m^n). Start at the row (F_1), and iterate (k)
times. The highest channel reached is (F_{k+1}), with coefficient
(1\cdot2\cdots k=k!), since reaching it requires the +1 step
at every move. Thus the functions (n\mapsto H(W^k(TT_n))),
(k\ge0), are linearly independent over the rationals: the
largest base (k+1) cannot cancel against earlier smaller bases.
There is no fixed-dimensional **linear** interface of arbitrary
functions of the middle that supports all those outputs by linear
readouts, hence no fixed matrix state update of that type.
This does not rule out nonlinear or variable-size encodings.
It also does not imply folded, repeated wraps cover the odd
spectrum. It gives an explicit mechanism for why they can escape
the single-wrap proof: the latter uses two coordinates, while
repetition calls for progressively larger color partitions.
