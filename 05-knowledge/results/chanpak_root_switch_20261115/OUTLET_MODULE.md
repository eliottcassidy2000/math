# Outlet module gluing and the cost of changing the rooted remainder

This extends the one-outlet identity to an arbitrary nonempty block, and gives a
closed transfer when the next root is either the current root or in the new block.
The scope of the low-switch obstruction is a fixed finite alphabet of such blocks;
it is not an impossibility for arbitrary tournament constructions.

## Four coordinates

For a rooted tournament `(W,u)` let `h=H(W)`, `g=H(W-u)`, and let `s,d`
be the total and odd-minus-even weighted counts of nonempty directed paths starting
at `u`, with residual weight `H(W-P)`; take `H(empty)=1`. Use
`z=s-g` (the weight of nontrivial such paths) and `b=h-s-d+2g`.
The signed endpoint identity gives `d=T_u(W)`, the ordinary full count ending at
`u`; see the signed specialization in the companion
[endpoint reciprocity note](../chanpak_ternary_emptying_20261115/ROOTED_SLOT_RECIPROCITY.md).
In particular `0<=d<=h`.

For vertex-disjoint nonempty R, sew `W -> R` except that every vertex of R beats
`u`. Write `a=H(R)`, `q=P_R(2)`, where
`P_R(2)=sum_{A subseteq V(R)} H(R[A])H(R-A)`. The result has path count

```
h' = a b + q z
   = a (h + s - d) + (q - 2a) z .
```

For `R=TT_l` this is `b+2^l z`, recovering `H+2E_u` at `l=1`.
For a general block the two coefficients are actual statistics of the same R,
not independently tunable integers.

**Proof.** Use the odd directed cycle-packing formula, with factor two for each
nontrivial cycle. A crossing cycle must enter W through u, hence must run on a
single nontrivial W-path starting there, and a single nonempty R-path. For any
nonempty R, optional R-path weight, including the empty choice, splits equally
by parity: assigning its odd packing blocks to the path or to the complement
gives zero signed total. When the W-path has even order, the R-path must have odd
order and contributes `q/2`. For an odd-order W-path of order at least 3 it
contributes `q/2-a` (the empty set cannot close a crossing cycle). The contribution
of no crossing cycle is `a h`. This proves the displayed expression. The omission
of the singleton W-path is essential: its final vertex is u and cannot go back
into R.

## Retain the root or choose a new one

Retaining u gives
```
s' = a g + q z,     d' = g' = a g,     b' = a b,     z' = q z.
```
Indeed a path starting at u either is the singleton or uses a nontrivial W-prefix
and an optional R-tail. The tail has zero signed sum, leaving only the singleton
for `d'`. The deletion `W-u -> R` is ordinal.

Now choose `r in R`. Let `s_i,d_i` for `i=1,2` be the rooted total and
odd-minus-even in R with residual weight `P_{R-P}(i)`, and let `f` be
the ordinary full count starting at r. Then
```
s' = (b+g) s_1 + 2 z s_2 + f (d-g)
d' = (b-g) d_1 + z d_2
g' = b H(R-r) + z P_(R-r)(2) + [|R|=1] (d-g).
```

To see the start formula, partition by whether the path enters W. A path staying
inside R has residual mass `b H(R-P)+z P_(R-P)(2)`, except when `P` exhausts
R, in which case the missing term is `d-g = h-(b+z)`. A path entering W must
enter through u. If its W segment is just u, the extra contribution is
`g s_1`; otherwise it may return to R to consume a second segment, so the
extra contribution is `z s_2`. This arbitrary return segment explains the
second copy of `z s_2` in the result. For the signed formula W-singleton
contribution has the opposite sign, and an optional tail of nonempty R-remnant
has zero signed sum. When the R-prefix exhausts R, entering W with a nontrivial
path cancels exactly the exceptional exhaustion term from staying inside R.
The deletion formula follows by applying the count formula to `R-r` with an
explicit empty convention.

For `R=TT_l`, index r from zero and put `e=l-r-1`. Direct subset selection
gives `s_1=2^e, s_2=2^r3^e, f=[r=0], d_1=[e=0], d_2=2^r`.
For example, if the new root is the last vertex of its block, then `b'=g'`
and `z'=g+2^(l-1) z`. Once the current state satisfies `b=g`, a next last-root
block of length at least two collapses to `b'=g'=z'`. Therefore repeated such
steps are genuinely multiplicative, not freely switchable binary digits:
on the rail `b=g=z`, each step multiplies the rail by `1+2^(l-1)`.

## A physical unit repair carried through any retained width

There are two strong order-five sources with identical `(s,d,g)=(15,3,3)`
and unequal path count. Outgoing rows indexed from zero, followed by the root:

* `(12,21,24,18,1)`, root 0, `h=13`, hence `(b,z)=(1,12)).
* `(14,20,24,18,1)`, root 4, `h=15`, hence `(b,z)=(3,12)`.

For **any** allowed appended module of profile `(a,q)), these give the
two counts `12q+a` and `12q+3a), preserving identical `s',d',g'`
when the old root is retained. In particular, adjoining `TT_l` gives
two adjacent odd counts `12*2^l+1` and `12*2^l+3` for every `l>=1`.
These are actual strong witnesses: the old strong W can reach every newly
added vertex through an old vertex other than u, and every new vertex beats u.
A target in this family is selected without a search: subtract its residue
1 or 3 modulo 12, divide by 12, and test for a power of two at least 2.

The repair is genuinely continuation-specific. Upon moving the root into r
in a new block R, the two states differ by `2*s_1` in s and by
`2*H(R-r)` in g; these are never zero. A sink r may preserve the z
difference at that move, but this does not preserve the full port state.
A route beyond the fixed-root dyadic pair needs to choose new physical
sources so their current weighted difference lands in a usable return
profile, not simply inherit the 2-offset after changing roots.

## A root-switch lower bound

Fix finitely many nonempty blocks and finitely many rooted seeds. A switch means
moving the root into a new block; retaining the root just glues on a block.
For any `k=o(log X)`, the number of distinct counts at most X attainable by
words with at most k switches is `X^{o(1)}`. Consequently, any positive lower
density construction in this fixed-library grammar requires at least
`c log X` switches for some witnesses of counts at most X, with c depending
on the library.

* Throughout, `h' >= a(h+s-d) >= max(h,s)`, because `q>=2a` and `d<=h`.
  Thus all pre-final z are bounded by X.
* With z positive, every preserve operation doubles it at least. A switch to
  a nonsingleton R likewise has `z' >= P_(R-r)(2) z >= 2z`: witness paths
  starting with r,u, a nontrivial W-path and an optional R-tail. Distinct choices
  are disjoint, with residual weights counted by the formula for two slots.
  A singleton switch obeys `z'=z+g, g'=h, h'>=s`. Any three consecutive such
  switches double z (the third adds at least the first h', itself at least
  the initial z). Therefore any three consecutive operations once z>0 double
  it. The number of operations in this part is O(log X).
* At z=0, u is a sink (any outneighbor would give a two-vertex path of positive
  residual weight). Thus `h=g=d=b`. A retained transitive block then has no
  state effect and can be dropped. Every other retained block multiplies h by
  at least 3. The first switch makes z positive. This proves the O(log X)
  count of effective operations from any seed.
* Retained operations within one nonempty run commute on the state, including
  the reset `d=g`. If the library has J blocks and L allowed new-root labels
  summed over blocks, and B effective operations, the word count modulo this
  commutation for exactly j switches is at most
  `L^j binom(B+J(j+1), J(j+1))` per seed: all retained runs are indexed by their
  multisets, and there are j+1 runs. As `B=O(log X)`, summing j≤k gives
  `X^{o(1)}` when `k=o(log X)`. Taking a sufficiently small positive
  multiple of `log X` instead of k gives an exponent strictly less than one.

There is no assertion that linearly many switches are sufficient, or that a
growing module library must satisfy the same obstruction. An exact-target
algorithm still has to choose a realizable sequence of root changes, or escape
this grammar. The reset makes the obstruction about physical continuation, not
just about fitting a congruence.

## Verification

`python replay.py` independently computes literal weighted path counts by
subset DP over six sizes for each of W and R, covering every available new root
and the retained root (2,430 comparisons, fixed seed). This is a check of the
closed formulas, not evidence for full spectrum coverage.

### The obstruction also covers every transitive width

The same `X^{o(1)}` conclusion for sublogarithmically many switches holds
even for the infinite library `{TT_l : l>=1}`. At positive z, retaining
`TT_l` multiplies z by `2^l), and switching into it increases z by at
least the factor `2^(l-1)`. By the preceding bound on the total number of
steps, the sum of lengths of all pre-final blocks after z becomes positive
is O(log X). The final block has the same bound from
`h'=b+2^l z >= (2^l-2) z`, since `b+2z=h+s-d>0`.
At z=0, all retained transitive moves are inert. On the first switch the
formula simplifies to `h'=g'=h`, `d'=0`,
`z'=(2^(e+1)-1)h`, where e is the number of vertices strictly after the
new root; the number before it can be normalized to zero. If it is not
the last step, e is O(log X) since pre-final z is bounded. If it is the
last step then its count is just h irrespective of e.

Within each retained run, concatenate all widths into one width; a
switch is determined by the two lengths before and after the new root.
Therefore the resulting finite count bound via nonnegative compositions
is `binom(O(log X)+3k+O(1), 3k+O(1))` times the number of root seeds.
This is subpolynomial when `k=o(log X)`. Allowing much longer individual
dyadic modules does not remove the cost of repeated root changes.

### More generally: arbitrary modules with bounded-size strong components

Allow an infinite module library whose strong components all belong to
a fixed finite set, with arbitrarily many components per module. The same
switch lower bound holds. Here are the extra sufficient statistics. In the
ordinal SCC chain decompose a switched module as `R=A -> C -> B`, with the
marked r in the single SCC C. Then
```
s_i(R,r) = P_A(i) s_i(C,r) P_B(i+1),
d_i(R,r) = P_A(i) d_i(C,r) P_B(i-1),
f(R,r) = [A is empty] f(C,r) H(B).
```
Path weight from A factors into the complement because r cannot reach it.
The optional last slot in B gives `i+1` colours for unsigned weights,
and contributes `-1` instead of `+1` for signed weights, by odd block
parity. Both the full and deleted profiles used by the transfer factor
over SCCs. Thus the required statistics of a switched module are exactly
determined by C and its root and by the two SCC multisets on either side.
Retained runs likewise compress to their SCC multisets.

When z is positive, a module with c SCCs charges at least c bits of growth
for a retain and at least `c-1` for a switch (delete r, then count the
remaining SCCs and use `P(2)>=2^(number of SCCs)`). The final count bounds c
by `h'>=(2^c-2)z`. As before all pre-final active operations number O(log X),
so total active SCC mass is O(log X). At z=0, first switches have
`h'=h H(R)`, `d'=0`, `g'=h H(R-r)`, `s'=2h s_1(R,r)`.
Arbitrarily many singleton SCCs in A are inert, hence can be deleted;
nontrivial SCCs in A are charged to `H(R)`. For a non-final first switch,
`z'>=h s_1(R,r)>=P_B(2)` charges all SCCs in B. For a final first switch,
singleton SCCs throughout R can simply be discarded when counting outputs.
Retained modules while z=0 similarly ignore all singleton SCCs.
Now replace the dimension `3k+O(1)` in the transitive composition bound by
`J(3k+1)`, for J allowed SCC types, and add at most an exponential factor
in k for marked-C choices. The resulting bound is still `X^{o(1)}` for
`k=o(log X)`.
