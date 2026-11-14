# Exact cut state and an observational test using two added vertices

Earlier [cut-vertex balance work](../../../01-canon/theorems/THM-4655-strong-tournament-balance-lemmas.md)
uses a pre/gap/mid/post transfer for distributions on full paths.
For the selective **count** and ordinal growth there is a different,
four-integer state with a short exact observation theorem.

Let B be nonempty, U the vertices beaten by the apex, and h=H(B),
N=H(B+v). Over all nonempty directed paths starting in U, sum
the Hamilton counts of their complements. Denote the even and odd
parts by E and O, and let S=E+O, d=O−E. These are weighted partial
paths, not just full Hamilton paths. Use (h,N,S,d) as the state.

Append an ordinal block R after all of B, with fixed mask V of
its vertices beaten by v. In R put a=H(R), q=P_R(2); let (s,e)
be its total and odd-minus-even start weights from V, and (t,f)
its total and odd-minus-even end weights into the complement
of V. Finally let n_R be the path count of R with the same apex
contact alone. Then the exact update is
[
 egin{aligned}
 h'&=ah,\
 N'&=aN+(n_R-a)h+tS+fd,\
 S'&=qS+hs,\
 d'&=he.
 end{aligned}
]
To prove the count update, in the odd-cycle packing expansion split
the apex cycle into three cases: internal to B, internal to R, or
crossing. The first two contributions together with the no-apex
contribution give aN+(n_R−a)h. A crossing cycle consists of one
nonempty path on each side of the seam. Their vertex parities must
agree to give an odd apex cycle. Its weight two turns the joint
parity indicator into 1+(-1)^(sum of path orders), giving tS+fd.
The residual tournament is an ordinal sum, so each complement
weight factors.

For the start update, paths starting in R contribute (hs,he).
A path starting in B has an optional unrestricted path in R.
Assign the odd packing blocks of R independently to the path or
the complement. The total is q and the signed sum is zero because
R is nonempty and each packing block is odd. This gives qS with
zero transported bias. In particular the balance d is reset, not
propagated through an arbitrary number of appended blocks.

## Four scalar observations recover exactly this entire state

Record the output of the marked B after each of the following
separate choices of suffix, with 1 for a vertex beaten by v:
* no suffix: N;
* one transitive vertex painted 0: N+S+d;
* two transitive vertices painted 00: N+3S+d;
* two transitive vertices painted 10: N+2h+2S.

For 00 the nonempty end total is 3 and bias 1. For 10 the only
eligible end is the second vertex; its end total is 2 and bias 0,
and the isolated apex block is a triangle. Solving these four
integer equations recovers N, S, d, h exactly. Therefore two
painted sources give the same scalar for **every ordinal painted
continuation** if and only if these four tiny observations agree.
This is a complete continuation invariant for this context class,
not a claim about selectively recontacting vertices in the old base.

## Inversion and independently switchable blocks

If two masks on the same appended block preserve (s,e,t,f), then
switching them changes only the output coordinate, by
h times their local output difference. Every later ordinal block
multiplies this change by its ordinary Hamilton count. Thus for a
fixed word of blocks, any number of such mask switches are exact
independent additions, with no overlap correction. If the left
and right surrounding bases have counts H_L,H_R, a local +2
becomes +2H_LH_R even with arbitrary coloring of those sides.
This is not an assertion that the resulting weighted sum has no
holes; in particular multiplicative amplification may make unit
moves impossible.

For a prescribed N, this state identifies what a suffix actually
has to supply. The exact correction is (n_R−a)h+tS+fd after
subtracting aN of the current prefix. A numerical guess for these
three reader coefficients is not a construction: they must come
jointly from one marked suffix. The dyadic inverse and the
independent switch examples give physical instances; the
quantitative sparsity result shows why holding the source too
simple cannot give universal target coverage.
