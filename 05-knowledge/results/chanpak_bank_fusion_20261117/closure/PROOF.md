# Positive recursive fusion does not by itself create a dense bank

The variable-height and rooted banks were designed to be *composed*.
Here is a genuinely nonlinear test of that idea, rather than just
the statement that products do not create new primes. A library
with a small positive weighted count remains thin under arbitrarily
deep ordinal and some cyclic compositions, including operations with
unbounded transitive clock lengths.

Use alpha=199/200. Count a physical witness history as one object, even
when another history has the same tournament or count. The weighted
sum of a family is the sum of H^(-alpha) over its histories. Any
bound on this sum also bounds the number of different counts.

## Two library inputs

For the [migration repair bank](../coupling/PROOF.md), the counting
proof bounds histories, not merely distinct outputs: it sums labeled
migration words, then the integer parameters and the four endpoints.
Consequently the history count is O(X^(99/100) log X). Summing over
dyadic ranges shows that its alpha-weighted sum converges. We can
choose a fixed cutoff beyond which the weight of its tail is at most
1/32. This only bounds a stated library, not all admissible migration
stages.

The original [variable back-edge bank](../../chanpak_joint_placement_20261117/physical/PROOF.md)
is small enough to include *without* such a cutoff. Use all even m>=2,
all dyadic G>=4 and 1<=z with |Delta|<=PG², counting all complementary
suffixes and four endpoints. Its exact bias satisfies B_r>=16PG²:
in its expansion the corrections beyond 16PG² are nonnegative in
aggregate, since A<=147P/256, C>=0 and P>=320h. The shift formula
forces y>=5PG²(z+1)/h. Both endpoints have H>=Puw>=256PG²y².
There are at most PG²/h possible y's for a fixed parameter quadruple.
At exponent 3/4 this gives the upper bound

    2·6400^(-3/4)
      · sum_(m even >=2) h^(1/2) P^(-5/4)
      · sum_(G dyadic >=4) G^(-3/2)
      · sum_(z>=1) (z+1)^(-3/2).

Use h<=7·2^m, P>=512·3^m, the dyadic sum <1/4, the z sum <=2,
and sqrt(2)/3^(5/4)<1/2. Even replacing the even m sum by all
m>=2 and weakening sqrt(7)<3, 512^(-5/4)<1/512, this is
less than 1/32. Raising exponent from 3/4 to alpha can only lower it.

Finally allow an explicit fixed seed H=129 (the extreme-back-edge
tournament of spine parameter 7). Its weight is less than 1/32 since
129^(3/4)>32. Hence the total library weight is below 3/32.

## A recursive mixed grammar

Leaves can be any library witness. Permit the following internal nodes,
in arbitrary combinations and with no bound on depth:

- Ordinal join of two children.
- Cyclic substitution C3[B,C,1] of two children and a singleton.
- Cyclic substitution C3[B,C,D] of three children.
- C3[B,C,TT_j] for any integer j>=8 (two children).
- C3[B,TT_j,TT_j] for any integer j>=8 (one child).

Each is a literal tournament. Ordinal join multiplies H. Cyclic
substitution has at least three times the product of the children's
counts, by following the three one-pass quotient Hamilton paths.
With one child slot TT_j, a further way is to split TT_j into an
ordered pair of nonempty paths and visit the quotient twice at
that slot. Every nontrivial subset choice works and there are
2^j−2 of them. Thus either clocked constructor has count at least
(2^j+1) times the product of its non-clock children. This is just
a lower bound; omitting other physical paths helps make an upper
bound on weight.

Let T_d be the total weight of represented histories of depth <=d.
The displayed inequalities give

    T_(d+1) <= 3/32 + lambda T_d
                 +(1+3^(-alpha)+lambda) T_d²+3^(-alpha) T_d³,
    lambda = sum_(j>=8) (2^j+1)^(-alpha) < 1/16.

The last bound follows by summing (3/5)^j, since 2^(-alpha)<3/5.
Also 3^(-alpha)<1/2. At T_d=1/4 the right side is strictly
below 1/4. Induction and monotone convergence give a finite
weighted sum for the whole grammar. Therefore the number of
different outputs up to X is O(X^(199/200)), even with the
nonlinear glues and unbounded clock lengths. It cannot realize
every admissible count or every odd prime.

This is a scope result: enlarging the small seed catalog or
using negative contributions, tagged coordinates, or more
permissive cuts can invalidate the weighted majorant. It does
not rule out an actual selector outside the stated grammar.
It identifies why composing positive readouts of the known
sparse banks is not a substitute for a new target-controlled
parameter.
