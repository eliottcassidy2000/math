# Exact mixed-return obstructions

Status: **PROVED** for the explicitly defined catalog and operation language.
Finite matrix certificates below are part of the proof, not extrapolation
from a few large scalar values.

## T1. Objects, coordinates and the complete finite return catalog

Inherit the actual marked bridge operation and row state v=(h,x,y) from
[G1-G4](../golden_hamilton_transfer_20261103/RECURRENCE.md). Here H=h,
E=y-x counts paths ending at the port, A=h+x-y=H-E, and C=x+y is the marked
two-cover count. Counting one-backward-step defects gives C>=H. Thus every
nonempty pointed source satisfies x>=0, y>=x, y<=h+x, x+y>=h, and h>=1.

For t=(2,-1,0), use the saturated quotient

    beta(v)=(Y,X)=(y,h+2x-y)=vB,
    B=[[0,1],[0,2],[1,-1]],    tB=0.

In particular every such source has

    Y>0,       0<=X<=3Y.                              (T1)

Indeed Y>=h/2>0; X=A+x>=0; and X<=3x<=3y. These are necessary physical
conditions, not a complete realization criterion.

The catalog consists of EVERY nonempty word of length at most six in the
fixed original65 alphabet for which tP=g t with nonzero integer g, plus
the supplied s78 atom. Distinct words remain distinct presentations.

| Length | Mixed-error words | Distinct error rows | Return words |
|---:|---:|---:|---:|
|1|24|20|2|
|2|227|113|6|
|3|1381|539|16|
|4|7651|2462|53|
|5|35650|9570|167|
|6|161538|35511|541|

The first six rows give785 return words; adding s78 gives786 blocks and423
distinct full matrices. The exhaustive generator discards a row only when
all its entries have the same weak sign. Nonnegative append matrices can
never turn such a row back into nonzero mixed-sign g t. An independent
row-multiplicity DP checks all displayed counts. Every saved full matrix
and correction identity is separately multiplied from the pinned atlas.

For each block, PB=BQ with Q a positive integer 2-by-2 matrix. These exact
identities are checked on all three rows. Arbitrary concatenations of the
catalog blocks preserve this quotient. They include all the earlier named
loops, D=L2 b L2, and s78^2. They do not include every possible longer
primitive return word in the full atom semigroup.

## T2. Projective images exclude the seven proposed exact returns

Write Q=[[a,b],[c,d]]. If beta is physical, its output ratio is

    X'/Y'=(b+d r)/(a+c r),       0<=r=X/Y<=3.

The derivative has constant sign, since det Q is nonzero. Its exact image
is the interval with endpoints b/a and (b+3d)/(a+3c). The union over the
786 blocks has25 connected intervals, recorded in `verification.json`.

All seven supplied states lie in the following complementary gaps:

| State h,x,y | X/Y | Gap between catalog images |
|---|---|---|
|(31,14,25), (95,41,75), (97,40,75)|34/25|(0,637/466)|
|(287,206,239)|460/239|(16498/9119,317/161)|
|(289,191,239)|432/239|(52493/29375,2876/1591)|
|(293,210,281)|432/281|(3244213/2121069,17/11)|
|(295,216,260)|467/260|(52493/29375,2876/1591)|

Consequently there is NO nonempty catalog word R, actual nonempty pointed
source w, and integer j such that

    v(w) P_R = 3v + j t                              (T2)

for any listed v. Project to beta: the target is 3 beta(v), which has the
same excluded ratio. Before the last catalog block there is still an actual
pointed source, so (T1) applies. The last block cannot produce the target
ratio. This argument has no word-length cutoff and permits arbitrary
initial tournaments, not just parents in our 66-symbol bank.

The same proof works with any positive slope in place of3. It applies only
when the correction is annihilated by B. Transverse integer carries change
the projective target and escape this hypothesis; CARRIES studies them.

## T3. Expansion versus determinant valuation rules out finite-bank cycles

The incoming [return-bank theorem R](../adaptive_state_induction_20261103/RETURN_BANKS.md)
handles quotient returns with determinants prime to3. The following criterion
also handles the nonunit loops in this catalog.

**General criterion.** Let Q_l be nonnegative nonsingular integer matrices.
Suppose their minimum row sums mu_l satisfy

    mu_l > p^v_p(det Q_l)                            (T3)

for a prime p. No nonempty product of them has a nonzero nonnegative left
eigenvector with eigenvalue p^L for any integer L>=0.

Proof: for beta>=0, its coordinate sum expands by at least mu_l at each
factor. If a product Q has total determinant valuation E, (T3) gives an
expansion strictly greater than p^E. An eigenvalue p^L would force L>E.
But an integer eigenvalue divides the determinant of an integer matrix:
factor its monic characteristic polynomial by X-p^L over Z. Thus L<=E,
a contradiction. No estimate of an average replaces this positive-vector
bound. Nonsingularity is needed for the finite determinant valuation.

For ALL786 catalog blocks, exact arithmetic certifies the stronger margin

    mu_l >= 5 * 3^v3(det Q_l).                       (T4)

The minimum margin is5, attained by b. An empty product has eigenvalue1
and cannot equal 3^L for L>=1.

Now consider any finite collection of actual seed states v_i with edges
i->j certified by

    v_j P_R = 3v_i + m t,                            (T5)

where R is any catalog word, possibly empty, and m is any integer. A cycle
of L edges makes beta(v_i) a nonnegative eigenvector of the reverse-ordered
product of the quotient returns with eigenvalue3^L. This contradicts T3,
or the empty-product observation. Hence **every such exact finite seed bank
is acyclic** and cannot provide unbounded recursive branching.

This excludes a specified representation of recursive closure, not the
unrestricted tournament spectrum, changing-port operations, transverse-carry
rules, or an infinite parameterized seed family. Unlike the previous p-unit
argument, it continues to apply after admitting s78 and D.

Strictness matters: Q=[[2,1],[1,2]] has det3 and minimum row sum3, with
(1,1)Q=3(1,1). This is an abstract equality-case control, not an asserted
tournament block. The all-input proof is new here relative to the recovered
packets; no claim of literature priority for this elementary criterion is made.
