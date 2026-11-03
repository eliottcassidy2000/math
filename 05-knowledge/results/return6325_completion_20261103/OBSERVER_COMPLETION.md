# Completion number and response-holonomy filtration

**Definitions and PROVED finite examples.** The general linear-observability
construction is inherited from
[Session171 framework](../session171_knots_information_20261103/framework/NOTE.md)
and [synthesis](../session171_knots_information_20261103/synthesis/NOTE.md).
The new content here specifies actual tournament contexts, distinguishes
state and law recovery, and gives exact completion numbers and an all-size
depth family. No general minimum native-complement carrier is asserted.

## OC1. Specify the task and the legal context bank

This specializes Session171's explicit observer-completion proposal to a
bank of actual tournament contexts and a desired task, including nonlinear
point identification as well as its original linear-law problem.

Let X be a source set, B a retained baseline, V the desired distinction,
and C a bank of actual contexts with responses R_c. Define the nonadaptive
completion number kappa_C(V|B) as the least size of a set D subset C such
that (B(x),(R_c(x))_(c in D)) determines V(x) for all x in X. Infinity is
allowed. Record maximum context depth and physical cost separately.

For finite X this is set cover on pairs x,y with B(x)=B(y), V(x)!=V(y):
context c covers the pair when R_c(x)!=R_c(y). This equivalence proves
minimum claims by a complete bank scan. For just one fixed pair, the number
is trivially 0,1 or infinity; the interesting problems concern a whole
state family, a weighted source law, or an operation-closed carrier.

An integer-valued response is not a one-bit query. Adaptively chosen
contexts are a different decision-tree problem. A decoder nonlinear in
responses may identify a point of X while failing to recover arbitrary
weights on X.

## OC2. Eight actual pair-orientation states

Use T=C3[TT2,TT2,TT2] on fixed names 0,...,5, with rows
`14,12,56,48,35,3`. Independently reverse the internal pair arcs, giving
eight states. The baseline retains the C3 quotient and the three unordered
pair fibres. Their internal bits are lost. All eight ordinary tournament
classes and unmarked cover profiles agree.

Allow one birth beaten by a nonempty proper subset of old vertices, and
observe H. There are exactly 62 cuts; every child stays strong. Exhausting
this bank gives kappa=2 for full state identification. No single context
has eight different responses. Two sufficient cuts, in bit-mask notation,
are:

|Cut|Responses in state order 000,001,010,011,100,101,110,111|
|---|---|
|1|95,109,95,109,95,109,95,109|
|20|123,123,131,131,139,139,151,151|

Their ordered pair determines all eight states. If a context may contain
two successive births, cuts 1 then 20 give

    375,417,415,467,407,463,453,523,

all distinct. Completion number is then 1, with attachment depth 2.
`completion.py` checks the complete one-birth bank and set-cover minimum;
`halving_response.py` counts each two-birth child directly.

Yet the full 62-by-8 one-birth response matrix has rational rank only 7.
Its kernel is the alternating parity vector. Uniform positive laws on the
four even and four odd states give identical expected single-birth
responses. The two-birth context has respective sums 1758 and 1762,
exposing that joint distinction. Seven independent one-birth rows plus this
row have determinant -143302656. This proves rational recovery of all
weights; it is not an integral unimodularity or all-moduli claim.

## OC3. Response-holonomy filtration

On V=Q^X let N_d be the common kernel of the baseline and all allowed
responses of depth at most d. If the context bank is closed under a
transport group G, then N_d is G-invariant. Define

    K_d = {g in G : (g-I)V subset N_d}.

Thus K_d is the kernel of the action on the observable quotient V/N_d;
the descending N_d track hidden laws, and the descending normal subgroups
K_d track invisible transport. In a linear controller language these are
the inherited kernels of all C A_w with |w|<=d. This algebraic definition
does not certify that an arbitrary matrix A_w has a realizable context.

For OC2 let G=(C2)^3 translate the pair bits. All fixed-cut banks are closed
under this action. At depths 0,1,2 the response ranks are 1,7,8, so the
blind dimensions are 7,1,0. The corresponding |K_d| are 8,1,1: a local
single-pair read detects every nontrivial translation already at depth 1.
The remaining parity line can still hide a distribution. Faithful
transport therefore does not imply full joint-law recovery.

This G is the explicit pair-flip action, not the nine-vertex reversal-frame
Bruck loop. Each pair arc reversal can be compiled into legal native
commands with helpers when its required path guard is supplied; the
definition of K_d does not charge that transport as a free physical move.

## OC4. What halving gives for free depends on its receipt

In C3[TT4,TT4,TT4], move a designated pole from the beginning to the end
of each four-vertex order, independently. If the baseline keeps only the
final C3 quotient and each unordered four-vertex fibre, one actual birth
with cut 275 distinguishes all eight states:

    1964273,2134907,2230415,2426313,
    2188715,2382133,2490921,2709603.

The intermediate halving receipt is cheaper. The local unordered pair
partitions are respectively {{0,1},{2,3}} and {{1,2},{0,3}}. They already
distinguish the two pole positions in each block. Retaining them gives
completion number zero. `completion.py` records all eight partitions.

This example deliberately changes the baseline: the one-query conclusion
omits intermediate fibres, while the zero-query conclusion retains them.
The older intrinsic decoder already requires those fibres; its interface
is not refuted. For TT2 fibres, by contrast, the unordered partition is
the same in both orientations, exactly the missing bit of OC2.

## OC5. Hostile source and unbounded depth

For the actual twelve-vertex 81/27 pair, H=462915 on both sources. One birth
with only old vertex 8 beating the new vertex (mask 256) gives H=1061775
and 1395279. Of all 4094 nontrivial cuts, 3584 distinguish the pair; a
singleton cut is the minimum cut cardinality. This is one additional
context for this named pair, not a universal classification result.
Both thirteen-vertex counts are independently checked by subset DP.
Automorphism orders and earlier joint-port ranks are inherited unchanged.

Uniform module substitution remains subject to the older hostile equality.
The displayed birth distinguishes individual old names and breaks that
uniformity. The full source rows, cuts and table hashes are retained.

[HB4](HALVING_BOUNDARY.md) supplies actual strong families whose positive
parity laws require exactly k attached vertices for every k. This supplies
a precise growing-memory mechanism after modulewise halving while keeping
two external port names. It concerns the birth-only menu; the completion
number under arbitrary selected complements and re-markings remains OPEN.
