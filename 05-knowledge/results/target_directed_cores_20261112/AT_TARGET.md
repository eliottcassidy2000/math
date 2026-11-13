# Crossing the first ready hole with a changed interface

**PROVED constructions / FINITE-EXACT marked obstructions.**
This note uses the complete exclusion in [TARGET_INVERSE](TARGET_INVERSE.md).
The distinction between target attainment, a root's positive defect, and
a physically reusable endpoint is essential throughout.

## 1. The least positive defect is two

Put N=53,084,231=405·2^17+71. The constructor in
[DEFECT_TWO](constructors/DEFECT_TWO.md) supplies a strong27-vertex graph
with an indegree-one root and state

    (H,s,q,F)=(53084231,16819383,16819385,53084233).

Here q=H(T−root), s counts full paths starting at the root, and F counts
ordered covers whose first slot is nonempty and ends at the root,
with the second slot possibly empty. Thus q−s=2
and F−H=2. Since no unit-ready source exists at N, this is the minimum
**positive** indegree-one defect. Defect zero is not excluded.

The two-entry core is (6,4,8,3), with entry set{1,2}, and outer integer
1133276. Its whole count is 5g(x)+43x+3(x mod2)+18. The two entry
channels require their joint correction; replacing them by two independent
one-port counts would give the wrong result.

The source also extends to an actual fixed-defect family, with its physical
root and predecessor retained:

    H_l=54397248·2^l−1313017,
    q_l=18132416·2^l−1313031,  s_l=q_l−2,  F_l=H_l+2,
    l≥0.

Its NR exit is balanced (3H_l+4,H_l,H_l,3H_l+4), always primitive because
gcd(H_l,4)=1. The inherited balanced bridges therefore enter indefinitely
reusable controllers. This is a positive family with full physical entry;
it does not change +4 to the requested +2.

## 2. Two adjacent attained counts have opposite unit-ready behavior

At l0 the preceding exit has H=159,252,697. The complete inverse and
the independent direct-word audit both exclude a unit-ready source at
that value. Their hereditary universe has16,075,330 cores; the independent
engine tests233,412,223 outer addresses. See the receipts and
[independent review](ports/COMPLETE_INVERSE_AUDIT.md).

In contrast, the desired value

    3N+2=159,252,695

has a ready source over the singleton core, with outer integer45157534,
that is the low-bit-first word

    01111001001100001000110101.

The full state is

    (159252695,68937626,68937627,159252696).

This source has29 vertices. The exact singleton-core formula is
Q(x)=3x+g(x)+(x mod2)=g(4x), and g(45157534)=23780093.
The full-count formula and unique-start bijection prove its state;
no29-vertex subset census is claimed.

## 3. A direct selective-port macro retains every old vertex

Use the inherited raw405-ray source at k17. Start with C3=(2,4,1),
root0, attach the funnel word1^17 0001, then append its ready root.
This is the unit-ready135·2^17+23 source in
[MODULE_RETURN](../cycle_opening_six_return_20261112/compiler/MODULE_RETURN.md).
Two newborn-root births produce the raw28-vertex source at N.
Its transitive module consists of labels3,...,19; its root and unique
predecessor are27 and26. Its full state is

    (405D+71,135D+24,243D+21,513D+68),  D=2^17.

Construct the target of §2 in append order: singleton0, tail1,...,26,
predecessor27, ready root28. Match the old module3,...,19 in increasing
order to tail1,...,17. Match old26→new27 and old27→new28. Match the
nine other old vertices, in increasing order, to the first nine unused
target labels. Insert the one remaining target vertex and give all
nonroot comparisons their declared target values.

The [literal certificate](gap_bridge.json) checks:

* all28 old vertices survive and exactly one vertex is inserted;
* the entire old root star, including its physical sole predecessor,
  survives;
* every internal arc of the original TT17 module survives;
* exactly75 old nonroot pairs reverse.

Some predecessor-to-module comparisons change selectively. Consequently
the old module need not remain homogeneous relative to that port.
This is exactly where the operation leaves the fixed-block obstruction's
domain. The complete endpoint has the desired target3N+2 and is ready.
Intermediate flips are not claimed count-neutral, and the macro is not
the standard two-birth operation. This one explicit bridge is not an
all-k replacement law.

## 4. The old defect-three source also has a primitive exit

The incoming [changing-core packet](../changing_core_returns_20261113/returns/PROOF.md)
already constructs N with state

    (N,5319648,5319651,N+3).

Keeping its final root and applying NR gives the primitive balanced state
(159252699,N,N,159252699), because gcd(N,6)=1. The inherited C3 guard
passes directly. This uses the source's declared root, avoiding the
different predecessor preparation whose particular guard failed there.
The source and its count are inherited; this is a short application of
the previous session's controller-entry mechanism.

The lower positive defect two supplies a different primitive exit, while
§3 supplies the exact +2 target. None of these facts asserts that a
unit-ready representative exists at N. Their separation is the mechanism
needed for a target-directed proof: choose both the count and the
interface appropriate to the next operation, with a physical conversion
between them when it is actually available.
