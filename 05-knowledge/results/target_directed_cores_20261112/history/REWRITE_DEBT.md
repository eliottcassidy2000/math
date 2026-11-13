# The 462 address rewrite has a positive inner-substitution debt

2026-11-12. **PROVED:** the exact debt and its equality boundary. This is
a boundary for one specified attempt to make the inherited finite chain
recursive; it excludes neither other address rewrites nor individual
clock coincidences.

The second and third cores of the inherited
[462 chain](../../multicore_metagraph_20261113/FAMILIES.md) have chronological
funnel addresses(11,1) and(1,11,empty). Replace the implicit singleton
at the start of both addresses by the same rooted tournament U. Write

    z_U=(h,e,s,p,f,t)=(H,E(1),S(2),P(2),E(2),S(3)).

The two resulting core states are

    z_2=M_1 M_11 z_U,
    z_3=M_empty M_11 M_1 z_U.

Their respective outer tails are1^k101 and1^k011, with k≥1. In the
unit-start case S_U(1)=1 these give actual unit-ready sources with

    A_2=6(s_2+p_2),       B_2+1=3h_2−f_2,
    A_3=7(s_3+p_3),       B_3+1=5h_3−f_3.

Direct multiplication of the inherited six-state matrices gives

    A_3−3A_2=16h−76s+312p−102f−462t,
    (B_3+1)−3(B_2+1)=18h−18s−14p−2f+30t.           (1)

The subscript3 in this display labels the third core, not the slot number;
t still means the source S_U(3). The checker verifies the six columns
individually, rather than fitting actual source values.

## A positive path-cover decomposition of the second row

Define three nonnegative quantities

    D=S(2)−H,
    J=S(3)−P(2)/2,
    K=P(2)/2−E(2).

They satisfy J≥D. To see this without a numerical inequality, split a
Hamilton path immediately before the root. This identifies H with the
two-slot covers(P,Q), P starting at the root, for which Q is empty or
ends at a predecessor of the root. Thus D counts the remaining covers,
whose nonempty Q ends at an out-neighbour of the root.

Likewise, forget the order of the two slots of an arbitrary P(2) cover;
the swap action is free because the vertex set is nonempty. There are
P(2)/2 unordered covers. Split the path containing the root immediately
before it, and put the root-starting suffix in slot1, its prefix in slot2,
and the other path in slot3. This identifies P(2)/2 with the S(3) covers
whose slot2 is empty or ends at a predecessor of the root. Consequently
J counts the complementary three-slot covers. The D-covers inject into
them by leaving slot3 empty. Hence J≥D≥0.

Finally K counts unordered two-slot covers whose root-containing path does
not end at the root. It vanishes iff the root is a global sink. If r→u,
the path(r,u), together with any Hamilton path on the remaining vertices,
is one of the omitted covers; the second path may be empty. Conversely a
sink must end its path in every cover.

The offset debt in(1) is therefore

    Delta_B=2(15J−9D+K)
           =12D+30(J−D)+2K≥0.                     (2)

Equality holds iff the root is a global sink: necessity follows from K=0;
for a sink, D=J=K=0 directly from the cover descriptions.

## Consequence and scope

Every nontrivial strong unique-start source has an out-neighbour at its
root. Thus substituting **any** such source U into both inherited addresses
makes Delta_B strictly positive. The original singleton is the only
eligible unique-start seed for which this fixed second-to-third rewrite
preserves the full all-k relation H_3=3H_2+2. One cannot make this finite
chain renewable simply by grafting an arbitrary recursive one-start core
into its implicit singleton.

For a general rooted U the same linear-row calculation applies to the
shifted quantity H+S_U(1), and the natural affine step is3H+2S_U(1).
We do not label that nonunit-defect graph unit-ready.

The slope debt in(1) can have either sign. Therefore(2) does not exclude
an isolated k with (A_3−3A_2)2^k+Delta_B=0. Nor does it exclude changing
the two interior sources differently. Its constructive output is the
exact nonnegative cover debt: a successful extension must compensate
12D+30(J−D)+2K while also matching the slope and the named port interface.
An actuator of the right numerical size at another baseline is insufficient.

The [periodic-extension theorem](PERIODIC_CORE_OBSTRUCTION.md) gives a
different boundary, allowing arbitrary repeated funnel words but requiring
append-only core evolution. The present result allows the specific
noncommuting interior rewrite and identifies precisely why the most direct
seed substitution fails.
