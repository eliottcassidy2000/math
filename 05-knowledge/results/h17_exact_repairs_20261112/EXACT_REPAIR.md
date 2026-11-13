# A two-arc repair at the same next target

**PROVED in note.** This generalizes the newly exposed one-step-hidden
H15 repair to arbitrary left factors, classifies its named-pair divisors,
and gives an infinite explicitly accepted return family. Chan–Pak and
arbitrary prescribed-target repair remain **OPEN**.

## 1. The physical five-vertex switch

Use outgoing rows and physical root0:

    K=(16,9,3,5,14),     K'=(18,24,3,5,12).

Reverse the two arcs1→0 and4→1. Both tournaments are strong, and their
rooted states are

    (H,s,q,F)(K,0)  =(15,3,3,21),
    (H,s,q,F)(K',0) =(15,3,5,21).                    (1)

Here s counts paths starting at the root, q counts paths after deleting
it, and F counts ordered spanning covers(P,Q) with P nonempty ending
at the root and Q possibly empty. Independent permutation/cut enumeration
certifies(1). K is an instance of the old cyclic triangle wrapper, rooted
at its sink-before-the-new-root; K' is the regular five-tournament.
The sources and ensuing balanced family were already used in
[RESET_JUNCTION J2–J4](../coverage_machine_induction_20261103/RESET_JUNCTION.md).
The present interface retains a different root and uses two literal arc
changes: both H15 and the next target H33 stay fixed.

The root-child move adds x with `root→x→all other old vertices`. R keeps
the root; N moves it to x. The inherited state maps are

    R(h,s,q,F)=(h+F−s,q,q,2F−q),
    N(h,s,q,F)=(h+F−s,h−s,h,h+F).                   (2)

Thus R sends(1) to `(33,3,3,39)` and `(33,5,5,37)`.
The new child has exactly the same incidence in both graphs. The physical
root and child are retained, while the old root's indegree is allowed to
change. This is not an indegree-one-preserving replacement at the old root.

## 2. Exact classification with any nonempty left factor

Let L be any nonempty tournament. Write

    h=H(L),       f=F_L(2),

where F_L(2) counts ordered two-slot spanning path covers allowing empty
slots. Make `T=L→K` or `T'=L→K'`, keeping the named root in the five-block.
Their states are exactly

    z=(15h,0,3h,21f),       z'=(15h,0,5h,21f).      (3)

Indeed H and deletion H multiply under ordinal join. A full Hamilton
path cannot start in K when L is nonempty. In the rooted two-slot cover,
the L parts are any ordered two-slot cover and independently the K parts
are a rooted cover ending at0. This proves F=21f without an inferred
unmarked-to-marked identification.

After one R, the two actual strong graphs have the SAME order and count

    M=15h+21f,
    (M,3h,3h,42f−3h),    (M,5h,5h,42f−5h).        (4)

They differ in just the same two old arcs, and retain the same root and
new child. Strongness is physical: the child reaches L, L reaches all
of K, K is strong and reaches the root, and the root reaches the child.

The inherited root-child interface theorem
[ROOTED_RESETS R3/R5](../two_port_return_induction_20261103/ROOTED_RESETS.md)
identifies the exact continuation divisor at the named final pair with
the gcd of the pre-birth rooted coordinates. Therefore the two divisors
are

    D_old=3 gcd(h,7f),          D_new=gcd(5h,21f).   (5)

The coordinate gcds in(4) are the same respective numbers. These are
exact formulas, even when the new interface is not primitive. In
particular the hypothesis

    gcd(5h,21f)=1                                  (6)

gives an actual exact-target repair `D_old=3 → D_new=1`.
It implies gcd(h,21f)=1 and5 does not divide f; no assertion that every
left factor meets(6) is made. If L is an ordinal product of factors L_i,
then h=product H(L_i) and f=product F_(L_i)(2), so(6) requires every
cross gcd `gcd(H(L_i),21F_(L_j)(2))=1` and `5∤F_(L_j)(2)`.
Individual primitivity tests without the cross terms are insufficient.

## 3. An explicit reusable continuation for every admissible L

Start from the repaired state z' in(3). Put

    d=21f−5h,       delta=21f−20h.

Both are positive: cutting each Hamilton path shows f≥2h, and in fact
f≥(|L|+1)h. For every k≥1, the inherited R^k N R theorem gives

    t=20h+(2^k−1)d,
    b=30h+(3·2^k−1)d,
    b−3t=2delta,
    gcd(b,t)=gcd(delta,2^k d+5h).                   (7)

Let g be the odd part of delta, and use period1 when g=1. Every positive
multiple k of ord_g(2) makes(7) primitive. In fact modulo g its second
gcd argument is `d+5h=21f`. Any odd common prime of delta and21f would
divide20h. Condition(6) excludes primes dividing h and excludes5, leaving
no such prime. The returned t is odd, so powers of two cause no defect.
This is an effective finite clock calculation, with no unproved
number-theoretic assumption or asymptotic prime assertion.

For every k, `b−2t ≡ 2^k h (mod3)`, nonzero under(6). Therefore, at a
primitive clock, the inherited C3 balanced return `(b,t)→(61b−2t,9b)`
remains primitive and has `3|t', 3∤b'`. If its projective ratio modulo11
is7, one further C3 makes the ratio3. A C5 then enters the inherited
eight-state repeatable controller. These are actual entered sources,
using [NATIVE_SELECTION section7](../source_selection_closure_20261103/NATIVE_SELECTION.md)
and the entry bridge in RESET_JUNCTION J4. No finite-field tuple is being
assumed physically realizable.

The first R in(7) is the already repaired common target in(4). All further
graphs are strong. The active root is kept during the R steps and moved
only by the named N command. Reusability refers to the entered balanced
controller, not to repeated neutrality of the old two-arc switch after
arbitrary later selective attachments.

## 4. A small explicit clock: the transitive prefix family

Take L=TT_a, a≥1. Then h=1, f=x=2^a and(6) always holds. The exact target
repaired in(4) is

    M_a=21·2^a+15                                  (8)

on a+6 vertices. The same two old arcs give D_old3 and D_new1 at the
same physical pair, at every a. The first targets are57,99,183,351.

One need not compute an order modulo the defect: k=a is always primitive.
Since5∤(21x−20), and

    xd+5 ≡5(3x+1) mod(21x−20),
    (21x−20)−7(3x+1)=−27,      3∤(3x+1),

the gcd in(7) is1. The completed return has2a+7 vertices and

    b=63x²−36x+35,       t=21x²−26x+25.             (9)

Its first counts are215,899,3779,15587. It has the accepted C3 entry just
proved. The unprepared source has an invariant factor3 at every R^k N R
clock, since both its defect21x−18 and its d=21x−3,q=3 are divisible by3.

The empty-prefix case is separate: s=3 rather than0. It yields the
H33 repair in section1 and the inherited balanced family
`b=48·2^k+11,t=16·2^k+1`. Substituting a=0 into(3) is invalid.

## 5. What this closes, and what it does not

This is an exact-target, retained-pair repair theorem for a specified
physical family, with a complete arbitrary-left-factor criterion and
an infinite accepted continuation. It goes beyond changing H17's
future count: here the first next count itself is unchanged.

It does not show that a prescribed arbitrary odd N admits
`N=15H(L)+21F_L(2)` with an admissible L, nor that every later balanced
source has the same repair module. All counts in(8) are multiples of3.
Exact positive target supply and unrestricted recursive repair remain
the Chan–Pak gaps. The all-word monotonicity obstruction in the history
lane explains why one cannot simply commute the earlier H17 preparation
past a birth; the pure deletion-coordinate change in(1) supplies a
different legal direction.

Run `python3 -B 05-knowledge/results/h17_exact_repairs_20261112/exact_repair.py`.
The receipt separates literal finite controls from the universal proofs.
