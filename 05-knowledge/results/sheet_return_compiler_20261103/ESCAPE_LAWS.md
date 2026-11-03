# Two literal all-size escape laws and the loss of the dyadic return

**PROVED with exact finite transfer certificates.** All vertices and arcs are
literal; no scalar-equivalent replacement of a supplied source is allowed.

## 1. Source, paths, newborns, and the preserved root

Use the rigid quotient Q=(2,4,9,19,39,15) from
[WINGS W1--W5](../defect_offset_dynamics_20261103/WINGS.md), with transitive
block sizes `(3,1,1,1,1,n)`, n>=2. A vertex has address (block,index).
Let T_n be this strong source. Its inherited count is

    t_n = 18*2^n + 17.

Fix r=(0,2), the old canonical root0 of the H161 control. Native complement
along a Hamilton path preserves its consecutive arcs and reverses all others.
Both paths below are legal for every n in the stated range.

* A_n: use the increasing block/address order as the path; complement; add z
  with incoming old vertices `{(0,0),(0,2),(2,0),(5,n-2)}`.
* B_n: use `(0,1)`, then blocks1,2,3,4, then block5 increasingly, then
  `(0,0),(0,2)`; complement; add z with incoming vertices
  `{(0,1),(1,0),(2,0)}`.

All other arcs incident to z point from z to the old vertex. Thus the old
root and every path endpoint are named; B revisits block0. At n=3 these are
exactly the published 485/487 escapes, transported from canonical labels.
The targets have n+8 vertices.

Write h_j for the inherited interval-tournament Hamilton count:

    h_0=h_1=h_2=1;   h_j=h_(j-1)+h_(j-2)+h_(j-3), j>=3.

The all-size formulas are

    H(A_n) = 97 h_(n+1),
    H(B_n) = 54 h_(n+1) + 62 h_n + 31 h_(n-1).       (1)

In particular both obey the cubic `f1(X)=X^3-X^2-X-1`.

| n | source t_n | A_n | B_n |
|---:|---:|---:|---:|
|2|89|291|255|
|3|161|485|487|
|4|305|873|889|
|5|593|1649|1631|
|6|1169|3007|3007|
|7|2321|5529|5527|

The equal scalar at n=6 does not equate the sources, their root responses,
or their next growth responses.

## 2. Why the transfer is exact

Reverse the order of block5 for counting. Its internal arcs are now all
forward except the backward arcs between consecutive vertices. Every
exterior vertex has uniform adjacency to its interior, except the explicitly
named first/last boundary vertices. Put exterior vertices pointing to a
generic interior vertex before the block, and the others after it. The first
two and last two block vertices contain every exceptional boundary incidence.
All subsequent interior insertions have the same paired-component transition.

Add a universal source and sink to turn an unrooted Hamilton-path count into
one fixed source-to-sink query. Process the resulting ordered tournament with
the inherited endpoint-expiry rule. A component retains its ordered head/tail
pair. A head with no future backward neighbour is rejected; an expired open
tail becomes anonymous. Never join both ends of the same component: that
would close a cycle instead of a Hamilton path.

For these two sources, after the prefix there are exactly four reachable
states. Write s for the fixed source, a for an anonymous tail, and v for the
most recently processed bulk vertex. Their component lists are

    ((s,a),(v,a)), ((s,a)), ((s,a),(v,v)), ((s,v)).

The bulk matrix, common terminal column, and initial rows at n=4 are

    B = [[1,1,0,0], [0,0,1,1], [1,1,0,0], [0,0,1,0]],
    beta = (3,3,3,2)^T,
    alpha_A = (163,97,31,0),
    alpha_B = (85,85,85,62).

Consequently `H(A_(4+j))=alpha_A B^j beta`, and similarly for B.
`bulk.json` records every prefix and suffix transition, including all
individual backward incidences. Directly from the displayed adjacency rules,
inserting a further interior vertex only repeats the one bulk transition;
this is why a finite matrix proves an all-size statement.

The first and third states have equal future transitions and equal terminal
weights. Merge their amplitudes. This gives the positive three-state matrix

    C = [[1,1,0], [1,0,1], [1,0,0]],
    a_A=(194,97,0), a_B=(170,85,62), b=(3,3,2)^T.

Its characteristic polynomial is f1. Equivalently the exact identity
`f1(B) beta=0` holds before merging. Three initial values prove (1).
The size-three Hankel determinants from n=4 are -3650692 and -171676,
respectively, so neither scalar recurrence has smaller rational order.
This is a polynomial-annihilation certificate, not recurrence fitting.

## 3. Strongness at every size

At n=3, directed Hamilton cycles are, in the natural vertex labels0,...,10,

    A: (7,8,10,1,2,3,0,4,5,9,6),
    B: (7,8,9,0,2,1,3,4,10,5,6).

Both use `(5,0)->(5,1)`. Insert the next block vertex between these two
addresses, shifting later indices. The old adjacent edge reverses, and the
cycle replaces it by the two forward edges through the new vertex. All other
old incidences are unchanged, including the tracked penultimate newborn
incidence on A. The cycles therefore extend inductively. The n=2 bases are
checked by literal cycle enumeration. This proves strongness for all n>=2.

This target-growth operation includes an old-edge reversal. It is not an
induced birth. Each target nevertheless has the explicitly supplied legal
native-complement-plus-birth realization from T_n. No helper deletion is
hidden in that realization or in the count formula.

## 4. Complete failure boundary for the two offset lifts

Among n>=2 and c in{2,4}, the identities `H(A_n)=3t_n+c` or
`H(B_n)=3t_n+c` occur exactly at `(A,3,2)` and `(B,3,4)`.

For a proof beyond finite testing, set `g_n=54*2^n+53-H(A_n)` or use B.
Both satisfy, for n>=5,

    g_n=g_(n-1)+g_(n-2)+g_(n-3)+54*2^(n-3)-106.

At n=4,5,6 the three g values are positive on both branches; the forcing
term is positive thereafter. Thus both outputs are strictly less than even
the smaller target3t_n+2 for all n>=4. Check n=2,3 using the table.

Mechanism: the retained source growth has a dyadic eigenvalue2; the actual
complement replaces it with a tribonacci boundary transfer for which
`f1(2)=1`. Indeed any sequence annihilated by f1 cannot agree with
`a*2^n+b`, a!=0, for five consecutive n: applying f1 on two consecutive
four-term windows would require both `a*2^n-2b=0` and
`a*2^(n+1)-2b=0`. A change of path/run policy or boundary operation is needed.

A_n is always composite since97 divides its count. For B_n, divisibility
by7 occurs exactly at `n=1,4,5,10 (mod16)`; for A_n it occurs exactly at
`n=8 (mod16)`. The certificate is the exact cubic clock modulo7, whose
sixteenth power is4I, followed by the sixteen residue checks. These are
divisibility classifications, not assertions of primality in other phases.

## 5. Both actual escape outputs have another pair of offset exits

The following **FINITE-EXACT constructive** programs continue the literal
eleven-vertex A_3/B_3 targets. Labels are their natural labels, not a freshly
selected isomorphism representative. A mask marks old vertices pointing to
the new vertex. Every final target is strong and has twelve vertices.

| Source | Target | Selected native path (or direct birth) | Birth mask | Intermediate H | Zero-guard root |
|---:|---:|---|---:|---:|---:|
|485|1457|(1,2,4,5,9,6,7,8,10,3,0)|194|303|6|
|485|1459|(3,4,5,9,6,7,8,10,1,2,0)|720|417|6|
|487|1463|(0,2,1,3,4,10,5,6,7,8,9)|2037|803|6|
|487|1465|direct birth|1823|487|1|

All rows, complete root-defect arrays, retained-root integer observers and
F49 phases are saved in `next_escapes.json`. Independent path enumeration
checks all four counts. The physical old root2 has respective phases
`(2,2,6),(1,2,1),(5,1,6),(3,5,4)` on the four rows. No phase neutrality
is asserted and no unnamed re-marking is hidden.

At each displayed zero-guard root, the ordered-two-bin count ending at that
root equals the full H. The inherited universal reset and canonical C/D
returns therefore continue every child h to all counts

    3^(a+1) 5^b h, a,b>=0,

with the prior literal guard programs. This credits the earlier
[canonical-return integration](../defect_offset_dynamics_20261103/INCOMING.md);
the new contribution is entry into those continuations from all four children.

This is a deeper source-preserving branch bank, not a reusable offset theorem.
For example the cheap487->1465 birth points to precisely the old boundary
triple `{(3,0),(4,0),(5,0)}` and loses to every other old vertex. Repeating
that same rule, with the first newborn now an old vertex, gives3745 rather
than4399=3*1465+4. The target needs a new controller or source-return rule.
