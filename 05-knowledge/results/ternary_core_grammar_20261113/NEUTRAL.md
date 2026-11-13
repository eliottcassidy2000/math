# A neutral core change that enters an unbounded renewal

**PROVED explicit macro and entry; inherited all-order renewal law.**
This establishes one new entry connection, not a universal source selector.
The core-changing macro preserves its completed count and physical names;
individual edits inside the macro need not preserve H.

Use z=(h,s,q,F), where s counts Hamilton paths starting at the root,
q=H(T−root), and F counts ordered two-path covers with the first path
nonempty and ending at the root. The readiness defect is F−h.

## N1. The count25 class obstruction becomes an endpoint preparation

The seven-vertex interval chain of two four-vertex blocks has count25 and
does not belong to the triangle-ear grammar at that order. Choose its
literal representative

    old rows = (2,8,3,69,15,31,55),   root r=6, port p=3.

Its profile is (25,5,5,25). Delete physical vertex5 and reverse precisely
the comparisons on the unordered pairs

    {0,1}, {0,3}, {1,3}, {1,6}, {2,3}.

The six retained physical names are0,1,2,3,4,6. The result is the rooted
count25 class isomorphic to rows (2,36,33,7,15,25), whose root is index5.
The root r and incoming arc p→r are retained. Both graphs are strong.
The new profile is

    (25,5,5,25)  ->  (25,3,3,27).                   (N1)

Thus h and q−s=0 survive, while the readiness defect changes0→2 and the
profile gcd changes5→1. The separate start/deletion counts are deliberately
changed. This is a nonisomorphic core replacement, not mere re-marking.
The [all-order ready classification](../two_entry_funnels_20261103/READY.md)
already contains the target class;
the contribution here is its concrete connection to the blocked count25
interval chain, with the physical edit and retained interface explicit.

Two births—first at r, then at its newborn, keeping that first newborn as
the final root—give the actual balanced profile

    (79,25,25,79).

This is the inherited identity (h,s,q,F)↦(h+2F,h,h,h+2F). In particular
the next construction starts with a literal root, not just a scalar79.

## N2. Entry remains proved indefinitely

The following renewal is inherited from
[NATIVE_SELECTION §§5–7](../source_selection_closure_20261103/NATIVE_SELECTION.md).
At a balanced source (b,t,t,b), birth at its root and mark the newborn.
Create three children there. They form a transitive homogeneous block.
Move the newest physical child to the oldest position by reversing its
two comparisons to the other children. This leaves the whole H count9b
unchanged. Birth at that selected child, mark its newborn, then birth once
more while keeping the marked root. The output is

    (b',t',t',b'),        b'=61b−2t,   t'=9b.        (N2)

This uses six permanent new vertices per cycle. The neutral preparation
has an actual final graph and root; the final graph again has precisely
the balanced interface required by the same recipe. No fixed dyadic chart
is being reused as a scalar surrogate.

For the entry (79,25), the first cycle gives (4769,711), which is coprime,
has3|711 and3∤4769. If gcd(b,t)=1,3|t and3∤b, then any common divisor of
61b−2t and9b divides9, while61b−2t≡b mod3. Hence the output is coprime
and has the same guard. This proves entry and primitivity at every cycle.
The counts begin

    79, 4769, 289487, 17572865, 1066733999, …

and obey b_(j+2)=61b_(j+1)−18b_j. Every completed cycle includes a neutral
endpoint preparation. The literal audit verifies the initial replacement,
the count79 graph and the entire first cycle through14 vertices.

The growing cores escape the finite dyadic-bank restriction. Their output
sequence is sparse and does not supply prescribed arbitrary counts. The
renewal was already known; discovering this entry does not close the
original405·2^k+71 family or the473 repair obligations.

## N3. A uniform family of contexts

Append any tournament Q after the changed core in SCC order, and set
m=H(Q), allowing the empty Q with m=1. Each of h,s,q,F at the retained
root is multiplied by m: after that root only the second cover path can
visit Q. Consequently N1 preserves25m, changes the defect0→2m, and its
two-birth entry has profile m(79,25,25,79). The same renewal yields
m(b_j,t_j,t_j,b_j) at every depth. Its profile gcd is exactly m.

This is an all-context statement for **ordinal suffixes** only. An
arbitrary uniform exterior can use several core runs and does not inherit
H-neutrality from the equality of the two core counts. The two cores even
have different orders, so their complete cover polynomials differ.
