# The native zero operation cannot restore this marked ready interface

**PROVED:** a forward-invariant obstruction for the restricted two-link
N/Z machine. **FINITE-EXACT:** independent literal controls below. The
universal Z operation is inherited from incoming
[NATIVE_ZERO.md](../../odd_tables_signed_carries_20261111/history/NATIVE_ZERO.md),
merged at a18596e598; this note adds the restricted-machine consequence.

Let the marked root r have sole in-neighbour a, and let a itself have
indegree one. Deleting a and r leaves a nonempty tournament X. The graph
is literally two successive births from X, because a loses only to its
predecessor in X and r loses only to a. The inherited birth formulas give

    h=H(X)+2F_X,       s=F_X,

where s counts full paths starting at r. Consequently

    h=H(X)+2s >= 2s+1.

The strict inequality uses H(X)>=1. The predecessor of a is distinct
from a,r, so X cannot be empty.

Write b=F-h. Operation N appends a vertex beaten only by r, followed by
one beaten only by that newborn, and marks the last newborn. Operation Z
instead makes the first newborn beaten only by a. Both preserve every
old arc and preserve the two-link condition. Their exact maps are

    N(h,b,s)=(3h+2b, h-s, h+b),
    Z(h,b,s)=(3h, 2s, h).

Suppose s>=2 and b>=0. Under N,

    b'=h-s >= s+1 >= 3,       s'=h+b >= 2s+1 >= 5.

Under Z,

    b'=2s >= 4,              s'=h >= 2s+1 >= 5.

For either operation, `h'-2s'=h>=1`. Thus the hypotheses persist, and
induction proves that **every nonempty word in N and Z has marked defect
at least 3**. It never returns this marked interface to defect 1 or 2.
The s>=2 boundary matters: the directed triangle with marked state
(h,b,s)=(3,0,1) has both N and Z outputs (9,2,3).

For the actual Q_k endpoint after its two standard births, the earlier
[target-return construction](../reset/TARGET_RETURN.md) gives, with
D=2^k and k>=1,

    (h,b,s)=(891D+215, 192D+47, 297D+72).

At k=1 this is (1997,431,666). Its literal next N state is
(6853,1331,2428), while Z gives (5991,1332,1997). The invariant applies
to every k and every subsequent N/Z word. Unlimited exact zero digits
therefore coexist with this marked readiness obstruction.

The restriction is essential. The marked root changes only to the last
newborn prescribed by the machine; Z consults its unique predecessor.
The result does not inspect all other roots of an output graph, change
orientation gauge, permit source/core replacement, or exclude a different
operation. In particular it does not invalidate the separately proved
same-count replacements or the finite defect-three escapes.

`incoming_machine.py` checks all rooted labeled X of orders 1 through 3,
rooted TT4 and S4, every N/Z word of lengths 1 and 2, and the actual Q_k
endpoints for k=1,2. It counts the resulting literal graphs and all marked
coordinates independently, verifies retained arcs and two-link conditions,
and includes the s=1 hostile control. Run with `python -B` or `python -O -B`;
`incoming_machine.json` retains the complete finite universe and source
hashes. The proof above, not the finite controls, supplies arbitrary word
length and arbitrary k.
