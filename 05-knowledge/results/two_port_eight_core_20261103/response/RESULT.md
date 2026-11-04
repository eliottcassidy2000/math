# Full twelve-vertex response: a parity hostile, not a collision

**FINITE-EXACT.** The specified twelve-vertex tournament has rows

    3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002.

Bit j of row i means i->j. Rows are retained explicitly: twelve vertices
have66 arc bits, so packing them into a single64-bit integer would lose
information. The second control reverses every arc while preserving names.
The first core has S=E modulo2 but28 unordered pairs with nonzero C skew
modulo4. It is one of the scout's near-self-converse controls, not a survivor
of both necessary conditions and not a reconstruction counterexample.

## Full data retained

For each of the two cores the artifact keeps every entry of:

- S and E: named starts and ends of Hamilton paths;
- M: joint start/end of a single Hamilton path;
- C: ordered two-path covers, with the first path ending at i and the
  second starting at j, whether or not i->j is an actual arc;
- U and V: both three-index two-path endpoint tensors;
- D_(i,j,k,l): three ordered nonempty paths, first ending at i, middle
  starting at j and ending at k, last starting at l. The middle singleton
  j=k is included. All12^4=20736 entries are saved for each core.

One implementation enumerates all12! permutations and all admissible
one- and two-cut positions. The other builds a squarefree path matrix and
multiplies over disjoint support partitions. All entries agree exactly.
Unsigned64-bit arithmetic is sufficient: the generous three-cover bound
55*12! and actual fourteen-vertex Hamilton bound14! are below2^64.

The converse checks are entrywise, not just equality of totals:

    S(T^op)=E(T),  M(T^op)=M(T)^t,  C(T^op)=C(T)^t,
    U_(i,j,k)(T)=V_(k,j,i)(T^op),
    D_(i,j,k,l)(T)=D_(l,k,j,i)(T^op).

They follow by reversing all three paths and their order. These identities
preserve the named endpoint data rather than treating separate marginal
histograms as the observer.

## Complete response comparison for this pair

For ports x->y let a_i mean i->x and c_i mean i->y. The full response has
bidegree at most(2,2), so masks of weight at most2 in each port form an
interpolation battery of(1+12+66)^2=6241 actual queries. The retained tensor
formula is evaluated on every battery entry and compared with independent
Hamilton DP on the actual fourteen-vertex extension. Nine further high-weight
queries per control check interpolation. There are12500 direct-query checks
in total. The complete two-port response on all2^24 mask pairs is determined
by each battery; it is not identified with its smaller H,S,E,C sidecar.

The two batteries differ in6240 of6241 positions. The first distinguishing
query has a=0 and c=1:

    H(T plus the same two ports)       =  627597,
    H(T^op plus the same two ports)    = 1139039.

Here x beats every old vertex and y, while only old vertex0 beats y.
In fact x is a universal source, so this witness is already a one-port
distinction after deleting x. It is recorded honestly as a hostile to
endpoint-parity sufficiency, not as a difficult full-response collision.

The universal formula also obeys
F_(T^op)(a,c)=F_T(1-c,1-a), by reversing the whole enlarged tournament and
exchanging its two port names. Equality at the **same** masks, required for
a named collision, is a different assertion and fails in this example.

## Reproduction and limits

Run `python3 run.py` and `python3 -O run.py` here. Both retain all arrays
and actual query values, with matching frozen receipts and source hashes.
The selected universe is exactly one core and its converse. The scout
contains more candidates but no complete order12 class enumeration is
claimed. A future joint-modular survivor still needs its full response
comparison; a necessary-condition survivor alone is not a collision.
