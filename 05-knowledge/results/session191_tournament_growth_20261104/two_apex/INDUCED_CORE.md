# The H=13 core imposes an exact birth cost and forbids strong multiplication by three

**PROVED with complete finite certificates and a cited all-order size cutoff.**
Let C be any strong tournament with H(C)=13.

1. Every one-vertex induced extension C+x has count in exactly

       {13,23,25,29,31,33,37,43,45}.                    (1)

   In particular neither 39=3H(C) nor 41=3H(C)+2 can be obtained with one
   birth, even with arbitrary incidences and no two-apex restriction.
2. Exactly two births are necessary and sufficient to obtain a strong
   H=41 induced extension of C. They can pass through a strong H=23 state.
3. **No strong tournament with H=39 contains any induced subtournament with
   H=13**, at any order. Strong H=39 tournaments nevertheless exist.

These statements concern induced containment: old vertices and old arcs are
retained. They do not exclude changing old arcs, choosing a new parent, or
using nonstrong products. In particular C->C_3 has H=39 and contains C, but
is nonstrong.

## 1. All-order reduction to one actual core

Busch's bound H(T)>=5^((n-1)/3) for a strong n-tournament gives

    H=13 implies n<=5, since 5^5 > 13^3;
    H=39 implies n<=7, since 5^7 > 39^3.                (2)

There is no H=13 tournament on fewer than five names; this can be checked
by the complete order-five census, because source padding would carry any
smaller counterexample to order five. Equivalently the small standard path
census gives the order-at-most-four support {1,3,5}.

Among all 1024 labelled five-vertex tournaments, exactly 120 have H=13.
They are precisely the 120 distinct relabelings of the explicit core

    C = (6,8,26,17,3),                                (3)

where rows are outgoing bitmasks on names 0,...,4. Thus C is the unique
isomorphism class and has trivial automorphism group. These are complete
finite statements, not an assumed uniqueness classification.

`verify_extension.py` checks each of all 1024 path counts by both full
vertex-subset DP and literal permutation enumeration. It independently
constructs all 120 relabelings of (3), checks their distinctness, and compares
them with the complete H=13 list. This establishes (3) without relying on an
external isomorphism bank.

## 2. Every one- and two-birth extension

Adjoin name 5 with outgoing mask a on the five old names, for every a=0,...,31.
These are all possible one-vertex extensions of the fixed physical core.
Their counts, independently checked by subset DP and permutations, give (1).

Next adjoin name 6 with outgoing mask b on names 0,...,5, for every b=0,...,63.
The resulting 32*64=2048 matrices exhaust all two-vertex extensions. Their
complete count histogram and the 39/41/43 target records are generated in
`extension_receipt.json`. **None has H=39.**

The second independent counting algorithm enumerates the 7! possible
Hamilton words. A word inconsistent with any fixed old-core arc is rejected.
Every other word fixes at most four of the eleven new-arc bits and contributes
one to every compatible completion of the remaining bits. This produces all
2048 path counts simultaneously, without subset DP or the two-apex transfer.
The verifier compares every entry with its independent subset-DP count.

Any strong H=39 tournament containing C would have at most seven vertices
by (2), so it would appear among C itself, its one-birth extensions, or these
two-birth extensions. None does. This proves the strong containment exclusion.
If an induced subtournament U has H(U)=13 but is nonstrong, SCC multiplicativity
and primeness of 13 give a strong induced component C with H(C)=13. Therefore
the exclusion also applies to arbitrary induced H=13 subtournaments, not just
strong ones.

## 3. Positive controls and the exact two-birth cost

Use the two-apex grammar of `TWO_APEX.md`. The literal chain

    120 -> 1201 -> 12011

has counts

    13 -> 23 -> 41.

Every state is strong, and each arrow adds exactly one vertex without changing
an old incidence. The explicit outgoing rows are

    C_0 = (6,8,26,17,3),
    C_1 = (38,8,58,49,35,2),
    C_2 = (102,8,122,113,99,66,2).

Thus two births suffice. The complete one-birth support (1) shows that one
cannot suffice, proving the exact cost. Relabeling transports the extension
to any strong H=13 parent.

The word 10222 is a strong seven-vertex tournament with H=39, independently
checked by subset DP and literal permutations. Therefore the induced-core
exclusion is not vacuous and is not a global missing spectrum value.

Finally the upper child 43=3*13+4 is already attained by one birth. For example
add name 5 with outgoing mask 5 to (3). This distinguishes the two arithmetic
branches: the lower child 41 needs two births, whereas the upper child 43
needs one. It does not assert a general branch-cost law for other H.

## 4. Reproduction and source scope

    python3 verify_extension.py
    python3 -O verify_extension.py

Both runs produce identical stdout and deterministic JSON, using only the
standard Python library. All checks use explicit exceptions.

The all-order inference uses only the complete finite extension universe,
SCC multiplicativity, and the cited Busch bound already used in this packet
(Arthur H. Busch, EJC 13 (2006), N3, https://doi.org/10.37236/1141).
No new independent verification of the source paper's full proof or general
literature-priority claim is made.
