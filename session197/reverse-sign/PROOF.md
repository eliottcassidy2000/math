# Two directed sign bridges on one marked ternary clock

## Precise result and limits

The same marked root can pass through three physical clone families with counts

`60*3^m-1 -> 120*3^m+1 -> 420*3^m-1`, for every m>=1.

Each arrow adds exactly one vertex. Five current ports suffice throughout; no old vertex is re-exposed. The root-cloning operation remains available at all three stages. Its returning scalar laws are respectively `3H+2`, `3H-2`, and `3H+2`.

The first arrow is `2H+3`. The second is `(7H-9)/2` on the middle family. With the same m at both ends, their composition is `7H+6`.

These statements concern the three specified physical families. The two-bridge word is not closed under repetition, and its scalar relation does not hold on all physical states. Explicit counterexamples appear below. This packet is separate from the fixed-five-old-port positive binary–ternary family: the first birth here makes a different forget choice to retain its newborn as a port.

## Core graphs and root

The root is old vertex7 in all three core arrays:

`T=[2,12,145,149,99,143,175,19]`;

`P=[2,12,145,149,99,399,431,19,159]`;

`N=[2,12,145,149,611,911,431,19,671,207]`.

P is obtained from T by adding A=8, which loses only to5,6 and beats all other old vertices. N is obtained from P by adding B=9, which loses only to4,5,A and beats every other old vertex. Thus A beats B; the original root loses to both.

For each core, replace root7 by a transitive module of m vertices. Every added gate has uniform incidence to that entire root module, so taking a root blowup and adding either gate commute as actual graph operations. This proves that the same two one-vertex extensions exist at every m, not only at m=1.

All three cores are strong. T contains the directed Hamilton cycle `(0,1,2,4,6,5,3,7,0)`. A has both an inneighbor and an outneighbor in T; B has both in P. A mixed one-vertex extension of a strong tournament is strong. Nonempty transitive blowups preserve strongness, so all three families are strong.

## Exact five-port protocol

Start with ordered ports `[0,1,4,5,6]`.

1. First gate: `L(mask24,drop4)`. The newborn A loses to port positions3,4, namely old5,6. Forget old6 as a port and retain A. The new list is `[0,1,4,5,A]`.
2. Second gate: `L(mask28,drop5)`. The newborn B loses to positions2,3,4, namely old4,5,A. Forget only B's temporary exposure. The list stays `[0,1,4,5,A]`.

Both steps retain five distinct ports. Vertex6 remains in the graph after it is forgotten, but is never re-exposed. A stays exposed through the second gate, and B is never retained as a port.

At every stage use `R(mask7,drop5)` to add a root clone. The newborn beats precisely the first three ports, which are always old0,1,4, and loses to every other old vertex. These are exactly the original root's incidences, including its losses to A and B when present. The same literal operation can therefore be repeated any number of times before, between, or after the sign gates.

If any root births occur, A and B need not have physical labels8 and9. The masks refer to current port positions, not guessed vertex labels. In the canonical order with all root births first, A has labelm+7 and B has labelm+8. The correspondence fixes old0,...,6, maps the m root vertices to original7 followed by root newborns in birth order, and maps A and B to their actual birth labels. Every incidence is preserved by this correspondence.

The two gates are used at most once in this protocol. No claim is made that the final marked graph is a new admissible starting state for the same two-gate law.

## All-order clone identities

For a marked root, let a_j count literal root skeletons with j root occurrences and no consecutive roots. Replacing the root by m transitive clones gives the exact identity

`H_m=sum_j a_j j! S(m,j)`.

A skeleton's j-1 internal fragments yield vertex-disjoint crossing edges from the root's outneighbor side to its inneighbor side. In each of T,P,N the fixed two-vertex cover `{1,4}` meets every such crossing edge. In T these edges are `1->2,1->3,4->5,4->6`. P adds none; N adds only `4->B`. Therefore j<=3 uniformly at every root multiplicity.

The verifier counts all finite core skeletons by a separate subset recurrence, giving

- T: `(a0,a1,a2,a3)=(59,179,180,60)`;
- P: `(a0,a1,a2,a3)=(121,361,360,120)`;
- N: `(a0,a1,a2,a3)=(419,1259,1260,420)`.

All higher entries vanish. The count a0 is the root-deleted Hamilton count and is not a term in the positive-clone formula. Since

`1!S(m,1)=1`,

`2!S(m,2)=2^m-2`,

`3!S(m,3)=3^m-3*2^m+3`,

the three run vectors prove respectively `60*3^m-1`, `120*3^m+1`, and `420*3^m-1` at every m>=1. This is a finite exact coefficient certificate in a proved three-mode span, not a fit extrapolated from small m.

It follows that each stage's R birth returns the stated ternary affine law, and that the sign gates have the stated scalar identities on their respective families. If r root births are inserted between the gates, the composite sends a boundary count H to

`7*3^r*(H+1)-1`.

Thus `7H+6` is the r=0 case, or the comparison between the two cores at equal root multiplicity.

## Exact scope controls

### Repeating the same five-port word fails

Let W be the exact position-word

`L(24,drop4), L(28,drop5)`.

Starting from T with `[0,1,4,5,6]`, one application gives H=1259. Applying W a second time is a legal five-port word, but its count is11071 rather than `7*1259+6=8819`. The retained ports are then `[0,1,4,5,10]` in the consecutive-birth labeling.

With the original root still marked, that second output has run vector `(3607,11071,11364,3900)`, and hence clone count

`3900*3^m-336*2^m+43`.

The unwanted second mode is physically present. Consequently these two arrows do not by themselves form an iterable sign-switching grammar.

### The composite is not a universal physical-state identity

Start instead from the transitive five-vertex tournament, with all its vertices as ports in order. The identical word W produces Hamilton count121, whereas `7*1+6=13`. Its output rows are

`[30,28,88,112,32,71,19]`.

Thus the relation `7H+6` cannot be promoted from the constructed ray to arbitrary actual five-port states.

### What is shared and what is not

The root and its three essential clone ports0,1,4 are shared physically. The five-port protocol explicitly tracks its other ports, and the root operation is closed within each of the three stages. The phase of the core determines which affine offset applies. There is no same-state choice between the two offsets, no old-port re-exposure, and no assertion that a scalar equality identifies two marked graphs.

## Reproduction

Run `python3 verify.py` and `python3 -O verify.py`; their deterministic receipts agree. The verifier reconstructs both gate arrays from literal births, checks all full root-run vectors, independently counts source Hamilton words by depth-first traversal, tests every distribution of up to three root births before/between/after the gates, checks graph relabeling and strongness, and reproduces both scope counterexamples. The positive fixed-old-port binary–ternary family is neither imported nor changed by this appendix.
