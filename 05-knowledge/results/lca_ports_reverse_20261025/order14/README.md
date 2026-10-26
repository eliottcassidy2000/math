# One relocation rectangle repairs the four order-fourteen exceptions

2026-10-25. **PROVED:** every completion in each of the four order14
failed source cells from
[LOW_ORDERS.md](../../port_residual_integral_20261025/LOW_ORDERS.md)
has three successive directed black K moves to its exact labelled target.
Use the inherited ledger unless its last word is blue; in that case
the single new ledger below repairs all four exceptions.

The audited [one-bit extension](GUARD_EXTENSION.md) strengthens this:
read only T(3,5), using the new ledger when it is0 and the inherited
ledger when it is1. The new ledger works on the entire zero-guard
half-cell, not just its unique symmetric exception.

The exceptional symmetry itself is structurally classified: each
28-dimensional source cell contains exactly one failed source.
**FINITE-EXACT:** the independent affine construction and literal
replays in check.py verify this classification and the inherited first
two words over the full cells. Full black-core connectivity is not claimed.

## 1. Inheritance, exact universe and concept board

The closest mechanism is the common-final symmetry argument in
[endpoint_three_repairs/theory.md section2](../../endpoint_three_repairs_20261020h/theory.md).
Its two equal cycles at the even critical order are retained here as an
actual automorphism, rather than used only for a degree obstruction.
The exact failed cells and the complete source-universe audit are in
[LOW_ORDERS_REVIEW.md](../../port_residual_integral_20261025/LOW_ORDERS_REVIEW.md).
The canonical hostile is a genuinely blue final word, not an illegal
path. The corrected near miss is interpreting such a ledger failure
as a disconnected black component. The least-used sidecar is the action
of the failed word's symmetry on a complete endpoint row.

The concept board is two seven-cycles / exact incoming ports / endpoint
row reconstruction / shared-vertex relocation squares / degree sums /
literal clone preservation. The small probe first tried the24 distinct
inherited base ledgers, which do not repair these four points. A menu
of at most three contiguous P blocks then found the common rectangle.
The proof below does not depend on that search or enumerate Hamiltonian paths.

Let T have vertex set0,...,13, with the directed closed blue word

    P=(0,1,...,13),       r_P(v)=13-v.

Closed means13->0. Blue means its actual positional reversal is anti:

    T(u,v)=T(r_P(v),r_P(u)).

Assume endpoint outdegrees(9,4), a closed middle12->1, and
N^-_T(0)={13} union I. The four assigned cells are

    I={3,5,10} or{3,7,10}, with2->12;
    I={5,9,10} or{7,9,10}, with12->2.

All other arcs remain arbitrary subject to those conditions. There is
no strongness, module, clone, overlap or isomorphism filter.
A K_Q move retains the arcs of the directed word Q and complements
every other unordered pair. Put B=K_P(T).

The inherited three-word ledgers have common last word

    F=(0,13,1,2,...,12).

For2->12 they are

    (12,1,13,0,2,...,11),
    (13,12,1,0,2,...,11),
    F.

For12->2 they are

    (12,2,...,11,1,13,0),
    (13,12,2,...,11,1,0),
    F.

Their first two words are black throughout all four affine cells;
their legality and literal target are already certified in the linked
low-order audit and rechecked here. F has the same blue status at both
ends of its own K edge. If F is black in B, the inherited ledger
therefore completes the task. Everything that follows treats F blue in B.

## 2. The failed symmetry determines exactly four points

P remains blue in B. The composition of its anti-reversal with that
of F is the automorphism

    g=r_F r_P
      =(0,11,9,7,5,3,1)(2,13,12,10,8,6,4).             (1)

Write A for the first orbit and C for the second. Each has seven
vertices and is transitive under g. In particular B[A] and B[C]
are regular tournaments of degree3. Also r_P exchanges A and C.

The known endpoint row is

    N^+_B(0)={1,13} union I,     d_B(0)=5.             (2)

Since the middle is closed, B(1,12)=1. Applying g gives

    B(0,10)=1,

so10 belongs to I. Regularity of B[A] then forces exactly two members
of I in A. Hence the other two ports are odd, and4,12 are not ports.
Inside the cyclic orbit A, opposite powers of g give opposite arc
directions from0. Its paired choices are

    {11,1}, {9,3}, {7,5}.

The first choice is1 because0->1 is a retained path arc. Therefore

    I={10,t,s},        t in{3,9}, s in{5,7},          (3)

which is exactly the four displayed port sets; in particular11 is
not a port.

The x,c chord is forced too. P-anti and g give

    B(2,12)=B(1,11)=B(0,9).

These pairs are nonpath pairs, so T(2,12)=1 exactly when3, rather
than9, belongs to I. Thus (3) reproduces the complete four-cell list,
including the branch bit.

For each I in(3), the entire B is uniquely determined. The row at0
and powers of g determine every ordered pair with its first vertex
in A. Reverse orientation determines pairs from C to A, and P-anti
maps any pair entirely in C to a pair entirely in A. Hence there
is at most one B, and thus at most one T, for each I. The explicit
carrier construction in check.py produces a consistent tournament
for every one of the four choices, proving existence as well.

This is the source/target map used in the repair: the literal failed
tournament maps to its two-cycle action plus its full row at0.
That carrier preserves every arc because P-anti restores the second
orbit. Retaining just the two cycle sizes would destroy the row data
and would not determine a tournament. No isomorphism contraction occurs.

Two consequences of the symmetry are all the new ledger needs:

    T(3,5)=0;       d_B(v)=8 for every v in C.        (4)

For the first, g^2 sends(3,5) to(0,1). Thus B(3,5)=B(0,1)=1,
and K_P complements the nonpath pair{3,5}. For the second,
P-anti gives d_B(13)=13-d_B(0)=8; transitivity gives that degree
at every vertex of C. In particular d_T(4)=d_T(12)=5:
these are interior vertices of P, so K_P sends degree d to13-d.

## 3. One common three-move repair

Use

    Q1=(0,1,2,13,3,4,5,6,7,8,9,10,11,12),
    Q2=(4,0,1,2,13,3,5,6,7,8,9,10,11,12),
    Q3=(4,0,1,2,3,5,6,7,8,9,10,11,12,13).          (5)

Legality is local. Q1 introduces2->13 and13->3. P-anti identifies
the former with0->11, true because11 is not a port; it identifies
the latter with10->0, true because10 is a port. Every other Q1
edge is an original P edge.

At S=K_Q1(T), Q2 introduces4->0 and3->5. Both pairs were omitted
by Q1. Their original directions were0->4 and5->3, respectively,
by(3) and(4), so both new joins are true. All other edges are retained.

Q3 restores2->3 and12->13, each omitted by both Q1 and Q2 and
therefore complemented twice. Its other edges are retained from Q2.
Every word contains every literal vertex once.

For the exact target, define unordered edge sets

    X={{2,3},{12,13},{2,13},{3,13}},
    Y={{3,4},{4,5},{0,4},{3,5}}.

Then

    E(Q1)=E(P) xor X,
    E(Q2)=E(P) xor X xor Y,
    E(Q3)=E(P) xor Y.                                (6)

The relocations share vertex3; they do not require disjoint port sets.
Their XOR is exactly E(P). Since three full complements leave one,
the three K moves end at the literal B=K_P(T). Every clone and its
adjacency remain present in this equality.

## 4. Blackness from three strict endpoint sums

If a directed word is blue in an n-vertex tournament, the sum of
the outdegrees of its endpoints is n-1. Here the required sum is13.

* Q1 has endpoints0,12 in T, with degrees9 and5. Their sum is14.
* Q2 has endpoints4,12 in S=K_Q1(T). Vertex4 is internal in Q1,
  so its degree becomes13-5=8. Vertex12 is last in Q1, so its
  degree becomes12-5=7. Their sum is15.
* Test Q3 at its other endpoint B. Its endpoints4,13 both belong
  to C, so their degrees are8 and8. Their sum is16.

Thus all three words are black. Testing the last word at B is
legitimate: its reversal preserves its own path-edge set, so K_Q3
preserves the word's blue/black status. Alternatively its degrees
at the actual pre-move source are6 and4, whose sum10 also fails13.

Combining this repair when F is blue with the inherited ledger
when F is black proves the full four-cell theorem. The choice of
ledger occurs at the original T; no prefix of one ledger is spliced
into the other.

## 5. A cheap hostile and the precise boundary

The new ledger is not asserted on every source in the four cells
without the F-blue condition. Start with the source for I={3,5,10}
and flip the two unordered pairs{3,5} and{8,10}. They are a full
P-reflection orbit, so the modified tournament stays in the same
source cell, with the same ports, path, closing arc and branch bits.
It now has T(3,5)=1. After Q1 its required Q2 join3->5 is false.
The receipt gives its exact rows and replays the inherited three-black
ledger successfully.

The first invalid implication after dropping symmetry would be
T(3,5)=0. The strongest survivor is the two-ledger rule, with the
actual F anti-reversal condition retained. A failed-word symmetry
can therefore serve as a source of forced joins, rather than being
only the obstruction encountered by an earlier proof.

## 6. Reproduction and consequence

Run from the repository root:

    python 05-knowledge/results/lca_ports_reverse_20261025/order14/check.py
    python -O 05-knowledge/results/lca_ports_reverse_20261025/order14/check.py

Both pass and produce identical normal.json and optimized.json.
The checker uses explicit exception gates and two paths:

1. Raw forward-arc Gaussian equations check every completion in the
   four original28-dimensional cells. They certify the old first
   two moves uniformly and that adding F-blue fixes all91 arc bits.
   A separate120-port-set pass, leaving the x,c chord unfixed,
   recovers exactly the four symmetry exceptions.
2. A direct g-orbit and P-anti construction reconstructs each
   tournament from its row0. Literal arc arithmetic checks the
   source predicates, all directed joins, blue flags, exact target,
   the degree certificates and the missing-symmetry hostile.

No Hamiltonian-path census or enumeration of2^28 completions occurs.
Together with the previously certified nonexceptional base cells
and residual I={h,11,12} repairs, this closes the closed endpoint-four
layer at order14. It remains a local marked-edge replacement result,
not a proof that every blue edge belongs to that layer or that the
entire black core is connected.

