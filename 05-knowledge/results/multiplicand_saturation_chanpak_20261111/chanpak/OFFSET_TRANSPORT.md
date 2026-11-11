# Physical offset transport and a two-parameter growth grid

**PROVED:** the three marked carriers below, their all-module substitution
formulas, and their commuting physical growth operations. **FINITE-EXACT:**
the stated order-seven discovery menu. **OPEN:** a uniform supply of such
carriers for arbitrary odd coefficients, and Chan–Pak. The complete
coefficient-35 grid is already proved in [PROOF.md](PROOF.md); this note
explains a further physical transport mechanism rather than extending that
coverage claim.

## 1. Retaining the offset makes transport possible

For a marked tournament (B,O), retain the entire state

    D=H(B)−W_B(O),     U=U_B(O).

Here U counts ordered spanning two-path covers whose second path is nonempty
and starts in O, permitting the first path to be empty. W requires both
paths to be nonempty, with first end and second start in O. A nonempty
transitive tail of length t and its final apex give an actual tournament with

    h=D+2^t U.                                            (1)

The arcs and the proof interface are those in [PROOF.md](PROOF.md), section2.
The old unit-offset requirement D=1 is unnecessary for (1).

The new question is physical: can a graph operation carry

    (D,U) → (2D−1,2U),                                   (2)

while retaining the selected set and a reusable port? The following three
carriers answer yes. All vertex indices are zero-based and rows are outgoing
neighbour masks. Replace the specified unselected vertex by TT_m, with
identical exterior neighbourhood, and retain the indicated original marks.

| name | base rows | O | port | H_m | W_m | U_m |
|---|---|---|---|---|---|---|
| 3 | (40,33,3,6,15,28) | {3} | 4 | 1+16·2^m | 0 | 24·2^m |
| 21 | (42,28,57,48,33,2) | {3,4} | 2 | 15·2^m−1 | 11·2^m−2 | 42·2^m |
| 35 | (86,124,40,113,100,65,4) | {6} | 1 | 1+32·2^m | 0 | 35·2^m |

Every displayed identity holds for **every m≥1**. The offsets are,
respectively,

    1+16·2^m,     1+4·2^m,     1+32·2^m.                 (3)

Thus adding one transitive twin inside the same module proves (2), for
every positive module length. No selected vertex is cloned or relocated.
The internal transitive module remains available after each operation.

The bases are strong. Replacing a vertex of a strong tournament by any
nonempty tournament preserves strong connectivity, by lifting quotient
paths. The marked tail output is also strong: its last vertex reaches a
nonempty selected set and hence the whole base; the base reaches the tail
and then its last vertex.

## 2. A finite complete proof, rather than extrapolation

We use the ordered-run theorem proved in
[Session260, constructors, section1](../../session260_quasiperfect_chanpak_20261111/constructors/RESEARCH.md).
With k exterior vertices and at most p ordered paths, compress every maximal
run of module vertices to a symbol. At most k+p run places can occur.
Endpoint restrictions refer only to exterior marked vertices, so the same
finite skeletons work for H,U,W. A skeleton with j canonically ordered run
places has exactly

    F_j(m)=j! S(m,j)

expansions in TT_m. Consequently each observable is a linear combination
of F_1,...,F_(k+p) for positive m. Its coefficients are uniquely determined
by evaluations m=1,...,k+p, because this evaluation matrix is triangular
with diagonal m!. H uses p=1; using p=2 for all three quantities gives the
common cutoff n+1 when the original base has n vertices.

The complete nonzero coefficients, starting at j=1, are:

| carrier | H coefficients | U coefficients | W coefficients | D coefficients |
|---|---|---|---|---|
| 3 | (33,16) | (48,24) | (0,0) | (33,16) |
| 21 | (29,15) | (84,42) | (20,11) | (9,4) |
| 35 | (65,32) | (70,35) | (0,0) | (65,32) |

Every remaining coefficient through the proved cutoff is zero: through7
for the first two carriers and through8 for the third. These integers are
computed directly from subset endpoint path counts on the literal graphs;
all evaluations and all coefficients are retained in
`offset_ladders_receipt.json`. Since F_1(m)=1 and F_2(m)=2^m−2 for m≥1,
the displayed formulas follow. The finite skeleton bound makes the
calculation an exact all-length certificate.

Deletion of the module, m=0, is outside these formulas. No absent-module
value or negative-length predecessor is used to extend their coverage.

## 3. Two commuting graph operations

After attaching the tail, keep the state (D,U,h), where h=D+2^tU. The
internal clone operation L and the tail-growth operation R are

    L(D,U,h)=(2D−1,2U,2h−1),
    R(D,U,h)=(D,U,2h−D).                                (4)

These are actual additions of a vertex in two disjoint transitive modules.
They commute, both physically and in the displayed state:

    L R(h)=R L(h)=4h−2D−1.

The offset D is the necessary sidecar. Dropping it would make R look like
an unsupported scalar doubling rule. Keeping it exposes a genuine growth
grid whose two directions may be taken in either order.

For N(d,r,s)=1+2^s+d·2^(r+s), r≥1, the parameter choices are:

| d | covered s | internal length m | tail length t |
|---|---|---|---|
| 3 | s≥5 | s−4 | r+1 |
| 21 | s≥3 | s−2 | r+1 |
| 35 | s≥6 | s−5 | r+5 |

The full3 and21 grids, and the35 quadrant s≥6, were already known by other
constructions. These are additional physical presentations with an exact
offset-advance operation. In particular this does not erase the finite
heads below each starting level. The new35 completion uses the separate
offset9,17,33 sources and product replacement in [PROOF.md](PROOF.md).

## 4. The stronger arbitrary-module interface

The same certificate is more general than a transitive-chain identity.
Let M be **any nonempty tournament**, let C_j(M) count ordered nonempty
spanning j-path covers, and put

    Q_2(M)=C_2(M)+2H(M).

This Q_2 counts ordered spanning two-path covers allowing either slot to
be empty. Replace the indicated port by M, leaving every mark exterior.
The very same skeleton coefficient multiplies C_j(M), as in
[Session261's marked module replacement theorem](../../session261_creative_chanpak_20261111/constructors/MODULE_REBASE.md).
The preceding table therefore proves the exact interfaces

    carrier3:  D=H(M)+16Q_2(M),   U=24Q_2(M),   W=0;
    carrier21: D=H(M)+4Q_2(M),    U=42Q_2(M),
               H=15Q_2(M)−H(M), W=11Q_2(M)−2H(M);
    carrier35: D=H(M)+32Q_2(M),   U=35Q_2(M),   W=0.       (5)

After the nonempty tail, these construct actual Hamilton counts

    H(M)+(a+b·2^t)Q_2(M),  t≥1,                         (6)

for (a,b)=(16,24),(4,42),(32,35), respectively. This is an exact additive
and multiplicative interface, with the two summands coupled through the
actual path-cover statistics of M. They are not independently selectable
integers. Taking M=TT_m gives H(M)=1 and Q_2(M)=2^m, recovering sections1–3.

Any module replacement preserving H and Q_2 preserves (D,U), hence all
future tail counts, in these three carriers. Their higher-run coefficients
vanish, so equality of every C_j is sufficient but unnecessary here.
Different internal root menus may still be retained as an additional
sidecar, just as in Session261. Conversely, composing arbitrary instances
of (6) requires the output's own Q_2; its Hamilton count alone is not a
closed continuation state. No universal interval saturation is claimed.

## 5. Scope, controls, and reproduction

The direct newborn menu on the three new35 offset sources found no cut
implementing (2) in one arbitrary birth while retaining the old marks
(with the newborn optionally selected). This finite failure is in
`offset_transport_receipt.json`; it motivates changing the source carrier
and does not contradict the proved repeated-clone constructions above.

The positive order-seven discovery menu used all456 ordinary tournament
representatives inherited from the literal n≤7 bank, every unselected port,
and every nonempty exterior selected set. It retained D>1 with D−1 a power
of2 and tested (2) at module lengths1,2,3. There are nine presentations:
four with odd part21 in U, four with odd part3, and one with odd part35.
This menu alone gives candidate identities, not all-length proofs. The
explicit carriers above were then proved by the complete finite coefficient
certificate. The order-six3 and21 carriers give smaller presentations than
the corresponding order-seven candidates.

Run the complete proof certificate:

    python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/offset_ladders.py

The optional discovery menu is reproduced separately by adding `--search7`.
It reads the bank at the inherited revision f92e070445, includes converses
when the bank records an ordinary multiplicity of2, and writes
`offset_ladders_search7.json`. It does not classify all larger carriers.

The remaining constructive gate is precise: obtain enough physical
offset/capacity carriers, with their actual starting thresholds, to cover
all residual odd coefficients and all below-threshold heads. The three
carriers prove that offset transport itself is possible; they do not yet
supply that uniform bank.
