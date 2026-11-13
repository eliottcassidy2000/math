# Target selection with shared owners and reusable module repairs

**PROVED IN NOTES / FINITE-EXACT; universal Chan–Pak and ordinary Collatz OPEN.**
2026-11-13. Starting snapshot `ed79a99437c3bf19509b2f9b7b696ab9349032d3`.
This packet answers the requested connection between the diagonal Fibonacci
coloring, oriented ancestry charts, and target-dependent tournament construction.
Its strongest constructive result is a literal handoff between changing
module frames, followed by a repair that works for every nonempty module.

## The physical advance

Write h=H(M), p=Q_M(2) for any nonempty tournament M. The
[chirality preparation](CHIRALITY.md) and [named handoff](HANDOFF.md) give

    3p−h → 9p−h → 9p−h → 15p−h.

These are actual graphs: two root-child births, a count-neutral module
relocation preserving the original root star, then a named helper birth.
The relocation changes the root response by
`(0,−(p−2h),−(p−2h),−2(p−2h))`; scalar equality alone was insufficient.
The original root and internal module survive, and the endpoint is in the
reusable repair frame. Examples are `11→35→35→59` and
`289→889→889→1489`, where 889=7·127 and 1489 is prime.
The relocation uses retained frame labels beyond the two-port language.

The [arbitrary-module repair](selection/PROOF.md) changes one old arc and
adds one vertex, connecting two strong tournaments with exactly the same
rooted state

    (H,s,q,F)=(15p−h,5h,5h,15p+5h).

The old root star is retained. A common fresh seal then gives equal
responses to every continuation in the declared six-command two-port
language. Complete two-run coefficient identities prove this for all M.
The repair commutes with an actual binary growth step, M↦TT1→M,
which induces `H↦2H+h`. This is a commuting graph diagram, not a fitted
scalar recurrence. It is not indefinite neutral pure growth at fixed H.

For the explicit module D with `(H,Q(2))=(11,100)`, ordinal products give
an unbounded family of primitive repair frames

    H(n,m)=15·100^n·2^m−11^n,   n,m≥0, n+m>0.

The [residue compiler](RESIDUE_SELECTOR.md) selects every unit modulo11^k
for each fixed n≥1 and can also impose `v_7(H)=1`. Its choices are
actual modules with the reusable ports. Formulae identify prime and7p
fixtures, without claiming infinitely many prime cofactors. All members
are4 modulo5, so this family cannot be the universal selector.

[Finite prefix/suffix control](PREFIX_SUFFIX.md) goes further: every
finite leading base11 word and unit trailing word can be imposed together,
even with exactly one factor7. The proof combines the retained-port
residue clock with an elementary irrational rotation. Arbitrarily many
middle digits remain unspecified; the O(log X) family is still sparse.

## What the older constructions contribute

| Recovered construction | Exact additional result | Boundary retained |
|---|---|---|
| Diagonal grid C→CA, A→B, B→CA | Every globally shared proper-expansion dictionary has a unique coherent stripe; visible choices have an exact owner map | Repeated occurrences of an owner cannot be selected independently |
| Additive pair counts in that stripe | A Boolean target polynomial in O(log N) owner bits, with reconstructible witnesses | At N=16 the image is{0,2}, not{0,1,2}; at N=18 it is{1,2,4,5} despite coefficient gcd1 |
| Signed ternary ancestry chart | Neutral cyclic reversal changes root starts by a weighted sum of response determinants | The chart's numeric minors and the response determinants are different objects |
| Colored hierarchy and marking | Weighted maximum-spanning-tree formula; exact adjacency/quota/residual-group rank tables; joint shared-owner selection | Greedy alone does not handle coupled quotas, and rank is not H |
| Physical marking | Accepted marks compile into root-child births with the requested actual global stabilizer | Symmetry breaking does not preserve H or imply unit readiness |
| Rooted returns | Arbitrary-module one-arc repair and the literal inter-frame handoff above | The pair(h,p) must still come from an actual source module |

The [history route](history/HISTORY.md) follows the requested grid and
neighboring sessions, including forced-zero elimination, assignment loads,
odd reflected Collatz branches, and boundary completion. The
[shared-dictionary proof](history/SHARED_DICTIONARY.md) resolves a
previously unclassified quantifier, rather than treating every occurrence
as a fresh bit. One legal occurrence-specific release repairs the first
pair-count gap; applying this idea to tournaments requires a physical
release with an exact response, as the module relocation now supplies in
one concrete setting.

The [chart selector](charts/PROOF.md) keeps ancestry, cyclic orientation,
physical palette, eligibility, and external role colors separate. Valuation
alone admits six sibling permutations; oriented leading digits retain
only rotations; actual child classes decide which rotations survive.
The weighted greedy theorem was independently recovered concurrently in
S313 and is credited there. Our extra interfaces are the decorated MST,
coupled constraints, shared-owner selection, and the stated birth protocol.

## The next precise proof obligation

Supply an operation on actual framed modules that escapes the current
absolute image while retaining a usable continuation interface. For a
target Y in the current repair frame, one needs `15p−h=Y` with
`(h,p)=(H(M),Q_M(2))`, not independently chosen odd/even integers.
The chirality move can change marked response while preserving H; the
handoff can change module role while preserving the root star; neither
yet makes all such pairs realizable.

Three tempting shortcuts have explicit hostile controls: independent
choice at repeated dictionary positions; replacing attainable target
sets by their min/max; and inferring exact integer coverage from dense
modular or rank images. The new construction removes a physical handoff
gap. It does not remove the joint-source selection gap, and it is not a
proof that every odd count other than7 and21 occurs.

## Incoming connections and reproducibility

[INCOMING](INCOMING.md) credits the concurrent actual-count charts,
additive ports, modular bypass, and scalar stabilization. In particular,
Q(2), the deletion response q, and a marked apex response are not
interchangeable merely because different notes use the same letter.
[AUDIT](AUDIT.md) records quantifiers, independent checks, hostile
examples, and what was not inferred. [VALIDATION.json](VALIDATION.json)
and the manifest record the exact local replay. There is no new theorem
ID reservation or universal coverage promotion.

Run the ten packet checkers from the repository root, as listed in
VALIDATION.json. The proof notes justify unbounded identities; the
receipts state their finite universes. No legacy large census was rerun.
