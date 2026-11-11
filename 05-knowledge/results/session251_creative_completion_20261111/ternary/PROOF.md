# A fixed-core arithmetic compiler and the next deletion threshold

2026-11-11. **PROVED IN THIS NOTE**: the general compiler, terminating
inverse, and reduction to five explicit cores for the stated grammar.
**FINITE-EXACT**: new bounded-support experiments and graph controls.
**OPEN**: the all-height covering conjecture. Chan–Pak remains open.

## Constructive result first

A single positive construction has a much smaller description than an
arbitrary one-vertex branch: put an explicit core C between two transitive
chains, then add one vertex x with an arbitrary cut. The two chain cuts
are binary integers A,B. Nine fixed core coordinates give H exactly by
formula (1) below. Every inverse solution constructs the actual tournament.
The inverse terminates at every target, including uniform core contact.

The five-core bank is

    K1       (0,)                      H=1
    C3       (2,4,1)                    H=3
    S4       (2,4,9,3)                  H=5
    S5a      (2,4,9,19,7)               H=9
    S5b      (2,4,25,3,11)              H=9.

Rows are out-neighbour bit masks. All flank lengths are arbitrary, all
core cuts are allowed, and the final vertex is the new apex. The compiler
remembers the core and cut, so this is a tournament grammar rather than a
formal recurrence on scalar H. The bank fixes the source-deletion count
at 1,3,5, or9 independently of the target size.

The exact new small-parent boundary is **173**: every H173 tournament has
every one-vertex-deletion count at least9, and9 is attained. This extends
the inherited H55 boundary constructively: H55 has a small H5-parent
certificate, while H173 has a small H9-parent certificate.

## Inheritance, live concepts, and scope

Closest mechanism: the H3-source two-flank compiler in
[decomposition16 attachment](../../decomposition16_odd_spectrum_20261110/attachment/PROOF.md).
Its corrected uniform-contact formula retains module reentry. The other
closest mechanism is Session250's seven-coordinate rooted ordinal gluing
law in [constructive/PROOF.md](../../session250_sixteen_20261110/constructive/PROOF.md).
The underused sidecar here is the polynomial of directed paths weighted by
H of the unused core, retaining start/end ownership. This complements,
rather than replaces, the ordered two-cover coordinates in Session250.
The inherited odd-cycle formula is THM-002,
`01-canon/theorems/THM-002-ocf.md`; the inherited one-apex recurrence is
THM-4875(D),
`01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`.
The latter already reports a much larger finite coverage range for a richer
one-arc-join grammar. **No new range record is claimed here.**

Anchor: ternary coverage. Niche: weighted residual-core paths. Wildcard:
replace growing parent h by a constant-size parent-core bank. The board is
actual tournament / root cut / binary flank words / residual H weight /
finite target inverse / small-parent cover. The canonical hostile is
uniform module reentry; the scalar H of a source does not determine its
attachment language. The inherited nine missing ternary edges also warn
against claiming any fixed input tournament has all desired updates.

## 1. The core coordinates

Let C be a nonempty tournament, O the set of core vertices beaten by x,
and I=V(C) minus O. Put d=H(C), with H(empty)=1. For each nonempty directed
simple path P in C, give it weight H(C minus V(P)). Define polynomials

    K(z) = sum_P weight(P) z^|P|,
    I(z) = sum_(P ending in I) weight(P) z^|P|,
    O(z) = sum_(P starting in O) weight(P) z^|P|,
    J(z) = sum_(P starting in O, ending in I) weight(P) z^|P|.

The polynomial names I,O are distinguished from the underlying vertex
sets by their arguments. Write k_±=K(±1), i_±=I(±1), o_±=O(±1),
j_±=J(±1). The nine-coordinate signature is

    (d,k_+,k_-,i_+,i_-,o_+,o_-,j_+,j_-).

The positive coordinates are literal counts. The negative evaluations
retain parity, not a positivity claim. They are computed once per fixed
coloured core and then used at every flank length.

Let the old source be TT_a -> C -> TT_b. Read the left x-cut in dominance
order as binary integer A. Read the complement of the right x-cut in
reverse dominance order as B. Put α=A mod2, β=B mod2. Empty words give0;
leading zeroes may be retained to recover named padding. Let

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2)  (q>=1).

This is the inherited H count for one apex over a transitive chain.

**Theorem 1.** The actual child has

    H = d(g(A)+g(B)-1) + (d+k_+)AB + (d+k_-)αβ
        + i_+ A - i_- α + o_+ B - o_- β + j_+ + j_-.       (1)

No strongness assumption on C or mixed-contact assumption on O is needed.

**Proof.** In the odd-cycle formula, collections avoiding x contribute d.
A collection containing x has one odd directed cycle through x. Removing
x gives an even nonempty directed path in the old source, beginning in
its out-set and ending in its in-set. The residual odd-cycle collections
sum to H of the unused old vertices, and the cycle through x supplies a
factor2. Hence, if F(z) is the generating polynomial of these paths with
residual H weights, the extra contribution is F(1)+F(-1).

Let L(z) enumerate left-chain paths beginning in its out-set, and R(z)
right-chain paths ending in its in-set. Their boundary values are
L(1)=A, L(-1)=-α, R(1)=B, R(-1)=-β. Internal chain paths give
`d(g(A)+g(B)-2)` after the two evaluations. A path meeting no core vertex
but both flanks gives d L(z)R(z). A path meeting the core gives exactly

    L(z) K(z) R(z) + L(z) I(z) + O(z) R(z) + J(z).

A path in an ordinal join meets each part in one segment. Its unused core
has exactly the residual weight in these four polynomials. This is a
classification of the cycle path, not a claim that a Hamiltonian path
cannot revisit a module. Substituting the boundary values, adding d, and
collecting terms proves (1). QED.

For C3 and |O|=1 the signature is
`(3,9,-3,6,-2,3,-1,2,0)`, so (1) is exactly the inherited H3 formula.
The arbitrary-core statement retains all the terms needed for its
uniform-contact repair.

## 2. A terminating inverse, including uniform contact

Put L=d+k_+>0 and

    m = min_(α,β in {0,1}) ((d+k_-)αβ-i_-α-o_-β+j_++j_-).

Since g(A),g(B)>=1, (1) gives

    H >= d+m+L AB+i_+ A+o_+ B.                          (2)

If O and I are both nonempty, singleton paths imply i_+>0 and o_+>0.
For target bound N, put cap=N-d-m. A negative cap rules out every pair.
Otherwise every solution satisfies

    0 <= A <= floor(cap/i_+),
    0 <= B <= floor((cap-i_+ A)/(L A+o_+)).              (3)

Enumerate these finite integer pairs and retain exactly those satisfying
(1). Binary expansions then give the explicit labelled child.

If O is empty, then o_±=j_±=0 and i_±=k_±. The axis A=0 has
`H=d*g(B)`. Every g-value has an even representative2q: repeatedly strip
a trailing binary1 using the recurrence, and use0 for the terminal value1.
Moreover g(2q)>=q+1, so `0<=q<=floor(N/d)-1` suffices for this axis.
For A>=1 the same bound on A and B in (3) is valid, since L A>0.

If I is empty, the axis B=0 similarly has `H=d*g(A)`. For B>=1,
(2) implies `A<=floor(cap/L)`, after which the B bound in (3) is valid
with i_+=0 and o_+>0. Thus all cuts have a finite target inverse. There
is no arbitrary source-order cutoff and no assumed absence of reentry.

The inverse is complete for its explicit fixed-core grammar. Distinct
pairs or coloured cores can yield isomorphic graphs; no class-count or
uniqueness assertion is made. Leading-zero padding changes names and
order, but not the possible H values. A witness retains the core rows,
O,A,B and can reconstruct every arc.

## 3. The complete H<=5 parent language and the sharp H173 certificate

Here completeness of the **small-parent** classification uses the
inherited strong minimum, as in
[decomposition16 induction](../../decomposition16_odd_spectrum_20261110/induction/PROOF.md).
A source of H1 is transitive. A source of H3 is TT_a -> C3 -> TT_b.
A source of H5 is TT_a -> S4 -> TT_b: SCC multiplication forces one
H5 strong component, the strong minimum at order5 is9, and the unique
strong order4 class is S4. Transitive padding accounts for everything else.
These facts imply that the first three core templates, with all cuts and
all flanks, are exactly the union of one-birth images of all H1,H3,H5
sources. The inverse in §2 makes fixed-target nonmembership an all-order
statement, not a graph census truncated in size.

`compiler.py` enumerates all core classes through5 and all cuts, computes
signatures from literal simple paths, and inverts through4095. Restricting
the parent count to1,3,5 gives first missing values

    7,21,173,335,663,1383,2373,2457.

This is the complete missing list through4095 for this small-parent
language. In particular every allowed odd below173 has a child with
one-vertex deletion H<=5. The inherited global omissions7,21 and oddness
then give: **every vertex deletion of any H173 tournament has H>=9**.

A sharp constructive certificate is recorded in `explore-results.json`.
Its H173 child has rows

    (126,68,72,82,230,78,128,47)

and deletion counts

    (77,71,51,27,51,71,27,9).

The new apex is last, and deleting it leaves an H9 source. Direct subset
DP recounts the selected certificates. The final experiment also retains
explicit H9-parent certificates for the other five small-parent misses
through4095.

For comparison, one H55 certificate is S4 with O={0,1,2}, A=0, B=2:
its signature is `(5,27,-5,6,-2,18,-4,3,1)`, and (1) gives55.
This supplies the next positive parent lane rather than treating55 as a
global obstruction.

## 4. The concrete all-height surrogate

**OPEN fixed-core covering conjecture.** Every positive odd N except7,21
is obtained from one of the five displayed cores by some O,A,B in (1).
Equivalently, every permitted target has a tournament with a vertex whose
deletion is TT_a -> C -> TT_b for one of these five C.

A proof would prove Chan–Pak without establishing the stronger-in-another-
direction branches h -> 3h+2 and3h+4. It asks for constant-size source
cores rather than a prescribed source count. It is not a theorem about
extending every existing realiser. It can be attacked as a finite union
of explicit binary arithmetic images, using g and nine integer coefficients
per cut, with an actual graph constructor attached to every image.

The first three cores do not suffice, by173. Both H9 cores are genuine
new choices; they are not obtained from a shared scalar H9 by a guaranteed
attachment transport. The sixth diagnostic template C3 -> C3 is kept
separate in the experiment. Even including that template does not include
all H9 sources: C3 -> TT_r -> C3 with arbitrary internal gap r is a further
family. No such internal gap is silently collapsed into the exterior
flanks.

## Reproduction and loss ledger

    python3 -B 05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py
    python3 -B 05-knowledge/results/session251_creative_completion_20261111/ternary/explore.py --bound 262143

`compiler.py` constructs all1099 labelled tournaments of orders1..5,
quotients by all vertex permutations (converse is not merged), obtaining
class counts1,1,2,4,12. It tests formula (1) on all cuts of these20 cores
and A,B in0..3: **7520 actual child counts** by subset endpoint DP. This is
a finite implementation control of the new formula, not an independent
reverification of old spectrum tables. It then records complete arithmetic
support through4095 using the proved inverse. The companion explore script
uses the five explicit cores and the diagnostic sixth template, with no
imported large census. Saved selected certificates are actual row masks.

Source: fixed coloured core, two named chain words, apex. Target: H and
its bounded inverse. Preserved: exact H and reconstructible arcs. Lost on
projection to the nine coordinates: graph identity, named endpoint paths,
and all response to a subsequent arbitrary birth. Required sidecar: the
actual core rows and selected set, always stored. Identical signatures
may be merged only for current H-value support, not for future attachment
semantics. The compiler does not itself prove finite-state coverage or a
ternary branch rule; those remain new obligations.

## 5. Compressed experiment and a nonredundant H9 choice

The five cores alone attain every allowed odd through262143. Adding the
gapless double-triangle diagnostic core supplies no new target in this
range. Across the six-core experiment the90+64 raw cuts have98 distinct
nine-coordinate signatures. This is **FINITE-EXACT support for a new small
source grammar**, not a theorem beyond the stated bound and not a range
record for the repository's much richer older constructions.

The exact list of targets through262143 not in the H1/H3/H5-parent image is

    173,335,663,1383,2373,2457,9495,9761,10271,13861,15033,
    19243,21901,30597,40685,42147,42223,42503,46363,46629,
    48541,53707,55555,62917,63165,71377,73355,77233,80423,
    86641,88177,90449,153741,156739,174993,178051,213447,215923.

The counts below2^k, k=8..18, are

    1,2,3,4,6,6,11,14,25,32,38.

`hole-repairs.json` contains every full prefix list, both core repair
searches, and the actual row masks of each chosen repair. The repair chosen
within each fixed core minimizes total flank bit length, then A+B and the
remaining tuple lexicographically. It is not a minimum-order certificate
among arbitrary tournament constructions.

**Both H9 core types matter, even though their H and their uncoloured
(k_+,k_-)=(75,-9) agree.** With the first three cores retained, removing
S5b leaves21901 uncovered, while removingS5a leaves86641 uncovered. The
inverse bounds make these two particular nonmemberships all-order results
for the restricted grammar. Each H9 type covers37 of the38 listed holes.
Concrete missing-coordinate repairs are:

* N=21901: C=S5b, O mask12, A=0, B=703; signature
  `(9,75,-9,50,-6,31,-3,21,3)`.
* N=86641: C=S5a, O mask9, A=9, B=109; signature
  `(9,75,-9,43,-5,32,-4,18,2)`.

The74 displayed repairs were counted directly by subset DP, not only by
(1); the largest chosen graph has19 vertices. Their minimal-order repair
cuts are diverse, and the smaller flank parameter can be43, so the data
do not support a single-cut or one-flank replacement. The status is a
positive two-chart construction, not two interchangeable H9 states.

## 6. Binary continuation: exact affine pullbacks, no claimed automaton

For fixed A and a fixed coloured core, (1) is

    F(B)=a g(B)+b B+c(B mod2)+e,

with explicit coefficients

    a=d, b=(d+k_+)A+o_+,
    c=(d+k_-)α-o_-,
    e=d(g(A)-1)+i_+A-i_-α+j_++j_-.

The recurrence for g gives exact binary pullbacks

    F(2q)   = a g(q)+(a+2b)q+a(q mod2)+e,
    F(2q+1) = a g(q)+2b q+(b+c+e).                     (4)

Thus the coefficient maps are

    (a,b,c,e) -> (a,a+2b,a,e),
    (a,b,c,e) -> (a,2b,0,b+c+e).

These are genuine algebraic digit rules for the actual graph compiler.
They are a useful starting point for a normalized carry-state proof, but
the slope b grows, so they do not themselves yield a finite-state
recognizer or establish that all target words are accepted. Any claimed
normalization must preserve target membership and retain the core/cut
chart and both flank addresses.

A cheap predefined test of the sparse hole list gives zero transitions
`h -> 2h+c` for each c in{-3,-1,1,3} among the32 holes whose candidate
children remain within262143. This rejects that particular carry-free
hole-closure guess; it says nothing against a more elaborate finite-state
language with carries. The next constructive problem is to use the two
H9 charts to close the complement of the first three charts under a proved
normalized version of (4), rather than extending finite ranges.

Additional reproduction:

    python3 -B 05-knowledge/results/session251_creative_completion_20261111/ternary/repair_holes.py

`python3 -B 05-knowledge/results/session251_creative_completion_20261111/ternary/reproduce.py`
runs all three commands in normal and optimized Python and fails on any
receipt difference. `replay.json` pins the three source-script hashes and
all three identical receipt hashes. This checks only the new constructions
and their declared finite arithmetic universe.
