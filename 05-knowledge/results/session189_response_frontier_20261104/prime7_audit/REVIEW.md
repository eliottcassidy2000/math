# Independent review of prime-seven ambient converse separation

**Status: approved for the stated quotient-substitution and fixed-ambient scopes.**
The two producer proofs are mathematically compatible and together cover every
prime seven-name quotient with arbitrary nonempty physical children and
profile-preserving target replacements. The sole remaining all-even carousel
branch in the profile packet is precisely covered by the separate modulo-sixteen
packet. No repair was needed. This does not assert unrestricted named-response
injectivity or separation between arbitrary quotient pairs.

Reviewed source pins:

- `profiles/PROOF.md`: SHA-256
  `ea7ca89df44ba732247291b00aadb6a70cd8f4c4fb13624ae51f63479502d9b8`.
- `prime7/PROOF.md`: SHA-256
  `e382ce566c061626cb309467c58928f2a7b4a4e387835166bde252ccce7f7672`.

The review read both full proofs and the relevant producer implementations.
The independently written programs in this directory import no producer code.
Both were run normally and with optimized Python; each pair of outputs was
byte-identical and reported PASS. All data inputs have SHA-256 pins in the
receipts. The additional physical controls retain complete physical row lists
and all first-five integer cover counts, not just a reported residue or hash.

## 1. Classification and physical names

`owner_check.cpp` independently checks all 2^21 labeled tournament encodings
against the complete retained owner array. It uses an incremental adjacency
update rather than the producer's per-code reconstruction. Exactly 973,680
labels are prime. It also checks the owner of all 197*7! = 992,880 representative
permutations. Each of the 197 input rows is separately verified to be a valid
prime tournament in Python. Therefore the retained orbits cover every prime
label and cannot overlap across representatives. Their automorphisms merely
repeat the same owner and do not cause overcounting.

The proof explicitly transports quotient labels, entire physical child sets,
and both masks together. Thus an orbit representative is not a license to
relabel source and target independently.

All 2,758 integer endpoint coordinates were independently recounted by 7!
permutation enumeration. The partition into 159 nonflat and 38 endpoint-flat
representatives, used by the inherited endpoint criterion, is consistent with
these full integer counts.

## 2. Singleton-sensitive anchor theorem

The proof correctly retains the exact singleton/non-singleton set V. A
non-singleton has normalized q2=2 modulo four; a singleton has q2=0. The
three-visit coefficient is even by swapping its two nonempty return segments,
including when the repeated child is the constrained start. This swap preserves
all run-word arcs and endpoints. Consequently the truncation to one repeated
child is valid modulo four in both the start-count and one-port calculations.

There is no missing normalization factor: the one-port converse contrast is
even, so the common odd product fixes it modulo four. Under p_i=0 the start
skew is even too, so doubling it removes the odd normalization ambiguity
modulo eight. Applying the exact inherited start identity gives the displayed
anchor formula with the stated signs. In particular, the implementation's
use of converse endpoint counts to obtain start-run counts has the correct sign.

A separate Python completion recurrence recomputed all 3,990 primitive anchor
integers and then all 4,864 pattern records. It verified all ten residual
patterns and the exact coverage of the six singleton-containing residuals.
As an additional literal control, `verify_residual_physical.py` counts the
12 actual singleton/pair extensions at their six witnesses. Every complete
first-five contrast is (2,0,0,0,0) modulo four. The largest such control has
13 physical vertices. These extra audit controls are distinct from the
producer's smaller-model interpolation route.

## 3. Odd-total ambient safety

The ambient argument was checked against the pinned Session188
`mod8_transfer/PROOF.md` and `quotient_cycle_parity/PROOF.md`, not inferred from
a Hamilton scalar gap alone. Both ports inherit the root's exterior edges,
so the augmented core is genuinely a common module. Its parity is the old
core parity because two ports are added.

For odd order, the needed identity is

    Delta q2 = Delta q3 = 2 Delta H - 4 Delta Z (mod 8).

The safe-query filter is applied to **actual query values**, not just to
polynomial coefficient positions. An independent XOR endpoint dynamic program
recounted the Hamilton-cycle parity of every one of the 6,728 bare quotient
extensions on the full 841-query basis, retaining the closing edge to the
unique root port. Every retained parity and safe index matched.

An independent reconstruction of the integer coefficient banks into actual
values, followed by independent binary rank calculations and enumeration of
all 2^14 type assignments for each of the four residual quotients, found:

- 445 safe queries, rank 14, no safe-blind type;
- 589 safe queries, rank 14, only all-E safe-blind;
- 421 safe queries, rank 13, only all-P and all-E safe-blind;
- 517 safe queries, rank 14, no safe-blind type.

Every safe-blind type has even total order. Therefore every odd-total type
has a safe actual gap-four query and a pure first-coordinate augmented-profile
gap. The anchor's odd-order argument is also valid: at query (e_i,0), every
Hamilton cycle must use x->y, leaving an old path ending in child i; its
converse parity difference is p_i=0.

## 4. Modulo-sixteen closure and full-cover interpolation

The even-child domain was checked against the pinned Session188
`even_cyclic_mod16/PROOF.md`. Its actual upper domain is exactly sufficient:
normalized coordinates are (1,14+8e,4+8f,8,8g) through index five, and all
later coordinates vanish modulo sixteen. The formal bits need not all be
realizable; proving the identities on the larger domain is valid.

The new reduction preserves cuts. In particular, a cut permits consecutive
occurrences of the same child, whereas an uncut adjacency does not. It does
not incorrectly apply a return-segment involution across separate paths.
Products of two perturbations vanish modulo sixteen, and every other child
must have one visit in a first-order term. The 120 baseline shapes plus seven
single-child five-visit shapes form a complete downset for the coefficient
functional. No large physical child is silently replaced by one with an
unmatched full profile.

All six stated identities were checked in every one of their 110 channels:
five baseline residues and 21 derivatives for each of five cover coordinates.
The two carousel queries have complementary first-coordinate residues
8(e3+e5) and 8(1+e3+e5). Their higher-coordinate contrasts are zero, so an
arbitrary ambient multiplies the selected residue eight by its odd contracted
Hamilton count. The common odd child normalization fixes eight as required.
For h>=6, sixteen divides h!, completing the infinite-profile statement.

The audit independently checked the complete 127-shape interpolation domain
and reevaluated all 1,524 retained physical instances from their exact
nonnegative kernels and integer Stirling profiles. A third physical
implementation then directly recounted 72 selected named graphs on at most
13 vertices. It uses a backward endpoint recurrence with the five exact
cover coefficients packed in base 2^64, rather than the producer's forward
vector recurrence. Every count matched; the bound
13!*binomial(12,4)<2^64 prevents coefficient carry. Full row lists and counts
are retained in `physical_controls.json`.

## 5. Scope and reproduction

The source-relative root-family extension in the modulo-sixteen proof uses
an additional, explicit inherited hypothesis: hereditary original-module
profile monotonicity. Its converse-complement dual pair gives opposite
integer baseline contrasts, and nonnegative substitution preserves a chosen
negative direction. This is not all-pairs orbit injectivity. Exterior
operations still need the stated highest-node restoration argument.

For reproduction, restore any losslessly packed producer data first, using
its manifest. The audit supports the published sibling layout `profiles/`,
`prime7/`, `prime7_audit/`, and the original workspace layout. Then run:

    python3 verify_audit.py > normal.json
    python3 -O verify_audit.py > optimized.json
    cmp normal.json optimized.json
    python3 verify_residual_physical.py > residual_normal.json
    python3 -O verify_residual_physical.py > residual_optimized.json
    cmp residual_normal.json residual_optimized.json

The installed C++17 compiler is used only for the independent owner census;
no installation is performed. These are exact finite checks and mathematical
proof review, not proof-assistant certification. No observation about arbitrary
prime-eight quotients or unrestricted response collisions follows from this
review alone.
