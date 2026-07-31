# Fragment decode: "the depth-2 tree is 3+3+1 = 7, not 9" (klein-S691, 2026-07-28)

**The fragment (owner puzzle bundle, truncated mid-word):**

> The depth-2 tree is 3+3+1 = 7, not 9: the μ₃-fixed branch P₀ has a
> degenerate fiber — two p

Precedent: the eq(27) owner-fragment episode (see auto-memory
`eq27-owner-fragment-decode`) — fragments arrive stripped of context and the
decisive move is reconstructing the *generating computation*, not
free-associating. Three candidate readings, each with a prediction and a
cheap disambiguation:

## Reading A — recursion tree of an arithmetic-Kakeya construction (preferred)

The Katz–Tao slice framework works with triples (g(0), g(1/2), g(1)) on a
line, i.e. coordinates (a, b, −a−b)-style with a natural **S₃ action whose
cyclic part μ₃ rotates the three slices**. A depth-2 ternary recursion
(3 branches per level) naively has 9 leaves; if the branch fixed by the μ₃
rotation (call it P₀) has a **degenerate fiber — two p[oints/preimages
coincide]** — the tree counts 3 + 3 + 1 = 7. In certificate-score terms
(see `04-computation/ak_forcing_engine.py`): a collision in the symmetric
branch is an *identification*. SIGN CORRECTION (late-session, empirically
settled): identification per se SHRINKS the denominator n − t, which alone
RAISES the score; its real value is that a merged class inherits all its
members' incident edges and becomes a species-rich junction hub, making
forcing cheaper — the numerator savings (m + r) outweigh the denominator
loss (measured exactly: merge worth 12 = 9+3 versus edge+seed 14 = 10+4 on
the twin [2,2,2] witnesses at 12/7 vs 7/4).
**Prediction:** the intended improved-score gadget is a depth-2, 3-branch
constructible graph whose symmetric branch carries an identification;
searching structures with dims [3,3] plus a forced two-vertex coincidence
(emulated by T-seeding or by label-sharing) should show a score drop.
**Test:** run `ak_strict_search.py` class D on [3,3]-shapes with an added
"identified pair" gadget; compare best scores with/without.

## Reading B — Ward-recurrence tree in vandermonde-snp

The vandermonde-snp Lean project verifies non-SNP-ness of Vandermonde powers
via **Ward recurrences** on Jack/Macdonald data indexed by partitions. A
depth-2 recursion over partition branches where the branch fixed by a ℤ/3
symmetry degenerates ("two p[artitions collide]") would likewise count
3+3+1 = 7 nodes. **Test:** read the paper's recursion section
(`paper/vandermonde-power.pdf` in that repo) and check whether a 7-node
depth-2 tree with a collapsed symmetric branch appears.

## Reading C — arboreal preimage tree of a cubic with a critical collision

For a degree-3 map, the depth-2 preimage tree of a base point has 9 leaves
generically; if the μ₃-symmetric branch passes through a **critical value
(two preimages collide)**, the count is 7. This is standard arboreal-Galois
bookkeeping. No current in-repo consumer; parked.

## Verdict and routing

A and B are not exclusive — the owner may be pointing at the *pattern*
(symmetric-branch degeneracy buys a count reduction) rather than one host.
The AK workbench treats it as a design heuristic (Reading A) because that is
where a 9→7 count *provably changes a score denominator*. Nothing here is a
claim; this note exists so the next session does not re-derive the candidate
readings.
