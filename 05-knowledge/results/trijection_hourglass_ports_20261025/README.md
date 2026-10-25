# Consecutive-port closure, affine cut certificates, and four-bit Fano codes

2026-10-25. **PROVED** scoped structural statements with **FINITE-EXACT**
controls and independent audits. No literature-priority claim. Complete
endpoint-four and all-order black-core connectivity remain **OPEN**.

## What changes the tournament frontier

1. [Uniform consecutive-port theorem](endpoint/GENERAL_TAIL.md): for every
   k>=4 and n>=2k+3, a closed blue marked path with incoming interior ports
   I={n-k,...,n-2} has a literal three-black detour. A single symbolic
   certificate handles the uniquely forced cyclic case at n=2k+3.
   The [terminal three-port theorem](endpoint/TAIL.md) is its endpoint-four
   case: for every n>=11,
   a closed blue marked path with incoming interior ports
   I={n-4,n-3,n-2} has a literal three-black detour. There is no
   closed-middle assumption. The proof treats the two directions of
   x versus4; one uniquely forced order-eleven source has an explicit
   exceptional certificate. This closes the former h=n-4 boundary.
2. [Arbitrary-port cyclic-prefix theorem](endpoint/THEOREM.md): for n>=11,
   I={h,y,c}, c->b, 3->c and (n-h)->x give another three-black detour.
   Its useful chord can be near the far end, beyond early-prefix guards.
3. [Overlapping cuts](trijection/README.md): the odd-prefix bound becomes
   n>=max(11,2ell+3), improving the inherited max(15,2ell+6) bound. Three overlapping
   positions are handled by an all-length permutation-orbit proof.
4. The same packet supplies an exact affine compiler and finite-menu
   coverage criterion. Every source filter, directed join, intermediate
   complement and positional reversal is retained. Inconsistent signed
   cycles certify coverage without enumerating all source tournaments.
   Two positional reversals can both be anti-automorphisms of some
   tournament iff their product has odd order; prescribed paths remain
   a separate constraint.

The historical order-ten target is already closed by the inherited
[endpoint-three affine certificate](../endpoint_three_order10_20261021/README.md)
and the all-layer census in [THM-4732's result packet](../black_detours_closed_blue_lines_klein_dictionary_20261020.md).
No new giant census is needed here. The new results keep exact vertex
labels; neither ports nor clone classes are contracted. Equal literal
target tournaments imply equal clone relations as well as equal arcs.

For n>=15, combine these results with the
[inherited endpoint-four frontier](../endpoint_four_frontier_20261021/README.md)
and [second-port closure](../ports_owners_f1_20261025/endpoint.md).
The terminal h=n-4 case is gone; the remaining closed-middle residual lies
in 3<=h<=n-5 with c->x and y->x, and must fail the applicable old and new
prefix/cyclic guards. This is a necessary description of the remainder,
not a claim that every such pattern is realizable or unrepaired. Complete
endpoint-four at orders11..14 and unrestricted connectivity remain open.

## Geometry and the recovered 6–4–6

[The code packet](codes/README.md) integrates the incoming Session60
surface trichotomy as recovered mathematics: 16 numerical Campedelli,
12 properly elliptic, and two Enriques minimal models. It independently
checks the branch-code and local singularity mechanism.

The new extension is a lossless four-bit model of the sixteen good codes.
With a labelled seed and opposite-pair order, the crossing-star group acts by

    (n,sigma):(v,epsilon) -> (n+sigma(v),epsilon+parity(sigma)),
    v in F2^3, epsilon in F2.

All48*16 transport identities are checked literally on Fano planes.
This group acts transitively with stabilizer order3. Consequently no
subgroup preserving the three crossings has orbit partition6+4+6.
The actual geometric dihedral subgroup has the inherited4+12 orbits.

[The incidence continuation](incidence/README.md) recovers the genuine
6–4–6 as dimensions of three line-incidence spaces, and proves a growth
theorem rather than identifying different sixteen-element objects.
The [common-kernel refinement](incidence/common_kernel.md) isolates the
additional modes caused by physical line identifications.
[The bridge ledger](bridges.md) retains the actual F1 carrier and the
Monster/Baby Monster fixed-point regimes, with their sources and losses.

## Inheritance, portfolio and concept board

The session began at clean82691d391b, after fetching the live shared
surface. Incomingf7041582a8 andc24c123d86 were read before finalizing
geometry; their stronger classifications supersede our initial open
minimal-model question. Those overlaps are corroboration, not new claims.
The final clean rebase also integrated the six/seven continuation through
e4d54b0131; its branch-cube criterion and resolved-reflection constructions
are linked in the bridge ledger, with their separate carriers retained.
The missing local.machine-id prevented the legacy message processor;
no registration or message state was invented.

Closest proved mechanism: literal path-edge XOR plus actual positional
anti-reversal. Canonical hostile: an order-seven pure-blue target, and
the order-eleven source where the generic tail word really stays blue.
Corrected near miss: a legal permutation/XOR ledger need not have black
moves. Least-used sidecar: signed constraints on unordered arc variables,
and the exact endpoint neighbourhood inside a permutation orbit.

**Anchor:** endpoint-four and black connectivity. **Niche:** branch-code
quotients with their symmetry action. **Wildcard:** the recovered
three-direction incidence kernel and modular line collisions. The broader
repository LRC(14) anchor remains open; this packet supplies no reduction.

| Concept | Source/target operation | Preserved predicate | Missing information or decisive test |
|---|---|---|---|
| Movable cuts | Three literal K operations | Directed words and exact target | Actual blackness at every intermediate state |
| Ordered graph chart | Orientations to pair bits | Every labelled arc | Observer and constant parity must stay attached |
| Reflection product | Two anti maps to an automorphism | Orbit restrictions | Odd order alone does not make prescribed joins legal |
| Fano code quotient | Rank-six deck to rank-three labels | Crossing stabilizers | Branch matroid and code kernel, not seven labels alone |
| Four-bit good-code carrier | W/C3 to F2^3 x F2 | Entire W action | Labelled seed and opposite-pair ordering |
| Physical incidence box | Grow2-by-L-by2 support | Points and adjacent cubes | New shared lines create new common-kernel modes |

The pull and each computation were compared against this board. The
successful common move is to retain the constraint that a quotient would
otherwise erase. This is a scoped methodological connection; no tournament
or arithmetic theorem follows merely from the surface counts.

## Verification and next proof target

[AUDIT.md](AUDIT.md) records the independent proof review and exact universes.
Run `python verify.py` and `python -O verify.py` from this directory for
the publication gates and source pins. Producer checks and the independent
bit-mask audit run in ordinary and optimized Python; their receipts agree.

The next tournament target is a finite legal-menu cover of the remaining
3<=h<=n-5 affine cells, particularly those with c->x,y->x but no usable
prefix chord. The signed-cycle criterion specifies exactly what a
structural certificate must prove. A surviving cell from one menu is not
a non-leaf bridge and not a counterexample to black connectivity.
