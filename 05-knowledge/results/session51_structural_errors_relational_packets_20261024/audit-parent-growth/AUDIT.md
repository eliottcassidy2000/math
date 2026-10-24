# Independent audit: parent-only growth and the order-ten prediction

Status: PROVED selection/action argument; FINITE-EXACT independent atlas and
trace certification. Production claims pass. No publication was performed by
this audit. Ordinary classes are used; converse classes are never merged.

## Results

The independently certified order-nine atlas contains 191,536 classes:
191,480 have a stable singleton and 56 have none. Its exact labelled mass is
68,719,476,736 = 2^36 and its refinement error is 758. Every supplied action
list is the complete automorphism group, not merely a closed subgroup.

The independent order-ten scan evaluates 97,213,472 pointed states and gives
U10 = 97,336,068. All 191,536 per-parent histograms agree with production.
The separately proved and replayed correction is E10 = 5,508, giving

    a10 = (97,336,068 - 5,508) / 10 = 9,733,056.

This is not an all-order scalar recurrence or a prediction from an unmarked
metagraph. The source retains the entire parent adjacency and permutation
representations. Neither the new atlas size nor a10 is an input to the
independent atlas engine. The incoming trace compares against production only
at the level of independently computed per-parent histograms.

## Structural proof audit

Parent subset orbits classify pointed extensions exactly: deleting the mark
recovers the parent, and changing its identification acts through the full
parent automorphism group. Stable refinement's ordered color labels are
intrinsic, so its first singleton is a unique fixed vertex. Selecting precisely
the state rooted there chooses one pointed class over every child with a
singleton. All child automorphisms fix this vertex, and restriction identifies
the full group with the attachment stabilizer in the parent group.

For singleton-free nine-vertex children, the only odd-cell patterns are 9 and
3+3+3. The deletion key (canonical eight-parent index, least attachment mask)
is a complete invariant of the pointed child, not merely a color signature.
The least re-rooted key is one automorphism orbit, hence one pointed state.
This selects exactly one representative of each remaining child.

For each root v attaining that key, the parent transporters h(S)=M_v give
child maps g(root)=v and g(u)=f_v^{-1}(h(u)). They preserve every arc, and any
child automorphism has a unique such image v and transporter h. Thus the
construction returns the full group once, without target-nine canonicalization.
Closure and labelled mass are additional checks, not substitutes for these
completeness and distinctness proofs.

The independent path intentionally DOES canonicalize order-nine audit objects.
This is solely an audit, not production computation, and is marked as such in
its receipt. It does not canonicalize order ten or enumerate its full atlas.

## Independent algorithms

`independent.cpp` uses incoming-edge equitable partition refinement, recursively
splitting each cell by incoming count vectors. It enumerates permutations within
its intrinsic ordered cells and minimizes a 36-bit upper-triangle encoding.
All minimizers recover the full automorphism group. This differs from the
producer's outgoing signatures, row-array key and eight-only deletion transport.

It independently checks all 6,880 order-eight parents and all 191,536 order-nine
outputs for pairwise nonisomorphism and exact full action-list equality. Summing
n!/|Aut| then proves coverage of the complete labelled universes, 2^28 and 2^36.
As a second construction, it enumerates all order-eight subset orbits, selects
an incoming-refinement singleton, and deduplicates the small residual by its
own canonical key. The reconstructed canonical set equals the production set.
This regenerates 1,716,608 pointed states, including 232 residual pointed states,
without assuming those counts. The run visits 36,941,580 canonical labellings.

`incoming_trace.cpp` works with incoming bit rows, starts from one color, packs
count signatures in integers, and groups by equality in first-occurrence order.
The numerical color names need not be intrinsic because only cell membership
and root-cell size are used. Retaining the old color prevents merging; unchanged
cell count therefore certifies stability. This differs from the producer's
outdegree initialization and lexicographic array sorting. Normal execution may
stop on a singleton; O3/NDEBUG/FULL_REFINEMENT disables that shortcut. Both modes
compare every parent separately, not merely the final scalar count.

The common histogram is:

    size 1: 97,153,964
    size 3: 58,118
    size 5: 1,160
    size 7: 60
    size 9: 170

## Correction and a further stratified check

The physically copied structural producer manifest is
905dd04456f9b4cf302e8b71c7bdd915db209ba081b4735b0f7a733378866f43.
Its independent structural audit manifest is
81381e8c41aac5bdb57d0b8c8a7744c45e0684e090ec03e854249ca91efd0068.
All 36 structural-audit payloads, including its nested 15-payload producer,
were byte-verified and physically copied. Its manifest-only replay reruns both
independent structural checkers and both producer checkers normally and with
Python -O. It confirms strong E10=3992, total E10=5508, 160 positive-error
classes and raw second moment 229416.

The independent full nine-atlas audit further gives strong error by stable-cell
size {3:72, 7:14, 9:630}, totaling 716. Its all-nine counterpart is
{3:72, 7:56, 9:630}, totaling 758. Applying the already-proved SCC placement
identity to these values and to the frozen strong-ten witnesses gives total
order-ten error by stable-cell size {3:720, 5:3360, 7:168, 9:1260}.
Subtracting within each root-size stratum yields root-orbit trace contributions
97,153,964; 173,634; 2,440; 252; 270, totaling 97,330,560.

One wording-only defect was corrected before freezing: four order-nine rigid
3+3+3 classes have two singleton paddings each, producing eight order-ten
classes. They are not eight order-nine classes. The numerical 32+8=40 and
size-three correction 720 were already correct.

## Replay, hostile controls, and provenance

The final production manifest is
161449ddf5c2a07b7afefd213438b85c74897850b89050dde66428819b0838d6.
All 29 payloads and the manifest itself are actual copied bytes in `producer/`.
The corrected source includes checked integer parsing before allocations.
The source pin chain authenticates the original certified JSON-eight atlas and
compares every text input row and action to that JSON, not only its class count.

Both production programs are compiled from this frozen copy with O2 and with
O3/NDEBUG; the latter trace also disables singleton early exits. The rebuilt
nine-atlas bytes and every trace TSV match their frozen originals exactly.
The two independent engines are replayed in the same two modes. Normal and
optimized independent atlas receipts are identical; trace mathematical fields
and complete per-parent outputs agree, while work counters differ as expected.

Five semantic mutations fail in both compilation modes: choosing the singleton
by labelled vertex rather than intrinsic color; omitting residual classes;
omitting root-moving transporter cosets; omitting an attachment; and corrupting
the independent root-cell size. Four more meaningful input tests check duplicate
same-mass classes and omitted full group actions in both modes. Eight malformed
input tests cover empty, nonnumeric, negative and truncated input in both modes.
These are 22 deliberate rejections. In particular, duplicating one rigid class
while removing another preserves labelled mass, but still fails nonisomorphism.

`manifest.json` lists the deliverable payloads and their exact SHA256 hashes.
`verify_audit.py` validates this envelope, the embedded envelopes and the
mathematical cross-checks. `replay.py` makes a new manifest-only physical copy,
then rebuilds both independent and production programs and reruns the structural
packet, mutation controls and stratified reconciliation. No binary is trusted
or shipped as evidence. Machine-dependent run times are stored in receipts,
not asserted as asymptotic improvements.
