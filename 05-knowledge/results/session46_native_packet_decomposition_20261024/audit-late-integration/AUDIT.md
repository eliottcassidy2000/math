# Independent audit of the late Session 46 integration

24 October 2026. **PASS within the stated bounded comparison scope.**
This addendum reviews incoming commit
`06306318651dfcddbdff9e9790794f11b0e39f94` and the late integration note.
It does not replace or modify the frozen original synthesis audit. Incoming
objects were read with read-only `git show`; no checkout or Git state was changed.

## Proof-level findings

### The same integer kernel does not imply the same legal routes

Incoming `SIX_TRADE_LATTICE.md` uses the same 28-source carrier and marked charge
map. It selects 27 actual four/six-slot columns with determinant 2048. The earlier
Session 46 witness instead gives integer combinations of 48 actual moves to
produce an index-2048 basis. Both prove the same containment-and-index argument
in the zero-margin lattice Z^27, not merely equality of rational ranks.

Its question about universal six-slot binary-fiber connectivity is answered
negatively by the earlier Session 46 complete cardinality-four census. This is
a legitimate comparison within the same authenticated supply, not a deduction
from matching group orders. The eight-slot hostile allows every one of the 28
slots to vary. No additional external source inventory is covered by that claim.

The integration note initially described an “enlargement” ambiguously. The
reviewed version states correctly that the move lattice enlarges at six slots,
while its rational annihilator shrinks from dimension five to the single
margin direction. Finite index then identifies the remaining paired-charge
obstruction. None of this establishes bounded-box connectivity.

### The incoming 18-source circuit is a different object

The incoming 45-clock relation has different source sets from both the 43-clock
and cost-optimal 37-clock Session 46 relations. Source 1923 occurs here, so its
extra target-1559 path must be excluded by the zero target margin rather than
ignored. The saved complete native proof does so. The 37-clock minimum remains
valid in its fixed supply; the existence of a different 45-clock circuit does
not change that optimization statement.

The earlier scalar b_1638+b_1974 is zero at both incoming endpoints and is
therefore blind. The incoming c vector annihilates all 27 four-slot moves and
has endpoint values -3 and -1. Its four-slot lower bound continues to hold when
all 28 supplied coordinates vary, without a one-active-helper restriction.
The “-6 at both” value in `CHARGED_RETURN.md` refers to a different concurrent
scalar and is not attributed to the Session 46 critical-pair invariant.

The source-1638 helper route uses four six-slot events with clock costs
14,15,15,13, totaling 57. Four events are minimal at slot budget six: a route
must change all 18 original slots, and an active restored binary helper adds
at least two changes, exceeding the capacity of three six-slot events. The
number 57 is the work of the displayed shortest-event route, not a proved
minimum cumulative clock cost. Its 42-state complete one-helper fiber and the
217-state one-active layer comparison remain distinct state spaces.

This helper is reversible borrowing, not strict temporal nonreturn. Retention
of the original all-dyadic bank is also not exact preservation of the unscaled
integer X current. The incoming endpoint debt is already -256 at U. The
source's explanation correctly retains diagonal degree terms when splitting
an intermediate current between the old slots and its active helper.

### The multi-target extension is sound with its full hypotheses

For t distinct equal-charge targets and sources outside that target set,
endpoint target-margin equality makes the changed-slot directed multigraph
balanced. A decomposition into edge-disjoint simple directed cycles retains
actual slot identities, uses at most t edges per cycle, and preserves every
integer target count. Same-target record changes are harmless one-edge loops
under the supplied local-record hypotheses.

For each cycle, the remaining target-owned charge vector sums to zero. It
belongs to the kernel of G^t -> G, identified with G^(t-1) by retaining t-1
coordinates and reconstructing the last. Equal target charge already makes
the source rows safe. The sum of cycle charges is zero because the endpoint
physical currents agree. A minimal nonempty zero-sum subcollection has at
most D(G^(t-1)) cycles; extract it and repeat. Thus each event uses at most

    t D(G^(t-1))

slots, and every differing slot changes once. For a finite abelian p-group,
the classical formula gives D(G^(t-1))=1+(t-1)(D(G)-1). Consequently D(G)=132
gives 264,789,1576 for t=2,3,4. The p-group formula is not extended to arbitrary
mixed-prime groups; the more general Davenport bound still applies there.

The essential realization hypothesis is that every mixture of the supplied
old/new records satisfies all fixed shared-coordinate, full-lift, local-guard,
resource and remainder conditions. Without it, a zero-charge cycle packet
need not be feasible. The general exact-incidence/Graver framework remains
available but does not retain the same charge-only bound. There is no claim
of native sharpness, arbitrary native supply, a bound independent of t, or an
order-independent constant.

The original marked label fiber is exactly {570,1965,2234}; it is not an extra
path menu. The complete retained 57-path inventory has targets 1965,2234 and
1559 only, so target 570 supplies no native choice under those pins.

The inherited all-source census strengthens this boundary only within its
explicit grammar. Its saved 700 menus, independently authenticated in Session
45 across 3,528 canonical merged sources and 1,102,392 actual paths, contain no
target-570 pair and exactly the same 28 menus for targets 1965/2234. Both the
original indexed-arc generator and the separate endpoint-fiber verifier were
read. They require four-word normalized black resources, no source-target loop,
two distinct equal-charge targets, the same actual first/last owners, and some
common internal indexed arc k,k+1 for 1<=k<=5, in the fixed canonical frame.

Consequently changing only the original source cannot supply a target-570
transition while retaining that complete menu grammar. At least one condition
or the ambient setting must change. This does not exclude standalone paths to
570, blue-resource moves, other pin geometries, overlapping source/target roles
or general native wiring. This audit scans the inherited exhaustive output;
it does not rerun the 1.1-million-path census.

The grammar-dependent target graph on {570,1965,2234} therefore has one edge
1965--2234, with 28 source-menu witnesses, and isolated vertex 570. This is not
the original metagraph: the original node table gives converse classes 732/736
and fibre 590, while direct reading of all original edge rows gives 213 distinct
neighbors, black weight 590, no blue contribution and no loop at 570. The absence
of eligible transitions is thus not absence of incident ports. The combined
grammar obstruction does not single out a particular pin requirement as its
cause.

The note's topology-sensitive corollary is also valid under exactly the same
realization hypotheses. In each active nontrivial target-graph component C,
the actual transition multigraph is balanced and its total target-owned charge
is zero separately. Its simple directed cycles have length at most ell_C, the
maximum simple cycle length in the bidirected component, with two-cycles
included. Zero-sum extraction in G^(|C|-1) therefore bounds component events by
ell_C D(G^(|C|-1)). Components can be processed separately; all-mixture
compatibility makes this scheduling legitimate. Same-target changes have
support one and are handled separately, including at isolated vertices.
For K2 plus isolated 570 the nontrivial component gives 264, not an available
three-target 789 construction. This extension was checked at proof level;
no new native supply, sharpness or efficient longest-cycle algorithm is asserted.

## Executed checks and explicit limits

`check.py` is a new independent implementation of cycle decomposition followed
by minimal zero-sum packet extraction. Normal and optimized runs each pass
**250,248 explicit checks**, with byte-identical output. They cover 28,454
accepted ordered endpoint pairs for C2/C3, two/three targets and up to four
slots, with every intermediate integer margin and target-owned charge checked.
A nine-slot control has three prescribed cycle debts forming a minimal zero-sum
sequence in C2^2. That control does not claim the underlying endpoint trade is
primitive or the nine-slot bound natively sharp.

`check_incoming.py` passes **156 explicit checks** in both modes, with identical
output. It reads exact commit objects and imports no incoming verifier. It:

- independently regenerates the complete 27/1,917 small-trade sets;
- computes determinant 2048 by rational Gaussian elimination, independently of
  the incoming fraction-free determinant method, and verifies all 27 columns
  are actual selected moves;
- recomputes the full charge image based at the first source, four-move rank 23,
  the incoming c invariant and the old invariant's blindness;
- checks the saved helper route's exact margins, marked group charges, literal
  endpoint clock costs and restored source-1638 coordinate;
- checks the complete three-target label fiber and absence of target-570 paths
  in the already authenticated common-clock domain, and scans the inherited
  complete 700-menu census for the stronger grammar-scoped exclusion;
- verifies the distinct original-metagraph and pin-compatible target-graph
  statistics at node 570;
- compares six incoming normal/optimized receipt pairs byte-for-byte.

The late comparison checker was separately replayed in both modes: **76 checks**,
with output identical to its supplied receipts. Its exhaustive incoming circuit
guard counts are 48,620 balanced, 696 mod-8 and two full-dyadic states. This replay
is labelled separately from the new algorithms above.

These controls do not enumerate fresh Hamiltonian paths, re-multiply the incoming
full Laplacian image witnesses, or rerun its entire route/transport/spectral
corpus. Those native claims retain the incoming proofs and audits as antecedents.
No stronger scope is inferred from receipt agreement.

## Reproduction and source boundary

    python3 -B check.py --out normal.json
    python3 -O -B check.py --out optimized.json
    python3 -B check_incoming.py --repo REPO --out incoming.normal.json
    python3 -O -B check_incoming.py --repo REPO --out incoming.optimized.json

`check_incoming.py` requires the named commit to be available locally and only
reads its objects. Replay the sibling late-integration checker using that
packet's documented arguments, writing outputs here if desired. Compare each
normal/optimized receipt pair. `SOURCE_PINS.json` also pins the reviewed final
README and integration note; the manifest lists only this audit's deliverables.

Preservation always means that every original harmonic guard remains valid
(kernel inclusion), not equality of harmonic kernels or critical groups.
Selected-packet signature retention is not a reclassification or rewiring of
the fully occupied canonical atlas. The final README keeps these boundaries
and the original frozen synthesis as a historical checkpoint, with the late
comparison linked separately.
