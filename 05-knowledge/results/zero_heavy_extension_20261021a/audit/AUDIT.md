# Independent audit: two-unit corridor and first owner correction

Date: 2026-10-21. This audit concerns only the separate extension frozen by
manifest SHA-256
`a6667f7a99d5fa186ba09b5b6df462c13d084599f1af2c0001aba9f9aff4c390`.
It does not modify or supersede the earlier main-packet audit. Its base
manifest remains
`a5490b17a42062db222676a44a1eb2ef1c3f0b55499b2389cc70691886f9112c`.

## Verdict and exact scope

**PASS:** the common binomial formula extends through `0<=delta<=2k-1`
for `k>=3`. At `delta=2k`, for `k>=4`, the exact weighted-class count is
`L_k(2k)+K_k-2H_k`, with the strong variant using the strong marked-edge
incidences. Generally this correction is `-3(k-1)(2k-3)!!`.

The `k=2,3` first-owner results are separately FINITE-EXACT, not an
unacknowledged extension of the at-most-two-unit-gap proof. The complete
`k=4` C/D check is owner-family-only, not a complete mass-layer census.
No claim for `delta>2k`, unrestricted tournament enumeration, or a new
palette-forcing theorem is included in this audit.

## 1. The corridor retains its grammar

For every good order, `delta<2k+1` forces exactly `k` backward edges and
unit endpoint weights. Thus `delta=(k+1)d+kh+U`. Below `2k`, the only
families are full matching gaps, a matching with one empty gap, and a
full-gap P3-plus-matching. The earlier stars-and-bars calculation really
is unchanged: the full-gap inclusion-exclusion has at most one violating
part; the missing-gap residual is at most `k-1`; the P3 residual is at
most `k-2`. The latter two residuals do not exceed the cap `k-1`.

The extension is therefore a quotient theorem, not a changed coefficient
formula. All weight-above-one vertices are forward in every good order,
so their transitive order is intrinsic. Their incidence signatures fix
ordered unit cells. Large weights may vary without weakening this order.

## 2. Every two-unit and mixed cell has the claimed owner

For a full matching, three unit gaps cost `3(k-1)>2k-1` for `k>=3`.
There are at most two. Separated exceptional cells can be handled without
circular partner identification: a cell without an internal matched pair
is transitive, hence ordered intrinsically; a cell with an internal pair
is C3 and has no external matching edge. If two cells are cross-paired,
they are both transitive, so both endpoint roles are fixed first.

Adjacent interior unit gaps create the five-position cell
`e1,u1,e2,u2,e3`. With no internal edge it is transitive. With an internal
edge, the remaining endpoint has one actual external partner, lying in a
fixed singleton cell. For internal pairs 12, 13, 23, that marked endpoint
has local outdegree 0, 2, 4 respectively. These are distinct pointed
objects. The middle case is rigid: the marked vertex splits the remaining
four positions, and the unique reversed cross-edge fixes their roles.
The other two cases have only the harmless internal C3 rotation.
At an exterior end the four-position object is transitive, or a source
joined to C3, or C3 joined to a sink at the opposite fixed end. No numerical
matching/gap record is changed by its local triangle rotation.

In the missing-gap family the residual budget allows at most one unit
gap. If separated, the missing-gap two-endpoint cell has a forward edge
after clone rejection, which fixes its roles before any external partner
is used. Thus an edge joining the two exceptional cells causes no circular
argument. If adjacent internally, the four-position words are
`e1,e2,u,e3` and `e1,u,e2,e3`. The four legal one-edge cases have marked
outdegrees 2, 3, 0, 1. The checker independently exhausts all pointed
permutations and verifies these four classes are distinct. Forgetting
the actual external port identifies the two extreme-edge-reversal cells;
retaining it separates them. This is an actual hostile to an unpointed
quotient, not a hypothetical warning.

With no internal edge the local order is transitive and the fixed external
partners retain the endpoint roles. At an exterior overlap, the allowed
internal edge makes C3 in `e1,u,e2`; the alternative consecutive endpoint
edge in `u,e1,e2` is exactly the clone rejection. Otherwise the local
transitive order and outside partners distinguish the two patterns.
Strong-end cells are the same interior four/five-position cells, not new
unexamined exterior cases.

The P3 family still has every gap at least two below `2k`, so it retains
the old rigid decoder. The three families have different total vertex
counts, so there is no cross-family collision in this corridor. These
arguments prove the enlarged all-order quotient, not merely its finite
controls.

## 3. At delta=2k, all and only A/B/C/D enter

The nonnegative identity permits exactly
`(d,h,U)=(0,0,2k),(0,1,k),(1,0,k-1),(0,2,0)`.
A lost backward edge or two shared endpoints still costs too much.
The four presentation counts in first_owner.md, equation (6), are correct:

- A adds the two-violating-parts inclusion-exclusion term `E binom(m,2)`.
- B subtracts the `m-1` forbidden residual parts of size `k`.
- C distributes `k-1` among its `m-1` occupied gaps with no cap violation.
- D chooses two allowed holes. A matching with `a` adjacent matched edges
  admits exactly `binom(m-a,2)` hole pairs. Expanding this expression
  yields `E binom(m,2)-(m-1)J+K`.

The sum equals `L_k(2k)+K`; the other correction terms cancel exactly.
This identity counts retained presentations only. It is not yet a class
count. In particular one must still subtract two for each nontrivial
owner fibre.

For two holes, the clone criterion remains complete. Separated holes
produce two forward endpoint cells unless a forbidden adjacent matched
pair occurs. Adjacent holes produce three endpoints: an adjacent internal
matching edge is the rejected case, an extreme internal edge gives C3,
and no internal edge gives a transitive triple with its outside partners
retained. There is no extra hidden clone contraction.

## 4. Complete C/D ownership, with the actual external port fixed

For `k>=4`, A has at most two unit gaps (`3(k-1)>2k`), and B at most
one (`2(k-1)>k`), so the previous cell analysis applies even when other
large gaps have been lowered. A and B each remain injective, and their
vertex counts differ from C/D.

C has no unit gap or exactly one. With no unit gap, its graph is rigidly
recovered and it has one more large vertex than D. With one unit gap,
all other gaps have weight exactly `k`. Its exceptional interior cell is
transitive, an unmarked C3 from an isolated matching component, or a C3
with exactly one marked vertex from a P3 edge. The mark is the P3 centre's
actual other backward edge to a fixed outside singleton. In the last
case the marked vertex can occupy either endpoint position, producing
exactly two C records of the same pointed C3.

D has no self-identifications. Separated holes create forward ordered
pairs. Adjacent holes create an ordered triple or a pointed C3. In that
C3, the vertex with the external matching edge must be the middle
endpoint. Thus there is exactly one D record. Exterior-hole overlaps
are ordered singletons or pairs and do not introduce another owner.

These are the only C/D collisions. Distinct unit-cell size patterns cannot
be identified by an isomorphism fixing the intrinsic large order. Where
the patterns match, a transitive D triple has one distinct external
matching edge at each vertex, whereas the transitive C cell has a gap
vertex with no external backward edge. Exterior two-cells similarly
have external-degree patterns `(1,1)` in D versus an actual gap and an
endpoint in C. These cannot agree under their unique transitive orders.
An unmarked C3 cannot become D's C3 with its external port retained.

Therefore every nontrivial fibre is exactly `(C,C,D)`, not an automorphism
factor and not a clone merger. The independent audit parses all eight
literal collision examples (`k=2,3,4,5`, general/strong) and reconstructs
their rows and weights. For each of the two C targets it explicitly finds
the unique isomorphism of the local triangle which fixes *every outside
vertex* and maps the actual distinguished port correctly. It does not
merely compare a scalar invariant or a port-erased canonical key.
Strongness is invariant under these literal isomorphisms, so the same
owner rule holds in the cut-crossing subfamily.

## 5. The marked-edge incidence corrections are exact

A D owner is specified by a marked edge `(i,i+2)`, with both intervening
gaps empty. Its middle endpoint must have an external partner because
its two neighboring endpoints are already matched. Conversely every such
marked edge supplies one valid pointed-C3 owner. Distinct D owners remain
distinct by the preceding decoder. Hence exactly H fibres each lose two
presentations.

Generally there are `2k-2` positions for the span-two edge and
`F_(k-1)` remaining matchings, so `H=(2k-2)F_(k-1)`. For K, choose two
nonoverlapping adjacent-edge slots in `binom(2k-2,2)` ways and match the
rest. This gives `K=(k-1)F_(k-1)`, and `K-2H=-3(k-1)F_(k-1)`.

For strong K, deleting the two marked adjacent pairs preserves cut
irreducibility twice. Conversely insert two indistinguishable marked
pair blocks into the `2k-5` interior gaps of an irreducible matching on
`2k-4` points. Both may use one gap; stars-and-bars gives
`binom(2k-4,2)I_(k-2)`. The marks are an unordered pair, so no extra
factor two belongs here. The exceptional `K^s_2=0` is correct.

For strong H, a marked span-two edge and its intermediate vertex belong
to one unique irreducible cut block. Arbitrary matching prefixes and
suffixes give the two unrestricted F factors. Coefficient extraction
gives precisely the stated convolution recurrence, including
`H^s_2,H^s_3,H^s_4,H^s_5 = 2,8,60,592`. This does not assume ordinary
connectedness of the matching graph and uses no fitted sequence.

## 6. Independent finite tests and provenance

The auditor imports only its own standard-library support from the prior
independent main audit, bundled as core_audit_support.py. No producer
checker or producer canonicalizer is imported or executed. The extension
canonicalizer uses weighted directed colour refinement followed by *all*
unresolved permutations, not the intrinsic-gap decoder being tested.
Every admissible record is checked for literal clones, mass, weighted
isomorphism ownership, and (where applicable) actual reachability.

The exact finite universe is:

- all declared corridor cases `(3,4),(3,5),(4,6)`, general and strong,
  plus the additional strong `(4,7)` layer;
- all four families at `(k,delta)=(2,4),(3,6)`, general and strong;
- all C/D records at `(4,8)`, general and strong, with full fibre histograms;
- all endpoint matchings through `k=7` for K/H and their strong analogues;
- 58 presentation-cancellation checks via independent bounded polynomial
  coefficients for `k=4,...,32`, both parameter rows;
- every pointed permutation of the crucial four/five-position cells,
  the port-erasure hostile, and all eight literal whole-graph collisions.

Incoming `zero_heavy_next_20261020h`, commit
`e265436678941418bc34a9eb3b697acd989524d0`, already proves deficit k,
corrects the shared-endpoint onset, and supplies the k=2 full-matching
two-unit-cell mechanism. The extension now credits and pins its README,
AUDIT, and receipt. That inheritance is appropriate. The extension's
further content is the wider corridor, mixed missing/unit cell distinction,
and first-owner correction. The older heavy-coordinate-deficiency owner
mechanism is also credited without identifying the two parameters.

A reproducibility issue was found and resolved transparently: incoming
work revised the old charge_palettes zero_heavy.md onset sentence, so the
current repository no longer matches the *unchanged* base manifest's hash
for that one dependency. The exact historical file was recovered read-only
from base commit `c4a7875652a5ea516e93a31c3c6410d2204773b2` and bundled as
frozen_base_zero_heavy.md. Its SHA-256 is
`34fc50fbc4c5d6c58d8f8a683a4b7d74b2c712347be80d063396a556ecc0eebf`.
The audit authenticates that original version for the old base, and the
new incoming files under their own new pins. It does not silently repin
or mutate the base, claim both versions coexist at one commit, or treat
source drift as a mathematical failure.

All gates are explicit exceptions and remain active under `python -O`.
Normal/optimized receipts and the portable verifier record exact results.
No repository or remote file was changed, and no post-freeze extension
edit is certified by this audit.
