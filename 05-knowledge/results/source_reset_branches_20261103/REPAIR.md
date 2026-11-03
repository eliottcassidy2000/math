# Two actual native moves repair the next branch, with a strong-core boundary

**FINITE-EXACT complete searches in the specified universes; PROVED operation
and observation identities. Full ternary closure remains OPEN.**

## R1. The stronger adaptive obstruction

Start at C3=(2,4,1). A birth at any old root r adds x with r->x and
x->every other old vertex. The [unique-ancestry theorem A1-A4](../adaptive_metagraph_reformulation_20261103/ANCESTRY.md)
makes the first-root-0 histories a disjoint list of ordinary classes.
Every birth strictly increases H. Consequently forward expansion, retaining
every child of H<=295 and rejecting only larger H, is complete at every
vertex order for this count cutoff. The resource guards in the program
are failures, never silent pruning criteria.

There are2820 classes, at orders3 through11 in counts

    1,1,4,20,120,750,1133,687,104.

All25196 candidate births were counted directly. No numerical observer
deduplication was used. For each source h=31,95,97, mask bit j denotes
the existence of some nonempty adaptive descendant of H=3h+2j, j=0,1,2.
Reverse propagation on this finite complete tree gives:

|mask|H31|H95|H97|
|---|---:|---:|---:|
|0|0|2|6|
|1|0|4|6|
|2|1|1|3|
|3|0|2|0|
|4|0|0|1|
|5|2|1|1|
|6|0|1|0|
|7|0|0|0|

Thus the three H31 classes do not have a common source for95 and97,
even with arbitrary old-root access. Their root histories are(0,3,1),
(0,3,0,5), and(0,3,4,5). This strengthens the inherited two-live-port
obstruction, not the global scalar spectrum. H95 has one source supporting
287 and289 in this language; H97 has none supporting293 and295.

The older three-command observer remains valid in its own operation
language. Allowing a source reset changes the physical graph and its
observer; merely renaming a previously forgotten root was insufficient.

## R2. Actual ears from the original31 source

Use the inherited marked source

    T=(2,4,41,115,103,67,7), root0,
    v(T,0)=(31,14,25), source word(a,a,r5).

For a cut I, add vertex7 with i->7 precisely when i belongs to I.
Exhausting all128 cuts yields exactly

    H95: I=40 or61;       H97: I=32 or63.

These are general actual ears, not the original single-bridge affine
states(95,41,75) and(97,40,75). Keeping the source relation records that
difference. A complete repeated-one-ear expansion, allowing only the two
target counts at each step, has13 labelled nodes and7296 evaluated cuts.
Every one of the four first children lacks at least one immediate next
branch. The tree terminates at depth3; the requested depth4 has no nodes.
It is therefore a complete failed four-level test, not a successful tree.

The all-cut count is obtained from the full ordered pair incidence:

    H(T+x)=sum_(b notin I) start_b(T)+sum_(a in I) end_a(T)
            +sum_(a in I,b notin I) W_ab(T),                (9)

where W_ab counts ordered two-path covers of T ending at a and starting
at b. This is the inherited mixed-insertion identity. The C++ program
computes W from separate subset start/end tables and evaluates all cuts;
every accepted child is independently recounted by ordinary Hamilton DP.
It does not replace ordered incidences with total cover mass.

## R3. Four exact reset certificates

Write K_P for actual path complement: retain the directed Hamilton spine P
and reverse every other arc. Physical names stay fixed throughout a reset.
The following two-spine programs return to the original H and admit both
subsequent target ears. The table gives cuts in the **unstripped** eight-
vertex endpoint; all paths and row masks are in `receipt.json`.

| Initial cut | H along reset | First path P | Second path Q |Next +2 cut|Next +4 cut|
|---|---|---|---|---:|---:|
|40|95,157,95|2,3,5,7,4,6,0,1|4,6,0,2,1,7,3,5|120|158|
|61|95,197,95|3,4,0,7,6,1,2,5|6,1,7,3,2,5,4,0|184|84|
|32|97,359,97|3,4,6,2,5,7,0,1|7,0,2,1,6,5,3,4|140|189|
|63|97,207,97|3,4,1,7,6,2,5,0|4,1,7,0,6,2,5,3|18|122|

The intermediate counts are deliberately included: the return is neutral
at its endpoints, not at every step. Each program uses two complements
and no birth, deletion, auxiliary vertex or free class-valued substitution.
All Hamilton paths were checked on the actual current graph.

The two H95 endpoints are strong. The H97 endpoint from cut32 has universal
source7; that from cut63 has universal source4. Removing that single name
leaves a strong seven-vertex core of the same H. Root0 survives in both
cases. The projected cuts are respectively(12,61) and(2,58), producing strong
eight-vertex293/295 children. The H95 children287/289 are strong of order9.

The deletion can be executed as a named operation, or the universal source
can be retained with a guard that it continues to beat every subsequent
new vertex. Both policies preserve the claimed counts. It must not be
silently forgotten if arbitrary subsequent native complements are allowed.
One of the deleted names is the original newborn7; hence this is not an
unqualified theorem preserving an arbitrary ordered pair of protected ports.
The golden port0 is explicitly protected. The full source and deletion
sidecar remain in the certificate.

**Exact small minimality.** For every one of the four inputs, neither zero
nor one native complement has a same-H endpoint with both immediate ears.
Two suffice, so the minimum is two for this precise reset predicate. For
each H97 input, exhausting all distinct endpoints at native radius<=2
(13944 labelled graphs) finds no successful **strong** endpoint. This does
not exclude three moves, another target operation, or another H97 source.

## R4. A universal source can be deleted without changing golden phase

**Proposition.** If z is a universal source of U, r!=z, and T=U-z, then

    (H(U),E(U,r),C(U,r))=(H(T),E(T,r),C(T,r)).             (10)

All Hamilton paths of U start at z; delete that prefix to get the first
two equalities. In a C cover, the marked path starts at r and therefore
cannot contain z. The other path starts at z. Delete it there; the other
path may become empty, which is exactly the allowed Hamilton-path term
in C(T,r). Conversely prefix z to that path. This is a bijection.

Consequently universal-source deletion preserves the entire three-coordinate
single-bridge observer, all five Q4 responses, and the F49/charge state.
This is also literal: under a subsequent single-bridge attachment at a
port other than z, z remains a universal source. The proposition does not
hold for deleting an arbitrary source component, a universal sink, or the
marked vertex. It does not commute with an unguarded native complement.

## R5. Transport is retained, not assumed neutral

For any actual rooted graph use v=(h,x,y) from the inherited golden
observer and encode

    alpha=4(h+2x-y)+y theta in F49,
    g=h-y in F7,        theta^2=4theta+2.

Record alpha as a pair of F7 coordinates, its theta exponent only when
nonzero, and g separately. Every endpoint and intermediate graph was
tested by literal attachments of **all five** q0,q1,q2,q4,q5 presentations;
the resulting full integer observers equal multiplication by the inherited
matrices. In particular a later q2 sends alpha to theta alpha and fixes g.

The short native resets themselves have nontrivial recorded transport:

|Initial cut|v before|v after, and after permitted source deletion|theta exponents before/after|charge before/after|
|---|---|---|---|---|
|40|(95,44,77)|(95,34,70)|16 ->0|4 ->4|
|61|(95,63,80)|(95,43,75)|14 ->20|1 ->6|
|32|(97,45,78)|(97,36,65)|35 ->44|5 ->4|
|63|(97,63,84)|(97,56,67)|40 ->3|6 ->2|

These are source-dependent changes, not a single claimed field automorphism.
The inherited Frobenius law theta^7=4-theta is unchanged. If a reset has
integer discrepancy Delta, every later single-bridge context transports
it as Delta P_word. The actual graph and legal spines are still necessary
to select and implement the reset; three H observations alone do not
certify a native operation.

**Stricter phase requirement: a decisive obstruction.** Exhaust all two-
native return endpoints of both H95 children, allow every vertex as the
new marked port, and require both immediate affine ears. No endpoint has
v congruent to the original affine interface(95,41,75) modulo7. Since
(alpha,g) is an invertible encoding of v modulo7, preserving that exact
joint clock is impossible in this template. This remains false after
allowing re-marking; it is not caused by a fixed choice of root0.

The H97 scans do have congruent marked endpoints, but the retained examples
have a universal sink and E=0 at the selected marks. They do not repair
the missing H95 half of a joint clock-preserving bank. E=0 is a real
realization datum: there is no Hamilton path ending at that port, although
the general single-bridge count identity still permits it. Do not turn
this endpoint guard into a false blanket prohibition on using the graph.

Thus the positive short reset retains and computes transport; it does not
solve the stronger demand for a phase-neutral ternary bank. The latter
needs a larger native language, a different source, or explicit transverse
phase corrections with a proved realization rule.

## R6. Why the small repair is not yet an infinite family

If two native paths have the same predecessor and successor at a vertex v
(with the same endpoint convention), replace v by k transitive clones and
expand that common run in both paths. Both moves remain legal and their
endpoint is the corresponding cloned target. This supplies an all-k
**operation** lift. H neutrality still requires equality of the clone
polynomials; the two desired branch identities require two additional
polynomial identities. These are independent obligations.

Using the inherited [Session155 clone-slot formula](../session155_spectrum_20261103/primes/README.md),
the complete two-step template scans for the inputs from cuts40 and32 give:

|Input|Two-path words|Neutral words|Words whose endpoints have both ears|Distinct eligible(target,v)|Neutral clone lifts|
|---|---:|---:|---:|---:|---:|
|95,cut40|17385|155|6|5|0|
|97,cut32|17857|156|10|6|2|

Repeated endpoints account for the difference between word and endpoint
counts. In both neutral clone lifts, H(k)=56*2^k-15, k>=1. Neither has
even one of its available k=1 ears transport as the required affine
identity for all k. Exhaustive clone coefficients certify this; no fitted
sequence is used. The scans cover only untouched single-vertex cloning
and cuts uniform on that clone module. They do not rule out partial cuts,
multiple modules, longer native words or different source resets.
