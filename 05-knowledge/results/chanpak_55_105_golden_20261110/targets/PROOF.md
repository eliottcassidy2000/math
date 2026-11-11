# Actual-source compatibility at 55 and 105

**PROVED:** the padding-compatibility reduction P1 below. **FINITE-EXACT,
with all-order consequences:** P2–P4 use the inherited audited complete
strong inventory below375, the cited strong minimum, and P1. **REFUTED:**
the proposed strengthening that one actual source can support both ternary
returns at every realised odd height at least95. Height105 is a permanent
counterexample. The weaker rule allowing different sources remains **OPEN**.

Local labels P1–P4 apply only to this file. Counts are ordered directed
Hamiltonian vertex words. Converse classes are not identified. A cut mask
has bit `i=1` exactly when the newly adjoined vertex points to old vertex i.

## Inheritance and the question changed by the experiment

The closest proved mechanism is J1 in
[the previous joint-deletion proof](../../decomposition16_odd_spectrum_20261110/induction/PROOF.md):
for an extension with source/child counts h,N, all old SCCs outside the
new vertex's SCC contribute a common factor dividing gcd(h,N). The
canonical hostile is the golden H55 source family: its scalar arrow
55→111 exists elsewhere, but no attachment to either golden source realizes
it. The corrected near miss is uniform module contact implying contiguous
Hamiltonian visits; the previous attachment packet supplies the explicit
27-versus15 counterexample. The least-used sidecar here is the **whole
source between its first and last nontrivial SCC**, rather than its count
or only the new vertex's SCC.

The board is actual source / strong child / exterior padding / joint cuts /
golden compiler / arithmetic factorization. The anchor is simultaneous
3h+2,3h+4; the niche is the exact source-class series; the wildcard is a
common-source lifting theorem for arbitrary finite sets of coprime targets.
The aim is procedural: reduce an infinite source family before computing.

## P1. Common sources are classified by a retained middle core

For a tournament S with H(S)>1, delete only the singleton SCCs before its
first nontrivial SCC and after its last nontrivial SCC. Call the remaining
tournament K(S). Thus uniquely, up to the usual isomorphism,

    S = TT_a → K(S) → TT_b.

Internal singleton SCCs are part of K(S) and must not be stripped. K(S)
starts and ends in nontrivial SCCs; it can have several such SCCs.

Fix h>1 and a finite list N_1,...,N_t>h, each coprime to h. For each N_i,
enumerate marked strong children (C,v) with H(C)=N_i and H(C−v)=h.
Normalize C−v as above, retaining a,b, an actual isomorphism to K, and
the complete cut of v. Let B_i be the set of resulting K isomorphism
classes. Then:

> There exists one actual source S of count h with a one-vertex extension
> to every N_i if and only if the intersection of the sets B_i is nonempty.

More precisely, a source TT_A→K→TT_B admits branch i iff some recorded
padding pair (a,b) for that K and branch has a≤A,b≤B. Consequently the
branch-ready padding set is an upper set in the product order on N².
Its minimal pairs, together with actual cuts and embeddings, are sufficient
to construct all the branches. For finitely many specified witnesses,
take A=max_i a_i and B=max_i b_i.

**Necessity.** Apply J1 to each child of the same S. Since gcd(h,N_i)=1,
its nontrivial part is a single strong child C_i; all vertices outside
C_i form exterior transitive prefix and suffix. In particular they are
singleton SCCs in S. Removing those exterior vertices changes neither
K(S) nor any internal singleton gap. The normalization of C_i−v_i is
therefore K(S), for every i.

**Sufficiency.** Choose witnesses with the same K and transport their
arcs by the recorded isomorphisms. Add A−a_i new universal sources before
the entire i-th child and B−b_i universal sinks after it. On deletion of
the distinguished vertex each resulting source is the same graph
TT_A→K→TT_B. Ordered joins with singleton SCCs preserve H, hence each
child still has count N_i. This gives all the cuts on one labelled source.
The coordinate inequalities give the sharper statement. QED.

**Scope and names.** The new vertices for the different branches are
alternative choices; this theorem does not put all of them simultaneously
in one child. Chosen named vertices of K are transported by the stored
isomorphism. Named exterior vertices require the ordered prefix/suffix
embedding, including where extra vertices were added. Other ports,
two-vertex attachment responses, F49/F343 phase, and the81/27 control are
not recovered from K or the padding pair alone. Preserve those sidecars
on the actual source and recompute them after any reset.

**Finite procedure.** The strong minimum tends to infinity, so every
B_i is finite and obtainable from a bounded strong census. For this
session max N_i=319, while the exact minimum at order12 is375 and is
increasing afterwards. Thus every strong child has order at most11,
and every reduced source C−v has order at most10. The inherited
[independent S50 audit](../../../../04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/README.md)
contains all required strong children. This is why absence of a shared
K here is an **all-order** obstruction, rather than a source-order cutoff.
The lower-bound import is
[THM-4744, interval chains and the strong minimum](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md).

The map is from marked strong children to source middle cores with padding
and cut embeddings. It preserves branch existence and the entire actual
source once the sidecar is retained. Dropping the cut loses future
attachment behavior; dropping K loses common-source compatibility. The
cheapest hostile is P3 below.

## P2. Complete source families and the 55 split

The inherited inventory has four strong classes with H55, all at order7;
H5 has one strong class at order4 and H11 one at order5. Ordered SCC
multiplication implies the all-order ordinary-class series

    F_55(z) = 4z^7/(1−z)^2 + 2z^9/(1−z)^3.

The first term is one strong55 core with exterior singleton runs. The
second term is either order of the strong5 and strong11 cores, with a
singleton run in each of three gaps. No other nontrivial factorization
of55 is possible. These are counts of source isomorphism classes, not
Hamiltonian paths or marked cuts.

Exactly two of the four strong classes admit55→167 and55→169. Each
supports both branches, without padding, and also admits55→111. The
other two support none of the three targets, even after arbitrary
exterior padding. No reducible H55 source supports any of the three.
In particular, the branch-ready source series for each of111,167,169
and for their simultaneous intersection is

    2z^7/(1−z)^2.

Here are explicit simultaneous ternary witnesses. The rows are outgoing
bitmasks in fixed labels0,...,6; the same row list is used by both cuts.

| Source rows | Cut to167 | Cut to169 |
|---|---:|---:|
|70,12,8,17,39,79,30|119|5|
|78,68,8,18,39,79,28|16|31|

The other two classes are precisely the literal golden outputs for
`s s q1 s` and `s q4 s s` from
[the golden55 calculation](../../golden_hamilton_transfer_20261103/CONNECTIONS.md).
This strengthens its source-language obstruction: exterior transitive
padding of either golden class still cannot supply111,167,169. The
golden scalar number55 does not encode which half of its source classes
the grammar reaches. An attachment repair must change that class.

## P3. Height105 permanently forbids a common ternary source

There are18 strong H105 classes: four at order7, eleven at order8,
and three at order9. The only other possible SCC factorization is
105=3·35: there are two strong35 classes, both at order7, and one C3.
The alternatives involving7 or21 are excluded by
[THM-1370, H-spectrum omissions](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Thus the complete source-class series is

    F_105(z) = (4z^7+11z^8+3z^9)/(1−z)^2
               + 4z^10/(1−z)^3.

The branch kernel sets have sizes3 for317 and1 for319 and are disjoint.
All their minimal padding pairs are(0,0). They are completely described
by the following fixed-label source matrices and cuts:

| Source rows | Target | All successful cuts |
|---|---:|---|
|82,76,9,17,38,79,28|317|68|
|102,84,8,19,37,78,28|317|93|
|94,20,8,18,36,79,158,63|317|15,60|
|158,140,24,16,34,79,159,60|319|26,31|

Their source-class series are respectively

    F_(105→317)(z) = (2z^7+z^8)/(1−z)^2,
    F_(105→319)(z) = z^8/(1−z)^2,
    F_(105→317 and319)(z) = 0.

By P1, **no tournament of any order with H=105 supports both one-vertex
returns317 and319**. Arbitrary cuts, relabelling, root switching inside
the same graph, and exterior transitive padding cannot fix the disjoint
source cores. A genuine source replacement is needed. In particular
the stronger form of the proposed threshold95 induction fails at105.
The previous threshold95 rule allowed different sources and is not
refuted. This result supplies a concrete boundary between those two
quantifier orders.

The four reducible3·35 source skeletons, and every internal/exterior
singleton padding of them, fail both branches. Thus105=3·35 is a valid
count realization but does not supply the necessary coprime return
ports. J1 forces the new vertex to fuse both nontrivial components into
one strong child; the complete child census shows that no such fusion
achieves317 or319. The factorization105=5·21 cannot supply a source at
all because21 is forbidden. Arithmetic factorization and a reusable
attachment construction are different predicates.

## P4. Consequence for an adaptive induction

A source-free theorem asserting both scalar edges is too weak to execute
a recursive construction on a fixed growing tournament. A universal
one-source-two-branches theorem is false. The surviving procedural
target is a source-sensitive switching system:

1. Retain the actual source core K, exterior names, and joint attachment
   data, together with the arithmetic target.
2. For each requested next branch, supply either a cut already valid on
   that source or a proved legal reset to a source in that branch's set.
3. Prove that resets preserve H and every advertised port/phase invariant,
   or explicitly account for their changes, and that the output supports
   the next recursive request.

P1 removes a spurious obstacle: different exterior padding requirements
can always be reconciled by coordinatewise maxima. P3 identifies the
real missing ingredient: incompatible **middle cores**. At105, any
reset serving both possible next digits must move between the displayed
317 and319 families. A mere re-marking or source/sink insertion cannot
be that reset. No bounded native reset between these classes is proved
here. Even a reset preserving H alone would not prove the required
joint-port or F49/F343 compatibility.

The historical105 atlas already warned against equating equal carrier
sizes and retained an actual source in its Rule105/free-input and carry
comparisons; see
[the mechanism atlas](../../number_105_mechanism_atlas_20261018q.md).
The precise transfer is that a recursive address language needs a
realizable, source-owned transition, not just an admissible output word.
No dynamical equivalence to cellular automata or a geometric equivalence
to square packings follows from the numeral105.

## Reproduction and hostile controls

Run from the repository root:

```bash
python3 05-knowledge/results/chanpak_55_105_golden_20261110/targets/check.py > /tmp/chanpak-targets.json
python3 -O 05-knowledge/results/chanpak_55_105_golden_20261110/targets/check.py > /tmp/chanpak-targets-optimized.json
cmp /tmp/chanpak-targets.json /tmp/chanpak-targets-optimized.json
```

The retained [receipt](receipt.json) and
[optimized receipt](receipt.optimized.json) are byte-identical. They
include the inherited bank SHA256, every relevant kernel's canonical
rows, all branch padding pairs, actual witnesses, and cut masks. The
universe is all strong children in the complete inherited H<375 bank
with counts107,109,111,167,169,317,319; all vertices are deleted, without
quotienting marked roots. Source canonicalization minimizes over every
actual Hamiltonian word, an exact isomorphism invariant for tournaments.
No automorphism-orbit count is inferred from the raw deletions.

Checks:300 selected child classes,2584 vertex deletions,343 independent
path-enumeration comparisons, all5376 cuts on all22 strong H55/H105
sources by a separate split-at-root endpoint/two-path-cover formula,
and12 exterior-padding lifts. The formula counts paths starting/ending
at the new vertex plus ordered two-path covers meeting the two cut
sides; it does not assume contiguous visits to a module. Permanent
35→107 and35→109 absences are negative controls. The two golden55
source graphs are reconstructed literally and all256 attachment counts
match the historical receipt. Bank completeness is inherited from its
audited ear enumeration, not rerun or claimed new here.
