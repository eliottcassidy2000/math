# Hamilton-path guards reconstruct every finite tournament

2026-10-26. **PROVED, with a CITED structural dependency:** G1-G3 below.
**FINITE-EXACT:** the accompanying reconstruction and graph controls.
No literature-priority claim. This resolves the all-order question left
open in [the preceding packet, C5](../strong_atom_continuation_vitali_20261026/README.md).
It does not give a small quotient for unrestricted named re-marking:
the minimum quotient is discrete.

**Incoming convergence.** While this packet was being checked, main advanced
to `a8fcb1c840` with [Session70's guard proof and audit](../session70_attachment_guard_alignment_20261026/session70-guards/REPORT.md).
It proves the same all-order reconstruction theorem from the same classical
structural mechanism, using vertex-deleted reachability. G1 below gives a
direct two-step-count decoder, O(n^3) once the trace support is supplied.
These are two formulations of the same reconstruction result, not separate
classification claims. The linked incoming packet also independently
obtains the graceful guard quotient discussed in CONNECTIONS.md.

## 1. Inheritance and the question

The closest proved mechanism is the preceding packet's operation-relative
continuation congruence C5 and its all-arity relative-alignment formula C1.
The hostile is C3's seven-vertex pair with equal rooted source, rooted target
and colour, and C4's independent amplification to 2^k presentations at
7k+2 vertices. The corrected near miss is to assume that every tournament
arc occurs in a Hamilton path. Even strong tournaments fail that test.
The neglected sidecar is the **support of the whole named path family**,
including the constraint imposed by the absent arcs of that support.

For a tournament T on a specified finite vertex set V, let P(T) be its set
of directed Hamilton words. Define its trace support D(T) by

    a -> b in D(T) iff a,b occur consecutively in some P in P(T).

Thus D(T) is an oriented spanning subgraph of T, obtained directly from
the positive answers to all legal re-marking guards. It need not be a
tournament. We show that D(T) nevertheless determines T.

## 2. The external theorem actually used

Arthur H. Busch, *Arc-Traceable Tournaments*, PhD thesis, University of
Colorado at Denver (2005), Theorem 4.12, printed pages45-47
(PDF pages53-55), is the cited dependency. The primary
[institutional PDF](https://digital.auraria.edu/downloads/ma7yx-pcr25/AA00005831_00001.pdf)
was read directly. Its source hash and the exact passages checked are in
[SOURCES.json](SOURCES.json). This is not Busch's different 2006 minimum-H
theorem; no minimum-path-count formula is used.

**CITED structural input.** If an arc u->v of a strong tournament S is in
no Hamilton path, there is z such that S-z has at least four ordered strong
components, with u in the first and v in the last. Only this part of
Theorem 4.12 is needed. The extra incidence condition in that theorem is
not needed here.

The thesis proves this using its Theorem 4.7 and a one-vertex separator.
The later paper *Extremal Results on Arc-Traceable Tournaments*, JCMCC63
(2007),3-15, has a [publisher record](https://combinatorialpress.com/jcmcc-articles/volume-063/extremal-results-on-arc-traceable-tournaments/).
That record confirms the non-arc-traceable setting, but the proof below
uses the specifically checked thesis statement. We make no claim that
this elementary reconstruction consequence is absent from prior literature.

## 3. G1: explicit all-order reconstruction

**Lemma (middle vertices).** Let S be strong and let S-z have ordered
strong components C1,...,Cr, r>=3. Every vertex in C2 union ... union C_(r-1)
can be deleted without destroying strong connectivity.

**Proof.** Strong connectivity supplies z->a for some a in C1 and b->z
for some b in Cr: no other component sends an arc into C1, and Cr sends
no arc into another component. Therefore C1 union Cr union {z} induces
a strong subtournament, since C1->Cr. Every intermediate vertex has an
incoming arc from C1 and an outgoing arc to Cr. Adjoining any selection of
these vertices preserves strong connectivity. In particular, deleting
one intermediate vertex preserves it. This proof also covers singleton
end components.

If S-w is strong, every arc incident with w is traceable. To realize u->w,
take a Hamilton path of S-w ending at u and append w; to realize w->v,
prepend w to one starting at v. Such prescribed-end paths exist by
cutting a Hamilton cycle of S-w. This is the elementary argument also
recorded as Busch's Proposition4.1.

**G1 (PROVED with the cited input).** Every finite tournament T is uniquely
reconstructed from D=D(T) as follows.

1. Compute the strong components of D and their reachability order.
   Between different components orient every pair in that order.
2. Retain every arc already in D.
3. For a missing pair {u,v} in one strong component S, count

       c_D(u,v) = |{w in S : u->w and w->v are arcs of D}|.

   Orient u->v exactly when c_D(u,v)>=2. Exactly one of c_D(u,v),
   c_D(v,u) is at least two.

**Proof of component recovery.** The strong components of T form an
ordered join. Every Hamilton path traverses them in that order and
restricts to a Hamilton path of each. Conversely every selection of such
local paths concatenates to a global one. Hence D restricted to a strong
component S is precisely D(S). A nontrivial S has a Hamilton cycle, all
of whose arcs occur in local Hamilton paths and therefore in D. Thus
S remains strong in D. No reverse cross-component arc has been added,
and adjacent components are connected by traceable cross arcs. The
components and their order are consequently recovered exactly.

**Proof of the missing-pair rule.** Suppose u->v is absent from D(S).
Apply the cited theorem: S-z has r>=4 ordered components with u first,
v last. Choose distinct middle vertices a,b from two intermediate
components. The lemma makes both non-cut vertices, so the four arcs
u->a, a->v, u->b, b->v all belong to D(S). Therefore c_D(u,v)>=2.
Every path from v to u in S must pass through z, because S-z orders
u's component before v's. In particular at most one vertex can be the
middle of a two-step path v->w->u. Thus c_D(v,u)<=1. This proves the
rule and its uniqueness. The singleton cases are vacuous.

The decoder uses **the original trace support**, not arcs inferred earlier
in a sequential closure. It runs in O(n^3) time once D is supplied.
Obtaining D from a full guard table can require reading factorially many
entries; this is not an O(n^3) oracle-query bound.

**Small hostile.** The five-vertex presentation with rows
`(2,4,9,19,7)` has nine Hamilton paths. Its only invisible arc is3->1.
The visible two-step routes are3->0->1 and3->4->1; the reverse has only
1->2->3. The decoder restores the missing arc. Strongness alone would
not justify taking D itself as T.

## 4. G2: the smallest closed named attachment carrier

For n>=1, let X_n be tournaments in the fixed Hamilton-path order0,...,n-1. Their
M=(n-1)(n-2)/2 nonconsecutive pairs are free bits. For every permutation p,
R_p is legal precisely when p is a directed Hamilton word; illegality
has a distinguished failure outcome, as in the preceding C5.

**G2 (PROVED).** The map from X_n to its complete named R_p guard signature
is injective for every n. Consequently every failure-aware operation-stable
quotient supporting these named re-markings is discrete on X_n. Its exact
number of states is 2^M and its optimal worst-case fixed-length encoding
uses M bits, with n supplied separately.

These are counts of **valid presentation states**. If partial operations
are totalized by adding an absorbing failure state, the total is 2^M+1;
the failure convention is not silently included in the M-bit count.

Indeed equal guard signatures give equal path sets, equal trace support,
and then equal tournaments by G1. If two states differ, a single named
guard distinguishes them; no long complement word is necessary. Conversely
the M-bit literal state implements re-marking, reflection, complement and
the earlier ordered-run assembly exactly. Complement flips all M bits.

This lower bound concerns the **whole literal carrier**. The earlier
conditional alignment-orbit index inside a specified pair of rooted
endpoint types is still the appropriate fibrewise discriminator. Counts,
fresh uniform re-marking kernels, unnamed operations and quotient actions
may admit smaller states: they are different protocols. The theorem does
not claim that every matrix entry must be physically stored, or that
variable-length coding cannot save on special families.

**G3 (PROVED, all-arity aligned forest recovery).** Suppose an ordered
spanning path cover has complemented relation B on the same named vertex
set as A, with its ordered endpoint list retained. Legal Hamilton-word
guards at A reconstruct A; guards after forest complement reconstruct B.
Their common directed arcs form exactly the path forest. The endpoint
list orders its components, retaining singleton coincidences. Thus these
depth-at-most-two observations recover the aligned ordered-cover state.
This statement uses the forest-complement operation of C1 at higher arity;
at arity one it is the actual metagraph complement. Without that operation
in the alphabet, the higher-arity conclusion is not asserted.

The independent multiplicity aggregation of Session69 P3 remains valid.
G2 says why that aggregation cannot serve as the full state for arbitrary
named continuation. Frobenius ancestry is an additional requested label,
not a substitute for the reconstructed relation; the phase transport F1
in the preceding packet is unchanged.

**The 81/27 hostile survives, and the guards separate it.** The independent
checker reconstructs both twelve-vertex members of the preceding packet's
control pair. They have the same 67 labelled cyclic triples and the same
462915 Hamilton paths by count, but different trace supports and therefore
different named path languages. Their automorphism orders 81 and 27 are
inherited from that packet's proof, not recomputed here. Equal triangle
data and even equal total path counts cannot replace the named guard data.

## 5. Mechanism and connected threads

[CONNECTIONS.md](CONNECTIONS.md) records the additional Vitali recovery,
the exact graceful factorial analogue, square-sum guard failure, and the
K5/Petersen/K3,3 maps. The key mechanism is **independently addressable
choices**, not a large symmetry order by itself. G1 adds a complementary
lesson: a missing coordinate in a support can still be forced by the
native class constraints. This forcing succeeds for tournaments and
fails for arbitrary simple graphs.

## 6. Reproduction and audit boundaries

Run from the repository root:

    python3 -B 05-knowledge/results/guard_reconstruction_20261026/verify.py

The producer reconstructs every fixed-path tournament through order seven
without a strong-source filter. The independent audit imports no producer
code: it uses subset dynamic programming for trace support and a different
component/closure implementation, checks all labelled tournaments through
five vertices, all fixed-path presentations through six, the full strong
representative atlas at seven, and specified larger controls. The graph
checks enumerate actual paths and literal edge replacements. Exact universes
and counts are in the receipts. Both programs are replayed normally and
with Python optimization; no assertion disappears under -O.

The proof has been self-audited against singleton boundaries, directed
versus undirected paths, labelled versus ordinary objects, the exact
external statement, and the non-traceable five-vertex hostile. Independent
code paths are computational corroboration, not an independent human or
agent proof review and not formal verification.
