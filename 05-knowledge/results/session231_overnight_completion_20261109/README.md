# Session 231: overnight constructive completion research

**Status: IN PROGRESS. Chan–Pak completeness remains OPEN.**

This ongoing packet contains independently audited partial results. It is not
a proof of completeness. The target is every positive odd tournament Hamilton-path
count other than 7 and 21. Base main is
`9c48b9ac1e2a9c22b0814686c76e3358ad56051d`.

## Inheritance and live concept board

1. **Physical descent.** Sessions228–229 give terminating physical-state and
   graft inverse procedures, but not a physically realizable inverse branch
   for every target. Revisit the older two-sided one-arc join and complete
   rooted triples rather than replacing states by scalar counts.
2. **Dense repair fibers.** Session229 gives context-stable adjacent pairs
   and triples, including unsaturated deletion cores. Repeated central
   switches scale by products and do not give unit intervals. Test changing
   deletion topologies and core positions, retaining the full endpoint data.
3. **Packing/carry operations.** Session230's positive cut kernel and S86's
   audited packing-polynomial work tie omission proofs to growth. Fixed
   sparse back-arc constructions, growing-budget variants and finite-width
   interfaces are separate candidate universes.
4. **Prime and creator structure.** A complete structural grammar is not an
   arithmetic coverage theorem. Birth fibers are not scalar representatives;
   W7 and the two63 creator orphans are inherited hostile examples.
5. **Independent verification.** Every promoted result needs a separate
   argument or implementation, explicit computational universe, hostile
   controls, and normal/optimized replay. Finite coverage remains finite.

## Fresh corrections applied to the research plan

The latest S86 audit (THM-4940, MISTAKE-739) and THM-4930 correction
(MISTAKE-738) were read before starting new derivations. Session229's
MISTAKE-736 separates core unsaturation from the scalar unit gap and from
unsaturation after closure. MISTAKE-737 separates count existence, birth
optimization and structural generation, and credits earlier directed-cut
and one-back-arc calculus. No original result is claimed in this checkpoint.

## Workflow

Parallel complementary construction rounds will be followed by independent
critical audits and nonforce mainline checkpoints. Exact scripts and proofs
will be added as they become ready. A finite computation, gcd-one family or
state-membership algorithm alone will not be described as completing the
conjecture.

## Audited checkpoint 2: two-sided physical descent

For a proper one-arc join of actual rooted triples (h,e,u), (k,f,v),

    N = hk - ef + uv,
    h,k <= floor((2N-1)/3) < N,
    uv <= N-1.

An exact finite divisor inverse uses the inherited physical-state oracle, so
no numerical cone point is admitted without a graph. Prime representations
normalize to strong outputs. This is a two-sided synthesis of inherited
counting, endpoint and membership mechanisms, not a new count formula or
coverage theorem. The two endpoint/response triples cannot be replaced by
scalar counts alone.

[Proof](descent/PROOF.md) · [Independent audit](descent_audit/AUDIT.md).
Both independent and producer implementations check all 1,207,801 pairs of
rooted labelled tournaments through order five; 5,625 smaller literal joins
are counted directly. The audit additionally checks strongness, prime
normalization and primitive interfaces. The enormous unrestricted membership
oracle is proved mathematically, not computationally exhausted.

The frozen audit input manifest describes the submitted candidate header.
Only its status and audit link were changed for publication, after PASS.
Subsequent construction rounds and independent audits remain active.

## Audited checkpoint 3: growing consecutive-odd intervals

**PROVED.** For every k >= 0, a single physical interface supports the k+1
consecutive odd counts D_k, D_k-2, ..., D_k-2k, where

    D_k = (245051264 * 4096^k + 1602510*k + 1118023) / 24843.

All witnesses are strong tournaments on 12k+13 vertices. The construction
moves one eight-vertex core through two fixed twelve-bit singleton words;
it does not repeat nontrivial cores. The common deletion core is saturated
(h=45 and gcd(h,q)=45), which is explicitly allowed here. The exact local
word identity and context law prove every k. Independent literal graph,
seven-coordinate and six-state path-cover checks agree. The first baselines
are 9909 and 40403039.

[Compact theorem](fibers/SHORT_INTERVALS.md) ·
[Standalone independent checker](fiber_audit/short12_audit.py).

This supplies unbounded consecutive repair width, the mechanism missing
from the earlier fixed-width fiber experiments. It does **not** select a
baseline over every target: these displayed intervals are separated by
exponential gaps. Universal Chan–Pak completeness remains OPEN.

Also audited in this checkpoint:

- [Sharper physical descent](descent/ROUND2.md): sole-exit zero-defect
  classification, equality cases and an improved two-sided descent bound.
  The independent census checks every labelled tournament through order six.
- [Sparse-budget and ribbon results](ribbon/REPORT.md): an explicit
  cutwidth-uniform sparse-back-arc kernel, all-order spectator obstruction,
  and a precisely finite repeat-block spectral exclusion. The qualitative
  bounded-budget sparsity mechanism is inherited, not new.

Further arithmetic and target-placement rounds continue. Earlier proof
candidates retained below their audited replacements are provenance, not
claims of a complete proof.

## Audited arithmetic and near-dyadic companions

- [Near-dyadic twins](dyadic/PROOF.md): for every t >= 0, set m=17+84t,
  P=(2^m-3672)/5096 and c=3672P+323f(P)+2539. Strong tournaments realize
  2^k-c and 2^k-c-2 for every k >= m+1, on k-3 vertices. Their source
  slope 2^m is optimal within the inherited positive-zero-tail dyadic route.
  At t=0 the deficits are 98861 and 98863. The independent audit proves
  the exact period, entry cost and limiting deficit ratio. Width remains two.
- [Source inverse](source/PROOF.md): a rigorous one-candidate-per-cylinder
  inverse for the inherited 2,732-label physical bank. Numerical running time
  is O(sqrt(N) log N) for the fixed bank, not polynomial in the bit length.
  This sharpens an inherited selection method but does not prove occupancy.
- [Two-clock rigidity](source/ROUND2/PROOF.md): an all-context obstruction for
  the specified equal-rate repeated-bit source model, with boundary and shift
  cases included. The positive centered-palindrome family is inherited;
  exact halfword coordinates and its square-root decoder are refinements.

Independent reports: [dyadic](dyadic_audit/AUDIT.md),
[source inverse](source_audit/REVIEW.md), and
[two-clock analysis](source_audit/ROUND2_REVIEW.md).

## Audited exact placement selector and graft boundary

[Placement proof](placement/PROOF.md) and [independent audit](placement_audit/AUDIT.md)
give a complete finite selector for this interval family with arbitrary
singleton prefix/suffix contexts. Any successful target N satisfies
5760*4096^k*2^m <= N-45, where m is the extra padding length. The search
returns actual graphs, but has not been shown to accept every target.

At the closing root, the interval is a pure-D four-state fiber. Iterating
arbitrary grafts there produces spacing 2H_*, where H_* is the product of
gadget Hamilton counts. For a nontrivial window, unit spacing forces every
gadget to be transitive and yields only existing zero-suffix placements.
The k=0 singleton boundary is explicitly excluded from that necessity claim.
Rerooting and other non-ordinal changes remain possible escape routes.

The main growing-width theorem is now routed as
[THM-4960](../../../01-canon/theorems/THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md).
Refinements to minimal loop length, quantitative diameter and alternative
root choices are still being audited; their local files keep their status.

The later [minimality and diameter audit](fiber_audit/MINIMALITY_PLACEMENT_AUDIT.md)
now passes. Twelve is the minimum unit-loop length for this specific moving
mechanism with physical cores of order at most nine. This is not a globally
shortest interval-constructor claim. For any fixed-(h,q) one-core bank with m
singleton positions, its equal-interface scalar diameter is at most
m(q+2h-2)+q+2h. The common baseline-selection problem remains open.
