# Session 231: moving-core intervals and the remaining placement problem

**Status: independently audited partial results. Chan–Pak completeness remains OPEN.**
All counts are raw directed Hamilton-path counts. This extended session pursued
complementary construction, physical descent, source arithmetic and ribbon
routes through several research and audit rounds. It did not complete the
proof that every positive odd count other than 7 and 21 occurs.

## Main positive result

[THM-4960](../../../01-canon/theorems/THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md)
gives genuinely unbounded consecutive-odd intervals using one moving core.

- **Common full interface:** for every k>=0, the k+1 counts D_k-2j,
  0<=j<=k, occur on strong tournaments with 12k+13 vertices, where
  D_k=(245051264*4096^k+1602510*k+1118023)/24843.
  The saturated deletion core has h=45 and gcd(h,q)=45. Its count is fixed,
  rather than a growing product of repeated nontrivial cores.
- **Smaller scalar witnesses:** deleting the protective flanks gives strong
  tournaments on 12k+9 vertices with counts Z_k-2j, where
  Z_k=(13617152*4096^k+1602510*k+717259)/24843.
  The first upper endpoints are 577, 2245227 and 9196067477.
  This smaller family does **not** have one full rooted interface: its last
  witness changes the root-ending coordinate by 14.

The proof is a finite exact twelve-bit word identity plus an all-order
context argument, not an extrapolation from finite k. The original 1188-bit
unit-loop construction remains valid but is less efficient.

[Full-interface proof](fibers/SHORT_INTERVALS.md) ·
[Independent audit](fiber_audit/SHORT12_THEOREM_AUDIT.md) ·
[Smaller-witness proof](rerooted/PROOF.md) ·
[Independent physical audit](rerooted_audit/AUDIT.md).

## What this does and does not solve

It closes the earlier missing **unbounded consecutive repair width** mechanism.
It does not place an interval over an arbitrary target. The displayed
baselines grow exponentially while widths grow linearly.

The final synthesis makes that gap precise:

1. [Exact selectors](placement/PROOF.md) enumerate every singleton-context
   placement of the protected family under explicit target-dependent bounds.
2. [Exterior rerooting](rerooted/PROOF.md) shows that common grafts confined
   to the closing root and its exterior suffix multiply interval spacing by
   the suffix Hamilton count. For a nontrivial unit window, this gives only
   existing transitive singleton padding. Protected-core grafts and other
   operations remain outside that obstruction.
3. [Absorption](bare_coverage/PROOF.md) proves that every binary-context
   placement of every bare moving interval is already a k=0 binary-context
   placement around the same eight-vertex core. Increasing k organizes nearby
   witnesses but adds no occupied targets to that complete context language.
4. The **complete target-bounded census of this fixed-core language** attains
   75,076 positive odd targets through 10^7 and 907,984 through 10^8, out of
   5,000,000 and 50,000,000 possible odds. These are construction-specific
   omissions, not missing tournament counts or an asymptotic density theorem.
   The largest internal gap has 854 absent odds between attained 987933 and
   989643. A separate transfer-tree implementation reproduces the entire union
   byte-for-byte, not just its size.

[Coverage audit](bare_coverage_audit/AUDIT.md). The exact remaining obligation
is a baseline selector, a more powerful core/interface family, or a different
operation that reaches every admissible target. Gcd-one, residue coverage,
large intervals and terminating inverses do not supply that obligation.

## Other audited outcomes

- **Physical descent:** [two-sided join inverse](descent/PROOF.md) and
  [equality refinement](descent/ROUND2.md). For a proper join,
  N=hk-ef+uv, h,k<=floor((2N-1)/3)<N and uv<=N-1. Sole-exit analysis sharpens
  the bound and classifies the equality cases. The inverse retains actual
  physical states; the full membership oracle is proved, not exhaustively run.
  [Audits](descent_audit/AUDIT.md), [round two](descent_audit/ROUND2_AUDIT.md).
- **Near-dyadic twins:** [an infinite optimally affordable source class](dyadic/PROOF.md).
  For m=17+84t, P=(2^m-3672)/5096 and c=3672P+323f(P)+2539, every k>=m+1
  has strong witnesses for 2^k-c and 2^k-c-2 on k-3 vertices. The source
  slope 2^m is minimal within the inherited positive-zero-tail dyadic route.
  This is width two, not coverage. [Audit](dyadic_audit/AUDIT.md).
- **Source arithmetic:** a [one-candidate-per-cylinder inverse](source/PROOF.md)
  for the inherited 2,732-label bank, and [arbitrary-context two-clock rigidity](source/ROUND2/PROOF.md).
  The positive centered-palindrome construction and cylinder strategy are
  inherited; their exact coordinate and decoder refinements do not establish
  occupancy. [Audits](source_audit/REVIEW.md), [round two](source_audit/ROUND2_REVIEW.md).
- **Sparse constructions:** [an explicit sparse-back-arc kernel](ribbon/REPORT.md),
  an all-order spectator obstruction, and precisely finite ribbon spectral
  exclusions. The qualitative fixed-budget sparsity mechanism is inherited.
  [Audit](ribbon_audit/AUDIT.md).
- **Carrier scope:** [independent minimality/diameter audit](fiber_audit/MINIMALITY_PLACEMENT_AUDIT.md).
  Twelve is the minimum unit-loop length in this particular moving mechanism
  with physical cores of order at most nine. For fixed h,q and m singleton
  positions, equal-interface scalar diameter is at most
  m(q+2h-2)+q+2h. Neither is a global minimum or all-core obstruction.

## Verification and reproducibility

The packet retains producer and independent proofs, exact witnesses, sources,
receipts, hostile controls, and normal/optimized comparisons. Major independent
checks include all 1,207,801 pairs of rooted labelled tournaments through order
five; every labelled tournament through order six for refined descent;
complete class/cut censuses through order nine; 35,567 tested block occurrences in the declared finite ribbon
universes; and full-subset DP for both 25-vertex named-port witnesses.

The order-nine fiber census has 191,536 classes and 193,541,720 nonempty erased
interfaces. Independent canonicalization/orbit mass and a separate sorted-record
replay certify its scope. The fixed-core coverage union has SHA-256
871cca7430a860219c3658c8672a18728cbecd01989f28e0ad3af4431bb99c77
in the documented 6,250,000-byte encoding.

[Reproduction instructions](REPRODUCE.md) separate compact checks, the large
physical subset checks, complete discovery censuses and generated artifacts.
Large class lists, union bitsets and binaries are regenerable and not committed.
Checker predicates remain active under Python -O; deliberate-failure controls
cover the repaired gates. MISTAKE-741 records the prepublication corrections.
The recorded preregistration is not claimed as independently timestamp-certified.

[Research synthesis and next obligations](SYNTHESIS.md).

## Inheritance and concurrent work

The session starts from main 9c48b9ac1e2a9c22b0814686c76e3358ad56051d and
builds on Sessions228–230, the older one-arc/F2 bank, Session217 cylinder
inversion, Session220's selector obstruction, and the earlier source-cost and
centered-fan packets. Counting identities, finite coverage cutoffs and classical
inputs are credited in each proof. No new unrestricted coverage cutoff is claimed.

Concurrent S87 and audited ternary-ribbon work were preserved through main
37f7d62b3f636807d8d67601edddc3cf62587804. Their base-three packing counts are
not raw Hamilton counts, and were not used as binary-completeness evidence.
The verified earlier checkpoints are 41329b268890302e9770212485777ada562547aa,
6092b657764d41f2428e088c5426368744ee9ccd,
b38a6a7d891eb699c9007f13925b27700b1cfb50 and
30127fb6e839b84eea01c78b30b1c526aae3694c.
