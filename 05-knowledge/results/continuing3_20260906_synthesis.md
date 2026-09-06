# Continuing synthesis: actual phase transport and three-atom rigidity

**PREVIOUS CHECKPOINT, 2026-09-06.** Current routing is the
[degree-eight geometry and retained-data synthesis](continuing4_20260906_synthesis.md). Four audited tracks extend the
[previous nine-result checkpoint](continuing1_20260906_synthesis.md).
This cycle inherited f0521b872 and integrated incoming work through
da2c452c3 before filing, then01c8b6887 and93a38412b during clean push rebases.
LRC(14), general actual Laurent two-rung
noncancellation, and full negativity in the anchored model remain **OPEN**.
No external priority or proof-assistant claim is made.

## What is proved, and what incoming work already settled

1. **PROVED relative to the proved/CITED real-root supplier, independently
   audited:** [Actual endpoint families](continuing3_20260906_laurent_endpoints.md),
   with [independent coefficient audit](continuing3_20260906_laurent_endpoints_audit.md).
   For h=6,7,8,9,10, N=6h+3, every integer g>=3h+2 coprime to N, and
   arbitrary nonzero complex alpha,beta,gamma,

       f(u)=alpha*u^(-N)+beta*u^(2g-N)+gamma*u^(3g-N)

   has first nonzero constant-term moment exactly g or 2g. Both values
   occur for every allowed g, and all h first-moment cancellation phases
   have nonzero doubled response. These are five specified families with
   unbounded positive endpoint, not all supports of width at most63.
   Incoming [endpoint39](long_frontier_sep06_endpoint39.md) independently
   settled h=6: all 258 coefficients and contents agree. Our extension
   beyond that incoming result is **endpoints45,51,57,63**, with 34 further
   characteristic polynomials and 2912 positive shifted coefficients.
   The full common certificate has 40 polynomials and 3170 coefficients.

2. **PROVED analytically, independently audited:**
   [Moment and entire-product characterization](continuing3_20260906_stability_near_minimizers.md),
   with [independent proof audit](continuing3_20260906_stability_near_minimizers_audit.md).
   Incoming [THM-4455 / three-atom minimizing-sequence rigidity](../../01-canon/theorems/THM-4455-three-atom-minimizing-sequence-rigidity.md)
   already proves the core unrestricted classification, using the sharp
   constant from [THM-4454 / signed-root duplication stability](../../01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md).
   Our derivation independently recovers it and its singular-boundary
   exclusions. The additions give the following equivalent conditions for
   arbitrary varying finite real lists with p1=p2=1 and positive energy:

   - The normalized stability quotient approaches its sharp constant K3.
   - The three largest positive roots tend to z=1/sqrt3, and all remaining
     square mass tends to zero.
   - M=sum r_i^2(r_i-z)^2 tends to zero.
   - p3 tends to z and p4 tends to1/3.
   - The actual products converge locally uniformly on the complex plane:
     product(1+r_i s) -> (1+zs)^3 exp((1-sqrt3)s).

   Writing d3^2=2-2z(a+b+c), the explicit bridge is
   `d3^2 <= (9600/49)M` when `M<5/6144`, and globally
   `M <= (1+z)^2 d3^2`. The sharp quotient-to-distance modulus remains open.
   Only the signed net dust first moment is fixed. A parameterized exact
   family allows arbitrary finite positive dust first-moment limits, or
   divergent positive and negative masses. The incoming theorem has the
   stronger explicit bound on the number of negative roots; retain it.

3. **PROVED analytically with complete exact certificates, independently
   audited:** [Smallest phase and effective tail](continuing3_20260906_laurent_finite_phase.md),
   with [independent eight-corner proof](continuing3_20260906_laurent_finite_phase_audit.md).
   In the specified nonnegative five-root model with sum13 and square
   sum59, the smallest positive original phase is unique and simple in
   (1/110,1/90). The complete carried response obeys sQ(-s)<-400 throughout
   [1/120,1/80], even on the containing coefficient box and without the
   two interlacers. The subsequent
   [two-residue floor](continuing3_20260906_residue_floor.md), with its
   [independent audit](continuing3_20260906_residue_floor_audit.md), proves
   **e4>161875/888583>9/50** when both prescribed weak interlacers hold.
   It combines a shifted C moment determinant, an ordinary D determinant,
   and the actual bounded-support inequality5e5<=M e4 with M<71/10.
   Retaining this floor and its slope linking e3 to e4 in the complete
   original-root envelope
   gives `Q(-s)/(s^7 e4^2)<-49000` at every original phase s>=**2500**.
   This improves the successive sufficient thresholds118163898523,75000 and4100;
   both earlier theorems remain valid. Any unresolved original phase now
   lies in **(1/80,2500)**, at most three per shape counted with multiplicity.
   This remains a continuous open problem. The later incoming
   [C-only residue proof](long_frontier_sep06_anchor.md) independently
   settles the same smallest branch under one interlacer. Our containing-box
   proof removes that hypothesis and strengthens the response margin. Its
   distinct residue envelope, described below, is useful for the other phases.

4. **PROVED elementary construction plus FINITE-EXACT controls,
   independently audited:** [Actual packet incompatibility](continuing3_20260906_lrc_packet_obstruction.md),
   with [full actual-entry audit](continuing3_20260906_lrc_packet_obstruction_audit.md).
   Set V=(325155600,429483600,498841200,2974571600,3212537328,5019589575)
   and U=(1,43,67,79,97,109,127). Every tV union U with
   `412164<=t<=25880486141010` and `gcd(t,305613336829)=1` is an actual
   primitive bounded6+7 equality entry retaining all current hereditary
   joint profiles. The primitive V first safe denominator is29. At
   t=412164 the physical first denominator is29, while at t=412183 it is31:
   the same two primitive packets have zero compatibility at29 in the
   latter row. The located scale residue is information that the component
   shapes and all these profiles omit. This refutes the stated small-clock
   and first-clock compatibility forcing, not adaptive choice of either
   component: U has a denominator2 phase. The inherited larger-unit theorem
   already proves this entire class safe, so this is a structural stopping
   object rather than a new LRC closure.

## Mechanisms recovered and the information each map retains

The inheritance pass began with the actual decoder's pair-gcd gates, the
complete carried factorial row and endpoint33, the sharp stability regional
envelopes, and the effective quadratic anchor. Hostile controls were the
actual balanced decoder residual, the failed derivative-window cone, the
one- and two-atom quotient boundaries, and the linearly anchored escaping
phase. The least-used sidecars were located CRT multipliers, principal-minor
degree bounds, signed dust's first moment, and the positive-term inventory
after eliminating the original zero.

| Source and target | Map and preserved predicate | Lost information, required sidecar, decisive test |
|---|---|---|
| Complete actual Laurent fibres -> finite quotient operator | Divide the full doubled row by the monic first row, retaining the inverse carry; the eigenvalues are same-root responses | A quotient can lose the raw complex scalar, whose nonvanishing is retained separately. Characteristic positivity excludes nonnegative real responses only after the real-root supplier applies. Full integer-Berkowitz reconstruction tests the map. |
| Signed root lists -> square-weighted measure and entire product | M detects concentration near0 and z; integer multiplicity forces exactly three macroscopic roots | Square mass alone forgets signed linear dust. Restore p1=1 to identify the exponential factor; a mixed-sign dust family tests the boundary. |
| Root/interlacer geometry -> coefficient box or eliminated response | The box preserves every response coefficient; original-root elimination preserves P(-s)=0 | The box loses beta real-rootedness and interlacing. Its sign proof works on the certified interval; the tail retains e4>1/100. An exact Newton/first-rooted hostile tests the failed wider relaxation. |
| Primitive cofactor cover -> actual scaled decoder entry | Scaling preserves primitive shapes and internal weighted graph, while uniform inequalities exclude every bounded crossing | Body gcd profiles omit both small internal divisor owners and the packet multiplier. The two actual scales give the cheapest direct test; the exact residue-labelled overlap is the sidecar. |

Two distinct certificates establish the smallest-phase box bound: all270
tensor Bernstein coefficients are below-400, and separate convexity reduces
the same continuous box to eight univariate corner polynomials, each settled
by exact Sturm counts. Agreement here is an independent proof mechanism,
not just a second implementation of the same basis change.

The exact relaxation hostile has (e3,e4,e5)=(104,50,37435088/3898125).
All four displayed Newton inequalities hold strictly, and the original
first polynomial has four simple positive phases. Yet Q(-15/2)>0. Its beta
polynomial has only one real root, so it is outside the model. First-row
real-rootedness does not reconstruct the beta-root geometry.

Incoming01c8b6887 supplies two concrete connections. The positive C/B residue
measure has moments(1,1,3,e3/3-16,16e3/3-373-4e4/7). Its moment matrix gives
`75<=e3<=135` and `e4<=(7/72)(e3-75)(135-e3)<=175/2`. This necessary
semialgebraic envelope is stronger than a box in that direction, but it is
not asserted equivalent to the root/interlacer predicate. Pursuing that
sidecar produced the stronger floor above: with delta=e3-75, C's shifted
third determinant and5e5<=M e4 force e4>(8750/8241)delta; D's ordinary
third determinant then gives the displayed rational floor. No previous
e4 floor is a dependency of this new argument.

The rational positive-response hostile passes both ordinary third and
shifted second moment-matrix tests for C and D. It nevertheless violates
mu4<=(71/10)mu3 by2159/105 for C and1091/42 for D. Thus restoring the
anchor-derived support endpoint already rejects it at moment order four;
C's shifted third determinant rejects it as well. Low moment positivity
alone had lost precisely that support coordinate. Higher bounded-support
residue constraints remain a concrete route into the surviving interval.

The incoming [fifth-clock family](long_frontier_sep06_phase.md), with its
[independent audit](long_frontier_sep06_phase_audit.md), gives384061 bounded
scales of an actual unitless5+8 entry, retaining all current profiles and
failing the listed native, endpoint, dual and whole-arc scalar gates. A
four-case signed affine phase word proves clearance greater than1/6 for
every member. At its actual hostile scale4912753 the nearest lift fails,
while a farther lift succeeds. This provides a positive counterpart to
our prescribed-clock incompatibility: the lost offset must be retained
and may need a different base point. The two results do not themselves
force an adaptive word for every actual decoder entry.

The final incoming [third-session board](third_20260906_board.md) changes
three next questions. Its audited [grid overlap reduction](third_20260906_grid.md)
proves weak safety for actual balanced entries with six-component scale
t>=97097. Its [rooted gcd-product theorem](third_20260906_decoder.md) also
closes every balanced equality entry with maxU<=28. Thus the balanced
residual search can retain t<97097 and maxU>=29. Our packet examples have
t>=412164: they remain valid obstructions to the stated unrestricted
forcing, but do not refute a forcing theorem restricted to this newly
smaller residual domain. This is a new hypothesis to retain in future probes.

The incoming [all-height carried factors](third_20260906_laurent.md) apply
to the exact same normalized quotient used by our five-family certificate.
Negative-integer parameter degenerations force explicit polynomial divisors
of every characteristic coefficient. For h=6 this compresses258 positive
coefficients to208 residual coefficients. The all-height divisibility does
not establish residual positivity. Deflating our h=7,...,10 certificates
is now a concrete first test before attempting a uniform recurrence.

The incoming [complete mixed(m,2,1) Smith law](third_20260906_smith.md)
extends the older(2,2,1) theorem using our complete-bank projective transport.
Clearing the growing multiplicity bank leaves exactly three residual rows,
and three generators control the second determinantal ideal for every m.
The next rank change should retain an additional intermediate ideal, rather
than reopen this now-classified multiplicity family.

The same incoming commit repairs the microscopic addendum of
[THM-739 / exact Bernoulli overlap](../../01-canon/theorems/THM-739-pairwise-coprime-bad-overlap-exact-bernoulli-closed-form.md):
tooth intersections require containment clipping and actual window placement.
The main full-circle identity survives. Read that repair before transporting
any overlap credit; different credits can only be summed with a valid
pointwise multiplicity bound.

## Six live concepts and the next decisive work

| Lane | Present object | Next test or stopping reason |
|---|---|---|
| Anchor: actual LRC equality | Incoming grid and rooted-product bounds leave t<97097,maxU>=29 in the balanced branch | Test adaptive signed-offset/farther-lift rules within this residual domain. Our existing prescribed-clock hostiles are outside it. A generic CRT compiler would repeat an existing theorem. |
| Niche: all-h carried response | Positive characteristic certificates through h=10 and incoming all-height divisors | Test the deflated h=7,...,10 coefficients, then seek a positivity-preserving operation in h or direct same-root proof. Divisibility alone does not provide positivity. |
| Niche: anchored finite phase | Smallest branch and tail closed; interval(1/80,2500) remains after the two-residue floor | Retain beta roots and both interlacers; test higher C/D residue moments. Challenge every necessary-inequality relaxation with the rational positive-response hostile. |
| Wildcard: sharp three-atom rigidity | Moment residual plus complete entire-product limit | Seek an effective quotient-to-M estimate near the three-atom boundary, retaining separate dust and imbalance scales. Do not reopen the proved scalar optimum or the core iff. |
| Projective higher jets | Incoming full(m,2,1) classification uses fixed residual rank three and complete-bank covariance | Explore the next residual rank only with every required intermediate ideal; read the incoming second round first. No additional Smith theorem is claimed by this cycle. |
| No-three-in-line and native operators | The recovered Guy-Kelly grid thread and exact Boolean sectors remain available | Require an explicit map preserving collinearity or the native operator. Analogy with packet counts, moments, or root signs currently supplies no such map, so no asymptotic transfer is claimed. |

The results change the board in specific ways: the LRC example rules out
profile-only clock selection; the endpoint proof bypasses a failed positive
response representation; the product theorem restores a coordinate lost by
square-mass compactification; and the finite-phase hostile identifies the
root predicate that a broader sign certificate must keep. None transfers
automatically to another lane. No new META-PATTERNS card is promoted solely
from these analogies.

## Reproduction and checkpoint discipline

The [manifest](continuing3_20260906_manifest.json) pins all33 proof, program,
output and JSON artifacts, recording pre-filing report identities where
status, credit or portable links changed. All source/output/certificate
bytes remain exactly those audited. The ten producer/referee programs
perform42535 always-active exact gates per combined normal run, and every
optimized output matches the corresponding normal and frozen output.
Universal claims rest on the audited analytic arguments or full polynomial
identities, not extrapolation from their finite controls.

The original shared dirty checkout is untouched. Continue in the isolated
worktree, read incoming mathematics at clean checkpoints, and explicitly
stage, commit and push audited changes before yielding. The active hourly
follow-up stays quiet without a substantive result, correction, completion,
failure or required user action.
