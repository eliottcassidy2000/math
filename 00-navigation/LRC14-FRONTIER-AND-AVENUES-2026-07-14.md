# LRC(14) Frontier and Research-Avenue Atlas — 2026-07-14

**Status:** current frontier audit, historical synthesis, and pull-card backlog.
**Owner:** codex-2026-07-14-S2.
**New mathematics in this session:** THM-760 and HYP-6785.
**Rebase note:** integrated the incoming S295--S297 Lean completion, THM-759
ratio bound, HYP-6800/HYP-6805, and proof-carrier automaton viewpoint.  The
sheet-dodge claim was renumbered from THM-759 to THM-760 after the collision.

## Executive assessment

For a set S of thirteen distinct positive speeds, write

M(S) = max over t in R/Z of min over v in S of ||vt||.

LRC(14) asks for M(S) at least 1/14.  The repository has made substantial
progress, but the honest frontier is narrower and less finite than several
recent theorem headlines suggest.

1. **The non-covering branch is settled.**  THM-366 gives the rational witness
   t=1/q whenever S omits a multiple of some q in {2,...,14}.  A strict
   counterexample must therefore be divisor-complete, or covering.

2. **The tight covering branch with at most three speeds above 14 is settled.**
   Such a family has at least ten speeds in {1,...,14}, so THM-738 applies.
   This includes the AP/Goddyn-Wong boundary atoms, the deep-well family, and
   every currently known low-margin family.

3. **The capped-envelope theorem is a genuine analytic advance, but its cutoff
   is per core.**  THM-755 proves the required discrepancy inequality above
   v*(P)=r_P/(pi |G'_P|).  THM-756 closes every bottom-core band and samples the
   assembled band successfully.  The paper proof is elementary.  After opus
   S295--S296, the geometric-to-spectral Lean chain is also complete:
   LRCClosedBudget.lean has 47 declarations, zero sorries, including the
   pair-overlap and family-assembly identities.

4. **The remaining far-count branch is not yet one uniformly bounded finite
   computation.**  THM-758 correctly proves its Claim A and gives strong
   evidence and per-core decision procedures for Claim B.  Its inference
   v*(P) is globally about 500 is false: HYP-6780 proves
   v*(cP)=c v*(P).  Therefore a finite band for each fixed P is not a finite
   raw-height band over all P.  The 120/120 small-good-period observation
   `q in [15,25]` was a useful sample, but THM-762/MISTAKE-143 now refute its
   proposed uniform extension with exact coherent and gcd-incoherent residuals.

5. **This session closes the largest elementary scale ray.**  THM-760 proves
   that for c at least 2 and gcd(c,w)=1,

   M(cP union {w}) >= min(M(P), 1/2 - 1/(2c)).

   Hence any primitive thirteen-speed family in which twelve speeds have a
   common divisor is already at least 1/13.  This closes the codimension-one
   scaled-core ray of HYP-6780 without a height cutoff, AP assumption, or
   enumeration.  The genuine scale residual has at least two exceptional
   residue classes, or has no twelve-speed common-factor core.

6. **The covering minimum 14/183 is a compelling candidate, not a proved
   global theorem.**  The family {1,...,12,182} has exact M=14/183, and the
   interval-core, tight, and many finite lanes point to its uniqueness.
   THM-724's body explicitly leaves a near-tight non-dilated residual
   empirical, despite its title; THM-726's multi-killer conclusion likewise
   relies on exact interval-core enumeration plus soft global routing.  These
   results are more than enough motivation, but they should not be cited as a
   completed global covering-min proof.

7. **The incoming tight-instance ratio theorem is proved but not
   closure-critical.**  THM-759 gives the elementary bound
   a_n <= n a_{n-1} for a tight n-set.  HYP-6800 uses it to localize the
   remaining tight-set classification to the non-extremal-core, Goddyn-Wong
   branch.  That sharpens rigidity and supplies finite ranges; it does not
   convert the unexecuted scale-quotient covering band into a theorem.

The shortest honest statement of the frontier is therefore:

> Prove LRC(14) for primitive covering families with at least four speeds
> above 14 after quotienting scale, excluding the one-exception common-factor
> lane closed by THM-760.  Equivalently, prove a uniform recursive
> scale-quotient certificate—small good period, empty pair-sum blocker edge,
> safe peel, or a finite normalized residue certificate.

There is no known counterexample to the conjecture, and the remaining
families are empirically loose.  What is missing is a uniform theorem that
turns that looseness into a certificate without silently bounding raw scale.

## Status ledger: what may safely be used

| Layer | Mathematical content | Status to cite |
|---|---|---|
| q-divisibility gate | omitted q in {2,...,14} gives t=1/q | proved, THM-366 |
| exact maximizer search | an optimum occurs at a pair-sum ruler q=v_i+v_j, retaining all multipliers | proved, THM-668 pair-sum ruler |
| near-AP ten-body tile | at least ten speeds in {1,...,14} implies LRC(14) | proved, THM-738 |
| safe peel | a runner safe at a core optimum can be removed without changing M | proved; irreducible-family tiling is empirical, THM-753 |
| capped envelope | discrepancy closes a peel for v>v*(P) | proved, THM-755 |
| bottom H-bands | all 91 bottom cores and their bands | exact finite proof, THM-756 |
| far-count Claim A | at most three far speeds | proved through THM-738, THM-758 |
| far-count Claim B | at least four far speeds | per-core decidable and sampled; not globally raw-bounded |
| q<=25 good periods | false uniformly; observed on 120 sampled rows | THM-762/764, MISTAKE-143; S312 sample only |
| scale behavior | |G'| fixed, component count and v* scale by c | proved, HYP-6780 note |
| one exceptional residue | cP plus coprime w | proved, THM-760 |
| tight-set ratio | a_n at most n times a_{n-1} | proved, THM-759; sporadic classification branch remains |
| tight-set height | primitive tight 12-tuples have sum at most 78^11 | proved, THM-763; finite is not empty |
| tight hereditary primitivity | every leave-one-out core of a primitive tight 12-set is primitive | proved, THM-765 |
| tight first-window geometry | a_11/a_1 at least 72/7; eleven sub-12 tooth cones | proved, THM-766 |
| tight zero-owner split | unique-maximum 13-multiple impossible; shallow packets are full residues; deep packets obey sheet capacity | proved, THM-768/769 |
| bounded shallow rigidity | through max 168, the only primitive tight full-residue row is {1,...,12} | finite-exact, THM-770 |
| global covering minimum | M at least 14/183, equality only at deep well | conjectural outside routed/finite lanes |

## How the repository's understanding sharpened

The history is best read recursively.  At each stage a quotient looked nearly
complete; an impersonating family survived because the quotient forgot one
coordinate; that forgotten coordinate became the next proof object.

### 1. From runner tournaments to labelled walls

Early tournament and observer-source models supplied a finite language for
dominance, wall switches, cycles, and Hamiltonian paths.  Their lasting
contribution is methodological: declare the vertices, pairwise observable,
switch, and tie path.  Their failure was magnitude blindness.  AP, GW, and
loose residue liars can share a raw tournament shadow.  The repair was to label
walls with exact endpoint coordinates and owners.

### 2. From endpoint counts to divisibility and pair-sum rulers

Endpoint/fold/pincer work showed that a safe interval is born at a collision
of tent-function walls.  THM-366 distilled one side into the small-denominator
divisibility sieve.  THM-668 later distilled the full maximization geometry:
one only needs endogenous rulers q=v_i+v_j, but every multiplier—including
nonunits—must remain.  Raw endpoint counts were replaced by labelled proof
obligations.

### 3. From positive measure to closed-threshold topology

Density and singular-series methods explained why generic, Sidon-like sets
are easy: additive relations suppress open safe mass.  They also exposed a
fundamental boundary correction.  Positive open mass proves M>1/14, but zero
open mass does not refute M>=1/14.  AP and GW live exactly on the closed
boundary.  Haar/Baire, Cech, and endpoint-owner packets were introduced to
carry open-versus-closed debt.

### 4. From covering arithmetic to a tight boundary stratum

The q-divisibility gate separated the problem into non-covering rational
witnesses and primitive covering sets.  AP/GW ceased to be candidate strict
counterexamples and became equality atoms at the boundary of the
non-covering side.  Covering is not invariant under dilation; this later made
the primitive/imprimitive distinction load-bearing.

### 5. From AP resemblance to state-lift labels

Residues alone could not separate AP, GW, 12-to-26 liars, and 12-to-36 K33
near-misses.  C27, unital, K33, petal, and state-lift packets kept the hidden
transfer data.  The useful conclusion was not that one finite geometry is the
object, but that any quotient must declare which LRC predicate and which
owner/lift information it preserves.

### 6. From a fixed denominator ladder to a dynamic blocker system

THM-566 killed every fixed external denominator window.  The twist-ladder
work HYP-2443/HYP-2972 repaired the idea dynamically: a failed rational
witness is not wasted; it identifies blockers and feeds a hypergraph or dual
certificate.  HYP-2978/HYP-2979 added exact-period and Ramanujan guardrails.
The missing move was to make the denominator set endogenous.  HYP-6785 does
that by attaching the twist blocker idea to every THM-668 pair-sum event.

### 7. From scalar density to moments, Fourier cancellation, and packet sheaves

Bonferroni, covariance, gK8, Fourier, Bernoulli, Dedekind, Fejer, and relation
lattice work produced real inequalities and powerful finite certificates.
Repeated failures of scalar monotonicity led HYP-2976/HYP-2980 to a
proof-object sheaf: exact scale, boundary status, owners, residues, lifts,
moments, and dual certificates travel together.  The July S312 computation
reinforces this lesson.  Low-order Bonferroni bounds become negative, while
absolute Fourier relation sums are enormous; the good-set value is made by
signed cancellation.

### 8. From global height bounds to scale-normal recursion

Near-AP, slow-fast, decorrelation, density-floor, safe-peel, and capped-envelope
work made the endgame look finite.  HYP-6780 found the missing coordinate:
the good-set component count scales with a circle cover, so the cutoff v*
scales too.  Raw height was the wrong compactness coordinate.  THM-760 then
turns the same cover into a positive tool—its witness sheets dodge one
exceptional runner.  The next proof must operate on normalized shape plus
residue data, not raw maximum speed.

This progression can be summarized as:

~~~text
runner set
  -> labelled endpoint walls
  -> q/Farey and pair-sum rulers
  -> open/closed boundary packet
  -> owner/residue/state-lift packet
  -> dynamic blocker hypergraph and dual weights
  -> scale-normal sheeted witness-obstruction complex.
~~~

The post-S296 Lean status removes a former implementation caveat from this
history: THM-755's overlap/Bernoulli bridge and family assembly are now
machine-checked.  What remains is not that analytic identity, but uniform
coverage of the normalized residual states and composition of their terminal
certificates.

## Viewpoint atlas: what each lens sees and forgets

| Lens | Natural object | Preserves | Destroys or hides | Best current use |
|---|---|---|---|---|
| Torus/geodesic | orbit tS in (R/Z)^13 | exact simultaneous phases | arithmetic provenance of walls | conceptual definition and compactness |
| Forbidden-arc cover | union of danger intervals | exact threshold topology | which collision creates an optimum | covering and Cech language |
| q-divisibility | multiples of 2,...,14 | cheap rational witnesses | witnesses at endogenous large q | first strict-counterexample gate |
| Farey/continued fractions | rational cells and neighbors | scale and boundary adjacency | simultaneous blocker identity unless labelled | local witness refinement |
| Pair-sum rulers | q=v_i+v_j, time p/q | complete exact maximizer set | global scale relations between rulers | exact M and HYP-6785 |
| Endpoint-owner graph | wall collisions with owners | protected runners and boundary debt | interior safe mass | equality/near-equality analysis |
| Haar/Baire/Cech | regular-open safe set plus boundary | open versus closed threshold | arithmetic labels | AP/GW boundary handling |
| Singular series | additive relation lattice | signed open-mass mechanism | closed equality and local owners | loose/generic intuition |
| Bonferroni/moments | intersection moments of danger arcs | finite lower certificates | high-order signed cancellation | near-AP finite trees |
| Fourier/Bernoulli | coefficients and autocorrelation | exact discrepancy and scale tails | combinatorial owner identity | THM-731/732/755 |
| Fejer/Toeplitz | positive trigonometric minorants | controlled spectral truncation | exact endpoint exceptions | possible signed-cancellation certificate |
| Ramanujan projectors | exact period/divisor channels | nonunit period information | cross-period interactions | denominator guardrail |
| Additive energy | repeated sums/differences | AP pressure and relation depth | which relation helps at threshold | structure-versus-spread split |
| BSG/Freiman | approximate additive structure | high-energy inverse structure | sharp 1/14 constants | coarse rigidity router |
| Schur triples/E3 | low-order additive incidences | third-order obstruction | higher relation lattice | sidecar, not scalar proof |
| AP/GW packet | tight boundary family | exact equality grammar | loose families | terminal equality atoms |
| C27/unital/K33 | state-lift incidence | hidden residue transfers | exact M without Farey labels | near-boundary packet routing |
| Twist ladder | ordered rational attempts | blocker evolution | completeness if denominators fixed | recursive primal/dual search |
| Blocker hypergraph | proof obligations versus unsafe runners | simultaneous obstruction coverage | continuous margins unless weighted | exact finite predicate at fixed S |
| LP/Farkas dual | weights on blocker obligations | obstruction pressure certificate | geometric intuition | local-to-global theorem target |
| Leader ledger | ownership changes along t | conservation and transitions | scale compactness | binders and irreducibility |
| Safe peel | deletion at a core optimum | exact M on safe branch | global classification of irreducibles | reduction to lower LRC |
| Slow-fast balance | perturb a resonant killer | quantitative near-witness | multiple coupled killers | deep-well lane |
| Ostrowski/three-gap | rotation gap words | low-denominator ordering | high-dimensional blocker coupling | interval-core and good-period work |
| Eisenstein/Phi6 | 13/14 doorway arithmetic | deep-well 183=13*14+1 structure | general loose families | extremal arithmetic |
| Two-scale/decorrelation | compact core plus far speeds | asymptotic independence | uniformity under core dilation | loose tail estimates |
| Capped envelope | measure cap plus jump envelope | per-core tail with explicit v* | global raw compactness | far-element peel |
| Scale sheets | fibers t=(t0+k)/c | whole core margin exactly | endpoint labels across several exceptions | THM-760 and next extension |
| Affine-slope suspension | V(c)=cP+R as slope-c fibers of a two-torus function | scale, offsets, and exact fiber nonemptiness | finite classification without a descent | HYP-6815 chart for multi-exception rays |
| Transverse tooth refinement | a high frequency cuts a fixed safe component | scale-free component growth and wall ownership | boundedness from maximal divisor scale | THM-772; peel-relative splice |
| Tropical/normal fan | active minima and wall cells | combinatorial type of optimum | metric clearance | finite chamber stratification |
| Matroid/circuit | minimal dependent proof obligations | irreducible obstruction support | phase and margin | candidate blocker-pressure abstraction |
| Tournament analysis | pairwise quotient of chosen objects | dominance fingerprints | higher-order intersections | diagnostics and route comparison |
| Formal certificates | rational witness/interval/LP proof terms | auditability and exactness | discovery flexibility | final closure and regression tests |

No single row in this table is the underlying object.  The useful object is
the smallest common refinement needed by the exact predicate.

## Proposed underlying object: the scale-normal endogenous witness-obstruction complex

HYP-6785 isolates an exact finite object for a fixed family S and threshold
delta.

For every unordered pair of runner indices i,j, set q=v_i+v_j.  For every
multiplier p modulo q, including nonunits:

1. retain the event t=p/q;
2. retain the protected endpoint set {i,j};
3. discard the event if a protected endpoint is already unsafe;
4. otherwise attach the hyperedge B(i,j,p) of all other unsafe runners.

By THM-668,

M(S) >= delta if and only if at least one retained hyperedge is empty.

This combines the exact endogenous witness list with the older twist-ladder
blocker hypergraph.  It preserves the LRC predicate exactly, including
simultaneous blocking and protected endpoints.  It also explains why a runner
tournament is inadequate: pairwise responsibility does not determine whether
one obligation is blocked simultaneously by zero, one, or several runners.

The object still needs three additional fibers for a uniform theorem:

- **scale fiber:** identify circle-cover copies cP while retaining exceptional
  residues and gcd strata;
- **margin/boundary fiber:** distinguish strict safe edges from equality
  atoms and record clearance;
- **reduction fiber:** record safe peels, sheet dodges, and normalized
  descendants.
- **fragmentation/peel fiber:** record `r_P/(v|G'_P|)`, the divisor-support
  profile, and the endpoint owners creating or merging safe components.

The incoming proof-carrier synthesis gives the same object an operational
form.  A state of the **scale-quotient peel certificate automaton** is

~~~text
normalized core shape
  + core scale and exceptional residue/offset word
  + peel-relative fragmentation load and wall-owner word
  + far count and peel/irreducibility state
  + blocker debt
  + terminal-certificate payload.
~~~

Its terminal labels are q-witness, near-AP/THM-738, safe peel, aligned tooth,
capped envelope, exact discrepancy/direct interval, good period, sheet dodge,
and named tight corner.  A Nerode-style quotient may merge two states only if
every hidden completion has the same terminal certificate.  This is the
controlled-forgetting rule missing from raw speed, runner-tournament, and
raw-height classifications.

HYP-6815 supplies a natural chart for the same state space.  An affine family
V(c)=cP+R is the slope-c geodesic fiber of
Phi(u,t)=min_i ||p_i u+r_i t|| on a two-torus.  The R=0 cylinder is exactly the
dilation law of HYP-6780; nonzero R is transverse residue holonomy and is the
right place to encode several exceptional offsets.  This representation is
exact, but it needs the blocker/certificate labels above before it becomes a
finite theorem.

The exactly-`f=4` family
`{1,...,9,15,110,N,1092N}` shows why the new fragmentation field is necessary.
It has no divisor packet of size seven, but the top-peeled core's component
count is unbounded because the `N`-runner inserts new teeth into a fixed safe
interval. Thus common scale and transverse frequency are independent
noncompact directions in the same four-far chart. THM-772 closes the entire
prime family by the capped peel, whose ratio normalizes component count against the named
peel; `c*` alone cannot support the automaton quotient.

THM-773 extracts the general mechanism at the state level. If the old good set
has mass `mu` and `r_B` components, a new frequency `N` leaves mass at least
`6mu/7-2r_B/(7N)` and creates at most `N+r_B` components. A proportional peel
`aN` is therefore eventually terminal under an explicit rational inequality.
One marked interval gives the simpler `L,sum(B)` corollary used by THM-772. The
genuinely open transverse face must include collapsing safe mass relative to
component and peel rates.

The incoming exact endpoint-sidecar audit independently confirms the payload
rule.  Runner and unweighted endpoint tournaments do not preserve covering,
exact M, capped-envelope status, or Bernoulli discrepancy.  The tested
formula-facing sidecars are respectively the divisor mask, projective cap
ratio v|G'|/r, signed endpoint phases, and exact peak witness.  These sidecars
should be fields of the automaton, not reconstructed from tournament scores.

A candidate counterexample is therefore not just a speed set.  It is a
projective arithmetic shape carrying a sheeted, weighted obstruction
hypergraph.  The likely proof is a local-to-global statement saying that such
a primitive covering complex cannot hit every protected obligation unless it
falls into a classified tight packet; outside that packet, an empty edge,
sheet, or small good period appears.

### First computational fingerprints

The companion exact atlas verifies THM-760 on 80 random 12-core instances and
computes blocker complexes for eight representative families.

- AP and GW have M=1/14 and repeated pair-sum rulers of multiplicity 5 or 6.
- The deep well has M=14/183 and multiplicity 6.
- The near-dilate c=26 rays have M=1/13 and multiplicity 6; THM-760 supplies
  their analytic witness.
- A spread family has pair-sum multiplicity only 2, thousands of empty
  obligations, and M=406/1669, about 0.243.
- Runner tournaments range from transitive to cyclic and do not track M:
  their value is diagnostic, while the hypergraph retains the proof predicate.

This is not yet an inverse theorem, but it suggests a precise split worth
testing:

> High pair-sum multiplicity should force an additive/sheet/ruler
> classification; low multiplicity should force blocker sparsity,
> decorrelation, or a bounded good period.

That statement is substantially sharper than “structured versus random”
because both sides have exact observables and certificate targets.

## New theorem: coprime exceptional-runner sheet dodge

THM-760 is elementary but removes an unbounded family from the frontier.
Choose a core optimum t0 for P.  At each lifted time

t_k=(t0+k)/c, for 0<=k<c,

every runner cp has exactly phase p t0, so all core clearances are unchanged.
If gcd(c,w)=1, the exceptional phases w t_k form a translate of the complete
c-grid.  One is within 1/(2c) of 1/2 and hence has clearance at least
1/2-1/(2c).  Therefore

M(cP union {w}) >= min(M(P),1/2-1/(2c)).

For a twelve-speed core, settled LRC(13) gives M(P)>=1/13, and the sheet
clearance is at least 1/4.  Thus the thirteen-speed family has M>=1/13.

The important conceptual change is the vertex choice: the natural vertices
for this lane are witness sheets, not runners, arcs, or pair-sum rulers.

## Ranked pull cards for future agents

Each card names a theorem-facing deliverable and a guardrail.  Priority A
cards attack the exact current residual.  Priority B cards build the inverse
structure.  Priority C cards revive niche threads only where they attach to
the retained object.

### Priority A — direct frontier attacks

**A1. Multi-exception sheet union bound.**  For
cP union {w_1,...,w_r}, count exactly how many sheets k are bad for each
exception at delta=1/14.  If every w_a is coprime to c, a translated c-grid
gives at most floor(c/7)+1 bad sheets per exception.  Determine the exact
small-c corrections and the largest r for which the union cannot cover all
sheets.  Deliver a theorem extending THM-760, not a random census.  Guardrail:
the exceptional phases share the same k, so they are coupled.

**A2. Gcd-stratified sheet dodge.**  Remove gcd(c,w)=1.  The exceptional
runner visits c/gcd(c,w) sheets with multiplicity.  Classify which gcd
profiles force a safe sheet and which descend to a smaller scale.  Target a
recursive lemma that either finds a witness or lowers c.

**A3. CRT sheet packet for two to six exceptions.**  Attach residue classes
of w_a modulo c to the sheet blocker hypergraph.  Search for a Hall,
covering-number, or circular-arc condition that prevents the union of their
bad-sheet sets from covering Z/cZ.  Stop if the condition is merely an
unlabelled union bound; record the exact residue counterexample instead.

**A4. Uniform good-period theorem or adversarial refutation.**  Reproduce the
S312 q in [15,25] observation with exact arithmetic on a normalized
scale-quotient generator, deliberately including cP plus several exceptional
residues and large c.  Either prove a q bound on the residual or produce the
first normalized family requiring q>25.  Do not infer a theorem from 120
floating-point samples.  This card is now owned by codex-2026-07-14-S3 in
05-knowledge/hypotheses/HYP-6820-q25-and-n12-uniformity-audit.md (renumbered
from HYP-6810 by opus-S299 after the collision with opus-S298's earlier
assembly-paper stub, per the first-pusher protocol).

**A5. Canonical scale quotient.**  Define a terminating normal form for a
primitive covering family: maximal common-factor cores, exceptional residue
word, normalized gap vector, and a descent measure.  Prove that each
nontrivial scale ray either invokes THM-760/A1-A3 or decreases the measure.
The deliverable is the missing replacement for the false raw bound.

**A6. Exact general band generator.**  Enumerate normalized shapes and gcd
profiles rather than speeds up to W.  Every emitted case must carry a proof
that all raw dilates are represented.  Route each to THM-755, THM-760, a good
period, safe peel, or an explicit unresolved packet.

**A7. Blocker-pressure theorem.**  In HYP-6785, prove that a primitive
covering complex at delta=1/14 cannot have every protected obligation hit
unless its repeated-ruler profile is near AP/GW.  Start with private blockers:
count obligations blocked by exactly one runner and double-count them against
protected endpoint pairs.

**A8. Weighted blocker dual.**  Form the set-cover LP whose rows are
pair-sum obligations and columns are runners, forbidding protected columns.
Search for a dual weighting derived from ruler multiplicity or Fourier
clearance.  A successful dual must imply an empty edge or classify equality;
an infeasible dual should output a small obstruction circuit for study.

**A9. Additive multiplicity dichotomy.**  Let r(q) count pairs with sum q.
Prove a quantitative alternative: max r(q) large implies a long symmetric/AP
core; max r(q) small implies enough independent obligations or small-period
witnesses.  Calibrate on AP, GW, deep well, near-dilates, and spread families.
Guardrail: raw additive energy alone was repeatedly nonmonotone.

**A10. Irreducible-to-blocker bridge.**  Relate THM-753 irreducibility to
HYP-6785: if every runner binds at every complement optimum, what minimum
degree or circuit condition follows in the blocker complex?  This could turn
the empirical “irreducibles are tiled” statement into a structural theorem.

**A11. Correct the far-count assembly.**  Rewrite THM-758 Claim B using a
scale-normal recursion.  State separately the proved capped-envelope tail,
completed bottom bands, sampled general bands, and q<=25 observation.  The
success condition is a well-founded cover of all primitive shapes, not a
maximum-speed estimate.  Coordinate prose with HYP-6810 (opus-S298/S299, which owns the
assembly-paper write-up; draft at 04-paper/lrc14-assembly.tex); this card owns
the missing uniform justification. THM-761 (opus-S299) now closes the
multi-exception coprime lane of that justification for r <= 6 exceptions above
explicit scale thresholds; see cards A1-A3.

**A12. Formalize THM-760 and compose the completed THM-755 bridge.**  Lean
should express the circle-cover sheet permutation and the 1/(2c) grid
approximation.  THM-755's circular-overlap and family-assembly chain is now
complete; reuse it as a terminal route rather than reopening it.

**A13. Affine-slope and transverse-tooth classification.** Execute HYP-6815
on the multi-exception residual `V(c)=cP+R`, but include scale-free paths such
as `{1,...,9,15,110,N,1092N}`. Stratify the two-torus strip arrangement,
attach HYP-6785 blocker edges and the audited divisor/cap/endpoint/peak
sidecars, and retain peel-relative fragmentation. Prove that each fiber either
has a threshold point or descends to a smaller normalized state. Guardrail:
the exact suspension and cap-normalized state are representations until a
finite or well-founded classification is proved.

### Priority B — structural inverse theorems

**B1. Pair-sum-ruler inverse theorem.**  Classify thirteen-point sets with
many pairs on one or a few sums.  Central symmetry, AP blocks, and defects
should be explicit outputs.  Feed each output to sheet or near-AP certificates.

**B2. Hypergraph circuit classification.**  Enumerate minimal blocker
complexes with no empty edge after quotienting runner relabeling, pair-sum
multiplicity, and protected endpoints.  Treat these as matroid-like circuits,
but retain the arithmetic realization question as a sidecar.

**B3. Realizability of abstract blocker circuits.**  Many abstract covers may
not arise from residues v_h p modulo v_i+v_j.  Derive congruence compatibility
conditions and use them to eliminate nonrealizable obstruction circuits.

**B4. Tropical normal fan of M.**  The function min_v ||vt|| is piecewise
linear.  Build the normal fan whose vertices are pair-sum events and whose
facets record active runners.  Relate fan circuits to blocker circuits and
safe peels.  Do not discard the rational ruler label.

**B5. Stability around the deep well.**  Prove that M near 14/183 forces a
twelve-speed near-interval core and one resonant outlier.  This would turn the
empirical covering-min rigidity into a true inverse theorem.  Keep the
large-s counterexamples to the unconditional rigidity in the test set.

The separate n=12 non-extremal-core/sporadic tight branch is also under the
S3 uniformity audit above.  THM-763 now makes it uniformly finite, THM-765
removes imprimitive peeled cores, and THM-766 gives a first-window tooth ladder,
while THM-768/769 split it into shallow full-residue and deep multiple-owner
binding packets.  THM-770 proves that the shallow sporadic slice is empty
through `max A<=168`; higher shallow lifts and every deep packet remain open.
Uniform emptiness therefore remains open.  It should remain
logically separate from the q25 covering certificate: one characterizes
equality, while the other was a proposed closure terminal and is now refuted.

The exact max-peel tooth atlas sharpens that residual.  Across all `2,453`
primitive escaping rows in its exhaustive `[1,20]` slice, cyclic nearest-tooth
winding is still one; `1,972` rows also have pure endpoint owners.  Explicit
liars survive both that quotient and a transitive component-phase tournament
with one Hamiltonian path.  The next theorem must therefore use the metric
slack `1/13-||wc||-wh` together with endpoint-owner splice incidence, not
topological winding alone.

**B6. Multi-killer balance as a convex program.**  Replace the one-direction
perturbation bound by a piecewise-linear max-min program over all killer
slopes.  Seek a global lower bound depending on the active core face, not only
the smallest killer.

**B7. Equality-kernel classification.**  Use closed-boundary owners,
pair-sum multiplicity, and exact M to classify M=1/14 and M=14/183 packets
separately.  AP/GW equality at 1/14 must not be conflated with the primitive
covering minimum.

**B8. Scale-aware compactness.**  Seek a projective compactification of speed
sets in which scale rays acquire boundary strata with residue labels.  Prove
that every limit stratum is discharged by a lower-runner theorem, sheet
dodge, or decorrelation certificate.

**B9. Good-period inverse theorem.**  If no a/q with q<=Q is lonely, translate
the residue exclusions into additive structure.  The goal is:
no small good period implies repeated rulers or a near-AP packet.  This is the
natural converse to S312.

**B10. Safe-peel confluence.**  Different peel orders can lead to different
cores.  Determine whether a canonical irreducible kernel exists, or attach a
diamond/branch certificate when it does not.  This is needed for a clean
recursive quotient.

### Priority C — niche historical threads with a precise re-entry point

**C1. Dynamic twist ladder.**  Revisit HYP-2443/HYP-2972 with endogenous
pair-sum rulers only.  Failed twists should add blocker edges; success should
emit a rational certificate.  Fixed external ladders remain forbidden by
THM-566.

**C2. Ramanujan exact-period projectors.**  Use HYP-2978/HYP-2979 to separate
primitive and imprimitive multipliers on each pair-sum ruler.  Test whether
nonunit obligations supply the missing pressure ignored by unit-only scans.

**C3. Fejer/Toeplitz packet certificate.**  HYP-2981's positive kernels may
control signed cancellation without the divergent absolute relation sum.
Attach the kernel to the good-set indicator and demand an explicit error below
the 1/14 margin.

**C4. Mollified discrepancy.**  Revive the HYP-6375 route only with a
Beurling-Selberg or Fejer minorant whose smoothing cost accounts for all
components of G'.  The old total-variation Erdos-Turan bound was about 700
times too weak.

**C5. Signed relation-lattice resummation.**  S312 refuted absolute-value
bounds.  Group relations by exact period, pair-sum ruler, or involution before
taking absolute values.  A valid proof must expose cancellation algebraically.

**C6. Higher-order E3/Gowers sidecar.**  Test whether third-order additive
relations distinguish repeated-ruler near-AP packets from low-multiplicity
spread packets after energy alone ties.  Do not use E3 as a standalone
monotone scalar.

**C7. BSG/Freiman router.**  Use high additive energy only to obtain a bounded
rank progression or symmetric core; then invoke a sharp LRC certificate.
The inverse theorem need not see 1/14, but the terminal route must.

**C8. Schur-triple dead-ruler analysis.**  Identify which Schur triples create
redundant or mutually blocking pair-sum obligations.  Preserve protected
endpoint labels; raw triple counts previously failed.

**C9. Ostrowski/three-gap witness words.**  For normalized exceptional residue
packets, encode sheet-bad sets and good periods by rotation words.  Seek a
bounded partial quotient forcing an uncovered sheet.

**C10. p-adic/Hensel scale descent.**  Factor the core scale prime by prime.
At each p, either lift a witness sheet or prove the exceptional residues
occupy a rigid p-adic pattern.  Reassemble by CRT with explicit coupling.

**C11. Automata for residue packets.**  Build a finite transducer whose state
is the exceptional residue word plus blocker debt under c -> pc.  A finite
state proof is acceptable only if state equivalence preserves the LRC
predicate.

**C12. C27/unital/K33 realization inside blocker complexes.**  Locate the old
state-lift packets as small labelled minors of HYP-6785.  This could turn a
large historical vocabulary into reusable obstruction motifs.

**C13. NORK/pinch boundary templates.**  Attach normal-cone or pinch
certificates to pair-sum events with nearly empty edges.  The target is a
local perturbation that empties the last blocker without losing protected
endpoints.

**C14. Cech barcode with arithmetic labels.**  Track how the danger-arc nerve
changes as delta moves from 1/14 to 14/183 or 1/13.  Require pair-sum and owner
labels on births/deaths; an unlabelled nerve cannot distinguish AP liars.

**C15. Tropical endpoint ownership.**  Treat active runners as a tropical
basis and blocker changes as wall crossings.  Search for a conserved signed
owner current stronger than the raw divergence that vanished on covering
rows.

**C16. Matroid proof obligations.**  Test submodularity, circuit elimination,
and deletion/contraction for protected blocker edges.  If a matroid axiom
fails, record the smallest arithmetic counterexample rather than forcing the
analogy.

**C17. Facility-location/game dual.**  View t as a facility maximizing the
minimum circular distance and runners as adversarial clients.  Seek a
minimax dual supported on pair-sum events.  It must reproduce the exact
THM-668 support theorem.

**C18. Topological parity.**  Borsuk-Ulam and winding arguments are best aimed
at proving an odd number of unblocked boundary events or a sheet parity
statement.  Pure topology only recovers non-covering/equality information.

**C19. Eisenstein/Phi6 arithmetic.**  Explain 183 and the deep-well resonance
inside the pair-sum and scale-sheet complex.  The deliverable is a rigidity
lemma, not numerology around 13*14+1.

**C20. Atkin-Lehner/cusp involutions.**  Re-enter only if an involution pairs
signed Fourier or Ramanujan channels on exact rulers.  Previous cusp analogies
without a preserved LRC observable should remain dormant.

**C21. Height descent.**  Design a Tao-style descent on maximal core scale,
exception count, and ruler complexity.  Each step must strictly decrease a
lexicographic measure and preserve covering or emit a rational witness.

**C22. Box-spline/polyhedral certificates.**  Encode simultaneous danger
constraints as a rational polyhedral complex in one time variable plus scale
parameters.  Export exact chamber certificates rather than a floating mesh.

**C23. Proof-carrying search.**  Every computationally closed family should
emit one of: rational t, empty blocker edge, sheet index, safe peel, exact
interval certificate, or LP dual.  This turns a census into a replayable
finite proof.

**C24. Alternate Tournament Analysis.**  Compare tournaments whose vertices
are runners, pair-sum rulers, sheets, blocker circuits, and proof routes.  For
each, report score histogram, directed cycles, SCCs, edge flips, and
Hamiltonian-path counts, but judge it only by which exact predicate it
preserves.

**C25. Adversarial frontier generator.**  Mutate normalized shapes to maximize
minimum blocker-edge size, smallest good-period denominator, peel
irreducibility, and low M simultaneously.  Include near dilates, multi-gcd
exceptions, and signed-relation-heavy families; uniform random sets are too
easy.

**C26. Nerode minimization of certificate states.**  Starting from the
incoming proof-carrier automaton, test whether two normalized states may be
merged by quantifying over exceptional-runner completions.  When a merge
fails, store the smallest distinguishing completion as the missing sidecar.
This imports the successful colored-gate/private-firewall discipline without
assuming runners are the automaton states.

## Guardrails: results that should not be rediscovered as false shortcuts

1. No fixed external denominator window can prove LRC(14); THM-566 supplies
   families killing every prescribed finite ladder.
2. A per-core finite band is not a uniform raw-height band; v*(cP)=c v*(P).
3. Covering is not dilation-invariant.  State primitivity and the order in
   which gcd reduction is used.
4. Pair-sum enumeration must include nonunit multipliers.
5. Positive open measure is sufficient but not necessary at the closed
   threshold; AP/GW have boundary witnesses.
6. Low-order Bonferroni truncations can be negative on loose band residuals.
7. The absolute relation-lattice sum diverges or is uselessly huge; signed
   cancellation is essential.
8. Raw additive energy, Schur counts, endpoint divergence, conductance, and
   unlabelled tournament class are diagnostics, not known monotone proofs.
9. Random samples underrepresent near-dilate and high-scale gcd families.
10. The proposed uniform q<=25 good period is false (THM-762/MISTAKE-143):
    retain S312 only as a 120-row sample and use the exact blocker deck.
11. THM-724 and THM-726 titles are stronger than the honest scope in their
    bodies; cite the remaining empirical/general lanes.
12. THM-753 proves the peel lemma and dichotomy, not the empirical assertion
    that every irreducible family is already tiled.
13. THM-756 proves the bottom-core battery; its assembly battery is sampled.
14. THM-758 Claim A is proved.  Claim B is a per-core procedure plus partial
    execution, not a completed globally finite census.
15. A runner tournament loses simultaneous hyperedge intersections and exact
    ruler labels.  It cannot replace HYP-6785.
16. Scalar compression should always state what LRC predicate it preserves,
    what data it destroys, and an explicit liar or guardrail.

## Recommended division of labor

The most efficient parallel program is:

~~~text
scale team:
  A1-A6 and A13, with THM-760 as base case

obstruction team:
  A7-A10 and B1-B4, with HYP-6785 as exact object

certificate team:
  A4, A12, C2-C5, C23

rigidity team:
  B5-B9 and C19

adversarial/formal team:
  A6, C11, C22-C26
~~~

The teams should exchange proof-carrying normalized packets, not raw lists of
large speeds.  A useful handoff record contains the scale quotient, exceptional
residue word, exact M or threshold certificate, smallest good period, safe
peel status, pair-sum multiplicity histogram, blocker-edge histogram, and the
information lost by every quotient tried.

## Primary entry points

- 01-canon/theorems/THM-366-lrc-small-denominator-divisibility-sieve.md
- 01-canon/theorems/THM-668-pair-sum-ruler-witness-structure.md
- 01-canon/theorems/THM-738-near-AP-three-slot-closure-all-1001-bodies-in-1-14.md
- 01-canon/theorems/THM-753-safe-peel-reduction-to-irreducible-cores.md
- 01-canon/theorems/THM-755-capped-envelope-H-band.md
- 01-canon/theorems/THM-756-H-band-closure-bottom-cores.md
- 01-canon/theorems/THM-758-far-count-split-covering-reduction.md
- 01-canon/theorems/THM-759-tight-instance-ratio-bound.md
- 01-canon/theorems/THM-760-coprime-exceptional-runner-sheet-dodge.md
- 05-knowledge/hypotheses/HYP-6780-lrc14-capped-band-scale-quotient.md
- 05-knowledge/hypotheses/HYP-6785-lrc14-endogenous-pair-sum-blocker-complex.md
- 00-navigation/LRC14-FRONTIER-2026-07-15.md
- 05-knowledge/hypotheses/HYP-2976-lrc14-holistic-lineage-synthesis.md
- 05-knowledge/hypotheses/HYP-2980-lrc14-holistic-route-atlas.md
- 05-knowledge/hypotheses/HYP-2443-lrc14-marked-ladder-support-gate.md
- 05-knowledge/hypotheses/HYP-2972-lrc14-twist-ladder-dual-certificate.md
- 05-knowledge/hypotheses/HYP-2978-lrc14-ramanujan-divisor-quotient-guardrails.md
- 05-knowledge/hypotheses/HYP-2979-lrc14-ramanujan-exact-period-projector.md
- 05-knowledge/hypotheses/HYP-2981-lrc14-robbins-robin-interval-fejer-packet-certificates.md
- 07-reflections/lrc14-holistic-history-and-current-proof-structure-codex-s160.md
- 04-computation/lrc14_endogenous_blocker_complex_codex_S2.py
- 05-knowledge/results/lrc14_endogenous_blocker_complex_codex_S2.out

## Bottom line

The repository has already learned that the LRC(14) object is not a runner
tournament, one denominator, one density, one Fourier norm, or one finite raw
box.  It is a scale-sensitive family of exact rational witness obligations
with boundary owners and simultaneous blockers.  THM-760 shows that choosing
the correct fiber—witness sheets—can turn the latest obstruction into a proof.
HYP-6785 supplies an exact common language in which the remaining scale,
additive, topological, Fourier, and computational routes can now meet.
