# LRC14 after the `q<=25` refutation: the decked threshold object

*codex-2026-07-14-S3.  A synthesis after THM-762--770, HYP-6820, the S2
four-far-cone/colored-threshold-sheaf reflection, and the exact blocker and
affine-suspension audits.*

## 1. The two starting goals did not fail in the same way

The request was to begin by proving two uniform statements.  Exact work
separated them.

The `q<=25` claim is false.  It fails inside the precise residual from which
it was extrapolated, not only on a remote scale ray.  It also fails for a
family with no large gcd packet and with very large lonely margin.  The right
response is not a larger sample or a slightly larger fixed denominator bank;
it is to identify the exact obstruction carried by each modulus.

The tight `n=12` sporadic-branch claim is still open.  Here the audit did not
find a counterexample.  Instead it proved that every deletion core is
primitive and that the entire primitive tight locus has explicit finite
height.  This changes the logical shape of the problem but does not yet give
the structural compression that empties the finite locus.

Concurrently, THM-761 proved a multi-exception sheet theorem: up to six
exceptions of arbitrary gcd type cannot cover all sheets above explicit
scale thresholds (uniformly from scale `43` in the coprime case), with exact
small-scale failure sets.  This matters here because the transparent `q<=25`
counterexample is dispatched immediately by a scale-native sheet theorem.
Fixed rational periods and witness sheets are complementary observers.

The honest outcome is therefore stronger than “no proof found” and weaker
than the requested closure:

```text
q<=25:          refuted, exact replacement object found;
n=12 sporadic:  uniformly reduced, exact remaining theorem isolated.
```

## 2. How the repository's view sharpened recursively

The long history is not a collection of unrelated analogies.  Each useful
view forgot something; the surviving liar named the next coordinate.

| View | What it made exact | What the next obstruction forced us to retain |
|---|---|---|
| torus orbit | simultaneous phases `tS` | arithmetic provenance of contacts |
| danger-arc cover | threshold as circle coverage | open/closed boundary and owners |
| divisor clock | cheap witnesses `1/q`, `q<=14` | endogenous larger denominators |
| Farey/continued-fraction cells | local rational refinement | simultaneous blocker identity |
| endpoint walls | cyclic event order and binders | rational lengths and tie blocks |
| pair-sum rulers | complete exact maximizer list | nonunit multipliers and cross-ruler compatibility |
| Haar/Cech/Baire | open mass versus equality boundary | arithmetic labels on components |
| moments/Bonferroni | finite intersection certificates | higher signed cancellation |
| Fourier/Bernoulli/Fejer | exact discrepancy and scale tails | owner incidence and coherent relations |
| additive energy/Freiman | structured-versus-spread routing | sharp threshold constants and realization |
| safe peel | deletion as an exact reduction | classification of irreducible cores |
| capped envelope | a rigorous tail for each fixed core | a scale-normal rather than raw-height quotient |
| scale sheets | common-factor rays as finite decks; up to six exceptions above explicit thresholds | small-scale and seven-or-more-exception tilings |
| affine suspension | all `cP+R` as slope fibers of one torus field | a descent on slope/offset states |
| blocker complex | exact simultaneous rational obstruction | metric component data and scale action |
| gcd deck | uniform obstruction from translation lifts | primitive-core structure |
| formal proof automaton | composable terminal certificates | a legal state quotient for all completions |

The progression is recursive in a precise sense.  A quotient is tried; two
families in one fiber disagree on the next theorem-facing predicate; the
smallest coordinate separating them becomes a sidecar.  Residues required
owners.  Owners required endpoint phases.  Endpoint order required metric
width.  Fixed-core widths required scale.  Scale required deletion and gcd
decks.  The `q<=25` failure now requires a modulus blocker deck even when all
gcd decks are trivial.

This is why there is no plausible scalar “master invariant.”  The object is
indexed by what we intend to do next.

## 3. The two exact decks discovered in this audit

### 3.1 Modulus deck

For a covering family and `15<=q<=28`, THM-762 replaces a search over all
`a` by

```text
Z_q(S) = q-divisible speeds,
B_q(S) = signed unit pairs represented by speeds modulo q.
```

There is an `a/q` witness exactly when `Z_q` is empty and `B_q` is incomplete.
The entire obstruction at one modulus is therefore a colored cover of
`(Z/qZ)^*/{+1,-1}` plus a zero bit.

The coherent counterexample `26*{1,...,12} union {339}` completes every
target deck by an inverse-pair mechanism.  More revealingly, the incoherent
counterexample

```text
{81,91,131,151,157,196,258,274,313,328,330,339,348}
```

does the same with all deletion gcds equal to one.  Its exact maximum is
`101/470`, so a complete short-period deck says almost nothing about the
continuous lonely margin.

The exact S105 replay reveals how the overclaim arose historically: `8,169`
of `8,260` capped rows do witness by `25`, but `91` do not; the last first
witness occurs at `q=38`.  A high success rate and a small observed tail were
compressed into a false universal sentence.  Keeping the certificate deck
per row prevents that quantifier loss.

### 3.2 Deletion/gcd deck

For a core `P`, a maximizer `t_0`, and `d=gcd(P)`, the lifts

```text
t_0, t_0+1/d, ..., t_0+(d-1)/d
```

preserve every core phase.  A completion `w` sees, after repetitions are
removed, a complete deck of order

```text
D = d/gcd(d,w).
```

The deck has a phase at distance at least `(D-1)/(2D)` from zero.  At the
`n=12` threshold this forces `D=1`; primitivity of the full set then forces
every deletion core primitive.

The proof is scale-uniform because increasing the common scale creates more
witness sheets rather than a larger unstructured integer.  This is exactly
the kind of conversion that HYP-6780 said was missing from raw cutoff
arguments.

THM-761 makes the same conversion for several completion speeds.  Its exact
bad-sheet sets form a discrete covering instance on `Z/cZ`; the tight sheet
case is a cyclic tiling problem.  The remaining scale frontier is therefore
not “large `c`” in general, but small-scale criterion failures, at least seven
coupled exceptions, gcd descent, or families with no coherent scale core.

The two decks are complementary.  The gcd deck controls lifts of one core
maximizer and becomes silent when the core is primitive.  The modulus deck
can remain complete when every gcd deck is trivial, but it forgets safe
component geometry and exact `M`.  Neither subsumes the other.

## 4. The exact residual in the tight twelve-speed problem

Let `A=P union {w}` be primitive, `|P|=11`, and suppose the sporadic
conditions

```text
M(P)>1/12,       M(A)=1/13.
```

For each component `J=(c-h,c+h)` of

```text
E(P)={t:min_{p in P}||pt||>1/13},
```

THM-765 says that the completion must satisfy

```text
||wc|| + wh <= 1/13.                                  (*)
```

This is not an estimate.  It is exactly equivalent to placing the closure of
`J` inside one `w`-tooth.  Thus full sporadic emptiness has become an
arithmetic interval anti-cover theorem:

> A primitive super-lonely eleven-core cannot have all of its strict safe
> components swallowed, one by one, by the teeth of a primitive completion
> satisfying the ratio and finite-height constraints.

The remaining data are all discrete or rational:

```text
hereditary gcd deck:       gcd(A\{a_i})=1 for every i,
height:                    sum a_i<=78^11,
max-peel ratio:            w<12*max(P),
active ruler:              13 divides a_i+a_j,
component stalk:           (c_J,h_J,endpoint owners),
completion obstruction:    sigma_J(w)=1/13-||wc_J||-w h_J.
```

A sporadic row is exactly one with `sigma_J(w)>=0` for every component.  The
unique tight row through maximum speed `30` is `{1,...,12}`, but a bounded
census does not control the enormous primitive domain supplied by the height
theorem.

## 5. The tooth-label winding word: exact checksum, not separator

Condition (*) is componentwise.  It has not yet used the global order of the
components around the circle.

If a completion covers every component, each component has a unique tooth
label in `Z/wZ`.  More generally, away from midpoint ties, assign the nearest
tooth to each cyclically ordered component `J_1,...,J_r` and record the cyclic
increments of the labels.  In a covered row the labels can be lifted in
cyclic order to a nondecreasing integer word

```text
k_1 <= k_2 <= ... <= k_r <= k_1+w.
```

After one circuit the total cyclic increment is an integer multiple of `w`;
that integer is the winding checksum.  The full metric word attaches to each
letter

```text
(k_j, c_j, h_j, left owner, right owner, active ruler, sigma_j).
```

Thus it retains

- cyclic nearest-tooth labels and increments;
- endpoint-owner transitions imposed by `P`;
- midpoint phase and half-width;
- the exact slack `sigma_J(w)`;
- deletion and scale actions on the same record.

The exact atlas carried out that falsification test.  It exhausts every
primitive twelve-subset `A` of `{1,...,20}` with `w=max(A)`,
`M(A\{w})>1/12`, and `M(A)<=1/10`, together with the eleven eligible tight
AP deletion controls.  The component criterion agrees with an independent
exact-maximum computation on all `2,464` rows.

All `2,453` escaping rows nevertheless have winding checksum one, and
`1,972` have pure endpoint ownership.  The first liar to “pure endpoints plus
winding one implies cover” is

```text
A={1,...,10,12,13},  w=13,  M(A)=1/11,
min_J sigma_J(w)=-11/26.
```

Requiring a transitive component-phase tournament with one Hamiltonian path
still fails for

```text
A={1,...,9,11,12,15},  w=15,  M(A)=1/10,
min_J sigma_J(w)=-43/91.
```

So winding is a consistency checksum, not the missing coordinate.  The old
one-component ratio proof sees one letter; safe-measure and component-count
bounds see aggregate width and word length; winding sees only topological
degree.  The theorem-bearing carrier must be the **metric labelled word**:
interval, phase, width, owners, and exact slack.  A viable global lemma must
transport or conserve these slacks across owner transitions, rather than
infer coverage from label topology.

## 6. Tournament threads: when the orientation carries mathematics

The repository has repeatedly used tournament fingerprints—score
histograms, directed triangles, SCCs, edge flips, and Hamiltonian-path
counts.  This audit gives one negative and one positive criterion for their
use.

### A negative example: modulus tournaments

Using `q=15,...,25` as vertices, two natural gauges compare pair retention
and shortest blocker cost.  For both counterexamples the switch flips
`47/55` edges, yet both tournaments remain transitive with one Hamiltonian
path and the exact witness verdict is unchanged.  The orientation is
telemetry.  The zero/sign-pair sidecar carries the proof.

### A positive example: deck-lift tournaments

Using gcd-deck lifts as vertices, the orientation comes from signed circular
phase displacement.  Odd decks give a cyclic tournament with one SCC and
directed cycles.  That wrap is the combinatorial shadow of the metric fact
that the deck cannot fit in a short arc.  The tournament still needs the
radius `(D-1)/(2D)`, but its cycles now arise from the exact group action
rather than an arbitrary ranking.

This suggests a legality test:

> A tournament quotient is theorem-bearing only when its orientation is
> functorially induced by the action whose non-linearity obstructs the target
> cover.  Otherwise it is a diagnostic ordering and must travel with the
> exact blocker/incidence sidecar.

For the tooth-label word, the exact atlas used safe components as vertices,
signed displacement of within-tooth midpoint phases as the pairwise
observable, clockwise half-circle orientation as the gauge, and chronological
component order as the tie Hamiltonian path.  The resulting `354` score/SCC
profiles are rich telemetry, but the explicit transitive liar above proves
that even `C3=0`, singleton SCCs, and `H(T)=1` do not preserve component
containment.  The orientation sees phase order and forgets the `w h_J` width
term.  Thus “all components are covered” remains a metric bipartite incidence
predicate; the full component intervals and slacks cannot be replaced by
their tournament.

Alternate vertex sets explicitly considered in this audit were runners,
gaps, fixed sections, section boundaries, wall crossings, residues, cover
arcs, Fourier modes, matroid circuits, rational obligations, moduli, deck
lifts, safe components, killer teeth, and proof states.  The correct choice
depends on which action and which predicate are being preserved.

## 7. The emerging common object

The S2 reflection described the LRC14 frontier as a four-far arithmetic cone
emitting a colored threshold sheaf.  The present audit adds the missing deck
structure.  A useful local record is

```text
base address:
  normalized core shape + scale/offset word + divisor mask

threshold stalks:
  safe components (midpoint, half-width, endpoint owners, boundary status)

arithmetic decks:
  deletion gcd lifts + modulus zero/sign-pair obligations

endogenous rulers:
  pair-sum event, multiplier, protected owners, blocker hyperedge

transitions:
  delete + complete + scale/shear + refine modulus

proof label:
  q-witness / safe peel / sheet dodge / capped envelope /
  component escape / exact finite certificate / named tight atom.
```

Call this an **observer-indexed decked threshold sheaf**.  “Sheaf” records
that local certificates must agree under deletion, scale, and slope changes;
“decked” records finite arithmetic lifts and modulus obligations; “observer-
indexed” records that the legal quotient depends on the next operation.

Its two simplest liar tests are now canonical:

- `S*` has complete short modulus decks but is continuously very loose;
- an imprimitive super-lonely core has wide safe stalks but its gcd deck
  forces a completion escape;
- `{1,...,9,11,12,15}` has winding one and a transitive component-phase
  tournament, but its negative metric slacks force escape.

Any proposed invariant that identifies either pair without an explicit
reconstruction map is too coarse for the corresponding theorem.

## 8. Reframing the two hard problems

The covering `f>=4` branch should no longer be phrased as “find a universal
small denominator” or “enumerate to a universal height.”  A sharper form is:

> Every primitive normalized affine-slope state either has an incomplete
> endogenous modulus/blocker deck, a safe sheet/component escape, a capped-
> envelope terminal, or a strictly smaller descendant in a well-founded
> decked-threshold order.

THM-761 supplies the safe-sheet terminal for coherent cores with at most six
exceptions above its exact thresholds.  Any proposed global descent should
use those thresholds as terminal faces and focus its new mathematics on the
discrete sheet-cover failures and the incoherent no-core stratum.

The tight `n=12` branch should no longer be phrased merely as “show the AP is
unique.”  Its exact hard kernel is:

> No hereditarily primitive super-lonely eleven-core under the MSS height
> bound admits a completion whose metric labelled component word has
> nonnegative slack at every letter.

THM-768--770 refine this into two recursive charts: shallow full-residue
packets above lift height twelve, and deep binding-scale packets whose
off-sheet runners must persistently own every lift of the quotient core's
loose set.  The component-word statement is the common geometric form; the
residue/sheet split supplies its arithmetic addresses.

The second wording exposes where the known theorems enter and what they do
not see.  It also gives exact counterexample output if the proposed local-to-
global principle is false.

## 9. Ranked next contributions

1. Prove the unbounded primitive descent trigger suggested by THM-770: a
   zero-defect full-residue packet is `{1,...,12}` or has nontrivial gcd.
2. Attack THM-769's exact `s=2` and `s=3` folded residuals as persistent
   colour covers of the quotient core's entire loose set, retaining endpoints.
3. Derive a global balance identity for `sum_J sigma_J(w)` refined by endpoint
   owner transitions; test whether it is monotone under deletion or scale
   descent.
4. Extend the exact atlas by structural strata rather than raw height:
   one endpoint-owner orbit, one pair-sum ruler orbit, fixed component count,
   and bounded relation rank.  Store the first liar to every quotient.
5. Prove a primitive metric component anti-cover lemma first for cores whose safe
   components form one orbit under a pair-sum or endpoint translation.
6. Couple THM-762 decks across endogenous pair-sum moduli and ask whether
   simultaneous completeness forces a bounded marked relation or repeated
   ruler.
7. Put the gcd deck and modulus deck on the affine suspension of HYP-6815;
   use a descent measure that counts unresolved deck orbits, not raw height.
8. Formalize the elementary deck and component equivalences before attempting
   the global inverse theorem; they are small, stable interfaces for the
   existing proof automaton.

## 10. What has and has not moved

The standard fourteen-runner conjecture remains open in the `f>=4` primitive
covering branch.  The global covering minimum remains conjectural outside its
proved and finite lanes.  The primitive tight twelve-speed sporadic branch
remains open despite the new finite height.

What did move is the structural frontier:

- a tempting finite-period closure was decisively removed;
- its exact residue obstruction was characterized;
- every imprimitive twelve-speed deletion core was eliminated uniformly;
- an infinite high-winding complete-residue family was eliminated uniformly;
- primitive tight twelve-speed instances received an explicit global height;
- the remaining equality theorem became an exact component-tooth anti-cover
  problem;
- the maximum cannot be the unique `13`-multiple, and every deep binding
  packet needs at least two off-sheet tighteners;
- the primitive shallow/full-residue sporadic branch is empty through
  `max A<=168`, by an exact `13^12` endpoint-owner classification;
- exact bounded liars eliminated winding, pure ownership, and transitive
  component-phase tournaments as sufficient quotients, isolating metric
  slack transport as the sharper target.

The durable lesson is not simply that more information is needed.  It is
which information, under which observer, and why.  The LRC object visible at
the current frontier is a cyclic metric incidence structure with arithmetic
decks and recursive transition maps—not a naked set of runners and not a
tournament without sidecars.

## 11. Live-mainline pull: endpoint current and normalized boundary intensity

Two concurrent advances sharpen this picture further.

First, THM-766 supplies a guaranteed first safe window and eleven projective
tooth cones at `n=12`.  The endpoint-splice audit then identifies the exact
local-to-global bridge.  For a core-safe component `J`, put

```text
e_w(J)=#pi_0(J intersect {t:||wt||>1/13}).
```

Then

```text
chi_13(P union {w}) = sum_{J in pi_0(G_P)} e_w(J),
sigma_J(w)>=0  iff  e_w(J)=0.
```

Thus negative slack is not merely a failed inequality: it emits endpoint
current, possibly with multiplicity.  The correct recursive carrier is an
endpoint-owned metric word whose letters carry `e_w(J)`, not its winding
checksum.  In the full-nonzero-residue height-one cube, the endpoint audit
finds one zero-defect row—the nonprimitive doubled AP—and defect at least two
on all `4,094` primitive rows.  Extending that splice coherence beyond the
height-one packet was open at that point; THM-770 now extends it exactly
through lift height twelve, while the unbounded trigger remains open.

Second, the newly proposed scale complementarity was falsifiable and false.
For

```text
P_N={1,...,11,N},
```

no divisor at least two divides seven members, so `c*(P_N)=1`, yet the fixed
safe interval `(1/14,13/154)` is cut by the `N`-comb into
`N/77-O(1)` distinct components.  Raw fragmentation is therefore not a
divisibility detector.  What survives is the endpoint bound

```text
r_P<=sum P<=|P|max(P)
```

and the scale-normal intensity

```text
kappa(P)=r_P/(max(P)|G'_P|),      v*(P)/max(P)=kappa(P)/pi.
```

Only `kappa>pi` can keep a genuinely larger peel below the capped-envelope
edge.  For `P_N`, `kappa->7/6`, so its exploding component count is benign.
This reframes the covering residual from “bound fragmentation by scale” to
“prevent endpoint-owned safe measure from collapsing relative to boundary
intensity.”

The concurrent sheet-Kirchhoff proposal met the same guardrail: seven
exception decks can maintain an overlapping cover of fourteen sheets with
total multiplicity `24`, so maintained coverage need not be an exact tiling.
THM-771 now makes the correction exact: the free-sheet defect is
`F=Q+Omega-sigma`, where capacity slack and overlap debt compete with
ramification surplus.  The original KCL necessity is withdrawn.  Sheets,
endpoint events, and safe components are related stalks, but their tournaments
and tiling graphs are not interchangeable without an explicit
predicate-preserving transition.

## 12. Live-mainline pull: binding scales and the bounded owner-CSP

THM-768 first separates the zero-owner cases.  If a prime `p` divides only
the maximum of a set of at most `p-1` speeds, a missing nonzero residue can be
gauged to `-1` and a height-weighted perturbation clears every runner strictly
above `1/p`.  At `p=13`, a tight twelve-set therefore cannot have its maximum
as the unique `13`-multiple.  Residues alone do not prove this; the decisive
sidecar is the owner's relative height.

THM-769 then exposes an endogenous scale at every binding maximum
`p/(13s)`.  The shallow chart `s=1` is exactly the full nonzero residue
transversal.  In the deep chart `s>=2`, the on-sheet quotient core is loose,
so off-sheet runners must cover all `s` lifts of one of its maximizers.  The
capacity observable is

```text
sum_(w in F) (floor(2D_w/13)+1)/D_w,   D_w=s/gcd(s,w),
```

and it must be at least one.  This kills one-exception packets, forces the
two-exception branch to `s=2`, and isolates an `s=3` equality edge for three
exceptions.  More importantly, it identifies the recursive object: a
persistent sheet-colour ownership over the **whole loose set** of a quotient
core.  A binding residue packet at one time is too weak.

THM-770 closes the shallow chart through lift height twelve.  It replaces
runner vertices by `24,008` atomic endpoint cells and 156 residue-labelled
speed options.  Unique-owner propagation enumerates every grouped cover among
`13^12` conceptual packets and finds exactly thirteen zero-defect leaves:

```text
c*{1,...,12},   c in {1,...,12,14}.
```

Only `c=1` is primitive.  Hence the bounded shallow sporadic branch is empty
for `max A<=168`.  The diagnostic tournament on the twelve residue-choice
obligations has all 66 burdens tied, a transitive tie path, no triangles, and
one Hamiltonian path: it is maximally uninformative while the incidence
hypergraph separates thirteen covers from more than `2.3*10^13` candidates.
This is the cleanest tournament guardrail in the audit.  Pairwise orientation
can summarize pressure; only the owner hyperedges preserve simultaneous
compatibility.

The bounded theorem does not combine with THM-763 into a practical closure:
the latter's height is vastly larger, and the deep multiple-owner chart is
outside the full-residue CSP.  The exact frontier is therefore no longer
“prove residue rigidity.”  It is to prove either an unbounded gcd descent in
the shallow chart or a failure of persistent sheet-colour ownership in every
deep chart, with the metric endpoint current retained in both.

## 13. Live-mainline pull: the seven-owner defect is an Euler current

THM-771 supplies the exact invariant missing from the seven-exception sheet
discussion.  For owner `a`, write

```text
g_a=gcd(w_a,c),  C_a=c/g_a,  u_a=w_a/g_a,
A_a=g_a ceil(C_a/7).
```

If `F` is the number of free sheets, `Q=sum(A_a-|B_a|)` is unused capacity,
`Omega=sum_k max(d_k-1,0)` is overlap multiplicity, and
`sigma=sum A_a-c` is ramification surplus, then

```text
F=Q+Omega-sigma.
```

This is an owner-labelled Euler defect, not a moment estimate.  In the
unramified layer `7|C_a` for every owner, `sigma=0`; every strict endpoint event
drops an owner's load by `g_a` and necessarily opens a free sheet.  The event
mesh is `1/u_a=g_a/w_a`, not `1/w_a`, which corrects the raw-winding scale
mistake.  A core-safe interval therefore closes the stratum once some reduced
winding satisfies `u_a>=7 max(P)`.

In ramified strata, surplus can be paid entirely by overlap.  The primitive
`c=21` exact row in THM-771 has `(F,Q,Omega,sigma)=(0,0,12,12)` on an open
chamber.  This explains both why a naked union bound stalls and why a static
exact-tiling/Kirchhoff analogy overreaches: equality is safe at a switch, while
ramification changes capacity away from the switch.

Tournament Analysis becomes precise here.  Owners may be oriented by private
sheet adjacency, with owner order as the tie path, but that quotient discards
all three scalars `Q,Omega,sigma`.  Its cycles and SCCs describe the private
word, not the free-sheet predicate.  The faithful vertex set is bipartite:
owner rows and sheet obligations, with multiplicity and endpoint state.  The
new identity is exactly the reconstruction sidecar a tournament would need.
