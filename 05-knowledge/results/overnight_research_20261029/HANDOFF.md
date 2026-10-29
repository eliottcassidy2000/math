# Overnight paired arithmetic research handoff

Status: the scheduled overnight sequence is complete on 2026-10-29.
The [morning synthesis](MORNING_REPORT.md) records the completed scoped
results and the remaining open obligations. The eight-run overnight
automation was deleted during the09:03 Denver close-out after the report
and the subsequent explicit user-requested checkpoint were published.
This file remains a continuation map for future requested research;
preserve proof status and do not restart completed work.

Latest explicit user follow-up: commit b39949c673,
[four-vertex gluing and Fermat clocks](../four_vertex_glue_fermat_20261029/README.md),
proves complete Smith forms for every positive power of the fixed T5
three-run response, the exact gluing sequence with its C5 clock, and the
composition-depth/precision threshold. Its normal and optimized checkers
agree on11628 gates, including81/27. Fixed-source matrix powers do not
establish native-complement closure. The next decisive test in that lane
is a mixed response word M_U M_T for an actual legal complement, retaining
the original named ports. Its [Fermat note](../four_vertex_glue_fermat_20261029/FERMAT.md)
supplies exact CRT squaring recurrences and an odd-order witness for F5;
Fermat-prime completeness remains OPEN. These results supplement the
overnight obligations below rather than discharging the unresolved
Selmer realization, ternary index3^10, or exceptional-integer problems.

Current worktree: `codex/session-sha_pairing_f1_20261029`, based on the
clean pushed prior packet `sha_marked_join_f1_20261029`.

## Completed checkpoints: inherit, do not restart

1. [TWO_PORTS](TWO_PORTS.md) proves the complete weighted classification
   H_s+H_(s n), s|gcd(d,e), and the bounded-weight Mobius/convolution
   recurrence. It exhausts all65,536 fixed native cuts and the separate
   simultaneous-C3 subcensus. At weight50 the packed square counts
   104,060,401 cuts. Minimum weights for H25,H27,H57 are4,4,6.
2. [MINIMAL_RATIONAL_PORTS](MINIMAL_RATIONAL_PORTS.md) proves that
   (det M,P^T M^-1 P) is the exact response quotient for retained-span
   determinant queries. It preserves re-marking and global sign reversal;
   integral cokernels and clock-restricted minimality remain separate.
3. [TERNARY_TERMINATION](TERNARY_TERMINATION.md) proves universal
   positive-integer termination with n(8/5)^Q. Incoming work gives stronger
   bounds, five increasing cylinders, periodic exclusions and an exact
   Cantor-cell recurrence. No positive aperiodic exception is excluded.
4. [FAMILY_HEIGHT](FAMILY_HEIGHT.md) proves degree24 for the record4239
   j-map and Theta(X^(1/12)) distinct curves at invariant height X. Its
   named section is exact. Generic rank17 is cited, not proved here.
5. [ARITHMETIC_INTEGRATION](ARITHMETIC_INTEGRATION.md) gives the explicit
   quadratic selector from15 matchings to6 K3,3 matchings, credits the
   incoming Frobenius/Hamming/fruit results and retains the12/13 clocks.
   [ALGORITHMS](ALGORITHMS.md) reads pinned primary techniques; the
   external manuscripts' full proof/build graphs have not been audited.

Six adjacent checkers pass normally and with Python -O, with identical
receipts. See [AUDIT](AUDIT.md) and mode_receipts.json. The earlier M4
Pfaffian sign was repaired in place and recorded in MISTAKES.md. The
August rank30 certificate and weighted81/27/A000568 controls were replayed.
After integrating incoming commit4d9c3ad1ce, the cyclic-kernel, exceptional
cell and all66 high-rank coverage checkers were also replayed successfully.

The first hourly continuation additionally completed
[FINITE_CLOCK_PORTS](FINITE_CLOCK_PORTS.md) and
[EISENSTEIN_CLOCK_JOINS](EISENSTEIN_CLOCK_JOINS.md). The latter includes the
minimal determinant/unmarked-group carriers for matched C3 cross joins,
the all-order norm/Mobius classification, and the full ternary precision
law K(3^k)=7*27^(k-1). It also explains exactly which extra internal probe
reads the forgotten source factors. These are scoped weighted-form results.
The 81/27 and A000568-through20 controls were replayed unchanged.

The second hourly continuation closes two more obligations:
[DYADIC_AND_FIELD_PORTS](DYADIC_AND_FIELD_PORTS.md) proves the complete
all-modulus response count, including K(2^k)=5*2^(3k-6) for k>=3.
It supplies coordinate equivalence, exact refinement branching and the
sharp two-extra-bit reconstruction bound. Do not repeat the now-closed
dyadic census as an open search. Over F64/F256 three field queries recover
all states, unlike the integer rings of the same sizes; the checker
decodes all262,144 and16,777,216 field triples.
[MARKED_CLOCK_PORTS](MARKED_CLOCK_PORTS.md) proves that the fully named
paired-group carrier is(d,e,a,b), including source signs. Eleven branch
observations suffice; there is no single-run survival claim at singular
probes. The induced clock on cokernel classes is T^(-T). These are exact
formal matrix-family results, not an arithmetic Selmer realization.

After integrating incoming4ff64c524c, both `two_run.py` and
`paley_pencil.py` replayed successfully without rewriting their receipts.
The joint-run81/27 separator is now modulo3 (ranks24/22), while the
single-run separator remains modulo9. Preserve both scopes. The all-prime
Paley-pencil Smith formula now includes ramified and singular weights;
the prime-power full Smith law was then open and is now closed by the
incoming PRIME_POWER packet recorded below. Its F27 rank8 hostile still
rejects blindly substituting27 into the prime-only rank formula.
[CYCLIC_COEFFICIENT_GLUE](CYCLIC_COEFFICIENT_GLUE.md) connects the actual
Paley7 coefficient triples(15,12,12),(-11,2,2) to phases(3,0),(-13,0).
Retain h and h-a+b=0 mod3 to invert the map; the clock is A^(-T).
Its image count m^3/gcd(m,3) is NOT the nonlinear determinant count.

The third continuation completes [THREE_RUN_GLUE](THREE_RUN_GLUE.md)
and [CLOCK_MINIMAL_PAIRING](CLOCK_MINIMAL_PAIRING.md). The S3 coefficient
lattice has index54 and an exact four-congruence inverse. T5's full
three-run Smith form is1^52,3^2,6^2,0^4; the unglued split adds four
spurious Z/6 factors. A complete order4 census proves the five-vertex
boundary for this specified split. Paley7's full210-factor Smith form
and six elementary odd-primary translation modules are certified.
At11,79,5861,392809 the module has no reciprocal constituents, so every
invariant bilinear form vanishes and hyperbolic doubling is optimal,
even for the full source automorphism group. This includes the square-
order, abstractly paired11-group: retaining clocks changes the minimum.
The raw Paley7 cokernel has squarefree order part79*392809 and is ruled
out as finite principally polarized Sha by cited Poonen--Stoll Corollary9.
No arithmetic realization or finiteness theorem follows from doubling.

Incoming e9e53cf3ed independently supplies the same S3 lattice and
native2/3-primary layers, plus a stronger T5 complement test: all15
Hamilton-path complements expose all four tetrahedron-boundary kernel
directions in the named weight space. It also closes the prime-power
Paley-pencil law in PRIME_POWER.md using a Jacobi carry recurrence.
Credit those shared/stronger results, and do not restart either census.
After rebase, both incoming `three_run.py` and `prime_power_pencil.py`
replayed their stored receipts unchanged. The latter verifies556 prime
specializations, actual F27/F243/F343 pencils and604 Jacobi sums.
Our documentation check and all100,190 normal/optimized gates pass.

The fourth continuation completes
[PRIMARY_CLOCK_COMPLETION](PRIMARY_CLOCK_COMPLETION.md). All2/3-primary
C7 isotypes now have reversible source charts and O-linear quotient
maps. At2 the partitions are trivial(2,2,1^10), cubic+(4,4,3,1,1),
cubic-(4,1^5). At3 the entire group has trivial C7 action, with partition
(4,2^6,1^2). The all-height minimum index recurrence is P2 of that note.
For Paley7 the C7 dyadic minimum is2^18, versus2^2 without a clock.
The direct-sum C7 minimum is2^30. This distinction is NOT a correction
of the earlier explicitly direct-sum statement: adjoining roots can
be cheaper than adding independent summands.

The next significant obstruction is now proved: an embedding that stays
injective under EVERY ordinary quotient reduction must be a direct
summand. Thus the cheaper finite-level minima cannot silently replace a
compatible quotient tower. Z/2+Z/4 inside(Z/4)^2 is the minimal hostile.
At2 the successive minimum log-index increments are3,6,6,3; at3 they
are1,1,1,1. Their numerical recurrence is not a transition-map theorem.

The same continuation also closes the extra clocks at5 and17. Literal
forms invariant under all126 elements of Aut(Paley7) x S3 are obtained
from one-dimensional alternating-form spaces, nondegenerate away from
zero. All six elementary odd-primary minima therefore retain the full
native vertex/run group: no addition at5/17; doubling at the other four.
Only the higher-primary2/3 extra clock actions remain open. Both new
checkers pass normally and optimized, with identical receipts and all
source maps retained; the81/27 slice controls are unchanged.

Incoming ef06067902 further closes SELECTED native complement recovery
for arbitrary three-port weights at every strong source through7.
COMPLEMENT_KERNEL.md supplies integral left inverses from coprime
minors, full396-class coverage, a B5 bad-path plane that survives vertex
invariance, and minimum unnamed orbit-count recurrences. These are not
the same objects as our paired cokernels. Its integral left inverse
is a useful comparison with P5's quotient-compatible splitting. The
remaining native question starts beyond this proved finite universe
or imposes actual attachment-realizability constraints; do not restart
the T5 or strong-through7 census. The larger8/9 examples are not complete
censuses. The incoming Paley count retains the Frobenius term14m^70.
The full incoming verifier replayed successfully after rebase, including
all396 source classes and its four larger examples; its receipts were
unchanged. The primary completion checkpoint is da03553dd0.

The fifth continuation is [FULL_NATIVE_COMPLETION](FULL_NATIVE_COMPLETION.md).
The2-primary arbitrary-injection minimum now retains ALL126 vertex/run
clocks and is exactly index2^40. The full group has native generators
C=translation*run-cycle and H=vertex-Frobenius*run-swap, with HCH^-1=C^2.
The ten C21 source orbits give partitions(2,2,1^4),(1^3),(3),(4,1),
(4,1),(1,1) in degrees1,2,3,3,6,6. Their full-clock minimum costs are
4,0,18,18 bits across the four reciprocal/self-dual sectors.

The first cost is a C4-squared block with involution[[3,2],[2,1]], whose
determinant-minus-one obstruction rules out every order64 paired ambient.
The cubic pair descends to C8 with involution3 and C2+C16 with involution
[[1,1],[8,5]]. All four index-two extensions fail; an index-four extension
C4+C32 with involution[[3,1],[8,5]] and character6x+y mod8 attains the
minimum. The exact all-size extension recurrence retains2x=v, Sx=x+w,
2w=(S-1)v, (S+1)w=0, with its character lift constraints. It enumerates
marked presentations, not unlabelled group classes. Sextic source bases
descend under actual Frobenius, giving a genuine F64 phase occurrence;
F256 has no primitive21 root and is not the same chart.

At3 the full invariant-form space has common radical(C3)^3+C9+C27,
order3^8. An explicit form attains that radical. The sixth continuation
below excludes completing at that index and raises the lower bound to
3^10; the hyperbolic upper bound is3^18. The actual minimum remains OPEN.
All six extra-C81 characters fail.
The new Heisenberg hostile proves why an optimal restriction radical
does not supply jointly compatible ambient clocks. The2-primary PURE
minimum also remains open: the index2^40 construction adjoins roots.

Named2/3 source projections, coordinate lifts and all four clock matrices are
saved in primary_native_clocks.json. Normal and optimized receipts agree;
an independent kernel audit checks every2x2 matrix over Z/4 and Z/9.
The81/27 singleton slices are replayed. Do not repeat the completed C21
or dyadic minimum calculation; use the strict ternary obstruction below.

The sixth continuation is
[TERNARY_EXTENSION_OBSTRUCTION](TERNARY_EXTENSION_OBSTRUCTION.md).
The index3^8 completion question is CLOSED NEGATIVELY, for ALL invariant
source forms, not only the stored optimal form. The complete first-step
recurrence retains3x=v, Fx=x+w_F, Rx=x+w_R, Sx=epsilon*x, all group
relations, and z(a)=b'(a,x). Sixteen of the18 sign-marked root classes
fail the group equations. With zero root there are81 carry/gauge classes
per sign. Every invariant source form forces z to vanish on the common
radical in all162 cases; an independent row-rank path checks all2916
literal carry pairs. No minimal coisotropic first step is possible.
Square order therefore raises the lower bound to3^10. The five possible
minimum indices are3^10,3^12,3^14,3^16,3^18.

ternary_root_extensions.json retains the literal source hash,16 left
annihilators and162 character-row certificates. Its checker also audits
all4825 affine systems of size2x2 with every RHS over Z/3 and Z/4,
and constructs the nonsplit positive chain C3^2 -> C3+C9 -> C9^2.
The separate ternary_obstruction_verify.py imports no extension solver,
re-enumerates all root classes, checks every literal certificate, rejects
two damaged certificates, and replays81/27. Both normal and optimized
runs agree. This is a same-agent alternative-path audit.

Next: try a native index3^10 completion. Starting from a source form with
radical3^8, it needs one radical-increasing and nine decreasing steps.
The seventh continuation below now excludes the other proposed branch:
there are NO invariant forms with radical3^10. A cheap exploratory pass of13 canonical
first carries with greedy continuations returned to radical3^8 and stopped;
this is NOT a census or an obstruction. No successful3^10 chain is saved.
The full paired minimum now has form2^40*3^a*11^6*79^3*5861^6*392809^3
with a in{10,12,14,16,18}. Keep original port projections and all clocks.

Incoming078be2b84d was integrated after that proof checkpoint. Its
[DEFECT_TRANSFER](../overnight_rank_return_20261029/DEFECT_TRANSFER.md)
gives an all-size joint-port evaluator linear in n for fixed k and
supplied backward-arc cutwidth c. Every first interval complement has
c<=2; component endpoint pairing survives expiry of vertex names. Its
full checker replayed here. The incoming full1320-by1320 three-run
81/27 matrices have ternary ranks430/435 and common mass80614764; this
larger census is INHERITED, not rerun in this continuation. Our own
singleton-slice hostile was rerun by the separate certificate verifier.

The incoming periodic two-run-length complement family has35 boundary
states, minimal rational scalar order18, and minimal ternary order12.
The latter is F3[X]/((X+1)^2*g^2), with g irreducible of degree5 and its
field root of order121. Periodic_defect.py replayed normally and under-O.
This actual F243 chart supplies a precise5/11 connection; the full clock
period is726, and forgetting five nilpotent coordinates preserves that
period while losing observable states. Keep this operator distinct from
our native126-element group action.

The new [PERIODIC_FROBENIUS_LIFT](PERIODIC_FROBENIUS_LIFT.md) constructs
an order-five algebra automorphism X->X^487, where487=3 mod242 and1 mod3.
It conjugates the full recurrence clock to its487th power and retains
all12 coordinates; naive cubing has rank6. The12-by12 Hankel chart H
transports it into ACTUAL future-response windows. Its initial window is
fixed, but its original scalar selector is not: b_2=0 and b_248=2 mod3.
Retain ell'=ell*Phi^-1 to preserve observations under source transport.
All726 orbit states and all144 basis products pass, with identical
normal/optimized receipts. This is NOT a constructed native attachment
or vertex symmetry. The next test is realization in the native operation
language with its selector transported. Do not identify the two ternary
modules, or infer arithmetic Selmer compatibility from either lift.

The seventh continuation closes the source-form census in
[TERNARY_FORM_STRATA](TERNARY_FORM_STRATA.md). The full invariant-form group
is C3^4+C9, exactly729 forms. Radical counts are324 at3^8,324 at3^12,
70 at3^14,10 at3^16, and one at3^18. In particular the radical3^10 branch
does not exist. Independent integer/local Smith calculations check every
form. Explicit clock-intertwining source maps reduce the324 optimal forms
to20 sufficient representatives (four orbits of9 and sixteen of18 under
the certified subgroup). These are not claimed inequivalent under the
entire centralizer. Keep the spanning-tree source maps and port projection.

Any index3^10 completion has a necessary extra symplectic F3 plane
P=V-perp/R. Its run cycle is trivial, its run swap is+I or-I, and its
vertex clock is identity or one of two unipotent symplectic conjugacy
classes. All18 literal actions and six classes are enumerated. The next
test is120 sufficient source-form/plane-action cases, with the extension
carries still unresolved: first construct L=R-perp of index9 over V,
then a coisotropic completion of L of index3^8. No candidate chain is
claimed complete. The overall lower bound3^10 has not increased.

[PERIODIC_PRECISION_OBSTRUCTION](PERIODIC_PRECISION_OBSTRUCTION.md) gives
an all-height, all-step recurrence for the single-selector weighted
observer: log_3|T^r W_h|=12h+5(h-2)_+ +(h-2-2r)_+. A unique contracting
Hensel root has valuation2; its source lattice amplitude has valuation2,
while the actual initial prefix has amplitude valuation4. Thus the
actual scalar sequence modulo3^h has exact preperiod ceil((h-4)_+/2)
and eventual period726*3^(h-1). These are symbolic consequences of the
literal projector and first-return valuations, not fits to finite data.

The full weighted carrier has |T W_5|=3^76 and |T^487 W_5|=3^75, so no
bijection can extend the clock conjugacy at modulus243. On the actual
source orbit alone, y_r->y_(487r) gives compatible bijections through
h=6 and fails from h=7 onward. At modulus2187, b_1=38 differs from
b_529255=1496, but the reindexed response windows coincide. The35 literal
boundary-state responses also cannot realize Frobenius by permutation:
27 move outside the column set. Pure bulk words are excluded by the
centralizer/normalizer distinction. The full operation-language question
is now closed negatively for the incoming bulk/name language below;
further actual complements and general vertex operations remain open.

Incomingbad77c3c2f closes the periodic observer under all six terminal
assignments: rational rank24, ternary rank18, and a sharp integral bank
of four selectors (two over Q or F3). A positive mass-two witness
separates under a terminal swap despite equality of all old observations.
Its full checker replayed here. The earlier single-selector precision
formula is explicitly scoped to the old language; do not silently apply
its12-dimensional formula to this larger arbitrary-weight carrier.

Incoming9c876584a1 gives reflection-fixed path recurrences and a joint
Frobenius chart. The reflection result is INHERITED. The joint chart
checker replayed normally and optimized. Its chosen lift changes names;
[JOINT_FROBENIUS_OBSTRUCTION](JOINT_FROBENIUS_OBSTRUCTION.md) now proves
that EVERY invertible F3-linear lift conjugating T to T^487 must change
the native name three-cycle. The native norm N=I+R+R^2 equals a recorded
degree11 polynomial f(T), while f(T^487)-f(T) has rank5. Any linear map
commuting with R would have to fix N, giving a contradiction. This is
a universal LINEAR obstruction, not a universal nonlinear residue claim.
The incoming positive lift remains valid with transported names/source.
All native name operations commute with bulk, so their entire word
language cannot realize the noncentral Frobenius clock exponent either.

[MORNING_REPORT](MORNING_REPORT.md) consolidates the overnight results,
evidence and remaining obligations. Resume only from these current
boundaries; no full completion at3^10, new A000568 term, or arithmetic
Sha finiteness claim has been made.

Incoming6b83fd8383 strengthens all396 native recovery certificates to
strong targets whenever a second observation is needed; C3 needs zero
moves. Its recovered March interval family has exact ordered-k-cover
recurrence of minimal order3k and a joint-port compiler, with a literal
strong complement that exits the family. The recurrence clock has
Smith(1,1,j^r), a precision-losing endomorphism distinct from our finite
group automorphisms. Both incoming `strong_target.py` and
`run_recurrence.py` replayed successfully after rebase. The core full-clock
completion checkpoint is6c45d7f245.

## Incoming results to inherit with their own scopes

Read the current Session97 README and the current ternary-integer packet,
then their targeted proofs. Recent Session97 additions include:

- controllers/EXCEPTIONAL_CANTOR_RECURRENCE.md: occupied-cell recurrence
  N_L=N_(L-1)+N_(L-6)+3N_(L-11), a dimension about0.441171 and an integer
  sparsity bound. This is not an emptiness or global Collatz theorem.
- supersingular/CYCLIC_KERNEL_RECURRENCE.md: all odd cycle lengths,
  reciprocal factors and odd ord_p(2), including Ogg intersection
  {7,23,31,47,71}. Retain the characteristic separately from cycle length.
- supersingular/ALL_M_AFFINE_PORT_RECONSTRUCTION.md: four-port affine
  reconstruction with its m=3 hostile; this is a different query language.
- ranks/phase2/README.md: all66 listed rank-at-least30 point banks have
  local certificates; two require F3 instead of F2. Three named sections
  and a common A3(2) lattice are proved for the newer families. The full
  classical rank17 certificate belongs to determinant948, not1152.
- halving/phase2/README.md: the two original binary gaps now have exact
  2-primary saturation quotients(C2)^13, index8192. We replayed125,459
  checks from the exact Git blobs at1128ee8635, preserving the original
  source hash and rejecting six hostiles. The raw verifier rejects a
  CRLF checkout; use replay_incoming_saturation.py without weakening it.
- overnight_rank_return_20261029/PRECISION.md has a separate integral
  Smith filtration for the fifteen-Lagrangian incidence map. Do not
  identify its precision counts with our nonlinear determinant quotient.
- overnight_rank_return_20261029/ATOM_PRECISION.md now gives actual Paley
  endpoint cokernels and an explicit Hamming/prime-bank map. Its checker
  was replayed here after incoming d87e5a8292. The81/27 pair first separates
  in its3-power endpoint channel at modulus9: the81 matrix is zero there,
  while the27 matrix has E[9,1]=3171=3 mod9. This is a sharper control to
  retain in future multi-run integer response experiments.
- fractions/FULL_IMAGE.md and DENSITY.md: an arithmetic rational-cover
  guard and one-third lift density, conditional on the explicitly cited
  source-generator theorem. Real positivity alone remains insufficient.

## Next bounded obligations, in priority order

1. The finite-clock unrestricted-integer subproblem is now closed in
   [FINITE_CLOCK_PORTS](FINITE_CLOCK_PORTS.md): orbit sums of pair probes
   are legal and recover every invariant Gram entry. A unipotent infinite
   clock is an all-word hostile. Next restrict the operation set to actual
   cross-block or native cuts; the orbit sums need not belong to it.
   The matched irreducible-C3 cross-block case is also closed for the
   stated determinant and UNMARKED paired-group outputs: retain phase(a,b),
   source product de and gcd(d,e) as specified in EISENSTEIN_CLOCK_JOINS.
   Named cokernel classes, induced clocks and every integer observation
   modulus are NOW CLOSED for the explicit matched C3 family; see the
   two second-continuation notes above. The next test is whether a native
   strong-atom attachment language realizes a separating probe family,
   with original ports, integral gluing and complemented target relation
   retained. Keep branching experiments separate from a supplied single
   run that terminates at singularity. No native minimality follows yet.
   The three-run S3 lattice and small native cases are now CLOSED in
   our continuation and the simultaneous incoming THREE_RUN packet.
   All2/3 clock matrices are NOW CLOSED, and the full dyadic injective
   minimum is index2^40. Ternary index3^8 is now excluded for all source
   forms. Next construct or exclude index3^10, with all joint clocks.
   Use the20 sufficient optimal source forms and six necessary extra-plane
   actions in the seventh continuation; the radical3^10 source branch is
   excluded. The supplied common radical type(1,1,1,2,3),
   nine-coordinate matrices, source maps and all-form first-step
   certificates must be retained. Separately classify
   pure full-clock completions. The root constructions do not retain
   ordinary quotient injectivity. Do not recompute the saved source maps.
   For native attachment minimality, extend the incoming strong-through7
   selected-complement theorem, or impose explicit realization equations
   on its arbitrary weight space. Retain target relation, chosen-path
   frame and paired responses. The B5 bad-path plane survives vertex
   invariance, so a vertex-phase or LCA-only repair cannot determine it.
2. Replace abstract paired-group realization by an actual source-fixed
   arithmetic map. For Sha finiteness one still needs finite prime support
   and bounded primary exponent (or another sufficient argument), with
   all Selmer levels compatible. A square-order count supplies neither.
   Specify the arithmetic transition maps before using a completion
   tower: ordinary quotients, torsion inclusions and actual Selmer maps
   are different. P5 now supplies a direct-sum obstruction for ordinary
   quotient injectivity, even when every finite-level minimum exists.
3. Develop adaptive local-character certificates: select cheap primes
   and observation moduli while keeping root labels, generators, curve
   transport and torsion guards. The two F2 gaps are now explained by
   actual joint halving and complete 2-saturation, in addition to the F3
   certificates. Use those exact source maps as the arithmetic anchor;
   do not infer a Sha class or bound from their subgroup index alone.
4. For the cheap controller, retain the ORIGINAL source, full carry and
   both clocks through resets. Try a decisive arithmetic exclusion of
   aperiodic positive points in the exact failure attractor; if it fails,
   record the missing coordinate. Do not substitute Haar probability.
5. Recover section banks for the actual newer K3 fibrations. Combine
   their exact family-height exponent with proved section/intersection
   data; no unproved rank-tail law may be promoted to an upper bound.

Each continuation should produce a bounded proof, exact experiment or
documented obstruction, then commit and push. Check incoming corrections
before extending these results. Do not repeat a completed census merely
to keep the session busy.

Anchor: marked paired groups and arithmetic realization. Niche: the
positive-integer parser certificate. Wildcard: high-rank families and
the two fifteen-element structures. Sha finiteness, Collatz and LRC(14)
remain open.

Do not submit anything to external leaderboards or message other people.
Commit and push coherent verified checkpoints. A final morning synthesis
must say exactly which obligations remain open.
