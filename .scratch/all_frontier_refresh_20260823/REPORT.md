# Bounded all-frontier refresh — 2026-08-23

**Status:** SCRATCH RESEARCH REPORT + FINITE-EXACT PROBES.  No canon file,
theorem namespace, hypothesis, navigation page, or message was edited.  The
four new probes below are diagnostics and bounded evidence, not promoted
results.

## Freshness and truth boundary

The assigned worktree was already dirty with 27 unrelated untracked paths;
all were preserved.  At startup its branch
`codex/session-all-frontiers-20260822` was 14 commits ahead of and 3 behind
`origin/main`, so no pull or rebase was attempted.  Repeated `git fetch
origin --prune` calls were read-only.  During the concurrent session the local
branch temporarily acquired the incoming ancestry and additional local
research commits, then the remote advanced again.  The last origin object
integrated into this report is `8a9ee9aff`; concurrent local promotions kept
the aggregate branch divergent from `origin/main`.  New remote commits were
inspected with
`git show origin/main:path`; this session itself did not rebase the dirty
shared worktree.

The important point is that the startup snapshot did not stay current.  The
complete incoming chain actually encountered was:

1. THM-3791's intrinsic universal Hensel root was promoted through a
   finite-etale independent audit (`f7093028b`);
2. THM-3794 and THM-3795 were promoted and independently audited
   (`4bbc62f67`, `a6c9e251c`, `1428ba19b`);
3. the THM-3797 confluent namespace moved from an honest reservation
   (`06734d90d`) to proved (`a58b7ea2f`) and then **PROVED + VERIFIED-EXACT +
   INDEPENDENTLY HOSTILE-AUDITED** (`dbf3385b3`);
4. THM-3798 was reserved as an empty common-AP step-three stub
   (`2b4c50453`); a separate scratch proof candidate exists, but is not a
   dependency; and
5. THM-3799 was promoted at `88f228619` and independently audited at
   `e225f67ca`, closing monomial `r e^m` repairs but leaving nonmonomial
   `r g(e)` and mixed corrections open; and
6. origin reserved five further lanes: THM-3800 sharp
   `{-6,1}` carrier/mate support, THM-3801 cubic normalization and companion
   sheets, THM-3802 contact-tree de Rham, THM-3803 affine-linear r-repair, and
   THM-3804 Rule 30 all-period Smith law.  A transient reservation race had
   given the last two files the same THM-3803 ID at `ac2a64986`; commit
   `791658bcc` repaired the Rule 30 file to THM-3804; and
7. THM-3802 was then promoted at `2f95a00fd` to **PROVED + VERIFIED-EXACT +
   PENDING INDEPENDENT HOSTILE AUDIT**, abstracting the actual plane-chart
   contact-tree residue law, and THM-3803 was promoted at `fb00b83c4` with the
   same status, closing every affine-linear `r(b_0+b_1e)` repair.  The live
   IDs are unique; and
8. on the current session branch, `f89c0f53a` promoted THM-3804 to **PROVED +
   VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**, giving the all-period
   amplitude image/kernel/cokernel/Smith law and ruling out physical finite
   scale-cycles at every declared period.  Every Rule 30 center-column prize
   remains open; and
9. the final origin wave corrected the reserved THM-3800 direction into the
   audited theorem that the sharp torus-escaping carrier has exactly fourteen
   reduced critical points (`f0140cd1e`, `b2a89a8a5`), promoted THM-3801's
   cubic-normalization/companion-sheet gate pending independent audit
   (`93eb1ac33`, sharpened by `e37585090`), independently audited THM-3802 and
   THM-3803 (`6348c09e1`, `7d0fe0049`), and reserved THM-3805 for the general
   quadratic r-repair obstruction (`8a9ee9aff`); and
10. the current session branch promoted THM-3798 at `ed11afa2e` to **PROVED +
    VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**, closing only exact pure
    common-AP 2x4/4x2 step-three support cells.

The live theorem surface used here is therefore:

- [THM-3790, cubic pseudoplane arm nodal-immersion gate](../../01-canon/theorems/THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate.md),
  [THM-3791, finite-etale intrinsic resonant-jet law](../../01-canon/theorems/THM-3791-moving-root-danielewski-resonant-jet-de-rham-law.md),
  [THM-3792, pure first-normal carrier obstruction](../../01-canon/theorems/THM-3792-pure-first-normal-nodal-carriers-have-critical-points.md),
  [THM-3793, inert primitive-sum singleton](../../01-canon/theorems/THM-3793-inert-prime-sum-all-scale-two-cube-singleton.md),
  [THM-3794, constant-unit quadratic-etale exclusion](../../01-canon/theorems/THM-3794-constant-unit-surfaces-have-no-quadratic-etale-plane-map.md),
  [THM-3795, r-independent quadratic-normal carrier obstruction](../../01-canon/theorems/THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points.md),
  and [THM-3796, first 2x4 support peel](../../01-canon/theorems/THM-3796-cubic-pseudoplane-first-two-by-four-support-peel.md)
  are **PROVED** with their displayed exact/independent-audit qualifiers;
- [THM-3797, confluent quadratic Hermite completion](../../01-canon/theorems/THM-3797-confluent-quadratic-hermite-jet-completion-no-go.md)
  and [THM-3799, monomial r-repair obstruction](../../01-canon/theorems/THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points.md)
  are **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** on
  `origin/main`; and
- [THM-3798, common-AP step-three support peel](../../01-canon/theorems/THM-3798-cubic-pseudoplane-common-ap-step-three-support-peel.md)
  is **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** on the
  current session branch at `ed11afa2e`; and
- `THM-3802-plane-chart-contact-tree-resonant-de-rham-law.md` and
  `THM-3803-affine-linear-r-repairs-of-nodal-carriers-have-critical-points.md`
  are **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** on
  `origin/main`;
- `THM-3800-sharp-torus-escaping-nodal-carrier-has-fourteen-critical-points.md`
  is **PROVED + VERIFIED-EXACT + CORRECTION OF RESERVED DIRECTION +
  INDEPENDENTLY HOSTILE-AUDITED**, and
  `THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate.md`
  is **PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT** on
  `origin/main`;
- `THM-3804-rule30-all-period-amplitude-lattice-smith-law.md` is **PROVED +
  VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** on this session branch at
  `f89c0f53a` (the `origin/main` copy is still the earlier empty reservation);
  and
- `THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points.md`
  is **RESERVED / UNPROVED EMPTY STUB** on `origin/main`.

The former `THM-3800-sharp-nodal-carrier-gap-seven-four-weight-mate-gate.md`
reservation is **SUPERSEDED** by the corrected THM-3800 slug above before its
support enumeration begins.  The only new empty namespace in the final wave
is THM-3805; it has no mathematical force.

The bounded coordination peek exposed broadcasts MSG-3061--3068.  The latest
usable mathematical signals were the THM-3743 flatness audit and its promoted
canon, the Rule 30 adaptive-observer report, and the THM-3756 ordinal audit.
The Rule 30 broadcast's bounded-gap/overflow result has no routed canon theorem
in this worktree, so it is treated below as a message-level lead, not a proved
dependency.

## Active concept board

| Object / representation | Predicate | Current operation | Lost coordinate | Cheapest decisive test |
|---|---|---|---|---|
| LRC rank-eleven star plus one `l1<=356` Graver row | exclude one of 165 rows | intersect bounded relation atlas with each star space | owner, sign partition, phase, arrival | exact per-star short-vector intersection, AP positive control |
| Labelled state with a root/owner | fixed-root response survives quotient | switching, blocker relabeling, reset descent | the covariant root or stabilizer-fixed choice | enumerate smallest same-shadow/different-root pair |
| Finite observer under its next native operation | observer equality is a congruence | Rule 30 section, FC reset, AMM continuation | carry, overflow, dynamic root, phase invoice | one more operation with all sidecars updated lawfully |
| Smooth surface plus normal/branch data | Darboux or etale plane entry | Hensel lift, affine modification, canonical normal correction, finite normalization | resonant coefficient, quadratic nonmonomial `r g(e)`, mixed support anchors, unramified branch sheet | adjacent-quadratic anchored support gate and first quadratic `g` resultant |
| Prime-local address | exact polygon or singleton fibre | residue lift / change of representation | exponent, gcd, root residue | structured prime-power hostile before extrapolation |
| Paid labelled forcing matrix `[B|D_rho B]` | legal Arithmetic-Kakeya forcing | distinguished-coloop peeling | grammar ties, slope labels, row provenance | exact same-`H` compiler or refuter at the `5/3` rung |

The operative method cards were “controlled forgetting requires a sidecar,”
“audit sections under their next native operation,” “fill operation columns,”
and “test filtration--observer commutation before scalarizing.”

## Topic-by-topic inheritance refresh

### 1. LRC(14) — **OPEN**

- **Closest proved mechanism.** [THM-3743, lonely-runner polyhedron
  Khinchin-flatness relation reduction](../../01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md)
  turns a hypothetical counterexample into one Graver speed relation of
  `l1<=356`, then joins it to the rank-eleven/rank-twelve split.  [THM-3731,
  valuation-owner equivariance](../../01-canon/theorems/THM-3731-lrc14-valuation-owner-equivariance-and-semantic-packet-boundary.md)
  repairs only the static tied-owner selector by retaining full positive
  coefficients.  No one of the 165 rows is excluded.
- **Canonical hostile.** The AP speeds `(1,...,13)` have the norm-four row
  `(1,-2,1)` and are safe at the exact boundary time `1/14`: short relation is
  necessary, not sufficient.  `26*{1,...,12} union {339}` separately refutes
  a uniform good-period bound `q<=25`.
- **Least-used sidecar.** THM-3731's joint table retaining canonical owner,
  predecessor owner root, translated deep root, endpoint mask, terminal word,
  and actual arrival.  A coefficient tie-break alone is not semantic drift.
- **Cheapest new exact probe.** Intersect the full `l1<=356` Graver atlas with
  each rank-eleven star space, stratify support and sign, and attach the actual
  coefficient-resolved owner before integrating.  The 5,855 two-cube addresses
  found below are a lossless index on one support-two subatlas, not an LRC
  certificate.

### 2. Jacobian / Dixmier — `JC(2), DC(2)` **OPEN**

- **Closest proved mechanism.** [THM-1300, explicit Jacobian counterexample in
  dimension three](../../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)
  settles the higher-dimensional boundary, not the plane.  [THM-2049, DC(2)
  Ore correction complex](../../01-canon/theorems/THM-2049-the-DC2-Ore-boundary-correction-complex-is-acyclic.md)
  gives an exact lift in the beta-completion, while polynomial termination is
  open.  THM-3791 now gives the intrinsic finite-etale resonant coefficient;
  THM-3794/3795 close quadratic-etale and every r-independent canonical
  nodal-carrier entrance.  Audited THM-3797 completes the
  first confluent quadratic Hermite surface by two affine-modification jets
  and still finds residue `(1,0,0)` modulo scalars, while THM-3799 closes every
  monomial `r e^m` repair.  Audited THM-3802 abstracts the common residue law
  to an actual finite plane-chart cover and audited THM-3803 closes every
  affine-linear `r(b_0+b_1e)` repair.  Pending independent audit, THM-3801
  forces any surviving constant-unit degree-three normalization to be
  nonmonogenic with transposition branch type and one visible companion
  sheet.  None classifies arbitrary planar Keller maps.
- **Canonical hostile.** [THM-2044, rank-two Poisson suspension](../../01-canon/theorems/THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension.md)
  does not quantize to `DC(2)`.  Completion solves every next graded equation
  without forcing a polynomial endpoint.  THM-3794's sharp hostile is
  `G_m x A1 -> A2`, `(u,v)->(u^2,v)`, when nonconstant units are allowed.
  THM-3797 also shows why stopping at the first singular hypersurface or its
  normalization is false: two further jets are needed before the smooth exact
  image is visible.  THM-3802 assumes an actual smooth affine plane-chart
  cover; a decorative contact tree alone constructs no surface.
- **Least-used sidecar.** Polynomial-support termination in the Ore lift and,
  on a surface factor, the actual affine chart cover and inverses,
  finite-etale root algebra/resonant coefficient, complete modification tree,
  unit group, and branch-sheet decomposition.  Stable equivalence does not
  preserve the displayed rank or formula.
- **Cheapest new exact probe.** Track support/degree growth for the next Ore
  correction.  In the surface lane, use THM-3803's affine-linear cell as the
  positive control, then test the first adjacent quadratic correction
  `r(g_1e+g_2e^2)` with the mandatory `e^2+f_0z` anchors retained.  Its support
  gate is executed below; the next positive test is the exact `(e,u=re)`
  critical resultant for a general quadratic `g`.  The cubic
  unramified-sheet hostile is a separate boundary control.

### 3. HFC(3) / FC(3) — **OPEN**

- **Closest proved mechanism.** [THM-3732, factorial cross-face reset relation
  census](../../01-canon/theorems/THM-3732-factorial-cross-face-shared-reset-edge-sign-separation.md)
  completely audits the named `F12/F13` 239-state, 568-edge graph.  [THM-3466,
  factorial face Stokes/Keller boundary current](../../01-canon/theorems/THM-3466-factorial-face-stokes-and-keller-boundary-current.md)
  supplies the separate HFC-null boundary primitive.
- **Canonical hostile.** Raw equal positive magnitude has only one non-full
  edge and never persists for two edges.  Identity rows instead retain a
  three-edge 12-row corridor, refuting a universal no-subatlas claim; the
  rooted-core survivor is only an `S`-half and omits canonical root 17.
- **Least-used sidecar.** Lawfully updated root multiplicity, native upset
  size, bank/provenance, and—on the HFC boundary—orientation plus a basepoint.
  Frozen root tags incorrectly kill valid paths.
- **Cheapest new exact probe.** Search for a target-preserving map from the two
  surviving dynamically rooted identity corridors to an oriented boundary
  primitive.  Reject it immediately if bank cell, sheet, or root 17 cannot be
  transported.

### 4. Factorial exact-support lane — finite boundary `r<=9998`; general lane **OPEN**

- **Closest proved mechanism.** [THM-3483, nondivisor residue digit polygons
  and pair compiler](../../01-canon/theorems/THM-3483-factorial-nondivisor-residue-digit-pair-compiler.md)
  upgrades raw Kummer--Legendre hulls with the finite residue `rho_p` at every
  extreme vertex, then uses adaptive primes to close exact-support quadratics
  through `d=10000` (`r<=9998`).
- **Canonical hostile.** At `d=4,p=5`, `rho=0` at the raw constant vertex and
  the actual polygon gains a negative edge.  At `d=6518`, prime 29 is needed,
  refuting every `p<=23` heuristic.
- **Least-used sidecar.** The complete extreme-vertex list and residue triple
  `(m mod p,j mod p,d mod p)`; when `rho=0`, the missing higher valuation and
  root residue remain unpaid.  Endpoint ancestry and positivity are not in the
  pair barcode.
- **Cheapest new exact probe.** Run the adaptive compiler at the first untested
  row `d=10001`, beginning with nondivisor primes through at least 29, and print
  the first unresolved extreme vertex rather than only a pass/fail count.

### 5. Rule 30 — every prize **OPEN**

- **Closest proved mechanism.** [THM-3511, orbit-signalizer gap
  renormalization](../../01-canon/theorems/THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile.md)
  retains the first-return section; [THM-3516, marked van der Put carry](../../01-canon/theorems/THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge.md)
  retains the moving digit and ordinary carry.  [THM-3778, odd-period finite
  scale-cycle no-go](../../01-canon/theorems/THM-3778-rule30-odd-period-finite-scale-cycle-projective-profile-no-go.md)
  linearizes and excludes only the saturated cyclic projective-profile ansatz.
  [THM-3804, all-period amplitude-lattice Smith law](../../01-canon/theorems/THM-3804-rule30-all-period-amplitude-lattice-smith-law.md)
  is now proved and independently audited on this session branch: for
  `n=2^a m`, `d=n/2^min(a,r)` and `s=r-min(a,r)`, the image of `T_n^r` is the
  period-`d` lattice whose block sum is `0 mod 2^s`, with Smith form
  `diag(1^(d-1),2^s,0^(n-d))`.  Together with THM-3778 it excludes physical
  exact finite scale-cycles at every declared period, but no center-column
  prize.
- **Canonical hostile.** Signalizers at scales 7 and 17 have the same
  depth-three portrait but gaps 8 and 5.  At depth six, XOR of visible digits
  misses an ordinary lower carry.  Even period two changes type because
  doubling is not a permutation and the amplitude operator is singular.
- **Least-used sidecar.** Marked child/phase owner, integral amplitude lift,
  ordinary carry, free even-period kernel coordinates, and an honest overflow
  state.  The latest broadcast reports
  a depth-`D+B` bounded-gap observer and a size-three selected-ray bank that is
  not operation-closed; without a routed canon theorem this remains a lead.
- **Cheapest new exact probe.** Close the first 23 selected-ray states under
  one further native section while retaining the integral total-sum carry and
  explicit overflow.  Test the `s_7/s_17` pair before minimizing the bank.

### 6. Tournaments, `H>=disc` — **OPEN**

- **Closest proved mechanism.** [THM-1950, reduction to strong
  tournaments](../../01-canon/theorems/THM-1950-h-ge-disc-reduced-to-strongly-connected.md)
  isolates the base.  [THM-3729, rooted Pfaffian response](../../01-canon/theorems/THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md)
  identifies the exact fixed-root odd energy and proves switching-class mean
  envelopes, not the pointwise all-one-root inequality.
- **Canonical hostile.** The five-vertex tournament in THM-3729 has unsigned
  odd induced-path square sum `272>240`; signed Pfaffian numerators are only
  `(32,128)`.  Mean control loses the fixed root.
- **Least-used sidecar.** The covariant root `u`, transformed with
  `K -> D K D`, plus the signed Pfaffian numerators.  Fixing `u=1` after
  switching is not a quotient invariant.
- **Cheapest new exact probe.** Enumerate strong switching classes by the joint
  `(H,disc,E_odd(1))` profile.  The new exact probe finds the first same-`H`,
  same-`disc`, different-root-response pair at order six; order seven is the
  next classification scale.

### 7. Integer sequences / positive two-cube support

- **Closest proved mechanism.** [THM-2000, support harmonic/Abel--Dini
  surface](../../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md)
  and [THM-2005, support Dirichlet atlas](../../01-canon/theorems/THM-2005-support-dirichlet-automatic-tournament-atlas.md)
  separate support from indexed multiplicity.  [THM-3730, positive distinct
  two-cube support abscissa](../../01-canon/theorems/THM-3730-positive-distinct-two-cube-support-abscissa.md)
  uses the THM-463 divisor split to prove singleton inert-prime rows and
  critical divergence.  [THM-3793, inert primitive-sum singleton](../../01-canon/theorems/THM-3793-inert-prime-sum-all-scale-two-cube-singleton.md)
  is now **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**: if all
  primes dividing `d=x+y` are inert and only the primitive quotient
  `s=d/gcd(x,y)` is cube-free, then the positive distinct fibre is a
  singleton.  Arbitrary inert powers may live in the common scale.
- **Canonical hostile.** `1729=1^3+12^3=9^3+10^3` is the first collision tax.
  The all-inert exponent cutoff fails at three:
  `515375=54^3+71^3=15^3+80^3`, with `54+71=5^3`.
- **Least-used sidecar.** The complete pair sum, its prime exponents, primitive
  gcd, and the good-divisor address—not merely support or a congruence class.
- **Cheapest new exact probe.** Audit the genuinely sharper scale seam omitted
  by the first scout: allow arbitrary inert exponents in the gcd while capping
  only the primitive quotient, and enumerate complete coordinate fibres.  The
  probe below predates that sharpening and corroborates the stricter subcase
  where the cap is imposed on all of `d`, including 57,829 nonprimitive pairs
  through pair sum 1000; it is not a replacement for the theorem's proof.

### 8. Planar source fibres and affine factor surfaces — **OPEN outside named cells**

- **Closest proved mechanism.** The older split/nonsplit atlas closes named
  degrees by response truncation, prefix gcd, infinity poles, and
  `q/h^3` contradictions.  The current chain is sharper: THM-3790 forces a
  noninjective immersed arm and exact normal Bezout law; audited THM-3791
  identifies the universal finite-etale resonant Hensel coefficient;
  THM-3792 closes pure first-normal lifts; THM-3795 closes the entire
  r-independent canonical cell `e^2+z f(e)+z^2h(e)`; and THM-3794 excludes
  surface degree two under constant units.  Audited THM-3797 closes the
  specific confluent quadratic Hermite completion after
  two modification jets, and THM-3799 closes every monomial `r e^m` repair.
  Audited THM-3802 extracts their common actual-chart residue law and audited
  THM-3803 closes all affine-linear `r(b_0+b_1e)` repairs.  THM-3800 corrects
  the proposed sharp-carrier support direction earlier still: the carrier
  itself already has fourteen reduced critical points.  Pending audit,
  THM-3801 supplies the degree-three normalization/companion-sheet gate.
  THM-3796 independently removes one-collision 2x4 cells and common AP steps
  one and two; audited THM-3798 now removes the exact pure common-AP step-three
  cell, and nothing with forgotten nodal anchors.
- **Canonical hostile.** Quadratic and higher nonmonomial `r g(e)` are the live
  carriers left by THM-3795/3799/3803; adding `z^2h(e)` creates a further
  mixed collision grammar.  The corrected THM-3800 sharp torus escape is the
  canonical near miss: it leaves the THM-3795 torus yet has fourteen critical
  points elsewhere.  On the completion side, root collision makes the naive
  hypersurface singular and normalization alone nonfinal (THM-3797), while a
  cubic finite cover can retain an unramified sheet over its branch divisor.
- **Least-used sidecar.** The universal finite-etale root algebra and resonant
  coefficient; the complete affine-modification tree at confluence; and, for
  support arguments, the mandatory `e^2` and `z f(e)` anchors, actual scalar
  bracket bucket, coefficients, and mate support.  The branch lane separately
  needs every ramified and unramified sheet.
- **Cheapest new exact probe.** Retain `r(g_0+g_1e)` as the THM-3803 positive
  control, then shift to the first still-live adjacent quadratic correction
  `r(g_1e+g_2e^2)`.  In both cases retain the minimal nodal core `e^2+f_0z`,
  pair it with a four-term step-three mate, and run the singleton-sign gate.
  The exact mixed-support probe below rejects all eight scalar rows.  The next
  positive test is the general quadratic `(e,u=re)` resultant, not an
  application of the pure AP grammar.

### 9. AMM 12592

- **Closest proved mechanism.** [THM-3340, cyclic single-donor pointwise
  floors](../../01-canon/theorems/THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors.md)
  closes the pointwise problem while the uniform frontier stays open.
  [THM-3648, `R=16384` terminal adjoint recovery](../../01-canon/theorems/THM-3648-amm-r16384-terminal-local-failure-adjoint-golden-phase-recovery.md)
  turns one failed fixed-Rule-A prefix into an exact departure deadline.
- **Canonical hostile.** At the same scale, moving `D0=400` to `D0=412`
  changes the horizon through three coordinates; scale or offset alone is not
  monotone.  A local adjacent failed/closed pair proves no global feasibility
  or value of `C*`.
- **Least-used sidecar.** The exact phase invoice `(m,h,ell)` together with the
  checkpointed adjoint support and an actual alternative continuation after
  departure.
- **Cheapest new exact probe.** Finish the independent optimized
  `R=32768,D0=854` replay routed by CURRENT-FRONTIER, then compute its adjoint
  wall before interpreting the adjacent `854:DIE@8246 / 855:CLOSED` pair.

### 10. Arithmetic-Kakeya forcing certificates — external benchmark **OPEN**

- **Closest proved mechanism.** [THM-2850, paid distinguished-coloop round
  budget](../../01-canon/theorems/THM-2850-paid-distinguished-coloop-round-budget-and-slope-grammar-rank-defect.md)
  identifies forcing with coloop peeling in `[B|D_rho B]`, gives the exact
  cycle/bicirculation defect, and forces `g/u>=1+1/q`.
- **Canonical hostile.** The literal loose rule lets one paid slot generate a
  complete bipartite row bank and drives the score toward one.  The quotient
  `9/5` witness is algebraically sound after identification but has no proved
  recursive same-`H` compilation; label-sharing can cancel forest monomials.
- **Least-used sidecar.** Paid row provenance, prefix/suffix slope grammar,
  loop/parallel multiplicity, and the separate topology, grammar-tie, and
  coefficient-cancellation defects.
- **Cheapest new exact probe.** Exact minimum-seed branch-and-bound at the
  `5/3` rung with two-arboricity filters, while separately proving or refuting
  same-`H` legality for the `9/5` quotient witness.  Never use the unsound
  rule-(1) engine as a positive control.

## Three orthogonal connection audits

### A. Equivariant roots are state, not decoration

This is a shared obstruction pattern, not a literal cross-domain theorem.

```text
tournament source: (K,u)
quotient:           forget u after K -> D K D
preserved:          switching class, disc; on the hostile, H as well
destroyed:          fixed-root odd Pfaffian response
required sidecar:   u -> D u

LRC source:         coefficient-labelled blocker triple
quotient:           sorted valuation orbit
preserved:          valuation profile and admissible minimum set
destroyed:          tied owner
required sidecar:   full coefficient or marked/set-valued owner fibre

FC source:          reset edge history with moving root multiplicity
quotient:           row relation without dynamic root
preserved:          one-edge identity/magnitude data
destroyed:          chronological congruence
required sidecar:   lawful root update plus bank/upset provenance
```

The exact probe exhausts labelled tournaments through order six and finds a
minimal same-`H`, same-`disc` separation at `n=6`.  Two strong switching-
equivalent tournaments both have `H=29` and `disc=2`, but fixed-root odd
energies `6` and `4`.  The switch is
`D=(1,1,-1,1,1,-1)`, and the covariant identity
`E_odd(K,D1)=E_odd(DKD,1)=4` passes.  A second implementation reconstructs
the energy from all odd-subset Pfaffian squares.

**Consequence:** any next-operation observer in these lanes must transport its
root/owner; a switching mean, valuation quotient, or one-edge relation is not
a congruence by itself.

### B. A two-cube value is a reversible address on one LRC Graver subatlas

This one is a genuine bounded map:

```text
source:      coprime 1<=a<b, a+b<=356 (THM-3743 support-two ratio)
restriction: every p|(a+b) has p=2 mod 3 and v_p(a+b)<=2
map:         (a,b) -> a^3+b^3
preserved:   ordered ratio, l1 mass, centred Christoffel address
destroyed:   scale, coordinate placement, eleven other speeds,
             phase, owner, arrival, loneliness
sidecars:    labelled runner pair, scale, full row, semantic packet
```

The independent complete-fibre probe reproduces:

- all `19,314` support-two ratios;
- `94` admissible sums and `5,855` admissible ratio addresses;
- `456,690` labelled coordinate-pair addresses;
- no collision for those values in the complete positive-pair universe up to
  their maximum `44,738,876`; and
- no failure among `57,829` nonprimitive admissible pairs through pair sum
  `1000`.

The same atlas has 28 collision values without the prime restriction, starting
at `1729`.  The exponent-three and split-prime hostiles are active.  This is
**FINITE-EXACT** evidence and a useful index.  It excludes no LRC row.  Because
the probe caps exponents on the whole sum rather than only the primitive
quotient, it is a bounded corroboration of a strict subcase of now-proved
THM-3793, not the source of that theorem and not an LRC consequence.

### C. Surface completion exposes distinct jet, anchor, and sheet coordinates

The THM-3791--3803 chain is best read as a sequence of typed sidecars, not one
scalar obstruction.

1. **Finite-etale Hensel resonance.**  For `n=3` and
   `Sigma(c,b)=b(b-1-lambda c^2)`, the `lambda=0,1` surfaces have identical
   root jets through order one, while their resonant vectors are `(0,0)` and
   `(0,1)`.  The exact Poisson-centrality and Hensel equations pass.  This is
   now a split positive control for audited THM-3791: the theorem's promotion
   makes the root jet intrinsic in the universal finite-etale algebra
   `E=k[b]/(Sigma_0)`, so no chosen splitting is part of the invariant.
2. **Confluence is a different type.**  THM-3797 began this session as a
   reserved confluent stub and was promoted, then independently audited, on
   incoming `origin/main`.  Its root
   collision is outside THM-3791's squarefree finite-etale hypothesis.  The
   singular hypersurface and first normalization are not final: adjoining
   two jets yields the smooth three-arm exact image, where the residue vector
   `(1,0,0)` modulo constants survives.  Thus the connection preserves the
   eventual moving-axis residue but loses finite-etale descent and requires
   the entire affine-modification tree; it is not a formal application of
   THM-3791.  Audited THM-3802 now proves that the shared
   residue calculation *does* extend to any supplied smooth affine
   plane-chart contact-tree atlas; crucially, it assumes the actual charts,
   inverses, overlaps, and separated surface rather than manufacturing them
   from a root diagram.
3. **R-independent normal corrections are now closed.**  For the canonical
   nodal carrier on `r^2e-z^3+r=0`, every `alpha!=0` in
   `A=e^2-z/3+alpha z^2` has the explicit critical point
   `e=0, z=1/(6alpha), r=1/(216alpha^3)`.  At `alpha=0`, the seven-point
   equation `8z^7+9=0` returns.  This probe is no longer an open THM-3792
   extension: it is a positive control inside audited THM-3795, which closes
   every r-independent `e^2+z f(e)+z^2h(e)`.
4. **Branch sheets do not extrapolate quadratically.**  The cubic cover
   `z^3-3tz+2=0` has discriminant `108(t^3-1)` and fibre
   `(z-1)^2(z+2)` at `t=1`.  The simple sheet `z=-2` has derivative `9`, so
   it remains etale above the branch divisor.  This is an exact hostile to
   extending THM-3794's “delete the whole branch fibre, hence obtain a unit”
   step from degree two to degree three.  Pending-audit THM-3801 now gives the
   honest degree-three replacement under the constant-unit surface
   hypotheses: the finite flat normalization is nonmonogenic, every
   codimension-one branch has `(2,1)` type, and one reduced companion sheet is
   visible.  The toy monogenic cover illustrates the sheet loss but is not a
   proof of that theorem.

The newest origin objects separate the next operations.  Corrected and
audited THM-3800 supersedes the proposed `{-6,1}` mate grammar because its
sharp torus-escaping carrier already has fourteen critical points.  THM-3801
proves the cubic normalization/companion-sheet gate pending audit.  Audited
THM-3802 proves the contact-tree de Rham abstraction suggested by
THM-3791/3797, and audited THM-3803 proves the affine-linear
`r(g_0+g_1e)` resultant.  Reserved THM-3805 names the next general quadratic
resultant.  The probes here are independent controls, not dependencies for
those theorems.

#### Honest THM-3795 to THM-3798 support map

For a live correction

```text
r g(e)=sum_j g_j r e^j,
```

the Euler weights are `3-3j`, a step-three progression.  Projecting a
two-term `g_0+g_1e` to its correction support therefore gives `{3,0}`, which
looks exactly like the two-term side of the pure common-AP 2x4 step-three cell
proved in THM-3798.

That projection is not a reduction.  A legal nodal carrier also contains the
irregular anchors `e^2` and `f_0z`, of weights `-6` and `1`.  The map

```text
source:       e^2 + f_0 z + r(g_0+g_1e), with a four-term step-3 mate
projection:   forget e^2 and f_0 z
target:       pure correction 2x4 common-AP step-3 support
preserved:    correction coefficient support and its step-three spacing
destroyed:    actual carrier support, bracket-bucket multiplicities,
              scalar placement, arm value and the Bezout datum f(0)
needed sidecar: both anchors, all coefficients, r-coupling, mate support,
                and the scalar bracket bucket
```

is therefore controlled forgetting without the needed sidecar.  THM-3798 is
now **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**, but its exact
hypothesis is the pure 2x4/4x2 cell.  It cannot consume the forgotten anchors
and therefore cannot be cited as a mixed-carrier reduction.

The cheapest exact mixed probe restores precisely those two anchors.  The
affine-linear control's
pure projected bucket histogram is

```text
[(0,1),(3,2),(6,2),(9,2),(12,1)],
```

whereas the actual carrier support `{-6,0,1,3}` against the four-term mate
has histogram

```text
[(-6,1),(-3,1),(0,2),(1,1),(3,3),(4,1),
 (6,2),(7,1),(9,2),(10,1),(12,1)].
```

The four possible scalar rows `(q,b)=(0,-2),(3,-5),(6,-8),(9,-11)` all fail
the proved singleton opposite-sign gate.  Their respective hostile singleton
contributors are `e^2[-6]` paired with mate weight `1`, then `f_0z[1]`
paired with mate weights `-5,-8,-11`.  This closes only that exact minimal
mixed 4x4 support configuration.  It neither follows from nor reproves
THM-3798, does not reduce the mixed lane to it, and does not close larger
nonmonomial `r g(e)` or mixed
`z^2h(e)` corrections; rather, it is an independent anchor-restored control
next to the proved pure theorem.  Incoming, independently audited THM-3799 closes every
one-monomial `r e^m` correction, and audited THM-3803 independently closes the
full affine-linear case.

The refreshed probe therefore also shifts the correction window once, to the
cheapest still-live adjacent quadratic pair `r(g_1e+g_2e^2)`.  Its projected
histogram is

```text
[(-3,1),(0,2),(3,2),(6,2),(9,1)],
```

while restoring the anchors gives support `{-6,-3,0,1}` and histogram

```text
[(-6,1),(-3,2),(0,3),(1,1),(3,3),
 (4,1),(6,2),(7,1),(9,1),(10,1)].
```

The scalar rows `(q,b)=(-3,1),(0,-2),(3,-5),(6,-8)` also all fail: the first
has the singleton `e^2[-6]` against mate weight `1`, and the remaining three
have `f_0z[1]` against mate weights `-2,-5,-8`.  Thus all eight rows across
the affine and adjacent-quadratic controls fail the sign gate.  The quadratic
result is only this exact 4x4 mate-support exclusion; a general quadratic
carrier critical-point theorem and larger mixed corrections remain open.

**Consequence:** the next planar-fibre search should keep
`finite-etale resonant jet -> confluence modification tree -> anchored
canonical support -> critical resultant -> branch-sheet/unit` as separate
operation columns.  None implies `JC(2)` or `DC(2)`.

## Prioritized next pulls

1. **LRC anchor:** intersect the `l1<=356` Graver atlas with every rank-eleven
   star space.  The 5,855 reversible two-cube addresses are a ready indexing
   sidecar for the support-two portion; retain the actual labelled pair and
   run the AP control.
2. **Planar-fibre niche:** after the minimal anchored support gate, compute
   the exact critical resultant for `e^2-z/3+r(g_1e+g_2e^2)`, then pass to a
   general quadratic `g` and add the smallest `z^2h(e)` term.  Keep THM-3794's
   unit hypothesis and degree-two scope, THM-3798's pure-support scope,
   THM-3805's reserved status, and THM-3801's pending audit explicit.
3. **Operation-closure wildcard:** use the order-six rooted tournament hostile
   as the regression test for any proposed root-free observer; then apply the
   same stabilizer audit separately to the Rule 30 selected-ray bank and the
   FC rooted corridor.  In Rule 30, carry THM-3804's Smith free part and cyclic
   `2^s` quotient rather than only a projective profile.
4. **Independent finite closures:** finish AMM `R=32768,D0=854`, factorial
   `d=10001`, and the Arithmetic-Kakeya `5/3` min-seed search.  These are
   orthogonal exact jobs and should not be rhetorically merged.

## Honest remaining frontier

LRC(14), `JC(2)`, `DC(2)`, `HFC(3)`, `FC(3)`, `H>=disc`, every Rule 30 prize,
the full positive two-cube support asymptotic/residue, arbitrary planar
fibres, the uniform AMM constant, and the Arithmetic-Kakeya `<=1.675` target
all remain open.  The refresh found no legal reduction among those problems.
It found one reversible bounded address, one exact cross-lane obstruction
pattern, two scoped mixed anchored-support exclusions, and two planar
boundary controls.  The `z^2` control is now subsumed by THM-3795 rather than
a new frontier result.

## Reproduction

Run from the repository root:

```powershell
python3 -B .scratch/all_frontier_refresh_20260823/rooted_switching_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/rooted_switching_probe.py
python3 -B .scratch/all_frontier_refresh_20260823/lrc_two_cube_address_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/lrc_two_cube_address_probe.py
python3 -B .scratch/all_frontier_refresh_20260823/planar_surface_sidecar_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/planar_surface_sidecar_probe.py
python3 -B .scratch/all_frontier_refresh_20260823/mixed_anchored_step3_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/mixed_anchored_step3_probe.py
```

Every normal/optimized pair matched, and every stream matched its frozen
`.out` file.  The rooted probe also uses two independent exact energy paths
(augmented determinants and odd-subset Pfaffian squares).  Root independently
replayed all eight normal/optimized stream comparisons and confirmed zero AST
assertions in the four companions.
