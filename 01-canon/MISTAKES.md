# Mistakes Log

**Purpose:** preserve failed implications, minimal witnesses, and repaired
statements so future agents do not repeat them. Search this ledger before using
a historical synthesis or promoting computational evidence.

Format per entry:
- What was assumed / done
- Why it was wrong
- The correct framing

## MISTAKE-441 (2026-08-22, first THM-3684 draft) -- three named generator subalgebras were broadened to every one-direction output

- **What failed:** the proof excludes a Jacobian mate for outputs in
  `C[b]`, `C[c]`, or `C[e]`, but its first concluding slogan said that each
  output must mix two target directions.
- **Minimal witness / first failed implication:** the Laurent argument uses
  the pure weights of exactly `b`, `c`, and `e`.  It says nothing about a
  different one-generated subalgebra such as `C[b+c]` or the stable algebra
  `C[w]`; calling all of these one "direction" erased the generator type.
- **Repair / strongest survivor:** THM-3684 now states precisely that neither
  output may lie in `C[b]`, `C[c]`, or `C[e]`, and arbitrary mixing only in
  the opposite output cannot repair that slice.  The master Laurent formula,
  all three leading coefficients, the univariate-composition extension, and
  the full-source no-mate result passed independent hostile audit unchanged.
- **Reusable rule:** a grading obstruction for named homogeneous generators
  transfers to their univariate polynomial algebras, not automatically to
  every rank-one or one-generated subalgebra.

## MISTAKE-440 (2026-08-21, first THM-3682 draft) -- moving a current-scale pure-word role inside a target-fixed delayed word

- **What failed:** the first role-to-target leak draft wrote the original
  THM-2365 successor row as `S(s,0)=int w_s(y)d(y-s/13)dy` by treating the
  redundant danger factor of a pure terminal word as a present factor.
- **Minimal witness / first failed implication:** the actual handoff has
  `F_(s,0)(x)=E_(s,0)(x)Q(Rx)` with `13|R`.  Its word factor transforms as
  `d(Rc_a x-Rs/13)=d(Rc_a x)`, so the target action fixes it exactly.  The
  pushforward contains `d(Ry)`, not `d(y-s/13)`.  Thus the forced
  factor-coloured role mode does not type the original delayed target row.
- **Repair / strongest survivor:** the normalized Fourier reverse-triangle
  tariff, `121 rho^2/2028` role floor, sharp cancellation hostile, and
  `1/2028`, `1/26364` energy invoices are correct.  THM-3682 now realizes
  them on an explicitly constructed current-scale present-`Q` auxiliary
  table.  Transport to the original handoff requires a target/word-action
  intertwiner or Bockstein sidecar.  The reusable warning is: a factor may be
  movable before transport but neutral after a clock divisible by the target
  modulus; always type the scale before composing character actions.

## MISTAKE-439 (2026-08-21, THM-3676 audit) -- disjointness from the blind space was mislabeled as the sharp detector premise

- **What failed:** the first THM-3676 draft called
  `supp(mu) intersection W=empty` the sharp cross-source phase-cone
  alternative.
- **Minimal witness / first failed implication:** a positive measure such as
  `delta_0+delta_(e_j)` meets the nonzero-sum hostile `{0}` but has mass
  outside it and is detected by the other source family.  Likewise
  `delta_1+delta_(e_j)` meets the zero-sum scalar line without being confined
  to it.
- **Repair / strongest survivor:** under a common strict phase cone, the
  exact premise is only `supp(mu) not subset W`, where `W={0}` or
  `W=span(1)`.  The cross-source intersections and eight-chart theorem were
  correct.  Source switching was also retyped as two packet families over
  one fixed labelled residue kernel.
- **Reusable rule:** when a contradiction says every support point lies in
  a hostile subspace, negate containment, not intersection.

## MISTAKE-438 (2026-08-21, THM-3674 audit) -- the zero three-site mask was left inside a Rayleigh quotient

- **What failed:** the first THM-3674 draft stated
  `Var(S)>=|Delta|^2/(kappa p^2)` for every pair of displacements, while its
  own collision table included `kappa=0` at `a=b=0`.
- **Minimal witness / first failed implication:** when `a=b=0`, the combined
  mask and defect are both zero.  The displayed right side is `0/0`, not a
  defined zero bound.
- **Repair / strongest survivor:** require the combined mask to be nonzero
  before dividing by its squared norm; the zero-mask case is vacuous.  The
  exact drift decomposition, sharp lawful constant, all nonzero-mask equality
  cases, and the distinct-site LRC tariff passed exhaustive hostile checks.
- **Reusable rule:** every norm or Rayleigh quotient needs an explicit
  nonzero-vector domain, even when the degenerate numerator vanishes too.

## MISTAKE-437 (2026-08-21, THM-3671 audit) -- termwise unit support was conflated with a nonzero grouped phase-cone measure

- **What failed:** the first THM-3671 draft described the blind-space theorem
  as ranging over all owner-pivot packets without fixing the selected source,
  allowed a formally zero measure, used an ambiguous “open half-plane,” and
  treated THM-2325/2331's termwise all-unit address bank as the grouped
  residue-measure support needed by the pushforward argument.
- **Minimal witness / first failed implication:** changing the selected source
  changes the blind source line, so two such fixed-source intersections may
  already meet only at zero.  The zero measure has every pushforward zero.
  Also `1,-1` lie in an affine open half-plane such as `Im(z)<1` but cancel;
  only a strict half-plane through the origin forbids cancellation.
- **Repair / strongest survivor:** for one fixed selected source and target
  pair, the exact all-packet intersection and four-chart reduction pass.  For
  a nonzero grouped measure whose nonzero coefficients satisfy
  `Re(exp(-i theta)z)>0`, all-coordinate-unit support forces a nonzero target
  pushforward off the scalar-sum-zero branch.  Neither nonzero grouped support
  nor this common phase cone is currently known for the LRC current.
- **Reusable rule:** a term bank is not its grouped coefficient measure, and
  a geometric cone prevents cancellation only when it is pointed at the
  origin.  State the grouping map, nonzero premise, and fixed ambient stratum
  before applying support-intersection arguments.

## MISTAKE-436 (2026-08-21, THM-3670 audit) -- a raw boundary jump was substituted for a signed successor-mass defect

- **What failed:** the first THM-3670 draft advertised its successor-mass
  gate for every owner-pivot packet and said that any boundary-jump
  discrepancy obstructed the rigid thirteen-number configuration.
- **Minimal witness / first failed implication:** the THM-2365 successor
  blocker must be one of the two target blockers; selecting it as the source
  leaves the theorem's two-dimensional target quotient.  Independently,
  `f_0=1_[0,1/2)`, `f_A=1_[1/2,1)`, and `f_B=1_[1/4,3/4)` have a nontrivial
  boundary-jump combination but
  `integral(f_0+f_A-2f_B)=0`.  Jumps can cancel after signed integration.
- **Repair / strongest survivor:** restrict to THM-2365 strata whose
  successor blocker lies in the target pair.  A boundary or valuation
  calculation is sufficient only if it proves a nonzero signed successor
  integral or directly contradicts the thirteen-number equations.  The core
  implication, the rank-11 system on thirteen masses, and the resulting
  two-dimensional recirculation locus all passed independent hostile audit.
- **Reusable rule:** when a detector is defined by an integral, local support
  or jump asymmetry is not itself the detector.  Push the local witness
  through every projection and signed integration before claiming survival.

## MISTAKE-435 (2026-08-21, THM-3667 audit) -- a frequency-dependent phase was mistaken for spectral conjugacy

- **What failed:** the first THM-3667 draft correctly proved that swapping the
  two positive weights preserves every Fourier-multiplier magnitude, then
  incorrectly transferred the thirteen eigenvalue collisions from one
  maximin orientation to the other and advertised a conditioning-versus-
  rigidity tradeoff.
- **Minimal witness / first failed implication:** the exact identity is
  `lambda_(1-a)(-u,v-u)=zeta^u lambda_a(u,v)`.  The factor `zeta^u` depends on
  the frequency.  It preserves singular values but not equality of complex
  eigenvalues and is not a conjugacy; already the two traces are `169a` and
  `169(1-a)`.  An exact `Phi_13` reduction gives 143 singletons plus thirteen
  doubles for `(a,b)=(a_*,b_*)`, but 169 singletons for the swapped optimizer.
- **Repair / strongest survivor:** the maximin weights, frame constants and
  magnitude equality cases survive.  The stronger repaired result is that
  optimal conditioning and full simple-spectrum rigidity coexist in the
  swapped orientation.  The exceptional real collision weights in `0<b<1`
  form the finite chord-ratio set
  `{sin(r*pi/13)/sin(s*pi/13):1<=r<s<=6}`.
- **Reusable rule:** an identity of the form `lambda'(phi(q))=c(q)lambda(q)`
  transports norm/frame data only.  Before transferring eigenvalue
  multiplicities, centralizers or conjugacy, prove that `c(q)` is constant or
  perform a second exact spectrum census.

## MISTAKE-434 (2026-08-21, Cohn repair-cycle transfer) -- parity was substituted for multiplier holonomy

- **What failed:** factorial rescaling turns the interior recurrence in the
  Cohn repair path into an alternating `d_i=-d_(i-1)` chain.  The first draft
  then closed this flattened chain cyclically and treated the usual odd/even
  kernel of `I+S_n` as if it were a cyclic Cohn repair.
- **Minimal witness / first failed implication:** rescaling is not periodic.
  For weighted cyclic equations `alpha_i c_i+c_(i-1)=0`, the determinant is
  `product(alpha_i)-(-1)^n`.  The Cohn seam retains weights
  `2,3,...,n+1`, so its determinant is `(n+1)!-(-1)^n`, nonzero for every
  `n`; even support does not close it.
- **Repair / strongest survivor:** odd/even rectangle rigidity remains a
  valid support filter.  A coefficient cycle must additionally carry signed
  multiplier holonomy `(-1)^n`; the smallest even rational positive control
  has reciprocal gains `2,1/2`.  THM-3653 and its companion separate these
  two predicates explicitly.
- **Reusable rule:** a gauge that flattens a transport equation on a path need
  not descend through a quotient seam.  Before cyclically identifying support
  endpoints, retain edge weights and test the exact product holonomy.

## MISTAKE-433 (2026-08-21, THM-3639 namespace race) -- a pre-push ID check did not survive a concurrent reservation

- **What failed:** the Berggren denominator-401 package checked that
  `THM-3639` was free and pushed a reservation, but another session reserved
  the same YAML identifier under a different filename during the subsequent
  fetch/rebase/push window.  Git could merge the two paths, so the mathematical
  files landed with a duplicate theorem ID.
- **Minimal witness / first failed implication:** commit `a25a2ee6d` reserved
  `THM-3639` for the Russell-cylinder universal stable debt before the cube
  reservation commit `34289e968` reached `origin/main`.  Filename uniqueness
  therefore did not imply YAML-ID uniqueness on the rebased tree.
- **Repair / strongest survivor:** the complete cube theorem, companion, and
  transcript were moved to the then-free `THM-3640`; the three obsolete cube
  paths carrying `3639` were removed in commit `4d49fa05a`.  Semantic digest
  `879a912b...c1836374` and every mathematical count are unchanged.  The
  Russell-cylinder `THM-3639` remains the sole owner of that identifier.
- **Reusable rule:** reserve scarce IDs in a tiny checkpoint and, after every
  rebase, recheck both filenames and YAML IDs on the exact candidate tree
  immediately before committing.  If the remote acquired the ID, renumber
  before pushing; an earlier availability check is not a concurrency lock.
- **Same-day recurrence / sparse repair:** the same race recurred while
  reserving the AMM threshold as `THM-3642`: root commit `b645cbbcd` acquired
  that ID between checks, while the already-started push of `dac6e6da8`
  completed despite a late interrupt.  The first renumber checkpoint added
  `THM-3644` but, because the stale path was skip-worktree in the sparse
  checkout, did not delete the duplicate.  After explicitly materializing the
  path, commit `fe8a7a4fe` removed it; the Russell-cylinder theorem is again
  the sole `THM-3642`.  Strengthened rule: after interrupting any push, fetch
  and inspect `origin/main` rather than assuming cancellation; in sparse
  worktrees, verify deletions in the index and finish with `git grep` for the
  YAML ID on the pushed tree.

## MISTAKE-432 (2026-08-21, THM-3592 support classification) -- simultaneous reversal was not a regularity symmetry

- **What failed:** the first THM-3592 proof classified three-point sumsets
  modulo simultaneous reversal, but audited only the canonical hook and
  Euclidean representatives.  Although reversal preserves the additive fibre
  pattern, it changes the absolute weights and therefore the divisibility
  floor `Sigma^ceil(-u/2)|f`; it is not a symmetry of Danielewski regularity.
- **Minimal witness / first failed implication:** the omitted hook
  `(d,d;2d,d)` has fibre word
  `00;10;01+20;02+11;12+21;22`.  Its `21=(1,-2)` arm gate survives every
  singleton test for all `d>=3`, so the canonical-hook census cannot be
  transported through reversal.
- **Repair / strongest survivor:** THM-3592 now audits the reflected hook and
  both reflected Euclidean families separately.  The hook survivor is killed
  by the exact arm-order mismatch `(2d-1)ell != (2d+2)ell-1`; the reflected
  Euclidean rows die by rectangle and singleton gates.  The companion now
  includes all three orientations and the universal seven-piece support floor
  survives, pending a second independent hostile audit.
- **Reusable rule:** quotient a finite support atlas only by symmetries that
  preserve every downstream side condition.  If the target uses signed
  weights, divisibility floors, or boundary regularity, replay every
  reversal-sensitive orbit representative before claiming exhaustiveness.

## MISTAKE-431 (2026-08-21, THM-3579 degree-operator lemma) -- the constant/constant zero boundary was omitted

- **What failed:** the first statement of THM-3579's degree-rigid operator
  lemma quantified over all nonzero polynomials `H,u` and then asserted
  `deg E_H(u)=deg J_H(u)=deg H+deg u-1` with nonzero leading coefficient.
- **Minimal witness / first failed implication:** `H=u=1` gives
  `E_H(u)=J_H(u)=0`, while the displayed degree is `-1` and there is no
  nonzero leading coefficient.
- **Repair / strongest survivor:** require `deg H+deg u>0`.  Every use in the
  proof has `deg H>=deg Sigma>=2` (or a nonzero arm-divisible companion), so
  all five scalar-row exclusions and the equal-step `3x3` theorem are
  unchanged.  An independent audit accepted every load-bearing valuation,
  support, first-integral, and coprimality step after this one-line repair.
- **Reusable rule:** before applying a leading-term degree formula for a
  differential operator, split off the constant/constant zero image even
  when later applications automatically have positive degree.

## MISTAKE-430 (2026-08-20, THM-3569 promotion suffix) -- the audited universal theorem retained its provisional scope block

- **What failed:** after THM-3569 was promoted to an independently audited
  theorem for every squarefree `Sigma` of degree at least two, Section 5
  still labeled that universal statement `PENDING AUDIT` and called only the
  old specialization `Sigma=b(b+4)` proved.
- **Repair / strongest survivor:** replace the stale block by the proved
  universal scope.  The support-collision proof, exact companion, hashes,
  open `3x3`/`2x4` boundary, and all mathematical claims are unchanged.
- **Reusable rule:** promotion audits must search scope tables and fenced
  status summaries as well as frontmatter and final paragraphs; stale truth
  labels can survive inside a theorem even when its proof was correctly
  generalized.

## MISTAKE-429 (2026-08-20, THM-3574 certificate bookkeeping) -- the intermediate quotient inherited the final quotient's term count

- **What failed:** the first THM-3574 proof said that the intermediate
  core-cofactor quotient `B_H` had 152 expanded terms.  An independent
  recomputation found 54; 152 is the term count of the final Bezout quotient
  `C_H=V(F)B_H/4`.
- **Repair / strongest survivor:** the companion now gates all three counts:
  8 for the linear-unit quotient, 54 for the core-cofactor quotient, and 152
  for the final inverse certificate.  Both divisibilities and the universal
  nonconstant-unit theorem are unchanged.  The proof also now spells out the
  smoothness, dominance, and generic-cubic argument that sends source
  reducibility into THM-3573's Pell compiler.
- **Reusable rule:** when one certificate is obtained by multiplying another,
  record and verify their sizes separately; a valid divisibility does not
  make its intermediate and final expanded quotients interchangeable.

## MISTAKE-428 (2026-08-20, THM-3571 irreducibility case split) -- vanishing quadratic `a` coefficients did not force `deg_a(phi)=0`

- **What failed:** the first THM-3571 proof said that `A=B=0` in the
  collision-compatible quadratic
  `phi=Aa^2+Bab+Cb^2+Da+Eb+D/4-A/16` forces `deg_a(phi)=0`.  It forgot the
  surviving linear term `Da`.
- **Minimal witness / first failed implication:**
  `phi=a+b^2+1/4` has `A=B=0`, is genuinely quadratic, contains the collision
  value, and has `deg_a(phi)=1`, not zero.
- **Repair / strongest survivor:** the Euler and discriminant proof was
  unaffected.  For irreducibility, split the residual cell by `D`.  If
  `D!=0`, THM-3565 would require `-2h(b)^3=D`; hence `h` is constant and the
  entire reducible row is affine, contradicting the `b^2` term.  If `D=0`,
  THM-3564 applies at `deg_a(phi)=0`.  More generally, when `A=0` and
  `(B,D)!=(0,0)`, comparison of `Bb+D=-2h^3` again forces `h` constant and
  excludes every genuinely quadratic row.  THM-3571 survives unchanged.
- **Reusable rule:** determine a partial degree from every coefficient of
  that variable, including lower-degree terms; a vanished top homogeneous
  block does not by itself put a polynomial on the degree-zero face.

## MISTAKE-427 (2026-08-18, degree-passport scope) -- the first survivor of a weaker sieve was called globally admissible

- **What failed:** the first synthesis called the four height-`105` cells
  surviving the Appelgate--Onishi/Nagata, Moh, and 2014 gcd filters the first
  “globally admissible” degree cells.  That wording omitted the stronger
  cited 2022 sub-`125` degree-pair classification.
- **Minimal witness / first failed implication:** the 2022 list leaves only
  `(72,108)` (and its transpose) below height `125`, so every height-`105`
  cell in the weaker screen is excluded.
- **Repair / strongest survivor:** retain the height-`105` enumeration only
  as a finite hostile validating the weaker filter stack.  The current cited
  passport has first reduced pair `(72,108)`.  Its consequences move the
  row-dominant torus boundary from `(84,105)` to the first not-yet-excluded
  `(100,125)` box and move the width-three mixed-Catalan coefficient cap from
  `102` to at least `105`.
- **Reusable rule:** every “first admissible” computation must name the exact
  inherited filter set.  When a stronger classification is imported, keep
  the old frontier only as a labeled hostile, never as the live global floor.

## MISTAKE-426 (2026-08-18, THM-3555 status suffix) -- a second promoted proof retained its reserved-stub disclaimer

- **What failed:** THM-3555 was promoted to `PROVED + VERIFIED-EXACT` with a
  full affine-equivalence proof and companion, but the old reservation suffix
  still claimed that no theorem, proof, or dependency was asserted.
- **Repair / strongest survivor:** remove the stale suffix.  The universal
  cubic marked-root normal form, discriminant calculation, and fixed-
  ramification-line no-go are unchanged.
- **Reusable rule:** the promotion audit in MISTAKE-425 is now triggered by a
  second independent instance; search the full body for reservation language
  before every status promotion, including the final paragraph after scope.

## MISTAKE-425 (2026-08-18, THM-3554 status suffix) -- a promoted proof retained its reserved-stub disclaimer

- **What failed:** THM-3554 was promoted to `PROVED + VERIFIED-EXACT` with a
  complete proof and replay, but its former reservation paragraph remained at
  the very end and incorrectly said that no theorem, proof, or dependency was
  asserted.
- **Repair / strongest survivor:** remove only the stale reservation suffix.
  The proved punctured-Kummer normal form, its exact quadratic filling
  obstruction, and all stated scope boundaries are unchanged.
- **Reusable rule:** when promoting a reserved theorem, search the entire body
  for `RESERVED`, `stub`, and old no-claim language; frontmatter promotion alone
  does not make the truth surface internally consistent.
## MISTAKE-424 (2026-08-18, THM-3554 image typing) -- a punctured Kummer cover's geometric image was called its untyped scheme/rational image

- **What failed:** the first promoted THM-3554 said without qualification
  that the image of the collision surface in the target plane was exactly
  `U=D(beta^2-16alpha)`.  This conflated the morphism's topological/geometric-
  point image with its scheme-theoretic image in `A^2` and its image on
  `k`-rational points.
- **Minimal witness / first failed implication:** over `Q`, the target point
  `(alpha,beta)=(-1/8,0)` lies in `U(Q)` because `delta=2`, but a preimage
  would require `s^2=delta/4=1/2`.  Thus it is not in the rational-point
  image.  Conversely, the scheme-theoretic image in the ambient `A^2` is all
  of `A^2`, because the topological image `U` is dense.
- **Repair / strongest survivor:** `S->U` is finite etale and surjective, so
  its topological image and its image on algebraically closed points are
  exactly `U`.  State the ambient scheme-theoretic closure and rational-point
  square-class qualification separately.  The Kummer normal form, deck
  collision, and no-etale-filling valuation proof are unchanged; the latter
  in fact holds over every field of characteristic different from two.
- **Reusable rule:** for a nonproper or punctured morphism, always type
  “image” as topological, scheme-theoretic, geometric-point, or rational-
  point image.  These need not agree.

## MISTAKE-423 (2026-08-18, provisional THM-3546 scope) -- rational-point graph containment was silently treated as a polynomial identity

- **What failed:** the first THM-3546 wording said that `F(Graph(h))` lies in
  `Graph(g)` over an arbitrary field, then inferred that the straightened
  normal component is divisible by the graph variable.  This is sound for
  graph-subscheme containment, but set-theoretic containment on `k`-rational
  points is weaker over finite fields.
- **Minimal witness / first failed implication:** over `F_q`,
  `F(x,z)=(z,x^q-x)` has determinant one and sends every rational point of
  `z=0` to `(0,0)`, yet its straightened normal polynomial is `x^q-x`, not a
  multiple of `z`; the induced rational-point map is constant.
- **Repair / strongest survivor:** require
  `F_w(x,h(x))=g(F_y(x,h(x)))` in `k[x]`, equivalently scheme-theoretic graph
  containment or containment on all geometric points.  Then the normal
  component is `rA`, the determinant factors as `det Jf * A(x,0)`, and both
  factors are units.  Over an infinite field, containment on all rational
  points already implies the identity.
- **Reusable rule:** over finite fields, equality on rational points is not a
  polynomial identity.  State subscheme containment or audit the vanishing
  ideal before using divisibility.

## MISTAKE-422 (2026-08-18, THM-3025 notation collision) -- the subleading Jacobian was mistaken for the Keller constant

- **What failed:** THM-3025 relabelled the multiplier in THM-3016's
  rigidity identity as `J = Jac(P,Q) in k*`.  In THM-3016 the displayed
  multiplier is instead the homogeneous subleading form
  `j = Jac(H,Q_(m-1))`, which may vanish.  Thus the original step
  `j Jac(W,H)=0 => Jac(W,H)=0` was unjustified, and the companion script's
  opening prose repeated the same collision.
- **First failed implication:** the zero branch `j=0` makes the product
  identity tautological, so it supplies no information about `W` by division.
  This is a proof-text error rather than a counterexample to the conclusion.
- **Repair / strongest survivor:** split on `j`.  If `j!=0`, THM-3016's
  common-form and coprime-degree argument gives `W=0` whenever `H` has at
  least two distinct roots.  If `j=0`, then `Q_(m-1)!=0` would make `H` and
  `Q_(m-1)` powers of a common form of degree dividing both `g` and `m-1`;
  since `g|m`, that degree is one and `H` has only one root.  Hence for
  `K>=2`, `Q_(m-1)=0`.  The subleading Keller equation then gives
  `Jac(P_(n-1),H)=0`, and the same `gcd(g,n-1)=1` argument forces
  `P_(n-1)=0`.  Therefore `W=0` survives in both branches.  THM-3025 and its
  replay now use distinct names `J0=Jac(P,Q)` and
  `j=Jac(H,Q_(m-1))`; no division by `j` occurs in the zero branch.
- **Reusable rule:** never reuse the full-Jacobian symbol for a homogeneous
  coefficient or response row.  Before dividing a named Jacobian, expand its
  arguments and degree and audit the zero branch explicitly.

## MISTAKE-421 (2026-08-16, THM-3540 proof gate) -- etaleness was used as a coordinate-projection derivative test

- **What failed:** the first THM-3540 proof said the inverse `x`-core
  derivative `d=E'(X)` was nonzero because the predecessor cover is etale.
  Etaleness separates full source points, but two distinct fibre points may
  still share one displayed coordinate, so it does not by itself make an
  `x`-eliminant squarefree at that coordinate.
- **Minimal witness / first failed implication:** any etale finite fibre with
  two points having the same `x`-projection defeats the implication.  No
  THM-3540 conclusion fails: its independent normalization-axis calculation
  gives `d|_(lambda=0)=4` exactly, proving that `d` is a nonzero rational
  function before the discriminant square-class division.
- **Repair / reusable rule:** distinguish separability of the full finite
  algebra from separability of a chosen coordinate eliminant.  Gate division
  by a coordinate derivative with a primitive-coordinate theorem or one
  exact nonidentity witness; here the axis witness supplies the latter.

## MISTAKE-420 (2026-08-16, provisional THM-3534 scope) -- relative rank-two H1 collapse was overextended past the endpoint-retaining coefficient space

- **What failed:** the provisional THM-3534 draft correctly computed
  `dim H^1(C4;V_rel,C)=1` for the two-dimensional relative response plane,
  then said the frozen data supplied no two-dimensional chamber-faithful
  cohomology.  That conclusion silently fixed the coefficient-space
  dimension while quantifying over all response spaces.
- **Minimal witness / first failed implication:** the already computed
  three-space `V_ext=S+I` is chamber-stable and satisfies
  `V_ext=I direct_sum E_6=S direct_sum E_6`.  Chamber swaps `i_+,i_-`
  and fixes `D_(6,*)`, so `rank(C_ext-I)=1` and
  `dim H^1(C4;V_ext,C_ext)=3-1=2`.  The two invariant directions are the
  relative line `i_++i_-` and the pure endpoint line `D_(6,*)`.
- **Repair / strongest survivor:** no **two-dimensional coefficient
  representation** carrying the descended chamber action has two-dimensional
  twisted H1; its invariant rank is one.  The unique minimal
  endpoint-retaining ambient has coefficient dimension three and formal H1
  dimension two.  One class is endpoint-only, the all-child Q10 core pairs
  only with the relative summand, and the same-copy closure edge is still
  absent, so neither formal group is a physical word-current or D5 flux.
- **Reusable rule:** when a cohomology dimension is used as a no-go, state
  both the base complex and coefficient representation.  Quotienting a
  sidecar can lower the invariant rank; restore the minimal extension and
  recompute before quantifying over all carriers.

## MISTAKE-419 (2026-08-16, THM-3532 promotion gate) -- target postcomposition was called conjugacy, and optimized replay erased the certificate

- **What failed:** historical W1/W2 prose called
  `W_i=T_i o F` a tame conjugate and treated the pulled coordinate-core
  discriminants as evidence that the fixed raw self-iterate tower should
  transport without moving its packet chart.  A target postcomposition has no
  compensating source factor `T_i^(-1)` and is not the self-map conjugate
  `T_i o F o T_i^(-1)`.  Separately, the first THM-3532 companion put every
  truth-bearing gate in a Python `assert`; `python -O` erased all of them, so
  matching normal and optimized stdout did not certify optimized execution.
- **Minimal witness / first failed implication:** exact evaluation gives
  `W_1^2(0,0,-1)=(1,-3,0)` but `T_1F^2(0,0,-1)=(0,0,-2)`, and
  `W_2^2(-1,0,0)=(62,4,0)` but `T_2F^2(-1,0,0)=(-2,0,0)`.
  The transformed boundaries also leave the standard packet chart: their
  five-extremum vectors are `(6,-1,-5,2,-8)` and
  `(1,-8,-16,8,-40)`, versus `(1,-1,-5,2,-8)` for `L`.
- **Repair / strongest survivor:** THM-3532 proves exact one-step covariance
  for every two-sided polynomial equivalence and full all-level covariance for
  every honest conjugate `phi o F o phi^(-1)`, with the five weights moved to
  the `phi^(-1)` coordinate chart.  W1/W2 retain their exact one-step boundary,
  norm, finite-divisor, discriminant, and literal-coordinate-core identities.
  Their honest controls `T_i o F o T_i^(-1)` inherit the complete tower.  The
  repaired companion uses explicit `require`/raise gates; normal and `-O`
  executions are byte-identical to each other and the stored transcript.
- **Reusable rule:** for `G=tau F sigma`, write
  `G^2=tau F(sigma tau)F sigma` before transporting an iteration theorem.
  Left/right equivalence transports one fibre square; conjugacy, or a separately
  proved intertwiner for the inserted `sigma tau`, is required for a tower.
  Never use an optimized replay as a validity gate until executable `assert`
  statements have been excluded from every truth-bearing path.

## MISTAKE-418 (2026-08-16, Fibonacci angle-language sidecar) -- raw Fibonacci slope was mistaken for the canonically normalized primitive Berggren slope

- **What failed:** the provisional Berggren angle-language reflection said
  every consecutive-Fibonacci slope `F_k/F_(k+1)` for `k>=2` gives a
  D-obtuse primitive Euclid triple because the raw ratio is at least `1/2`.
  This skipped the THM-3339 odd/odd normalization.
- **Minimal witness / first failed implication:** at `k=4`, the raw pair is
  `(F_4,F_5)=(3,5)`, but both entries are odd.  The primitive normalization
  is `T(3,5)=((5-3)/2,(5+3)/2)=(1,4)`, whose slope `1/4` lies below the
  U-obtuse wall `alpha=(sqrt(145)-9)/8`.  Normalization is not
  slope-preserving.  More generally it is required exactly when
  `k=1 mod 3`, and the normalized slope `(n-m)/(n+m)<=1/3<alpha`.
- **Repair / strongest survivor:** the canonical Fibonacci three-ray locus
  has exact chamber word `D,D,U` periodically for `k=2,3,...`; it is sparse
  in both obtuse angle languages.  The wall cycles, CDF returns, count
  recurrences, `16/41,9/41,16/41` densities, regular languages, harmonic
  asymptotics, and signed-`C4` theorem are unchanged.  A raw Fibonacci
  projective cycle may still be used only when explicitly typed as the raw
  recurrence quotient rather than the normalized Berggren node.
- **Reusable rule:** before transporting a recurrence on Euclid parameters
  into the primitive Pythagorean/Berggren tree, apply the parity-content
  normalization and recompute every slope-dependent predicate.

## MISTAKE-417 (2026-08-16, U_full owner-node common-base spectrum) -- full Fourier support of a delta-cell lift was mistaken for genuine two-coordinate mixing

- **What failed:** the first checkpoint of the owner-node source/endpoint
  coupling reported `(91,1,6,12,72)` raw support and `(72,0,0,0,72)` after
  output ANOVA as though the actual `U_full` integrand had retained a
  load-bearing seven-cell coordinate.
- **Minimal witness / first failed implication:** `U_full` contains the owner
  factor `||13t||<1/14`.  Under the exact desheeting
  `t=(y+u)/13`, it becomes `||y||<1/14`, which is precisely
  `cell_0`.  All 2,197 character rows therefore have cell support
  `(2197,0,0,0,0,0,0)`, and the inverse table has exact form
  `delta_0(ell)R(t)` and rank one.  Its ANOVA is the separable outer product
  `(delta_0-1/7)(R-mean R)`.  A delta vector has all seven Fourier modes, so
  full mixed *support* does not imply nonseparable cell/residue interaction.
- **Repair / strongest survivor:** the common-base product, literal guard
  controls, pointwise-zero same-root hostile, thirteen-class profile, and
  nonzero `(1,0,6)` role bridge survive.  The seven-cell spectral-closure
  interpretation is withdrawn.  A future candidate must show at least two
  occupied cells and matrix rank at least two, or use a refiner not implied by
  the endpoint owner factor, before a mixed-mode census is evidentiary.
- **Reusable rule:** before reading multidimensional Fourier support, compute
  coordinate support and tensor rank.  Centering a separable boundary delta
  can populate every formal mixed frequency without creating a coupled
  observable.

## MISTAKE-416 (2026-08-16, historical outside Keller family) -- an unstored prose description was promoted as an explicit verified family

- **What failed:** the Keller file
  `THM-1605-infinite-family-extent-vs-mechanism.md` labelled the historically
  reported maps `E_m` verified and asserted that `m=2` was the fixed sporadic
  map, even though the repository retained no coordinate definition, source
  artifact, or higher-member replay.  Its companion checked only the fixed
  map and the arithmetic identity `2m-1=1+2(m-1)`.
- **Minimal witness / first failed implication:** targeted searches of the
  current canon, computations, results, reflections, and retained source
  references find no literal `E_m` formula.  The nearest lawful explicit odd
  family is THM-3448's cyclic weighted subfamily, reindexed by `ell=2m-3` in
  THM-3517.  Its first-coordinate degrees are exactly `7,17,27,...`, which
  disagree with the old THM-1605 prose string `7,13,26,43,64,...`.  Therefore
  the numerical `m=2` overlap cannot identify the higher constructions.
- **Repair / strongest survivor:** THM-1605 is now a partially verified
  historical record.  The fixed map's degree, determinant, three-point fibre,
  and THM-1350's conditional `1+2k` involution-orbit mechanism survive.
  THM-3517 supplies a disjoint explicit `m=3` test with full formula,
  `S5` monodromy, three coordinate quintics, and exact Jelonek components;
  it makes no claim to recover the unstored family.
- **Reusable rule:** a numerical degree/fibre synopsis is not a mathematical
  object.  Before extending a named family, require a literal formula or a
  citable source, instantiate the first new parameter, and compare an
  invariant finer than the overlapping base case.

## MISTAKE-415 (2026-08-16, fixed Keller norm-face extrapolation) -- the exposed pole pair is not a closed scalar recurrence

- **What failed:** the pinned values `e=1,7,43` and exposed factors
  `m=0,3,15` suggested `e_next=6e+1` and `m_next=2e+1`, hence the first
  untested pair `(259,87)` for `G=L^43N(J)`.  This treated the visible
  `(e,m)` data as though the norm read only the exposed face.
- **Minimal witness / first failed implication:** exact extraction from the
  canonical `66,146`-term `J` ledger gives
  `min(i-k)=-43` with face `y^99z^43` and
  `min(i-j-2k)=-185` with face
  `y^99z^28(y^2+27z)^10(y^2+108z)^5`.  The inverse-chart linear and
  quadratic edges therefore give
  `e_next=43-(-43)-(-185)=271` and `m_next=99`, not `(259,87)`.
- **Repair / strongest survivor:** THM-3506 proves the five-face one-step
  transform `(e,m)->(7e-2m,3e-2m)`.  It is derived from the inverse chart and
  resultant geometry, and it gives the exact fixed-map pair `(271,99)` plus
  `v_L(N(G))=-271`.  All-level iteration remains conditional on renewal of
  the `z`-top and `gamma=i-j-5k` faces.
- **Reusable rule:** a Newton face observed on divergent sheets need not be a
  state variable closed under norm.  Before extrapolating, include every
  inverse-branch valuation, test each complete initial form in its residual
  algebra, and separately test the finite sheet.

## MISTAKE-414 (2026-08-16, provisional THM-3497 pole terminology) -- a normalized polar coefficient is not an analytic residue

- **What failed:** the provisional fixed-drift Berggren note called
  `lim_(x->1/3)(1-3x)F(x)=17/96` and
  `lim_(x->-1/3)(1+3x)F(x)=-1/96` the residues at the two dominant poles.
- **Minimal witness / first failed implication:** if
  `F(x)=A/(1-lambda*x)`, then the displayed normalized polar coefficient is
  `A`, whereas the analytic residue at `x=1/lambda` is `-A/lambda`.  Here the
  actual residues are `-17/288` at `1/3` and `-1/288` at `-1/3`, not the two
  numbers printed in the provisional note.
- **Repair / strongest survivor:** THM-3497 names the two defining limits
  normalized polar coefficients.  Their coefficient asymptotic
  `(17/96)3^n-(1/96)(-3)^n+O(3^(n/2))`, the parity densities, and logarithmic
  density `17/96` are unchanged.
- **Reusable rule:** whenever a generating-function pole is summarized by a
  coefficient, print the defining limit.  Convert to an analytic residue only
  after dividing by the derivative of the vanishing denominator factor.

## MISTAKE-413 (2026-08-16, provisional THM-3494 discriminant unit) -- an odd branch divisor does not determine the full field square class up to constants

- **What failed:** after proving `D_n=u_n B_n` with `u_n in k^*`, the
  provisional statement wrote `[D_n]=[B_n]` in `K^*/K^{*2}`.  This silently
  killed the constant square class `[u_n]`; only divisor parity is insensitive
  to a nonzero constant unit.
- **Minimal witness / first failed implication:** the theorem's own cubic row
  has `D_3=-B_3`.  Over `Q`, `-1` is not a square, so
  `[D_3]=[-B_3]!=[B_3]`, even though `B_3` is exactly the unique reduced odd
  branch divisor.  The quartic scalar `-1/16` has the same negative unit
  class; the fifth-degree displayed scalar is a rational square.
- **Repair / strongest survivor:** every primitive coordinate eliminant has
  square class `[D_n]=[u_nB_n]`, because its trace-form index and raw
  normalization contribute squares.  The unique odd divisor carrier remains
  `B_n`.  THM-3494 now distinguishes the full field square class from its
  image after quotienting constant units and was repaired before promotion.
  THM-3498 supplies a second load-bearing check: the complete level-four
  ledger gives `[Delta_4]=[2G]`, not `[G]`; the factor `[2]` survives both the
  cubic norm and old-`L` cancellation.
- **Reusable rule:** a UFD factorization controls valuations, not unit square
  classes.  Whenever a discriminant is normalized "up to scalar," retain that
  scalar before making a claim in `K^*/K^{*2}`; discard it only in an explicitly
  unit-quotiented divisor group.

## MISTAKE-412 (2026-08-16, provisional THM-3494 XOR and hash typing) -- an exact volume coboundary need not be edgewise square-trivial

- **What failed:** the provisional THM-3494 tournament paragraph labelled
  `g_ij=v_j/v_i` and then said that passage to square class, or to the bit
  recording square versus nonsquare, sends every edge to zero.  It also
  labelled the primary output hash as LF-normalized although the pinned value
  was the raw Windows-CRLF replay hash.
- **Minimal witness / first failed implication:** over `Q`, vertex volumes
  `(1,2,6)` give exact edge ratios `(g_12,g_23,g_13)=(2,3,6)`.  All three are
  nonsquares, so their square classes are nonzero and the naive bits
  `(1,1,1)` violate `b_12 XOR b_23=b_13`, even though
  `g_12 g_23=g_13`.  The square/nonsquare indicator is not a character on a
  square-class group of rank greater than one.  Separately, the raw CRLF hash
  `723798b4...83f4` normalizes to the stored LF hash
  `6402ddf6...fb1f`.
- **Repair / strongest survivor:** the unsquared classes obey the exact
  coboundary law `[g_ij]=[v_j]-[v_i]` in `K^*/K^{*2}` and may be nonzero
  edgewise.  Every chosen character `K^*/K^{*2}->F_2` gives a lawful XOR
  coboundary, while the discriminant ratios `g_ij^2` are individually
  square-trivial.  Thus the common discriminant square class, zero graph
  `H^1` class, and unsquared index sidecar all survive.  THM-3494 was repaired
  before promotion, and its frontmatter now pins LF-normalized bytes.
- **Reusable rule:** distinguish an exact multiplicative cochain, its full
  square-class image, a selected binary character, and the non-homomorphic
  predicate "is nonsquare."  Normalize line endings before naming a hash
  basis rather than inferring the basis from replay equality.

## MISTAKE-411 (2026-08-15, provisional THM-3481 Walsh scope) -- odd support gives full Walsh spectrum only beyond the one-variable cube

- **What failed:** the first provisional summary applied THM-3481's full
  Walsh-support conclusion to every odd innovative terminal profile.
- **Minimal witness / first failed implication:** at depth `k=3`, the phase
  cycle has size `P_k=2`, so the innovation cube has one variable.  Its exact
  terminal row has weight one and maximal ANF degree one, but one of its two
  Walsh coefficients is zero.  Odd truth-table support makes all Walsh
  coefficients `2 mod 4` only when the cube size is divisible by four.
- **Repair / strongest survivor:** every odd innovative terminal profile has
  maximal ANF degree and uses every innovation variable.  Full nonzero Walsh
  support and the `2 mod 4` law hold when `log_2(P_k)>=2`; `k=3` is the sharp
  boundary.  The YAML summary, theorem statement, and hostile scope were
  repaired before promotion.
- **Reusable rule:** when a parity argument subtracts twice an odd signed sum
  from the ambient cube size, audit the smallest cube separately before
  turning a stable congruence into a universal spectral claim.

## MISTAKE-409 (2026-08-15, 7x13 private-support sidecar) -- zero absolute H1 flux was mistaken for weighted spectral singularity

- **What failed:** the first private-support reflection correctly proved that
  every owner-potential edge current is a coboundary, but its closing sentence
  then said spectral closure required a non-coboundary weighting.  The graph
  cohomology quotient and the signed matrix-tree determinant are different
  target predicates; the latter does not descend to `H^1`.
- **Minimal witness / first failed implication:** at `k=1`, the canonical
  private-count potential is `(4,2,4,4,8,4,4,4)`.  Its edge gradient pairs
  zero with all six cycles, while the two tetrahedral tree sums and bridge
  weight are `-64`, `-136`, and `-4`, giving determinant `-34816`.  Thus
  “coboundary, therefore spectrally singular” is false.
- **Repair / strongest survivor:** THM-3482 derives all three
  residue-class determinant factorizations and proves they are nonzero for
  every `k>=1`.  Holonomy remains necessary for a nonzero **absolute H1
  class**, not for this owner-order weighted determinant.  The original
  `7x13` incidence rank, six-dimensional cycle space, and coboundary no-go all
  remain valid.
- **Audit correction:** the first THM-3482 draft wrote
  `<c,w_k>=<Bc,f_k>` with `Bc in Q^7` but the displayed `f_k in Q^8`.  The
  typed identity uses the hub-gauged representative
  `fbar_k=(f_(k,i)-f_(k,5))_(i!=5)` and reads
  `<c,w_k>=<Bc,fbar_k>=0`; equivalently one may pair `B_full c` with `f_k`.
  This repairs notation only, not the cycle-vanishing or determinant claims.
- **Reusable rule:** after quotienting out exact objects, do not assume a
  different nonlinear observable factors through that quotient.  Write the
  actual target map and hostile-test it on a canonical exact element before
  promoting a quotient obstruction to a spectral obstruction.

## MISTAKE-408 (2026-08-15, concurrent theorem reservations) -- a stale namespace check allowed duplicate theorem IDs

- **What failed:** commit `958234d4b` first published `THM-3479` as the
  relation-current two-transplant reservation.  The later Rule-30 commit
  `a0a5fe066`, prepared against a stale view of `origin/main`, independently
  created another `THM-3479`.  Both were honest empty stubs, but the second
  publication made the YAML theorem identifier nonunique.
- **Minimal witness / first failed implication:** the next fetched tree
  contained both
  `THM-3479-literal-half-twist-relation-current-two-transplant-certificate.md`
  and `THM-3479-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum.md`.
  A negative scan of one stale checkout was therefore not a lock across the
  concurrent fetch/push window.
- **Repair / strongest survivor:** the earlier relation-current reservation
  keeps `THM-3479`; the Rule 30 cyclic-arc candidate is renumbered to
  `THM-3481` before any theorem promotion or proved dependency is added.  Its
  companion `THM-3480` reservation is unchanged; no mathematical claim is
  affected.
- **Second exact recurrence:** commit `01f394823` published the private-count
  gradient `THM-3482` reservation first; the later factorial commit
  `b79304f74` independently published a second empty `THM-3482`.  The former
  keeps its ID and the still-empty factorial stub moves to the globally
  scanned `THM-3483`.  This recurrence occurred even after a full remote-tip
  pre-scan: a scan is an observation, not a lock, and uniqueness must be
  repeated after every expansion rebase as well as every reservation rebase.
- **Third exact recurrence:** commit `bb86eb0a2` published the Rule-30 dyadic
  wrap reservation as `THM-3493`; `7aa9908e72` then published a weighted-lift
  `THM-3493`, and `2e2dc346bd` published a third empty `THM-3493` after its
  pre-push command printed the collision but failed to gate the following
  rebase and push.  The first keeps `3493`; the weighted atlas moves to `3494`
  and the fixed-map Keller stub to `3495` before proof attachment.  A printed
  diagnostic is not a safety gate: the command must exit on a collision hit.
- **Reusable rule:** immediately re-fetch and repeat both YAML-ID and filename
  searches after any reservation or expansion rebase and before pushing.
  If two reservations still race, publication order keeps the earlier ID and
  the later empty stub moves.

## MISTAKE-407 (2026-08-15, THM-3472 layer transport) -- owner doubling was falsely called augmented-primitive and even nonbijectivity was mistaken for a cover boundary

- **What failed:** the provisional THM-3472 proof observed
  `gcd(Q,2S)=gcd(Q,S)` for odd `Q` and concluded that doubling a primitive
  fixed-zero owner family produced a primitive half-twist family.  THM-3415's
  half-layer gate is instead the augmented condition
  `gcd(2Q,r_1,...,r_k)=1`; an all-even doubled family fails that gate.  The
  proof also treated failure of `ell -> 2ell+1` to be bijective at even `Q`
  as a boundary for cover transport rather than only for conjugacy.
- **Minimal witness / first failed implication:** at `Q=15`, the primitive
  fixed-zero cover `S=(1,2,3,4,5,7)` has `gcd(15,S)=1`, but its valid doubled
  literal half cover `(2,4,6,8,10,14)` has
  `gcd(30,2S)=2`.  Thus “preserves primitivity” was the first false
  implication.  At `Q=8`, owner `s=4` occupies only even fixed-zero sheets
  and doubles to the empty self-opposite half owner; this destroys a mask
  conjugacy but it cannot help cover the odd image sheets.
- **Repair / strongest survivor:** the exact one-way identity
  `B_Q(2s)(ell)=Z_Q(s)(2ell+1 mod Q)` holds for every `Q`.  A fixed-zero
  family covering all sheets therefore covers the entire image of the affine
  map; when `Q` is even, the self-opposite `s=Q/2` is empty on that image and
  may be deleted.  The result is an ordinary transverse literal half cover
  with no more owners, which is all the rank comparison needs.  Deriving the
  active-gcd divisor interface explicitly from THM-3405 (and matching
  THM-3415's divisor minimum) then proves the stronger equality
  `rho_ZMC(q)=rho_H(q)` for every `q>=2`.  Augmented primitivity and even-sheet
  invertibility are not preserved.
- **Reusable rule:** distinguish a cover transport from an isomorphism of
  typed primitive presentations.  Before declaring a parity obstruction,
  identify the actual image, check whether discarded masks meet it, and ask
  whether the target theorem needs bijectivity, cardinality equality, or only
  a cover of no larger size.

## MISTAKE-406 (2026-08-15, THM-3469 rank typing) -- literal half-twist rank was silently used as full zero-cochain rank

- **What failed:** the provisional THM-3469 proof used THM-3455's cap-seven
  literal half-twist rank `rho_H` as a lower bound for the full zero-mode-
  cochain rank `rho_ZMC`, which also ranges over the fixed-zero layer and
  proper-divisor ancestry.  The numerical rank word was correct, but the type
  conversion needed a proof.
- **First unsupported implication:** “THM-3455's cap-seven atom sieve makes
  the eight-owner upper bound exact” does not follow until a hypothetical
  fixed-zero cover on every odd divisor is transported into the classified
  half layer.
- **Repair / strongest survivor:** for odd `Q`, the sheet map
  `phi(ell)=2ell+1 mod Q` is bijective and
  `B_Q(2s)(ell)=Z_Q(s)(phi(ell))`, because
  `dist_(2Q)(2a,0)=2 dist_Q(a,0)`.  After canonical sign representatives are
  chosen, every fixed-zero `k`-cover becomes a half-twist `k`-cover.  Every
  divisor of the odd family `42k-3` is odd, so THM-3405 divisor descent plus
  THM-3455's atom sieve gives the required full `rho_ZMC` lower bound.  The
  exact rank-4/6/7/8 law survives.
- **Reusable rule:** never identify a layer-specific rank with a global rank
  merely because all displayed witnesses use that layer.  Give an explicit
  conjugacy or audit every omitted layer and every divisor.  The oddness of
  `Q` makes this transport a literal sheet conjugacy.  MISTAKE-407 later
  strengthens the rank theorem: multiplication by two is not a sheet
  permutation at even modulus, but one-way cover transport still survives.

## MISTAKE-405 (2026-08-15, THM-3464 q=123 ancestry) -- a primitive mixed-order witness was mistaken for noninheritance

- **What failed:** the provisional THM-3464 package correctly found a
  primitive full-modulus rank-eight witness at `q=123`, with quotient-order
  profile `3^1,41^3,123^4`, but then said that rank eight was not inherited
  from a proper divisor.  Primitivity of one realization does not exclude a
  second realization with nontrivial active gcd.
- **Minimal witness / first failed implication:** at `Q=41`,
  `(3,5,11,19,28,33,37,39)` is an exact primitive half-twist eight-cover,
  while the independently audited cap-seven exclusion proves it is minimal.
  Its threefold pullback
  `(9,15,33,57,84,99,111,117)` covers at `q=123` with active gcd three.
  The first false implication was “this displayed q=123 witness is primitive
  and mixes divisor orders, therefore the rank cannot also descend from
  q=41.”
- **Repair / strongest survivor:** `rho_ZMC(41)=rho_ZMC(123)=8`.  The q=123
  grade has **coexisting** inherited and primitive realizations: the pullback
  from q=41 and the distinct mixed-order active-gcd-one packet
  `(1,40,42,81,82,83,117,122)`.  All q123/q227 lower bounds, witnesses,
  centres, widths, multiplicities, the first-seven U-spine rank word, and the
  square-root hostile survive unchanged.
- **Reusable rule:** a primitive representative proves existence in the full
  divisor layer, not uniqueness of ancestry.  Before declaring a grade
  noninherited, search every proper divisor at the same rank and pull each
  positive witness back.  Record “primitive versus inherited” as realization
  types that may coexist, not mutually exclusive rank labels.

## MISTAKE-404 (2026-08-15, THM-3453/3455 half-twist scope) -- a fixed common mode centre was discarded with the arbitrary-centre problem

- **What failed:** THM-3453 correctly classified the literal half-twist masks
  `B_(q,r)`, but its status, statement, companion transcript, and closing scope
  said that the witnesses provided no common physical time.  THM-3455
  inherited that wording and said its positive spine labels realized neither a
  common time nor a zero-mode cochain.  Those sentences conflated failure to
  classify arbitrary centres with failure of the centre already built into
  the literal half twist.
- **Minimal witness / first failed implication:** directly from the two
  definitions,
  `B_(q,r)=D_(q,r)(1/(2q))`.  At this source centre the containing THM-3398
  mode has `h=r mod 2q`; taking `n=0` in its centre lattice gives
  `x=(n+h/(2q))/r=1/(2q)`.  Thus every labelled literal cover has one common
  physical centre and identically zero complete mode cochain.  The `q=11`
  witness `(1,2,3,5,7,9)` is already an exact six-owner partition at that
  centre.  The first failed implication was “the bare rank word forgets the
  owner/mode sidecar, therefore no common time exists.”
- **Repair / strongest survivor:** THM-3453 now states the fixed-half-twist
  realization explicitly and distinguishes it from the still-unclassified
  arbitrary-centre problem.  THM-3455 says that each positive rank grade has
  such a realization after a witness is reattached, while its compressed
  periodic word still forgets the owners and widths.  The projective wedge at
  this centre collapses: `A_i=r_i` and `P_ij=0`.  Hence no nonzero current,
  endpoint relation, decrement, spectral closure, or LRC(14) conclusion is
  gained.  All cap-seven support, rank, density, and Fibonacci calculations
  remain unchanged.
- **Reusable rule:** when a quotient object is defined by evaluation at a
  named phase, separate “the phase is not retained by the compressed
  invariant” from “no realizing phase exists.”  Before denying a physical
  lift, substitute the defining phase and then check the selected-mode centre
  lattice; only after that ask which current or endpoint data were lost.

## MISTAKE-402 (2026-08-15, THM-3454 spine-index/depth and replay typing) -- a shifted branch index was called rooted depth and raw newline bytes were overclaimed

- **What failed:** THM-3454 repeatedly called the parameter `t` in
  `P_t=U^(t-1)(3,4,5)` a rooted depth, and it said normal, optimized, and
  stored transcripts agreed byte for byte without naming newline
  normalization.
- **Minimal witnesses / first failed implications:** `P_1=(3,4,5)` is the
  root at depth zero, so `t` is the spine index and its rooted depth is `t-1`.
  Consequently the marked index seam `c_12=x_0` becomes `c_12=h_0+1` for
  rooted depths `h_i=x_i-1`, and the homogeneous recurrence becomes
  `h_2=h_0+h_1+1`.  Independently, Windows stdout uses CRLF while the stored
  artifact uses LF, so their raw bytes need not agree.
- **Repair / strongest survivor:** all Fibonacci selections are now typed as
  spine indices `F_n`, at rooted depths `F_n-1`.  The theorem records both
  the homogeneous index recurrence and its affine rooted-depth form.  Every
  metric identity, six-cost level, Pell classification, recurrence mode, and
  semantic output is unchanged because differences erase the shift.  Replay
  equality is asserted only after LF normalization.
- **Reusable rule:** whenever a ray is parameterized by
  `P_t=M^(t-t_0)P_(t_0)`, declare the index/depth shift before transporting an
  absolute recurrence.  Difference invariants survive affine shifts;
  homogeneous recurrences generally do not.  State the exact normalization
  basis for every byte-equality claim.
## MISTAKE-403 (2026-08-15, THM-3456 Rule 30 source scope) -- a citation was scoped below its actual use and an active prize listing was paraphrased as a literal openness statement

- **What failed:** the provisional THM-3456 frontmatter said Wolfram's 2019
  announcement was used only for the prize statements and local rule, while
  the theorem also cited its sideways two-column discussion.  The companion
  Rule 30 reflection separately cited the announcement's reported finite-size
  observations.  The theorem then said the official page itself stated that
  the problems remained open, although its literal evidence is an active prize
  listing and an ongoing submission rule.
- **Minimal witness / first failed implication:** THM-3456 Sections 2 and 7
  attribute the sideways mechanism to the announcement, and the companion
  reflection's finite-size paragraph makes the separate historical use.  On
  2026-08-15 the official page listed all three prizes
  and accepted submissions until a satisfactory solution is achieved, but
  did not contain the literal sentence “these problems are open.”
- **Repair / strongest survivor:** the two sources now have separate entries.
  The 2019 announcement is cited for the questions, Boolean rule, sideways
  discussion, and historical finite observations; the current page is cited
  for dated prize status only.  The repository treats the questions as open
  on that evidence and says explicitly that this is an inference.  No
  mathematical claim, exact artifact, or no-go implication changed.
- **Reusable rule:** a citation's declared import scope must cover every
  mechanism later attributed to it.  When a current web page supplies status
  by an active listing rather than a declarative sentence, record the check
  date and distinguish the page's literal content from the repository's
  inference.
## MISTAKE-400 (2026-08-15, provisional theorem status routing) -- an unrecognized provisional token exposed an audit-pending theorem as established canon

- **What failed:** the packaged THM-3452 frontmatter began `PROVISIONAL /
  AUDIT REQUIRED + VERIFIED-EXACT COMPANION`.  The bounded-startup router does
  not recognize `PROVISIONAL` as a status token, so it skipped that word,
  selected the later recognized token `VERIFIED`, and placed the unpromoted
  theorem under `Established canon`.
- **Minimal witness / first failed implication:** an exact
  `agents/start_session.py --topic "THM-3452 unequal-depth noncommuting Hensel"`
  replay printed `[VERIFIED]` in the established-canon group even though the
  theorem body and results index both said that independent audit and explicit
  promotion were still required.  Provisional prose did not itself keep the
  file outside the proof graph.
- **Repair / strongest survivor:** the independent immutable-file audit found
  the mathematics, controls, replays, dependencies, hashes, and scope clean,
  so THM-3452 was explicitly promoted and now begins with the recognized token
  `PROVED`.  Its proof body, script, output, and frozen hashes did not change.
  Had promotion not been authorized, the correct audit-pending status would
  have begun `RESERVED / PROVISIONAL ...`.
- **Reusable rule:** every audit-pending theorem status must begin with the
  recognized token `RESERVED`; never place an unrecognized provisional label
  before `VERIFIED-EXACT`.  Before and after promotion, run a topic-specific
  bounded-startup replay and inspect the actual canon/candidate group, not
  only the human-readable status sentence.
## MISTAKE-401 (2026-08-15, refined exact-six transcript portability) -- native path separators defeated byte-for-byte replay

- **What failed:** the refined exact-six SCC sidecar serialized dependency
  paths with `str(Path)`.  Its stored transcript used POSIX `/`, while a fresh
  Windows replay printed `\\`, contradicting the stated cross-platform
  byte-for-byte reproduction even though every graph count agreed.
- **Minimal witness / first failed implication:** on Windows the fresh and
  stored LF hashes were respectively
  `2930d926be17346497d6371d5238726d1130d06c6b44b7cdd8435c5468a269b3`
  and
  `cebb631596ca3cb04d95a39fa7a6edfdd8dea6bf89be96d1f751d5392b812496`;
  replacing only printed `\\` by `/` made the transcripts identical.
- **Repair / strongest survivor:** dependency paths now use `.as_posix()`,
  the stored transcript is refreshed, and the verdict explicitly says
  `k2_body_quotient`.  The exact weights, SCCs, engine hash, and semantic hash
  `d3be3507...bed31` are unchanged.
- **Reusable rule:** deterministic outputs must serialize paths with an
  explicit platform-independent convention.  LF normalization alone does not
  normalize directory separators, and semantic agreement does not justify a
  byte-identical replay claim.
## MISTAKE-399 (2026-08-15, concurrent theorem reservation race) -- filename merges did not protect a shared semantic ID

- **What failed:** THM-3448 was absent when the weighted-Keller boundary
  theorem was first checked and locally reserved, but another clean session
  reserved the same YAML ID for the noncommuting Hensel theorem before the
  next fetch.  Because the filenames differed, rebase had no textual conflict
  and the duplicate IDs were briefly pushed together.
- **Minimal witness / first failed implication:** after commit `8504e472e`,
  both `THM-3448-noncommuting-smooth-hensel-heisenberg-orbit-law.md` and
  `THM-3448-weighted-keller-cyclic-jelonek-inertia-family.md` declared
  `id: THM-3448`.  “Free before local editing” did not imply “free at push.”
- **Repair / strongest survivor:** after both sessions' corrective commits
  were inspected together, the already-pushed Keller reservation keeps
  THM-3448 and the Hensel reservation takes THM-3449.  The mathematical
  content is unaffected, and Git history preserves the collision lineage.
- **Reusable rule:** immediately before pushing a reservation, fetch and
  recheck the intended YAML ID against the fetched remote tree—not merely the
  working tree—and explicitly audit duplicate IDs after every conflict-free
  rebase, since filename-level merging cannot detect semantic ID collisions.
  Once a collision is shared, inspect the other session's corrective commit
  before choosing a successor; two independent “move to the next ID” repairs
  can collide again.

## MISTAKE-398 (2026-08-15, THM-3446 group and alignment typing) -- a direct product was called a free product, and a sharp universal bound was called every packet's exact level

- **What failed:** the provisional THM-3446 truth surfaces called the finite
  abelian exponent group a "free product," imposed `r<=d` even though the
  displayed dependent hostile used two generators on `A^1`, and called
  `M+1` the unqualified exact detection level.
- **Minimal witnesses / first failed implications:** `C_3*C_3` is infinite,
  whereas the equal-depth two-generator fibre uses `C_3 x C_3` of order
  nine.  The `A^1` hostile has `r=2>d=1`.  Conversely, at depths `(1,1,3)`
  two identical shallow translations already give a stabilizer at level
  three, before the universal alignment level `M+1=4`.
- **Repair / strongest survivor:** the group is the direct product
  `prod_i Z/p^(a-c_i)Z`, acting freely under pointwise independence.  The
  ambient statement allows arbitrary `r>=1`; independence itself forces
  `r<=d`.  Every dependence produces a stabilizer by `M+1`, and the
  `(1,2)` delayed hostile shows that bound is sharp, but some relations appear
  earlier.  The mixed-depth carry proof, orbit invoice, transitivity test,
  and exact artifacts are unchanged.
- **Reusable rule:** distinguish algebraic coproducts from freely acting
  groups, and distinguish a sharp universal detection bound from the first
  failure time of each individual packet.

## MISTAKE-397 (2026-08-15, THM-3437 Prüfer-limit semantics) -- the Tate module was called the Prüfer arm

- **What failed:** the promoted synthesis described the inverse limit of the
  selected `Tor_1` packets as recovering the divisible Prüfer channel.  Under
  the chain maps induced by `R_(q+1)->R_q`, however, the transition on kernels
  is multiplication by `lambda`, not the canonical inclusion.
- **Minimal witness / first failed implication:** for
  `Pr_lambda=A[lambda^(-1)]/A`, the inverse system
  `Pr_lambda[lambda^(q+1)] --lambda--> Pr_lambda[lambda^q]` has limit
  `K'[[lambda]]`, a torsion-free complete Tate module.  It cannot be the
  divisible torsion module `Pr_lambda`.  The first failed implication was
  identifying two different transition systems merely because their finite
  stages are the same kernels.
- **Repair / strongest survivor:** the inverse limit recovers the selected
  arm's Tate module and hence its presence bit.  The canonical inclusions
  `Pr_lambda[lambda^q] -> Pr_lambda[lambda^(q+1)]` instead have direct limit
  `Pr_lambda`.  The finite Tor formulas, Euler cancellation, Mittag--Leffler
  statement, and filtration-loss boundary are unchanged; THM-3437 now records
  both limits explicitly.
- **Reusable rule:** a tower is its objects plus its arrows.  Before naming a
  limit, write the induced transition map and compare algebraic type (torsion,
  divisibility, completeness), not just finite-stage ranks.

## MISTAKE-396 (2026-08-15, Keller degree-spectrum inheritance) -- a public all-degree family was ignored by later open verdicts

- **What failed:** THM-1330, THM-2465, HYP-9027, and HYP-9030 continued to say
  that degrees `4` through `7`, global G1, and even-degree Keller maps were
  open or conjecturally impossible.  Yet THM-1300's maintained attribution
  amendment already named the public weighted-lift family and its explicit
  quartic witness.
- **Minimal witness / first failed implication:** the polynomial 2-jet map
  `G` in THM-3438 has `det JG=-6`, generic degree four, and
  `G(1,0,0)=G(-1,0,2)`.  Thus `4 in KDeg(3)` immediately refutes
  `KDeg(3)={3^k}`, the claimed global-open G1 verdict, and the proposed
  all-shape odd-degree lens.  The failure was inheritance/search discipline,
  not a subtle implication inside the older local exclusions.
- **Repair / strongest survivor:** the weighted-lift inverse equation proves
  `KDeg(m)={1,3,4,5,...}` for every `m>=3`; this classifies degree values only.
  Its two-root incidence proves `S_n` monodromy and an explicit atom in every
  grade, while compositions populate exactly grades `ab` with `a,b>=3`.
  The monoid law, no-degree-two theorem, ternary `F` subfamily, local G1
  exclusions, and z-affine order-`{1,3}` problem survive.  Classification of
  arbitrary maps in mixed grades still requires monodromy/intermediate fields.
- **Reusable rule:** before declaring a realization gap, search maintained
  theorem amendments as well as filenames and hypotheses.  Separate
  classification of numerical degree values from classification of maps in a
  degree grade.

## MISTAKE-395 (2026-08-15, exact-six truncated minimum) -- a depth-six sentinel was labelled as pool-14 infeasibility

- **What failed:** the first exact-six mutation reflection said that the
  companion computed the unrestricted minimum over all subsets of
  `{1,...,14}` and labelled the `None` histogram bucket “no pool-14 cover.”
  The solver deliberately searches only depths zero through six.
- **Minimal witness / first failed implication:** for
  `F=(1,2,4,6,9,10)` and `D=1260`, the strict target has no cover by at most
  six pool clocks, but `(1,2,3,5,8,9,10)` covers it with seven.  Thus this row
  lies in the old `None` bucket although its exact pool-14 minimum is seven.
- **Repair / strongest survivor:** read `None` as “no cover by at most six
  (`>6` or uncovered).”  The companion now prints that scope and checks the
  seven-clock hostile exactly.  The depth-one through depth-six counts, every
  exact-six completion, the full mutation relation, both SCCs, and the typed
  `7+6=13` stopping boundary are unchanged.
- **Reusable rule:** a capped-search sentinel records failure within its
  searched budget, not global infeasibility.  Freeze a first-outside-budget
  positive control whenever a finite depth cap is load-bearing.

## MISTAKE-394 (2026-08-15, Fibonacci--Berggren 17-adic torsor scope) -- a parameter norm-square was misnamed and the tied root entered a `T6` support claim

- **What failed:** the local-`T4` gate called `m^2+n^2` the “hypotenuse
  squared,” although it is the squared Euclidean norm of the parameter and
  the hypotenuse itself of its Euclid triple.  Separately, the periodic-support
  paragraph allowed an unspecified Fibonacci base index and then said that a
  labelled `T6` was frozen along every resulting support.
- **Minimal witnesses / first failed implications:** for `(m,n)=(1,2)`, the
  Euclid triple is `(3,4,5)`, so `m^2+n^2=5`, not the squared hypotenuse `25`.
  For the support claim, take `n_0=2`, `a=1`, and `J={0}`.  The support
  contains the root window `W_2=(1,1,2,3)`, whose edge products tie at `2`
  and `3`; the reflection itself correctly declares that no `T6` exists
  there.
- **Repair / strongest survivor:** call the `T4` observable the parameter
  norm-square, equivalently the Euclid-triple hypotenuse.  For the pure
  periodic realization choose `n_0` in `{3,...,11}`; these nine indices
  represent every line of the square base cycle, lie in the first period,
  and are all in the tie-free `T6` domain.  Changing to any other compatible
  base lift is exactly the already-recorded translation gauge.  The general
  Hensel lemma, both Legendre cycles, all finite rows, affine hostile, Boolean
  identities, harmonic density, tournament loss, hashes, and no-LRC boundary
  are unchanged.
- **Reusable rule:** distinguish a parameter's norm-square from the square
  of the derived object's norm, and intersect a periodic carrier with the
  declared tie-free domain before claiming that its tournament state is
  constant.

## MISTAKE-393 (2026-08-15, THM-3435 covering degree/sign scope) -- a componentwise bijection was assigned to the grid union

- **What failed:** the proof said that the degree-`d` circle map sends all `d`
  small arcs *bijectively* to the widened arc.  It also quantified the partner
  law modulo `4R` while excluding only `r=R` from its empty sign class.
- **Minimal witnesses / first failed implications:** already at `d=2`, both
  inverse-image arcs map bijectively onto the target, so the map on their
  union is two-to-one rather than one-to-one.  For `Q=2R`, the residue
  `r=3R` is sign-equivalent to the excluded residue `R`; it is empty and
  satisfies `2R-r=r modulo 4R`, so it is not a complementary pair of two
  owners.  The fibre criterion and partner identities themselves remain
  true.
- **Repair / strongest survivor:** each inverse-image component maps
  bijectively and the whole union maps with degree `d`; a fixed grid orbit
  therefore has a unique selected point exactly under the stated widened-arc
  condition.  In the two-sheet paragraph, choose the canonical odd sign
  representative `0<r<2R` and exclude `r=R`.  Equations `(4)--(11)`, every
  endpoint/carry/Boolean consequence, both companions, and the census survive.
- **Reusable rule:** record the degree when passing from several inverse-image
  components to one target arc, and state whether coefficient quantifiers are
  on residues or on sign orbits before calling a partner pair distinct.

## MISTAKE-392 (2026-08-15, refined dyadic auxiliary capacity) -- a full-order half-block bound was applied to arbitrary transverse pullbacks

- **What failed:** the first exact-six mutation reflection called
  `ceil(3822/7)=546` the one-half-clock capacity on the two residual k=2
  dyadic rows.  That scalar bounds a full-order interval block; it does not
  bound a lower-order mask pulled back to `3822` sheets.
- **Minimal witness / first failed implication:** at `Q=3822`, residue
  `r=2548` has quotient order `Q/gcd(Q,r)=3`.  Its strict half-twist block is
  the pullback of the nonempty order-three block and has `1274` sheets, so
  `1274>546`.  The first failed implication was “the ambient modulus controls
  every block size” instead of descending to the owner's quotient order.
- **Repair / strongest survivor:** exact enumeration of every transverse
  residue modulo `2Q` gives maximum block size `1274`, attained by the stated
  order-three residue.  The two unsupported targets have sizes `1530` and
  `1560`, so neither can be covered by one auxiliary block even under the
  repaired sharp capacity.  The refined six-pool counts, zero dyadic-hit
  verdict, and absence of any LRC(14) terminal are unchanged.  The exact
  repair is reproduced by
  `lrc14_refined_six_pool_dyadic_stopping_audit_20260815.py`.
- **Reusable rule:** for a half-twist residue of quotient order `m`, charge
  the order-`m` block and its `Q/m` pullback multiplicity.  Never substitute
  `ceil(Q/7)` for an arbitrary transverse capacity without proving full
  order.

## MISTAKE-391 (2026-08-15, zero-cochain rank artifact) -- the indexed output hash did not name the committed transcript

- **What failed:** the results index recorded LF-normalized output SHA-256
  `582b901b...4085e6` for
  `lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.out`, but that hash
  does not equal the committed file.
- **Minimal witness / first failed implication:** direct hashing of the
  committed `6,745` LF bytes gives
  `52684d84bba6076c760285937e29cfd4a81c998324d6a9019e2919d4f764ab5d`.
  Fresh normal and optimized replays are byte-identical after LF
  normalization and give that same hash.
- **Repair / strongest survivor:** replace only the stale index hash with
  `52684d84...764ab5d`.  The source hash, semantic digest
  `233c092a...30c7e`, theorem statements, witnesses, and every audit result
  are unchanged.
- **Reusable rule:** never copy an output hash from working notes.  Hash the
  committed path and an independent fresh replay in both normal and
  optimized modes before publishing the artifact tuple.

## MISTAKE-390 (2026-08-15, zero-cochain divisor ancestry) -- two lift lemmas were stated beyond their used cover scope

- **What failed:** the first divisor-ancestry reflection said without a
  cardinality qualifier that residues `r_1,...,r_s mod M` have positive
  gcd-one lifts iff `gcd(M,r_1,...,r_s)=1`, and called
  `gcd(V)=1` equivalent to `lcm_i(Q/gcd(Q,v_i))=Q` for fixed literal owners.
- **Minimal witnesses / first failed implications:** for the singleton
  `(M,r)=(5,2)`, the augmented gcd is one but every positive lift is
  `2 mod 5` and no one-element family has gcd one.  For fixed
  `Q=5,V=(2,4)`, both quotient orders are five and their lcm is five, while
  `gcd(V)=2`.  Thus both reverse implications fail as literally written.
- **Repair / strongest survivor:** every transverse strict full cover uses at
  least two nonfull blocks.  In that `s>=2` scope, the CRT avoidance proof
  gives gcd-one positive lifts exactly when the augmented gcd is one.
  Primitive literal gcd one always **implies** the quotient-order lcm is `Q`;
  that is the only direction used by the rank floor and rank-four
  classification.  The finite Boolean gate, q15--28 ranks, universal floor
  four, and `rho_ZMC(q)=4 iff 8|q or 9|q` are unchanged.  An independent
  rare-coordinate branch-and-bound census through `Q=500` found no rank at
  most three and only the primitive half-twist positives `Q=8,9` at rank
  four.
- **Reusable rule:** distinguish literal integer gcd, gcd modulo the finite
  owner modulus, and existence of a gcd-one lift.  State the number of
  selected types before promoting a modular gcd condition to an iff.

## MISTAKE-389 (2026-08-15, all-owner divisor-chart probe) -- a synchronized half-grid physical time was mistaken for a common mode centre

- **What failed:** the first all-owner divisor-chart draft correctly derived
  the necessary condition `2quc in Z` for every owner at a common THM-3398
  mode centre, but then treated that condition as sufficient.  Its exact
  `direct_mask` computation evaluates danger sets at the physical time `c`;
  it does not prove that `c` is the centre of each selected consecutive mode.
  Consequently the draft falsely labelled its ranks as zero-mode-cochain
  ranks and compared them with the capped mobile common-centre atlas.
- **Minimal witness / first failed implication:** at `q=15,c=1/150`, owners
  `(5,40,50)` have exact danger sets equal to the three residue classes modulo
  three, so they partition all sheets and satisfy `2quc in Z`.  Their active
  gcd is `d=5`; writing `a=2qdc` and `g=gcd(q,d)` gives `(a,g)=(1,5)`.
  THM-3405 proves that a zero mode cochain requires `g|a`, so this physical
  half-grid partition is not a common-mode-centre certificate.  Equivalently,
  owner five has half-grid residue `h=1`, while the THM-3398 mode formula
  requires `gcd(15,5)=5` to divide `h`.  The first failed implication was
  “`2quc` integral implies `c` belongs to the mode-centre lattice.”  The
  unique containing mode centres are `(0,1/120,1/150)` and their exact
  THM-3398 pair cochain is `(-50,-50,100)`, with norms `(L1,Linf)=(200,100)`;
  the failure is positive drift, not an endpoint ambiguity.
- **Repair / strongest survivor:** THM-3405 supplies the missing mode
  divisibility and proves the genuine zero-cochain gauge has at most the two
  classes `a=0,g mod 2g`.  The divisor-chart artifact is renamed and typed as
  a **synchronized half-grid physical-time** theorem.  Its affine
  normalization, exact q15--28 ranks
  `(3,2,8,2,9,2,3,2,6,2,5,2,3,2)`, capacity bounds, and infinite exact
  half-grid families of ranks two/three/five survive.  The claimed
  zero-cochain ranks, comparison with the owner-14 mobile ranks, finite
  “saving” support, and corresponding harmonic weights are withdrawn.  No
  LRC row or ledger count ever followed.  A subsequent exact reconstruction
  strengthens the surviving separation: on every even `q>=8`, the same
  even/odd rank-two partition at `c=a/q^2`, for odd `1<=a<q/7`, has unique
  containing-mode cochain `P=-a q^2/2` and THM-3405 scalar/gauge
  `(a,q/2)`.  Thus the quotient can hide a whole quadratic positive-drift
  ladder, not merely the isolated q15 hostile.  Applying the corrected gauge
  and retaining the primitive owner-gcd sidecar gives the genuine
  unrestricted positive-transverse zero-cochain ranks
  `(6,4,8,4,9,5,8,6,6,4,6,7,4,7)` for q15--28.  Their proof is the divisor
  minimum over primitive fixed/half-twist covers, not a reinterpretation of
  the physical chart minima.
- **Reusable rule:** an interval containing a physical time is not centred at
  that time.  After deriving a half-grid integrality condition, separately
  verify the mode-residue divisibility (or the full centre-lattice formula)
  before setting the affine mode cochain to zero.

## MISTAKE-384 (2026-08-15, THM-3401 scope prose) -- fixed source centre zero was identified with the entire zero-cochain locus

- **What failed:** THM-3401's precise statement and proof correctly compute
  the cover rank at physical source time `t=0`, but its status and final scope
  sentence called this “equivalently the zero-cochain slice.”  In THM-3398,
  `p_ij=0` says that the selected centre lifts are equal to one common
  rational `c`; it does not force the surviving additive gauge `c` to vanish.
- **Minimal witness / first failed implication:** at `q=16`, owners
  `(2,6,10,14)` have selected blocks
  `(0,7,8,15)`, `(2,5,10,13)`, `(1,6,9,14)`, and `(3,4,11,12)` at common
  centre `c=1/32`.  They partition all sheets and every affine cochain value
  is zero, so the mobile zero-cochain rank is four.  THM-3401 proves the
  fixed-zero rank is five.  Translating `c` to zero would shift a common
  sheet label by `qc=1/2`, which is not a permutation of `Z/16Z`.
- **Repair / strongest survivor:** THM-3401 is now consistently labelled a
  **fixed-source-centre-zero** theorem, a proper sub-slice of the zero-cochain
  locus.  Its statement, classification, ranks, boundary theorems, script,
  output, and semantic digest are unchanged.  The independent mobile
  common-centre atlas on owners `1,...,14` gives ranks
  `(6,4,8,4,9,6,8,6,6,6,7,8,9,8)` for `q=15,...,28`, strictly below the
  fixed-zero ranks exactly at `{16,18,22,23,25,27}`.
- **Reusable rule:** a vanishing difference cochain kills relative
  coordinates, not a common additive gauge.  Before normalizing that gauge,
  verify that the induced translation acts on the retained labelled fibre;
  here it is a common cyclic sheet relabelling only when `qc` is integral.

## MISTAKE-388 (2026-08-15, THM-451 skew-tower Smith status) -- a theorem for the ambient matrix class was left labelled as a tower conjecture

- **What failed:** THM-451 correctly computed the flat Smith form of the
  skew-doubling tower at orders 16 and 32, but labelled its continuation
  conjectural beyond order 32.  The classification step searched the tower
  internally and missed a theorem for every skew-Hadamard matrix.
- **First failed implication:** THM-447 already proves that every tower level
  is skew-Hadamard.  Michael--Wallis (1998), reproved as Theorem 2.4 of
  Hacioglu--Keman (2014), proves that every skew-Hadamard matrix of order
  `4m` has Smith form
  `(1,2^(2m-1),(2m)^(2m-1),4m)`.  Thus the displayed THM-451 formula was
  true, but its `CONJECTURED` status was false immediately after its own
  ambient-class identification.
- **Repair / strongest survivor:** THM-451 now marks the Smith law
  **CITED/PROVED** at every tower order.  Its binary rank law is likewise
  all-level through THM-482's `d^+` code theorem.  The finite Hadamard
  equivalence, chirality, Hall-class, density, and transform computations are
  unchanged.  The HYP-2361 half-life question must seek a different invariant;
  neither Smith form nor binary rank can be the order-128 failure.
- **Reusable rule:** after recognizing an object as a standard ambient class,
  search the class theorem before promoting a few computed instances to a
  conjecture.  A correct formula can still carry a false epistemic status.

## MISTAKE-383 (2026-08-15, q=8--15 finite-mode probes) -- a rank-bounded edge list was used to print a global-looking profile

- **What failed:** the q=8 and q=8--14 probes enumerated minimal physical
  covers only through rank five, and the q=15 probe only through rank six,
  but their output fields said `minimal_edges` and `independence_profile`
  without the cutoff.  The latter profiles were those generated by the
  truncated edge lists, not the true full physical independence profiles.
- **Minimal witness / first failed implication:** at q=8,
  `(1,3,5,11,13,14)` covers all eight sheets at source `17/504` and every
  five-subset is noncovering.  Thus q=8 has six additional rank-six minimal
  edges and true `I_6=1217`, not the truncated profile's `1223`.  At q=11
  there are 23 rank-six edges despite the rank-at-most-five list being empty;
  q=15 also has 16 rank-seven and six rank-eight edges beyond its 157
  rank-six edges.
- **Repair / strongest survivor:** the original outputs now explicitly say
  `minimal_edges_through_rank5` / `rank6` and
  `profile_generated_by_rank_le5_edges` / `rank_le6_edges`.  An independent
  exact boundary-and-open-cell bitset audit of every literal subset computes
  the full q=8--15 clutters.  All intended low-rank counts, every `I_5`
  value, all no-rescue statements, the 155 q=15 triphase-required rank-six
  edges, and the finite-mode/cochain formulas survive unchanged.
- **Reusable rule:** an independence profile is global only when every rank
  has been searched, or when a proved rank bound excludes later minimal
  edges.  Otherwise name both the edge cutoff and the fact that the displayed
  profile is generated by that truncated clutter.

## MISTAKE-382 (2026-08-14, THM-3387 canonical grid iff) -- openness does not remove endpoint-only sheet covers

- **What failed:** THM-3387 correctly proved the pointwise image identity
  and the universal aligned-grid identity, but then claimed that the core
  clocks cover the unsupported open cells iff every transverse full cover is
  pointwise core-dangerous.  The inference used only that the two danger
  unions are open and silently discarded the removed grid.
- **Minimal witness / first failed implication:** for
  `F=(4,5,8,9,10,18)`, `q=2`, `C=(2,4,5,9)`, `U=(5,9)`, and
  `(L,D)=(5040,2520)`, the only core-safe transverse full covers are
  `{3/14,11/14}`, both on the removed `D`-grid.  The core therefore covers
  every unsupported open cell although the pointwise safe image omits those
  two base points.  Open sets can differ at a point when another open union
  covers both punctured sides by an endpoint handoff.
- **Repair / strongest survivor:** universally, core-cell completion is
  equivalent to `B_q(U) minus A_C` being contained in the `D`-grid, whereas
  pointwise image equality is equivalent to `B_q(U) subset A_C`.  The literal
  `F subset {1,...,14}` atlas has no endpoint-only exceptions: its core clocks
  are at most seven, and an opposite strict-boundary handoff by distinct
  clocks `c,d` would force `14 | c+d`, impossible when `c+d<=13`.  All
  `15,393` atlas rows, the twelve core rescues, every sector/occurrence count,
  the S172 key identification, the q=2 gcd graph, and the harmonic lattice
  formulas survive unchanged and were independently replayed.
- **Reusable rule:** after deleting a finite boundary grid, distinguish
  equality of open-cell words from equality of the underlying point sets.
  The exact sidecar is the residual set on the deleted grid; openness alone
  does not reconstruct it.

## MISTAKE-381 (2026-08-14, THM-3366 artifact pins) -- escaped text was normalized instead of line endings

- **What failed:** four THM-3366 companion-script hashes were advertised as
  LF-normalized, but the shell recipe replaced the literal four-byte source
  text `\r\n` by `\n` instead of replacing byte pair `13,10` by byte
  `10`.  The resulting values therefore described modified Python source,
  not an end-of-line normalization.
- **Minimal witness / first failed implication:** the refined `k=3` source
  has LF line endings already, so true LF normalization leaves its hash
  `27e4bff5...` unchanged; the advertised `127ef53b...` appears exactly when
  literal backslash-`r`-backslash-`n` inside the source is rewritten.  The
  same mechanism produced the stale primary, `k=1`, and `k=2` pins.
- **Repair / strongest survivor:** the four true LF hashes are now pinned by
  prefixes `372f1b0d`, `65f1e598`, `414b3777`, and `27e4bff5`.  The scripts,
  stored outputs, semantic decisions, and THM-3366 mathematics were
  unchanged; independent all-sector and refined-`k=3` audits replayed the
  exact claims and diagnosed the metadata-only failure.
- **Reusable rule:** normalize actual bytes (`0d 0a -> 0a`) or use a checked
  in-script byte hash.  A quoting-sensitive shell literal is not a valid
  artifact pin unless a hostile file containing both CRLF and literal
  `\r\n` has distinguished the two operations.

## MISTAKE-380 (2026-08-14, THM-101 surplus interpretation) -- a cycle kernel was identified with homology after omitting the next boundary

- **What failed:** THM-101 and three exploratory `beta2_*` companions said
  that `s=dim(Omega_3)-dim(Z_2)` becomes `beta_3` when `beta_2=0`.
  They silently replaced `dim ker(d_3)` by `beta_3` and omitted
  `im(d_4)`.
- **Minimal witness / first failed implication:** at `n=5`, tournament
  `bits=0` has `dim Omega_3=5`, `dim Z_2=4`, `rk(d_4)=1`, and
  `beta_2=beta_3=0`.  Hence `s=1`, not `beta_3=0`; the surplus kernel is
  filled from dimension four.
- **Repair / strongest survivor:** every chain complex satisfies
  `s=rk(d_4)+beta_3-beta_2`.  Under THM-101's verified `beta_2=0`
  conclusion, `s=dim ker(d_3)=rk(d_4)+beta_3`; thus `s=0` still implies
  injectivity and `beta_3=0`, while `s=k` gives only
  `beta_3=k-rk(d_4)`.  The DT+cancellation filling computation and the
  `beta_2=0` conclusion are unchanged.
- **Reusable rule:** a kernel is homology only after quotienting the next
  boundary.  In dimension `p`, retain `rk(d_(p+1))` whenever a chain-rank
  statistic is translated into `beta_p`.

## MISTAKE-379 (2026-08-14, THM-300 quadratic signature) -- the first untested staircase layer contributes two negative directions

- **What failed:** finite evidence through `n=8` was extrapolated to the
  conjecture that the quadratic coefficient matrix has exactly `n-2`
  negative eigenvalues for every `n>=5`.  The same statement also called
  the matrix full-rank at `n=5` despite its displayed nullity one.
- **Minimal witness / first failed implication:** exact one-/two-flip
  Hamiltonian-path DPs and rational symmetric congruence give inertia
  `(20+,8-,0)` at `n=9`, versus the claimed seven negative directions.
  In the THM-299 old/new block split, the seven-tile hypotenuse Schur
  complement has inertia `(5+,2-)`; the preceding `n=7,8` layers each
  contributed only one negative direction.
- **Repair / strongest survivor:** mark THM-300 `REFUTED`.  Exact finite
  evidence retains full rank for `6<=n<=12`, the old negative count for
  `5<=n<=8`, and the new pattern `negative=n-1` only for `9<=n<=12`.
  The latter remains OPEN beyond the audited range.
- **Reusable rule:** when a nested matrix family suggests constant inertia
  increment, audit the exact new-layer Schur complement at the first
  untested size; boundary concentration of eigenvectors does not fix its
  inertia.

## MISTAKE-378 (2026-08-14, THM-309 Paley cycle design) -- ordered 2-transitivity and trace moments were both overclaimed

- **What failed:** THM-309's proof called the full affine group an
  automorphism group and invoked 2-transitivity.  Nonsquare multipliers
  reverse the Paley orientation.  A nearby S24e artifact also used
  `tr(A^k)/k` as the simple-cycle count beyond its valid short-length range,
  thereby counting repeated-vertex closed walks as cycles (for example it
  reported `318` instead of `24` simple 7-cycles in `P_7`).
- **First failed implication:** the theorem asks for equal incidence of
  **unordered** vertex pairs, so ordered-pair 2-transitivity was unnecessary.
  Separately, a closed walk ceases to be forced simple once its length permits
  repeated directed subcycles; trace data cannot be inserted as `c_k` there.
- **Repair / strongest survivor:** the square-affine subgroup
  `{x -> ax+b : a square}` has order `p(p-1)/2` and acts sharply
  transitively on unordered pairs because exactly one of slopes `a,-a` is a
  square for `p=3 mod 4`.  It preserves every simple directed cycle, so the
  cycle multiset is a `2-(p,k,lambda_k)` design for **all** `3<=k<=p`, even
  as well as odd, with `lambda_k=c_k C(k,2)/C(p,2)`.  The old 5-cycle counts
  and 3-incidence analysis survive; the S24e long trace-count conservation
  table is quarantined.
- **Reusable rule:** match the group action to the incidence being balanced,
  and distinguish closed walks, primitive necklaces, and simple cycles before
  transferring a trace formula.

## MISTAKE-377 (2026-08-12, disconnected affine-ray workload) -- a many-turn skip omitted its short-rotation hypothesis

- **What failed:** the first affine-ray quotient audit and provisional carrier
  scanner applied the Dirichlet many-turn inequality whenever its displayed
  turn statistic was at least five.  The proved estimate also requires the
  short-rotation hypothesis `9|c|<=p`; neither analytic skip tested it.
- **Minimal witness / first failed implication:** at physical pair
  `(p,q)=(264,302)` and body-safe context `(L,j,e,f)=(168,90,2,1)`, the
  lawful witness `(d,a,c)=(7,1,2)` has fewer than five turns.  The unrelated
  witness `(8,1,40)` passes the provisional turn test but violates
  `9|c|<=p`, so it caused an uncertified skip.  The exact physical mass
  `4591428/225009725` is above `Dmax/5`; the error is routing, not a physical
  counterexample.
- **Repair / strongest survivor:** require `9|c|<=p` at both the universal
  and contextwise gates.  The primitive quotient `22890 -> 14168`, carrier
  chambers, and strict cutoff `14913` are unchanged.  The honest residual
  occurrence count is `8,079,264`, not `8,013,156`.  Grouping all lawful
  witnesses by physical `(p,q)` remains valid.
- **Reusable rule:** when several affine witnesses represent one physical
  pair, a witness may route that pair through an analytic theorem only after
  satisfying every hypothesis of that theorem.  A stronger-looking statistic
  cannot import a missing scale regime from another witness.

## MISTAKE-376 (2026-08-12, disconnected-low frontier routing) -- navigation advertised two reductions absent from its proof package

- **What failed:** the incoming LRC proof-map headline said the primitive
  `3:5` lane was closed at every dilation and that a repaired Dirichlet
  reduction left exactly `22,890` affine rays. Its attached theorem note and
  scripts instead exclude `3:5`, prove only the non-`3:5` `g>=4` cone, and
  contain no compiler or derivation of the `22,890` count.
- **First failed implication:** a result anticipated by concurrent work is not
  inherited evidence. The proof router cannot outrank the exact artifacts it
  routes to, especially when their declared universe expressly omits the lane.
- **Repair / strongest survivor:** retain the proved `36,520` component-profile
  count, complete-multipartite five-edge tree, `q>=8p` physical floor, and
  non-`3:5` moderate-ratio `g>=4` closure. Restore `3:5` at arbitrary dilation
  and all other small-ruler `g<=3` channels to OPEN. Treat the proposed
  `p>=264`/`22,890` reduction as exploratory until its proof and compiler land.
- **Reusable rule:** a navigation merge must be audited against the universe
  and exclusions in the routed artifact; never let an anticipated sequel
  enter the current headline before its certificate.
- **Superseding evidence (later 2026-08-12):** commit `365f6b8bd` subsequently
  supplied the two artifacts that were absent at this checkpoint.  The exact
  symbolic compiler now proves the `3:5` lane at every dilation on all 2,530
  small-ruler contexts, and the generalized Dirichlet companion proves the
  `p>=264` reduction and enumerates the 22,890 nonzero-resonance affine cover.
  This does not make the earlier routing lawful retroactively, and it does not
  certify the remaining affine rays or finite raw head; it supersedes only
  the two "artifact absent" status statements above.

## MISTAKE-375 (2026-08-12, THM-3352 integration) -- stale dependency pins made both exact replays fail before their audits

- **What failed:** the collision-safe integration of THM-3350/3352 retained
  the pre-integration SHA-256
  `32587f0b965de7da1096e0f817cee46429ed2842495790cb9d489a21d2ed24c4`
  in both THM-3352 replay scripts, although the integrated THM-3350 tail
  script, both theorem front matters, and the results index consistently use
  `78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9`.
  The argmin replay also retained the pre-integration reference-engine hash
  `b125427d...` instead of the declared integrated hash `da941a42...`.
- **Minimal witness / first failed implication:** invoking
  `lrc14_connected_low_all_heads_universal_forest_thm3352.py --limit 1`
  raised `RuntimeError(('tail hash', ...))` before constructing even one
  context.  The independent argmin replay failed first at that tail pin and,
  after it was repaired, at the stale reference-engine pin.  This was a
  reproducibility regression, not a counterexample to any overlap or forest
  inequality.
- **Repair / strongest survivor:** update only the three dependency pins to
  the already-declared LF-normalized hashes.  The repaired one-channel compiler
  reaches its partial semantic output, and the argmin replay reaches its
  canonical inventory `(4044,261254,4148)`.  The mass engine, expected
  semantic digests, theorem statements, and mathematical dependencies are
  unchanged.
- **Reusable rule:** after integrating a proof family across branches, audit
  reverse dependency pins, not only each file's own front matter.  A theorem
  can have internally consistent declared hashes while downstream replay
  scripts still pin an obsolete parent blob.

## MISTAKE-374 (2026-08-12, provisional THM-3354 response row) -- localization vanishing was mistaken for an integral polynomial mate

- **What failed:** the provisional D5 comparison record placed the unit
  observer `theta=[1]` on THM-3348's generic de Rham site while asserting the
  integral biconditional `polynomial mate iff theta=0`. This silently treated
  localization as conservative.
- **Minimal witness / first failed implication:** for
  `P=x+x^2z`,
  `D_P=(1+2xz)partial_z-x^2 partial_x`. In the localization at `x`, the
  rational primitive `Q=1/x` gives `D_P(Q)=1`, so the generic image
  `theta_gen` vanishes. But THM-3348's exact one-root formula gives
  `Ann_(K[P])(theta_int)=(P)`, hence `theta_int!=0` in
  `K[x,z]/D_P(K[x,z])` and no polynomial mate exists.
- **Repair / strongest survivor:** split the response into the integral record
  `(K[P] subset R,C_P,theta_int)`, where the mate biconditional is exact, and
  the generic record `(Spec K(P)[x,g^(-1)],H^1_dR,theta_gen)`, where a mate
  implies vanishing but the converse fails. The direct coefficient-map no-go
  of THM-3354 survives; its comparison cospan is explicitly definitional.
- **Reusable rule:** after localization, distinguish a necessary generic
  obstruction from an integral existence criterion. Record the kernel of the
  localization map--here vertical `K[P]`-torsion--before transferring an iff.

## MISTAKE-373 (2026-08-12, THM-3353 raw Gaussian toggle) -- a parent-torsor edge was stated on gauge-dependent lifts

- **What failed:** the provisional THM-3353 candidate said its rational
  Gaussian rotation changed exactly the local factor above `p` on the
  displayed signed Gaussian lift.
- **Minimal witness / first failed implication:** at `p=5,X=A,s=1`,
  `(2+i)(16+5i)=27+26i` is sent to
  `(2+i)(16-5i)=37+6i`.  Conjugating only the `5`-factor instead gives
  `(2-i)(16+5i)=37-6i`, the global conjugate of the displayed target.
- **Repair / strongest survivor:** on
  `X_C~=F_2^omega(C)/<1>`, the two outputs define exactly the same
  `p`-coordinate edge; raw lifts may differ by unit and global-conjugation
  gauges, which complement every allocation bit.  The compiler, exact
  addresses, valuation one, complementary roots, fixed-`p` transducers, and
  arbitrary-rank/dispersion theorem all survive with this quotient typing.
- **Reusable rule:** a one-prime allocation statement belongs on the parent
  Boolean torsor unless a signed associate and global-conjugation section have
  both been fixed.  Test a literal factorization before promoting a raw-lift
  claim.

## MISTAKE-372 (2026-08-12, projected `z1=216` two-high closure) -- a safe cell address was mistaken for a safe projected section

- **What failed:** an exploratory terminal counted cells whose distinguished
  grid point was safe for two high drifts and inferred that the whole aligned
  projected completion was obstructed.  A second shortcut used
  `(Z/2Z)^*={1}` to identify two denominator-two safe-cell sets.
- **Minimal witness / first failed implication:** every cell address `c` has
  projected coordinate `phi_L(c/L)=0`, and `0` belongs to every aligned danger
  set.  Thus safety of the address in the centered drift coordinate is fully
  compatible with projected containment.  For denominator two the unit is
  fixed, but the local slope `h+1/2` still depends on the ray height, so two
  safe sets need not coincide.
- **Repair / strongest survivor:** retain the local coordinate `y`.  For a
  fixed-safe multiset `C`, maximize weighted cell multiplicity over every unit
  and every translated open cyclic interval of length `d/7`; if
  `sum_i m_C^tr(d_i)>(r-1)|C|`, inclusion--exclusion supplies a cell safe from
  all `r` highs at each `y`, hence the full projected section.  The three
  denominator-two equality cases instead use a measure bound: one high removes
  at most `2/7`, so two leave at least `3/7>36/91`.  The repaired exact theorem
  is THM-3351; no point-only artifact was promoted.
- **Reusable rule:** projection statements require a witness on every local
  fibre, not a favorable section origin.  A unit-orbit collapse fixes residues,
  not affine slopes; retain translations and multiplicities until the target
  predicate has been proved.

## MISTAKE-371 (2026-08-12, THM-3346 modular/full content scope) -- one sufficient fixed-grade condition was stated as necessary

- **What failed:** the first promoted form of
  [THM-3346](theorems/THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy.md)
  said its selected `N`-primary channels could equal the full Gaussian
  contents only in the literal specialization `N=C_r` with `N|C_s`.
- **Minimal witness / first failed implication:** take `N=5,r=3,s=6`.
  Then `(C_r,C_s)=(25,85)`, while both selected and full channel pairs are
  `(1,5)` even though `N!=C_r`.  The smallest actual modular/full mismatch is
  `N=5,r=3,s=21`, where the pairs are `(1,5)` and `(1,25)`.
- **Repair / strongest survivor:** for every `r,s in R_N`, each signed
  channel satisfies `delta_pm=gcd(N,d_pm)`.  Because the two channel pairs are coprime and have
  products `N` and `gcd(C_r,C_s)`, respectively, they agree coordinatewise
  iff `gcd(C_r,C_s)=N`.  The fixed-grade condition is sufficient, not
  necessary.  The root cube, Hensel clocks, conjugation monodromy, metric,
  Pell compiler, and every other conclusion of THM-3346 are unchanged.
- **Reusable rule:** when a theorem identifies the primary part of an
  invariant, equality with the full invariant is controlled by absence of
  extra common factors, not by equality of one ambient object with the
  selected modulus.

## MISTAKE-370 (2026-08-12, THM-3336 primitive matrix extension) -- primitive-entry matrices were treated as a composition class

- **What failed:** the first promoted Section 9 introduced `d_A,mu_A` under
  the hypothesis that the four entries of `A` have gcd one, then stated the
  matrix composition cocycle without declaring its larger codomain.  Read as
  closure of the declared class, this is false; it also invites the false
  inference that the scalar-free effective determinant grades composition.
- **Minimal witness / first failed implication:** with
  `A=B=[1 -1;1 1]`, both factors have entry-content one and determinant two,
  but `AB=[0 -2;2 0]` has entry-content two.  After scalar removal the
  effective degrees are `Delta(A)=Delta(B)=2` and `Delta(AB)=1`, not four.
- **Repair / strongest survivor:** define `d_M,mu_M` on every nonsingular
  integral matrix.  Then
  `d_(ML)(u)=d_L(u)d_M(mu_L(u))` and
  `mu_M(mu_L(u))=mu_(ML)(u)` are exact.  Reserve entry-content one for the
  Smith/range conclusions; for `g=cont(M)`, primitive normalization obeys
  `mu_M=mu_(M/g)` and has objectwise effective degree
  `Delta(M)=|det M|/g^2`.  All Gaussian, Farey-face, Vieta, Boolean-groupoid,
  and determinant-gate conclusions of THM-3336 survive with this typing.
- **Reusable rule:** after dividing every output by a content, test both
  closure of the normalized class and multiplicativity of the proposed
  degree.  A content cocycle can preserve the action while destroying a
  monoid grading.

## MISTAKE-369 (2026-08-12, THM-3341 U-spine Pell synthesis) -- an unoriented Markov branch and shared Pell field were overcompressed

- **What failed:** the first promoted form of
  [THM-3341](theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md)
  said every ordered positive fixed-two Markov solution had the displayed
  positive inverse, called five negative intersection rows “four,” and let
  its typed-carrier language obscure that the negative-Pell, Pell-8, and
  norm-17 equations all live in `Q(sqrt(2))`.
- **Minimal witnesses / first failed implication:** `(2,5,1)` is an ordered
  positive Markov solution but the unoriented formula gives `R=-1`;
  `(2,1,1)` gives the algebraic boundary `R=N=0`, not a positive
  square-triangular row.  Also `sqrt(8)=2sqrt(2)` and
  `17+6sqrt(8)=(3+2sqrt(2))^2` directly refute any different-field reading.
- **Repair / strongest survivor:** modulo swapping the two non-two Markov
  coordinates, orient uniquely as `1<=M_-<=M_+`; equality `(1,1)` is the
  boundary and strict inequality is exactly the positive selector.  The
  three recurrences share `Q(sqrt(2))` and its unit, but occupy different
  norm/parity cosets and preserve different typed predicates.  Every infinite
  classification, Gaussian branch transplant, Boolean-fibre theorem, and
  norm-17 orbit in THM-3341 survives unchanged.
- **Reusable rule:** normalize unordered Diophantine coordinates and isolate
  zero/degenerate seeds before writing an inverse.  Distinct norm equations
  or typed orbits inside one quadratic field are not distinct fields.

## MISTAKE-368 (2026-08-12, AMM 12592 endpoint rigidity) -- nonattainment at slope one was confused with a strict gap above the infimum

- **What failed:** the Szegő endpoint reflection and several downstream
  syntheses inferred `C*>1` from the fact that no exactly fair extractor has
  a bounded-additive deadline `T(n)<=n+D`.  The archived Lane B draft also
  proved the stronger per-extractor statement `T(n)-n != o(n)`, but this was
  cited under a nonexistent `THM-2967` and was not on the canonical truth
  surface.
- **Why it was wrong:** nonattainment of an endpoint does not separate an
  infimum from that endpoint.  Different extractors could have deadlines
  `T_epsilon(n)<=(1+epsilon)n+D_epsilon` for every `epsilon>0` even though no
  single extractor has sublinear excess.  The toy feasible-slope set
  `(1,2]` already has infimum one without containing one.  No compactness or
  uniform positive lower bound was supplied to exchange these quantifiers.
- **Repair / strongest survivor:**
  [THM-3342](theorems/THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors.md)
  now canonizes the full valid Lane B result: every fixed feasible deadline
  obeys `limsup (T(n)-n)/n>0`.  Its Pólya--Carlson/Fatou/Kronecker proof is
  stronger than the bounded-alphabet Szegő proof, but still gives no positive
  constant uniform over extractors.  Thus the general-class infimum `C*`
  remains open with only `C*>=1`; the slope-one endpoint and all
  `n+o(n)` envelopes are unattainable.  Balanced-block lower bounds in
  THM-3009 are unaffected.
- **Reusable rule:** before turning an impossibility at a boundary parameter
  into a strict infimum inequality, exhibit either a uniform quantitative
  gap or a compactness/closure theorem.  Always write which object is fixed
  before taking a limit.

## MISTAKE-367 (2026-08-12, THM-1880 Chebyshev--Pell frame) -- the coupled recurrence transposed its even and odd outputs

- **What failed:** [THM-1880, the a/b functional frame](theorems/THM-1880-the-a-b-functional-frame-chebyshev-pell-companions.md)
  stated
  `E_n=E_(n-1)+x O_(n-1)` and
  `O_n=O_(n-1)+x E_(n-1)`, and its status incorrectly said that this
  displayed recurrence had zero symbolic residual.
- **Minimal witness / first failed implication:** from the defining forms,
  `E_1=x`, `O_1=1`, and `E_2=x^2+1`.  The old first recurrence instead gives
  `E_2=2x`.  Thus the printed formula fails already at `n=2`; a historical
  verification claim did not match the expression canonized.
- **Repair / strongest survivor:** direct multiplication of
  `E_n+O_n=(x+1)^n` and `E_n-O_n=(x-1)^n` gives the corrected crossed system

  ```text
  E_n=x E_(n-1)+O_(n-1),
  O_n=E_(n-1)+x O_(n-1).
  ```

  [THM-2142, the half-angle bridge](theorems/THM-2142-the-half-angle-bridge-ab-monoid-is-the-ctu-cyclotomic-skeleton.md)
  had already recorded this correction; THM-1880 now agrees with it.  The
  defining closed forms, Pell identity, cotangent roots, and triangular
  coefficient remain valid.
- **Reusable rule:** when a later theorem explicitly corrects a live proved
  theorem, repair the original truth surface as well.  A verifier label such
  as “residual zero” is not evidence unless the frozen expression is the one
  printed in the theorem.

## MISTAKE-366 (2026-08-03, THM-3321 Hesse-torus normalization) -- noncovariance of one named torus was promoted to a classification of all continuous normalizations

- **What failed:** THM-3321 correctly showed that the formal Hesse torus
  `z -> tz, zbar -> t^(-1)zbar` gives different weights to nonzero pure terms
  of `M_3`, but then said that only global projective scaling was a lawful
  continuous coefficient normalization.
- **Why it was wrong:** the weight check refutes covariance of that particular
  `C^*` action.  It does not classify every continuous automorphism of the
  moment ideal or rule out an unrelated symmetry.  Failure of a proposed
  action is not a universal symmetry classification.
- **Repair:** THM-3321 now says exactly that the Hesse torus supplies no second
  normalization.  Projective scaling is established and used, without a
  classification of all continuous moment symmetries.  The five homogeneous
  Macaulay certificates, projective support-`<=4`
  exclusion, and all exact computations are unchanged.
- **Reusable rule:** state a symmetry negative with the same quantifier as the
  action tested.  To classify all normalizations, compute the automorphism
  object rather than extrapolating from one failed carrier action.

## MISTAKE-365 (2026-08-03, THM-3319 base/deck type collision) -- a completed base section was identified with its rank-two quadratic cover

- **What failed:** THM-3319 reused `R_2` for both the localized gradient cubic
  and the constant coefficient of a quadratic subresultant, then said that
  the quadratic cover's completion was the formal section `(12)`.  Its prose
  also placed a first-order resultant tangent check next to an all-orders
  persistence claim without displaying the scheme-theoretic identity.
- **Why it was wrong:** `(12)` parameterizes the completed base germ
  `Spf A_i[[h_d,h_k]]`; the deck is a separate rank-two algebra obtained by
  adjoining a root of the surviving quadratic.  Vanishing of two first
  derivatives alone does not imply that a resultant vanishes identically on
  the germ.
- **Repair:** `R_2` is reserved for the gradient cubic and the quadratic
  constant is now `R_2^(0)`.  The theorem distinguishes the base completion
  from `R_i[t]/(S_2(-t))` and inserts the unit-scaled subresultant identity

  ```text
  P_2b^2-Q_2ab+R_2^(0)a^2=u Res_y(R_1,R_2).
  ```

  On `a=b=0` this proves all-order resultant vanishing; the unit-leading
  nonzero quadratic row then gives gcd degree exactly two.  The etale germ,
  connected cover, gradient obstruction, status, and scope are unchanged.
- **Reusable rule:** never identify a parameterizing base section with a
  finite cover over it.  Track algebra rank and type explicitly, and use an
  all-orders identity rather than tangent vanishing for scheme-theoretic
  persistence.

## MISTAKE-364 (2026-08-03, cyclic-quartic support-five Macaulay target) -- a dimensionally impossible degree-21 full-rank certificate was proposed as the next proof test

- **What failed:** the first THM-3321 boundary paragraph and its synthesis
  proposed full column rank of the five-variable degree-21 Macaulay map for
  `M_3,M_6,...,M_21` as the next support-five emptiness certificate.  The map
  has `13,972` rows and `12,650` columns, which made full rank look plausible.
- **Why it was wrong:** row count ignores universal syzygies.  The first five
  forms have degrees `3,6,9,12,15`; the complete-intersection Hilbert
  coefficient in degree `21` is `1705`, a lower bound for the quotient by any
  five such forms in five variables.  `M_18 R_3` removes at most `34`, not
  `35`, dimensions because the `M_3 M_18` Koszul overlap is already zero in
  the quotient, and `M_21` removes at most one.  Thus every such map has
  quotient dimension at least `1670` and rank at most `10980`.  Full degree-21
  rank is impossible independently of the Hesse coefficients.
- **Repair:** THM-3321's support-`<=4` theorem is unchanged, but its open
  support-five boundary now records the rank bound.  The formal product has
  coefficient `39` at degree `28` and first becomes nonpositive at degree
  `29`; the raw degree-29 map is `66,486 x 40,920`.  Degree 29 is only the
  first candidate not excluded by this count, not a generic-Hilbert-series
  theorem.  Projective emptiness needs a sufficient-degree certificate or a
  different saturation/chart argument.
- **Reusable rule:** before budgeting a Macaulay rank computation, compute the
  Hilbert/Koszul lower bound.  More rows than columns do not imply that full
  column rank is algebraically possible.

## MISTAKE-363 (2026-08-03, simplex moments mod p) -- a modular reduction silently corrupted every moment because the prime divided the factorial denominator

- **What failed:** while searching the degree-four cyclic eigenspace on
  `Delta_2` for a homogeneous-Factorial-Conjecture counterexample, an
  exploratory cascade reduced the uniform-simplex moment
  `<s_1^a s_2^b s_3^c> = 2 a! b! c! / (D+2)!` modulo small primes `p = 7, 13`.
  It reported `743` surviving points of `P^4(F_7)` at the very first condition
  and identical survivor counts at `m = 3, 6, 9`, which is not how a cascade of
  independent conditions behaves.
- **Why it was wrong:** the test polynomials have barycentric degree up to 4,
  so `<g^m>` involves monomials of total degree `D = 4m` and denominators
  `(4m+2)!`.  For `p = 7` already `(4*3+2)! = 14!` is divisible by `7`, so the
  denominator is `0 mod p`, the naive `pow(den, p-2, p)` returns `0`, and every
  moment silently becomes garbage.  Sixteen distinct moment denominators
  through degree 12 die at `p = 7`.  No exception is raised; the run simply
  produces confident nonsense.
- **Repair:** the reduction is sound only when `p > D + 2`, i.e. `p > 4m + 2`.
  `THM-3310` enforces this with an explicit `require`, refuses to invert zero,
  and records the thresholds `m = 3, 6, 9, 12, 15, 18 -> p >= 19, 31, 43, 61,
  67, 79`.  The published modular certificate in `THM-3300` was audited against
  this and is unaffected: it uses barycentric degree 3 with `m <= 9`, so
  denominators are at most `29!` against `p = 10^9 + 9`.
- **Genus, and where else it can recur:** any reduction of a rational
  combinatorial constant -- simplex or sphere moments, Beta and multinomial
  weights, Bernstein coefficients -- must check that the prime exceeds every
  factorial in the denominator.  Prefer a prime far larger than the maximal
  degree, and make the inverse routine fail loudly on zero rather than return
  it.

## MISTAKE-362 (2026-08-03, exceptional quadratic geometric count) -- a relative two-point fibre was reported as the total geometric base change

- **What failed:** the first frozen output and reflection for the affine-`c`
  exceptional quadratic printed `geometric_points=2` and said that the
  degree-72 closed point becomes two points over an algebraic closure, without
  specifying the base field.
- **Why it was wrong:** the quadratic algebra `B_i/A_i` has relative degree
  two, but `A_i/K_i` has degree `36`.  After fixing one geometric embedding of
  `A_i`, the relative fibre has two directions exchanged by `C_2`; after base
  change from `K_i` to an algebraic closure, all `36` base embeddings split,
  giving `72` geometric points in `36` conjugate pairs.
- **Repair:** the transcript now records both
  `relative_geometric_fibre_points=2` and
  `total_geometric_points_over_K=72`; the reflection and synthesis use
  fibrewise deck language.  The irreducibility, nonsquare norm certificates,
  degree `2/72`, lack of an `A_i`-rational direction, and first-normal
  nonstationarity are unchanged.  Reusable rule: every geometric-point count
  must name the field being algebraically closed.

## MISTAKE-361 (2026-08-03, THM-3024 general-class floor promotion) -- a truncation-edge Hall cut was read as a cross-shell floor

- **What failed:** THM-3024 promoted the balanced-block archimedean floor
  `C*_block >= log_5(5 phi^2)` to ALL exactly fair AMM 12592 extractors via a
  forward-routing transportation model, citing (G1) an aggregate Hall-cut sign
  flip at golden and (G3) numeric equality of degree-resolved and per-shell
  cuts. The cited script computed only per-shell continuum (ARCH) margins in
  floats; no cross-shell cut was ever computed. In exact arithmetic, within
  the theorem's own model (forward routing at preserved absolute degree),
  every tail cut with a deeper shell available is satisfied for ANY
  `gamma > 0`: at fixed degree `d`, demand `binom(m-1,d)` is outrun by supply
  exponentially. Independently confirmed twice (audit agent + orchestrator):
  at `gamma = 71/125 < gamma*`, the genuine per-shell deficit `2^242` at
  `(m,d) = (256,155)` is absorbed by shell `512` with `~2^128` room. The
  reported binding cuts were truncation-edge artifacts — the deepest
  windowed shell's own per-shell constraint, i.e. the balanced-block floor
  read back circularly.
- **Why it was wrong:** an aggregated Hall cut is only as strong as the
  model's routing bound. With an unbounded forward window the relaxation has
  no floor at all, so the "most generous possible degree mobility" reading
  inverts the logic: generosity destroys the cut rather than validating it.
- **Correct framing:** the general-class `C*` floor is OPEN again. Surviving:
  balanced-block `C* > 1.5970` exact through `m = 4096` and the asymptotic
  block/checkpoint barrier `log_5(5 phi^2)` modulo one unwritten Stirling
  transfer lemma (THM-3009 sec 10.3); THM-3027's tangency collapse, now
  scan-free in the floor direction (concavity + gamma-monotonicity upgrades).
  The missing ingredient for a general floor is a DEADLINE-BOUNDED routing
  window derived from the extractor axioms — the pathwise deadline is exactly
  what the transportation relaxation forgot. Audit with exact certificates:
  `05-knowledge/results/amm12592-golden-floor-audit-boxeph.md`.

## MISTAKE-360 (2026-08-02, THM-3214 offset-six application) -- canonical PRS rows were identified literally with unnormalized iterates

- **What failed:** the first promoted wording of THM-3214 equation `(25)`
  identified the raw second pivot and connection literally with
  `P_2(r,a)` and `P_1(s,r)`.  The fraction-free operator is homogeneous:
  rescaling an earlier row propagates accumulated unit powers into later
  iterates.  The displayed `r,s` in THM-3192 are canonical renormalized PRS
  rows, not the unnormalized iterates of THM-3214 equation `(7)`.
- **Why it did not invalidate the theorem:** the accumulated factors are
  `p`-units on the stated `p>=197` range.  They neither change the `H,J,K`
  chart ideals nor the `2,4,5` jet-locality budget, and the universal locality
  and Catalan sharpness statements never use the mistaken literal equality.
- **Repair:** state `(25)` after the standard `p`-unit row renormalizations and
  display the raw scales
  `rho_2=U_H^2P_2(r,a)` and `chi_2=U_H^3U_JP_1(s,r)`, then identify the
  normalized coordinates by their generated ideals and jet orders.  Equality
  of principal ideals is not equality of raw representatives.  Catalan
  sharpness is universal pseudo-division sharpness, not necessity inside the
  thinner factorial-moment family.

## MISTAKE-359 (2026-08-02, THM-3183 top-jet evidence) -- an exact symbolic helper admitted a float and constructed two truncated rows without their predecessors

- **What failed:** THM-3183's maintained companion described its offset-six
  top-jet calculation as exact, but the codimension-zero branch formed
  `1/1` with Python division and therefore introduced `Float(1.0)`.  The same
  helper constructed the unused lower entries `R[0],S[0],S[1]` after replacing
  unavailable predecessor coefficients by zero.  For example, at `p=5` its
  displayed internal `R[0]` was `1428382771200`, whereas direct multinomial
  expansion gives `-36838416384000`.
- **Why it did not invalidate the theorem:** the asserted and independently
  audited THM-3183 identities use only the leading entries `R[3]` and `S[2]`.
  Those entries require no omitted predecessor, agree with direct integer
  expansion, and still factor by the stated `H` and `J`.  No proved matrix,
  Smith, continuant, wall, or PRS-leading formula changes.
- **Repair:** use an explicit SymPy `Rational` factorial normalization, allow
  exact negative factorial shifts relative to `(2p)!`, extend `A,B,R` far
  enough that every constructed row has its true predecessor, and compare all
  66 resulting `A/B/R/S` entries directly with integer multinomial sums at
  `p=5,7,11`.  The reusable rule is that an unused symbolic row is still part
  of the evidence surface: either construct it lawfully or do not construct
  it, and never infer exactness merely because a `1.0` coefficient is
  mathematically integral.

## MISTAKE-358 (2026-08-02, THM-3169 stale post-QED status) -- a promoted theorem retained its candidate disclaimer

- **What failed:** THM-3169's frontmatter and audit record correctly marked the
  depth-six certificate `PROVED + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED`, but its final line still said
  `QED (candidate pending independent audit)`.
- **Why:** the promotion changed the maintained theorem status without removing
  the provisional EOF suffix.  That contradiction is especially harmful here
  because the exact depth-seven successor uses THM-3169 as its immediate
  dependency.
- **Repair:** replace the stale suffix by plain `QED.`.  No mathematical claim,
  dependency, script, output, or hash changes.

## MISTAKE-357 (2026-08-02, THM-3163 stale post-QED status) -- a promoted theorem retained its candidate disclaimer

- **What failed:** THM-3163's frontmatter and headline were correctly promoted
  to `PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`, but its final
  line still said `QED (candidate pending independent audit)`.
- **Why:** the promotion patch updated the maintained status surfaces but
  missed the provisional suffix at EOF.  The contradiction could make later
  consumers treat an audited theorem as unavailable even though its proof and
  exact evidence were already accepted.
- **Repair:** replace the stale suffix by plain `QED.`.  No mathematical claim,
  dependency, script, output, or hash changes.

## MISTAKE-356 (2026-08-02, THM-3162/3164/3167 live reservation races) -- non-atomic repairs can collide again

- **What failed:** the depth-six selector non-resurrection session reserved
  `THM-3162` from a fresh snapshot, but the unrelated falling-factorial
  order-join session had already added its `THM-3162` file to `main` 51 seconds
  earlier.  The crossed pushes left two distinct files with the same YAML ID.
- **Why:** a fetch-and-scan reservation is not atomic.  The selector stub was
  mathematically honest but chronologically second, so its namespace claim
  was invalid under the first-on-`main` rule.
- **Failed first repair:** both sessions independently reacted to the
  `THM-3162` collision by choosing `THM-3164`.  The proved order-join packet's
  repair landed first; the selector repair arrived 12 minutes later and
  recreated the duplicate under a new number.
- **Second failed repair:** a later fresh scan made `THM-3167` appear free, but
  an already-pushed inverse-different repair owned that ID in newer remote
  history.  The selector packet had not yet been pushed under this name.
- **Final repair:** the already audited order-join theorem keeps `THM-3164`,
  and the inverse-different theorem keeps `THM-3167`.  The later selector
  packet is coherently renamed to `THM-3169`, including theorem, script,
  output, hashes, and headings, before its proof is promoted.  This is a
  namespace repair only; none of the theorems' mathematics changes.

## MISTAKE-355 (2026-08-02, THM-3159 candidate audit) -- losing the odd reflected-pole sign turned the zero-face system into an unrelated quartic gcd

- **What failed:** the first THM-3159 candidate converted a partial fraction
  written with denominator `(beta*t-1)^3` to the normalized pole
  `(1-beta*t)^3` without the required minus sign.  Solving that incorrectly
  signed `F_2=F_3=0` system produced a quartic numerator `N` and the
  degree-`q+2` test polynomial `x^2(x-1)U(x)-N(x)`.
- **Why it was load-bearing:** exact coprimality of that quartic transform did
  not imply coprimality of the actual two endpoint faces.  Re-expanding the
  partial fractions gives
  `(beta*t-1)^(-3)=-(1-beta*t)^(-3)` and changes the solved values to
  `U(alpha)=-1/(2alpha)` and
  `U(1-alpha)=-1/[2(1-alpha)]`.  Thus the lawful transformed polynomial is
  `P(x)=2xU(x)+1`, of degree `q`, not the old quartic expression.
- **Repair / strongest survivor:** at `q=249727`, a fresh exact FLINT gcd and
  extended-gcd verification give `gcd(P(x),P(1-x))=1`; the singular charts
  `alpha=0,1,1/2` are checked separately.  THM-3159 therefore survives, but
  its theorem text, script, transcript, and hashes were replaced before
  promotion.  The old candidate SHA `67450bdb4` is not a proved dependency.
- **Reusable rule:** when normalizing partial fractions, parity is data.
  Convert every pole denominator before collecting residues, and audit the
  solved linear system against a direct finite-field face evaluation before
  launching a large resultant or gcd computation.

## MISTAKE-354 (2026-08-02, THM-3155 selector-barcode scope audit) -- treating abstract terminal-law realization as a substantive sequential obstruction

- **What failed:** THM-3155's scope said that an arbitrary law on legal
  unordered pole prefixes need not be the stopping distribution of a
  sequential pole-removal process.  Without further restrictions on the
  transition kernel, that statement is false.
- **Why:** for any finite terminal law, sample the terminal labelled subset and
  reveal a uniformly random ordering.  Conditional probabilities given the
  current prefix depend only on that prefix, so they define a state-dependent
  one-letter-at-a-time Markov chain with exactly the prescribed terminal law.
  Equal-value symmetrization descends to unordered multiplicity states.
- **Repair:** THM-3155 now distinguishes automatic abstract Markov realization
  from the missing value-only, prescribed-hazard, selector-current-compatible,
  or original-response-compatible transport.  THM-3163 records the exact
  posterior kernel, its proof, and the labelled lift of THM-3158's law.  No
  barcode, Hasse-positivity, NC2, or GMC conclusion changed.

## MISTAKE-353 (2026-08-02, reflected two-star replay) -- a semantic-preserving upstream normalization silently invalidated raw dependency pins

- **What failed:** the all-`649`-body upper-median two-star referee stopped at
  import time before checking any of its `55,606,320` assignments.  It pinned
  raw hashes of the low-phase and robust-edge-11 referees, but their later
  dependency refresh and line-ending normalization changed those bytes while
  preserving their proved statements and semantic digests.
- **First failed implication / exact witness:** a previously stored `PASS`
  transcript does not imply that the current dependency graph replays.  The
  exact witness was the stale low-phase pin
  `b2418d...`, whereas the current audited source is `416c36...`; refreshing
  that one pin exposed a second stale edge-11 source/output pair.
- **Repair:** refresh all three dependency pins together, regenerate the
  transcript, and rerun both ordinary and optimized Python.  Both modes again
  exhaust the declared `55,606,320` assignments with zero failures and the
  unchanged semantic digest `aa7f4927...`.  This is a reproducibility repair,
  not new mathematical evidence and not a change to the finite-exact scope.
- **Reusable rule:** dependency normalization is still a dependency change.
  A downstream exact proof must pin the whole current source/output chain and
  be replayed after any byte-level refresh; matching an old semantic payload
  is a conclusion of the replay, never a substitute for it.

## MISTAKE-352 (2026-08-02, modular free-factor synthesis) -- an invariant cyclic substitution was called a `C3` action, and the four-point torsor was mistaken for a faithful modular carrier

- **What was claimed:** the first version of the modular free-factor reflection
  said that THM-3121/3134 cyclic substitution supplied an actual order-three
  action on the full endpoint jet, and treated a pointed involution on the
  THM-2950 `V4` torsor as a candidate faithful `C2*C3` carrier.
- **First failed implication / exact witness:** the genuine `C3` action is
  rotation of the **labelled input triple**; the numerical profile/jet output
  is invariant under it and forgets the rooted slot (even though equal-factor
  substitutions have a genuine cyclic-wreath automorphism).  On the four-point torsor a
  chosen pair flip `S` and oriented pair-cycle `R` satisfy the extra relation
  `(SR)^3=1` and generate `V4 semidirect C3 = A4 = PSL2(F3)`.  In
  `PSL2(Z)`, `(SR)^3=T^3` is nontrivial; for example it sends slope `(0,1)`
  to `(3,1)`, although those slopes coincide modulo `3`.
- **Repair / strongest survivor:** the exact intrinsic object is the
  level-three congruence shadow, not a faithful modular action.  A genuine
  lift needs both a rooted `C3` block and a `Gamma(3)`/Farey sidecar carrying
  the lost `T^3` coordinate, plus a common physical atom.  The corrected
  reflection now states this and scopes `C2*C3` only as a search grammar.

## MISTAKE-351 (2026-08-02, live theorem-reservation races) -- a fresh fetch is only a snapshot, and a reported collision is useless unless it gates the commit

- **What happened:** two independently valid results acquired `THM-3130`
  while their reservation pushes crossed.  The later prime-resonance
  factorial candidate was repaired to `THM-3131`.  A second race then put the
  tournament endpoint-jet theorem and the fixed-reference Hasse no-go at
  `THM-3134`.  Before the later owner completed its repair, the tournament
  theorem moved to the apparently free `THM-3135`; the other session made the
  same move from its own snapshot, creating a second collision.  A third race
  later put the first-on-main quartic congruence reservation and a later
  prime-power factorial stub at `THM-3141`.  The reservation command actually
  printed the collision after fetching, but `grep ... || true` was chained to
  the commit, so the discovery did not stop the write.
- **First failed mechanism:** `fetch; search; choose next; push` is not an
  atomic reservation.  It cannot see a reservation pushed after the search,
  and a diagnostic search is not a safety check when its result is ignored.
  Once a collision exists, “both choose the next free integer” has the same
  race and is not a repair protocol.  The mathematics and file contents were
  sound; the YAML namespace was not.
- **Repair / authoritative allocation:** first-on-`main` provenance decides
  ownership.  The current distinct allocation is `THM-3130` for the
  divisor-antichain response theorem, `THM-3131` for prime-resonance
  factorial closure, `THM-3134` for the tournament endpoint-jet theorem, and
  `THM-3136` for the later fixed-reference Hasse no-go.  (`THM-3135` was free
  at that checkpoint and now names the directed-cycle LRC theorem.)  The later
  prime-power factorial stub moved from `THM-3141` to `THM-3142`, leaving the
  first-on-main quartic reservation fixed; it was later promoted as the
  congruence-shadow theorem.  Later allocations include promoted `THM-3143`
  for the two-step-prime Euclidean--Newton lane and reservations `THM-3144`
  and `THM-3146` for their declared lanes.  After every
  reservation push, fetch again and search YAML IDs on the resulting remote
  history.  If a collision appears, inspect the two add commits, keep the
  first-on-`main` ID fixed, and move only the later reservation.  Repeat the
  global ID check after the repair lands; do not move both sides speculatively.
  Run the collision check as a separate hard gate: never append `|| true` and
  then continue to commit in the same shell command.

## MISTAKE-350 (2026-08-02, factorial-conjecture type audit) -- indexing SFC by slot count and restricting FC to homogeneous polynomials

- **What was claimed:** THM-3018 defined `FC(m)` only for homogeneous
  polynomials and promoted its simplex polar reduction as an equivalence for
  the full Factorial Conjecture.  A separate univariate lane then wrote
  `SFC(N)` when `N` was the number of monomial slots, so in particular a
  three-term polynomial in one variable was repeatedly called `SFC(3)`.
- **First failed implication / minimal typed witnesses:** in Edo--van den
  Essen, arXiv:1304.3956v2, Definition 2.1 and Conjecture 2.4 quantify
  `FC(m)` over **every** polynomial in `m` variables.  The legal input
  `1+x` has no single homogeneous degree and hence no single radial Gamma
  factor, so the homogeneous polar identity cannot be a global equivalence.
  Definition 2.7 and Conjecture 2.8 likewise use `m` for ambient variable
  dimension and `N(f)` separately for the number of nonzero monomials and
  the length of each shifted window.  Thus `1+x+x^2` is a three-slot slice
  of `SFC(1)`, while `1+x+y+z` is an ambient `SFC(3)` input with a
  four-moment window.
- **Repair / strongest survivors:** THM-3018 now states only the proved
  exponential-integral identity and the exact **homogeneous-subclass**
  simplex/Bernstein dictionary; its claimed all-homogeneous Laplace closure
  remains `AUDIT-REQUIRED` because of the already recorded oscillatory-saddle
  gap.  The univariate exact results are retyped as `N`-slot restrictions of
  `SFC(1)`; their polynomial systems, Macaulay certificates, resultants, and
  finite ranges are unchanged.  Historical filenames and output labels such
  as `sfc3_*` are provenance only and must be read through this correction.
  Full `FC(3)` and ambient `SFC(3)` remain open.
- **Rule:** always record the pair `(ambient dimension m, support size N(f))`.
  Never use the same index for both, and never pass from a homogeneous polar
  slice to the full nonhomogeneous conjecture without a map that preserves
  all powers of the factorial functional.  When a short label is needed,
  write `SlotSFC_1(N)` for the `N`-monomial restriction inside ambient
  `SFC(1)`.

## MISTAKE-349 (2026-08-02, THM-3101 post-promotion module audit) -- treating a normal-variable lower quotient as a finite base algebra

- **What was claimed:** after quotienting the remotely perturbed lower moment
  forms, the first repaired proof called the result a finite real algebra
  `A_C` and the subsequent upper quotient a finite-free `A_C`-algebra.  It
  then described multiplication by `h_m` as a scalar extension of an
  `o(tau)` operator on `A_C`.
- **First failed implication:** the perturbed lower forms contain the normal
  variables `v`.  Their quotient is finite free over `R[v]`, not finite over
  `R`; imposing the upper equations is a quotient, not a scalar extension.
  Thus the displayed finite-free tower did not literally justify the
  operator estimate, even though the scheme-zero mechanism remained
  available.
- **Repair / strongest survivor:** keep the lower quotient as a finite-free
  `R[v]`-module.  The exact scheme-zero identity
  `e_0 h_m=sum a_r H_r` becomes
  `e_0 h_m=-sum a_r delta_(r,C)=o(tau)` after the lower deformation.
  THM-3093 bounds the upper roots on finitely many compact composition
  charts, so fixed monomial reduction preserves this estimate after the
  upper quotient.  Pure upper perturbations cannot create `h_m`, because
  the identity already vanishes before any upper relation is imposed.  The
  An independent hostile re-audit verified the `R[v]` finite-free bases,
  bounded upper reduction, and the spectral projector formed only after the
  finite upper quotient.  The repaired theorem is promoted.

## MISTAKE-348 (2026-08-01, THM-3097 post-promotion scope audit) -- pointwise good supports were given the remote family's uniform condition-number bound

- **What was claimed:** after proving that every sufficiently translated
  fixed-width support has a normalized resultant uniformly bounded away from
  zero, Section 7 said that the normalized Bezout map has a width-only
  smallest-singular-value floor on every resultant-good support.
- **First failed implication:** a good extension may approach a gap face over
  a bad child.  Every finite extension can have nonzero resultant while the
  limiting lower-child resultant, and hence the limiting singular value,
  vanishes.  Pointwise surjectivity does not make the full good-support family
  a compact subset of the surjective locus.
- **Repair / strongest survivor:** the width-only right-inverse bound is scoped
  to the remote translated family `N>=N_tr(t)`, whose compactified closure was
  actually proved positive.  Every other good support still has the literal
  complete-intersection certificates, but no uniform norm bound.  The same
  repair checkpoint displays both row and column potential multipliers in the
  block line matrix and describes a gap face by its lower-child zero-locus
  support with the correct resultant multiplicities.

## MISTAKE-347 (2026-08-01, reflected cone split-tail orientation) -- reversing both the ratio interval and the ordered label pair covered one ordering twice

- **What was done:** the reflected one-cone referee split the harmonic-body
  tail into a low-ratio lane with ordered pair `(1,0)` and a high-ratio lane
  with ordered pair `(0,1)`.  The cap-`7/3`, cap-`5/2`, and proposed cap-`3`
  refinements inherited the same pattern and treated the two lanes as a cover
  of both relative-level orientations.
- **Minimal witness / first failed implication:** for an ordered pair `(i,j)`
  the script's ratio is `Q/P=q_j/q_i`.  The low lane `(1,0)` with `Q/P<=1`
  says `q_0<=q_1`; the high lane `(0,1)` with `Q/P>=1` says the same thing.
  Thus the admissible assignment `q_0=2,q_1=1` lies in the cap-two interval
  but in neither asserted lane.  Exhausting every channel *inside* each lane
  did not prove that the lanes cover every assignment.
- **Repair / strongest survivor:** the cap-`3` half-cone claim is quarantined,
  and the promoted one, three-quarter, and two-thirds refinements are
  `AUDIT-REQUIRED` until an assignment-level orientation cover is proved.
  The last unaffected assembled cone is the full-interval cap-`7/4` theorem
  `D>=6, 3m>=4D`; hence the current proved reflected certificate-failure wedge
  is `561` bodies with `D>=6,1<=m<4D/3`.  Future split tails must assert a
  disjunction on the original level assignment, not merely enumerate two
  oriented channel lists.  THM-3135 later supplies exactly that disjunction
  for `H2` at caps `7/3` and `5/2` by keeping one pair fixed and using both
  halves; it also proves that the eligible standard single-pair uniform lanes
  for `H` form a DAG, so the global cone claims remain `AUDIT-REQUIRED`.

## MISTAKE-346 (2026-08-01, concurrent THM-3063 reservations) -- two distinct results acquired the same new theorem ID

- **What happened:** the audited quartic cofactor-blindness theorem was frozen
  as THM-3063 on a topic branch while a concurrent session reserved THM-3063
  on `origin/main` for the five-slot terminal-suspension theorem.  Both landed
  in one merge before the cross-branch identifier check was repeated.
- **Why it was wrong:** filename checks in one worktree do not see an unmerged
  remote topic reservation, and a clean fetch cannot prevent a reservation
  pushed after that fetch.  The mathematical statements were distinct, but
  the shared YAML ID made theorem references ambiguous.
- **Repair / rule:** the cofactor-blindness theorem, companion, output, and
  current routes move together to THM-3066; the empty THM-3063 suspension stub
  retains its first-on-main namespace.  Before final promotion, recheck IDs
  against both `origin/main` and all fetched topic refs, and repeat after the
  last merge; a reservation is not complete until its ID is globally unique.

## MISTAKE-344 (2026-08-01, THM-3052 pre-promotion evidence audit) -- a final wording repair changed the frozen transcript length without updating the theorem's byte count

- **What was done:** THM-3052's status-neutral evidence repair changed the
  canonical output wording and refreshed the output artifact and hash, but
  the theorem retained the earlier claim that the LF transcript had `22,683`
  bytes.
- **Why it was wrong:** the immutable stored transcript, both fresh disjoint
  ordinary/optimized rebuilds, and the declared SHA all agree on `22,664` LF
  bytes and `96` lines.  The stale count differed by exactly the final
  wording repair's `19`-byte contraction; no mathematical census changed.
- **Correct framing / rule:** the theorem now records `22,664`.  Whenever a
  frozen transcript is regenerated, refresh and independently compare its
  byte count, line count, and content hash as one evidence tuple; a matching
  hash does not excuse stale human-readable metadata.

## MISTAKE-343 (2026-08-01, THM-3043 evidence audit) -- zero measure was still labelled COVERED, and a composite-denominator witness family was overstated

- **What was done:** THM-3043's prose correctly repaired the implication
  `mu(Safe)=0 => Safe=empty`, but its stored transcript still labelled all
  zero-measure rows `COVERED`.  The committed 58-line companion generated
  only the first census block and could not reproduce the stored quantisation
  transcript.  The theorem also wrote the `(1,...,7)` witness set as `a/8`,
  suggesting every nonzero numerator.
- **Minimal witness / first failed implication:** at `t=2/8=1/4`, speed `4`
  has distance zero, so it is not safe.  The exact safe set is only
  `{1/8,3/8,5/8,7/8}`.  It still has measure zero and is nonempty, so this is
  precisely the MAX-not-MEAN boundary rather than a counterexample to LRC.
- **Exact repair / strongest survivor:** the replacement companion performs
  the full finite rational census, explicitly enumerates a nonempty safe-point
  set for every zero-measure row, checks all six displayed witness sets, and
  reproduces the quantisation samples in normal and optimized modes.  Every
  transcript label is now `TIGHT_NONEMPTY`.  The invariant surviving the
  numerator repair is sharper: every displayed witness has **reduced**
  denominator exactly `n+1`.
- **Rule:** a measure computation cannot certify pointwise covering.  Store an
  explicit witness/emptiness certificate, and expand shorthand `a/q` whenever
  a composite denominator makes the allowed numerator set load-bearing.
## MISTAKE-345 (2026-08-01, reflected four-thirds cone proof candidate) -- monkeypatching the projective cap did not widen a hard-coded primitive-channel generator

- **What was done:** the first four-thirds referee set the imported engine's
  `RATIO_CAP` from `5/3` to `7/4`, then described its primitive bank as complete
  on `[4/7,7/4]`.  But `primitive_universe` still selected rows with the
  literal predicate `3Q<=5P`.  Its `2,492`-row bank therefore omitted every
  unordered channel in `(5/3,7/4]`.
- **First failed implication / minimal witness:** the reduced endpoint
  `(P,Q)=(4,7)` satisfies `Q/P=7/4` and lies far below the computed product
  bound, yet it was absent.  Changing a global consumed by the CSP embedding
  did not change a separate literal inside the alphabet generator.  Thus the
  original `530/31` verdict was not complete evidence even though it happened
  to be numerically correct.
- **Exact repair / strongest survivor:** the repaired generator accepts the
  cap explicitly, checks every emitted row against it, and asserts inclusion
  of the rational endpoint whenever its product is within the exact bound.
  The complete bank has `2,728` rows.  Re-running both CSP search orders gives
  the same `530` closed bodies and same `31` traps; all `158` finite policies,
  `314` tail-head controls, and direct reflected comparisons still pass.
- **Rule:** a monkeypatched geometric/search parameter is not evidence that
  every upstream enumerator consumes it.  Thread changed bounds as explicit
  arguments and freeze a boundary witness (here `(4,7)`) before trusting an
  enlarged-domain census.

## MISTAKE-340 (2026-08-01, THM-1254 Lean full-invoice consumer) -- integer tooth addresses were generalized to rationals across a discrete gap step

- **What was done:**
  `LRCCoherentBlockerChronology.binary_speed_descent_same_edge_full_invoice`
  typed the carrier address `c`, gap address `k`, and binary relative digit
  `delta` as rationals.  Its reflected branch then tried to infer
  `0<=c-k-delta` from `0<=k<c` and `delta in {0,1}`.  The paper theorem and
  the preceding Lean address lemmas all use integer tooth addresses.
- **Why it was wrong / minimal witness:** over the rationals take
  `(c,k,delta)=(1/2,0,1)`, marked positions `(a,aNext)=(0,1)`, and reflected
  data `n0R=nrR=DeltaR=s0R=1`, `residualR=1/2`.  Complete the unused original
  data by `n0=nr=Delta=s0=1`, `residual=2`.  Then every identity and positivity
  hypothesis holds, with the reflected identity following from
  `residualR=1+1*(1/2-0-1)`, but the asserted invoice would require
  `1<=1/2`.  The original disjunct is unavailable because `aNext<a` is false.
  Thus the theorem statement itself, not merely its attempted tactic proof,
  was false at the generalized rational type.
- **Exact repair / strongest survivor:** only `c,k,delta` in the full-invoice
  consumer are restored to `Int`; their two nonnegative factors are then cast
  into the rational residual identities.  All speeds, drifts, and residuals
  remain rational.  There are no downstream calls to the theorem, direct Lean
  elaboration passes, and its axiom report is exactly within
  `{propext, Classical.choice, Quot.sound}`.  THM-1254's paper statement and
  mathematical consequence are unchanged.
- **Rule:** when a proof uses `k<c => k+1<=c`, the discreteness is a load-bearing
  hypothesis.  Keep addresses in `Int` (or state the unit-gap hypothesis
  explicitly) until after that step; casting an address identity into a field
  does not license field-typing its order argument.

## MISTAKE-341 (2026-08-01, THM-3001 promotion audit) -- an asymptotic two-end expansion was read as an exact finite curvature sign

- **What was assumed:** from
  `log(R_k/R_(k-1))=C(mu_d)d^-2+O(d^-3)` and its reversed analogue, the
  candidate asserted the finite-degree necessary screen
  `C(mu_d)>=0>=C(mu_d*)`.  It then said any positive reciprocal curvature at
  one width refutes global no-return.  The same section slid from the proved
  conclusion "the ratio sequence is constant" to Newton equality `R_k=1`.
- **First failed implication / sharp hostile:** multiplying a nonnegative
  circuit by `d^2` gives only `C(mu_d)>=-O(1/d)`; at the other end it gives
  `C(mu_d*)<=O(1/d)`.  A curvature of size `1/(2d)` can be positive while an
  allowed `d^-3` remainder keeps the reversed circuit nonnegative.  Separately,
  the positive coefficient polynomial with normalized coefficients
  `h_k=2^(-k(k-1)/2)` has the constant ratio `R_k=2`, so a constant ratio need
  not mean equality in every Newton inequality.  Constant-ratio and
  Newton-equality are different conclusions; the former is all that the
  reversal-closed class argument supplies.
- **Exact repair / strongest survivor:** audited THM-3001 now states the
  quantitative screen and only concludes `liminf C(mu_d)>=0` and
  `limsup C(mu_d*)<=0`, or exact limiting signs under a fixed margin.  A
  reciprocal curvature bounded below by a positive constant still refutes
  eventual no-return.  The audit also proved the exact converse fixed-locus
  law: ratio palindromy is equivalent to
  `N*(x)=A^-1 N(x/B)`, its circuit is antipalindromic, and every odd-degree
  reversal-equivariant path hits the central circuit wall.
- **Rule:** an `O(d^(-r-1))` remainder leaves an `O(1/d)` ambiguity after
  rescaling the leading invariant.  State finite signs only with a remainder
  bound or a fixed margin; otherwise state liminf/limsup.  Do not silently
  upgrade a constant Newton-ratio sequence to the equality sequence `R=1`.

## MISTAKE-342 (2026-08-01, THM-3030 hostile audit) -- finite slot fitting was promoted as an all-order law, and its next test was assigned to the terminal slot

- **What was assumed:** the `j<=8` corner tables were said to prove the unique
  all-order odd-slot sign/constant laws.  The canonical companion also described
  two disjoint interpolation grids and out-of-sample checks that it never runs,
  and the theorem named `P_9` as the next test of `c_5`.
- **First failed implication / evidence defect:** finitely many values never
  uniquely determine an unrestricted sequence.  The current script only reloads
  frozen tables; neither grid engine nor its referenced pickle is stored.  More
  sharply, `c_m` occupies slot `k=2m-1` under the explicit hypothesis `k<j`.
  For `m=5`, `k=9` is the exceptional terminal slot of `P_9`, so it supplies no
  test at all; the first nonterminal occurrence is `P_10`.
- **Exact repair / strongest survivor:** audited THM-3030 is scoped
  `FINITE-EXACT TABLE-INTERNAL` through `j=8`.  A new independent companion
  checks all eight table hashes, `48` visible odd slots, `36` even-zero slots,
  and the corrected `(-1)^(j+m)` sign.  It also discovers the exact finite
  identity
  `(-1)^(j-1)jC_j=46 sum_(s=1)^(M-1)s^j+20M^j+K_j`, so Faulhaber proves on the
  visible range `c_m^C=46|B_(2m)|/(2m)!`.  This closes the observed denominators
  `6,360,15120,604800` without pretending to prove the continuation.
- **New decisive tests:** the Bernoulli continuation preregisters
  `c_5^C=23/23950080` at `P_10`.  At `P_12` it predicts
  `c_6^C=15893/653837184000`, with `15893=23*691` from
  `B_12=-691/2730`; this is the first test that separates the Bernoulli law from
  the observed reduced-numerator-`23` extrapolation.
- **Rule:** distinguish an exact finite atlas from an all-order theorem; keep
  reported build controls separate from executable evidence; and compute the
  first legal index from every strict slot inequality before advertising the
  next experiment.

## MISTAKE-339 (2026-08-01, pre-promotion THM-3000/3003 audit) -- the leading third jet was charged as a remainder and an asymptotic threshold was reported as an exact finite bound

- **What was assumed:** the first THM-3000 candidate imposed
  `m_j/m_1^j=o(d^(j-3))` starting at `j=3`, described uniform boundedness as
  its `o(1)` case, and promoted a decimal finite-width threshold for the third
  edge as though it were uniform on the whole `(x,z)` box. THM-3003 then
  translated that invoice to the spread exponent
  `kappa=o(d^(1-3/(k+1)))` and called `2/5` the fourth-edge exponent.
- **Minimal witness / first failed implication:** `J_3` is part of the leading
  curvature `3J_2^2-2J_3-d^(-2)`, so it is not a remainder jet. At the exact
  box point
  `(d,x,z,w)=(701,129/100,39/20,-149/20)`, the advertised decimal condition
  applies, but the exact third-edge numerator is
  `-114191274399994230172453/10000000000<0`. Downstream, for edge `k=4` the
  old spread claim permits `kappa=d^(3/10)`, yet it does not imply the required
  `q_4=o(d)`. This is realized inside the positive-coefficient universe by
  `N_n(t)=((t+1)^2+n^6)^(n^10)`: here `d=2n^10`,
  `kappa=sqrt(1+n^6)=o(d^(2/5))`, but
  `q_4=1-6n^6+n^12` and `q_4/d~n^2/2`.
- **Exact repair / strongest survivor:** the graded remainder condition starts
  at `j=4`. Uniformly on curvature at least `923/10000`, the exact third-edge
  formula gives
  `G_3/d^6=C+6w/d+O((1+|w|/d)/d)`, so the safe boundary is the strict
  `liminf w/d>-923/60000`; in particular `w=o(d)` suffices. A labelled-polymer
  cluster expansion proves the all-order degree and first-occurrence law, so
  `q_j=o(d^(j-3))` for `4<=j<=k+1` is correct. Under
  `|q_j|<=kappa^j`, the binding exponent is always `j=4`, hence the single
  sufficient condition is `kappa=o(d^(1/4))` for **every fixed edge**.
  THM-3003 sharpens this with the cancellation tax
  `chi=mean|r|/|mean r|`: bounded `chi` improves `1/4` to `1/3`.
- **Rule:** separate jets already present in the leading invariant from true
  remainder jets. Never turn a leading-order asymptotic threshold into a
  finite non-strict bound without an exact monotone remainder estimate. When
  collapsing a family of exponent inequalities to one spread exponent, take
  the minimum over the actual jet range rather than substituting only its
  largest index.

## MISTAKE-338 (2026-07-30, pre-promotion THM-2980 audit) -- a positive-ray cutoff omitted zero and negative suffix rays

- **What was assumed:** the first THM-2980 candidate enumerated only suffix
  rays with positive scalar numerator. On three exceptional rows it split the
  remaining universe into four positive low labels or one positive high label
  plus three positive lows, then claimed this exhausted every packet.
- **Minimal witnesses / first failed implication:** on
  `E=(1,8,10,12,13,14)`, the packet
  `(1612,1836,2004,2340,20384)` is scalar-admissible and the last numerator is
  zero. On `E=(1,10,11,12,13,14)`, the same happens for
  `(1612,1736,1800,2340,210210)`. A negative ray `A/z` increases toward zero
  as height increases, so positive-ray monotone truncation reverses direction:
  if three positive lows have strict surplus, every negative residue ray is
  eventually admissible. The candidate also called the `z_1=1650` row
  non-finite although its positive cutoff is exact and positive.
- **Exact repair / strongest survivor:** the positive finite, zero-high, and
  one-high censuses survive unchanged. Exact gaps show that a packet with one
  nonpositive suffix must have exactly three positive companions, all below
  the high floor. There are only `18+9` such low triples on the first two
  rows and none on the third. The repaired verifier keeps every nonpositive
  residue/unit. For first representatives `r<L`, the correct translated-band
  capacity `kappa(d)=ceil(d/7)` closes all `3,861` fixed-carrier denominator
  pairs; the `5,159,799` exact unit instances independently regress that
  quotient. The smaller centered capacity is forbidden here by MISTAKE-334.
  If `z>=L` and that carrier contains one complete cell, `a=z/L>=1` supplies
  at least `floor(a)` full phase turns; the punctured safe image has mass at
  least `(6/7)floor(a)/a>=3/7>25/91`, so every later height and residue zero
  close.
- **Rule:** before compactifying a rational ray, audit the sign. Positive
  `A/z` decreases and admits an upper cutoff; negative `A/z` increases and
  generally creates an infinite tail. A finite residue quotient must retain
  zero and negative rays whenever the other slots already have strict scalar
  surplus.
## MISTAKE-337 (2026-07-31, THM-3001 section 6 classifier census) -- a 42/42 census held the failing axis fixed

- **What was done:** THM-3001 section 6 proposed that the two end curvatures
  `(sign C(mu), sign C(mu*))` classify the global shape of the Newton-ratio
  sequence, and supported it with an exact-rational census reporting `42/42`
  agreement across two-cluster, three-cluster and geometric families.
- **Why it was wrong:** the census varied the root ratios, the cluster count and
  the degree, but every three-cluster row used **equal** cluster sizes
  (`d//3` each).  Unequal multiplicities are exactly where the classifier
  breaks.  Restricted to equal sizes the failure rate is `0/30`; over all
  three-cluster configurations with `d=6..12` it is `51/2100 = 2.43%`.
- **Minimal witness (THM-3004):** `N(n)=(n+1)^2(n+3)^2(n+8)`, degree `5`,
  `R=(256/215, 1849/1600, 10000/8643, 4489/4000)` -- down, up, down, so two
  circuit sign changes, while both end curvatures are positive.  All roots real
  and positive, so the witness is PF-infinity, Hurwitz and strictly ULC: it is
  interior to every class in the lane, not a boundary artefact.
- **Exact repair / strongest survivor:** the classifier is true exactly for at
  most two clusters (exhaustive, `936` configurations, zero violations).  The
  correct general law is a **cluster count**: `m` well-separated clusters give up
  to `2m-3` sign changes, attained for every `m<=6`.  THM-3001's proved
  quantitative necessary condition
  `C(mu_d)>=-O(1/d), C(mu_d*)<=O(1/d)` survives (MISTAKE-341); it is not
  sufficient, and no bounded set of moments can be, since the sign-change count
  is a property of the support structure.
- **AMENDMENT (same day, klein-S428, after an independent adversarial pass).**
  The diagnosis above is the *dominant* mechanism but not the only one, and one
  circulating description of it is too strong.  Exact split over the `51` failing
  three-cluster configurations: the census's own `shape_of` classifier would have
  **disagreed with the prediction on 46 of them** (it returns `MIXED`), and only
  `5` slipped through.  So the census was **not** vacuous, contrary to an
  external audit note claiming its agreement was "logically compatible with
  arbitrary interior oscillation".  What is true is that `shape_of` decides
  `INTERIOR-MAX`/`INTERIOR-MIN` from `R[2]>R[1]` and `R[d-1]<R[d-2]` alone -- the
  two **end** circuits -- so it has a genuine blind spot on exactly the W-shaped
  palindromes (`c_2>0`, `c_(d-1)<0`, oscillating interior), e.g.
  `(n+1)^2(n+3)^2(n+9)^2`.  Two independent defects, therefore: the pinned
  cluster-size axis (which meant no failing configuration was ever presented) and
  the end-only shape branch (which would have mislabelled `5/51` if one had
  been).  Second rule below is amended accordingly.
- **Rule:** a census is evidence only about the coordinates it actually varies.
  Before quoting an `n/n` census, enumerate the coordinates of the configuration
  space and mark which were moved and which were pinned; report the pinned ones
  next to the score.  Sample size never substitutes for an un-varied axis.
  Second rule: when a mechanism is available (here the multipole/step-function
  picture of THM-3003), derive the predicted failure mode and search for it
  directly instead of sampling.

## MISTAKE-336 (2026-07-31, merge `f737bbe22922`) -- unresolved conflict markers were committed into PROVED canon

- **What was done:** the merge commit `f737bbe22922` ("Merge branch 'main'")
  was committed with **unresolved** conflict markers still in five tracked
  files: `THM-2596`, `THM-2597`, `THM-2598` and the companions
  `04-computation/jacobian_quartic_v4_resolvent_thm2598.py`,
  `04-computation/modular_farey_gram_owner_cocycle_thm2596.py`
  (`35` conflict blocks total).  It was then carried forward through
  `e4ee2e93710a` onto `origin/main`.
- **Why it was wrong:** the three theorem files are PROVED canon.  Every reader
  and every downstream citation between the merge and the repair saw
  interleaved `<<<<<<<`/`=======`/`>>>>>>>` text with two contradictory status
  lines in the same frontmatter, and `agents/check_docs.py` failed globally
  ("tracked files contain merge-conflict markers").  A conflicted PROVED file
  is worse than a `RESERVED` stub: it *looks* citable.
- **Minimal witness:** `git grep -lE "^<<<<<<< "` on `origin/main` returned five
  files; `THM-2598` carried `13` blocks including two different `status:`
  values (`PROOF-COMPLETE CANDIDATE ... AWAITING INDEPENDENT HOSTILE AUDIT`
  versus `PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`).
- **Exact repair (klein-S428):** the two merge parents were both clean.  Parent
  `949488f90070` holds the later, audited promotions of exactly these
  theorems; parent `db42ecbe246e` holds their pre-audit candidate copies, and
  every line unique to it is superseded status text or prose that the audit
  rewrote -- no unique mathematics.  The five paths were restored to
  `949488f90070`, `check_docs.py` was re-run green, and each restored
  `.md`/`.py` pair was re-hashed against its recorded `script_sha256`
  (both match).  No other path in the merge was touched.
- **Rule:** never commit a merge without `git grep -lE "^<<<<<<< "` returning
  empty, and run `agents/check_docs.py` **before** committing, not after.  When
  repairing someone else's conflicted canon, do not adjudicate mathematics:
  check whether one parent is the audited successor of the other, verify no
  unique content is lost, restore that parent verbatim, re-verify the recorded
  hashes, and report the exact commands.

## MISTAKE-335 (2026-07-30, first THM-2991 promotion) -- a late directional turn was called a global leading-edge return

- **What was done:** the first promoted THM-2991 used
  `(x+1)^n(x+B)^n`.  Its exact reciprocal symmetry gives an arbitrarily long
  ladder `R_1<...<R_n` followed by `R_(n+1)=R_(n-1)<R_n`.  The theorem then
  said this disproved global no-return from the leading edge.
- **First failed implication / minimal witness:** `R_(n+1)<R_n` only says the
  ratio path changes direction.  For `n>2`, the same ladder and symmetry give
  `R_(n+1)=R_(n-1)>R_1`; the later circuit has not returned below the leading
  edge.  The proof and exact transcript were correct, but the claimed
  consequence conflated directional monotonicity with the stronger property
  `R_j>=R_1` for all `j`.
- **Exact repair / strongest survivor:** the two-cluster family remains a
  sharp directional-turn control.  Repaired THM-2991 uses
  `(x+1)^n(x+3)(x+C)^n`: as `C` grows, `R_1<...<R_n`, while reciprocal-edge
  asymptotics give `R_(2n)<R_1`.  This is the required arbitrarily delayed
  global return inside PF-infinity, Hurwitz stability, and strict ULC.
- **Rule:** define a path obstruction by the actual comparison relation.
  A late negative discrete derivative does not imply return below an earlier
  baseline; freeze both a directional witness and a baseline-crossing witness.

## MISTAKE-334 (2026-07-30, projected k3 unit-count simplification) -- a centered bad-band count was applied after translation

- **What was done:** a proposed simplification of the `z_1=270..247`
  terminal closure used THM-2984's centered bad-band cardinality
  `beta(d)=2 floor((d-1)/14)+1` after the three aligned coordinates had left
  an arbitrary local phase translation.  The proposal would have replaced a
  pair certificate by the test `|S|>beta(d)`.
- **Minimal witness / first failed implication:** for `d=28`, `beta(28)=3`,
  but the four residues `S={0,1,2,3}` all lie in the translated open circular
  interval `(-1/2,7/2)` of length `d/7=4`.  Rotation preserves arc length but
  not which strict endpoints meet the residue lattice.  The first invalid
  step was silently recentering a translated open band at the integer zero.
- **Exact repair / strongest survivor:** an arbitrary translated danger band
  has the different sharp lattice capacity
  `kappa(d)=ceil(d/7)`.  THM-2984 now proves this separately and preserves
  `beta(d)` only for the absolute centered phase.  The projected terminal
  simplification still works because its stronger deterministic count is
  `|S|>alpha`, where `alpha=d/R`, `R|d`, and `R<=7`; hence
  `alpha>=ceil(d/7)`.  No pair-selection policy is required, but the replay
  must track `kappa`, not `beta`.
- **Rule:** every circular-band count must declare its center/translation
  sidecar.  Centered strict-open endpoint savings do not survive arbitrary
  rotation; use `beta(d)` only at an absolute phase and `ceil(d/7)` when the
  local phase is translated.

## MISTAKE-333 (2026-07-30, THM-2941 `z_1=297` torsion descent) -- the repaired LP-dual digest bug recurred in a new descendant

- **What was done:** the first `z_1=297` located-torsion verifier exactly
  checked all `830` Farkas certificates, but hashed each complete status
  witness, including the arbitrary solver-selected dual basis.  It also put
  those complete witnesses and the basis-dependent minimum contradiction in
  its semantic payload.  The first `z_1=294..276` descendant imported that
  payload wholesale.
- **First failed implication:** exact validity of a returned Farkas
  certificate was again treated as uniqueness of its representation.  The
  closure counts `1172=271+830+71`, the `73` located torsion cases, and the
  later `1549=659+882+8` descent are canonical; the chosen HiGHS dual vectors
  are not.  Thus the advertised byte/semantic replay could change under an
  equally valid optimizer basis without changing any infeasible instance or
  proof conclusion.  This is the precise mechanism already recorded in
  MISTAKE-331.
- **Exact repair / strongest survivor:** every solver certificate is still
  rebuilt and verified over exact rationals.  Reproducibility records now bind
  only `(denominators,q,M,marginals,capacity-set,load-histogram)` and the number
  of verified instances; the certificate itself and its basis-dependent
  contradiction value are excluded.  The repaired `z_1=297` artifact is then
  repinned through the lower descent.  The later `z_1=275..272` septimal
  descendant likewise imports only this certificate-free instance ledger;
  all `1,844` duals are still checked exactly.  All ray, high-wall,
  finite-low-pair, torsion-density, cap, and ledger conclusions survive
  unchanged.  A post-promotion audit also found that this descendant declared
  an LF-normalized dependency-hash basis while its helper hashed raw checkout
  bytes.  The helper now normalizes CRLF and bare CR to LF before hashing;
  the stored output, profile, semantic digest, cap, and ledger are unchanged.
- **Rule:** a correction to a proof-search digest is a dependency-graph
  invariant, not a one-file patch.  New descendants of a nondeterministic
  solver must reuse the canonical instance ledger explicitly; never copy the
  old witness tuple merely because its exact checker passes.

## MISTAKE-332 (2026-07-30, first `z_1=1736` hybrid replay) -- a repaired status rule was not inherited by a duplicated solver

- **What was done:** commit `a90c87e293cd` independently duplicated the K5
  common-status solver for the fifteen-row `z_1=1736` closure after
  MISTAKE-331, but its stage digests again hashed the floating solver's full
  `(alpha,z)` Farkas representative.  The resulting closure is mathematically
  valid, yet its advertised byte-level evidence is noncanonical.
- **Minimal witness / first failed implication:** on
  `E=(1,2,8,10,12,14)`, all counts remain `46 -> 0+16+30`, while the raw-dual
  status digest `9e7a9bb7...46732` becomes the canonical infeasible-instance
  digest `f63ad05b...7f4b`.  Both runs exactly verify their returned rational
  certificates; equality of a solver-selected basis is the first invalid
  reproducibility demand.
- **Exact repair / strongest survivor:** retain exact rational verification of
  every dual, but hash each deterministic status-kill row with
  `witness[:-1]`, omitting only the noncanonical certificate representative.
  The repaired all-label replay preserves the complete fifteen-row closure
  and strengthens the exceptional body to all `749` unrestricted packets,
  with minimum projected margin `121/18109`.
- **Rule:** copying a proof engine also copies its repaired evidence boundary.
  Every duplicate status solver must be audited against MISTAKE-331 before it
  can become a canonical replay dependency.

## MISTAKE-331 (2026-07-30, THM-2941 ray/status replay) -- a semantic digest bound a noncanonical LP dual basis or magnitude

- **What was done:** the three THM-2941 residue-ray/status companions verified
  every rational Farkas certificate exactly, but then hashed the literal
  certificate representatives returned after a floating HiGHS search.  The
  theorem consequently advertised byte-identical replay through a quantity
  that was not mathematically canonical.  The first repair removed the raw
  certificates but still retained a minimum exact-contradiction field in the
  `z_1=324`, `z_1=328`, and `z_1=312` records, semantic hashes, and outputs.
  That derived magnitude is equally basis- and scale-dependent.
- **Minimal witness / first failed implication:** an unchanged replay of the
  `z_1=250` companion preserved all `35,990 -> 1,965 -> 1,266 -> 0` counts,
  the decisive-modulus histogram, and the displayed first exact witness, but
  HiGHS selected another valid dual basis.  The raw certificate digest changed
  from `f017d330...046611` to `bcb55c09...324bfc`, so the frozen semantic hash
  failed even though every returned certificate passed the exact inequalities.
- **Exact repair / strongest survivor:** exact verification of every returned
  certificate remains load-bearing.  Replay hashes now bind the deterministic
  infeasible instances, counts, first witness, exact-check count, and the
  basis-invariant fact that every replayed contradiction is negative, never
  the solver-selected dual representative or its magnitude.  This final
  magnitude repair covers `z_1=324`, `z_1=328`, and `z_1=312`; the same narrow
  basis repair is applied to `z_1=378`, downstream `z_1=306/302/298`, and
  `k=2` companions because they shared the latent mechanism.  None of the ray,
  capacity, Hunter, Farkas, cap, or ledger conclusions changes.  The later
  `lrc14_j7_k2_z1736_hybrid_closure_thm2941.py` and
  `lrc14_j7_k2_z1736_exact_ray_status_projected_closure_thm2941.py`
  status-stage replays regressed to hashing raw HiGHS certificates; their
  mathematical checks may survive, but their transcript, profile, and
  semantic pins are noncanonical until the same repair is applied or a
  deterministic closure supersedes them.
- **K2 follow-through audit:** the later canonical descent supersedes those two
  `z_1=1736/1732` replays.  A targeted duplicate-solver scan also repaired the
  independent `z_1=1784` referee and the redundant `z_1=1836 -> 1824 -> 1812`
  side chain.  Those programs still reconstruct and verify every rational
  certificate, but their stage, profile, semantic, and displayed replay data
  now exclude the certificate basis and every magnitude derived from its
  arbitrary normalization.  The cap-`1656` lineage uses this canonical
  boundary throughout.
- **Projected-k3 follow-through audit (2026-08-02):** THM-3078 appended an
  exact-audit digest as worker field `row[21]`, but that digest itself hashed
  the solver-selected certificate and contradiction magnitude.  THM-3102,
  THM-3106, THM-3109, THM-3111, and the first THM-3113 candidate then hashed
  or imported the complete worker rows.  Their state, Farkas, terminal,
  carrier, ledger, and cap mathematics survives.  The repaired descendant
  chain still rebuilds and exactly checks every returned dual, but screen and
  semantic digests bind only `tuple(row[:19])` plus the basis-invariant
  direct/legacy branch counts.  THM-3078's historical semantic hash remains a
  noncanonical evidence artifact and is never inherited by the repaired
  descendants.
- **Rule:** a proof search may use a noncanonical optimizer witness, but a
  reproducibility digest must bind canonical problem data or a separately
  canonicalized certificate, not an arbitrary valid basis or a scalar derived
  from its normalization.

## MISTAKE-330 (2026-07-30, chained k=2 frontier verifiers) -- descendant evidence pinned superseded atlas bytes

- **What was done:** the standalone `z_1=1824` verifier was produced against
  the first `1810..1835` scalar-atlas transcript.  The atlas was then repaired
  to state the scalar-empty `z_1=1811` slice explicitly, but the verifier kept
  the old source/output hashes and old exact-line sentinel.  The later
  `z_1=1812` verifier inherited those stale atlas pins and also pinned the
  pre-repair `z_1=1824` source/output pair.
- **First failed implication:** mathematical equality of the surviving
  `z_1=1824` and `z_1=1812` rows was treated as sufficient dependency
  stability.  Replaying either descendant at the current canonical tip
  stopped at its hash gate before reaching the unchanged exact quotient.
  Thus these were reproducibility failures, not counterexamples to either
  row closure.
- **Exact repair / strongest survivor:** the `z_1=1824` verifier now pins the
  repaired atlas pair and full empty-slice sentinel, with regenerated profile
  and semantic digests; the `z_1=1812` verifier then pins that repaired parent
  and the current atlas pair.  Ordinary and optimized replays are
  byte-identical.  Their `38=15+23` and `11=4+7` crude/status closures survive,
  and the independent combined ten-row descent already proves the stronger
  current cap `z_1<=1799` from current atlas bytes.
- **Rule:** after repairing a truth-bearing artifact, traverse its descendant
  hash graph in topological order and replay every pinned verifier at the
  fetched tip.  An unchanged mathematical row does not make obsolete evidence
  pins current.
## MISTAKE-329 (2026-07-29, first mixed Lorenz/activity sidecar) -- an ambient half-cap was substituted for the exact reflected parity fibre

- **What was done:** the first four-aligned/three-drift Lorenz/activity
  sidecar handled repeated denominators `(2,2,d_3)` with the coarse test
  `|S_D|>2C_3`, where `C_3` is the third needle's ambient capacity.  It
  reported a residual of `29,221`.
- **Minimal witnesses / first failed implication:** reflection
  `r -> D-1-r` makes the two parity loads of `S_D` exactly equal, but the
  third needle's capacity inside one parity is the sharp `q=2` fibre value,
  not half of its full ambient capacity.  The correct kill is
  `|S_D|/2>fibre_cap(D,d_3,2)`.  The coarse replacement misses exactly
  two rows:

  ```text
  F=(1,4,5,7,9,11),  D=194040, S=55392,
  (d_1,d_2,d_3)=(2,2,194040),  27696>13860;

  F=(1,5,7,8,9,11),  D=388080, S=109044,
  (d_1,d_2,d_3)=(2,2,388080),  54522>27720.
  ```

- **Exact repair / strongest survivor:** the sidecar now reconstructs the
  reflected parity histogram and applies the exact `q=2` fibre cap.  Its
  corrected chain is `544,571 -> 419,364 -> 29,219`; ordinary and optimized
  replays agree.  The former scratch `29,221` TSV/digest is superseded.
  The load-bearing closure
  `544,571 -> 419,511 -> 29,364 -> 19 -> 0` never depended on that TSV and
  is unchanged.
- **Rule:** a symmetry-balanced target load must be compared with the exact
  capacity in the corresponding symmetry fibre.  A full ambient capacity,
  even when divided by the orbit size, can lose gcd/remainder information
  and is only a relaxation.

## MISTAKE-328 (2026-07-29, THM-2908 infinity replay) -- a SymPy structural comparison was used as a polynomial identity gate

- **What was done:** after repairing the removable `n+3` content in the first
  projective-infinity invariant, the THM-2908 companion compared the reduced
  second numerator to its displayed factorization with Python/SymPy `==`.
- **First failed implication / minimal witness:** the exact reconstruction
  reduces the second invariant to
  `-2(n+2)(4n+5)/(3(3n+4)(3n+5))`, so its numerator differs from the stated
  `-2(n+2)(4n+5)` by the zero polynomial.  Nevertheless the long ordinary and
  optimized replays both stopped at that structural comparison.  SymPy `==`
  is not a general algebraic-identity decision procedure across differently
  normalized expression trees.
- **Exact repair / strongest survivor:** every infinity comparison is now
  made after cancelling the difference and testing it against zero.  A
  separate direct reconstruction confirms both displayed reduced invariants;
  the degree-`2804` finite factorization and all preceding certificates were
  unchanged.
- **Rule:** truth-bearing symbolic gates compare a normalized difference with
  zero (or use an explicit polynomial identity), never raw expression-tree
  equality after a long exact computation.

## MISTAKE-327 (2026-07-29, THM-594 pair overlap) -- a one-wrap Fourier branch was used as an all-pair formula

- **What was claimed:** THM-594(B) asserted that for every coprime
  `p<q` with `(p+q)/14>1`,
  `mu(D_p intersect D_q)=(q+2p-14)/(7pq)` at danger radius `1/14`.
- **Minimal witness / first failed implication:** for `(p,q)=(50,51)`, that
  expression is `137/17850`, while the exact finite trapezoid sum of the
  later LEM-042 is `73/3570`; their difference is `38/2975`.  The
  product-to-sum quadratic was continued after the first wrap without
  retaining subsequent periodic wraps, so the purported second branch was
  not global.
- **Strongest survivor:** THM-594(A)'s Fourier support, THM-594(C)'s
  divisor-minimal no-exact-tiling argument, and the Parseval floor in
  THM-594(D)--(E) do not use the false branch and survive.  LEM-042
  supersedes all pair-mass calculations with
  `g sum_j min((W-|j| delta)_+,w_min)` and records both below- and
  above-`1/49` families.
- **Rule:** pair overlap past the first Farey cell is a periodic trapezoid
  sum, not a two-cell formula.  Preserve gcd, reduced residues, and all
  wraps; for inherited carriers also retain the positioned-window
  discrepancy, since full-circle mass does not localize.

## MISTAKE-326 (2026-07-29, concurrent THM-2926 reservations) -- a successful push does not prevent a later session from reserving the same namespace

- **What happened:** the seven-wall session reserved `THM-2926` at
  `9578144779f6`.  Two concurrently based GMC sessions then used the same
  identifier; one was repaired to `THM-2927`, while the consecutive
  four-slot stub remained `THM-2926`.
- **First failed implication:** a pushed reservation was treated as globally
  exclusive without every later session rechecking the fetched remote tip.
  The result was two filenames and two YAML records with `id: THM-2926`.
- **Exact repair:** before either stub became a proved dependency, the
  seven-wall namespace was moved to the freshly fetched and independently
  checked `THM-2928`.  The mathematical statement and provenance are
  unchanged; no theorem ever depended on the ambiguous ID.
- **Rule:** immediately before *each* reservation commit, fetch and check the
  remote filenames, YAML IDs, and recent history.  After the next pull, scan
  the claimed ID again: reservation is a protocol, not a lock.

## MISTAKE-325 (2026-07-29, THM-2908 projective-infinity replay) -- `together` was mistaken for reduced rational form

- **What was done:** the candidate THM-2908 companion compared the numerator
  returned by `together(I_1).as_numer_denom()` at the missing infinity point
  with the reduced cubic invariant stated in the theorem.
- **First failed implication / minimal witness:** `together` combines a rational
  expression over one denominator but need not cancel its polynomial gcd.  In
  this cell its raw first numerator is
  `-(n+2)^2(n+3)(28n^2+87n+66)` and the raw denominator also contains
  `n+3`.  The candidate therefore aborted after its long resultant replay even
  though the reduced invariant is exactly the theorem's
  `-(n+2)^2(28n^2+87n+66)`.  The second invariant has no such removable
  content.
- **Exact repair / strongest survivor:** apply `cancel` before extracting the
  numerator, and retain an explicit hostile gate proving that the raw
  numerator and denominator each differ from the reduced pair by exactly
  `n+3`.  The projective-infinity nonvanishing statement and every finite-chart
  certificate are unchanged.
- **Rule:** when a proof identifies a rational invariant rather than a chosen
  clearing, normalize by gcd before comparing numerators.  If removable
  content is mathematically relevant to replay portability, freeze both the
  raw and reduced pairs explicitly.

## MISTAKE-324 (2026-07-29, THM-2920 depth-two working audit) -- a pair statistic was mistaken for a two-slot residual

- **What was inferred:** after an ordered five-cover branch fixed centres
  `x` and `y`, a preliminary audit compared the twice-subtracted carrier
  with its global top two singleton coverages and reported the branch
  closed.
- **First failed implication / minimal witness:** fixing `x` leaves four
  labels and fixing `y` leaves **three**, not two.  The fact that the next
  certificate computes a pair cap does not consume a second label.  On the
  five live THM-2920 tips, the top-two inequalities are true but insufficient:
  a third allowed danger comb is still available.
- **Exact repair / strongest survivor:** retain the correct three-slot type.
  The exact global pair cap combined with the leading singleton closes three
  of the five tips.  The lawful three-slot Hunter envelope closes all five,
  with strict margins
  `4211/280280`, `2459243/348107760`, `52873/2802800`,
  `143831/12892880`, and `7009/420420`.  Thus the intended two-H3-row
  closure survives, but not by the discarded top-two argument.
- **Rule:** every recursive toothpick node carries an explicit remaining-slot
  counter.  Subtracting one chosen centre decrements it by exactly one;
  computing a pair cap, matching bound, or tree credit changes the bound,
  not the arity.  Before promoting a terminal inequality, reconstruct the
  family size `fixed labels + remaining slots` directly.

## MISTAKE-323 (2026-07-29, first THM-2921 scratch certificate) -- a divided-power coefficient table was convolved as an ordinary polynomial

- **What was done:** the first diameter-four Macaulay probe stored one
  symmetric order-`m` tensor entry for each exponent triple and inserted
  those values directly as coefficients of ordinary monomials.  The
  resulting fixed determinant had an attractive degree-`196`
  Gregory--Newton certificate.
- **First failed implication / minimal witness:** a symmetric tensor table
  is a divided-power coefficient table.  Ordinary multiplication requires
  the multinomial multiplicity `m!/(alpha_0!alpha_1!alpha_2!)`; equivalently,
  divided powers multiply with binomial structure constants.  Already at
  depth zero on offsets `(0,1,2,4)`, the correct mixed quadratic coefficient
  is twice the stored tensor entry.  A direct original-form determinant
  therefore failed the claimed scaling identity.  The entire first minor
  and its apparent factorization were invalid and were retracted before any
  theorem reservation or promotion.
- **Exact repair / strongest survivor:** THM-2921 multiplies every tensor
  entry by its multinomial count, proves the common denominators by exact
  division in `Z[n]`, and recomputes the minor.  A genuinely separate direct
  four-variable multinomial expansion reproduces the repaired forms and
  selected minors modulo `1000003` at `21` family/depth controls.  The
  fixed-chart plus Gregory--Newton strategy survives and closes all three
  nonconsecutive diameter-four families, but none of the false determinant's
  coefficients or factors is retained.
- **Rule:** declare the polynomial basis before turning a symmetric tensor
  into a coefficient dictionary.  For an ordinary monomial basis, audit one
  mixed coefficient against an ordered-word or direct multinomial expansion
  before computing any resultant, determinant, or positivity certificate.

## MISTAKE-322 (2026-07-29, first THM-2901 promotion) -- the full-arity same-cap deadlock was applied to a lower-arity pair flag

- **What was written:** the first promotion of THM-2901 said that reusing a
  parent cap after selecting an `H_3` pair `L` was “precisely the same-cap
  deadlock” of THM-2893.
- **First failed implication:** THM-2893's algebraic deadlock assumes
  `ell=s`, so the selected flag itself is heavy and
  `U_C(L)+B>=h`.  Here `(s,ell)=(3,2)`: the pair need not be heavy, and an
  inherited parent-carrier certificate can succeed.  On THM-2900's exact
  control it succeeds on `204/784` pair flags, which directly refutes the
  claimed impossibility.
- **Exact repair / stronger survivor:** every extension
  `L union {w}` inside a hypothetical cover is a heavy triple.  On the
  literal residual `C_L`, this gives the exact link condition
  `c_(C_L)(w)>=h_L-B_2` for **each** of the three remaining labels.
  THM-2893 now records the general heavy-link recursion, and THM-2901 routes
  its `12,919` pair branches as link-constrained children.  Parent caps are
  optional sufficient tests; literal child caps retain the link and
  forbidden-prefix sidecars.
- **Rule:** before invoking a no-go theorem, audit its arity equality.  A
  lower-arity selected flag retains nontrivial links into the child; it is
  not a full-arity heavy flag and need not inherit the same obstruction.

## MISTAKE-321 (2026-07-29, concurrent THM-2889 reservation) -- a local filename scan was not replayed after the reservation fetch

- **What happened:** the GMC low-sector session found no `THM-2889` in its
  freshly fetched tree and reserved that identifier.  The immediately
  preceding remote commit had already reserved `THM-2889` for the dicyclic
  reverse-action carrier, but the filename check was run before that fetched
  commit was made visible in the worktree.
- **First failed implication:** “absent from the checked tree and prior
  history” was treated as “still absent after the concurrent fetch/rebase.”
  The dicyclic reservation at `b01066a1c736` precedes the GMC reservation at
  `0a614c4e623c`, so chronology gives `THM-2889` to the dicyclic theorem.
- **Exact repair:** the GMC stub is moved to the independently rechecked free
  namespace `THM-2890`; no proved dependency ever used either empty stub.
- **Rule:** after the fetch that immediately precedes a reservation commit,
  rerun the filename, YAML-ID, index, and remote-history checks against the
  fetched remote tip—not only against the pre-fetch worktree.

## MISTAKE-320 (2026-07-29, first THM-2878 promotion) -- a zero-address census is not global uniqueness

- **What was claimed:** the first promoted wording described `(guard,u1)` as
  the unique full factor square and `(u3,D->S)` as the unique one of the
  eighteen oriented factor toggles without retaining the zero-address
  qualifier.  The companion had exhausted only the distinguished
  `ell=0` word orbit for those two uniqueness claims.
- **Minimal witnesses / first failed implication:** at canonical address
  `(1,0)`, the unique full pair is `(u1,u5)`, not `(guard,u1)`.  At address
  `(7,0)`, both `(u1,D->S)` and `(u3,D->S)` count
  `floor((q+h)/13)` on every one of the `169` positive edges.  Therefore a
  section-level exhaustive census cannot be promoted by silently
  quantifying over the quotient address.
- **Exact repair:** THM-2878 now audits all `169` canonical addresses on both
  source and target word orbits.  Full-pair counts have histogram
  `{0:54,1:64,2:38,3:10,4:3}`, with `38` distinct full-pair sets.  The
  marker census has `121` addresses with only `u3`, `44` with `u3` plus
  one shifted `u1/u2` event, and `4` with `u3` plus two.  Thus `(guard,u1)`
  is uniquely full at zero address, while `(u3,D->S)` is the unique
  **address-uniform** carry marker: it works at all `169` addresses, whereas
  each extra orientation works at only `13`.
- **Strongest survivor:** the carry theorem, all `169` edge identities, all
  `2,197` composition identities, and the `omega^3` seam survive unchanged.
  The full word supplies the carry transition but still not its initial
  ancestry state or a physical `QA->QAB` current.
- **Rule:** whenever a computation selects a distinguished section, audit
  the claimed property on every quotient fibre/address before asserting
  global uniqueness.  Record separately “unique in this chart,” “pointwise
  unique,” and “unique object present uniformly across charts.”

## MISTAKE-319 (2026-07-28, provisional THM-2846 geometry) -- cubic divisibility cancellation is not orientation cancellation

- **What was claimed:** the provisional positive-cone hostile said that two
  opposed transport contributions “cancel exactly” on the locus
  `I1=I2=0`, suggesting that the THM-2830 orientation `D(U,V)` vanished
  there.
- **First failed implication:** THM-2824's exact factorization is
  `I2=-g11 D-t111 g22^2`.  Hence `I2=0` forces
  `D=-t111 g22^2/g11<0`, because the Gram entries are positive and
  THM-2841 gives `t111=L(U^3)>0`.  For the displayed witness the numerical
  value is about `-4101.29`, emphatically not zero.
- **Strongest survivor / repair:** the construction itself survives.  The
  two **quadratic-division remainders** vanish, so the positive plane really
  contains a nonzero factorial moment-one/two/three null line.  Its
  orientation is strictly negative, its fourth moment is nonzero, and its
  first variance jet survives.  THM-2846 now states this corrected
  transverse geometry and the exact companion certifies the complete
  resultant factor and the algebraic point's rational rectangle.
- **Rule:** when several derived quantities encode the same cubic tensor,
  name the exact vanishing object.  Cancellation of Euclidean-division
  remainders does not imply cancellation of an oriented determinant;
  substitute the defining relation before giving geometric prose.

## MISTAKE-318 (2026-07-28, THM-2834 first Picard-lattice proof) -- conjugate split curves are complementary, not numerically equal

- **What was claimed:** for the fourteen curves on
  `X_14 subset P(2,2,7,7)`, the first proof wrote `C_i ≡ D_i` and
  `C_i.C_j=delta_ij/2`, then treated Frobenius as the ordinary seven-cycle
  on the `C_i`.
- **First failed implication:** the weight-two form `ell_i` cuts
  `C_i union D_i`, so its divisor gives the unavoidable identity
  `C_i+D_i=2H`.  Thus `D_i=2H-C_i`, not `C_i`.  The old positive-definite
  matrix also violated the Hodge-index signature.  On the tame local covers,
  distinct `C_i,C_j` meet at the common `1/7(1,1)` point with intersection
  `1/7`, while `C_i,D_i` meet at an `A_1` point with intersection `1/2`.
  Therefore
  `C_i^2=-5/14`, `C_i.C_j=1/7` for `i!=j`, and `C_i.D_j=0` for `i!=j`.
- **Strongest survivor / repair:** the theorem's ranks and point counts
  survive.  The corrected `C_1,...,C_7` Gram determinant is `1/128`, so the
  seven classes are independent.  Since Frobenius interchanges the two
  factors, it acts by
  `C_i -> D_(i+1)=2H-C_(i+1)`, or by the **negative** seven-cycle
  `v_i=C_i-H -> -v_(i+1)` on `sum v_i=0`.  Its fixed space is exactly `Q.H`,
  proving `rho(X/F_3)=1`; its eigenvalues
  `{1} union {-zeta_7^j:1<=j<=6}` still give traces `2,0,2,0`, matching the
  direct counts.
- **Rule:** whenever a reducible coordinate divisor produces two conjugate
  curve families, write every principal-divisor relation before identifying
  numerical classes, and test the proposed Gram matrix against the
  Hodge-index signature.  A transitive orbit on geometric curves does not
  mean Frobenius acts by the corresponding permutation after quotienting by
  divisor relations.

## MISTAKE-317 (2026-07-28, provisional THM-2827 all-pole extrapolation) -- a unique Newton face need not miss the prescribed response valuation

- **What was assumed:** after the `nu=3` and `nu=4` local Faber obstructions,
  the same polar/pure/regular trichotomy was claimed provisionally for every
  `nu=ord(V)>=3`.  In the pure-`q` cone, exact controls verified that the
  `Phi/H/Psi` top face is unique and strictly below every retained lower row.
- **First failed implication / minimal resonance:** for `R=3k+2`, put
  `delta=nu-2a>0`.  The unique `H` face has valuation `-2k delta`, while
  `K=TH` prescribes `v(H)=nu+1-3a=1-a+delta`.  These are equal exactly when
  `a=1+(2k+1)delta`, equivalently
  `nu=2+(4k+3)delta`.  The first relevant ray is
  `R=8,k=2,nu=13,a=6,delta=1` (with `b>=7`).  Face uniqueness prevents
  coefficient cancellation but does not contradict an equal required
  valuation.
- **Repair / strongest survivor:** no false proved theorem was promoted.
  The pushed empty THM-2827 reservation was immediately narrowed and is being
  rebuilt as a nonresonance atlas: all `R!=2 mod 3` pole orders and all
  off-ray `R=3k+2` pole orders remain candidates for the uniform argument;
  `nu=4` is unconditionally off every ray.  Any resonant closure must use
  the leading coefficient in `A_src K=lambda M` or a later response, not
  pure-face valuation alone.

## MISTAKE-316 (2026-07-28, provisional THM-2784 passport wording) -- a distinguished degree-one cycle is not nontrivial inertia

- **What was written:** the balanced square-potential passport was described
  as having an involution over zero and one nontrivial cycle over the third
  value.
- **Minimal boundary witness / first failed implication:** the sharp map
  `F=4x/(x-1)`, `V=x(x-1)^3` has degree `N=1`, with `e=0` and `r=1`.
  Every permutation is the identity.  A displayed length-`r` distinguished
  cycle need not be nontrivial, and `2^e1^s` may likewise be identity.
- **Repair / survivor:** the product-one passport is correct with “involution
  or identity” over zero and “at most one nontrivial cycle” over the third
  distinguished value.  The Riemann--Hurwitz ledger, transitivity, and all
  divisor classifications survive; the independent `3185`-gate audit checks
  the repaired degree-one boundary explicitly.

## MISTAKE-315 (2026-07-28, THM-2781 exact companion) -- degree four does not certify non-fourth-power status

- **What was done:** the unreduced-exponent hostile `f=(1+z^2)^2`, displayed
  with exponent `2/4`, was correctly described as a square that is not a
  fourth power, but its original exact gate checked only `deg(f)=4`.
- **First failed implication:** a polynomial of degree four can be a fourth
  power of a linear polynomial.  The gate verified a compatible intermediate
  statistic, not the advertised consequence.
- **Repair / strongest survivor:** a constant-one fourth root would have to
  be `1+uz`; its linear coefficient forces `u=0`, while `f` has quadratic
  coefficient `2`.  The primary companion now tests the unique linear
  candidate explicitly.  An independent coefficient-root engine exhausts
  `3408` small polynomials and verifies the theorem, the repaired hostile,
  the all-degree sharp family, and a characteristic-two failure.  The theorem
  statement and example were correct; only the truth-bearing gate was weak.

## MISTAKE-314 (2026-07-28, concurrent THM-2759 reservation) -- a clean local scan was not replayed after the reservation rebase

- **What happened:** the exact-prefix even Faber flux-gcd theorem was first
  reserved as `THM-2759`.  Concurrent mainline work had already claimed that
  identifier for the split degree-six componentwise monicization closure.
  Replaying the reservation above the incoming commit temporarily left two
  distinct files with the same YAML theorem ID.
- **First failed implication:** a namespace check made before a concurrent
  reservation does not remain valid after rebase, even when the textual rebase
  is conflict-free.  Git detects path conflicts, not duplicated identifiers in
  different paths.
- **Repair / survivor:** the first-pushed split degree-six theorem keeps
  `THM-2759`; the exact-prefix gcd theorem and all its citations moved
  coherently to the freshly checked `THM-2760` namespace.  No mathematical
  statement changed.  After every reservation rebase, search both filenames
  and YAML IDs on the rebased tree before pushing; a clean merge is not a clean
  namespace audit.

## MISTAKE-313 (2026-07-28, THM-2749 Section 5 / provisional THM-2751) -- a clock-blind natural-sheet carrier was subtracted from a fully marked common section

- **What was assumed:** the legacy helper
  `lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py::restrict_e3_and_sheet`
  was treated as the fully marked natural `U_(0,3)` source/target sheet.  Its
  coefficients were compared with THM-2749's two-sided common section, leading
  to provisional THM-2751's symmetric weighted wings and normalized gains
  `12` (common), `2` (formal wing), and `7` (folded natural sheet).
- **Why it was wrong:** the helper inserts `E3` and the four shifted
  `q1/c2/q2/c3` gates but never intersects the source with the source-one clock
  comb `d_(1,ell)`.  THM-2749 constructs its common section through
  `source_present_section(...,source_clock=1,...)`.  The delayed prefix banks
  are identical, so the omitted clock factor is the first and only typing
  mismatch; exact replay of the legacy script merely certified the wrong
  carrier.  A hostile audit that checked the candidate only against its pinned
  legacy constructor reproduced the same error until the factor lists were
  compared directly.
- **Repair / survivor:** the former gain-`2` THM-2751 body was retracted before
  promotion.  The current file at that slug is a separately rebuilt
  `RESERVED PROOF-COMPLETE` fixed-clock candidate awaiting hostile audit.  With
  the actual clocked constructor, the natural source coefficient equals the common
  coefficient `C=339633525654239542165440`, so its physically nonempty left
  wing is coefficient-null.  At `t=3` the target coefficient is
  `345341652135823400016960`; the right-wing coefficient is
  `5708126481583857851520`.  After content/root normalization the source,
  target, and right-wing profiles are `9,8,4`, respectively, and the actual
  folded target/source ratio is `11`, not `7`.  The source-clocked one-sided
  bank extends gain `11` uniformly across all `81` canonical endpoint labels;
  it is not THM-2749's two-sided fibre product.  THM-2749's two-sided common
  clutch, all-rail raw equality, normalized sign `-1`, target window/unit,
  open cylinder, and primitive target-character profiles are unaffected.
  General rule: before subtracting two exact interval tables, compare the
  complete factor constructors, not merely their prefix banks and output
  hashes.
- **Independent retyping check:** if the two legacy outer sheets are instead
  treated honestly as full unclocked unions, their literal intersection is
  `M=disjoint-union_e M_e`, not `M_1`; the same-clock `e=2,3` pieces are
  nonzero and all cross-clock pieces are empty.  The resulting wings have
  `v_13(L)=1`, `v_13(R)=2`, with source profile `(0,0,0,12,0,0,0)` and target
  profile `(0,9,2,2,0,0,0)`.  Both physical-present profiles are `Phi_7`
  units: their determinants are `1` and `11`, and the exact coefficient-ring
  ratio is
  `g=2z+2z^2+2z^3+2z^4+6z^5`, with `det(g)=11` and `g s=p`.
  The target nevertheless augments to zero, so folding to the repeated
  delayed-clock scalar gives the zero profile and no augmented scalar gain.
  The ratio `g` is not a physical carrier map.  The clock-blind whole-sheet
  scalar ratio `7` belongs to the unclocked union, not to the fully marked
  fixed-clock carrier.  Thus the two valid repairs are distinct and neither
  supplies a physical `L -> R` map, relative phase, endpoint current, or LRC
  exclusion.
- **All-label sharpening:** the forgotten-`e` target augmentation is zero at
  the displayed `(s,t)=(0,3)` cell, not uniformly on the `81`-label bank.
  Exact fixed-`j` reconstruction gives a source unit precisely on a Cartesian
  `7 x 7` rectangle.  On its `49` cells the augmented gain census is
  `0:10, 2:19, 11:20`; the target physical-clock profile is a `Phi_7` unit on
  all `81` cells.  There are four coefficient-ring ratios, so neither the
  chosen determinant-`11` ratio nor the revived gain `2` is a global action.
  The strongest survivor is the exact rectangle/chamber law recorded in
  `lrc14-physical-wing-rectangle-and-target-chamber-20260728.md`.
- **Endpoint/address hostile:** on the fixed THM-2334 triangle, exhaustive
  search finds no scalar-times-character-times-quotient-shift covariance for
  either endpoint wing or even the two common carriers.  Transporting every
  present endpoint factor does restore an exact scalar, paired with its
  reciprocal in the reverse chart.  The reason is structural: left/right
  carrier harmonics `k,l` satisfy `r.W=l-k`, so exact address lives on the
  extended lattice `Lambda_tilde={(r,k,l):r.W=l-k}`; the old `Lambda` is only
  its diagonal `k=l`.  Interval collapse remembers merely `pi(r mod 13)`.
  A Bezout correction is noncanonical and changes the address.  Hence the
  factor-level harmonic origin (or a declared Bezout section) is an essential
  sidecar, and the formal gain is one moment rather than an addresswise map.

## MISTAKE-312 (2026-07-28, root-zero clutch dependency pins) -- LF evidence hashes were checked against raw platform bytes

- **What was recorded:** the finite-exact root-zero overlap companion pinned
  six audited dependencies by their LF SHA256 values but hashed each dependency's
  raw worktree bytes at runtime.  It passed on an LF checkout and failed before
  mathematics on a Windows CRLF checkout.
- **Why it was wrong:** a declared LF evidence address is portable only if the
  verifier applies the same normalization.  Replaying the main script alone
  had not exercised this cross-platform branch.
- **Repair / survivor:** dependency bytes are now normalized from CRLF to LF
  before hashing.  Normal and optimized Windows runs again byte-match the
  stored transcript; all six dependency pins and every finite-exact overlap
  count survive.  The script digest was refreshed after the repair.

## MISTAKE-310 (2026-07-28, relative-present root-zero clutch) -- a forbidden root label was mistaken for empty physical support

- **What was assumed:** the relative-present scout correctly found that every
  edge-preserving translation from the residue-7 right-root-12 chart lands in
  the right-root-0 chart.  Because root `0` is omitted from the private-root
  label bank, the accompanying reflection then declared the induced
  private-root lift relation empty and treated an unconstructed edge switch as
  the terminal obstruction.
- **Why it was wrong:** the private half-tooth labels form overlapping charts,
  not a partition.  In the canonical `182`-grid,
  `H_0^R=(0,13)/182` and `H_1^L=(1,14)/182`, so their open overlap is
  `(1,13)/182`.  Translation by `7/13^6` sends right root `12` to right root
  `0`, but a point in that overlap is simultaneously a lawful left-root-`1`
  point.  The exact adjacent witness
  `47850889647341/100360982066072 ->
  47851035194197/100360982066072` lies strictly inside the overlap and has
  common rail, clock, relative-present, semantic, and full-target support.
- **Repair / survivor:** landing at the zero **label** proves only that the
  full edge-preserving labelled row does not intertwine.  It does not prove
  physical emptiness after a chart change.  THM-2744 now proves the corrected
  partial Cech clutch: the restricted raw vectors agree and are units, while
  their root-normalized classes differ by `-1`.  The relative coefficient
  census, aggregate unit tables, endpoint counts, midpoint full-target
  coexistence, and empty unrestricted fixed-edge nonzero-root label relation
  survive.  Whole-cylinder target stability is proved only for the displayed
  rail-`8` pair; a global target action, endpoint current, row exclusion, and
  LRC(14) still do not follow.

## MISTAKE-311 (2026-07-28, THM-2727 evidence and provenance) -- a reserved candidate carried a stale script hash and cited an empty namespace for a proved formula

- **What was recorded:** reserved/provisional THM-2727 declared a script
  SHA-256 beginning `1e95f268` and called its private-root equation the
  “THM-2717 formula.”
- **Why it was wrong:** the tracked LF-normalized script hashes to
  `a6368190...`; normal and optimized executions still byte-match the stored
  transcript.  At the time of the bad attribution THM-2717 was still an
  `UNPROVED EMPTY STUB`; after its independent promotion, the formula
  `r=2c+floor(d/13)+1_(edge=left) mod13` remains sourced directly to the
  earlier proved THM-2640.  The same
  synthesis also blurred fourteen coefficient-only constant matches with the
  physical signature census, whose intersection is actually zero.
- **Repair / survivor:** THM-2727 now depends explicitly on THM-2640 and the
  subsequently proved THM-2717 bank, declares the tracked hash, and
  separates fourteen non-source constant coefficient matches from zero
  same-clock physical signature matches.  Its exact fixed-rail Frobenius
  no-go survives unchanged and is now independently audited and proved; it does not
  constrain THM-2742, prove THM-2744, or decide the joint `C13` coefficient.

## MISTAKE-308 (2026-07-28, THM-2636/2692 Kummer field typing) -- a spectral nonsquare was identified with its physical square before pullback

- **What was recorded:** the first promoted THM-2692 text wrote
  `mathscr H=rho^3*zeta/t^3=Z=T^2` inside the spectral-curve discussion and
  immediately called `mathscr H` nonsquare.  THM-2636 used the same compressed
  equality while separately describing the intended function-field embedding.
- **Why it was wrong:** `mathscr H` lives first in the abstract spectral
  function field, whereas `Z,T` live in the physical rational function field.
  A literal equality in one field would make `mathscr H` a square and
  contradict the odd-valuation argument.  The missing symbol was the
  hypothetical trajectory embedding, not a missing mathematical hypothesis.
- **Repair / survivor:** define `mathscr H=rho^3*zeta/t^3` on the spectral
  curve.  Its five valuations `-3` make it nonsquare and its connected double
  cover has genus at least two.  Only after a nonconstant physical trajectory
  gives `phi:K(C)->C(x)` do the equations assert
  `phi(mathscr H)=Z=T^2`; the physical `T` then lifts `phi` to that cover, giving
  the contradiction.  THM-2636 and THM-2692 now state this typed pullback
  explicitly.  Their exact computations and closure conclusions are unchanged.

## MISTAKE-306 (2026-07-28, concurrent THM-2693 reservation) -- a later stub duplicated a live theorem namespace

- **What was recorded:** commit `2e3a42992` first reserved THM-2693 for the
  odometer delayed-tail theorem.  Commit `1fee894c1` later added a distinct
  mixed-dilation/slope-seven empty stub with the same YAML ID and theorem
  number before observing the first push.
- **Why it was wrong:** a filename reservation is shared only after its push;
  concurrent local scans do not serialize namespace allocation.  Two current
  files declaring the same theorem ID make every bare citation ambiguous even
  when both are honest empty stubs.
- **Repair / survivor:** the first claimant remains THM-2693 and has since
  been proved.  The later mixed-word stub is coherently moved to the freshly
  checked THM-2694 namespace; its mathematical intent and unproved status are
  unchanged.  Re-fetch immediately before reserving, and yield to the first
  pushed claimant after a race.

## MISTAKE-307 (2026-07-28, THM-2588 separation constants) -- the threshold for rho+1 was reported as the threshold for rho

- **What was recorded:** with descent ratio `rho=h/w`, THM-2588 called
  `134,99,...,15` the sharp integer ratio thresholds and called `533/4`
  the binding real threshold in the one-fold case.
- **Why it was wrong:** the snap loss is `1/[2(rho+1)]`.  Solving
  `1/[2(rho+1)] <= 1/13-3/41=2/533` gives
  `rho+1>=533/4`, hence `rho>=529/4`.  The script correctly subtracted one
  internally and then added it back when printing.  The same display error
  shifted every row of the table by one.
- **Strongest survivor / repair:** the proof strengthens.  The exact real
  one-fold threshold is `529/4`; the legal integer thresholds for
  `k=1,...,12` are `133,98,84,74,65,57,50,42,35,28,21,14`.  Uniform
  separation `rho>=133` closes every legal depth, and every gap family obeys
  `v_max/v_(2)<529/4`.  The fold identity, snap inequality, and all existing
  certificates remain valid.

## MISTAKE-305 (2026-07-28, THM-2588 fold-cascade boundary) -- a thirteenth fold was assigned to a 13-speed family

- **What was recorded:** THM-2588 and its referee printed a constants table
  for `k=1,...,13` and claimed the uniform separated-tower conclusion over
  that range.  The final row used the formal expression
  `1/(14-13)=1` as an empty-body floor.
- **Why it was wrong:** a 13-speed family has only twelve consecutive ratios
  and only twelve legal snap-fold removals, ending at a singleton.  There is
  no thirteenth partner, no thirteenth fold modulus, and no empty-family
  lonely-runner floor.  The script's `k=13` row was not used by any
  certificate and was internally inconsistent with its own special-case
  comment.
- **Strongest survivor / repair:** the fold identity, snap lemma, all
  certificates, T-A/T-B, and every constants row `k=1,...,12` are unchanged.
  THM-2588, its companion, and the stored transcript now quantify only over
  this legal range.  When iterating a peel/fold recurrence, count actual
  edges of the ordered speed list, not vertices or the stationary runner.

## MISTAKE-304 (2026-07-28, THM-2688 hardening) -- a concurrent script/output hardening left stale declared hashes

- **What was recorded:** the hardened THM-2688 frontmatter declared primary
  script/output SHA-256 values beginning `4013f471` / `3f1c7ad9`.
- **Why it was wrong:** the same incoming commit changed the primary companion
  and transcript after those digest values had been calculated.  The committed
  LF-normalized bytes instead hash to `16685599...` / `fecc6762...`.  This was
  a reproducibility-metadata race, not a mathematical or transcript mismatch:
  normal and optimized runs both reproduce the committed output byte for byte.
- **Repair / rule:** THM-2688 now records the hashes of the actual hardened
  files, and its independent referee has separate frozen hashes.  After a
  concurrent rebase or any hardening edit, recompute digests from `HEAD` and
  replay both interpreter modes; never carry pre-edit hashes across the
  checkpoint merely because the visible theorem statement is unchanged.
## MISTAKE-309 (2026-07-28, THM-2683 evidence portability) -- checkout bytes and Python's integer-print cap leaked into exact replay

- **What failed:** THM-2683 declared LF-normalized hashes for its THM-2636
  and THM-2671 executable dependencies but called `sha256(read_bytes())`.
  The promoted companion therefore aborted before its mathematics on a
  Windows CRLF checkout, even though both normalized dependency hashes were
  exactly the declared values.
- **Repair:** dependency text is normalized from CRLF to LF before hashing;
  the theorem now describes LF-normalized transcript/dependency equality.
  A retained optimized replay then exposed a second platform boundary: Python
  3.11+ refused to stringify the already-computed BCDE certificate because
  one exact integer exceeded its default 4,300-digit display cap.  The trusted
  companion now explicitly disables that display cap before forming its
  declared digest.  These are evidence-portability defects only: no algebraic
  certificate, stored output, or theorem conclusion changed.
- **Rule:** a declared normalized hash basis must also be implemented at every
  dependency lock, not merely when frontmatter hashes are computed; a digest
  serializer must also declare any runtime size limit that its exact objects
  exceed.
## MISTAKE-303 (2026-07-28, THM-2648 edge-thinning sharpness) -- a minimum relative to the affine chart was promoted as the unrestricted minimum

- **What was claimed:** the first repaired THM-2648 proved that no
  two-edge change of the matched-wall affine rainbow remains rainbow, found a
  three-edge repair, and called its fourteen-edge union globally sharp.
- **First failed implication:** the hostile fixed the affine chart before
  minimizing.  Two nonlinear rainbow charts can differ in two positions even
  though a transposition of the affine chart cannot.  In the parallel normal
  form the target vectors
  `(2,3,4,6,8,5,10,7,12,9,11)` and
  `(2,12,4,6,8,5,10,7,3,9,11)` have disjoint holes `{0,2}` and `{6,9}`, share
  nine edges, and have union thirteen.  Subtracting one from every target
  gives the antiparallel witness, and affine transport covers every matched
  endpoint pair.
- **Strongest survivor / repair:** fourteen is sharp **conditional on
  retaining the affine chart**; its three-edge alternating six-cycle is the
  local `K_(3,3)`/`C3 semidirect C2` mechanism.  Thirteen is the unrestricted
  sharp value, since distinct bijections cannot differ at exactly one source;
  it is attained by the displayed nonlinear pair through an alternating
  four-cycle.  The resulting binary/ternary fork is structural: the binary
  switch is globally thinner but loses the affine reference, while the
  ternary switch pays one edge to retain it.  Always state which baseline or
  symmetry stratum is fixed before calling a local repair minimal.

## MISTAKE-302 (2026-07-28, THM-2648 rainbow atlas census) -- a selected nonlinear atlas was reported as the universe of nonlinear charts

- **What was claimed:** promoted THM-2648 said that there are exhaustively
  `11,154` affine and `1,014` nonlinear rainbow charts.  Its companion had in
  fact enumerated every affine chart but had only *constructed one chosen
  nonlinear repair per matched-step pair*.
- **First failed implication:** counting the members selected by a complete
  covering atlas does not count every chart that could occupy one atlas slot.
  In each matched normal form there are `50,283` rainbow matchings and
  `43,122` nonlinear second charts compatible with the affine hole pair.  A
  minimal same-hole witness is the pair of target vectors
  `(2,4,11,5,6,7,8,9,3,10,12)` and
  `(2,10,3,5,6,7,8,9,11,4,12)`, both with holes `{6,9}`.
- **Strongest survivor / repair:** the exact `11,154` affine census, existence
  of a two-chart cover for every pair, two-chart minimality, `1^4 2^9`
  incidence, and charged energy all survive.  THM-2648 now calls `1,014` the
  chosen nonlinear-atlas size.  Its selected matched-wall chart was sharpened
  to an affine-anchored alternating six-cycle: three moved edges give a sharp
  fourteen-edge union relative to that anchor, and the affine chart plus the two inverse repairs are
  the three one-factors of `K_(3,3)` with local `C_3 semidirect C_2=S_3`
  symmetry.  The occurrence-level chart bit, not carry support or residual
  reflection, is the missing physical selector.  Never label a constructed
  atlas count as a universe count without enumerating the ambient objects.

## MISTAKE-301 (2026-07-28, THM-2621 trace--Liouville sidecar) -- closed was mistaken for possibly nonexact after forgetting the polynomial source primitive

- **What was assumed:** the first promoted form of THM-2621 correctly proved
  that `omega_F=Tr(x dy)-4 kappa^(-1)u dv` is a closed rational one-form, then
  proposed that a polynomial degree-four Keller map could have a nonzero
  trace--Liouville class or divisor residue along the Jelonek boundary.
- **First failed implication:** the trace was taken before using the stronger
  source fact.  On `A^2`, the polynomial form
  `theta_F=x dy-kappa^(-1)P dQ` is closed by the Keller equation and hence has
  a polynomial primitive `H_F`.  Trace commutes with differentiation, so
  `omega_F=d Tr(H_F)` is exact.  Every divisor residue is therefore zero for
  every actual polynomial Keller realization, regardless of monodromy.
- **Strongest survivor / repair:** the inverse-quartic congruence, resultant
  pole ledger, depression warning, affine exact-change formula, and curvature
  calculation all survive.  The corrected THM-2621 proves that the trace
  primitive is regular off the Jelonek set and has only second-kind derivative
  poles there.  The missing coordinate is branchwise: every normalization
  branch obeys `Res(x dy)=sum_j j x_(-j)y_j=0`, while in the `D_4` lane the
  anti-invariant potential `H_F-tau(H_F)` survives the trace.  For an abstract
  rational inverse-spectral pair, a nonzero traced residue is an obstruction
  to polynomial realization, never a feature of a polynomial Keller map.
  The repaired theorem strengthens this to a localized power-sum potential
  and, in the `D4` lane, an exact opposite-pair potential before the final
  quadratic trace.  A rational `D4` one-sheet-defect family realizes arbitrary
  residue despite the PDE, pole law, and polynomial companion coefficients;
  a second `D4` family has exact base trace but cancelling nonzero pair
  residues.  These are minimal hostiles to recovering polynomial realization
  from the abstract inverse-spectral data.  Always pull a quotient-side closed form back to the claimed
  affine source before promoting it to an independent boundary invariant.

## MISTAKE-300 (2026-07-28, THM-2618 candidate scope) -- the canonical selected source was conflated with a distinct target-informed head

- **What was assumed:** the first candidate form of THM-2618 proved a lawful
  product-torus orbit for THM-2537's canonical selected source and then said
  that this unfreezes the selector inside THM-2573's `w_(N,h)`, closing its
  target-action type gate.
- **First failed implication:** both objects choose a marked root, but they are
  not the same Boolean event.  THM-2537 selects the occupied tail
  `s_tau(e)` from the complete old-carrier root word.  THM-2565/2569 instead
  use the target-informed head `A=T_delta`; its marker comes from a separate
  `k_a` failure mask and slope stratum, and `w_(N,h)=A_h(x)A_h(T^N x)` also
  contains a genuinely later occurrence.  Covariance of the first selector
  says nothing formal about covariance of the second.
- **Strongest survivor / repair:** THM-2618's product-torus restriction,
  Boolean Mobius factorization, dipole covariance, proper-profile hole,
  complete function-level source spectrum, and free profile-atom orbit all
  survive independent hostile audit.  They form a lawful source-side model
  and ingredient only.  Closing THM-2573 still requires a labelled
  whole-event lift of `A_h(x)A_h(T^N x)`, including every failure-mask,
  old-danger/safe, selector, and later-occurrence factor before any target DFT
  or Abel smoothing.  Shared vocabulary such as "selector" is never a map;
  equality of the underlying truth functions must be proved before transport.

## MISTAKE-299 (2026-07-28, THM-2603 proof-graph overwrite) -- a new theorem reused an already proved canonical ID

- **What happened:** a concurrent PSL2(F13) transition audit replaced the
  proved and independently audited THM-2603 Hurwitz root-owner atlas, together
  with its exact companion, under the same theorem ID.  The incoming result is
  mathematically distinct: it studies retained-target normalizers, natural
  chart holonomy, and a Bruhat-square target-zero wall.
- **First failed implication:** thematic continuity and reuse of the same
  matrices do not make two theorem statements the same proof-graph node.  A
  later, stronger-looking obstruction cannot silently discard the earlier
  proved projective atlas, its normalized root-scale guardrail, or its exact
  reproducibility paths.
- **Repair:** the proved Hurwitz/root-owner atlas and its original companion
  remain `THM-2603`.  The incoming PSL2(F13) retained-target theorem, script,
  and output are preserved under the freshly all-ref-checked `THM-2619`
  namespace, still at candidate status.  Cross-links may be added after audit;
  neither result is treated as a replacement for the other.

## MISTAKE-298 (2026-07-28, concurrent THM-2611/2615/2617/2618 reservations) -- a local reservation broadcast was mistaken for a shared proof-graph reservation

- **What happened:** the endpoint-pair parabolic-transvection theorem reserved
  `THM-2611` after a fresh `origin/main` and all-remote-ref scan, and broadcast
  that choice.  Its empty stub was committed locally but could not yet be
  pushed while stronger concurrent theorem edits were being integrated.
  Before the next fetch, `origin/main` assigned `THM-2611` to the distinct
  principal-deck holotopy theorem.  The same failure then recurred after the
  endpoint theorem was moved locally to `THM-2615`: before that move could be
  pushed, the shared branch assigned `THM-2615` to the physical diagonal
  toric-kernel theorem and `THM-2616` to a Keller reduction.  A third local
  move to `THM-2617` again lost a push race when the Keller reduction was
  renumbered there and `THM-2616` was reused for the cross-time diagonal.  A
  fourth move to `THM-2618` met two independently pushed stubs bearing that
  same ID: the selected-source product-torus lift and the affine holotopy
  spectrum.
- **First failed implication:** a local commit plus coordination broadcast is
  not a proof-graph reservation.  Only a pushed, current shared stub occupies
  the namespace; an earlier clean scan says nothing about assignments made
  before the local branch rejoins `origin/main`.
- **Repair:** the principal-deck and physical diagonal-kernel theorems retain
  `THM-2611` and `THM-2615`, respectively, while the two later incoming
  theorems retain `THM-2616` and `THM-2617`; the shared `THM-2618` collision is
  left visible for its authors to repair.  The unpublished endpoint theorem,
  script, and output are coherently moved to the freshly all-ref-checked
  `THM-2620` namespace before promotion.  Mathematics is unchanged.  Future
  sessions that cannot push a reservation immediately must treat its ID as
  provisional, repeat the full remote/YAML/history scan immediately before
  integration, and prefer an immediate small pushed reservation checkpoint
  over another local-only rename.

## MISTAKE-297 (2026-07-28, degree-four Keller monodromy audit) -- Smith theory was applied across the omitted Jelonek divisor

- **What was assumed:** T1549 and the first THM-1375/1440/2465 formulations
  argued that the deck group of the finite etale cover over
  `A^n minus A_F` acts freely on the contractible source `A^n`.  Smith
  fixed-point theory was then used to require a self-normalizing point
  stabilizer, leaving only `A4/S4` at degree four and hence no proper
  intermediate field.
- **First failed implication:** a deck transformation is regular and free on
  the open source above `A^n minus A_F`; it need not extend regularly, let
  alone freely, across the omitted divisor.  The toy map `z -> z^2` on
  `C^*` already has a free deck involution which extends to `C` with the
  fixed boundary point `0`.  Contractibility of the completed source
  therefore cannot be fed to Smith without a separate extension theorem.
- **Minimal algebraic hostile:** `T^4-2` has transitive Galois group `D4`.
  Its quartic root field contains `Q(sqrt(2))`, while its squared-pair
  resolvent `W(W^2+8)` exposes the different quadratic
  `Q(sqrt(-2))`.  Thus degree four can have a proper intermediate field and
  the resolvent quadratic need not be that field.
- **Strongest survivor / repair:** Campbell's Galois-Keller theorem still
  excludes degree two and the regular degree-four groups `C4,V4`.  A wild
  degree-four Keller cover is restricted to `D4,A4,S4`; maximal
  point-stabilizer arguments apply only in the `A4/S4` lanes.  The useful
  primitive-coordinate conclusion survives in `D4` for a different reason:
  a quartic `D4` root field has only one proper intermediate, so if
  `K_4=K(x,y)` then at least one of `x,y` generates `K_4`.  The `S4/A4`
  tower statements of THM-2465/2598 remain valid when lane-scoped, while
  `D4` is restored as an open branch with a reducible `1+2` resolvent
  algebra.  Historical Smith-selection tables are provenance, not current
  canon, unless an extension-across-Jelonek theorem is supplied.
  More precisely, the generic involution extends over the finite Zariski-main
  normalization and is polynomial exactly when it preserves the open affine
  source.  Since normalization ramification is invariant, any surviving `D4`
  map needs an unramified missing boundary divisor exchanged with an included
  divisor over the same target component; normalizer data alone cannot see
  that ownership sidecar.
- **Later superseding repair (THM-2633):** restoring `D4` was the correct
  response to the invalid Smith proof, but it was not the final classification.
  THM-2633 excludes `D4` by a different mechanism: etale openness and the
  absence of nonconstant units on affine space give a fixed affine sheet for
  every boundary inertia group, while Kummer support forbids a nonzero
  monodromy quotient that kills every point stabilizer.  Equivalently, a
  Keller point stabilizer normally generates `G` (hence surjects onto
  `G^ab`); the point stabilizer of `D4` normally generates only its proper
  source-deck kernel.
  No deck transformation is extended across the Jelonek divisor, so this does
  not revive the failed implication recorded here.  The current degree-four
  list is `A4,S4`; the `D4` ledgers remain conditional hostile controls.

## MISTAKE-296 (2026-07-28, concurrent theorem reservation) -- a later empty stub reused the live THM-2604 accessibility ID

- **What happened:** the independently audited future-root accessibility
  theorem was reserved and proved as
  `THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary.md`
  on a public session branch.  Before that branch was merged to `main`, a
  later session reserved the distinct alternative-rail clock-collapse
  candidate under the same numeric ID.  The filenames did not conflict, so
  both YAML IDs briefly appeared on `main` after integration.
- **Why it was wrong:** theorem IDs are proof-graph keys.  An empty reserved
  stub cannot share the key of an earlier public, proved, audited theorem;
  checking only `origin/main` misses concurrently published session branches.
- **Repair:** the proved accessibility theorem retains `THM-2604` and its
  full slug.  The later `RESERVED / UNPROVED EMPTY STUB` is moved without any
  mathematical promotion to
  `THM-2608-alternative-rail-clock-collapse-and-missing-transition-index.md`.
  Future reservations must scan all remote refs immediately before commit;
  current citations continue to use theorem ID plus slug.

## MISTAKE-295 (2026-07-28, THM-2594 hostile control) -- constant-column erasure was mislabeled as fixed-root slaving

- **What was claimed:** because replacing every deep column by the same
  marginal made the centred contraction vanish, the first THM-2594 draft
  said that its nonzero value was carried specifically by the affine coupling
  `theta=t-2u`.
- **First failed implication:** the coded control copied
  `sum_theta A(ell,theta)` into every column.  This erases the deep label; it
  is not the genuine `u`-independent experiment in which one fixes an
  absolute root `t` and, on the `u` fibre, reads column `theta=t-2u`.
- **Minimal hostile:** for the exact same table, define
  `A_abs(ell,t)=sum_(u,q)N(u,q,ell,t-2u)`.  Its doubly centred interaction is
  nonzero and `Psi_(1,1)(1,1)` is nonzero exactly in `Q(zeta_91)`; the tested
  primitive companions are nonzero as well.
- **Strongest survivor / repair:** THM-2594 still proves that its defined
  rational `theta`-indexed table is non-replica and hence fires all `5,184`
  primitive coefficients.  The beta-zero and constant-column controls remain
  valid, but they do not identify the cause of nonvanishing.  The fixed-root
  hostile is now computed explicitly, and every statement that the signal is
  uniquely carried by slaving is removed.

## MISTAKE-294 (2026-07-28, concurrent theorem reservation) -- the audited theta contraction collided a third time at THM-2590

- **What happened:** after the second namespace repair, the theta contraction
  was provisionally moved to `THM-2590` on the session branch.  While its
  stalk semantics were being audited, `origin/main` independently assigned
  `THM-2590` through `THM-2593` to four earlier-pushed results.
- **Repair:** a fresh remote scan found `THM-2594` free, so the contraction is
  now `THM-2594-realized-theta-slaved-contraction-at-the-r5-window.md`.
  The earlier numbers retain their first-pushed meanings.  This third repair
  is also why current citations use ID plus slug rather than ID alone.

## MISTAKE-293 (2026-07-28, THM-2594 stalk typing) -- one ancestry fibre was described as one circle point, and an exact return exponent as an eventual lawful clock

- **What was written:** the first theta-slaved contraction draft said that
  all THM-2449 response factors were evaluated at the same stalk point and
  called `k=6` the smallest lawful clock in its congruence class.
- **Why it was wrong:** in the exact Boolean expansion, the owner cell is
  evaluated at the collision-chart node `w_u`, the deep and delayed-word
  factors at its descendant `X_(u,a)`, and the source factor at the linked
  node `Y_(q,e')`.  They share one base `y` and one finite ancestry fibre,
  but not one circle point.  Also, THM-2449 only licenses its eventual clock
  law for `k>=k_0`; the residue calculation `k=6 mod 16380` does not prove
  `6>=k_0`.
- **Strongest survivor / repair:** the `k=6` table and its nonzero
  cyclotomic contraction are computed directly and exactly, so they need no
  eventual-overlap inference.  Sheet-independence lets every linked-node
  factor be pulled through the finite Boolean ancestry sum, proving a genuine
  same-fibre product before marginalization.  THM-2594 is therefore an exact
  auxiliary multilevel-stalk contraction, not THM-2449's one-point lawful
  response table and not a physical target current.

## MISTAKE-292 (2026-07-28, concurrent theorem reservation) -- the first repair of the theta-slaved contraction collided again at THM-2588

- **What happened:** MISTAKE-291 repaired the theta-slaved contraction's
  original duplicate `THM-2581` by moving it to `THM-2588` on a session
  branch.  Before that branch reached `origin/main`, another live session
  assigned `THM-2588` to
  `THM-2588-fold-cascade-separated-towers-gap-empty-all-k.md`; the shared
  hypothesis index also explicitly reserved `THM-2589`.
- **Why it was wrong:** checking the local tree before a concurrent repair was
  not enough.  A theorem-ID repair must be rebased onto the current shared
  namespace before its replacement ID is treated as canonical.
- **Repair (superseded by MISTAKE-294):** after fetching and rebasing onto the
  then-live shared surface, the theta-slaved contraction was renamed
  `THM-2590`; a third concurrent collision required the final `THM-2594`
  name.
  `THM-2588` remains the fold-cascade theorem and `THM-2589` remains reserved
  by HYP-9075.  Historical messages and the intermediate MISTAKE-291 wording
  remain provenance; current references use `THM-2594`.

## MISTAKE-291 (2026-07-28, concurrent theorem reservation) -- the realized theta-slaved contraction reused live ID THM-2581

- **What happened:** the candidate theorem for the realized depth-five
  theta-slaved contraction was committed as a second `THM-2581`, after
  `THM-2581-b-word-depth-five-owner-clock-host-and-reflection-breaking.md`
  had already been proved, independently audited, and routed under that ID.
- **Why it was wrong:** theorem IDs are namespace keys, not mathematical
  analogies.  The new contraction depends on the existing depth-five host;
  sharing its ID makes dependencies and status searches ambiguous even though
  the filenames differ.
- **Repair (superseded by MISTAKE-292):** the first repair renamed the later
  candidate `THM-2588`; a second concurrent collision required the
  intermediate `THM-2590` name, followed by the final `THM-2594` repair in
  MISTAKE-294.
  Historical broadcast messages retain their original text
  as provenance only.  The original `THM-2581` remains the canonical
  owner-clock host, while `THM-2594` remains a separately audited candidate
  contraction.

## MISTAKE-290 (2026-07-28, THM-2576 resultant convention) -- a PRS degree reorder lost the odd-by-odd resultant swap sign

- **What was written:** THM-2576 equation (14) displayed the standard
  resultant as `Res_X(E,Q)=-a^8 c^18 S^8 H`, with `deg_X(E)=3` and
  `deg_X(Q)=15`.
- **First failed implication:** the computation called SymPy's PRS
  `resultant(E,Q,X)` with the lower-degree polynomial first.  That routine
  reordered the inputs but did not restore the factor
  `(-1)^(3*15)=-1`.  At the exact target `(a,b,c)=(1,1,1)`, the PRS value is
  `-198662708168284531608320`, while the standard `18 x 18` Sylvester
  determinant and root-product resultant are
  `+198662708168284531608320`.
- **Strongest survivor:** the primitive polynomial `H`, its coefficient
  ledger, irreducibility, image hypersurface `V(H)`, and
  `S_(F o F)=V(LH)` are unchanged: multiplying a defining equation by `-1`
  never changed the ideal or zero set.  The independently derived norm law
  `prod_(q in F^-1(t)) L(q)=H(t)/(64L(t))` has the positive sign and remains
  correct.
- **Repair:** compute the higher-degree-first PRS value and restore
  `Res(E,Q)=(-1)^45 Res(Q,E)` explicitly.  In the standard convention the
  exact identity is `Res_X(E,Q)=+a^8 c^18 S^8 H`.  Whenever resultant input
  degrees are reordered, audit both the swap parity and one Sylvester or
  root-product control before using constants or square classes.

## MISTAKE-289 (2026-07-28, recurrent concurrent THM-2570 collision) -- a later theorem reused an ID after the earlier collision repair

- **What was done:** the LRC word-depth/owner-clock theorem was committed as
  `THM-2570` after the Jelonek cusp-cylinder normalization had already retained
  that ID in the repair recorded by MISTAKE-288.  The live tree again contained
  two distinct theorem files with one YAML ID.
- **First failed implication:** a prior namespace repair is not a reservation
  service for later concurrent writers.  Choosing an ID from an earlier local
  view, without fetching and checking the complete remote filename/YAML
  namespace immediately before the commit, did not establish uniqueness.
- **Repair:** the earlier Jelonek theorem remains `THM-2570`; the already
  reserved Cayley and deep-energy results remain `THM-2571` and `THM-2572`;
  the pending Abel-normal endpoint result owns `THM-2573`; the oriented-tooth
  component-holonomy reservation owns `THM-2574`; and the later
  word-depth/owner-clock theorem is renamed `THM-2575`.  Its scripts, outputs,
  mathematical statement, and exact constants are unchanged.  Future
  reservations must fetch and scan the prefix `id: THM-` and matching
  filenames immediately before both reservation and push.

## MISTAKE-288 (2026-07-27, concurrent THM-2570 reservation collision) -- an end-anchored ID scan missed a CRLF frontmatter line

- **What was done:** the LRC deep-augmentation energy theorem was reserved as
  `THM-2570` after a filename/history scan appeared clear.  The earlier
  `THM-2570-jelonek-cusp-cylinder-normalization-and-conductor.md` reservation
  was already present, so the live tree temporarily contained two distinct
  files with the same YAML ID.
- **First failed implication:** the decisive YAML search used the pattern
  `^id: THM-2570$`.  On this mixed-line-ending worktree, the earlier file's
  CRLF carriage return prevented the end anchor from matching.  A clean
  command result therefore did not imply a clean namespace.
- **Repair:** the later LRC theorem and both exact artifacts were coherently
  renamed to `THM-2572`; the first Jelonek claimant keeps `THM-2570`.
  Namespace audits must search filenames, remote history, and YAML prefixes
  without a line-ending-sensitive end anchor (or normalize line endings
  first).  No mathematical statement or computed constant changed.

## MISTAKE-287 (2026-07-27, HYP-9033/HYP-9065 cusp and elliptic pullbacks) -- square coefficient factors were discarded in global zero-set claims

- **What was claimed:** from the exact identity
  `S^2-T^3=27c^2L`, HYP-9033 said that the Jelonek hypersurface
  `V(L)` was the full inverse image of the cusp `S^2=T^3` under
  `(a,b,c) |-> (S,T)`.  HYP-9065 then made the parallel inference from
  `-4p^3-27q^2=-2^4 3^6 a^2L`: it identified the global singular-
  Weierstrass locus with `V(L)`, identified `p=q=0` only with the omitted
  curve, and promoted the resulting `3/1/0 = smooth/nodal/cuspidal` count
  dictionary globally.
- **First failed implication:** multiplication by `c^2` is invertible only on
  the chart `c!=0`.  Globally the scheme-theoretic inverse image is
  `V(c^2L)`, its effective divisor is `2V(c)+V(L)`, and its underlying set is
  `V(c) union V(L)`.  The affine target `(a,b,c)=(1,0,0)` is the minimal clean
  hostile: it maps to the smooth cusp point `(S,T)=(8,4)` while `L=16`.
  Likewise the Weierstrass degeneration set is `V(a) union V(L)`.  More
  sharply, `V(p,q)=V(a,b) union E`: every target `(0,0,c)` gives a cuspidal
  Weierstrass cubic but THM-2546 gives one affine `F`-preimage `(c/2,0,0)`,
  not zero.  Hence the global count dictionary also fails.
- **Strongest survivor:** after localization at `c`, or equivalently saturation
  by `c`, the claimed identification is exact.  The cusp-point fibre
  `(S,T)=(0,0)` is still exactly the omitted rational curve because `T=0`
  forces `bc=4/3` and hence `c!=0`.  The plane `c=0` is an honest affine target
  plane with an explicit finite section, not a compactification boundary.
  On `a!=0`, the elliptic statements are repaired exactly: degeneration is
  equivalent to `L=0`, `p=q=0` is equivalent to `E`, and the `3/1/0`
  fibre-count dictionary does match smooth/nodal/cuspidal type.
- **Repair:** THM-2566 records the saturated `c`-chart together with the
  companion `a`-chart from THM-1335, their parasitic planes, the punctured
  two-chart cover, and the exceptional origin.  Never pass from
  `h^m f=0` to `f=0` globally without localizing at or saturating by `h`.

## MISTAKE-286 (2026-07-27, THM-2560 bridge typing) -- a specific `r=5` slaving ansatz was attributed to THM-2512's generic bridge test

- **What was written:** THM-2560 called the proposed
  `theta=t-2u` toothpick contraction “THM-2512 Section 5's `u`-slaved
  bridge” and said the bridge required collision index `r=5`.
- **Why it was wrong:** THM-2512 Section 5 asks generically whether a signed
  toothpick contraction can be realized on one common Boolean
  owner/deep-ancestry fibre.  It contains no collision index, `theta`, or
  `r=5` hypothesis.  Those coordinates belong to THM-2471's deep-root sidecar
  and to THM-2560's particular candidate realization.
- **Repair:** the exact `r=3` computation proves only that this particular
  `u`-slaved candidate is unavailable because the deep phase is base-only.
  It does not refute THM-2512's bridge test, a base-only coupling, or another
  common-ancestry realization.

## MISTAKE-285 (2026-07-27, THM-2560 phase-cone verifier) -- a cyclotomic sign error bound used coefficients from a different representative

- **What was done:** the phase-cone companion evaluated the original
  `Z[zeta_13]` coefficient vector against a fixed sine/cosine table, but bounded
  its rounding error by the `l1` norm of the conjugate-difference or
  conjugate-sum vector used for the algebraic norm separation.  Those norms do
  not control an arbitrary representative modulo `1+zeta+...+zeta^12=0`.
- **First failed implication:** a small conjugate-sum/difference vector does not
  bound the numerical error contributed by a large cyclotomic-kernel or
  symmetry component of the unreduced input.  The original trigonometric
  self-checks also tested consistency, not a rigorous approximation interval.
- **Strongest survivor:** every zero test and all `168` phase-cone verdicts are
  unchanged.  The repaired ordinary and optimized runs again give `168/168`
  strict failures and the same first cone-breaking witness.
- **Repair:** evaluate the conjugate-difference/sum vector itself, certify all
  thirteen fixed sine/cosine values by outward-rounded rational intervals from
  Machin plus alternating Taylor bounds, and use the matching `B` both for the
  table-error budget and the degree-twelve algebraic norm separation.  This is
  now enforced inside the exact companion before any angular comparison.

## MISTAKE-284 (2026-07-27, concurrent THM-2560 reservation collision) -- an empty reservation and a proved candidate shared one theorem ID

- **What happened:** the empty canonical-duty-commutator reservation and the
  independently developed second-rung diagnostics were both pushed with
  frontmatter ID `THM-2560`.
- **Why it was wrong:** theorem IDs are a shared namespace even when one file
  is only a `RESERVED / UNPROVED EMPTY STUB`; status separation does not make a
  duplicate identifier safe.
- **Repair:** the already substantive diagnostics retain `THM-2560`; the empty
  reservation is renamed `THM-2562` after a fresh fetch and repository-wide ID
  check. Cite theorem ID plus slug, and recheck the remote namespace immediately
  before every reservation push.

## MISTAKE-283 (2026-07-27, THM-2550 consequence scope) -- two exact nondegeneracy computations were treated as one common-root packet

- **What was claimed:** after separately computing positive owner-loop drift
  and a non-replica lawful response table on the same canonical typed row and
  word label, THM-2550 said the drift fed THM-2466 (32) and THM-2471 (48)
  directly in the same semantic-word-filtered packet.
- **First failed implication:** Part (A) uses the `k=2` owner packet, while
  Part (B)'s lawful ANOVA law uses clocks `k=6 mod 16380`.  A common row and
  word label do not identify their clocks, packets, marked triangles, or
  ancestry atoms.  Independently, THM-2466 assumes a supplied common oriented
  root base and positive service density, and THM-2471 does not identify its
  first-collision root with THM-2365's endpoint/deep root.  None of these
  hypotheses is created by drift positivity.
- **Strongest survivor:** the independent hostile audit accepts both exact
  computations: both drift tensors are positive; both large-clock ANOVA
  defects are nonzero in all `91` cells and ignite all `5,184` primitive cut
  modes.  THM-2365 branch 1 applies to Part (A), and THM-2557 conditionally
  transports each Part (B) signed interaction without cancellation.  These
  are separate typed non-cover controls, not one physical current.
- **Repair rule:** equality of row, owner, or word labels is not a fibre
  product.  Before composing two positive/nonzero certificates, name and
  match their clock, packet, marked triangle, common base, root gauge, and
  required semantic roles.  If any coordinate is absent, retain the results
  componentwise and state the missing intertwiner.

## MISTAKE-281 (2026-07-27, THM-2547 cut/root pairing) -- a product of unrelated Fourier controls was promoted as a physical current

- **What was claimed:** THM-2547 multiplied the synthetic weighted
  root-correlation control from THM-2471's companion by THM-2506's canonical
  two-row high-septimal defect, found all `432/432` external pairings nonzero,
  and called the result the first physical cut-character current on one
  owner-retaining ancestry cospan. It treated THM-2478's printed `k=alpha`
  label equality as a typed bridge and read `K_{0,beta}=0` as target-root
  exclusion.
- **First failed implications:** the two inputs have no exhibited common
  Boolean base, atomwise product, or physical morphism. THM-2478 (21a)-(21b)
  aligns a future collision colour with one fixed old target-component label;
  it does not identify that temporal sheet with the static cut-row character.
  Moreover, `alpha=0` is the trivial Fourier character, whereas `C(0)` is a
  physical root-coordinate value. Owner and deep margins were calculated
  separately and never occurred in the pairing.
- **Minimal hostile mechanism:** one predecessor packet at root `0` and one
  disjoint arrival packet at root `2` give the Boolean correlation
  `C=delta_2`. All twelve nontrivial root colours are nonzero, and the same
  defect has all `5,184` primitive modes nonzero, but the external dual sum
  vanishes in `108/432` placements (in particular all 36 `(a,beta)` at
  `tau=4`). Here `R_{4,a,c}(2)=0` for all 42 `(a,c)`, so character
  orthogonality forces the whole alpha-sum to zero. Nonzero factors do not
  prevent cross-alpha cancellation.
- **Strongest survivor and repair:** for the specifically displayed synthetic
  `C=(6/25)(delta_3+delta_6+delta_10+delta_12)` and canonical `d`, both dual
  and untwisted external pairings really are nonzero `432/432` in
  `Q(zeta_91)`. A separate direct implementation agrees at the five split
  primes `547,911,1093,2003,2549`; covariance and both torsor collapses also
  survive as facts about this external array. THM-2547 is downgraded to
  `FINITE-EXACT EXTERNAL CONTROL + REFUTED PHYSICAL INTERPRETATION` and gains
  the exact `delta_2` non-universality theorem. No physical current, target or
  owner retention, typed bridge, live-row transplant, or LRC(14) exclusion
  follows.
- **Rule / missing quantifier:** multiply cut and collision data on one common
  ancestry base before marginalization, DFT, and integration. A live-branch
  theorem must quantify over every one of the 165 rows, construct same-atom
  owner/word/deep/cut data and a lawful typed intertwiner, choose a
  row-dependent but ancestry-point-independent cut probe, and prove the
  resulting integrated sum nonzero. Full support alone is never a
  noncancellation theorem.

## MISTAKE-280 (2026-07-27, THM-2545 ID collision) -- a theorem ID was claimed without re-pulling at write time

- What was done: opus canonized "positive cut-character ancestry pairing"
  as THM-2545 while [arrival-hall]'s reservation of THM-2545
  (word-stratified hall arrival criterion) had already reached
  origin/main; the rebase interleaved both files under one ID.
- Why it was wrong: IDs are reserved by first push (AGENTS.md concurrent
  protocol); the free-ID check was done before, not at, commit time.
- Correct framing: opus's theorem is renumbered THM-2547
  (01-canon/theorems/THM-2547-positive-cut-character-ancestry-pairing.md);
  [arrival-hall] keeps THM-2545. Any reference to "THM-2545" in
  opus-2026-07-27 letters/broadcast (MSG to mac-mini/kind-pasteur, the
  close-out broadcast) means THM-2547. Re-pull and re-grep the ID at the
  moment of commit, not merely at drafting time.

## MISTAKE-279 (2026-07-27, THM-2533 Galois gain propagation) -- an absolute-square character integral was called rational

- **What was written:** THM-2533 Section 6 said that expanding the owner-
  weighted absolute-square energy `Gamma(kappa,b,a)` before its character
  sums made every scalar integral rational, and used rationality to explain
  the Galois propagation of three positive gain slices.
- **Why it was wrong:** rational step-cell measures do not make the squared
  modulus of a cyclotomic character sum rational.  Two equal rational branch
  cells in target slots `0,1`, for example, contribute the nonrational real
  cyclotomic factor `|1+zeta_13|^2=2+zeta_13+zeta_13^(-1)`.  The exact referee
  checked the covariance and nonvanishing consequences, not this prose-only
  rationality assertion.
- **Strongest survivor and repair:** the `216` phase ladders and every later
  consequence survive unchanged.  Galois covariance sends a positive seed
  `Gamma` to the corresponding gain entry.  A field automorphism preserves
  nonzero, while that image entry is independently an integral of an
  absolute square against a nonnegative owner, hence is nonnegative in the
  distinguished real embedding.  It is therefore strictly positive.  The
  proof now uses nonzero preservation plus independent nonnegativity, never
  order preservation or rationality.
- **Rule:** a rational partition supporting a character sum does not make its
  quadratic Fourier energy rational.  When propagating positivity by Galois,
  separate the algebraic step (conjugates of nonzero elements stay nonzero)
  from the analytic step (each conjugate has its own manifest absolute-square
  realization and is nonnegative).

## MISTAKE-278 (2026-07-27, THM-2521 complex Hilbert norm proof) -- the real part was omitted from a complex cross term

- **What was written:** THM-2521 equation (15) asserted
  `2 sum_(x<y)<p_x,p_y>=-sum_x||p_x||^2` for real or complex
  Hilbert-valued centred potentials.
- **Why it was wrong:** with a sesquilinear complex inner product the
  unordered sum of cross terms need not be real.  For example
  `(p_1,p_2,p_3)=(1,i,-1-i)` is centred, but the unordered cross-term sum
  has a nonzero imaginary part (whose sign depends only on the inner-product
  convention).  The norm expansion retains twice its **real part**, not the
  complex sum itself.
- **Strongest survivor and repair:** the claimed norm identity
  `||Sp||^2=12||p||^2`, and hence every downstream constant `156` and
  `1872`, was already correct: its proof uses
  `2 Re sum_(x<y)<p_x,p_y>=-sum_x||p_x||^2`.  Equation (15) now states that
  formula explicitly and notes that `Re` may be omitted over a real Hilbert
  space.  The same audit sharpened the aligned Radon map to
  `R+R^*=13I` and `||Rp||>=(13/2)||p||`, with exact companion checks.
- **Rule:** whenever a real quadratic norm identity is extended
  coefficientwise to a complex Hilbert space, write either an ordered
  sesquilinear sum or the real part of the unordered cross terms; do not
  copy the real bilinear expansion verbatim.

## MISTAKE-277 (2026-07-27, concurrent THM-2510 reservation collision) -- remote uniqueness was not rechecked after the final rebase

- **What happened:** the quadratic affine-cut energy theorem was checked
  against the then-current remote namespace and found `THM-2510` free.  While
  its exploratory computation was being developed, a concurrent session
  reserved the unrelated truncated-Radon tight-frame theorem as `THM-2510`
  in commit `dc69ebaa800`.  This session then rebased over that commit but did
  not repeat the YAML-ID search before pushing its own energy reservation in
  commit `dd9e4964336`, leaving two distinct theorem slugs with the same ID.
- **Why it was wrong:** the earlier uniqueness check had expired, and a clean
  rebase does not detect semantic collisions between different filenames.
  The bounded reservation protocol requires checking filename, YAML ID,
  indexes, and fresh remote history *after* the last fetch/rebase and
  immediately before the reservation commit.
- **Repair:** the quadratic energy reservation, companion, and all pending
  citations were moved to the freshly rechecked free namespace
  `THM-2511-affine-cut-quadratic-root-service-and-pair-space-boundary`.
  `THM-2510` now names only the truncated-Radon tight-frame theorem.  No
  proved mathematical statement depended on the transient duplicate.
- **Rule:** after every rebase performed during an ID reservation, rerun the
  remote filename and YAML-ID search before committing, even when an earlier
  check in the same turn found the number free.

## MISTAKE-276 (2026-07-27, THM-2478 promotion audit) -- a factor-wise target action was written as one common physical translation

- **What was written:** after proving the delayed owner graft target-neutral,
  THM-2478's promoted time/charge paragraph represented a lawful THM-2365
  target co-shift as `x -> x-theta/13` on one physical variable.  It then
  claimed that the resulting sheet displacement was the digit read by every
  old deep probe.
- **First failed implication:** THM-2365 acts independently on the atomic
  phases `y_i=w_i x`, with role-dependent `theta_i`; there need not be one
  physical translation of `x`.  Worse, the displayed common shift
  `13^(L-1)theta` is zero modulo the claimed deep quotient
  `13^(L-lambda)` whenever `lambda>=1`, so that formula could not move the
  advertised digit.
- **Strongest survivor and repair:** the graft, all-colour mixing, and sharp
  sheet threshold survive.  On each atomic natural extension,
  `a_i -> a_i-13^(L-1)theta_i (mod 13^L)`.  For the deep role
  `C=13^lambda u`, this induces
  `a -> a-13^(L-lambda-1)u^(-1)theta_C`
  modulo `13^(L-lambda)`, which really moves the top digit of the essential
  sheet.  THM-2478 now states this factor-wise cocycle, makes the Boolean-leaf
  target-neutrality check explicit, and gives a root-uniform two-sheet
  witness.
- **Rule:** never turn a quotient-dual action on several atomic phases into a
  common translation of the orbit variable unless the shift vector is proved
  proportional to the coefficient vector.  Test any claimed moving digit in
  the exact quotient where it is supposed to move before using its semantics.

## MISTAKE-275 (2026-07-27, concurrent THM-2479 reservation collision) -- a pushed reservation did not keep the theorem ID unique

- **What happened:** the degree-twenty-two `B,C` plane closure reserved and
  pushed `THM-2479` in commit `7b15601fb81`.  Before its proof was promoted,
  a concurrent session canonized the unrelated two-colour trichotomy under
  the same YAML ID in commit `9942471dbe3`.  The two filenames differed, so
  Git rebased cleanly while the theorem namespace became ambiguous.
- **Why it was wrong:** a reservation is coordination evidence, not a lock.
  Checking the filename, YAML ID, indexes, and remote history only at initial
  reservation does not protect against a later concurrent reuse.  Any final
  citation to bare `THM-2479` would have had two possible targets.
- **Repair:** the `B,C` reservation and exact companion were renamed to the
  freshly rechecked free namespace
  `THM-2480-degree-twenty-two-BC-plane-hensel-ramification-closure` before any
  proved promotion.  `THM-2479` now refers only to the canonized two-colour
  trichotomy.  The `B,C` mathematics and exact output were unaffected; only
  the namespace changed.
- **Rule:** immediately before every status promotion, fetch and recheck the
  candidate number across filenames, YAML IDs, navigation, and remote history,
  even when the same session already pushed a reservation.  Cite theorem ID
  plus slug so a transient duplicate is detectable rather than silently
  resolved by context.

## MISTAKE-274 (2026-07-26, reserved THM-2440 strict two-comb radius) -- a closed/a.e. handoff radius was assigned to the literal open component

- **What was proposed:** two radius-`1/14` integer danger combs were said to
  cover a centred pullback window of radius at most `15/182`, with equality
  for `{n,13n}`, using the literal strict-open convention needed for a
  pointwise component.
- **Minimal witness and first failed implication:** after gcd normalization,
  take `{1,13}` and the handoff `x=1/14`.  Then
  `||x||=||13x||=1/14`, so both strict combs miss the point.  Their open union
  covers the window only almost everywhere and its literal centred component
  stops at `1/14`.  Passing from an a.e. union identity to a connected open
  cover silently filled the internal seam.
- **Strongest survivor:** the normalized closed radius, equivalently the open
  almost-everywhere radius, is sharply `15/182`, with equality only for
  `{1,13}`.  The literal open connected radius has the different sharp value
  `15/196`, with equality only for `{1,14}`; its first positive tooth strictly
  overlaps the central `1`-tooth.
- **Repair:** THM-2440 proves both statements separately and supplies a
  dependency-free exact referee.  The valid a.e. theorem can support a
  `16/182>15/182` window contradiction only after a lawful argument produces
  coverage of that whole normalized window.  It does not by itself decrement
  a scalar row or close the surviving noncirculant graft.
- **Rule:** for every interval-comb handoff, record three radii separately:
  closed, open almost everywhere, and literal open connected.  Audit every
  equality extremizer at its internal tooth endpoints before transferring a
  component bound.

## MISTAKE-273 (2026-07-26, proposed THM-2434 endpoint cage) -- an almost-everywhere exact tiling does not cover its open handoff seams

- **What was proposed:** in a `b=0,t=5` THM-2427 residual, THM-2430 makes
  the six top words an exact one-fold common-root tiling.  It was then
  proposed that every top-word endpoint must lie in a blocker: otherwise
  another open top word covering the endpoint would overlap the exiting word
  on one side.  Endpoint counting would have implied
  `sum_j gcd(u,c_j)>=2u/3` for every top speed.
- **First failed implication:** the root tiling is an almost-everywhere
  identity.  An exiting open top tooth and an entering open top tooth may
  abut at one common endpoint.  Neither contains that point, while the
  punctured neighborhood is exactly one-fold.  Thus there need not be
  "another word covering the endpoint."  The local hostile is the handoff at
  `x=1/14` from guard speed `H=2` to ordinary speed `q=27`; both words miss
  the common point and all other displayed masks can be absent nearby.
- **Strongest survivor:** away from the high-blocker danger union, every top
  endpoint is either high-blocker-caught or paired with exactly one
  opposite-signed top endpoint.  For an ordinary or guard speed `u` and a
  blocker `c`, with `g=gcd(u,c)` and `A=u/g`, each endpoint-sign branch has
  at most `g ceil(A/7)` blocker hits.  The top-to-top handoff instead obeys
  the exact seam divisibility law of THM-1156.  A valid invoice must retain
  both signed event types.
- **Repair:** THM-2434 remains a superseded unproved stub.  THM-2431 closes
  the intended `b=0,t=5` branch by a parent-mass rounding argument which does
  not require pointwise endpoint coverage.  Any future endpoint theorem must
  use the signed handoff/blocker event current, or explicitly add a literal
  pointwise cover with the correct closed-guard convention.
- **Rule:** an a.e. partition of open interval masks controls punctured
  chambers, not common endpoints.  Before converting a tiling into an
  endpoint cage, state which mask contains the seam point and audit the
  strict/closed convention.

## MISTAKE-272 (2026-07-26, THM-2411 fixed-coefficient projection) -- eliminating `E` turned a constant-field coefficient into a moving function

- **What was done:** the first degree-twenty-two pole-divisor reduction solved
  the residual first flux for `E`, projected from
  `(B,C,D,E,W;Z,y)` to `(B,C,D,W;Z,y)`, and then treated the resulting
  quadratic/square-class family as the live Keller trajectory carrier.
- **First failed implication:** on a trajectory, `E` is one fixed Faber
  coefficient in the algebraically closed constant field. The projection
  preserves the formal quadratic `F_2(Z,y)=0` but forgets the equation
  `P_5(y)=2108304E`. Since `P_5` has fixed leading coefficient `945`, that
  equation makes `y` constant; the second flux then makes `Z,T,q` constant
  and contradicts the genuine deck.
- **Strongest survivor:** the quadratic `F_2`, its sextic discriminant
  `R_6`, the perfect-square classification, and the explicit
  `H_2S_2^2` and squarefree controls are exact algebra on the projected
  universal coefficient cone. They are not live trajectories.
- **Repair:** hostile-audited
  `THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction`
  retains fixed `E` and proves only the full `mathcal A=0` chart empty.
  Its full-mate sidecar retains all lower Faber rows. The complementary
  `mathcal A!=0` chart remains separate.
- **Rule:** solving for a constant coefficient does not authorize forgetting
  that it is constant. Before projecting a coefficient variety, name the
  next consumer and retain every fixed-field coordinate it tests.

## MISTAKE-271 (2026-07-26, THM-2418 fixed-source amendment audit) -- a nonwrapping prefix-block proof was stated with wrapped circular endpoints, and its verifier checked the wrong source

- **What was written:** the fixed-source amendment described
  `[A/B,(A+p)/B)` as a circular source without restricting the displayed
  lift, called `J diag(q)` rank one for every terminal profile, and cited a
  variation check on an unrelated three-cell tuple as verification of the
  thirteen-cell source `[3/13,10/13)`.
- **Minimal witnesses:** for `(B,p,A)=(5,3,4)`, wrapping the cells
  `{4,0,1}` gives residue counts `(1,2,0)`, not one of each.  For `q=0`,
  `J diag(q)` has rank zero.  The unrelated tuple did not test the claimed
  interval at all.
- **Repair:** the universal law now assumes the nonwrapping lift
  `0<=A<A+p<=B` (and analogously `0<=A<A+7t<=13^s`), calls the centred
  half-open interval inversion-even only almost everywhere, and states rank
  at most one, with equality exactly for a nonzero terminal profile.  The
  companion constructs the actual thirteen-cell block and checks its two
  circular jumps, mass, and carry counts through eight clocks.  A separate
  denominator-39 Boolean control realizes the claimed flat/all-six period
  alternation.
- **Rule:** a consecutive integer block is not automatically a consecutive
  block after circle reduction.  State the chosen nonwrapping lift, include
  zero-profile boundaries in rank claims, and make a verifier construct the
  exact object named by the theorem.

## MISTAKE-270 (2026-07-26, THM-2419 affine-shell audit) -- an affine torsor was used without exposing primitivity, Fourier normalization, or the amplitude sidecar

- **What was written:** the first promoted affine-sideband theorem applied
  the full-kernel residue count without explicitly normalizing the selected
  speed row to primitive form, suppressed the quotient-character
  normalization, and described a same-shell residue-zero address as
  sufficient without separately requiring a nonzero labelled cross-packet
  coefficient and phase.
- **Minimal witness / first failed implication:** for `w=(2,2),X=2`, the
  shell residues modulo two are only `(0,1),(1,0)`, while
  `ker(w mod2)` is all four residues.  An address in a residue-zero reference
  fibre also carries no amplitude by itself.  Summing all `M^n` ambient
  twists repeats every quotient character `M` times.
- **Repair:** THM-2419 now assumes primitive `w` (or performs a lawful
  common-gcd normalization), writes
  `C(u)=M^(1-n) sum_[ell] chi_ell(-u)H(ell)` over
  `(Z/M)^n/<w>` and the equivalent `M^-n` ambient formula, and scopes the
  same-shell reference as address-level sufficiency only.  Its exact
  companion tests signed primitive fibres, the nonprimitive hostile,
  quotient-versus-ambient Fourier multiplicity, signed `X`, and the CRT
  reductions.  A physical conclusion still needs a nonzero cross-reference
  coefficient with a common Abel gauge and phase.
- **Rule:** affine-shell counting needs primitivity; quotient Fourier sums
  need their repetition factor; and a surviving address is not a surviving
  current until amplitude and relative phase are transported.

## MISTAKE-267 (2026-07-26, THM-2412 verification and quotient audit) -- optimized transcript equality was not optimized verification

- **What was written:** the promoted THM-2412 companion used sixteen
  truth-bearing Python `assert` statements and cited equality of its normal,
  optimized, and stored transcripts as exact optimized verification.  Its
  prose also displayed the lowering identities without an explicit `k>=1`
  range, called two nested subset half-spaces a partition, and called the raw
  central sign-vector count the gauge-class count.
- **First failed implication:** `python -O` removes every `assert`, so the
  optimized transcript exercised no truth-bearing checks.  Separately,
  `P_n={|S|<=n}` contains `M_n={|S|<n}`; only the reflected strict upper copy
  of `M_n`, together with the fixed equator, gives the numerical partition.
  The involution `x~-x` also halves the balanced sign-vector count for
  `v=2n>=2`.
- **Repair:** every check now uses an explicit `require` guard that remains
  active under optimization; normal, `-O`, and stored transcripts were
  replayed exactly and the LF script hash was updated.  Equations (1) and (5)
  now state `k>=1` and handle constants separately, the central split is
  called a reflected numerical partition, and `C(2n,n)/2` is recorded as the
  post-gauge balanced class count.
- **Rule:** transcript equality is evidence only when both executions retain
  the predicates being checked.  For Python companions, use explicit raises
  (or another optimization-stable guard) for every truth-bearing assertion,
  and distinguish a numerical reflected decomposition from a literal
  disjoint set partition or quotient count.

## MISTAKE-269 (2026-07-26, THM-2418 carry-cocycle quantifier) -- a base-thirteen prefix formula omitted the integer part of an unrestricted real phase

- **What was written:** candidate THM-2418 introduced an arbitrary real
  `y`, then identified `floor(13^k y)` with the alternating reduction of
  the first `k` digits of the fractional part.
- **Minimal witness:** for `y=1,k=1`, the displayed digit is zero while
  `floor(13y)=13=6 mod 7`.
- **First failed implication:** the digit expansion reconstructs
  `floor(13^k {y})`; for unrestricted `y` it omits
  `13^k floor(y)`, whose residue is `(-1)^k floor(y) mod 7`.
- **Repair:** the theorem now works on the canonical torus phase
  `y in [0,1)`, exactly the domain used by its Haar kernel, physical
  application, and verifier. Equivalently, an unrestricted-real version
  must add the omitted integer-part term. The terminal all-colour
  consequence also retains its necessary rational-profile hypothesis.
- **Rule:** before identifying a radix prefix with a floor, either reduce
  the phase to its canonical fundamental domain or carry the integer-part
  cocycle explicitly.

## MISTAKE-268 (2026-07-26, HYP-1994 twin-center necklace typing) -- the exceptional center 4 was put in a carrier asserted to lie in `6Z`

- **What was written:** HYP-1994 used one symbol `C` both for all twin-prime
  centers, including the center `4` of `(3,5)`, and for a carrier declared to
  satisfy `C subset 6Z`; it then formed `K=C/6`.
- **First failed implication:** `4` is not divisible by six, so `C subset 6Z`
  and the integral typing of `K` already fail at the first member. This is a
  carrier error, not evidence against the finite necklace computation.
- **Strongest survivor:** every nonexceptional twin center is divisible by
  six, and the finite centered-Goldbach packets checked in the hypothesis
  retain their stated exact scope; the exceptional pair is checked
  separately in that range.
- **Repair:** write `C_all=A014574`, `C_6=C_all minus {4}`, and
  `K=C_6/6`. All scaled sumset statements now use `C_6`; statements about the
  full twin-center sequence use `C_all`.
- **Rule:** when an eventually congruence-restricted sequence has startup
  exceptions, type the restricted subcarrier before scaling it. Do not hide
  an exceptional element inside quotient notation.

## MISTAKE-266 (2026-07-26, THM-2403 candidate typing audit) -- a clean-cell projector froze the blocker gate that the target action was moving

- **What was written:** the first proof candidate for THM-2403 restricted
  integration to a fixed clean owner cell `S_j`. That cell already contained
  the unshifted present-blocker safety factor `g_0`, while the proposed target
  table multiplied it by the shifted gate `g_s`. The ensuing clean incidence
  count was then read as a lawful THM-2365 target orbit.
- **First failed implication:** writing `1_(S_j)=f g_0` and translating only
  the displayed blocker factor produces `f g_0 g_s`, not the physical shifted
  layer `f_s g_s`. Scalar-cover redundancy at the untwisted terminal word does
  not make `S_j` invariant under the present-factor target action. The table
  therefore froze or doubled a moving gate and was not the claimed canonical
  orbit.
- **Strongest survivor:** the finite clean-cell incidence inequality remained
  exact and useful as a lower bound. It could not itself type the analytic
  packet.
- **Repair:** THM-2403 now defines the global lawful packet with all six
  unit/guard safeties and all three present blocker safeties shifted together,
  inserts only the target-neutral transported terminal word, and uses the
  clean set solely inside a nonnegative global integral. Scalar cover gives
  the genuine anchored zero globally, while the surviving local incidence
  yields the exact `18 rho_R/13` mass gap and all twelve target colours.
- **Rule:** if a localization predicate mentions a factor moved by the
  action, do not insert it as a fixed projector. Globalize the orbit first and
  restrict only the estimate, unless the entire predicate is transported
  covariantly or separately proved target-neutral.

## MISTAKE-265 (2026-07-26, THM-2396 orbit-localization audit) -- an empty finite orbit universe was used only at global zero mass

- **What was written:** the first promoted THM-2396 used the empty
  relaxed `Z/49Z` word universe only to contradict `delta=0`. It then
  said that coefficient-dependent chamber widths prevented any uniform
  lower bound for the clean-hole mass.
- **First failed implication:** the THM-2394 trichotomy does not require
  global clean-set emptiness. On one fixed high-safe forty-nine-root
  orbit, absence of a clean root supplies (9c) on each of its six
  `q_*`-safe seven-root rows; (9a), (9b), and the THM-2391 top-row
  partition remain pointwise. Thus an individually clean-free orbit
  already embeds into the same empty finite universe.
- **Repair:** for
  `H_0=D_V^c intersection D_(13V)^c`, every generic `z in H_0` has at
  least one clean root among `(z+j)/49`. Since `mu(H_0)=66/91` and
  `integral sum_j 1_S((z+j)/49) dz=49mu(S)`, one gets the uniform
  common-core floor `delta>=66/4459`. Combined with THM-2393, the whole
  last septimal lane has `delta>=1/26754`.
- **Rule:** whenever a finite-orbit contradiction is proved from
  pointwise local constraints, test the orbitwise contrapositive before
  recording only qualitative nonemptiness. Root disintegration can turn
  one forced witness per base orbit into an exact uniform mass floor.

## MISTAKE-264 (2026-07-25, THM-2393 base/root typing audit) -- an undivided high-safe point event was used as a seven-root base event

- **What was written:** THM-2393 Section 3 disintegrated
  `x_r=(y+r)/7` and then wrote `y in X_0`, where
  `X_0=D_(q_*)^c intersection D_(C_3)^c intersection D_(c_3)^c`.
  It said those three high conditions were consequently safe at every
  root `x_r`.
- **First failed implication:** for a speed `v` divisible by seven,
  `x_r in D_v` is equivalent to `y in D_(v/7)`, not to
  `y in D_v`. Invariance among the seven roots does not identify the
  base point with any root. This is the same base-versus-root typing
  hazard isolated in THM-2382.
- **Repair:** define the base event
  `Y_0=D_(q_*/7)^c intersection D_(C_3/7)^c intersection
  D_(c_3/7)^c`. Then `X_0=T_7^(-1)(Y_0)`, so Haar invariance preserves
  the exact mass `396/637`, and `y in Y_0` gives high safety at every
  root. All THM-2393 inequalities and the THM-2394 transversal proof
  are unchanged after this typing repair.
- **Rule:** after a root disintegration, divide every coefficient by
  the root degree when defining the base event. Fibre invariance says
  the root values agree; it does not say that they equal the undivided
  value at the base coordinate.

## MISTAKE-263 (2026-07-25, owner-pivot repair audit) -- a pure primal relation class was treated as its target covector

- **What was assumed:** for the THM-2309 owner pivot, the balanced
  omitted-unit/helper relation
  `h=w_k e_(u_0)-w_(u_0)e_k` has a pure nonzero class in `K/L`.
  A proposed repair table therefore shifted just those two factors and
  called its nonzero colour the corresponding lawful pure target.
- **First failed implication:** target classes are primal elements of
  `K/L`, while lawful target co-shifts are covectors in
  `(K/L)^*=L^perp/<w>`. For normalized
  `g=e_(u_0)-(w_(u_0)/w_k)e_k` and any other ungrafted unit row
  `ell=w_(k')e_(u_0)-w_(u_0)e_(k') in L`,
  one has `g.ell=w_(k')!=0`. Hence `g notin L^perp`, even though its
  primal class is pure. Adding a gauge multiple of `w` cannot help
  because `w.ell=0`.
- **Strongest survivor:** the two-factor table still has two zero axes,
  total mass at least `(11-2L)rho`, and a mixed duplicate-probe
  coefficient at least `rho/2704` for an ordinary failed unit or
  `7rho/24336` for a failed guard. These coefficients are positive but
  are not canonical owner-target currents. The true pure dual is the
  target-blocker/helper dipole `e_a-e_k`, which has different factor
  geometry.
- **Repair:** THM-2384 separates the primal obstruction, true dual
  representative, and surviving analytic corner. It does not claim a
  target landing or row exclusion.
- **Rule:** before turning a relation address into a Fourier shift,
  test it against every row of the chosen relation packet. Membership
  in `K` or purity in `K/L` is not membership in `L^perp`; quotient
  vectors and quotient characters have opposite variance.

## MISTAKE-262 (2026-07-25, post-promotion scope audit of THM-2377) -- a contrapositive retained the distinctness hypothesis it was negating

- **What happened:** THM-2377 proved circulancy from two hypotheses:
  every lower valuation lies strictly below the graft/target layers, and
  the lower valuations are pairwise distinct. Its promoted prose then
  stated the collision contrapositive “under (7),” where (7) still
  included pairwise distinctness, while concluding that two valuations
  coincide.
- **First failed implication:** the proof established
  `hierarchy + distinctness -> circulant`, whose useful contrapositive is
  `hierarchy + noncirculant -> not distinct`. Reusing the conjunction as
  an ambient hypothesis made the displayed conclusion contradictory and
  the statement vacuous, even though the spectral proof was sound.
- **Repair:** the hierarchy is now (7a), distinctness is (7b), and the
  collision gate retains only (7a). Independent hostile checks of the
  translated-factor exclusion, exact overlaps, same-layer counts,
  `chi_7` transfer, Bockstein carries, transcript, and hashes otherwise
  passed unchanged.
- **Rule:** when taking a contrapositive of `A and B -> C`, state the
  retained hypothesis and the negated conclusion separately; never cite
  the original conjunction while asserting `not B`.

## MISTAKE-260 (2026-07-25, post-promotion hostile audit of THM-2362) -- an unordered inverse-root count was promoted to a nonzero-mode sum without fixing the anchor sheet

- **What happened:** THM-2362 correctly separated the target-factor shift
  `d(y+s/13)` from the inverse-root collection `d((y+s)/13)`, but Section
  4 then claimed inverse-root nonzero-mode sums `12rho/13` on danger
  support and `2rho/13` on complement support. The companion checked only
  the arithmetic tautologies `rho-rho/13` and `rho-11rho/13`; it never
  checked the actual `s=0` root sheet.
- **First failed implication:** the count
  `sum_s d((y+s)/13)=2-d(y)` determines only the zero-character average.
  Fourier inversion gives `sum_(k!=0) Rtilde(k)=R_0-Rtilde(0)`, and
  `R_0=int w(y)d(y/13)dy` depends on a chosen lift/section which the
  unordered count does not contain. Minimal hostiles are `y=99/100`,
  where danger count is one but the chosen anchor is zero and the
  nonzero-mode sum is `-1/13`, and `y=1/2`, where complement count is
  eleven but the anchor is zero and the sum is `-11/13`.
- **Strongest survivor:** THM-2362's target-shift identities, successor
  statistic, all danger/complement mode floors in Sections 2--3,
  pure-word redundancy, and fork boundary are unchanged. THM-2364 uses
  only those target-shift identities and is unaffected.
- **Repair:** Section 4 now retains only the almost-everywhere
  inverse-root count and zero-character average, names the missing
  measurable section/anchor hypothesis, and records both negative
  hostiles. The exact companion now checks the actual anchor entries.
- **Rule:** a multiset or unordered fibre count determines the zero
  character only. Before inferring nonzero Fourier sums, explicitly
  identify the inversion anchor (`s=0`), prove that it is canonical or
  supply a section, and make the executable measure that anchor rather
  than a rearranged arithmetic consequence.

## MISTAKE-261 (2026-07-25, planar-graph tomography audit) -- a refined joint target was conflated with the coarse THM-2334 target

- **What was written:** the first THM-2356 candidate correctly proved that
  planar chirp intensities reconstruct every off-diagonal Gram entry inside
  a retained graph
  `Z_c(q)=A(q,q^2/2+c)`. It then claimed that a graph supported at two
  targets forces a nonzero THM-2334 target fibre and called the remaining
  obstruction only the one-sparse graph-location boundary. Commit
  `e0839fe1e` strengthened this to a vertical-tensor classification, and
  `59c933aae` promoted it after an independent audit; both retained the
  same refined-to-coarse category error. The audit's planar-tomography and
  refined-support checks remain valid, but its coarse interpretation is
  retracted.
- **First failure:** a nonzero refined coefficient `A(q,z)` with `q!=0`
  does not imply that the coarser target coefficient
  `C(q)=sum_z A(q,z)=sum_c Z_c(q)` is nonzero. The first forgotten
  operation is summation over graph/jet labels, and distinct graph
  coefficients at the same target can cancel. Over `F_169`, the three
  entries
  `A(0,0)=1`, `A(a,a^2/2)=1`,
  `A(a,a^2/2+c)=-1` for `a,c!=0` give a two-supported graph with a
  nonzero off-diagonal Gram entry but `C=delta_0`. More strongly, changing
  one whole graph by a global sign preserves every graph chirp intensity
  and labelled singleton energy. The two-graph pairs
  `(delta_a,delta_0-delta_a)` and
  `(delta_a,-delta_0+delta_a)` also preserve the total scalar current while
  changing the coarse target Dirichlet energy from zero to positive.
  More structurally, the exact coarse zero-only space has dimension
  `169+168^2=28,393`; vertical tensors occupy only `169` dimensions.
- **Repair:** THM-2356 now retains its planar tomography theorem and the
  valid conclusion that support two forces a nonzero **joint** target/jet
  fibre with nonzero target coordinate. It explicitly records one
  independent `U(1)` phase per active graph, the exact row-sum kernel, the
  missing cross-graph terms in the coarse Dirichlet energy, both exact
  hostiles, and the repaired residual. Averaging target Dirichlet energy over first-jet characters
  does recover a positive edge defect whenever some refined `q!=0` fibre
  survives, but only for an exact factor-coloured THM-2337 probe; a
  nontrivial jet character is not the original THM-2334 coarse current,
  and the defect need not occur for the trivial jet character. Reaching
  that current therefore still requires cross-graph phase transport, jet
  depolarization, a target-preserving physical realization, a target row
  supported on one graph, a familywise cone, or another
  coefficient-sensitive sidecar. A second hostile audit accepted this
  repaired statement.
- **Rule:** whenever a nonlinear decomposition refines a target quotient,
  write the map back to the original target and test its fibrewise sum.
  Complete reconstruction inside every refined fibre does not compare
  phases between fibres.

## MISTAKE-259 (2026-07-25, concurrent namespace and rebase audit) -- THM-2356 crossed and a second conflict block escaped inspection

- **What happened:** a session fetched and verified that `THM-2356` was
  free, then committed
  `THM-2356-degree-eighteen-perfect-quartic-wall-closure.md`. During the
  interval before push, the earlier
  `THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing.md`
  reservation landed. The later reservation rebased across it without a
  post-rebase namespace search. In the same rebase, a status-word conflict
  in `THM-2347` was resolved after inspecting only the frontmatter and first
  visible marker block; a second conflict block in the audit appendix was
  accidentally committed.
- **First failure:** the session repeated MISTAKE-258's forbidden combined
  rebase-and-push pattern. It also treated one resolved conflict hunk as
  evidence that the whole file was marker-free instead of running a global
  marker search and `agents/check_docs.py` before push.
- **Repair:** the earlier chirp reservation keeps `THM-2356`; the still-empty
  perfect-quartic stub moved coherently to the next free namespace
  `THM-2359`. The two independently valid `THM-2347` audit summaries were
  merged into one marker-free appendix without changing its proved claim.
- **Rule:** after every reservation rebase, search the rebased tree again for
  the exact YAML ID and filename before a separate push. After every manual
  conflict resolution, search every edited file for all three conflict
  markers and run the documentation gate before committing.

## MISTAKE-258 (2026-07-25, concurrent namespace audit) -- THM-2350 crossed during the reservation rebase

- **What happened:** after fetching `origin/main` at `6b23df8f4`, one
  session verified that `THM-2350` was free and committed the q-adic
  prefix-residue collision-spectrum reservation. During that interval,
  earlier commit `4890b8d25` landed the distinct
  `THM-2350-owner-pivot-dual-dipole-normal-form.md` reservation. Rebasing
  the later commit exposed the collision, but the combined pull/push command
  pushed it before a post-rebase YAML/filename search was run.
- **First failure:** the pre-reservation search was correct but stale after
  the rebase. A successful automatic rebase does not validate a scarce ID,
  and combining rebase with push removes the manual collision gate required
  by MISTAKE-257.
- **Repair:** the earlier owner-pivot reservation keeps `THM-2350`. The
  later carry-spectrum stub first moved to `THM-2351`, but before that repair
  could push, earlier commit `341e689d9` reserved `THM-2351` for the distinct
  self-target ANOVA ledger. The carry stub therefore moved again to
  `THM-2352`, still before receiving proof content. None of the reservations
  had entered the proved dependency graph.
- **Rule:** never combine the final reservation rebase and push. Rebase,
  search the resulting tree again for the exact YAML ID and filename, and
  only then push.

## MISTAKE-257 (2026-07-25, concurrent namespace audit) -- THM-2333 crossed between fetch and push

- **What happened:** one session fetched `origin/main`, verified that no
  `THM-2333` filename or YAML ID was present in the fetched tree, and reserved
  `THM-2333-relation-residue-current-and-character-twist-pushforward.md`.
  During the short interval before its push, earlier commit `d40de3b406`
  landed the distinct
  `THM-2333-abel-target-fibre-sum-landing-and-zero-fibre-boundary.md`
  reservation. Rebasing the later reservation therefore created a duplicate
  ID on `origin/main`.
- **First failure:** a clean fetch-and-search reduces but does not eliminate
  the final compare-and-swap race between concurrent reservations. A
  successful rebase is not sufficient: the rebased tree must be searched
  again for both the exact YAML ID and filename before pushing.
- **Repair:** the earlier Abel target-fibre reservation keeps `THM-2333`. The
  later relation-residue pushforward stub moved coherently to `THM-2334`
  before receiving proof content. Neither stub is in the proved dependency
  graph.
- **Rule:** after every reservation rebase, rerun the namespace search on the
  actual rebased tree. If a collision appeared, preserve the earlier remote
  owner and renumber or delete the later reservation before push.

## MISTAKE-256 (2026-07-25, concurrent namespace audit) -- a fetched but unreplayed THM-2325 reservation was missed

- **What happened:** commit `ba9d0cd284` reserved
  `THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank.md`.
  A second session fetched that commit while it was still one commit behind,
  searched only its checked-out files plus commit subjects, and saw no
  `THM-2325` string because the reservation commit subject did not contain the
  numeric ID. It then rebased but failed to repeat the filename/YAML search
  before reserving the distinct cyclotomic-chromatic candidate as THM-2325 in
  later commit `dced8e4253`.
- **First failure:** a remote fetch is not a namespace check until the fetched
  tree has either been inspected directly or replayed into the worktree. A
  subject-only history grep is not a substitute for searching YAML IDs and
  filenames in the fetched commit.
- **Repair:** the earlier prescribed-target-gain theorem keeps THM-2325. The
  later cyclotomic-chromatic file was deleted rather than renumbered because
  its narrower divisor-splitting argument was independently superseded before
  promotion by THM-2323's nonnegative-sector strengthening. No proved
  dependency ever cited the duplicate file.
- **Rule:** after the final fetch/rebase immediately preceding a reservation,
  rerun all three checks on the resulting tree: filename, exact YAML ID, and
  remote history. If inspecting without rebasing, search `origin/main^{tree}`
  directly; do not rely on the current worktree or commit subject.
- **Same-day recurrence:** two sessions independently found the same
  two-colour marked-edge triangle and crossed in flight at `THM-2327`.
  Commit `681c74c8b` reserved
  `THM-2327-two-colour-marked-unit-c3-triangle.md` first; the later
  `211f569a9` reservation used a different slug but the identical target.
  The later stub was deleted rather than renumbered, because there is only
  one mathematical theorem to promote. The proof and audit were already
  broadcast to the earlier owner; no proved dependency cited either stub.

---

## MISTAKE-255 (2026-07-25, quantitative handoff-gluing audit) -- a positive-value floor was used in the wrong direction

- **What was written:** for the arrival density
  `g=P_beta^k 1_F`, every positive value is at least `beta^(-k)`;
  the first THM-2310 candidate then cited this floor as the reason that
  `measure({g>0})>=measure(F)`.
- **First failure:** a positive lower bound on `g` gives an upper bound on
  the measure of its support from its integral, not the claimed lower bound.
  The displayed conclusion was true, but that implication did not prove it.
- **Repair:** use the other pointwise bound, `0<=g<=1`, together with
  `integral g=measure(F)`. Then
  `measure(F)=integral_{g>0}g<=measure({g>0})`. Retain
  `g>=beta^(-k)` on `{g>0}` only for the later lower bound on the mass of a
  lifted overlap. The independent audit caught the error before THM-2310 was
  promoted out of candidate status.
- **Rule:** for a nonnegative density, an upper pointwise bound controls
  support measure from below; a positive lower pointwise bound controls it
  from above. Record which direction is being used before transporting a
  mass statement through a Perron image.

---

## MISTAKE-254 (2026-07-25, concurrent namespace audit) -- THM-2304 was reserved twice in flight

- **What happened:** one session pushed the canonical blocker-word handoff
  hypergraph reservation as `THM-2304`. Before seeing that remote commit, a
  concurrent session pushed the deepest-boundary cyclotomic-current
  reservation under the same YAML ID. The filenames were distinct and
  neither empty stub had become a proved dependency, but the bare identifier
  was temporarily ambiguous on `origin/main`.
- **First failure:** checking the next ID before the first reservation did not
  prevent a second reservation made from a stale remote view. The collision
  was detected only after both in-flight commits landed.
- **Repair:** the deepest-boundary cyclotomic-current reservation keeps
  `THM-2304`. The earlier handoff-hypergraph stub is deleted and its theorem,
  companion, transcript, hashes, and reproduction commands move coherently to
  `THM-2305` before promotion. No mathematical claim or proved dependency
  changes.
- **Rule:** an ID check is a compare-and-swap operation, not a one-time read.
  Fetch and search YAML IDs immediately before the reservation push, then
  fetch and verify uniqueness immediately afterward. During a collision
  repair, cite the slug and keep both files unproved until uniqueness is
  restored.
- **Same-day recurrence:** the dual rank-six reconstruction stub and the
  mirror-double stable-mixture stub crossed in flight at `THM-2307` despite
  both owners checking a then-current remote. The earlier dual-rank
  reservation keeps `THM-2307`; the later mirror-double reservation moves
  immediately and coherently to `THM-2308`. Both were still unproved stubs,
  so no mathematical statement or proved dependency changed.
- **Second same-day recurrence:** the quantitative handoff-gluing reservation
  at `09:47:56` and the sparse-root bispectrum reservation at `09:49:32`
  crossed in flight at `THM-2310`. The earlier handoff theorem keeps
  `THM-2310`; the later unproved bispectrum stub moved to `THM-2312` before
  promotion.
- **Third same-day recurrence:** the biprime collision-frontier reservation
  at `10:04:29` and the marked target-gain corolla reservation at `10:08:31`
  crossed in flight at `THM-2313`. The earlier biprime reservation keeps
  `THM-2313`; the later corolla file, companion names, transcript name, and
  reproduction commands move coherently to `THM-2315` before promotion.
  The corolla was still unproved when the collision was detected, so no
  proved dependency changed.

---

## MISTAKE-253 (2026-07-25, concurrent namespace audit) -- two theorem identifiers were allocated twice

- **What happened:** concurrent sessions promoted both the expiration
  owner-absorber cut and the mixed scalar relation theorem as `THM-2271`.
  They also allocated `THM-2272` to both the persistent knot-interaction
  theorem and a later empty common-anchor reservation. Filename slugs kept
  the mathematical statements distinguishable, but bare-ID citations became
  ambiguous.
- **First failure:** reservations checked their intended filenames but did
  not re-fetch and check the full YAML-ID namespace immediately before the
  reservation commit. Concurrently allocated IDs then crossed in flight.
- **Repair:** preserve the earlier allocations by creation time:
  `THM-2271` remains the proved expiration owner-absorber cut and
  `THM-2272` remains the proved persistent knot-interaction theorem.
  Concurrent reservations had already assigned `THM-2274` to the mixed-scalar
  relative-rank candidate and `THM-2276` to the shallow-owner residue-aligned
  crossing candidate. The later proved mixed-scalar theorem, its companion,
  and its transcript therefore move coherently to `THM-2275`; the later
  unproved common-anchor stub moves to `THM-2277`. An intermediate repair
  that sent the stub to already-occupied `THM-2276` was itself corrected
  before either file became a proved dependency. All hashes and internal
  reproduction commands were regenerated.
- **Rule:** before reserving, fetch and check filename, frontmatter `id`,
  current remote history, and the next several candidate IDs. Cite every
  theorem by ID plus slug; a bare number never resolves historical collisions.

---

## MISTAKE-252 (2026-07-25, integrity repair for THM-2243) -- correct mathematics was paired with stale frontmatter hashes

- **What was recorded:** THM-2243's frontmatter claimed script hash
  `756b85f...` and output hash `b7cc692...` on working-tree LF bytes.
- **Why it was wrong:** direct `sha256sum` on the tracked companion and
  transcript gives respectively `81e51c5...` and `6b44e1c...`. The stale
  metadata made a faithful replay look like artifact drift even though the
  computation itself remained reproducible.
- **Correct framing:** the frontmatter now records the actual full hashes.
  Normal and optimized-mode runs are byte-identical to one another and to
  the stored transcript, including the hostile `(3,4,5)` value
  `8907541/62377224`. This is an integrity-metadata repair, not a change to
  THM-2243's statement, proof, constants, or status.

---

## MISTAKE-250 (2026-07-25, correction of the post-THM-2250 all-equal hostile controls) -- a 169-step kernel was assigned the wrong 13-step stochastic square root

- **What was claimed:** a session-level continuation of THM-2250 used a
  three-state unit-13 chain to evaluate the all-equal `(3,4,5)` whole-clause
  cylinder. It reported `354710/2599051<961/6930` for `u=H` and
  `1057220/7797153<961/6930` for `u=6H`, apparently closing both hostile
  relations. These values were not promoted into a theorem.
- **Why it is wrong:** the calculation began from the known unit-169 kernel
  and chose the parity-twisted square root `JP`. The exact unit-13 transition
  is `P=(-I+14Pi)/13`. Both `P` and `JP` have the same square, but `JP`
  reverses the one-step eigen-sign of the antisymmetric annular observable
  `1_A-1_B`. It also fails the known binary danger lump:
  its `A->A` probability is `2/13`, whereas exact thirteen-root counting
  gives `1/13`. Squaring destroyed precisely this midpoint `Z/2` sign.
- **Correct framing:** THM-2254 derives `P` pointwise from all thirteen roots
  and records the exact midpoint bridge. The corrected bounds are
  `28460/199927>961/6930` for `u=H` and
  `734515/5198102>961/6930` for `u=6H`. Thus neither relation is closed by
  this cylinder. THM-2250's reduction of `(3,4,5)` to equal normalized cores
  is unaffected. **Subsequent status:** THM-2257 closes the all-equal branch
  by a different exact `169`-image sieve; the parity-alias correction remains
  valid, while the `165` first-depth-one scalar rows and LRC(14) remain open.

---

## MISTAKE-251 (2026-07-25, correction of the depth-one odd-transfer reflection) -- the 165-profile family was split as 151+14 instead of 150+15

- **What was written:** the hostile Bellman classification in
  `07-reflections/lrc14-depth-one-odd-transfer-cone-and-private-stratum-no-go-codex-20260725.md`
  said that `151` profiles have exactly one depth-one blocker and `14`
  profiles have two.
- **Why it is wrong:** the current family is
  `(1,b,c)` with `5<=c<=19` and `1<=b<c`. The repeated-first-depth subfamily
  is `(1,1,c)` for all fifteen values `c=5,...,19`; the remaining sum is
  `3+4+...+17=150`.
- **Correct framing:** the exact split is `150+15=165`. The two hostile
  Bellman constructions and every bound in the reflection are unchanged;
  only the class counts were wrong. THM-2255 uses and freezes the corrected
  split.

---

## MISTAKE-249 (2026-07-24, correction of THM-2176) -- right continuation was called a two-sided congruence for an arbitrary monoid

- **What was claimed:** the universal continuation theorem began with an
  arbitrary monoid and called equality of the right-continuation profiles
  `z |-> f(x+z)` the largest monoid congruence contained in `ker(f)`.
- **Why it is wrong:** the displayed proof establishes stability under right
  multiplication. In a noncommutative monoid it gives the largest right
  congruence, not necessarily a two-sided congruence.
- **Correct framing:** THM-2176 now states the theorem for a commutative
  monoid, which is exactly the connected-sum setting used there. In the
  noncommutative generalization one must either retain both-sided contexts or
  say “largest right congruence.” The min-plus kernel theorem is unaffected.

---

## MISTAKE-248 (2026-07-24, correction of THM-2158) -- one quartic pole congruence was presented as two independent conditions

- **What was claimed:** polynomiality of the canonical quartic approximate
  root was said to require two independent conditions,
  `V|beta` and `V^3|(4 gamma V^2-beta^2)`, and the examples were described
  as showing both directions of independence.
- **Why it is wrong:** the second congruence implies the first in the UFD
  `C[x]`. Modulo `V^2` it gives `V^2|beta^2`; irreducible valuations then
  give `V|beta`. There is no example satisfying the cubic congruence while
  failing linear-coefficient regularity.
- **Correct framing:** THM-2158 now uses the single condition
  `V^3|(4 gamma V^2-beta^2)`. Its failure defines one effective pole
  divisor. Whenever the constant coefficient has a pole, it is strictly
  deeper than the linear coefficient, so polynomial fibre translations
  cannot cancel it. The divisor, not an independent pair of jet flags, is
  the faithful remaining quartic carrier.

---

## MISTAKE-247 (2026-07-24, audit of THM-1795) -- generic recurrence order was promoted to a uniform nullcone cutoff

- **What was claimed:** four computed cells and a creative-telescoping
  heuristic were promoted to the exact formula
  `D(M,N,d)=(M+N)(2d+1)`, and this recurrence order was declared to be the
  number of moments cutting out every bounded coefficient stratum.
- **Why it is wrong:** the general `s`-degree bound, telescoper order, and
  minimality were not proved. Even a valid P-recurrence propagates initial
  zeros only where its leading coefficient is nonzero; that coefficient may
  vanish at small indices or on parameter loci. A generic specialized order
  therefore does not by itself give a uniform algebraic detection depth.
- **Correct framing:** retain the formula as a conjecture verified on four
  generic cells. To turn it into a cutoff, construct a uniform symbolic
  recurrence, prove its order and leading coefficient, control every
  exceptional locus, and then prove zero propagation on the required initial
  interval. THM-1740 supplies bounded-stratum finiteness without this explicit
  size, while THM-1790 supplies the rigorous unbounded linear lower scale.

---

## MISTAKE-246 (2026-07-24, audit and repair of THM-1790) -- affine equation counting was mistaken for projective nonemptiness

- **What was claimed:** the old proof of THM-1790 parameterized a monic
  degree-`d` polynomial by `d` affine coefficients and asserted that the `d`
  equations `L(B^j)=0`, `1<=j<=d`, must have a solution. It then advertised
  exact survival depth `d` for every exact degree.
- **Why it is wrong:** `d` equations in `d` affine variables may cut out the
  empty set or only a forbidden point; equation counting alone supplies no
  solution. The exact computations through `d=4` did not prove the all-degree
  statement.
- **Correct framing:** work in the `(d+1)`-dimensional vector space
  `C[s]_{\le d}`. The `d` moment equations are homogeneous. If their only
  common affine zero were the origin, the radical of their ideal would have
  height `d+1`, contradicting Krull's height theorem. Hence a nonzero
  degree-at-most-`d` witness exists. This proves a degree-cap detection lower
  bound, not an exact-degree or monic witness and not exact survival `d`.
  Applying the same argument in `s*C[s]_{\le d+1}` gives the repaired
  two-charge lower bound in THM-1790. A matching effective upper bound is the
  one-variable boundary of the Strong Factorial Conjecture and remains open.

---

## MISTAKE-245 (2026-07-24, audit of mac-mini S171/P14) -- a bounded local census was promoted to the global second OPEN-Q-108 measure

- **What was claimed:** an exhaustive census on twelve-subsets of
  `{1,...,18}`, supplemented by random rows at larger height, found
  `426/35035` (the AP drop-12 value) as the second measure above the drop-6
  minimum `7/858`; this was called a global two-level stability law.
- **Why it is wrong:** the structured one-swap row
  `{1,2,3,4,5,7,8,9,11,12,13,20}` lies just outside the exhaustive bank.
  Two independent exact interval algorithms give
  `3859/420420=7/858+1/980<426/35035`. Random sampling did not visit the
  thin structured locus. This is the same failure mode as reading a random
  search away from a dilation or relation subvariety as an extremal census.
- **Correct framing:** THM-2162 proves that the drop-6 row is nevertheless a
  strict local minimizer under every one-element replacement. Its exact local
  neighbor floor is `3859/420420`, uniquely at replacement `10->20`.
  No global second-value classification is currently proved.

---

## MISTAKE-244 (2026-07-22, codex audit of HYP-8950/S109 and HYP-8955/S110) -- a Newton-face Wronskian was called DvdK, and a linear edge polynomial was confused with a weighted-linear coordinate

- **What was claimed:** the resonant face equation `1 in im J(F,-)` was called
  a one-variable DvdK constant-term problem; `xy` was called nonprimitive, and
  a linear collapsed edge polynomial was treated as a weighted-linear
  coordinate and a terminating descent.
- **Why it is wrong:** at terminal weight the bracket is the finite Laurent
  Wronskian `[x]F[y]G-[y]F[x]G`, with no supplied DvdK map. The edge polynomial
  of `x^3+y^2` is linear at weights `(2,3)`, but its degree is `6>5` and it
  has no mate; `xy` is primitive. Coordinates `x+(y+x^2)^m` also have
  repeated proper-power top faces.
- **Correct framing:** THM-2102 proves power-free top-face descent directly.
  Proper-power faces obey a first-defect Hamiltonian equation whose live data
  are the approximate-root quotient and resonant target-shear class. S109's
  Hamiltonian-cokernel viewpoint survives; no GMC/DvdK-to-JC implication is
  known.

---

## MISTAKE-243 (2026-07-22, codex dependency-aware Lean audit of THM-2101) -- an unbuilt additive-residue source was called kernel-checked

- **What was claimed:** the additive DvdK checkpoint described Check A, the
  additive orbit contradiction, and the full-root Lagrange identity as
  kernel-checked and printed only the standard Mathlib axioms.
- **Why it is wrong:** `lake build TournamentH7.GMC2LaurentShiftCheckA`
  reaches the source only after building its missing dependencies, then fails
  on Laurent coefficient application and ambiguous `C`, an unsolved `nsmul`
  commutation goal, and a failed Lagrange rewrite. The subsequent `#print
  axioms` lines include `sorryAx`; a direct `lake env lean` had stopped earlier
  at a missing dependency `.olean` and therefore certified nothing.
- **Correct framing:** the original formal claim required a dependency-aware
  rebuild and a genuine analytic bridge. The core was subsequently repaired
  and root-imported. THM-2101 is proved on paper by three product-free routes:
  additive monodromy, transcendental specialization, and a purely t-adic
  Newton-packet proof. Their global wrappers remain unformalized in Lean.
  Later `GMC2DvdKOmegaWiring.singlePolyCrux_holds` closed the independent
  small-root route, so `GMC2Main.gmc2` is now unconditional; this does not
  retroactively validate the failed build described above.

---

## MISTAKE-242 (2026-07-22, codex referee of THM-2098) -- a pure-transverse collision budget was transported to mixed lanes and an unguarded depth-zero core

- **What was claimed:** the first pushed THM-2098 consumer said every live
  rank-eight-through-eleven coefficient-plane cover fell into either seven
  transverse bands governed by the exact budget `5(n-7)/49`, or seven
  guard-proportional bands. It treated the high vertical cover pointwise.
- **Why it is wrong:** the budget uses both that every band has mass `5/49`
  and that the transverse family itself covers the guard complement. Neither
  holds in a low mixed row such as `(n,r,t)=(8,1,7)`. Live nontrivial towers
  have guarded terminal sizes only `7..10`; the depth-zero size-eleven core
  has no preceding guard. Finally, containment gives an almost-everywhere
  two-torus cover, so a single uncovered boundary fiber proves nothing.
- **Correct framing:** for a rank-two character image, `r=0` has the exact
  collision/tree budget; `1<=r<=n-7` is only count-isolated; and `r>=7`
  forces the vertical bands to cover the one-dimensional guard complement
  almost everywhere, by a positive-measure leftover plus Fubini argument.
  This applies to live guarded sizes `8,9,10`; rank-one images keep their
  freeze route and the depth-zero size-eleven core remains outside the split.

---

## MISTAKE-241 (2026-07-22, codex semantic audit of HYP-8935/S106) -- floating root asymptotics were promoted to a formal-series/Hensel reduction

- **What was claimed:** HYP-8935 described THM-2067 as four
  “Mathlib-ready” pieces plus one unramified-Hensel gap, called an elementary
  formal-log derivation and unramified descent verified, and inherited
  HYP-8931's unique-face bypass as a working alternate consumer.
- **Why it is wrong:** the NumPy checks compare finite floating evaluations.
  Check B reports the small-root product approximately equal to its leading
  term even for examples whose displayed higher moments are nonzero; a ratio
  tending to one cannot exclude fractional powers or identify a formal series.
  No proof selects the Hensel roots inside the rational splitting field,
  descends through the roots of unity, or establishes the required local/global
  subset product. HYP-8931's consumer is independently void by MISTAKE-240.
- **Correct framing:** HYP-8935 is an open dependency map. The later
  `GMC2OrbitProduct.lean` kernel-checks the abstract transitive orbit-product
  identity and fixed-product valuation-zero lemma, and
  `GMC2PhiIrreducible.lean` now checks irreducibility. The claim audited here
  remains invalid, but later exact Weierstrass, frame, transpose, and Omega
  modules independently close the small-root route and prove
  `GMC2Main.gmc2`. THM-2101's additive wrappers remain optional Lean work.

---

## MISTAKE-240 (2026-07-22, codex semantic audit of HYP-8931/S230) -- the lowest-face unique-channel bypass has an inconsistent class predicate

- **What was claimed:** HYP-8931 called
  `exists_nonzero_lowest_face_seed_of_uniqueChannel` a nonvacuous replacement
  for the DvdK seed on an asserted `98/116` or “84%” unique-channel class, and
  described the remaining work as a one-line connection to NC2.
- **Why it is wrong:** `LowestFaceUniqueChannel P` quantifies over every
  `lambda`, `delta`, and finite set satisfying only an exact-level-set
  equivalence. Take `lambda=0`, `delta=-1`, and `F=empty`. Radial
  exponents are nonnegative, so this is a valid empty level set for every
  `P`. The predicate then demands a positive-mass balanced composition on
  the empty subtype, which cannot exist. Thus no polynomial satisfies the
  premise. Lean correctly checks the implication, but the theorem is vacuous.
  The `98/116` number is only a bounded support census, not a theorem or a
  universal proportion.
- **Correct framing:** HYP-8930's fixed-support theorem
  `ct_ne_zero_of_unique_balanced` and its `dvdk1_of_uniqueChannel`
  corollary are substantive: a specified unique balanced composition prevents
  coefficient cancellation. HYP-8931 is quarantined until its class predicate
  is restricted to a genuine nonempty/straddling lower face and a resulting
  seed is wired into the NC2 descent. Independently, the later root-imported
  `GMC2Main.gmc2` closes GMC(2) unconditionally through
  `GMC2DvdKOmegaWiring.singlePolyCrux_holds`; HYP-8931 is neither used nor
  rehabilitated by that closure.

---

## MISTAKE-239 (2026-07-22, codex THM-2080 direction audit) -- terminal guard containment was reversed into a cover of the guard instead of its complement

This correction initially landed under the already-used `MISTAKE-231`; 239 is
the canonical ID. MISTAKE-231 remains the observable-relative entropy
correction.

- **What was claimed:** the first pushed version of the mixed-radius fold
  argument said that `G_Q subset E_h` implies
  `E_h subset union_(q in Q) D_q`. It derived a low reduced-ratio invoice for
  rank six and described `ab<=36` as the remaining projective conductor.
- **Why it is wrong:** points of `G_Q` are allowed to lie inside `E_h`; that is
  exactly what the containment permits. The implication runs on the other
  side. A point outside `E_h` that lies in no `D_q` would be in `G_Q` and
  violate containment, so the correct statement is

  ```text
  E_h^c subset union_(q in Q) D_q.
  ```

  The exact fold formula itself was correct; only its set-cover consumer was
  reversed.
- **Correct framing:** THM-2080 proves the mixed overlap floor
  `measure(D_q intersect E_h)>=1/42` for odd `h`, with equality only at
  `q=6h`. Hence each danger comb has outside-guard capacity at most `5/42`.
  Six distinct combs have total capacity strictly below `5/7`, so they cannot
  cover `E_h^c`, whose measure is `5/7`. The corrected argument eliminates
  terminal rank six altogether and sharpens tower depth to four; there is no
  surviving `ab<=36` lane.

---

## MISTAKE-238 (2026-07-21, codex audit of HYP-8920/S227) -- the empty safe set of the full dyadic counterexample was transported through a homeomorphism that starts only at its nonempty quotient core

This correction initially landed under the already-used `MISTAKE-230`; 238 is
the canonical ID. MISTAKE-230 remains the inverse-form-class outcome-counting
correction.

- **What was claimed:** HYP-8920 combined mirror parity with THM-2075 and
  asserted that a strict dyadic-seam counterexample `S=2C union {x,y}` has
  `chi(G_S)=0`, that this zero Euler characteristic descends through the
  dyadic tower, and hence that the hereditarily primitive terminal core has
  `chi(G_(Q_r))=0`. It then called the dyadic wall equivalent to exclusion of
  such zero-component terminals.
- **Why it is wrong:** THM-2075's homeomorphisms are

  ```text
  G_C=G_(Q_0) -> G_(Q_1) -> ... -> G_(Q_r).
  ```

  They do **not** start at `G_S`. At the inner guard levels, exactly one safe
  child survives and doubling is bijective. At the final outer step from
  `C` to `S`, the two original odd tails kill both lifts over every point of
  `G_C`; there are zero safe children. Thus `G_S=empty` while `G_C` is
  nonempty, and no homeomorphism connects them. Transporting `chi(G_S)=0`
  across THM-2075 is a domain error.
- **Correct framing:** THM-2075 gives

  ```text
  chi(G_C)=chi(G_(Q_r))>0.
  ```

  The strict dyadic obstruction is instead THM-2061's folded containment of
  the nonempty `G_C` in the two-tail danger locus. Mirror parity remains
  valid: every divisor-complete terminal contains an even speed, reversal
  acts freely on its nonempty safe set, and therefore its component count is
  even (hence at least two). That fact does not stop the original tails from
  killing the two outer lifts. The live carrier remains the terminal
  component word, its THM-2075 dyadic addresses, and the original tail-owner
  constraints.
- **Unaffected work:** the two-charge DvdK Lean theorems in the same S227
  session are independent of this LRC level mismatch and remain valid.

---

## MISTAKE-237 (2026-07-21, codex audit of HYP-8905/S225/S103) -- a valid binary symmetric-Hessian subcase and three analogous descent programs were promoted to an NC2-to-JC bridge and equivalent formulations of JC(2)

This correction initially landed under the already-used `MISTAKE-229`; 237 is
the canonical ID. MISTAKE-229 remains the tangent-disk/Heegner correction.

- **What was claimed:** S103 wrote `NC2 => GMC2 => ... => JC(2)`, identified
  two-dimensional Zhao/Laplacian vanishing with the Gaussian nullcone, and
  called the rank-one/rank-two nilpotent-matrix divide the true/false Jacobian
  boundary. S225 called VC(4), leading-form descent, and a proposed
  Lame-for-polygons bound three equivalent forms of the JC(2) obstruction.
- **Why it is wrong:** no arrow from GMC(2), let alone NC2, to JC(2) is proved
  in THM-1830 or elsewhere in the repository. The symmetric reduction relevant
  to a general planar Keller map increases dimension, landing on the
  four-variable symmetric/Vanishing-Conjecture problem; solving the symmetric
  problem in two variables does not solve the image of that reduction.
  Laplacian iterates and complex Gaussian moments are different functionals on
  different spaces until a predicate-preserving transform is supplied.
  Finally, the existence of rank-two nilpotent matrices in dimension three
  neither constructs nor classifies Keller collisions. The leading-form,
  Jelonek, VC(4), and continued-fraction programs may share a termination
  motif, but no equivalences among those programs were proved.
- **Correct framing:** retain the exact binary homogeneous calculation for
  `d>=2`: if `P=A z^d+B zbar^d` is harmonic, then
  `det Hess(P)=-4d^2(d-1)^2AB|z|^(2d-4)`, so nilpotence forces one side and
  THM-2063 makes the resulting gradient map tame. This is a classified
  two-dimensional symmetric subcase. Treat VC(4), planar leading forms,
  inverse-Jelonek data, and Lame-for-polygons as separate live routes. Any
  transfer among them must state its map, dimension, predicate, and loss.

---

## MISTAKE-236 (2026-07-21, codex audit of THM-1330) -- necessary Keller-monoid anatomy was titled and described as a classification of all Jacobian counterexamples

This correction initially landed under the already-used `MISTAKE-228`; 236 is
the canonical ID. MISTAKE-228 remains the Paley/LRC-role correction.

- **What was claimed:** THM-1330 was titled "the exact picture of the set of
  all counterexamples" and called the Jelonek hypersurface, complement cover,
  monodromy, and degeneration package a "covering-triple classification."
- **Why it is wrong:** the same theorem's honest ledger says that the inverse
  Jelonek realization problem, realizable generic degrees, irreducible seeds,
  unique factorization, and `JC(2)` are open. A package recorded by every map
  is a necessary invariant, not a classification, unless one proves both
  realization and completeness up to the stated equivalence. The phrase
  "subject to polynomial-with-constant-Jacobian gluing" simply puts the hard
  problem back into the definition. Monoid ideality and degree factorization
  likewise organize unknown seeds without enumerating them.
- **Correct framing:** THM-1330 is a necessary-structure atlas: Keller
  nonunits have function-field degree at least three, form a two-sided ideal,
  factor finitely under composition, are etale and nonproper, and impose
  covering/monodromy restrictions on the Jelonek set. The verified explicit
  families occupy cells of this atlas. A complete classification would still
  have to solve realization and equivalence for the covering data and classify
  irreducible seeds. THM-2063 legitimately classifies one empty planar stratum;
  it does not turn the atlas into a global classification.

---

## MISTAKE-235 (2026-07-21, audit of S102 / HYP-8879) -- a shared kernel-sum syntax and four truncated examples were promoted to an LRC-to-GMC equivalence and an AP-core reduction

This correction subsumes the narrower strict-measure repair that initially
landed under a colliding `MISTAKE-227` in commit `c270f8aaa`. Canonical
`MISTAKE-227` remains the AP-chain saturation error below; the modular-cusp
bridge is `MISTAKE-233`, and additive-energy/augmentation is `MISTAKE-226`.

- **What was claimed:** S102 wrote the Fourier expansion of the LRC lonely-set
  measure as a sum over the integer relation lattice, compared it with a
  balanced-channel expansion of a Gaussian moment, and called them the same
  noncancellation problem. It identified the zero-frequency term `(6/7)^13`
  with THM-878's clock floor, then used four supports of sizes four and five
  and a frequency cutoff `|k_i|<=9` to announce that LRC(14) reduces to
  maximal-resonance or AP-neighborhood cores.
- **Why it is wrong:** “sum over an integer kernel” is a schema, not an
  equality of typed objects. The LRC sum ranges over all vectors in
  `ker_Z(v)` with sinc weights from discontinuous interval indicators; a fixed
  GMC moment ranges over mass-constrained balanced channels with multinomial,
  factorial, radial, and coefficient weights. No map preserves either
  nonvanishing predicate. THM-878's `6/7` is the minimum of a pair-overlap
  functional averaged over primitive residue classes, not the coefficient
  `hat g(0)=6/7`; a repeated number does not identify the functionals. The
  script computes only the box `[-9,9]^n`, supplies no Fourier-tail bound, and
  its own output does **not** reproduce the direct integral (`0.0145` versus
  `0.0000` on `{1,2,3,4}`). Four low-dimensional examples cannot imply a
  uniform thirteen-speed theorem. More strongly, the known tight non-AP
  Goddyn--Wong row `{1,...,11,13,24}` has strict lonely measure zero, so its
  full nonzero-frequency sum cancels the main term. Thus the zero-measure locus
  is already non-AP; “AP-neighborhood” was never given a metric, radius, or
  theorem. Even a unique primitive relation does not save the transfer:
  `S={1,2}` at threshold `1/3` has kernel generated up to sign by `(-2,1)`,
  yet all harmonic multiples cancel the `1/9` main term and the strict lonely
  measure is zero. Sparse shortest relations do not control the signed sum of
  all higher relations. Moreover measure zero conflates a valid tight weak
  witness `M(S)=delta` with failure `M(S)<delta`; LRC needs the boundary data.
- **Correct framing:** the exact LRC relation-lattice expansion and zero mode
  are valid through THM-2047's Fejer-regularized limit (with earlier forms in
  THM-501/503/515/1061). They suggest diagnostics—relation height, support,
  Fourier decay, and initial channels—but transfer to GMC or tournament zeta
  requires an explicit weight-preserving map or a separately proved tail
  inequality; isolated weak witnesses additionally require phase-height/Euler
  boundary data. The S102 computation is a finite illustration only. It proves
  neither an AP-neighborhood reduction nor a new LRC theorem. THM-2051 is the
  proved sparse-relation localization; its residual contains every genuine
  support-`3..5` relation plane, not just AP cores.

## MISTAKE-234 (2026-07-21, THM-2070 audit) -- quotient data were promoted to weighted noncancellation

- **What was claimed:** S223 correctly encoded the return lengths
  `R={m:0 in mS}` of a two-sided Laurent support as an additive semigroup, then
  claimed that for mixed/complex coefficients only finitely many reachable
  lengths can cancel. It therefore called the Frobenius conductor an effective
  DvdK bound and announced a complete elementary replacement for the imported
  one-variable constant-term theorem.
- **Why it is wrong:** a support semigroup records which coefficient monomials
  exist, not whether their signed sum vanishes. The aperiodic support
  `S={-2,-1,1,2}` already has return lengths `2` and `3`, hence every
  `m>=2`. But for

  ```text
  f(z)=z-z^(-1)+z^2-z^(-2)
  ```

  one has `f(z^(-1))=-f(z)`. Constant-term invariance under inversion gives
  `CT(f^m)=(-1)^m CT(f^m)`, so every odd power vanishes. Thus cancellations
  persist infinitely far beyond the support conductor even though the support
  gaps have gcd one. S222/HYP-8890 explicitly leaves the general complex
  dominant-saddle/noncancellation step open, so citing it cannot close this
  gap.
- **The same lost-coordinate error in THM-1840:** the original functional
  corollary called `F` "diagonal by charge" (or said that it depended only on
  charge) and then used `F(h)=F(1)CT(h)`. Diagonality does not imply projection.
  For `Lambda=u+u^(-1)`, take `F(1)=1`,
  `F(u^2)=F(u^(-2))=-1`, and `F(u^k)=0` otherwise. This is a perfectly
  charge-labelled linear functional, but `F(Lambda^2)=0` while
  `F(1)CT(Lambda^2)=2`. The missing coordinate is the functional's weight on
  every nonzero charge, just as the semigroup quotient forgets every channel
  coefficient and phase.
- **Correct framing:** the return semigroup, pair-period law, and conductor are
  exact **reachability** data. For positive coefficients,
  `CT(f^m)!=0 iff m in R`; for signed or complex coefficients, retain channel
  phases and use a separate noncancellation theorem. HYP-8878's unique-minimum
  channel is one valid elementary sufficient case. HYP-8895 is therefore a
  useful support-sidecar construction, not a DvdK replacement or an effective
  coefficient-uniform bound. Likewise THM-1840's functional corollary requires
  the charge-**projecting** rule `F(u^k)=0` for every `k!=0`, with `F(1)!=0`;
  the two-charge constant-term formula itself remains proved. THM-2070 is the
  repaired theorem and supplies both exact counterexamples.

## MISTAKE-233 (2026-07-21, audit of S219/S220 / HYP-8880) -- classical theta and modular-curve facts were promoted to an LRC cusp-form obstruction without a map

- **What was claimed:** S219 called THM-515's sinc-weighted relation-lattice
  sum a modular theta series and transferred the binary-form decomposition
  `theta=Eisenstein+cusp` to LRC. S220 retracted that particular attachment but
  replaced it with the stronger claim that LRC scaling is `Gamma_0` level,
  clock divisors are modular cusps, a covering-min second moment is a form on
  `X_0(14)`, and the level-14 newform `f_14=14a` is the LRC(14) obstruction.
  It further identified genus with S218's hidden entropy and asserted that the
  period field of `f_14` is `Q(sqrt(-7))`.
- **Why it is wrong:** neither session constructs a map from an LRC row,
  good-set functional, or floor moment to a modular form, proves a modular
  transformation law, or shows that the LRC predicate is preserved. A cusp
  of `X_0(N)` is a `Gamma_0(N)`-orbit in `P^1(Q)`, not a subgroup of the
  additive clock `Z/NZ`; equality of two divisor counts is not an
  identification of objects. Dilation gives the exact LRC symmetry
  `M(cS)=M(S)`, but it cannot by itself be level structure: scaling the
  `12`-clock to `24` leaves `M` unchanged while `X_0(12)` has genus zero and
  `X_0(24)` has genus one. The proposed modular invariant therefore changes
  under an operation that preserves the alleged source object.

  The discriminant bridge is also false. The normalized level-14 elliptic
  newform has coefficient field `Q`; its elliptic-curve isogeny class `14a`
  is non-CM. In particular the standard curve `14a1` has
  `j=9938375/21952`, not an algebraic-integer CM `j`, so it has no CM period
  field `Q(sqrt(-7))`. Its coefficients `a_2=-1` and `a_7=1` are Fourier
  coefficients at bad primes, not labels for two cusps. `X_0(14)` is only the
  first positive-genus curve in the selected subfamily `X_0(2p)` with
  `p=3,5,7,...`; `X_0(11)` already has genus one. Finally, genus is a vector-
  space dimension, not the observable-fiber quantity refuted in MISTAKE-231,
  and class number greater than one does not alone force a cuspidal theta
  component: genus characters can account for the class variation.
- **Correct framing:** the finite identities in S219 survive as independent
  classical examples: for `x^2+xy+2y^2`,
  `r(n)=2 sum_(d|n)(d/7)`; for discriminant `-23`, the stated class average is
  Eisenstein and the principal-minus-nonprincipal theta series is a nonzero
  weight-one CM cusp form. The cusp-count formula, genera of `X_0(12)` and
  `X_0(14)`, the displayed coefficients and Atkin--Lehner signs of `f_14`,
  and the standard Rankin--Selberg factorization are likewise classical facts.
  They have no demonstrated LRC consequence. HYP-8880 and the LRC portions of
  S219/S220 are **REFUTED**; retain the artifacts only as historical records
  and do not route agents to them as an LRC strategy without a typed,
  predicate-preserving construction.

## MISTAKE-232 (2026-07-21, audit of S99 / HYP-8876) -- a scale-and-clock word analogy was called an exact proof-shape bridge, and a composite modulus was assigned a Paley spectrum

- **What was claimed:** S99 identified GMC(2) dilation/Frobenius with
  THM-2057's scaled clock witnesses, wrote
  `Z/pZ <-> Z/12aZ,Z/14aZ` and `Frobenius <-> modular-orbit periodicity`, and
  said that each of the clock moduli `{7,13,14}` carries a Paley
  `sqrt(p)` spectrum. Its LRC row also used `84a=12a*7a` and invoked the
  transitive-core/Paley/class-number entropy syntheses as explanation.
- **Why it is wrong:** the script constructs Paley objects only at the prime
  moduli `7` and `13`; `14` is not a prime power, has no finite field
  `F_14`, and no Paley graph or tournament was constructed. Printing
  `14=2*7` is not a spectral calculation. The displayed divisibility identity
  already fails at `a=2`: `84a=168`, while `(12a)(7a)=336`. The exact
  statement is `lcm(12a,14a)=84a`. More fundamentally, Frobenius is a ring
  endomorphism in characteristic `p`, whereas the composite LRC clocks supply
  only additive residue orbits. For example on `Z/12Z`, the proposed
  seventh-power analogue is not additive:
  `(1+1)^7=8 (mod 12)` but `1^7+1^7=2 (mod 12)`. No source-to-target map,
  preserved predicate, or sidecar turns the shared words “scale” and “clock”
  into a proof transfer. MISTAKES-227, -228, -230, and -231 independently
  block the transitive-core, Paley-LRC, and zero-entropy explanations.
- **Correct framing:** two exact, independent spectral facts survive. For the
  Paley tournament at `7`, its skew matrix satisfies
  `S^2=J-7I`, hence has spectrum `0,(+i sqrt(7))^3,(-i sqrt(7))^3`. For the
  Paley graph at `13`, `A^2+A=3I+3J`, hence the spectrum is
  `6,((-1+sqrt(13))/2)^6,((-1-sqrt(13))/2)^6`. THM-2057's scaled missing-
  clock proof also survives independently. “Scale, then close on a clock” may
  be kept only as a heuristic phrase; HYP-8876 is not verified and yields no
  theorem or transfer.

## MISTAKE-231 (2026-07-21, audit of HYP-8875 / S218) -- four unrelated fiber sizes were called one repo-wide entropy invariant

- **What was claimed:** S218 defined
  `H_arith(X|L)=log2|{X':L(X')=L(X)}|` and identified it with class-group
  depth, tournament reconstruction, continued-fraction partial-quotient size,
  and LRC/GMC moment depth. It then asserted that every rigid extremum has zero
  arithmetic entropy and every hard object has positive entropy.
- **Why it is wrong:** the displayed quantity is at most a finite-fiber
  Hartley information after fixing an ambient universe, equivalence relation,
  and observable. It is not Shannon entropy without a probability law, and
  changing `L` can force any chosen object to have fiber size one or the whole
  universe. For the two isomorphism classes of three-vertex tournaments, a
  constant `L` gives the transitive class a fiber of size two, while the
  identity observable gives every class a singleton fiber. Thus the proposed
  rigid-versus-hard law can be reversed by changing the observable. The four
  examples do not use one map or even one type of fiber.
  In particular:
  - the tournament claim is already false at `n=5`: the regular score
    `(2,2,2,2,2)` has one isomorphism class, whereas the score
    `(1,2,2,2,3)` has three. Those three classes are not cospectral; their
    adjacency characteristic polynomials are
    `x^5-4x^2-4x-1`,
    `(x^2+x+1)(x^3-x^2-3)`, and
    `(x^2+x+2)(x^3-x^2-x-1)`. MISTAKE-227 also blocks identifying a
    translated transitive score vector with the LRC AP relation carrier;
  - a finite continued-fraction partial-quotient geometric mean is not an
    entropy and has no maximum (`[0;N]` makes it arbitrarily large). Every
    nonterminal finite CF prefix has continuum many real extensions, while a
    complete rational expansion has a singleton fiber, so the displayed
    Hartley formula cannot distinguish golden prefixes from `[0;13,14]`;
  - bounded LRC danger count gives exact inclusion--exclusion by depth `13`,
    while `B5>0` is only a sufficient certificate, not a uniform depth-five
    classifier. Explicitly, the distributions with weights
    `(1,15,15,1)/32` on `{0,2,4,6}` and `(6,20,6)/32` on `{1,3,5}` have the
    same moments through degree five, by the sixth finite-difference identity,
    but masses `1/32` and `0` at zero. THM-1790 gives an unbounded-in-degree
    GMC moment-depth floor, not an infinite entropy; and
  - the genus routine labels its formula a heuristic and silently falls back
    when divisibility fails. It is a finite table, not a proof of genus theory
    or of invisibility to every congruence test. Its prose also calls the
    within-genus part “odd,” but the fundamental discriminant `D=-39` has
    `h=4`, two genera, and even within-genus fiber size `2`.
- **Correct framing:** for a fixed finite universe and quotient, define the
  explicit Hartley fiber size `H_0=log2|L^{-1}(L(X))|`; add a probability law
  before using Shannon language. The exact survivors are the `n<=5`
  tournament fiber census, transitive-score uniqueness, elementary
  score-distribution entropy, the displayed reduced-form tables, and the
  qualitative contrast between finite LRC inclusion--exclusion and no
  degree-uniform GMC cutoff. These are separate examples of information loss,
  not one invariant and not an LRC theorem. HYP-8875 is retracted as a unified
  mathematical claim.

## MISTAKE-230 (2026-07-21, audit of HYP-8870 / S217) -- inverse form classes were counted as distinguishable rational-prime outcomes

- **What was claimed:** S217 treated the `h(D)` proper classes of positive
  binary quadratic forms as `log2 h(D)` hidden bits beyond the Legendre symbol.
  For `D=-23` its script assigned split rational primes among the three forms
  with counts `[2,9,0]`, then transferred class-number-one rigidity to the
  THM-2053 LRC gate.
- **Why it is wrong:** proper inverse classes `(a,b,c)` and `(a,-b,c)` represent
  exactly the same rational integers by `y -> -y`. The two nonprincipal
  `D=-23` forms therefore represent the same displayed primes; choosing
  `reps[0]` manufactures the asymmetric `[2,9,0]` count. A rational split
  prime determines an inversion orbit of ideal classes, not an oriented class.
  Moreover `log2 h` is not an entropy until a random variable and distribution
  are specified. For oriented prime ideals one can formulate an Artin-class
  distribution and invoke Chebotarev; for rational primes the inversion
  quotient is different. “Invisible to any local test” also overstates
  “not determined by the quadratic Legendre symbol”: Frobenius in the Hilbert
  class field is itself local prime data. The finite search bound is evidence,
  not a proof of the global representation criterion. Finally, MISTAKE-229
  already shows that THM-2053 has no discriminant `-7` form to which any of
  this could transfer.
- **Correct framing:** the reduced-form and class-number tables are useful
  classical checks. With a precisely oriented prime-ideal sample, class-group
  size can support a Hartley/Shannon analogy. It has no demonstrated LRC
  consequence. The S217 use of `HYP-8870` also collided with the normal-fan
  program now correctly numbered HYP-8871; the S217 claim has no canonical
  hypothesis ID and is retracted.

## MISTAKE-229 (2026-07-21, audit of HYP-8865 / S216) -- the tangent-disk gate was assigned a nonexistent Heegner discriminant

- **What was claimed:** S216 called THM-2053's residual the short vectors of
  one binary quadratic form, inferred discriminant `-7` from `LRC(14)=2*7`,
  used `h(-7)=1` to declare every residual plane a single rigid class, and
  identified rank-versus-Euler with isotropic-versus-anisotropic directions.
- **Why it is wrong:** THM-2053 compares thirteen linear determinant forms
  `a z_i-b u_i` with the chosen-coordinate Euclidean norm `a^2+b^2`. The
  failure locus is a **union of 26 tangent disks**, not a representation set
  for one quadratic form. The Euclidean form has discriminant `-4` and changes
  covariantly with the chosen saturated basis; neither the number `14` nor a
  Paley spectral factor supplies a plane discriminant `-7`. At the prime `p`,
  the discriminant `-p` has Legendre symbol zero, so the script's computation
  of `(-1/p)` did not test its stated local criterion. An extra bounded
  relation is a linear annihilator condition on all columns, not an isotropic
  parameter direction, and failure of the sufficient gate does not imply
  either unsafety or Euler survival.
- **Correct framing:** retain the exact wedge identity and the classical
  class-number table as separate facts. The useful LRC geometry is the
  intersection of a primitive lattice and positivity cone with 26 explicit
  tangent disks, considered up to saturated `GL_2(Z)` basis changes. THM-2055
  compresses the determinant sidecar exactly to the signed column polygon's
  normal fan: hull vertices suffice for the maximum, and each owner cone meets
  one tangent disk. Non-hull runner data, pair sums, and endpoint ownership
  remain necessary sidecars. Any
  arithmetic-form compression must first construct an actual invariant form
  from the full column configuration and prove that it preserves the gate;
  Heegner or Paley language alone supplies no such map. See MISTAKE-230 for the
  independent inverse-class error in S217's entropy continuation.

---
## MISTAKE-228 (2026-07-21, audit of HYP-8860 / S215) -- the odd-prime Paley dichotomy was promoted to an LRC periodic table

- **What was claimed:** S215 assigned primes `2,3,5,7,11` intrinsic LRC roles
  from Paley spectra, identified adjacency eigenvalues with Gauss sums, and
  treated `p mod 4` as a tight/slack law. It called `2` part of the same Paley
  theorem, merged the `Phi_3` and `Phi_6` factors, and linked the Paley-7 field
  directly to the real cubic cap field.
- **Why it is wrong:** the verified classical statement is for odd primes.
  For `p=3 mod 4` the nonprincipal Paley-tournament eigenvalues are
  `(-1+-i sqrt(p))/2`, while the Gauss sum is `i sqrt(p)`; for `p=1 mod 4`
  the Paley graph also has principal eigenvalue `(p-1)/2`. The script does not
  test `p=2`. `Phi_3=x^2+x+1` and `Phi_6=x^2-x+1` need a sign twist, and
  `Q(sqrt(-7))` is not the real cubic field of discriminant `49`. The roles
  “5 = Fibonacci foil” and “11 = rank” supply no map, and multiples of `11`
  are not a forced speed coordinate.
- **Correct framing:** retain the quadratic-residue Cayley relation and the
  odd-prime Paley graph/tournament spectral formulas as tournament background.
  Treat the small-prime table as a modulus-selection mnemonic only. It proves
  no LRC tightness, slackness, rank, or transport statement.

---
## MISTAKE-227 (2026-07-21, audit of HYP-8855 / S214 and HYP-4552 / S25) -- a finite-index AP chain frame was called a saturated relation-lattice basis

- **What was claimed:** the rows
  `d_k=(k+1)e_k-k e_(k+1)` were called a `Z`-basis of
  `L(AP)=ker_Z(1,2,...,12)`, their tridiagonal Gram matrix was assigned to the
  whole lattice, and they were identified with THM-2052's eleven private
  anchor-star relations. The resulting rank coincidence was promoted to a
  common AP/tournament nullcone and to the claim that the rank and Euler
  branches meet at the AP. The 60 minimal vectors were also called additive
  energy.
- **Why it is wrong:** the `d_k` span only a sublattice `D` of index
  `11!=39,916,800`. Exactly,
  `det Gram(D)=1,035,678,099,456,000,000`, while
  `disc L(AP)=650`; their square-root ratio is `11!`. The saturation contains
  the norm-three Schur rows. THM-2052 instead concerns a bounded signed code
  on thirteen coordinates and its private rows depend on two chosen anchors;
  no map to this twelve-coordinate lattice was given. The kernel rank is
  `|S|-1=11`, not `|S|-2`. The 60 vectors are the signed off-diagonal Schur
  triples, not four-term additive energy. Tournament scores `0,...,11` are a
  translation of speeds, and LRC is not translation invariant. Reversal of
  phase is universal, not an AP-fixed configuration theorem. THM-2053 also
  shows rank twelve is not the only finite terminal.
- **Correct framing:** `D` is a useful finite-index path/Jacobi frame and
  `L(AP)=sat(D)` has discriminant `650`; saturation is exactly what restores
  the short Schur relations that the path frame misses. The AP, transitive
  order, palindromic pair sums, and Schur census remain diagnostic views, not
  one object or an LRC descent. Before promoting a rank match, compute the
  ambient type, explicit map, saturation index, and discriminant.

---
## MISTAKE-226 (2026-07-21, codex audit of HYP-8840) -- unanchored additive energy, LRC Fourier volume, and GMC radial moments were promoted to one constant-term bridge

- **What was claimed:** S211 promoted nine finite AP energy maxima to an
  all-order Wall-A theorem, called THM-730 its `m=2` case, identified the
  unanchored additive-energy lattice with the full LRC relation lattice, used
  separated-digit tensorization as a rank reduction for arbitrary cores, and
  treated a geometric model weight as the GMC radial kernel. It then called
  the one-dimensional LRC maximum a minimax saddle linking GMC volume to Euler
  positivity.
- **Why it is wrong:** `CT[P^m Pbar^m]` counts only affine relations
  `k dot v=0` with augmentation `sum k_i=0`; the LRC Fourier formula uses the
  full linear lattice. THM-730's Schur relation `a+b=c` has augmentation one,
  so it is not an unanchored two-energy relation. The script checked only nine
  `(N,k,m)` cases. Translation preserves every unanchored energy but changes
  Schur data: `{1,...,5}`, `{2,...,6}`, and `{1,3,5,7,9}` have the same
  unanchored energy at **every** order, yet their LRC maxima are `1/6,1/4,1/2`.
  Adjoining observer zero gives anchored `E_2=146,130,106` and Schur counts
  `10,6,0`. Separated-digit factorization
  requires a no-carry scale and therefore removes the carry relations that
  make an AP extremal; it does not factor an arbitrary 13-speed core. The
  geometric weight was invented, and a raw box-truncated sinc sum is not the
  Fejer-regularized measure in THM-2047.
- **Correct framing:** retain the universal identity
  `E_m(S)=CT[P^m Pbar^m]`, the single-character CT seed, conditional no-carry
  tensorization, the AP cyclotomic factorization, and the exact volume ceiling:
  positive Haar measure certifies a strict exit, while Euler characteristic
  also sees a tight isolated witness. A faithful energy comparison must adjoin
  the stationary observer `0`; then a linear relation lifts by
  `k_0=-sum k_i`, and THM-730 is the `a+b=c+0` sector of anchored `E_2`.
  Equivalently, retain the mixed table
  `M_(r,s)=CT[P(z)^r P(z^(-1))^s]`: it selects augmentation `r-s`, THM-730 is
  exactly `M_(2,1)`, and every full LRC relation appears in some bidegree.
  Whether this augmentation sidecar helps Wall A is OPEN. Do not call the
  maximum of a lower envelope a game/Morse saddle without an actual
  two-player or variational correspondence.

---
## MISTAKE-225 (2026-07-21, codex audit of HYP-8835) -- skew parity, tournament games, torus topology, and three sign actions were promoted to one antisymmetry theorem

- **What was claimed:** S210 asserted that antisymmetry alone forces odd
  tournament-game support; that a pure optimal strategy exists iff the
  tournament is transitive; that intransitivity is precisely toroidal
  recurrence; that the solid-bagel cutting deficit is a reduced Euler term of
  `T^2`; and that matrix complement, torus inversion, and Vandermonde parity
  are one involution.
- **Why it is wrong:** odd optimal support also uses the tournament block's
  mod-2 form `I+J`, which is nonsingular in even order. A four-vertex
  tournament with vertex `0` beating all and `1->2->3->1` is intransitive but
  has pure optimum `e_0`; replicator dynamics even satisfies
  `x_0'=x_0(1-x_0)`. Antisymmetry does not force zero column sums. The regular
  three-cycle's conserved product comes from a positive kernel vector, and its
  positive levels in `Delta^2` are circles `S^1`, not a surface `T^2`. The
  cutting bagel is `D^2 x S^1`, while the `1,2,1` Morse count belongs to its
  boundary. MISTAKE-222 already blocks the sequence/Euler jump. Finally the
  three `Z/2` actions have no supplied equivariant map. The RK4 check was
  nonsymplectic and rounded a nonzero return distance to `0.0000`.
- **Correct framing:** Fisher--Ryan's tournament-game support is unique and
  odd; the parity proof needs both skew singularity and the tournament mod-2
  block. Pure optimum is equivalent to a Condorcet winner. If `Mp=0`, then
  `d/dt sum p_i log x_i=0`, so a positive equilibrium vector yields the
  monomial invariant `product x_i^(p_i)`. Keep game, Morse, and dynamical
  saddles distinct; keep solid and boundary tori distinct; and require an
  equivariant map before transferring any sign representation.

---

## MISTAKE-224 (2026-07-21, codex audit of HYP-8830) -- a Fourier relation lattice was identified with a toric layer poset, and a thickened safe set with an arrangement complement

- **What was claimed:** S209 called `{k in Z^n:k.v=0}` the layers of a finite
  toric arrangement, identified `|G_delta|` with ordinary toric-complement
  volume and its Fourier series with an arithmetic-Mobius sum, and promoted a
  `B=2` relation count to Betti/Mobius mass. It also called LRCMod a finite-
  field characteristic-polynomial engine and braid cohomology a per-tournament
  invariant.
- **Why it is wrong:** the lattice indexes Fourier characters surviving
  integration over the orbit subtorus; toric layers are connected components
  of intersections of a **finite** character list. THM-1820's series is
  infinite, sinc-weighted, and `delta`-dependent. For `S={1}`,
  `|G_delta|=1-2delta`, whereas the ordinary circle complement of the
  character kernel has Haar measure `1`. The script computed no layer poset,
  Mobius function, arithmetic Tutte polynomial, or LRCMod equivalence. Its AP
  count is cutoff-dependent over exactly 72 triples; the natural circle-
  hypertorus union has four points for `(1,2,3)` but six for `(2,3,4)`.
- **Correct framing:** `L_v={k:k.v=0}` is the annihilator of
  `H_v=image(t |-> (v_jt))`, and `G_delta=phi_v^(-1)(F_delta^n)`. The full
  embedded `L_v` determines primitive `v` up to scale/relabeling; only an
  abstract, truncated, or summarized lattice is automatically lossy. The
  useful finite toric object has characters
  `X_S={(v,+1),(v,-1)}` on `(t,delta)`, with walls
  `vt +/- delta in Z`. Its selected inequality side, height, owner, sign, and
  paired deletion are essential. Opposite-sign wall determinants are `v+w`,
  recovering the pair-sum ruler. THM-2047 now proves this exact carrier and
  paired-deletion operation; localization forcing an AP core remains open and
  is not a proved Wall-A consequence.

---

## MISTAKE-223 (2026-07-21, codex audit of HYP-8825) -- a braid-flat factorization and common determinant syntax were promoted to hyper-Bessel and Euler-characteristic factorizations

- **What was claimed:** S208 called cake, bagel, full g-bonacci, and transitive-
  tournament counts one Zaslavsky/Mobius inversion; identified
  `bagel(n)-cake(n)=T_n-1` and a finite-shadow first deficit as the same reduced
  Euler characteristic; called the Vandermonde's within-block vanishing order
  a flat codimension; and asserted that braid-arrangement localization factors
  the NC2 hyper-Bessel boundary into block factors. It also called general
  `Phi_(p,q)` Laguerre--Polya open and described a duplicate recurrence
  comparison as a companion-determinant/topology verification.
- **Why it is wrong:** the `n!` chambers are components of the **real** braid
  complement; the complex complement used in the note is connected. For a
  `3+2` coincidence partition, the ordinary Vandermonde vanishes to order
  `C(3,2)+C(2,2)=4`, while the flat codimension is `(3-1)+(2-1)=3`.
  Localization gives within-block braid factors and a nonzero cross-block unit;
  no map from that polynomial expansion to THM-2017's univariate channel limit
  `Phi_(p,q)` was stated. THM-2023 had already proved every positive-integer
  `Phi_(p,q)` is type-I Laguerre--Polya by Gauss multiplication and the
  Baricz--Singh theorem. A bagel is a solid torus, not `T^3`; no common complex,
  Mobius function, or valuation was supplied for its regions and the shadow
  lattice. The script compared two implementations of the same recurrence and
  computed neither a companion determinant nor an Euler characteristic.
- **Correct framing:** retain four separate exact shadows: the Vandermonde is
  the braid defining polynomial; its **real** chambers are labelled total
  orders/transitive tournaments; near a coincidence flat it has a leading
  within-block/transverse product with exponent equal to vanishing order; and
  an explicit companion matrix satisfies
  `det(I-x*M_g)=1-x-x^(g+1)`. THM-2033 applies to its special moment matrix, not
  every NC2 moment. The cake--bagel identity has the more direct decomposition
  `bagel(n)=cake(n+1)-2`, so
  `bagel(n)-cake(n)=(T_n+1)-2`. The possible topology bridge remains an open
  source-to-target valuation test under MISTAKE-222, not a theorem.

---

## MISTAKE-222 (2026-07-21, codex audit of HYP-8820) -- a shared binomial ambient array and a matching minus-one offset were promoted to a geometric and cross-problem mechanism

- **What was claimed:** cake, bagel, Moser, and Fibonacci were described as
  ordinary truncated rows or diagonals of “one Pascal triangle”; the identity
  `bagel-cake=T_n-1` was said to prove that the torus hole **is** the deficit-one
  boundary in a g-bonacci shadow lattice and to place JC(2) and LRC(14) on the
  same scaffold. A broken polygonal-skip calculation was also printed as if it
  agreed with the cited indexing.
- **Why it is wrong:** caterer and cake are prefixes of the `n`th Pascal row.
  Moser is the even-column functional `C(n,0)+C(n,2)+C(n,4)`, equivalently the
  first five terms of row `n-1`; bagel is the weighted functional
  `C(n,3)+2C(n,2)+2C(n,1)`, equivalently `cake(n+1)-2`. Thus they are not all
  the same unshifted prefix operation. A common ambient array is a useful
  coordinate system, not a map between the counted objects. The same scalar
  offset `-1` supplies no bijection, operation, or preserved predicate, and the
  JC/LRC analogy supplies none either. The removed polygonal-skip code produced
  negative values and contradicted its own “first deviation” label.
- **Correct framing:** the four binomial formulas, Fibonacci diagonal identity,
  full-rank gap-diagonal g-bonacci identity in klein-S313's explicit indexing,
  and
  for `n>=1`, `bagel(n)=cake(n+1)-2` and
  `bagel(n)-cake(n)=C(n,2)+C(n,1)-1=T_n-1` are exact. Treat them as a binomial-
  reading atlas. Finite-rank shadows depart at their first deficit. The
  shadow-lattice/torus link remains a useful test only if both defects are
  derived from one explicit boundary-cell or valuation operation; JC and LRC
  links additionally require a preserved predicate and loss ledger.

---

## MISTAKE-221 (2026-07-21, codex audit of HYP-8815) -- a finite LRC witness scan and an AP/Fibonacci contrast were promoted to a necessary autocorrelation characterization

- **What was claimed:** HYP-8815 called a scan over rationals `a/q` with
  `q<=Qmax` an exact computation of `M`; inferred that every hypothetical
  LRC(14) counterexample must be near-AP, anti-golden, organized around one
  CF-blocking far element, and beat the AP by a higher-order autocorrelation
  mechanism; and used smaller THM-731 discrepancy in the wrong direction.
- **Why it is wrong:** a denominator-truncated scan computes only
  `L_Q(S)<=M(S)` unless a completeness theorem covers every maximizer.
  The exact strict witness is
  `L_1200({1,...,12,5460})=92/1197 < M=420/5461`.
  THM-730 proves only the AP's unique Schur-triple maximum. THM-731 is
  peel-dependent and gives
  `L_cert=(6/7)|G'|-sqrt((6/49)disc_v)`, so **smaller** discrepancy
  strengthens safety. A few loose Fibonacci packets prove no uniform
  exclusion, near-AP necessity, continued-fraction condition, or ownership of
  every covering modulus by one far speed.
- **Exact computational repair:** the pair-sum theorem
  [THM-1002](theorems/THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case.md)
  puts every maximizer at `t=p/(v_i+v_j)`, including the self-pair case.
  Enumerating every numerator on every pair-sum ruler computes `M` exactly;
  alternatively `Q>=2 max(S)` is a sufficient complete denominator cutoff.
  Do not filter to coprime numerators unless all divisor rulers are also
  present: `S={1,4,5}` attains its maximum at `2/6=1/3`, while `3` is not
  a pair sum. The repaired script validates inputs and computes all fifteen
  displayed rows exactly; all are safe. This finite bank is not a global
  minimizer theorem.
- **Correct structural framing:** after gcd normalization, a counterexample is
  primitive Cover14 with `M<1/14`. THM-1017 excludes a maximum-deletion core
  `d{1,...,12}`, so THM-730 gives that core Schur count `T<=65`. Since a
  counterexample has lonely-time measure zero, THM-731 forces
  `disc_v>=6|G'_{~v}|^2` for every peel. Anti-golden, near-AP,
  Fibonacci-foil, one-far-blocker, and joint-order autocorrelation language is
  a hostile-control program, not a characterization.

  `L_1200=92/1197 < M=420/5461` at `t=420/5461`.

  More starkly, for `Q>=14`, dilate the deep well by `lcm(2,...,Q)`. Every old
  sampled phase then has value zero while dilation invariance leaves
  `M=14/183`. A bounded lower bound above `1/14` does safely exclude one row,
  but a finite candidate list cannot prove a global minimizer or ranking.

  Primitivity is a WLOG normalization, not a literal condition on every
  counterexample. The old random generator did not validate cardinality and
  distinctness, although its six displayed seed-2026 rows happen to be valid.
  THM-730 counts Schur triples, not four-term additive energy; Fibonacci is not
  minimal for either statistic. THM-731 is peel-dependent and gives
  `L_cert=(6/7)|G'|-sqrt((6/49)disc_v)`, so **smaller** discrepancy strengthens
  safety—the claimed direction was backwards. Cover14 distributes divisibility
  among all runners, not necessarily a twelve-runner core plus one far blocker;
  THM-1017 excludes only a maximum-deletion core `d{1,...,12}` in the relevant
  normalized regime.
- **Correct framing:** exact `M` computation uses every numerator, including
  non-coprime ones, on the THM-1002 pair-sum rulers. The corrected script
  validates its inputs and computes all fifteen rows exactly; every one is safe,
  and none is a disproof. The rigorous necessary kernel is: after gcd
  normalization a counterexample is primitive Cover14 with `M<1/14`; its
  maximum-deletion core is not `d{1,...,12}`, so THM-730 gives Schur count
  `T<=65`; and `L=0` with THM-731 forces
  `disc_v>=6|G'_{~v}|^2` for every peel. THM-2048 now refines the latter by
  the exact fiber-occupancy tax
  `{7vmu}(1-{7vmu})/(7v^2)` and its THM-732 integer obstruction. Anti-golden,
  near-AP, Fibonacci-foil,
  one-far-blocker, and joint-order-autocorrelation language remains a useful
  hostile-control program, not a characterization.
---
## MISTAKE-215 (2026-07-21, codex audit of HYP-8795 / THM-2040) -- a prime-local minimum-face normalization was promoted to a global common-factorial/Vandermonde factorization

- **What was claimed:** dividing every NC2 moment by a purported common
  factorial `(p*A0)!` was said to leave a bounded residue and to make general
  NC2 noncancellation equivalent to one (possibly confluent) Vandermonde.
- **Why it is wrong:** `(p*A0)!` in THM-2022 is attached to one deliberately
  chosen amplified moment and its lowest balanced face. It divides all channel
  factorials there, but off-face quotients remain nontrivial; it does not peel
  every general moment to a universal determinant. The accompanying S91
  computation itself gives `E[P^m]/m! = 2, 7/2, 17/3, ..., 1441729/40320`
  for `P=Z+(1+s)+Zbar`, an increasing sequence through the displayed range,
  contradicting its printed "bounded" conclusion. A moment is a scalar sum of
  channel monomials, whereas the Vandermonde in THM-2033 is a determinant of a
  special moment matrix; division by one scalar factorial cannot identify the
  two in general. The claimed positive-zero / Laguerre--Polya / Paley iff chain
  also has no logical basis and repeats the type error in MISTAKE-214.
- **Correct framing:** retain two valid statements. On the exact symmetric
  monomial wall `P=a Z^p + beta*s^(p/2) + c Zbar^p` (necessarily even `p`),
  every balanced channel really has radial degree `mp/2`, so the special
  central-trinomial factorization is exact. In full NC2, THM-2022 provides a
  **prime-local initial-form normalization** only: at order `p*m0`, division by
  `(p*A0)!` kills higher/carry layers modulo a good prime and leaves the whole
  lowest-face residue `Q^p`. That residue is generally a face constant term,
  not a Vandermonde and not a Paley spectrum.

---

## MISTAKE-214 (2026-07-21, codex audit of HYP-8785/8790) -- repeated Vandermonde nodes were identified with repeated tournament scores and hence with a regular/Paley tournament

- **What was claimed:** HYP-8785 asserted that equal radial channel degrees on
  the NC2 central wall are equal scores in the tournament expansion of the
  Vandermonde, hence a regular, doubly regular, or Paley tournament. HYP-8790
  then promoted the true central-trinomial channel-weight identity to claimed
  equivalences among NC2 noncancellation, Laguerre--Polya zero geometry, and the
  Paley critical-line spectrum.
- **Why it is wrong:** in
  `V(x)=sum_T sgn(T)*x^score(T)`, the `x_i` are Vandermonde **node variables**
  and `score_i(T)` are their **exponents**. Setting two node values equal does
  not set two score exponents equal; it causes cancellation among the evaluated
  monomials. The type mismatch already fails at the smallest tied core: for
  symmetric `m=2` there are two equal-degree channels, while no regular
  tournament exists on two vertices (the scores are `0,1`). More generally
  tied channel counts have arbitrary parity, whereas regular tournaments need
  odd order, doubly regular tournaments require order `3 mod 4`, and Paley
  tournaments exist only at special prime-power orders. The identity
  `sum_i m!/(i!^2(m-2i)!)=[x^0](1+x+x^-1)^m` is valid, but it does not imply the
  stated Paley/free-probability/zero-set equivalences.
- **Correct framing:** retain THM-2033's determinant/Vandermonde identity and
  HYP-8790's central-trinomial coefficient identity. Confluence means repeated
  **nodes** and calls for divided differences or derivative rows; regularity is
  a property of a tournament's **score vector** and is not forced. Paley and
  free-probability parallels may be used as analogies, not iff reductions or an
  identification of the NC2 wall. THM-2022 supplies the actual tied-face
  certificate: after common-factorial normalization the whole face survives as
  `Q^p`, with no tournament-score interpretation required.

---

## MISTAKE-213 (2026-07-21, codex proportional-NC2 audit) -- finite toral-zero recurrence was treated as necessary to combine an exact factorization with EMP

- **What was claimed:** the first THM-2021/HYP-8771 write-up said that
  no-consecutive-zero for the toral factors `A_m(kappa)` was too weak for NC2,
  because it did not make `A_m(kappa)` eventually nonzero and therefore could
  not force `L(b^m)=0` eventually from
  `E[P^m]=A_m(kappa)L(b^m)=0`. It promoted finite recurrence of a fixed toral
  zero to the load-bearing proportional-slice target.
- **Why it is wrong:** the supports were combined in the wrong direction. If
  `b!=0`, EMP says `L(b^m)!=0` for **every sufficiently large** `m`. Moment
  nullity would therefore force `A_m(kappa)=0` for every sufficiently large
  `m`. Any theorem showing merely that `A_m(kappa)` is nonzero arbitrarily far
  out contradicts this. In particular, no-consecutive-zero is already more
  than enough. For the symmetric Legendre factors the three-term recurrence
  forbids consecutive zeros (with `kappa=1/4` handled directly), so it gives a
  self-contained full symmetric proportional closure. THM-2018's root-of-unity
  EGF proves the needed non-eventual-zero statement for every charge pair.
- **Correct framing:** for a factorization `M_m=A_m L_m` with `L_m` eventually
  nonzero, NC2 needs only that `A_m` is **not eventually zero**. Finite zero
  recurrence is a strictly stronger sequence-profile theorem, interesting in
  HYP-8771 and announced for Legendre polynomials, but not an NC2 dependency.
  Keep three distinct predicates separate: no consecutive zeros; unbounded
  nonzero support; finite zero support. The first implies the second, and the
  second is sufficient here. See THM-2018, corrected THM-2021, and HYP-8771.

---

## MISTAKE-212 (2026-07-21, codex audit of HYP-8772) -- channel-tournament transitivity was promoted from a sufficient dominance certificate to an iff characterization of NC2 noncancellation

- **What was claimed:** the reservation version of HYP-8772 said that NC2
  noncancellation is equivalent to transitivity/source-dominance of the
  channel-degree tournament, and that the moment nullcone *is* a tournament
  nullcone on its channels. Resonance was identified with a tied/regular
  tournament and therefore with the only possible cancellation wall.
- **Why it is wrong:** THM-2021 gives exact noncancellation inside a fully tied
  central channel family. Take `b=s(s-2)`, `a=s-2`, `c=s(s-2)`, and
  `P=Z*a(s)+b(s)+Zbar*c(s)`. Then `h=s*a*c=b^2`, so every return channel in
  `E[P^m]` has the identical radial word `b^m` and identical factorial degree:
  there is no asymptotic source. Nevertheless `E[P]=L(b)=0` and
  `E[P^2]=S_2(1)L(b^2)=3*8=24`, so this two-sided polynomial is immediately
  non-null. More generally THM-2021 proves NC2 on this tied proportional slice
  outside a countable exceptional set. Also, a dominance relation with exact
  ties is not a tournament until an extrinsic tie-breaker is imposed.
- **Correct framing:** a strict channel source/transitive dominance order is a
  **sufficient certificate** for noncancellation and explains why THM-2017's
  degree-gap argument works. Its failure identifies where that certificate
  stops, not where NC2 fails or can fail. In a tied core, phase-sensitive
  algebraic information can still prove noncancellation (Legendre recurrence
  in THM-2021, the all-charge EGF in THM-2018, and Sheffer/resultants in
  HYP-8769). The channel
  tournament remains a useful regime classifier, but the claimed iff and
  "moment-nullcone = tournament-nullcone" statements are withdrawn. See
  THM-2021, HYP-8771, and the correction banner in the HYP-8772 reflection.

---

## MISTAKE-218 (2026-07-21, codex audit of THM-1979/2013) -- score variance was promoted from a cyclicity coordinate to the whole structural axis

- **What was claimed:** THM-1979 put `σ²=Var(scores)` in
  `[0,(n²−1)/12]` for every `n`, called zero the regular center, and said
  fiber size, strong fraction, and modular-prime fraction vary monotonically
  opposite to `σ²`. THM-2013 consequently wrote
  `c₃/c₃,max=1−σ²/σ²_max` uniformly in `n`. The limit discussion
  identified tournamenton space with its degree functions.
- **Why it is wrong:** for even `n`, integral scores with half-integral mean
  force minimum variance `1/4`, not zero, and the maximum cyclic count is
  `n(n²−4)/24`. More fundamentally, the frozen S198 census refutes the
  monotonic structural claim: at `n=6`, variance `5/4` supports both size-3
  all-strong score fibers and size-1 reducible fibers; at `n=7`, maximum fiber
  size 47 occurs at both variances `4/7` and `10/7`, while the zero-variance
  fiber has size 3. A tournamenton is not determined by its degree function;
  many structured regular tournamentons share `d_W=1/2` with the quasirandom
  tournamenton.
- **Correct framing:** set `ε_n=0` for odd `n` and `1/4` for even `n`, and
  `σ²_tr=(n²−1)/12`. Then
  `c₃=n(n²−1)/24−(n/2)σ²` and
  `τ=c₃/c₃,max=(σ²_tr−σ²)/(σ²_tr−ε_n)`.
  The variance determines the cyclic shell exactly but forgets score shape and
  within-fiber structure. The finite score map and limiting degree map are
  projections with rich fibers, not complete parameterizations. THM-1979,
  THM-2013, their scripts/results/reflections, the hypothesis index, and the
  session log are corrected.

---

## MISTAKE-217 (2026-07-21, codex audit of THM-2013/2016) -- signed Redei data was labeled as the invariant magnitude `|R|`

- **What was done:** the continuum-coordinate scripts grouped tournament
  classes by `signed_redei(A,n)` while their output and THM-2013/2016 called
  the coordinate `|R|`.  THM-2016 then attributed signed-key resolution counts
  `36→44` and `41→50` to the absolute invariant.
- **Why it is wrong:** the sign changes under the parity gauge of a relabeling;
  the isomorphism invariant in THM-1966 is its magnitude.  Signs are mixed in
  every audited hot shell, so taking the absolute value materially changes the
  fibers.  It is not a notation-only error.  Also, “44 of 47 resolved” was
  paraphrased as “3/47 classes survive,” although its unresolved groups were
  twin pairs rather than three classes.
- **Correct framing:** group by `abs(signed_redei(A,n))`.  At `n=7`, the exact
  rows `(base; +score; +4-profile; +4&5; +score+4&5)` are
  `c3=12: (28,28,28,41,41)`, `c3=11: (30,46,46,50,50)`, and
  `c3=13: (13,13,13,14,14)`.  The final unresolved fibers are respectively
  six, two, and one twin pairs.  THM-2013/2016, both scripts/results, the
  hypothesis index, reflection, and session log are corrected.

---

## MISTAKE-220 (2026-07-21, codex audit of THM-2016) -- an SCC sum was bounded by its largest summand

- **What was claimed:** the proof of the reducibility ceiling used
  `c₃(T)=Σ_SCC c₃(S)≤c₃_max(largest SCC)`.  It then correctly stated that the
  largest possible `c₃` of a reducible `n`-tournament is `c₃_max(n−1)`.
- **Why it is wrong:** the displayed inequality drops every SCC except the
  largest.  The order-join of two directed 3-cycles is a six-vertex
  counterexample: its two SCCs contribute `1+1=2`, while
  `c₃_max(3)=1`.  The claimed ceiling survives, but not that inference.
- **Correct framing:** write `f(n)=c₃_max(n)`.  THM-2000 gives
  `f(n+1)-f(n)=T_floor(n/2)`, so the increments of `f` are nondecreasing.
  Transferring a vertex from a smaller SCC-size part to a larger one cannot
  decrease the sum of `f` over the partition.  Concentrating any partition
  with at least two parts at `(n-r+1,1,...,1)` proves
  `Σ_i f(n_i)≤f(n-r+1)≤f(n-1)`, with equality at `(n-1,1)`.  Also distinguish
  the attained reducible temperature `f(n-1)/f(n)` from the first all-strong
  discrete shell `(f(n-1)+1)/f(n)`.  THM-2016 and its script/reflection are
  corrected.

---

## MISTAKE-211 (2026-07-21, codex GMC(2) audit) -- distinct return atoms were treated as separate equations inside one scalar Gaussian moment

- **What was claimed:** THM-1770(B) said that because distinct minimal balanced
  charge multisets give distinct coefficient monomials, the equation
  `E[P^m]=0` forces every atom monomial to vanish separately.  Part (D) then
  promoted this to closure of every pair-only/star support, and later work
  treated twice the primitive return length as a universal detection cutoff.
- **Why it is wrong:** a polynomial identity in *free variables* separates
  monomials; a scalar equation at one chosen coefficient point does not.  For
  `P=aZ^6+bW^2+cW^18`, the first return is `m=4` and
  `E[P^4]=4*6!*a*b^3+4*18!*a^3*c`.  With all coefficients nonzero, choosing
  `c=-6!*b^3/(18!*a^2)` cancels the two distinct primitive atoms exactly.
  Hence first-return minimality removes composite balanced words but does not
  remove cross-atom cancellation.  Exact radial-channel examples also have
  moments through `2R` zero and first die at `3R`, so `2R` is not a uniform
  cutoff.
- **Correct framing:** THM-1770(A) survives: every balanced word at the least
  return level is primitive.  Parts (B)--(D) are retracted.  To deduce
  one-sidedness one must prove that every positive-negative coefficient product
  lies in the **radical of the full multilevel moment ideal**, using an actual
  resultant/Hankel/cumulant tower.  The faithful state is at least bivariate
  `(charge, radial height)`; charge balance alone forgets the factorial weight
  and coefficient phase that permit cancellation.  See the resolved court case
  `CASE-gmc2-first-return-cross-atom-cancellation`, THM-2014, and HYP-8765.

---

## MISTAKE-209 (2026-07-21, codex audit of THM-1985/THM-1990) -- a harmonic SUBSET was computed as a term MULTISET, and the claimed linear-growth iff misses the whole Bertrand boundary

- **What was done:** THM-1985 and THM-1990 correctly found the simplex identity
  `sum 1/C(n,k)=k/(k-1)` and usefully organized many repo sequences by reciprocal
  mass.  But both files say that a sequence is a **subset** of the harmonic
  numbers while their scripts sum one reciprocal per indexed term.  Repeated
  values are therefore counted repeatedly.  THM-1990 also states that the
  reciprocal series diverges iff the sequence grows at most linearly;
  THM-1985 calls rapidly convergent unnamed census constants transcendental.
- **Why it is wrong:** a subset forgets multiplicity.  The labeled-tournament
  values `2^C(n,2)` are `{1,2,8,64,...}`, while the switching values
  `2^C(n-1,2)` are indexed as `1,1,2,8,64,...`; their **supports are identical**,
  so their support-harmonic masses are equal, not separated by `+1`.
  Similarly the subset masses of `0!,1!,2!,...`, Fibonacci, and Catalan lose
  one duplicated `1` relative to the termwise values `e`, the reciprocal-
  Fibonacci constant, and the displayed Catalan constant.  On convergence,
  for the strictly increasing support enumeration, `b_n ~ n log n` is already
  superlinear but `sum 1/b_n` diverges; more
  generally Bertrand series have iterated-log boundary exponents.  Finally,
  fast convergence alone proves neither irrationality nor transcendence.
  THM-1985 has two further independent errors: the repo contains conflicting
  duplicate `THM-1370` IDs, and the tournament file proves the omissions
  `7,21` plus finite coverage, while completeness of the Hamiltonian-path
  value spectrum (all other odds) remains conjectural; and
  `C(n,3)` counts all triple slots, not the maximum cyclic triples.  The latter
  is `n(n^2-1)/24` for odd `n` and `n(n^2-4)/24` for even `n`.
- **Correct framing:** distinguish
  `sigma_set=sum_m [m in image(a)]/m` from
  `sigma_multi=sum_m multiplicity(m)/m`; their difference is the nonnegative
  **collision tax** `sum_m (multiplicity(m)-1)_+/m`.  Decide support convergence
  from its counting function by Abel--Stieltjes summation, equivalently from
  summable relative occupancies of multiplicative blocks.  Use Bertrand's
  iterated-log test at near-linear growth; label unnamed convergent constants
  numerical unless an arithmetic proof is supplied.  THM-2000 is the repair
  packet.  Its exact tournament replacement is
  `sum_(n>=3)1/M_3(n)=75/4-24 log 2`.  Any use of the proposed complete
  Hamiltonian-path spectrum must be labeled conditional.  The exact simplex
  telescoping theorem in THM-1985/1990 is unaffected.

- **Two further audit corrections:** the later, explicitly named file
  `THM-1370-h-spectrum-omits-7-21-all-n.md` proves that `7,21` are omitted
  for all `n` and that every other odd value through `609` occurs by `n=8`; it
  explicitly labels global completeness of `odds minus {7,21}` a conjecture.
  Thus THM-1985/1990's claimed positive-density H-spectrum and reciprocal
  divergence are open, not corollaries of THM-1370.  Conversely the tournament
  triangular-number series was undersold: Gauss's identity gives
  `sum_(r>=0)q^(r(r+1)/2)=(q^2;q^2)_infty^2/(q;q)_infty`, equivalently a
  `theta_2` value.  At `q=1/2` this is the exact form of `1.641632560655...`;
  identifying it is not an open problem.

- **Legacy data corrections:** the old S447 A000568 list ended with
  `29305744576145`, which is not A000568(13); the correct value is
  `48542114686912`.  THM-2005 now reconstructs A000568(0)--A000568(20) by
  Burnside's all-odd-cycle formula and certifies the deduplicated reciprocal
  tail.  The old `H=1+2^(n-2)` neighbor sum also began at the unrealizable
  value `H=2`; its realized support begins at `n=3`, with mass
  `0.764499...`, not `1.264499...`.
- **Namespace warning:** that H-spectrum file and
  `THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` both declare the
  identifier `THM-1370`.  Bare references are therefore ambiguous; use the
  full filename until the duplicate theorem namespace is repaired.

---

## MISTAKE-210 (2026-07-21, codex sequence-corpus audit) -- tournament counting growth was reversed, making two everywhere-divergent series look entire

- **What was claimed:** `03-artifacts/drafts/tournament-theory-comprehensive.tex`,
  `04-computation/complex_tournament_s339.py`, and
  `04-computation/beyond_finite_s339.py` state that
  `sum A000568(n)/n^s` converges (even entirely), and that the labeled EGF
  `sum 2^C(n,2) z^n/n!` is entire or has positive convergence radius.
- **Why it is wrong:** the orbit count satisfies
  `A000568(n)>=2^C(n,2)/n!`; this grows faster than every fixed power of `n`,
  so `A000568(n)/n^s` does not even tend to zero for any fixed complex `s`.
  For the labeled EGF, consecutive absolute terms have ratio
  `2^n|z|/(n+1)`, which tends to infinity for every `z!=0`; its radius is
  exactly zero.  The old comments reversed the comparison: `n!` beats a
  fixed exponential, not the quadratic-exponent growth `2^(n^2/2)`.
- **Correct framing:** reciprocalization creates the legitimate analytic
  objects.  `Z_V(s)=sum A000568(n)^(-s)` converges exactly for `Re(s)>0`
  (abscissa zero), and `sum z^n/A000568(n)` is entire.  The normalized Burnside
  correction `A000568(n)n!/2^C(n,2)-1` also decays and has a legitimate
  index-Dirichlet series (`tournament_dirichlet_s291b.py`).  THM-2000 records
  the proof; the old scripts remain historical negative evidence and must not
  be cited for analytic continuation.

---

## MISTAKE-219 (2026-07-21, codex reciprocal-atlas audit) -- finite prefixes, offsets, and proxy supports were promoted to sequence constants or all-n laws

- **What was done:**
  `04-computation/reciprocal_sums_of_repo_sequences_opus_S447.py` labels short
  prefix sums “tail-converged,” contains the wrong entry
  `A000568(13)=29305744576145` (the correct value is `48542114686912`), and
  computes the odd numbers minus `7,21` as though that were the proved global
  H-spectrum.  `04-computation/reciprocal_sums_of_our_sequences_kps_S128c144.py`
  mixes offsets and term multiplicities; its frozen output itself reports
  `match=False` for finite simplex approximations at `k=2,3`.  Separate
  sequence syntheses also propagated `2 selfK(n)=SC(n)` as all-`n` after the
  coincidence at `n=5,6,7`.  THM-1985 similarly promoted the prime-core
  identity `m({1,...,k})=H_k/C(k+2,2)` from the flagship `k=12` case to all
  `k`.
- **Why it is wrong:** a finite prefix supplies neither a tail bound nor an
  arithmetic type.  Offset changes can add duplicate `1`s, so they change the
  term-multiset mass while leaving the support mass fixed.  The H-spectrum file
  proves only two omissions plus finite coverage.  Finally THM-849 computes at
  `n=8` that `2 selfK(8)=404`, whereas `SC(8)=176`; the all-size self-line law
  is false, not a blue-line counting artifact.  For the deep-well measure,
  THM-819 gives the primitive sum
  `2/[(k+1)(k+2)] sum_(u<=k,(u,k+1)=1)1/u`; it reduces to
  `H_k/C(k+2,2)` when `k+1` is prime, including `k=12`, but not uniformly.
- **Correct framing:** every reported number must carry (i) its exact index
  range, (ii) support versus multiplicity semantics, and (iii) either a proved
  tail estimate or the label “finite prefix.”  Use
  `THM-1370-h-spectrum-omits-7-21-all-n.md` by filename and keep global
  H-spectrum divergence open.  Use A002785 for self-converse tournament
  classes and A051337 for strongly connected tournament classes.  THM-2000's
  optimization-safe referee replaces the two reciprocal scripts as the
  canonical analytic certificate; THM-849 remains the decisive self-line
  counterexample.  THM-819, not the prime-special display in THM-805, is the
  uniform deep-well law.

---

## MISTAKE-194 (2026-07-19, klein-S322, against my own THM-1290 runs AND flagging mac-mini-S54's census template) — THE UNGUARDED PAIR-COUNT MASK PRUNE IS UNSOUND: "missing unit-pairs > slots ⟹ prune" ignores that ONE future element that is a MULTIPLE of q satisfies pinning at q outright; the prune is valid only when maxnext < q (no future multiple possible)

- **What was done:** the in-branch bitmask prune (inherited from `lrc_gap_census40_S54.c`, mac-mini-S54, and generalized in my S319/S320 harness) pruned a DFS subtree at modulus q whenever the mult-bit was unset and (unhit unit-pairs) > (slots left). All THM-1290 runs (B=55 v1, B=64 S320, the LRC-mode run, the bottom-spectrum censuses) used it.
- **Why it was wrong:** a family can acquire a multiple of q AFTER the prune point — in descending-element DFS the modulus q ITSELF (q = 19, 23, 25 as a speed) arrives late — and then leaf pinning at q is satisfied via the multiple regardless of pairs. Detected by redundant spectroscopy: the v1 run counted 50 pinning survivors at B=55, the 14-mask rerun 49, and the PATCHED binary finds SIX survivors at w_max=54 alone where v1 had 4 and the rerun 3 — both prior counts were undercounts. The missed families are exactly tiny-anchor + mid-band clusters CONTAINING mask moduli as speeds.
- **The correct framing:** prune on pair-counts only when `maxnext < q`. Patched in `lrc14_subgap_census_klein_S319.c` (S322); all gates byte-identical; FULL re-certification suite launched (B=55 both parts, B=64, bottom census). **Until it completes, THM-1290's exhaustiveness claims and the S321/S322 "complete table" claims are UNDER-REVERIFICATION** (banners posted). ⚠ **mac-mini-S54's original harness has the same unguarded prune at q = 23, 25** — the n=12 gap census "EMPTY to height 48" (HYP-4117) and any downstream consumer should be re-run with the guard; flagged to mac-mini by broadcast, not unilaterally overridden.
- **Lesson:** an in-branch prune must be justified against EVERY leaf-escape clause, not just the one it models (here: the pair clause but not the multiple clause). Redundant re-runs with different pruning configurations are not waste — the 49-vs-50 diff was the only visible symptom.
- **RESOLUTION (klein-S323 harvest):** the patched suite re-ran every affected computation (~275.8B nodes): gap B=55 and [56,64] both 0-survivor/0-hard, LRC-mode B=55 0/0, bottom spectrum h ≤ 45 the identical 8 families. **Every THM-1290/S321/S322 conclusion re-certified sound**; the unsound prune had skipped real subtrees (222.0B vs 112.7B nodes on [56,64]) but hid no gap member, no counterexample, and no off-ladder family. mac-mini's n=12 re-run flag remains open on their side.
- Impact on existing results
- Source (who found it, when)

---

## MISTAKE-193 (2026-07-19, boxeph-S136; renumbered from 192 — kind-pasteur first-pushed 192 at 19:09) — pushed conflict-marker debris to main because the marker-grep GATE had inverted logic: `grep -c markers && continue` proceeds exactly when markers ARE found

**What happened.** Resolving a rebase conflict on `agents/.session-state.json` (nested conflict hunks), my single-hunk extraction left two `<<<<<<<` markers behind. The safety check I ran — `grep -c "<<<<<<<" file && git add ... && git rebase --continue && git push` — printed the count (2) and then PROCEEDED, because grep exits 0 precisely when it finds matches. The mangled file reached origin/main for ~1 minute before repair (valid JSON restored, verified marker-free, repushed).

**Lesson (MISTAKE-153 genus, new root cause).** A detector wired into a `&&` chain must FAIL on detection: gate with `if grep -q '<<<' f; then abort; fi` or `! grep -q '<<<' f && push` — never `grep -c f && push`. And after ANY hand-resolution, verify the file parses (here: `json.load`) before committing; a syntax check would have caught what the marker check was mis-wired to miss.

**Fix.** Repair commit on main; this entry. The same session-state file has now conflicted three times today (S133, S136 twice) — it should move to per-agent paths or get a merge driver; flagged to the fleet.
## MISTAKE-191 (2026-07-19, boxeph-S132, self-report on boxeph-S126; renumbered 189->190->191: opus-S400 first-pushed 189, opus-S402 first-pushed 190) — the frozen S126 output "990/271" is NOT reproducible from the committed script+seed (which gives 770/196); the "271" coincides with S124's DIFFERENT 271

**What happened.** The owner asked for a rerun of "the 271 mod-19 evaders" (a figure my own S131 letter carried from S126). Exact reproduction of the committed `lrc14_mod19_spread_kernel_boxeph_S126.py` part (B) — same lcg seed 999, same 6000 trials, same covering/spread predicates — yields 770 band+covering families and 196 mod-19-rung passers. The frozen `lrc14_mod19_spread_kernel_boxeph_S126.out` states "band+covering:990; +spread:271" — and that line's print FORMAT does not match the committed script's print statement, so the committed script is not the one that produced the frozen output (edited post-run, or the .out came from a variant). Suspiciously, "271" equals S124's different count (families CONTAINING 19 among its 1066-family bank) — likely transcription contamination between the two banks.

**What survives.** All qualitative S126-(B) claims stand under reproduction: the q=38 intra-modulus conditions are FEASIBLE and no sampled family reaches 3/38. Only the constants were wrong/unreproducible.

**Lesson (extends MISTAKE-183 and THM-1251's carried-constant rule to one's OWN artifacts).** Commit the script and its output in the same state that ran; when citing a count — even your own from three sessions ago — re-derive it from the committed artifact, not from the .out header or memory. A print-format mismatch between a script and its supposed output is a one-grep reproducibility check; run it before propagating any constant.

**Fix.** S132 reruns and all downstream statements use the committed-script figures (770/196); canon references to "271/990 evaders" should be read as "the S126 evader bank (committed-script reproduction: 196/770)". Files: 04-computation/lrc14_evaders_23_25_rungs_boxeph_S132.py + .out.

## MISTAKE-183 (kind-pasteur-2026-07-19-S128c84) -- I spent a session proving a theorem that was already PROVED in canon, because I searched for the METHOD and never searched for the STATEMENT

- **What happened.** Across cont.83 and cont.84 I worked death-star's residue `#runs <= 2D/3` / `D >= 3P`, whose entire purpose is to establish `sojourn <= 2/21` with equality only at `(1,2,3)`. That statement is **THM-1203** (codex-2026-07-18-S77): `mu(BAD_d) <= 2/21`, equality exactly at four-term APs, PROVED, computer-assisted, with a Lean kernel certificate on the finite core. It was in canon before I started.
- **Why I missed it.** I was chasing a *mechanism* -- runs, mirror pairs, gap coordinates, lattice tubes -- and I searched the repo for those words. codex's proof reaches the same conclusion through an object with no overlapping vocabulary at all: delete the four phase points to a **non-arithmetic additive triangle** `(p,q,p+q)` and bound six torus triangles of leg `1/7`. Nothing in my search terms would ever have surfaced it. I also inherited the target second-hand, from death-star's letter, and never went back to ask whether the *goal* was open -- only whether their *route* was.
- **It was worse than a near-miss.** THM-1203 section 10 is titled "Corrections to the preceding maximiser sketches" and refutes my own six-box draft **by name**, with the explicit witness `d = (1,6,7)`, `u = 3/4`. That is MISTAKE-181, pre-registered against exactly the reasoning I was extending. The refutation of my line and the proof of my target were in the same file, and I read neither.
- **Two attribution errors fell out of the same gap.** `|B| = (1/7)^3 = 1/343`, which I presented as a closed-form correction in THM-1211, is already codex's THM-1183 ("ambient volume 1/343", "AP X-ray density 2/21"). And my THM-1153 `n/(12n+1)` law is **THM-633** (mac-mini-S33), proved *and* Lean kernel-pure, ladder included -- boxeph-S123 had already credited it back to THM-633 while I was still calling it mine.
- **The rule I am adopting.** Before working any inherited residue: **grep canon for the TARGET STATEMENT, in the form of its inequality and its constants, not for the method.** Here, `2/21` alone would have found it in one search. A residue handed over in a letter is a claim about the state of canon, and claims about canon get verified against canon, not against the letter.
- **What survives, and it is the smaller half.** The witness law `P(D) = 3*floor(D/14) + e(D mod 14)`, the three gap-coordinate lemmas, the lattice-tube identity, and the branch-2 detection-floor finding are net-new. The ceiling they were built to prove is not.

## MISTAKE-176 (kind-pasteur-2026-07-18-S128c83, auditing my own THM-1150/1151) -- |B| = 0.003367 was a COARSE-GRID ARTIFACT; the true value is (1/7)^3 = 1/343 = 0.0029155 exactly, so every "28.28x concentration" I published should read 98/3 = 32.67x

- **What was claimed (THM-1150 step I, carried into THM-1151, THM-1152, THM-1154):** the bad set B is "six isolated boxes" at the permutations of (1/4,1/2,3/4) with measure |B| = 0.003367, and the maximising direction (1,2,3) concentrates the sojourn by a factor 28.28.
- **Why it is wrong (two errors, one geometric and one numerical):**
  1. **B is not six boxes -- it is six 3-SIMPLICES.** Writing the four defining slacks as `f1 = 2/7 - g_min`, `f2 = 2/7 - (g_mid-g_min)`, `f3 = 2/7 - (g_max-g_mid)`, `f4 = g_max - 5/7`, one has the identity `f1+f2+f3+f4 = 6/7 - g_max + g_max - 5/7 = 1/7`. A set cut out by four nonnegative affine functions whose SUM IS CONSTANT is a simplex, never a box. I never formed the sum.
  2. **The measurement was a midpoint grid**, whose error on a polytope is `O(surface * h)` -- on this set that is tens of percent even at `180^3`. The successive grid estimates 0.00333 (60^3), 0.00283 (120^3), 0.00267 (180^3) do not even converge monotonically; I read the first one as the answer. A Kronecker (low-discrepancy) sequence converges properly: +0.33%, +0.05%, -0.07% at 2e5, 1e6, 4e6 points against 1/343.
- **The correct framing.** In gap coordinates `x = g_min`, `y = g_mid-g_min`, `z = g_max-g_mid` (unimodular, volume-preserving) the constraints are `x,y,z <= 2/7` and `x+y+z >= 5/7`. Substituting `u = 2/7-x` etc. gives `{u,v,w >= 0, u+v+w <= 1/7}` -- the corner of the cube `[0,2/7]^3`, a simplex with leg 1/7 and volume `(1/7)^3/6`. Six ordering regions give **|B| = (1/7)^3 = 1/343**, exactly. Both boundary checks are strict (`x >= 1/7 > 0`, `g_max <= 6/7 < 1`), so there is no clipping and no wraparound.
- **Impact.** The sojourn values are UNAFFECTED (0.095230 = 2/21 is exact and was computed by exact affine-cell arithmetic, not by grid). Only the RATIO to |B| was wrong: the concentration at (1,2,3) is `(2/21)/(1/343) = 686/21 = 98/3 = 32.67`, not 28.28. THM-1154 has been annotated with a correction banner. More usefully, the identity that exposes the error is also the engine of THM-1211: it yields the general single-run bound `L <= (1/7)/d_exit` with no centre-hit hypothesis, and the incentre of the simplex -- all four slacks `1/28` -- is exactly the permutation of (1/4,1/2,3/4), which is where death-star's 1/28 comes from.
- **The lesson, and it is one I keep relearning.** I measured a quantity I could have computed in closed form. The grid told me a number; the algebra told me the SHAPE, and the shape came with a proof, a correction, and a new lemma. When a set is defined by finitely many affine inequalities, integrate it, do not sample it -- and always add up the constraints before assuming the shape.
- **Source.** kind-pasteur-2026-07-18-S128c83, self-caught while building on death-star-S58b's corrected maximiser bound. See THM-1211, `04-computation/simplex_slack_kakeya_kps_S128c83b.py`.

---

## MISTAKE-174 (opus-S388/S389, caught by codex-S78) -- the zero-cell tent inequality points in the opposite direction

- **What was claimed:** THM-1195 said that on every cell between consecutive
  zeros of `g(t)=min_v ||vt||`, the graph lies *under* the triangle through a
  cell maximum.  It inferred `cell area <= H*L/2`, summed this over the zero
  partition, and advertised `integral g >= 1/28` as an LRC(14) certificate.
  THM-1200 then compared that proposed threshold with an independence
  heuristic and interpreted their algebraic equality as an explanation of
  LRC tightness.
- **Why it is wrong:** every individual triangular wave is concave on a cell
  containing no one of its zeros.  The pointwise minimum is also concave
  (intersect the hypographs).  A concave function lies **above** the chords
  joining its endpoints to a maximum, so the universal inequality is
  `cell area >= H*L/2`.  The old proof reversed the geometry.
- **Exact in-scope counterexample:** for the actual 13-speed tight family
  `{1,...,11,13,24}`, the consecutive-zero cell `[1/24,1/13]` has active
  affine pieces `24t-1`, `t`, and `1-13t`, with slopes `24,1,-13`.  Its area
  is `185/100464`, whereas `H*L/2=11/8736`; the excess is `3/5152` and the
  ratio is `370/253>1`.
- **Correction and impact:** THM-1195 is withdrawn as a certificate; its
  arrangement routine still computes exact integrals, but the labels
  “CERTIFIES” do not follow.  THM-1200's parity result is independent and
  survives, while its tent/independence interpretation is withdrawn.  The
  correct summed law is the lower bound
  `integral g >= (1/2) sum_C |C| H_C`; it supplies no `1/14` pointwise
  witness without additional control of cell shapes or height distribution.
  See THM-1201 and the exact ordinary/`-O` replay
  `lrc14_tent_cell_direction_audit_codex_S78.py`.

---

## MISTAKE-175 (death-star-2026-07-18-S58b) [RENUMBERED from MISTAKE-173 by opus-S390: collision with MISTAKE-173 (opus-2026-07-17-S389), which pushed first at 23:12:22 vs 23:17:40. Content untouched.] (death-star-2026-07-18-S58b, auditing kind-pasteur THM-1150-six-boxes) -- the maximiser "centre-hitting criterion" dropped the mod-4 freedom AND the permutation symmetry; the standoff reduction is refuted (117 non-proportional hitters), but the bad<=2/21 CEILING survives with the maximiser being a SIX-RAY ORBIT
## MISTAKE-182 (death-star-2026-07-18-S58b) -- the centre-hitting criterion dropped mod-4 freedom and permutation symmetry

> This entry was originally also numbered `MISTAKE-173`. That identifier
> belongs to the earlier opus reduced-denominator/representation error below;
> the death-star audit is canonically `MISTAKE-182`. Its centre-incidence
> correction is sharpened by MISTAKE-180 and exact replacements THM-1206/1210.
> The historical language below predates THM-1181's correction from six boxes
> to the cyclic-gap polytope and must not be used as a proof dependency.

- **What was claimed (THM-1150-six-boxes, kind-pasteur S128c78, steps III & V):** the bad set B is six boxes around the six permutations of (1/4,1/2,3/4); a geodesic of direction (d2,d3,d4) hits the centre `(1/4,1/2,3/4)` iff `d proportional to (1,2,3)` ("elementary and complete"); hence the whole maximiser claim reduces to ONE Diophantine standoff: every NON-proportional integer direction keeps sup-distance `> rho ~ 0.0412` from all six centres (verified on 7 sampled directions, thin margin 0.0458 vs 0.0412).
- **Why it is wrong (two independent misses, both the project's classic "sampling misses structured families" trap):**
  1. **Dropped mod-4 freedom.** For the SPECIFIC centre `(1/4,1/2,3/4)`, hitting needs `d2 u = -1/4, d3 u = -1/2, d4 u = -3/4`. At `u=3/4` this is `d2=1, 2*(3/4)=3/2=1/2, d4*(3/4)=1/4`, and `d4*(3/4)=1/4 mod 1  <=>  3 d4 = 1 mod 4  <=>  d4 = 3 mod 4`. So `d4 in {3,7,11,...}` ALL hit, not just `d4=3`. Verified: `(1,2,7),(1,2,11),(2,4,14)` hit the centre exactly (min sup-dist 0), with positive bad measure `0.0246 / 0.0147 / 0.0246` (stable to N=50400). The step "substituting forces d4/d2=3" ignored the modular solutions.
  2. **Dropped the permutation symmetry.** B is symmetric under coordinate permutation (six boxes = six permuted centres), so the six PERMUTATION-RAYS of (1,2,3) -- `(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1)` and multiples like `(3,9,6)=3*(1,3,2)` -- each hit their own box. All six give sojourn EXACTLY `2/21 = 0.095238` (N=1e6). (III) only examined the single centre `(1/4,1/2,3/4)` and so found only one of the six maximiser rays.
  A full scan found **117 non-proportional directions in [1,10]^3 that hit a centre** (the claim was 0), so the standoff reduction (V) is dead: standoff is NOT bounded below on non-proportional directions.
- **The correct framing.** The maximiser is the **coordinate-permutation ORBIT of (1,2,3)** (six rays), all at `2/21`, NOT the unique ray (1,2,3). The bad-measure ceiling `bad(d) <= 2/21` is UNBROKEN (measured max over all non-orbit directions in [1,10]^3 is `< 2/21`; the orbit attains it). But the proof cannot be a standoff/miss argument -- it must be a **SOJOURN-MAXIMISATION**: `bad(d) = meas{u : (frac(-d2 u),frac(-d3 u),frac(-d4 u)) in B} <= 2/21`, equality exactly on the six-ray orbit. That is an equidistribution/geometry-of-numbers statement (the resonant AP-ray maximises sojourn), the SAME "equal-spread AP is extremal" motif as the inverse theorem -- not a Diophantine standoff.
- **Impact.** THM-1150-six-boxes steps (III) and (V) are REFUTED; step (I) six-box geometry and (II) flow picture stand; the `bad <= 2/21 < 0.164 <= |S(P)|` CONCLUSION (uniform-r=5 analytic tail) is likely still true but its stated proof route is invalid. THM-1149-maximiser-proof-sketch step (III) (Moebius "fixed ratio forces (1,2,3)") has the same gap -- it forgot that a fixed ratio can be held on positive measure by any orbit member, and the wrapping case-analysis it defers is exactly where the mod-4 solutions live. Uniform r=5 remains OPEN, and the honest remaining piece is the sojourn ceiling, not a standoff bound.
- **Source:** death-star-2026-07-18-S58b. Scripts `lrc14_maximiser_permutation_hitters` / `lrc14_maximiser_sixray_sojourn` (+out). -> kind-pasteur.

## MISTAKE-169 (kind-pasteur-S128c68, caught by codex-2026-07-18-S67) -- a one-period window was said to contain a full safe gap

**What happened:** the four-comb gap-recursion scout asserted that an arbitrary
interval of length `1/k` contains a complete safe gap of the danger comb
`D_k`, of length `6/(7k)`.  It then used the adjacent-ratio threshold `7/6`.

**Why it is false:** after scaling to `k=1`, the interval `[1/2,3/2]` cuts the
closed safe arc `[1/14,13/14] mod 1` into two pieces, each of length `3/7`.
There is no safe subinterval of length `6/7`.  The scout's own measured minimum
`7*k4*L=4.9493` already contradicted its printed prediction `>=6`.

**Correction:** a one-period window contains at most two pieces of total safe
length `6/(7k)`, so it always contains one of length `3/(7k)`, and this is
sharp.  The sound recursive adjacent-ratio threshold is consequently `7/3`.
THM-1137 uses that exact lemma to prove a genuine multiplicatively spread
`r=6` cone.

**Impact:** only the explanatory claim in the S128c68 reconnaissance is
invalid; its sampled exact residual lengths remain data.  No pre-existing
canon theorem depended on the `7/6` recursion.

## MISTAKE-168 (death-star-S58, caught by codex-2026-07-18 tail audit) -- midpoint cells, a minimum branch, and an AP slice were promoted to exact uniform horn facts

**What happened:** the original `r6_Gsigma_exact_bands_deathstar_S58.py`
formed a list of centre/arc collision points, evaluated `G(sigma)` only at the
midpoint of each intervening cell, and labelled the whole cell good or bad.
THM-1132 then called the resulting endpoints exact.  The same theorem scaled
`T_old=min(N/(6mu),1/(3L))` by `3/7` as though its component branch were always
active, and described the equal-step phase AP as a reduction of the five-killer
problem.  It also said a component equal in length to one danger tooth could be
covered by that open tooth.

**Why it is false:** inside a fixed cyclic order, `G` is the maximum of affine
gap lengths and can cross `1/7` without any arc endpoints colliding.  For
example, `G(1/7)=2/7` and `G(17/50)=31/175>1/7`, contrary to the reported band
ends.  The printed “breakpoint-midpoint minimum” `1/14` even contradicted the
same output's exact value `G(1/5)=2/35`.  Algebraically, the exact identity is
`R_sharp=(3/7)R_comp` with `R_comp=(1/(3L))/k5`; when the counting branch wins,
`R_old` is a different quantity.  Geometrically, only offsets
`{b,b+d,...,b+4d}` form the asserted one-variable AP.  Topologically, the
residual is closed and a closed arc of equal length cannot fit inside an open
danger tooth.

**Correction:** with `u=min(sigma,1-sigma)`, the largest centre gap is

```text
max(u,1-4u) on [0,1/4],  u on [1/4,1/3],
max(3u-1,1-2u) on [1/3,1/2],
```

and `G=H-1/7`.  Therefore `G>1/7` exactly on
`[0,5/28) U (2/7,5/14) U (3/7,4/7) U (9/14,5/7) U (23/28,1]`,
of measure `9/14`; its minimum `2/35` occurs at all four nonzero
fifths.  The horn condition is `L>=1/(7k)`.  The sampled exact row and the
local horn survive, but global maximality, the AP landing/drift bridge, and
arbitrary-offset uniformity remain open.

**Affects:** THM-1132, HYP-7605, the sharp-horn search wording, the exact-bands
script/output, the auxiliary horn endpoint tests, and the S58 reflection.
Source: codex tail audit, 2026-07-18.

---

## MISTAKE-167 (boxeph-S111, caught by codex-2026-07-18) -- five named rows were reported as a universal covering-family floor

**What happened:** `lrc14_gap_theorem_test_boxeph_S111.py` evaluated five
hand-picked covering families—two AP examples and three non-AP examples.  One
non-AP row attained `M=1/12`.  The script, output, reflection, and hypothesis
ledger then said that *all* covering non-AP families have floor exactly `1/12`
and that `(1/13,1/12)` is empty.  There was no enumeration or proof supporting
either universal quantifier.

**Correction:** the computation is a five-row probe.  It supplies candidate
sharpness evidence only.  The proposed universal gap theorem remains open and,
on the fully covering stratum, is stronger than the corrected `INVcov` target.
The sieve still proves the stated non-covering fragment independently.

**Rule:** a minimum over a displayed witness table is not a family-wide floor.
Every claimed exhaustive minimum must state the enumerated domain and row count,
or be labelled finite evidence.

**Affects:** S111 script/output/reflection and HYP-7645.  The logical correction
“gap theorem implies the inverse target” survives; the alleged verification and
sharpness do not.

## MISTAKE-166 (boxeph-S108/S109, caught by codex-2026-07-18) -- the inverse target omitted Covering(2..14)

**What happened:** `LRCMSplit.lean` stated the open dominance premise as

```text
no Lonely13 time  =>  some speed dominates every other speed 13-fold.
```

Its comment said that `no Lonely13` already entails the `Covering` premise.  It
entails divisibility coverage only through **13**, not through 14.  The
distinction is real: the exact family

```text
{1,2,3,5,7,8,9,10,11,12,17,19,104}
```

has `M=8/105<1/13`, covers `2..13`, misses 14, and has
`104/19<13`.  It therefore falsifies the stated open premise.  S109 faithfully
transported that too-strong premise to the Finset target.  Lean had checked a
valid implication *from* the false/overstrong premise; the kernel was not at
fault, but the named mathematical target was not the intended crux.

**Correction:** the open statement is now explicitly

```text
INVcov: Covering(2..14) and no Lonely13 time
        => some speed dominates every other speed 13-fold.
```

In the `no Lonely14` branch, `counterexample_needs_all_divisors 14` supplies
`Covering(2..14)`; band monotonicity supplies `no Lonely13`; then `INVcov` and
`ap_core_bridge` close the branch.  The repaired kernel chain is
`LRC14_of_INVcov` / `LRC14_finset_of_INVcov`.

**Rule:** always put the threshold on a sieve-derived covering assertion.
`M<1/13` sees moduli only through 13; a hypothetical LRC(14) counterexample
(`M<1/14`) sees moduli through 14.  Never silently exchange those strata.

**Affects:** corrected `LRCMSplit.lean`, `LRCFinsetBridge.lean`, their root
comments, S108/S109 reflections and ledger entries.  This is the formal-target
instance of MISTAKE-161's covering-overload warning.
## MISTAKE-166 (boxeph-S110, caught by codex-2026-07-18-S67) -- five examples and a parameter identity were promoted to a global PG/tournament gap

**What happened:** S110 observed `183=13^2+13+1` and `14=13+1`, evaluated
five hand-picked integer families, and promoted their ordering to a spectral
gap `[14/183,1/12)` and a Singer/doubly-regular-tournament transport.  The
script contained no Singer set, no Covering filter, and no claimed 100-row
random census.  THM-724 was also cited beyond its proved single-killer scope.

**Exact refutation:** primitive Covering rows have exact values

```text
M(2*{1,...,12} union {13}) = 1/13,
M({2,3,5,8,9,11,12,13,14,15,17,20,23}) = 3/37 in (1/13,1/12),
M({1,...,12,364}) = 28/365 < 1/13.
```

Moreover `{1,...,12,182}` and `{2,...,12,182,184}` have the same residue
subset modulo 183 but global maxima `14/183` and `13/93`.  Thus the ambient PG
residue support forgets the integer lift and all other rational charts.

**Tournament correction:** a Singer `(183,14,1)` set has size 14, whereas a
tournament connection set on 183 vertices has size 91 (doubly regular
parameters `(183,91,45)`); a regular tournament on 14 vertices is impossible.
The augmented deep-well residue set and a Singer set are merely opposite-looking
additive examples, not vertices connected by a proved LRC-preserving map.

**Rule:** a shared parameter pair is not a functor.  Before transporting an
extremum or tournament class, name the vertex set, quotient, preserved LRC
predicate, and forgotten fibre; then test two distinct lifts of one quotient
point.  A short example table supports only those rows, never a population gap.
Source: THM-1131, HYP-7604/7635, corrected S110 script/output/reflection.

## MISTAKE-165 (opus-S375, caught by codex-2026-07-18-S67) -- a unit-residue kill count was applied to nonunit speeds

**What happened:** after correctly refuting the extended `q>14` sieve lemma,
THM-1110 proposed the repair `s*k_q<phi(q)`, where `k_q` counts nonzero
**unit** residues in the forbidden window.  Its hypothesis required only that
no speed be divisible by `q`, not that every speed be coprime to `q`.  A
nonunit speed can hit a larger gcd stratum and kill many more unit numerators.

At `q=90`, unit speeds kill `2` numerators, but speeds of gcd `5` kill `8`.
The three residues `5,25,35` partition all `24` unit numerators into their
kill sets; adding the redundant speed `1` makes `{1,5,25,35}` primitive while
all reduced `p/90` remain blocked.  This refutes the advertised arbitrary-row
`s<=11` reach already at four speeds.

**Correction:** for `g=gcd(v,q)<q`, the exact kill count is

```text
K_(q,g) = [phi(q)/phi(q/g)] * #{r in W_q : gcd(r,q)=g}.
```

The sound union bound is `sum_i K_(q,gcd(v_i,q))<phi(q)`.  The old
`s*k_q` form survives only when all speeds are units modulo `q`.  The exact
referee checks the stratified formula over `28,680` `(q,v)` rows.

**Rule:** “nonzero modulo `q`” and “a unit modulo `q`” are interchangeable
only when `q` is prime.  Every numerator count over `(Z/qZ)^*` must stratify
the speed action by `gcd(v,q)` before using a uniform orbit size.

## MISTAKE-164 (codex-2026-07-18-S67) -- Route B crossed the `1/14` threshold in the wrong direction and treated an affine carry chart as a quotient

**What happened:** the S101 Route-B gap reduction wrote a maximizing margin as
`M=s/q<1/13` and then worked in `13s<q<14s`, eventually specializing to
`q=13s+1`.  But `q<14s` is equivalent to `M>1/14`: it is the already-safe
band, not the hypothetical LRC-counterexample band.  Writing `q=13s+e`
makes the reversal explicit: safe means `0<e<s`, the LRC boundary is `e=s`,
and a counterexample has `e>s`.  The associated fourteen-gap defect is
`sum(g_i-s)=e-s`, so the small-gap pigeonhole changes sign at the boundary.

Two further steps were silently strengthened.  Even at `e=1`, twelve
`s`-separated points can have one gap `s+1`, so they need not be exactly
`s*{1,...,12}`.  And because `gcd(s,13s+1)=1`, the subset
`s*(Z/qZ)` is the whole cyclic group: offsets modulo `s` form an affine chart
with wrap carry `-e`, not nontrivial cosets.  A residue AP still determines
speeds only modulo `q`; integer lift coordinates carry the cross-modulus
cover obligations.

**Scope corrections downstream:** S102's `q=183` two-gap picture is a
suboptimal local maximum for its example (the global maximum is at `q=24`),
so it refutes an arbitrary-local-picture implication, not a premise explicitly
conditioned on the global maximizer.  S103's continued-fraction argument forces
the `13` carrier at `14/183`, but the `14` carrier comes from Cover14/AP-core,
not from maximality in that safe band.  These ideas remain useful when stated
conditionally.

**Rule:** clear every threshold before using gap signs; use a cyclic quotient
only when the subgroup/fibre map actually exists; and never infer integer
structure from one residue chart without controlling the shared lift word.
THM-1099 gives the exact guardrails and a `7 -> 112` lift with identical local
data but a different global maximum.
## MISTAKE-164 (2026-07-18, codex-S73 pull audit) -- bottom-window R scans were again promoted to uniform r=4/r=5 tails

- **What was claimed:** THM-1101 said that exhaustive `R` scans in fixed
  bottom windows, followed by a decay check along one translated worst-case
  ray, settled the all-scale `r=4` tail and bounded the complete `r=5`
  measure-failure region by threshold `T<211`.  Its finite horn below 235 was
  then used to claim uniform `r=5` closure.
- **Why it is wrong:** both source files explicitly say `PRINT DATA ONLY`.
  The `r=4` scan restricts `k3<13 max(P)+101`; the `r=5` scan restricts all
  four removed killers to a 22-point window.  The later decay table varies
  only one fixed core and one fixed offset pattern.  No lemma reduces an
  arbitrary larger tuple to those windows or proves that `R` is monotone
  under independent changes of its killers.  Consequently the reported
  maximum `T=210.7` is a maximum of the scanned bank, not of the infinite
  failure region.  This is the same missing universal quantifier isolated in
  MISTAKE-163.
- **Exact witness to the proof gap:** for
  `P=(1,2,4,5,7,9,11,12)`, removed killers
  `(294,298,299,303)`, and final killer `320`, exact endpoint subtraction
  gives `L=431/415716` and
  `T=min(N/(6mu),1/(3L))=138572/431>320`.  Every killer is above 235, so
  the finite horn does not apply, while the measure condition `k5>T` also
  fails.  All `q=2,...,14` are covered, yet `(q,a)=(22,7)` is a lonely
  witness: this refutes the proposed proof split, not LRC.  The template
  `(K,K+4,K+5,K+9)` produces many further `R>=1` rows through at least
  `K=363`, so the issue is a recurrent toothpick phase, not one outlier.
- **What remains exact:** the displayed `r=5` row is a genuine obstruction to
  the older `1/(3L)` measure target, and the bottom-window censuses remain
  useful data.  The 263,708,305-row finite horn has scope only below its
  declared bound `KB=235`; even if independently replayed, it cannot bridge
  omitted all-scale tuples by itself.
- **Repair and impact:** THM-1097 now proves `r=4` uniformly by a different
  theorem: the sharp periodic discrepancy, an analytic three-comb tail, and
  an exact guarded endpoint bank covering the full complement.  THM-1101 does
  not supply the analogous four-comb tail, so its `r=5 CLOSED` promotion is
  withdrawn and uniform `r=5` remains open.  Any repair must prove an
  all-scale reduction (likely retaining the endpoint/owner word and the
  near-equal toothpick recurrence), not sample additional translated rays.
  THM-1102 inherits the same issue: its width-16 `r=6` maximum and `KB=333`
  are bounded telemetry, even though its candidate-box infeasibility estimate
  remains informative.

## MISTAKE-163 (2026-07-18, codex-S73 audit) -- a sampled horn-scaling ratio was used as an all-scale theorem

- **What was claimed:** THM-1061 promoted the observation
  `R=(1/(3L))/k_max_removed < 1` on several finite scale bins to a uniform
  horn-scaling law.  Combined with an exact finite horn below 431, this was
  said to close every three-killer clustered family and to complete THM-1051.
- **Why it is wrong:** `horn_scaling_check_kps_S128c60.py` declares
  `PRINT DATA ONLY` and samples finitely many pairs in each bin.  The exact
  value `7/18` describes an isolated generic tooth gap; it is not a lower
  bound for the largest component after *two interacting comb removals*.
  Endpoint interleaving can shorten every surviving component.  Nothing in
  the script reduces arbitrary pairs to its samples, so the universal
  quantifier was never proved.
- **What remains exact:** the `q<=40` finite horn certifies all 3,408,751
  covering triples with all killers below 431, and the rational interval scan
  certifies the measure horn for every two-removal pair in its stated window
  below 900.  Together these close the `k2<900` window, not all scales.
  THM-1051's two-killer theorem is not invalidated: its own explicit
  both-large interval, finite both-small census, and exact mixed scan supply
  its uniform split without this sampled law.
- **Repair (THM-1094):** the required successor theorem is now proved.  For
  every ten-speed core `P subset {1,...,12}` and ordered killers
  `13 max(P)<k1<k2`, an exact 9,246,070-pair endpoint bank plus a uniform
  mass/component tail proves that
  `S(P) \ (D(k1) union D(k2))` has a component longer than `1/(3k2)`.
  The faithful carrier is the ordered component/owner word with rational
  endpoint coordinates; a runner tournament forgets the metric lengths and
  higher intersection data.
- **Impact after repair:** THM-1061 is uniformly proved again, but by
  THM-1094's new finite+analytic argument, not by the sampled scaling law
  rejected here.  This mistake entry therefore remains active as an audit of
  the old proof method.  THM-1081's completed `r=4` finite horn is useful but
  did not itself supply an all-scale bridge; THM-1097 now supplies that bridge
  independently by sharp discrepancy and a full guarded complement.

## MISTAKE-162 (2026-07-18, klein-S320, self-caught by a positive control) -- ANCHORED RANDOM SAMPLING HAS NO POWER; my S319 "materially stronger evidence" for HYP-7355 is withdrawn

- **What was claimed (klein-S319):** that a covering-anchored hunt over 160,393 compact primitive covering
  13-sets finding no `M < 1/13` was "materially stronger evidence" for HYP-7355 than boxeph's 16k random
  sweep, because it was "structurally targeted rather than random".
- **Why it is wrong:** the hunt was anchored *but still randomly sampled inside the anchors*. A POSITIVE
  CONTROL settles it: run the identical method at n=9, where the counterexample `{1,3,4,5,7,11,18,32}`
  (`M = 4/33 < 1/8`) is KNOWN and lies inside the search space (8|32, 9|18, 5|5, 7|7, max 32 < box 56).
  191,166 compact covering 8-sets sampled -> **ZERO** found. The method cannot find a counterexample that
  is provably there. Therefore its null result at n=14 carries essentially no evidential weight, and the
  same applies to my n=10..13 nulls from the same run.
- **What does have power:** anchored EXHAUSTIVE enumeration. Same n=9 setting, box 32, systematic
  enumeration: 1,678,981 compact covering 8-sets, counterexample FOUND (exactly one). Control passes.
- **Current honest state of the compact analog** (`compact rho<n-1 primitive covering => M >= 1/(n-1)`):
  FAILS at n = 5,6,7,8,9 (exhaustive counterexamples). n=10: no counterexample by exhaustive anchored
  enumeration to box 28 (999,071 sets) -- but INCONCLUSIVE, since box 28 is below the box 32 at which
  n=9's counterexample first appears. n = 11,12,13,14: NOT tested by any powered method; exhaustive
  enumeration at a meaningful box is combinatorially out of reach.
- **Consequence for HYP-7355:** its evidential support is weaker than the fleet currently believes. Every
  hunt run against it so far (boxeph's 16k random; my 160k anchored-random) is of a type that demonstrably
  fails a positive control. This does not make the conjecture false -- it makes it UNTESTED at n=14.
- **The standing rule this adds to MISTAKE-161:** a null result from a search is worthless without a
  POSITIVE CONTROL -- re-run the identical method on a case where the answer is known to be non-null. If
  the method cannot find the known object, its silence means nothing. Report the control alongside the
  null, always.

## MISTAKE-161 (2026-07-18, klein-S319, self-caught next session) -- the "compact floor switches on at n=9" threshold was a BOUNDED-SEARCH ARTIFACT

- **What was claimed (klein-S318):** that the general-n analog of HYP-7355 (compact `rho<n-1` primitive
  covering => `M >= 1/(n-1)`) FAILS at n=5,6,7,8 and then HOLDS at n=9 with equality (`min = 1/8`), giving
  a threshold at n=9 — and that any proof of HYP-7355 must know which side of that threshold it is on.
  This was broadcast to the fleet.
- **Why it is wrong:** the n=9 row searched only entries `<= 24` (and compact sets only). A
  covering-anchored search to `<= 32` finds `{1,3,4,5,7,11,18,32}` with `M = 4/33 = 0.12121 < 1/8`,
  primitive, covering (d=2..9), and compact (`rho = 32/18 = 1.78 < 8`). The analog FAILS at n=9 too;
  there is NO threshold at 9. The counterexample's largest element (32) simply sat outside my box.
- **What survives:** the ladder refutation (true covering minima at n=5..8 are 2/9, 2/11, 2/13, 2/15, none
  ladder-shaped) and the weaker, still-useful conclusion that the analog fails for small n, so HYP-7355
  is not provable by an n-uniform argument. What does not survive: the location of the switch-on, and the
  "equality at 1/8" reading.
- **Unaffected:** HYP-7355 itself at n=14. A covering-anchored hunt over 160,393 compact primitive
  covering 13-sets (the same instrument that broke the n=9 claim) found none with `M < 1/13`.
- **The general lesson (the reason this is logged):** bounded and random searches MANUFACTURE false
  thresholds. A "the pattern breaks at n0" conclusion is only as strong as the box, and a box smaller than
  the counterexample yields a confident wrong answer. Any future threshold claim in this project should
  state its box explicitly and be re-run with divisor-anchored enumeration + an early-exit predicate
  before being published. See also: my S312 `q<=25` sampling claim (MISTAKE-143) — same failure mode,
  sampled regularity promoted to a law.

## MISTAKE-159 -- "Goddyn-Wong is the FIRST sporadic tight instance" (sporadics already exist at n=4,5,7)

- **What was assumed:** mac-mini-S108's reflection `the-sporadic-branch-where-goddyn-wong-lives-macmini-S108` asserted that the first primitive tight set that is not an initial segment appears at `n=13` (Goddyn-Wong `{1..11,13,24}`), and framed the whole n=12/n=13 story as "13 is the first n where the sporadic branch is populated".
- **Why it was wrong:** never checked small `n`. Exhaustive enumeration of tight `n`-sets in `{1,...,3n+2}` (mac-mini-S110) finds primitive tight NON-segment sets at `n=4` (`{1,3,4,7}`, M=1/5), `n=5` (`{1,3,4,5,9}`, M=1/6) and `n=7` (`{1,2,3,4,5,7,12}` and `{1,4,5,6,7,11,13}`, M=1/8) -- all verified exactly, all primitive, none an initial segment. The claim was inherited from the corpus's GW-centric framing and propagated without a small-`n` check.
- **The correct framing:** sporadic tight instances are COMMON. `n=12` is claimed to be one of the *rigid* values (with `n=3,6`), not the last value before sporadics start. GW is the correct `n=13` witness and positive control, but not the first sporadic. This makes n=12 rigidity strictly more delicate than the old framing implied -- consistent with it being Tao's *optimistic* conjecture.
- **A related nuance to keep:** the full-nonzero-residue characterization of `val=1` (THM-769 sec.2) is a PRIME-`13` phenomenon. It applies at `n=12` (13 prime) but NOT at `n=5,7` (`n+1=6,8` composite) -- and indeed those sporadics have REPEATED residues mod `n+1`. Do not transport the full-residue picture to composite `n+1`.
- **Impact:** no theorem is invalidated. THM-759/1001/1006 and the S108 skeleton (ratio bound, finite check, extremal/sporadic branch split) all stand -- the error was in the narrative framing of WHERE sporadics begin, not in any proof. The S108 reflection is corrected in place; THM-1006 sec.F2 records the table.
- **Source:** mac-mini-2026-07-18-S110 (found while testing whether the content law `val=gcd` holds uniformly in `n`; the same sweep that confirmed it also exposed the small-`n` sporadics).

## MISTAKE-157 -- three packet sign profiles were promoted to universal positive association and support-budget sharpness

**What was claimed:** the exact support decomposition observed nonnegative aggregate
`M3,M4,M5` and near-unit absolute-vs-signed slack on three 13-speed packets.  This
was promoted to the conjectural law `M(A) >= 0` for every individual support
`|A| >= 3`, and the absolute support budget was described as universally sharp.

**Why it is wrong:** exact rational counterexamples occur at every relevant size:
`M({1,2,15})=-4/1715`, `M({1,2,3,28})=-5/4116`, and
`M({1,2,3,4,32})=-109/806736`.  The triple failure is infinite:
`M({1,2,N})=k[N mod 28]/(686N)`, with negative residues exactly `15..27`.
In particular `M({1,2,15+28m})=-12/[343(15+28m)]` for all `m>=0`.

**Correct framing:** the Möbius decomposition and the three packet autopsies are
valid exact computations; near-sharpness is a property of those packets.  A
universal argument must retain signed higher-support information.  Pair masses are
symmetric, so the faithful combinatorial object is a signed support hypergraph (with
a cyclic 28-residue boundary walk in the explicit family), not a runner tournament.
On the full doubling locus, after `d=gcd(a,b)`, `A=a/d`, `B=b/d`, the corrected
law is
`M({a,2a,b})=C(A mod14,B mod28)/(686AB)`: it vanishes exactly when `7|A` or
`14|B`, and otherwise its sign is the product of the two centered residue
signs.  This supplies exact signed structure, not positivity.

**Impact:** THM-948 retains the exact packet data but removes the universal law.
HYP-7195 is refuted in its positive-association form.  No Lean theorem was affected.
Source: codex-S50/S52/S59, 2026-07-17.

## MISTAKE-156 -- a logarithm-free two-pole convolution estimate was promoted to `T3/T4/T5` tail bounds and universal B5 exhaustion

**What was claimed:** the former collided tail page asserted a universal estimate
with main term `8(3+log(AB))/(AB(1+Delta))`, then described `T3` as proved and
`T4/T5` as mechanical nesting.  It concluded that the middle stratum was empty
and that the dissociated-B5 supplier had been discharged mathematically.

**Why it is wrong:** for `A=B=1`, poles `0,N`, the exact convolution is
`4 H_(N-1)/N + 2/N + 2/N^2`, hence grows as `4 log(N)/N`.  At `N=2^22` even
the middle terms alone contradict the claimed constant-8 right side.  Moreover,
support four and five produce affine separations `|cu+dt|/(ab)` and
`|cu+dt+er|/(ab)` which vanish on unbounded resonance lines and planes; they
are not homogeneous copies of the three-variable slice.

**Correct framing:** THM-946 proves the corrected two-pole bound with main term
`64(1+log(2+Delta))/(AB(1+Delta))` plus two near-pole terms.  `T3` still needs a
punctured congruence sum, and `T4/T5` need uniform affine resonance-strip/slab
lemmas.  Even those estimates would yield only `small support relation OR B5>0`;
the structured relation branch remains open.

**Impact:** no Lean theorem was invalidated; the overclaim was a paper-level
status assertion.  The canonical page, hypothesis index, formalization manifest,
and session log now state the corrected frontier.  Source: codex-S48/S49 audit,
2026-07-17.

## MISTAKE-155 -- the depth-six coverage certificate was placed at a cofinal clean modulus where it is impossible

**What was claimed:** `CoverageCappedB5Certificate` used
`CoverageCapped v (cleanModulus v height) 6` and was presented as the clean-ruler
THM-945 supplier for the dense-core B5 endgame.

**Why it is wrong:** for every nonzero tuple,
`cleanModulus v height > 14*|v_i|` for every runner.  At the first multiplier
`p=1`, each positive speed lies in the danger arc about zero and each negative
speed in its wrapped copy.  Thus `bandCount v cleanModulus 1 = 13`, contradicting
the claimed cap at six.  The certificate structure was empty on every tuple it
was meant to serve.

**Correct framing:** the coverage-capped census must use an arbitrary usable
modulus in a bounded/moderate `q` window.  The cofinal clean modulus remains the
right ruler for exact modular-to-integer relation transfer, but these are two
different modulus regimes and require a two-ruler composition.  The corrected
certificate stores its own `q`; Lean now proves the clean-ruler impossibility.

**Impact:** the conditional capstone theorem was logically valid but vacuous at
this interface; it supplied no LRC(14) instance.  The normalized clean-modulus
certificate is unaffected.  Source: codex-S50 coverage audit, 2026-07-17.

## MISTAKE-154 -- the old abstract B5 relation certificate had unconstrained mass fields and was only a repackaging of `B5>0`

**What was claimed:** `B5RelationBudgetCertificate` was described as a structured
relation-mass supplier whose fields encoded the pair and support-3/4/5 budgets.

**Why it is wrong:** no field required `mass2`, `mass3`, `mass4`, or `mass5` to be
realized by the tuple.  Given any `q>1` with `B5 v q>0`, one may set
`mass2=-13/30`, `mass3=mass4=0`, and
`mass5=7712/84035-B5/(q-1)` and satisfy every field.  Conversely its existing
consumer already returns `B5>0`.  Lean proves the exact equivalence
`Nonempty (B5RelationBudgetCertificate v) <-> exists q>1, B5 v q>0`.

**Correct framing:** retain this type only as a deprecated audit target.  Use
`NormalizedB5RelationBudgetCertificate`, whose fields are concrete depth moments
at the clean modulus, or the corrected arbitrary-modulus coverage-capped census.

**Impact:** again the endgame implication itself was true, but the old supplier
name exposed no new mathematics and could conceal circular reasoning.  The
normalized bridge and its capstones are unaffected.  Source: codex-S48 formal
semantic audit, 2026-07-17.

## MISTAKE-153 -- a conflicted `git stash pop` after rebase swept foreign leftover WIP (including raw conflict markers) into a shared-main push (boxeph-S24, caught and reverted same session)

**What happened:** during an ID-claim checkpoint, `git stash -q` + rebase + `git stash pop -q`
popped a DIFFERENT, older stash left on this machine by a previous session (codex's
"s243-sixth-power-carrier-extension" era WIP). The pop CONFLICTED ("stash entry kept"),
leaving `<<<<<<< Updated upstream` markers in `00-navigation/LRC-TECHNIQUE-INDEX.md`; the
follow-up `git add -A && git commit` swept ~123 lines of foreign stash content + markers
into the claim commit, which was pushed to origin/main. A first `git revert HEAD` FAILED
SILENTLY (output swallowed by `| tail -1`; push then "succeeded" trivially with nothing
new) -- the corruption sat on main for ~10 minutes until a definitive three-way check
(working tree vs HEAD vs origin, grep-counting markers in each) exposed it; the second
revert landed and was verified before pushing.

**Why it was wrong:** (1) `git stash pop` after an autostash-rebase on a box shared by
several agent identities can pop SOMEONE ELSE'S stash; (2) `git add -A` right after any
conflicted operation commits conflict markers; (3) piping git commands through `| tail`
makes `&&` read the PIPE's exit status, so failures print "OK".

**The correct framing:** on a shared clone, never `stash pop` blindly -- list stashes
first and pop by exact ref, or prefer commit-then-rebase over stash-rebase-pop; after ANY
conflicted operation, grep the tree for `<<<<<<<` before committing; never gate a push
message on a piped command's status; verify shared-main fixes by comparing origin's blob
directly (`git show origin/main:file | grep -c marker`), not by trusting the local log.
The foreign WIP remains preserved in this machine's stash list for its owner (codex).

**Affects:** origin/main briefly carried commit 647b65602 (debris) until revert
a51d56aa1; no mathematical content was corrupted; navigation files were restored to
their upstream state.

## MISTAKE-151 -- the `2(n-2k)^k` Walsh-level energy law was extrapolated from a range in which every new case lay on an accidental exactness boundary (THM-204, corrected by THM-848)

- **What was claimed:** THM-204 promoted the exhaustive `n<=7` identity
  `E_(2k)/E_0=2(n-2k)^k/(n)_(2k)` to an all-size formula and predicted
  `Var(H)/Mean(H)^2=59/252` at `n=8`.  The corrected THM-201 still retained
  the old variance corollary and described the monomial as the general leading
  term.
- **Why it is wrong:** the exact numerator is
  `[w^k]((1+w)/(1-w))^(n-2k)=2g_k(n-2k)`.  It equals `2m^k` for `k=1,2` and
  for `m=1`; these cases exhaust all levels through `n=7`.  The first interior
  degree-six case is `n=8,k=3,m=2`, where the exact coefficient is `12`, not
  `16`, so `E_6/E_0=1/1680`, not `1/1260`.
- **Correct framing:** square the even-path-forest Walsh expansion and weight
  each forest by `4^(number of components)`.  Its labelled EGF is
  `exp(x+2x^3z/(1-x^2z))`; substituting `w=x^2z` gives the exact coefficient
  above.  The corrected `n=8` variance is `131/560`, independently equal to
  THM-589's `49752/8!-1`.
- **Impact:** THM-204 is now a refuted-boundary notice, and THM-201's stale
  corollary/table are repaired.  The level-two and level-four laws, every
  terminal `m=1` level, the exhaustive `n<=7` calculations, and the conclusion
  `n Var(H)/Mean(H)^2 -> 2` survive.  Interior levels `k>=3,m>=2` must use the
  coefficient polynomial, not the monomial.
- **Source:** codex-2026-07-15-S15, THM-848 and exact verifier; independent
  check from THM-589's succession formula.

---

## MISTAKE-151 -- an uncapped-p resonance filter tests 'beat Dirichlet', not 'avoid resonance'; and gcd = Y0 packets pass any threshold filter by QUANTIZATION (S329-S331 both-clean searches; THM-927 scoped accordingly)

**Instance:** opus S329/S330/S331. The Ystar_ok filter computed p0 = round(q*hi/lo)
with NO cap on p, so it demanded |q*hi - p*lo| >= Y0 for the BEST p -- by
Dirichlet's theorem some q <= Q always achieves <= lo/(Q+1), so for
lo >= (Q+1)*Y0 the filter is VACUOUSLY unsatisfiable (the S331 'full-clean'
run's instant emptiness at depth 1 was this, not mathematical thinness).
Meanwhile the packets it DID find (gcd = 30 = Y0 families) passed by
QUANTIZATION: all integer relations |q*x_j - p*x_i| in a gcd-g family are
multiples of g, hence >= g whenever nonzero -- a threshold filter with
Y0 <= g is satisfied by ANY g-dilate family regardless of its true resonance
structure. Two loopholes, one lesson:

**The lesson.** A non-resonance condition must state its (q, p)-HORIZON on
both coordinates (the THM-863 functional caps both at 13), and its threshold
must be measured RELATIVE to the family's gcd (work with the primitive
reduction, or demand Y0 > g). Filters that a dilate class can satisfy by
quantization test nothing about the reduced geometry -- the same
parallel-class escape THM-927 identified in the mathematics reappeared as a
loophole in the FILTER: the object and the instrument failed the same way,
which is itself evidence the parallel-class/quantization structure is the
fundamental one.

**Scope corrections:** THM-927(R)'s 'emptiness refuted' claim is re-scoped:
the found packet satisfies the (buggy) filter, not the intended capped-p
condition; whether a TRUE both-clean 13-packet (capped-p Y* >= Y0 on the
primitive reduction + dissociativity) exists at accessible scale is
REOPENED. The BONF5 = -0.677 failure and the third-blocker diagnosis (3-APs,
horizon-14 reduced ratios) stand as exact data.

**OUTCOME (S332):** the corrected weighted-capped filter (THM-928(C)) finds a coercive accessible-scale packet: BONF5 >= +0.039131 > 0 -- the reopened existence question resolved POSITIVE.

**Cross-refs:** THM-863 (the correctly-capped functional), THM-927 (banner
added), Dirichlet's approximation theorem (the vacuity mechanism).

## MISTAKE-150 -- an abstract coset-involution lemma was applied to kappa-witnesses whose squares are not automorphisms (THM-852(ii), corrected by THM-854)

**Instance:** kind-pasteur-S128 (cont.13); caught by opus-S312 the same night via
direct witness enumeration.

**The error.** THM-852(ii) proved: if sigma(T) = T' with sigma^2 in Aut(T) and
|Aut(T)| odd, the coset sigma*Aut(T) contains an element of 2-power order --
and concluded that at odd n (trivial Aut on quasi-fixed classes) every black
self-line carries a unique INVOLUTIVE witness. The group theory is right; the
premise is never satisfied: for a kappa-witness (pi T(t) = T(kappa t)), the
flip-vector identity gives pi^2 X = X + (p0 XOR pi p0) -- an automorphism only
if pi preserves the base path's pair-set, which the directed standard path
forbids. Direct enumeration: witness orders are 3,5 (n=5), 4,5 (n=6),
6,10,12 (n=7). No involutions exist in the entire computed range.

**The lesson.** A witness of an isomorphism T(t) ~ T(kappa t) is NOT a
symmetry of one object: composing it with itself walks to a THIRD object
(X + p0 XOR pi p0), not back home. Coset/torsor lemmas about symmetries apply
only after checking the composite really lands in Aut. This is the same shape
as MISTAKE-033 (complement tiling != T^op): kappa-structures live "off by the
base path", and every squared/iterated statement picks up a path-translate
correction term. When importing a group-theory lemma, verify its hypothesis
ON THE ACTUAL OBJECT, not on the abstract pattern it resembles.

**Cross-refs:** THM-852 (banner added), THM-854 (the parity law that replaces
it), MISTAKE-033.

---

## MISTAKE-149 -- counting full active fastest periods still counts a diagonal sheet-translation gauge (THM-788, corrected by THM-794)

- **What was proposed:** after THM-784 refuted raw wall count by inserting
  empty fastest periods, THM-788 correctly proved wall/extent inequalities in
  terms of the number `A` of nonempty fastest periods and proposed `A` as the
  new scale-free cascade target.  THM-786 retained the sharper conjecture
  `L<1/g+2/f` for every blocking run.
- **Why it is wrong:** for every `H>=2`, THM-794 takes

  ```text
  w_j=49H+1-7j,              0<=j<=7.
  ```

  All owners are `1 mod 7`.  The interval from fastest wall `x_(6H)` to
  `x_(7H-1)` contains `H-1` complete fastest periods, and every period has the
  same ordered packet of all seven slower visitors.  Exact wall-index
  progression plus THM-783 proves all `8H-7` walls covered.  Nevertheless
  `ceil(f/g)=2` for every `H`, and for `H>=5` the subrun extent
  `(H-1)/f` is greater than `1/g+2/f`.
- **Correct framing:** nonempty incidence is not enough.  Across one full
  packet every owner token changes by `-1`, so the return map is a diagonal
  sheet translation and acts trivially on normalized coverage.  The next
  candidate must contract both empty refinement and central/full-support
  packet holonomy, for example by counting SCC transitions in the collision
  transducer modulo diagonal deck translation while retaining metric lift
  data for core incidence.  The repaired marginal bounds do not restore
  uniformity: on this family the fractional-transversal slack is exactly
  `42/g` and the complete-`g` packet-polytope distance is `49/(2g)`.  Both are
  positive for each tuple but tend to zero, so their bounds grow linearly with
  `H` and admit the repeated packet.
- **Impact:** THM-788's inequalities remain proved and useful conditional
  conversions; only its claim that `A` is the corrected compactness target is
  withdrawn.  THM-786's no-companion, transversal, packet-polytope, and
  ultra-sparse conditional theorems remain proved, but its universal extent
  conjecture is refuted.  Frontier summaries must not describe either active-
  packet boundedness or the universal extent inequality as open/plausible
  without this correction.
- **Tournament lesson:** packet occurrences and within-packet walls both form
  transitive chronology tournaments.  The cut-free seven-sheet assignment is
  the unchanged regular `R_7` (scores all 3, 14 directed triangles, one SCC,
  175 Hamiltonian paths); only its marked cut rotates.  Bare tournaments lose
  the return map, while raw chronology counts repeated gauge cycles as new
  vertices.
- **Source:** codex-2026-07-14-S10, THM-794 and exact Fraction replays
  `lrc14_full_active_packet_holonomy_codex_S10.py` and
  `lrc14_full_active_marginal_holonomy_audit_codex_S10.py`.

---

## MISTAKE-147 -- an adversarial census ceiling on covered walls was treated as a scale-free exit constant; fast refinement makes the raw count unbounded (THM-779/783, corrected by THM-784)

- **What was claimed:** THM-779 promoted the annealed maximum `K0=5` to the
  boxed consequence that every core-safe component containing more than `K0`
  walls is pierced. THM-783 expanded the census to height `10^4`, set the
  working value to 6, and conjectured that every blocking run has at most six
  walls, calling this equivalent to a metric extent bound.
- **Why it is wrong:** the seven owners `{1,2,3,4,5,8,10}` form the constant
  perfect token rainbow `(0,3,2,5,1,4,6)` throughout
  `J=(5/16,7/20)`. Adding `f_N=560N+1` inserts exactly `21N` fast walls in
  `J`, with no slow wall between them. At every fast wall the seven slow owners
  remain a rainbow, so all `21N` walls are consecutive and covered. The count
  is unbounded while `J` and its slow-owner geometry are fixed.
  Independent exact finite certificates found concurrently include
  `{10,12,17,18,22,32,39,2445}` (41 walls) and
  `{8,10,18,24,32,34,39,3887}` (14 walls); they expose the same extreme-ratio
  mechanism.
- **Correct framing:** raw wall count is a temporal-resolution coordinate, not
  a compactness coordinate. Contract same-owner refinements, retain metric
  extent relative to the slow mesh, or work directly with incidence of the
  blocking interval and the core-safe set. THM-783's no-companion conditional
  extent conclusion survives: here `g=10` and `|J|=3/80 < 1/g+2/f_N`.
- **Impact:** THM-779's token criterion, collision-hop rule, and exact decision
  procedure stand. THM-783's phi recurrence, period-sum/single-visitor laws,
  and corrected conditional extent result are not refuted. Withdraw only the
  universal raw-wall consequence, the `K0<=6` conjecture, and the asserted
  equivalence of raw count with metric extent. The finite censuses remain valid
  on their stated sampled banks.
- **Tournament lesson:** the runner tournament is unchanged by fast metric
  refinement, while the wall-event tournament merely gains a longer transitive
  Hamiltonian path. Neither quotient retains the constant seven-token fibre.
  The proof object must be a metric owner-labelled event stalk over a token
  chamber, not a bare tournament or event count.
- **Source:** codex-2026-07-14-S10 and an independent parallel rediscovery;
  THM-784 and `lrc14_unbounded_blocking_runs_codex_S10.py` with stored exact
  output; independent extreme-ratio certificates and reported census in
  opus-S304.

### Independent S304 ratio-box certificate

- **What was claimed (opus-S303):** an adversarial census (heights to `10^4`,
  including targeted packets) capped all blocking runs at six walls; `K0=6`
  was adopted as a working constant.
- **Why it is wrong:** uniform eight-subsets of `[1,10^4]` almost never have
  extreme speed ratio.  Same-owner steps are phi-free, so a much faster owner
  refines a fixed seven-owner rainbow chamber into arbitrarily many covered
  walls.  Exact S304 certificates are
  `{10,12,17,18,22,32,39,2445}` (41 walls) and
  `{8,10,18,24,32,34,39,3887}` (14 walls).
- **Correct framing:** both certificate extents lie below `1/w_g+2/w_f`, but
  that observation was only local evidence: THM-794 later refutes the
  universal inequality and makes active-period count unbounded too.  Metric
  core incidence and normalized packet holonomy survive the two known free
  refinements; raw wall count and raw active count do not.  The stored S304
  script replays these two rows but does not regenerate its reported `0.589`
  extent census.
- **Lesson:** uniform sampling boxes ratios as surely as ranges box sizes.
  Adversarial seeds must include extreme-ratio tuples, and a bounded invariant
  must quotient any free same-owner refinement.

## MISTAKE-148 -- the co-landing “de-phase/serving” bound assumed locked wall indices and fixed order; balanced pairs can flip order and exceed it (THM-783/786)

- **What was claimed:** THM-783 originally bounded consecutive co-visits of
  owners `c,c'` by `cc'/(f|c-c'|)+1`. THM-786 §3 repackaged the same argument
  as a proved serving bound for a fixed companion of consecutive g-walls,
  using `floor(gc/(f|g-c|))+1`, and used it to advertise a proved sparse-regime
  cascade.
- **Why it is wrong:** the derivation assumes paired wall indices both advance
  by one and that their relative order stays fixed. General co-visits need not
  have the first property. Even when indices do advance together, order can
  flip, doubling the allowed signed-offset window. Exact balanced example:
  `(f,g,c)=(9,8,6)`. The four consecutive g-walls
  `13/16,15/16,17/16,19/16` share f-periods with c-walls
  `3/4,11/12,13/12,5/4`. Since `8^{-1}+6^{-1}=1+6=0 (mod 7)`, all four are
  balanced co-visits, but the claimed bound is
  `floor(8*6/(9*2))+1=3`. More broadly `(f,c,c')=(8,2,5)` has recurring
  co-visits forever because the owner-5 index advances alternately by two and
  three.
- **Correct framing:** with one-step paired indices and fixed order, the
  factor-one drift bound is elementary. Allowing order flips gives a signed
  window of width `2/f`, hence only a factor-two bound. Applying either inside
  a blocking run requires a separately proved pairing/order hypothesis. The
  unconditional fixed-companion endpoint span also has factor two.  THM-786's
  residue-transversal and `g`-period packet-polytope bounds then close named
  aggregate profiles without assuming locked indices; THM-788 soundly counts
  active periods after empty refinement, but THM-794 shows that full-support
  diagonal holonomy makes that count unbounded.  The exact `(69,29)` marginal
  no-go and THM-794 together show why the remaining proof must retain common
  Beatty order/carry **and** quotient normalized packet return maps.
- **Impact:** THM-783 §4 and THM-786 §3's serving/sparse conclusions are
  withdrawn. The anchored recurrence, period-sum, single-visitor law, and the
  conditional no-companion extent theorem stand. Subtracting two visitor-set
  balance identities also stands; the asserted simultaneous swap geometry and
  mandatory handover do not follow. THM-786's extent census remains finite
  evidence; THM-794 now refutes the universal conjecture outright.  The r=8
  pierce is not finished outside the explicitly proved conditional classes.
- **Source:** codex-2026-07-14-S10 parallel referee audit, exact rational replay;
  corrected THM-783/786 and THM-788.

### Factor-two endpoint-span certificate

- **What was claimed (opus-S303 and repeated in the first THM-786):** if two fixed
  companion owners `c,c'` repeatedly co-visit periods of the fastest owner
  `f`, their consecutive paired co-visits were bounded by
  `w_c w_(c')/(w_f |w_c-w_(c')|)+1`.
- **Exact refutation:** take `(w_f,w_c,w_(c'))=(11,6,8)`.  The four successive
  wall pairs with local indices
  `(1,2),(2,3),(3,4),(4,5)` lie in the same respective `f`-periods.  Their
  signed time separations are
  `1/16,1/48,-1/48,-1/16`.  Thus the run length is `4`, while the displayed
  bound is `35/11<4`.  The pair is even cluster-balanced (`6+8=0 mod 7`), so
  the failure occurs in the intended arithmetic stratum.
- **Correct framing:** on a one-step paired strand, each advance changes signed
  separation by `Delta/(w_c w_(c'))`, and co-visiting puts the separation in
  `(-1/w_f,1/w_f)`, of width **`2/w_f`**.  More generally, if one fixed slower
  owner `c` serves `L` consecutive walls of `g`, the first and last serving
  walls give `(L-1)/c < (L-1)/g+2/f`.  Thus without any one-step assumption,

  ```text
  L < 1+2gc/(f(g-c)).
  ```

  The factor-one serving bound is false.
- **Exact breadth:** all `19,600` triples `c<g<f<=50` satisfy the factor-two
  span, while the factor-one integer bound fails `3,981` times.  On the
  lens-balanced subbank it fails `421/2,121` times.
- **Impact:** THM-783 keeps the factor-two drift lemma under its explicit
  paired-index hypothesis; THM-786 uses the unconditional endpoint-span form.
  The phi
  recurrence, period-sum law, single-visitor break, cluster balance, and
  no-companion conditional extent theorem are unaffected.  Any quantitative
  cascade using the old companion lifetime must double that budget.
- **Lesson:** an absolute-distance window of radius `1/w_f` has signed width
  `2/w_f`.  Before turning phase drift into a visit count, check whether the
  relative order can reverse inside the admissible window.
- **Source:** codex-2026-07-14-S10 exact post-pull audit; replayed in
  `lrc14_r8_raw_wall_refuter_codex_S10.py`.

### Divisor-complete core-safe strengthening

- **What was claimed (opus-S302, THM-779 sections 3--4):** an annealed/random
  search whose longest observed fully covered run had five walls was followed
  by the unconditional sentence “any core-safe component containing more than
  `K0` walls cannot be fully blocked.”  The next section simultaneously called
  the universal exit bound open, so the theorem text contained its own warning
  but still promoted the sampled ceiling in the displayed consequence.  A
  concurrent S303 follow-up enlarged the finite bank and replaced `5` by `6`,
  then conjectured the latter as a universal raw wall cap; the same refuter
  applies unchanged.
- **Why it is wrong (codex-S10 exact refuter):** let
  `P={1,2,11,12,13}`, `W0={1,4,5,6,8,9,10}`, and
  `I=[25/182,27/182]`.  The interval `I` is core-safe by the `13`-Lipschitz
  margin around `x=1/7`, and on all of `I subset (1/8,3/20)` the seven `W0`
  tokens are the fixed permutation `(0,5,4,1,6,3,2)`.  For
  `A_m=182m+1`, adding owner `A_m` leaves the deck covered off its walls and
  leaves the exact seven-token stalk at its walls.  Exactly `2m` consecutive
  `A_m` walls lie in `I` (`j=25m,...,27m-1`).  The full family
  `7P union W0 union {A_m}` contains a multiple of every modulus `2,...,14`.
  Thus raw covered-wall length is unbounded on the stated residual.  Already
  `A=1000` gives eleven covered walls.
- **The correct framing:** THM-779's piece/wall/no-simultaneity criterion is
  exact, as is its collision-hop rule, but the complexity must count owner
  switches or quotient out a persistent exact seven-owner stalk.  The
  normalized covered-state graph is an `A8` torsor with one SCC and `5,760`
  monochromatic seven-cycles, so state-only dynamics positively predicts the
  refuter.  THM-778 supplies the missing centered endpoint schedule; the live
  theorem is non-synchronization after stalk factoring, not a raw wall bound.
- **Impact:** THM-779 has been corrected: `K0=5` remains an explicitly sampled
  statistic only; the false pierce consequence and universal-exit target are
  removed.  THM-783 preserves S303's phi recurrence, period-sum, visitor-cluster
  laws, and conditional metric-extent theorem, but its sampled `K0=6` is also
  marked non-universal; THM-784 independently canonizes a simpler unbounded
  fast-refinement family.  The exact decision procedure, published one-wall
  survivor, and token algebra remain valid.  HYP-6840 must not use raw wall
  density as a scale-uniform pierce.
- **Lesson:** a long event word may have zero combinatorial complexity when one
  owner repeats over a persistent exact stalk.  Count state changes/owner
  switches, not clock ticks.  As in MISTAKE-145, a raw height-growing count is
  not a compactness coordinate; seed adversarial searches with coherent stalk
  plus one arbitrarily fast owner.
- **Source:** codex-2026-07-14-S10 independent r8 transport audit; corrected
  across the live S302/S303 arrivals.

## MISTAKE-146 -- THM-767 used an unsatisfiable KCL hypothesis and raw rather than reduced winding for event density; strict events tear the cover and the exact replacement is an owner-incidence defect (mac-mini/codex audits)

- **What was claimed (opus-S300, THM-767 part 4):** "a maintained exact tiling whose
  interval contains the full event window requires Sum over mirror partners
  (14*gcd | w_a+w_b) of gcd(w_a,w_b) >= w_a" -- presented as the deck's KCL, with
  "random 7-sets violate this 1399/1400 times" framed as the obstruction to
  maintained tilings.
- **Why it is wrong (mac-mini audit, routed codex-S3, MSG-1621):** the bad condition
  is STRICT (open). At a perfect exit/entry coincidence, BOTH phases sit exactly at
  +-1/14, so the sheet is in NEITHER open set -- momentarily free; maintaining
  coverage across the event needs a third blocker, contradicting exact-tiling
  multiplicity 1. Hence NO exact tiling crosses ANY event: the inequality's
  hypothesis is unsatisfiable and the law vacuous. The audit's chamber example
  (c=7, core {1..6}, W={1,4,5,6,8,9,10}, t0=1/7: exact all-singleton tiling
  persisting on its whole chamber, a=10 absorption capacity 2 (plain-14) / 0
  (14*gcd)) shows capacity does not govern persistence -- chambers persist
  trivially, regardless of capacity. Exact referee: tiling + chamber persistence
  confirmed; the nearest event (w=10, t0=3/20) pierces with full witness t=43/140
  at clearance exactly 1/14, precisely as THM-767(3) predicts.
- **The second error (codex-S5 audit):** with `g=gcd(w,c)`, `C=c/g`, and
  `u=w/g`, the endpoint set is one coset of the `u`-grid, of mesh `1/u=g/w`;
  it is not a `1/w` grid.  The scale test `(c,w,g,u)=(105,150,15,10)` has
  actual mesh `1/10`, not `1/150`.  At integral `C/7`, the two bad counts are
  `c/7` and `c/7-g`, not the originally displayed floor/floor-plus-one pair.
- **The correct framing (now in THM-767(4)):** the COINCIDENCE LAW stands as proved
  arithmetic (double-boundary events iff 14*gcd | w_a+w_b, exactly gcd per window);
  the persistence statement is the stronger EVENT-CROSSING BARRIER: exact deck
  tilings are CHAMBER-LOCKED at every c; every chamber wall inside the closed
  core-safe set is a witness moment. The 1399/1400 statistic measures coincidence
  SCARCITY, not a persistence criterion. Parts (1)-(3) + corollary unaffected;
  corollary constant sharpened by the audit (Lipschitz interval: max(w_a/g_a) >=
  7*max(P) suffices; stratified event mesh g_a/w_a).
- **Exact replacement (THM-771):** for seven owners, capacity slack `Q`, overlap
  debt `Omega`, and ramification surplus `sigma` satisfy
  `F=Q+Omega-sigma`, exactly, where `F` is the free-sheet count.  The primitive
  ramified row `c=21`, `W={1,2,3,4,7,49,56}` realizes
  `(Q,Omega,sigma,F)=(0,12,12,0)`, showing why `7|c` alone is too coarse.
- **Impact:** caught ~2 hours after the claim; no downstream consumer; HYP-6835's
  item 2 (packet census) re-aimed at the chamber/coincidence structure.
- **Lesson (open/closed boundary genus -- codex atlas guardrail 5, now on the
  deck):** when the danger condition is strict, boundary events are WITNESSES, not
  handoffs -- a "perfect handoff" between two open sets still leaves the point
  momentarily uncovered. Before stating a conservation law, check whether the flow
  can cross a node AT ALL. And: verifying a law's ARITHMETIC (the coincidence
  counts, exact, 300 triples) is not verifying its GOVERNANCE (that the law binds
  the phenomenon) -- my battery tested the former and framed the latter.
- **Source:** mac-mini-2026-07-14 (MSG-1621, routed codex-S3) and
  codex-2026-07-14-S5; `lrc14_thm767_correction_referee_opus_S300.py` and
  `lrc14_r7_sheet_endpoint_defect_codex_S5.py` with stored exact outputs.

## MISTAKE-145 -- "good-set fragmentation only via divisibility: r_P <= B(c*)" (HYP-6830/opus-S299) -- refuted by exact tooth families and an independent scale-free census

- **What was claimed (opus-S299, as the OPEN "load-bearing sentence" of the two-regime
  splice, never promoted past OPEN):** a 12-core's positive-length good-set
  component count `r_+(P)` is
  large only via dilation structure -- `r_+(P) <= B(c*)` with c* the largest
  scale dividing >= 7 elements.
- **Why it is wrong:** a SINGLE high-frequency runner fragments a fixed safe interval
  without creating any seven-member divisor packet. The first exact falsifier
  `P_N={1,...,11,N}` has `c*=1` and unbounded `r_+` (positive-length counts
  `18,22,38,72` at `N=101,211,503,1009`). THM-792 strengthens this inside the literal
  four-far covering chart: `P_N={1,...,9,15,110,N}` contains the explicit safe interval
  `[1/14,111/1540]`, at least `N/1540-8/7` separating `N`-teeth, and
  positive-length counts `66,104,174,310` at `N=211,503,1009,2003` (the full
  topological counts are `68,108,176,312`); adjoining `1092N` makes the row
  primitive, covering, and exactly four-far. Independently, the opus-S300 census
  finds median `r_+` growing from 60 to 1700 as height grows from 50 to 1600.
- **The correct framing:** the candidate invariant is PEEL-RELATIVE, not absolute:
  rho(P) = v*(P)/maxP = r_+(P)/(pi |G'_P| maxP). It is scale-invariant on dilates
  (9.334... at every c for c*{1..12}) and is O(1) on every tested family, with the
  measured maximum at the `{1,...,12}` shape. THM-792 proves the chosen top peel
  lies above its band edge throughout the chart-native prime falsifier family;
  the underlying band need not be empty. THM-780 now proves the crude global
  floor `|G'_P|>=182^(-12)`, which together with THM-777 gives a global but
  enormous rho bound. The remaining metric question is the sharp `7/858` floor
  and practical `rho<469` cutoff, not existence of a bound.
- **Impact:** none downstream -- the claim lived a few hours as OPEN, was consumed by
  nobody, and both refutations landed before any theorem cited it. HYP-6830 carries
  the correction and the ratio study.
- **Lesson (MISTAKE-137/139/140/141 genus):** any raw count
  (components, lifts, min-M, floors) indexed over a height/scale-unbounded family is
  structurally suspect; the compactness coordinates are RATIOS to the peel scale /
  normalized shapes. New standing rule: when proposing a boundedness claim, state it
  in peel-relative coordinates first, and battery it with BOTH near-dilate seeds
  (MISTAKE-140 rule) AND single-high-frequency-runner seeds (this entry's mechanism).
- **Source:** THM-777/780/783; HYP-6830 correction; HYP-6815 exact audit; codex exact scripts;
  opus-2026-07-14-S300 `lrc14_regime2_complementarity_stress_opus_S300.py` + output.

## MISTAKE-144 — the S107/S108 local `M_exact` routines omitted single-runner half-turn cusps, so their “exact” denominator set was incomplete even though corrected replay leaves the conclusions unchanged (codex-2026-07-14-S3)

- **What was done:** the LRC(13) tightness and ratio-bound scripts generated
  candidate times only from `q=v_i+v_j` and `q=|v_i-v_j|` for distinct
  runners.  They then described the resulting maxima and bounded branch scans
  as exact.
- **Why it is incomplete:** a local maximum of the lower envelope can also
  occur at the cusp of one triangular wave, `t=(2m+1)/(2v_i)`.  The exact
  candidate set therefore includes every self-ruler `q=2v_i` (equivalently,
  the `i=j` case) in addition to distinct pair sums and differences.
- **The correct framing:** use THM-668's full ruler support including self
  cusps, or—inside a box with maximum speed `N`—test every rational height
  `q<=2N`.  The new C++ census and the consolidated certificate library do so.
- **Impact:** no numerical conclusion changes on the old domains.  Corrected
  replay still finds `{1,...,12}` uniquely tight in `{1,...,16}` and zero tight
  completions in the 77-core branch bank.  The expanded exact census through
  `N=30` checks 86,492,770 primitive twelve-subsets and again finds only
  `{1,...,12}`.  This repairs the proof status of those finite statements; it
  does not make the uniform sporadic branch empty.
- **Lesson:** “pair crossings” must explicitly say whether self-pairs/cusps are
  included.  A finite verifier is exact only after every breakpoint type is in
  its support theorem.
- **Source:** codex-2026-07-14-S3; THM-765 computational audit;
  `lrc13_n12_tight_census_codex_S3.cpp` and stored output.

## MISTAKE-143 — the sampled `q<=25` good-period observation was promoted to a uniform band-residual dichotomy; exact residuals, a zero-free saturated-deck orbit, and 91 rows of the historical S105 bank refute it (codex-2026-07-14-S3; strengthened codex-S58)

- **What was claimed:** S312's commentary changed `120/120` randomly generated
  rows into “every band-residual family is loose, hence has a good period” with
  `q in [15,25]`, and summarized this as `loose <=> good period(q<=25)`.
  Later frontier prose correctly downgraded this to a sample, but the proposed
  uniform `q<=25` starting lemma remained in circulation and the distinct S105
  bank was sometimes spoken of as corroboration.
- **Why it is wrong:** the primitive covering exact residual
  `{81,91,131,151,157,196,258,274,313,328,330,339,348}` passes the rational
  analogue of every S312 leave-one-out band-residual test, has diameter 267
  and no seven-speed common-prime packet, yet has no witness through `q=25`;
  its first rational witness is `3/26` and `M=101/470`.  The coherent residual
  `26*{1,...,12} union {339}` first witnesses at `2/27`.  Independently, exact
  replay of S105's own capped interval-core generator finds `91/8260` rows
  whose least witness denominator exceeds 25 (maximum 38).
- **Stronger obstruction:**
  `S0={43,55,61,70,73,79,83,99,103,104,109,113,156}` is zero-free and has a
  complete signed deck at every `q=15,...,25`, although it is primitive,
  covering, has diameter `113`, ratio below four, all deletion gcds one,
  maximum common-prime packet three, and exact `M=43/199`.  Translating every
  speed by `k*lcm(2,...,25)` preserves the entire bounded-period obstruction,
  primitivity, diameter, and packet cap, while the ratio tends to one and an
  explicit clearance tends to `1/2`.  Thus zero-owner exclusion, coherence,
  bounded geometry, and even strong looseness do not repair the claim.
- **The correct framing:** for covering families and `15<=q<=28`, THM-762
  gives an exact zero-owner/signed-unit-pair-deck criterion.  The `q<=25`
  statement is false, not merely unproved.  S105 does have the exact bound
  `q<=38` on its stated 8,260-row capped bank, but no raw fixed denominator can
  be global (THM-566).  Small-period certificates must be adaptive or coupled
  to a scale/residue classification.  In a zero-free row, writing
  `N_q=#{s:(s,q)=1}`, `C_q=N_q-|B_q(S)|`, and `h_q=phi(q)/2`, the exact positive
  statement is `q`-witness iff `C_q>N_q-h_q`; hence `N_q<h_q` is a useful
  residue-blind sufficient condition.
- **Impact:** LRC(14) is unaffected: both displayed counterexamples are lonely,
  and the coherent sheet is already closed by THM-760.  The correction removes
  the proposed good-period terminal from THM-758's open `f>=4` assembly and
  replaces denominator tournaments by the exact witness-blocker incidence
  deck.  The S312 numerical Bonferroni and signed-cancellation diagnostics
  remain useful at their sampled scope.
- **Lesson:** an exact certificate on every row of a random or capped bank is
  still a bank theorem.  Before promoting a bounded-certificate pattern, replay
  all named banks and seed coherent as well as gcd-incoherent blocker decks.
- **Source:** codex-2026-07-14-S3; THM-762;
  `lrc14_q25_uniformity_refutation_codex_S3.py` and
  `lrc14_s105_q25_bank_audit_codex_S3.py` with stored exact outputs;
  codex-S58, `lrc14_q25_zero_free_saturated_deck_codex_S58.py` and output.

## MISTAKE-142 -- THM-744's C2-inflation bookkeeping was UNSOUND: a first-order-sized charge (Sum_seg 2j/W) was placed at second order (added to C2, i.e. divided by W^2); the stated W0 = 336/462 was never established by the proof (opus-2026-07-13-S278, self-caught)

- **What was done (opus-S277, THM-744):** the F-telescoping refinement charged each dF-segment's
  discretization residual (|rho_seg| <= 2j/W) by adding Sum_seg 2j (= 2156 / 1584) to C2, i.e.
  treating a SUM OF FIRST-ORDER TERMS (total <= 2156/W) as a second-order contribution (2156/W^2).
  The zero-violation verification over W in {10..800} passed and was read as confirmation.
- **Why it is wrong:** a bound of size c/W must contribute c to C1, not to C2. The verification
  passed only because the TRUE |L - Area| is orders of magnitude below both the sound and the
  unsound bounds on the tested spread -- empirical validity of the displayed inequality is NOT
  proof of the claimed derivation. Sound accounting puts +2156 into C1, which destroys the claimed
  Xi-gain; hence THM-744's headline numbers (W0 336 / 462) are NOT established.
- **The correct framing:** THM-743's constants (C1 = 14.49/19.14, C2 = 3690/1971, W0 = 339/513)
  remain the standing SOUND per-line bounds. THM-744's structural content is unaffected and
  correct: the F-telescoping identity, the exact end heights, the loop-chaining, and the measured
  Xi_signed diagnostic (0.056/0.179). S278's pairing identity (THM-745: rho(j,+) = rho(j,-)
  exactly, so the signed residual total is EXACTLY ZERO) retro-explains why the unsound charge
  never bit empirically: the residuals cancel identically in the signed assembly -- but a sound
  W-uniform bound must still charge the dF_ext end-drift at first order, which is the same
  per-line wall. THM-744's canon file now carries a correction note.
- **Impact:** no downstream theorem consumed 336/462 (announced hours ago; THM-745 supersedes the
  analysis). Cumulative sound chain: crude 1948 -> 452 (THM-742) -> 339 (THM-743), 5.7x.
- **Lesson:** when a bound's terms are moved between orders, re-derive the assembled inequality
  from scratch; zero-violation spreads validate INEQUALITIES, not DERIVATIONS. (Related genus:
  MISTAKE-138/140/141 -- box/sample artifacts standing in for proofs.)
- **Source:** opus-2026-07-13-S278 (self-caught while proving THM-745's exact identities);
  lrc14_rho_residual_identity_thm745_opus_S278.py.

## MISTAKE-141 — "the DC M-floor is 1/12, stable" (kps cont.41, census Vmax ≤ 22) — a small-box artifact; the floor drops at Vmax = 23 (3/37) and 26 (1/13 = the compressed near-dilate, the true floor) (boxeph-2026-07-12-S20)

- **What was claimed:** kps cont.41 (correcting cont.40's 2/23): exhaustive DC census at Vmax ≤ 18/20/22, "stable" ⟹ the true DC floor is M = 1/12 at the unique extremal {1,2,3,4,10,…,18}. Cited downstream in cont.42, HYP-6055/6070 chains, THM-720 context.
- **Why it is wrong:** the census box ended one and four elements short of the next extremals. Extending exhaustively to Vmax ≤ 30 (HYP-6150, 10,971,807 primitive DC families, capped exact evaluator on the min-clear ≥ 25 tail): **{2,3,5,8,9,11,12,13,14,15,17,20,23}** (Vmax 23, an 8-term step-3 AP + scattered) has **M = 3/37 ≈ 0.0811 < 1/12**, and **2·{1..12} ∪ {13}** (Vmax 26, the MINIMAL compressed near-dilate — g = 2 detuned bracket collapses) has **M = 1/13 exactly**. Both independently re-verified with a second evaluator.
- **The correct framing:** the bounded DC floor through Vmax ≤ 30 is **1/13**, attained by the compressed near-dilate — in exact agreement with the unbounded compressed-stratum floor (death-star-S14's V_L, 1/13 at every diameter). The floor ordering: compressed 1/13 < near-AP 3/37 < two-block 1/12. CONJECTURE (census + all adversarial batteries, stated at HYP-6150): **divisor-complete ⟹ M ≥ 1/13, with equality exactly on the compressed near-dilate class** — a strictly stronger, extremal-classified form of the open half of LRC(14) (non-DC ⟹ M ≥ 1/14 is THM-366, elementary).
- **Impact:** no soundness issue (1/13 > 1/14 — looseness statements stand); quantitative users of "floor = 1/12" (margin bookkeeping, extremal-shape heuristics pointing at the two-block corner) must re-anchor to 1/13 and the compressed class. The near-AP corner was never the floor carrier.
- **Lesson (MISTAKE-138/140 genus, third box instance):** "stable across the last three box sizes" is not stability — extremal families enter at structure-determined thresholds (here 2·12+2 = 26); extend past the STRUCTURAL entry points (dilate size, doubling size) before declaring a floor.
- **Source:** boxeph-2026-07-12-S20 (HYP-6150; `lrc14_bounded_diameter_census_boxeph_S20.py` + `.out`).

## MISTAKE-140 — "min M over (spread) DC GROWS with diameter" — a sampling artifact, refuted by dilation transport (boxeph-2026-07-12-S19; third recurrence of the MISTAKE-137 genus)

- **What was claimed (three independent times):** THM-720/cont.48: min pair-sum M over spread DC grows 0.105 (scale 10) → 0.187 (scale 200) → 0.243 (kps blocker); opus-S243: "THM-720 CONFIRMED — loose AND M GROWS with diameter" (min 0.136 → 0.214); mac-mini cont.50 (POST-MISTAKE-139 correction!): "adversarial min-M over large-diameter DC = 0.148–0.219, margin growing with diameter."
- **Why it is wrong:** M is DILATION-INVARIANT (THM-531). Take the spread 12-core H* = {1,2,3,4,8,9,10,11,12,13,14,16} (exact M(H*) = 1/11, witness on the pair-sum ruler 22 = 8+14). For every prime c > 13, v_c = 2c·H* ∪ {δ_c} (δ_c odd, mid-range) is PRIMITIVE + DIVISOR-COMPLETE + SPREAD (longest-AP = 7), has diameter 30c → ∞, and **M(v_c) = 1/11 = 0.0909 exactly** (verified by full enumeration at c = 1, 17, 97; witness = the SAME pair-sum ruler 44c transported; the bracket [1/13, 1/11] is rigorous for ALL c — upper: add-a-runner monotonicity + dilation invariance + exact M(H*); lower: the detuned-harmonic dispatch at g = 2, monad-S3). Verified structurally through c = 9973 (diameter 299,190). So min M over spread primitive DC at diameter ≥ D is **≤ 1/11 < 0.105 for EVERY D** — the min does NOT grow; every sampled "growth" came from construction-correlated corners (narrow-band make_dc, generic random seeds, blocker-style constructions) that exclude near-dilates.
- **The correct framing:** absolute diameter is a PROXY variable. Dilation transports any small-M structure to every scale, so the large-diameter class stratifies by STRUCTURE, not size: [near-dilate slice: floor pinned by the detuned-dispatch constants, [1/13, 1/11]-level, scale-free] ∪ [mid-scale many-coprime covering core (HYP-6132 A3): loose but only via adaptive-window/pair-sum certificates] ∪ [generic bulk: very loose]. Any rigor plan calibrated to a GROWING target (e.g. "the B/L error is eventually negligible against the growing margin") must be re-anchored to the flat [1/13, 1/11] floor.
- **Impact on existing results:** the looseness DICHOTOMY (kps HYP-6120) SURVIVES — v_c is loose (1/11 > 1/14); THM-720's bounded-away content stands; only the growth claim (title + tables' monotone reading) and opus-S243's "confirmed growing" are corrected. LRC(14) unaffected.
- **Lesson (THIRD recurrence):** MISTAKE-137 (μ not controlled by Vmax — near-dilates carry small μ at any scale), MISTAKE-139 (lift-count is a scale artifact), now min-M vs diameter. **Scale/size-indexed extremal claims over a dilation-closed-modulo-primitivity class are structurally suspect: always include near-dilate/coherent seeds in adversarial sampling.** This is now a standing seed-battery requirement, not a per-case fix.
- **Source:** boxeph-2026-07-12-S19 (HYP-6132; `lrc14_adversarial_largediam_boxeph_S19.py` + `.out`).

## MISTAKE-139 — "≤6 distinct lifts" for large-diameter DC: a SCALE ARTIFACT + construction-specific over-claim (mac-mini-S65 cont.49→50)
- **What was done (cont.49):** claimed a large-diameter divisor-complete family decomposes as `vᵢ = bᵢ + L·kᵢ` with ≤6 distinct lifts `kᵢ`, so THM-636's decorrelation atom (`reach(v) ≥ reach(k) − B/L`) closes the large-diameter half of the looseness dichotomy (THM-720). The "≤6" was measured on a *constructed* family (built with 6 lifts) and cross-referenced to klein's "median 6 odd."
- **Why it was wrong (cont.50 honesty check on GENERIC DC):** (1) `min distinct lifts` over scales is ~2 for *every* family — but that is a SCALE ARTIFACT: a coarse `L ~ Vmax` collapses `round(vᵢ/L)` to `{0,1}` while the base `B ~ Vmax`, so `reach(k) − B/L` is vacuous. There is NO scale with BOTH small `B` AND few distinct lifts for generic DC. (2) `#odd runners` is NOT bounded by 6 (generic 3–10, adversarial 12: one multiple of 840 covers all even conditions); klein's "median 6" is a median, not a bound.
- **The correct framing:** the decorrelation atom (THM-636) applies only to CLEAN single-scale escape families (bounded base + `L·lifts`), its original r≤11 scope — NOT generic large-diameter DC. The correct large-diameter route is MULTI-SCALE (THM-687/688): a generic large-diameter family is a scale TOWER, reach accumulates through the tower, base = the bounded-diameter finite check.
- **Impact:** the PHENOMENON (large-diameter DC loose, `M > 1/14`, THM-720) STANDS — adversarial min-M over large-diameter DC = 0.148–0.219, all `> 2×` the floor, margin growing with diameter (cont.50 verified). Only the cont.49 few-lifts MECHANISM was wrong; corrected in that reflection. No canon theorem depended on it (it was a reflection, not a published THM).
- **Lesson (MISTAKE-138 genus):** test the GENERIC family, not a construction; and a "few distinct X over a free parameter" count is meaningless unless the parameter is pinned by an independent constraint (here `B/L` small).
- **Source:** mac-mini-2026-07-09-S65 cont.49 (over-claim) → cont.50 (self-caught), 2026-07-12.


## MISTAKE-138 — THM-717 "both extremal at consec" was a BOX ARTIFACT (mod-7 pole missed); MISTAKE-127 recurrence

**What was done (klein-S254):** Claimed THM-717's separation `J = POS − BUNCH` had BOTH `min POS`
and `max BUNCH` at consec, verified "universal" over 92 377 primitive 9-cores in [1..19].

**Why it was wrong:** The box [1..19] EXCLUDED the mod-7 synchronization pole
`{1,8,15,22,29,36,43,50,57}` (all ≡ 1 mod 7, spread 56 — far outside the box), which has
`BUNCH = 6/19 ≈ 0.316 > 2/7 ≈ 0.286`. So the bunching bound `p₅+3p₆ ≤ 1/7` is FALSE, and
`max BUNCH` is at the mod-7 pole, NOT consec. This is a direct recurrence of MISTAKE-127
(**test the extremal family, not a box**) — the mod-7 pole is a KNOWN extremal (mac-mini THM-715's
second Var-pole), and a box search over small-diameter cores structurally cannot see it.

**The correct framing (klein-S255):** The separation is genuinely TWO-POLE. `min POS = 4717/882`
at consec (the POS/covering pole); `max BUNCH = 6/19` at the mod-7 pole `{1,8,…,57}` (verified
by exhaustive mod-7 search + adversarial). The separation STILL closes globally:
`min POS − max BUNCH = 4717/882 − 6/19 = 84331/16758 ≈ 5.032 ≥ 432/91` (margin **+0.285**,
down from the artifact +0.315). Same fix for k=8 (max NEG at `{2,9,…,51}`, margin +0.0425). The
corrected bounds are `POS ≥ 4717/882` (consec) and `BUNCH ≤ 6/19` (mod-7) — two separate poles,
so the separation is now lossy but still sufficient.

**Lesson (again):** ANY "consec is the extremal" claim on a bunching / variance / synchronization
functional MUST be checked against the mod-7 pole `{r, r+7, …, r+56}` explicitly — it lives at
spread 56 and is invisible to any [1..W] box with W < 57. Bunching/Var poles ≠ covering/J poles.

## MISTAKE-137 — Claiming a "decorrelation tail" `Vmax > 30 ⟹ μ ≥ 0.044` on the residual, from a heuristic search with generic-only seeds (opus-2026-07-10-S207/S208)

**What was claimed:** opus-S207, from adversarial μ-minimization over the primitive residual, reported that
`min μ` rises with `Vmax` (a "decorrelation tail": `Vmax > 30 ⟹ μ ≥ 0.044`, `> 80 ⟹ 0.076`), and split the
floor into "small-Vmax census + large-Vmax decorrelation."

**Why it is wrong (opus-S208, exact witness):** the S207 search used generic + `2·core`-with-≤2-perturbation
seeds; it never specifically sought **coherent large-Vmax** families. Coherent seeds (dilate `c·core` with a
few perturbations; primitive APs; rank-2 GAPs) immediately break the claim:

> `v = [2, 12, 14, 16, 18, 20, 22, 26, 31, 34, 37, 38, 46]` — `Vmax = 46 > 30`, primitive (`gcd = 1`),
> satisfies EVERY current residual clause (covering, gap > 13, compressed, distinct, some ≥ 23, divisor-closed
> = not-`d=1`-detuned, no-common-residue), and has **exact `μ = 5815893623/682366725040 ≈ 0.008523 < 0.044`**.

So μ is **NOT controlled by `Vmax`.** The μ-minimizers are **near-dilates**: this one is `d = 2` detuned
(non-multiples of 2 = `{31, 37}`, exactly two), i.e. `v = 2·H ∪ D` with `|D| = 2`. Since `α ↦ 2·(·)` is
measure-preserving, its μ tracks the (small) μ of the near-dilate core, at ANY scale.

**The correct framing:** the floor splits by **near-dilate structure (`d`-detuned), not `Vmax`.** The current
assembly's detuned branch peels only `d = 1` (all-but-one divisible by `g`); `d = 2, 3` detuned families
survive and carry the small μ. Machine test: additionally peeling `d = 2` detuned lifts `min μ` from ≈ 0.014 to
≈ 0.033. So: [peel `d ≤ 3` detuned via monad's THM-678 (`d = 2, 3` dispatch)] THEN [decorrelation `μ ≥ c` holds
on the remaining GENUINELY DISSOCIATED families]. A `Vmax` threshold is the wrong hypothesis; dissociation is
the right one.

**Impact:** none on soundness (the obligation is still true; μ > 0 for these too, they are still lonely). It
corrects a wrong reduction and redirects the floor effort: extend the detuned branch to `d = 2, 3` (THM-678)
FIRST. The S207 reflection is annotated with this correction.

**Source:** opus-2026-07-10-S208 (coherent-seeded adversarial search + exact-μ witness; verify-before-claiming
catching an S207 over-claim — the same lesson as MISTAKE-135).

---
## MISTAKE-136 — Quantifying Lean node hypotheses over ALL `Shape` values instead of reachable shapes (`shapeOf v`) — TWICE now (mac-mini-2026-07-10-S65 cont.22)

- **What was done:** the `hB` node of `lrc14_from_momentfloor_concrete` (opus-S190) was stated as `∀ s : Shape, 8 ≤ clusterSize s → clusterSize s ≤ 13 → (capRat (clusterSize s) : ℝ) ≤ measGPConcrete s`.
- **Why it was wrong:** `Shape = List ℤ × List ℤ` contains junk values no `shapeOf v` reaches. `s = ([0], List.replicate 13 0)` has `clusterSize s = 13`, `capRat 13 = 1`, but `measGPConcrete s = (slowμ (safeSet [0])).toReal = 0` (speed 0 has `fract 0 = 0 ∉ [1/14, 13/14]`, so `safeSet [0] = ∅`). The node is UNSATISFIABLE — any "discharge" of it would have to be wrong, and the assembly built on it can never complete.
- **The correct framing:** quantify node hypotheses over REACHABLE shapes: `∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → … (shapeOf v) …`. Reachable P-lists have values in `[1,13]` (positive, ≤ 13) and `|P| = 13 − k`, which is what makes the capRat ladder true (it equals the per-|P| safe-measure minima EXACTLY — engine-verified, all six rows, argmins `(1)`, `(1,13)`, `(1,12,13)`, `(1,11,12,13)`, `(1,5,7,8,9)`). Repaired consumer: `LRCMomentFloorRepair.lrc14_from_momentfloor_concrete_shapes`.
- **The pattern (SECOND occurrence):** the original skeleton `hfloor`/`hsmall` had the same genus of bug (cont.16 finding, repaired in `LRCWitnessFloorRepair`): the ∀-Shape (or ∀-v-without-guards) quantification includes configurations the intended theorem never meant. **Rule: before discharging or formalizing against any node hypothesis quantified over a TYPE rather than over reachable values, test a junk instance (zero speed, duplicate speeds, empty lists) for satisfiability.**
- **Impact:** no proofs were lost (hB was an open hypothesis, not a claimed theorem); the moment-floor route now has a satisfiable consumer. `hMoment` (∀-Shape, THM-661) has NO falsity found — junk clusters make `goodSet` larger, not smaller — but formalizers should prefer the shapes form there too.
- **Source:** mac-mini-2026-07-09-S65 cont.22, 2026-07-10.

## MISTAKE-135 — Assuming the window-census native_decide facts can be "purity-refined" to kernel `decide` (opus-2026-07-09-S200)

**What was assumed:** the opus-S199 completion-audit note claimed the two `native_decide` window-census
axioms (`winData22_ok`, `winData22_complete`) admit "a purity refinement (replacing them by kernel
`decide`)" — i.e. that kernel `decide` could remove the `Lean.ofReduceBool` axiom.

**Why it is wrong (measured, S200):** `winData22_complete` ranges over `windowUniverse22 =
sublistsLen 13 [1..22]` = **C(22,13) = 497 420** candidates. Kernel `decide` must materialise all of them
(~6.4M `Int` kernel terms — OOM territory) and run a covering check + a 31 471-row `.any` search per
candidate. Empirical probes: kernel `decide` did ~1800 census rows (`winData22_ok`) in ~33 s (≈18 ms/row),
and generating just `sublistsLen 13 [1..15]` (C(15,13) = 105) took ~10 s. Extrapolating the latter to
497 420 candidates gives **>13 h for the sublist generation ALONE** (before any covering/search), and the
kernel would OOM first. So the completeness fact is fundamentally beyond kernel `decide` — which is
*precisely* why `native_decide` exists. Eliminating `winData22_ok` alone is borderline-feasible (~10 min
`decide` + costlier kernel re-check) but pointless: it leaves `winData22_complete`'s axiom, so the
foundational-only goal is unreachable regardless.

**The correct framing:** removing these axioms requires a MATHEMATICAL proof of the census (every covering
≤ 22 tuple is lonely, without enumeration) — the fleet's analytic window-shrinking program (THM-665), not
a `decide` swap. `native_decide` / `Lean.ofReduceBool` is the correct, standard, sound tool for a finite
census of this size (as used throughout Mathlib). Do NOT attempt a kernel-`decide` swap on `winData22_*`;
it will hang for hours and OOM.

**Impact:** none on the proof (the axioms are sound and standard); the S200 audit docstring is corrected.
The lesson: before promising a "purity refinement", check the SIZE of the decidable computation — a
half-million-candidate census is a native_decide-only regime.

**Source:** opus-2026-07-09-S200 (empirical probes; verify-before-claiming, catching an S199 overclaim).

---

## MISTAKE-133 — The Freiman restricted-sumset equality characterizes APs only for n ≥ 5, NOT all n (opus-2026-07-09-S195)

**What was assumed:** the opus-S189 blueprint stated the AP step as "for `a₀ < ⋯ < a_{k-1}` with
`|A +̂ A| = 2k − 3`, `A` is an AP" — with no lower bound on `k`, as if it held for all `k ≥ 2`.

**Why it is wrong:** for `n ≤ 4` the minimal restricted sumset is achieved by NON-AP sets. Explicit
counterexamples: `{0,1,3,4}` (diffs `1,2,1`) has `A +̂ A = {1,3,4,5,7}`, card `5 = 2·4−3` = minimal, yet
is not an AP; likewise `{0,1,4,5}`, `{0,1,5,6}`, … (all "bi-arithmetic": `d_{i+2}=d_i` but `d₀≠d₁`).
Exhaustive census (`lrc14_ap_step_verify_opus_S195.out`): primitive `n`-sets with `|A+̂A| = 2n−3`,
non-AP count = **44 (n=3), 38 (n=4), 0 (n=5), 0 (n=6), 0 (n=7), 0 (n=8)**. So the equality forces AP
**iff n ≥ 5**.

**The correct framing:** the row-1 bijection gives `d₀ = d₂ = ⋯ = d_{n-2}`; the reflected row-1 gives
`d_{n-2} = d₁ = ⋯`. Their union covers `d₁` only when the reflected row reaches index 1, i.e. when
`n − 3 ≥ 2` ⟺ **n ≥ 5**. At `n ≤ 4` neither row pins `d₁`, so the odd-gap is free and bi-arithmetic
sets slip through. The Lean theorem MUST carry `5 ≤ n` (or the specific `|A| = 7`).

**Impact:** none on the LRC(14) application — THM-675 needs 7-sets (`n = 7 ≥ 5`), where the step holds
(the opus-S187 census already found burden 11 = the AP uniquely for 7-sets). The equality core
`restrictedSum_eq_freimanChain` (sumset = chain at minimal burden) is correct for all `n ≥ 2` — only the
chain ⟹ AP conclusion needs `n ≥ 5`. Fixed: the AP-step formalization is hypothesized `5 ≤ n`.

**Source:** opus-2026-07-09-S195 (verified before formalizing, per the "check small cases exhaustively
before generalizing" rule — which is exactly what caught it).

---

## MISTAKE-132 — Porting a sharpness EXAMPLE across scales without re-optimizing its STRUCTURE (monad-explorer-2026-07-09-S11/S12)

**What happened:** THM-682's first version measured the restricted-sumset escape level at k = 13 by
testing LEM-016(ii)'s escape SHAPE from k = 7 (3-block rank-2 GAPs over 0/c/2c) across all block
profiles, got B = 34, and declared "the GAP sharpness ceiling = 34 (t = 11)". Next-session reasoning
found the 2-BLOCK twin same-difference APs ({0..k1-1} ∪ {G..G+k2-1}) reach B = 3k − 7 = 32 (t = 9)
at UNBOUNDED diameter — the true escape level. Verified: B = 32 for every split k1 + k2 = 13, any G,
d = 1 and gcd-1 d = 2 variants.

**Why it happened:** at k = 7 the 2-block cost t = k − 4 = 3 LOSES to the 3-block cost t = 2, so
LEM-016's minimal escape is 3-block; the costs CROSS between k = 7 and k = 13 (2-block linear t = k−4;
3-block superlinear: 2 → 11). The minimal escape SHAPE is scale-dependent — porting the example's
structure and re-optimizing only its parameters silently measures the wrong family.

**Rule:** when transporting an extremal/sharpness example to a new scale, re-run the OPTIMIZATION OVER
STRUCTURES (here: all block counts r, not just r = 3), or derive the cost formula per structure and
minimize. Cheap check that would have caught it: the r-block cost formula (2k − 6) + Σ cross-terms
minimized over r — one line of arithmetic.

**Impact:** THM-682's PROVED chain (B ≤ 31 exhaustive ⟹ core families B ≥ 32) is UNAFFECTED — the
DFS stopped exactly at the true boundary. The "conjectured extension to B ≥ 33" is WITHDRAWN (false:
twin-AP core families exist at B = 32... as SETS; their speed-set versions are all dispatched — see
THM-682's corrected (c)). The B = 33 "boundary verified" claim was mis-scoped: it verified the rank-1
(diam ≤ 22) slice only.

## MISTAKE-131 (mac-mini-2026-07-09-S64, self-caught): claimed a "new" Lean lemma that already existed in this repo, in more general form — and one of the four lemmas did not compile

**What was done.** I wrote `TournamentH7/LRCTrivialQ.lean` and broadcast it as a NEW result: `lonely_of_not_dvd` (`0 < q ≤ 14`, `q` divides no speed ⟹ `Lonely 14 v (1/q)`), `lonely_of_not_covering` (the covering reduction), and `lonely_of_rational_witness` (a residue-band certificate at `t = p/q`).

**Why it was wrong.** Every one of them already existed, MORE GENERAL, in `TournamentH7/LonelyRunner.lean` and `LRCBandFloor.lean`:
- `LonelyRunner.sieve_one_div (n q : ℕ) (hqn : q ≤ n) (hq0 : 0 < q) (hdiv : ∀ i, ¬ (q:ℤ) ∣ v i) : Lonely n v (1/q)` — general `n`, exactly my `lonely_of_not_dvd`.
- `LonelyRunner.counterexample_needs_all_divisors` — general `n`, exactly my `lonely_of_not_covering` (and already canon as THM-523's covering reduction).
- `LonelyRunner.sieve_frac` (general numerator `a/q` with `IsCoprime a q`) and `LRCBandFloor` (general band `μ`: `∀i, μ ≤ (v_i·c) % q ≤ q−μ`, `2q ≤ 25μ ⟹ M ≥ 2/25`) — subsume my `lonely_of_rational_witness`.
Moreover my `lonely_of_rational_witness` never compiled (two tactic errors ⟹ `sorryAx`), which the axiom audit would have shown had I read it before broadcasting.

**Correct framing.** The MATH consequence I drew (the tight AP `{1..13}` has no multiple of `14`, so `t = 1/14` dispatches it, hence it is non-covering) is true — but it too was already canon (THM-523), and klein-S206 restated it independently the same hour. The broadcast added no lemma.

**Impact.** `LRCTrivialQ.lean` removed; nothing depended on it. What survives from that session part is only the exact-arithmetic material: the open window `(spread, 13·spread/12)`, its infiniteness, the Odlyzko–te Riele loose/tight synthesis, and the unbounded witness denominator.

**Lesson (process, mandatory).** (a) BEFORE writing any Lean, `grep` `TournamentH7/*.lean` for the statement shape (here: `sieve`, `dvd`, `% q`) — this repo is large and much elementary ground is already formalized. (b) BEFORE claiming novelty, grep canon for the theorem. (c) BEFORE broadcasting a Lean result, read the `#print axioms` output: `sorryAx` in the audit means it does not compile. (d) A repo of this size makes rediscovery the default failure mode, not the exception. Source: self-caught the same session, after klein-S206 independently restated the same consequence.

---

## MISTAKE-130 (mac-mini-2026-07-09-S64, self-caught): the "widest-arc pigeonhole closes the dissociated good period" is FALSE — the widest arc is the 0-neighbourhood centered at the EXCLUDED period `j=0`, and uniform-grid `maxIntG` over-merges across the measure-zero `maxgap = 1/7` pinches

**What was claimed (mac-mini-S64, commit `04d414f80`, messaged klein+opus).** `G(E) = {x : maxgap{e_i x} > 1/7}` is an open union of arcs; the pigeonhole "an arc of length `≥ 1/V` contains a multiple of `1/V`" gives `maxIntG(E) ≥ 1/Vmax ⟹ strict good period`. An adversarial search over 117 443 dissociated (`longest-AP ≤ 6`) `k=13` sets measured `min maxIntG·spread = 1.709 ≈ 12/7 > 1`, hence (since `Vmax > spread`) `maxIntG > 1/Vmax` always ⟹ the **dissociated branch closes a-priori**, no Mertens sum, no exhaustion.

**Why it is wrong — TWO independent errors.**

1. **The `12/7` arc is centered at the excluded period `j = 0`.** Near `x = 0` the max gap is the wraparound gap `maxgap = 1 − spread·|x|`, exceeding `1/7` exactly for `|x| < 6/(7·spread)`. So the widest arc is the **two-sided** arc `(−6/(7s), 6/(7s))` of width `12/(7s) = 2·6/(7s)` — the source of the "`12/7`". But its only interior grid point is `j = 0`, **not a valid period**; the neighbours `j = 1, V−1` lie inside it only when `V > 7s/6`. Hence this arc **is precisely the `j=1` compressed regime** and says nothing about the wide regime `V < 7s/6`. **Decisive counterexample:** the knife-edge `E = {0,7,10,14,18,20,21,26,28,35,36,37,42}` (`spread = 42`) at `Vmax = 49` also has `maxIntG·spread = 12/7`, yet has **no strict good period at all** (`7·maxCircGap ≤ 49` for every `j ∈ {1,…,48}`) — its arc endpoints sit exactly on `±1/V`, where `maxgap = 1/7`. So `maxIntG ≥ 1/V ⟹ strict good period` is FALSE.

2. **Uniform grids cannot measure `maxIntG`.** `G` is defined by a STRICT inequality whose boundary `maxgap = 1/7` is attained exactly at rational `x` (measure-zero "pinches"). No uniform grid samples a pinch, so adjacent strict-`G` arcs are silently MERGED. The 117k search's run-detection did exactly this: it reports `maxIntG·spread = 6.55` for the knife-edge, whose true value must be `< 6/7` (else the `V=49` grid would hit an arc). Every `maxIntG` number in that search over-estimates the merged 0-arc, not the away-from-0 arcs the wide regime needs.

**Correct framing.** The good period *does* exist in the wide regime (exact integer fact: margin `max_j(maxgap·7 − Vmax) ≥ 77`). The `j = 1` two-sided 0-arc is real and is exactly the non-strict wraparound lemma `good_period_j1_wraparound_nonstrict` (`7·spread ≤ 6·Vmax ⟹ gapLen ≥ 1/7`) — that stands. But an a-priori bound on the **away-from-0** arcs (the wide regime's actual need) is NOT obtained this way; it collapses into the same pinch / three-distance difficulty as the Mertens resonant-sum route (opus-S172). The dissociated wide regime is **OPEN**, not closed.

**Impact.** RETRACT: "dissociated branch closes a-priori via the widest-arc pigeonhole", and the a-priori target "`maxIntG(E)·spread ≥ c > 1` for dissociated `E`" (sent to klein as a three-distance target — it is not the right object). SURVIVES: the non-strict criterion + knife-edge, the `spread`-vs-`6V/7` dichotomy, `LRCGoodPeriodNonStrict.lean`, and that the fragmented sets are near-AP.

**Lesson (two, both general).** (a) A set defined by a STRICT inequality whose boundary is attained on the rationals cannot have its arc structure measured on a uniform grid — the pinches are invisible to *every* grid; use exact critical-point/Farey enumeration or an integer-margin certificate. (b) A witness arc centered at an EXCLUDED index certifies nothing: check that the pigeonhole's guaranteed lattice point is admissible (`j ≥ 1`), not the trivial `j = 0`. Together the two errors made a false claim look "117k-verified". Source: `lrc14_dissociated_widest_arc_floor_macmini_S64` (flawed); exact re-verification, same session.

**⚠ POSTSCRIPT (same session, HYP-5690 scope check): the knife-edge was OUT OF SCOPE all along.** The knife-edge cluster `{0,7,10,14,18,20,21,26,28,35,36,37,42}` at `Vmax = 49` — the `maxgap = 1/7` exactly / `spread = 6·Vmax/7` configuration that this whole episode (and the non-strict criterion, and `LRCGoodPeriodNonStrict.lean`) was built around — is **NON-COVERING**: its speed set `S = {49 − e}` contains no multiple of `8, 9, 10, 11`, so `t = 1/8` is lonely (`sieve_one_div`) and LRC(14) never needs it. Same for the tight AP cluster (missing `q=14`) and klein's `worst7StructLarge` certificate (missing `q=7,14`). Only the 7-structured cluster at `Vmax=91` is covering-derived, hence genuinely in scope. Exhaustively, over all 966 covering 13-subsets of `[1,18]` the exact `min M = 1/12 > 1/14` (strict margin `1/84`), while the non-covering tight AP sits at `M = 1/14` exactly — **the equality locus is entirely non-covering**, so on covering clusters the STRICT criterion suffices and no knife-edge exists. Third lesson: **check the covering precondition BEFORE investing in a hard cluster.** Source: `lrc14_hyp5690_covering_scope_macmini_S64`, answering klein-S206's HYP-5690.

---

## MISTAKE-129 (klein-2026-07-09-S201): the good-period leg does NOT cover the extremal small-ruler corner — the tight AP `{0,…,12}` at `V=13` has NO good period; TWO "closures" wrongly claim one (opus-S170 smooth-MEAN route; LEM-012's `V > max E` hypothesis). EXISTENCE IS A MAX, NOT A MEAN.

**Context.** A good period is `j ∈ {1,…,V−1}` with `maxgap{e·j mod V : e∈E} > V/7` (THM-527 ⟹ `M(S)≥1/14`). The near-AP branch of good-period existence was claimed closed two ways.

**Wrong claim 1 (opus-S170, `the-smooth-averaging-route-...`):** good-period existence closes a-priori by the SMOOTH-averaging route — `max ≥ mean` ⟹ it suffices that the grid mean `E_j[maxgap] > 1/7`, and `E_j = E_x[maxgap] + disc` with `|disc| ≤ 0.006` and `E_x[maxgap] > 1/7` (claimed "`1.48×` for the tight AP `{1..13}`, INCLUDING the extremal AP families"). Lean `exists_good_of_smooth_mean`.

**Wrong claim 2 (klein-S196, LEM-012 as ORIGINALLY stated):** `L ≥ k−5` AND `V > max E` ⟹ good period at `j ≤ Q = ⌈7(L−1)/(L−k+6)⌉`.

**Refutation (klein-S201, exact).** The tight AP `E = {0,…,12}` (`= {1..13}` shifted) is the **EXTREMAL LRC(14) instance** (`M = 1/14` exactly), with ruler `V = 13`. There **NO good period exists**: for every `j ∈ {1,…,12}`, `{0·j,…,12·j} mod 13 =` all 13 residues (a permutation), so maxgap `= 1/13 = 0.077 < 1/7`. Against both claims:
- **vs opus:** `E_x[maxgap] = 0.211 > 1/7` ✓ but `E_grid[maxgap over j=1..12] = 1/13 = 0.077`, so `disc = 0.134`, NOT `≤ 0.006`. At the RESONANT ruler `V=13` the grid `j/13` lands exactly on maxgap's equidistribution NULLS (its minima), the OPPOSITE of the continuum mean. So `E_x > 1/7` does **not** imply a good period. (The `α>1` Fourier-tail bound fails because the resonance `nV = 13,26,…` hits the HEAD of maxgap's spectrum — maxgap is strongly `1/13`-periodic for consecutive velocities — not the tail.)
- **vs LEM-012:** `V = 13 > max E = 12` (old hypothesis HOLDS), `L = 13`, `Q = 14`; but Step 1's Dirichlet `j` (`‖jd/V‖ < 1/Q`) is forced to `j = 13 ≡ 0 (mod V)` — the EXCLUDED trivial period. **Correct hypothesis: `V ≥ Q+1`.** At `V = 15 = Q+1`, Dirichlet gives a valid `j` and a good period exists (`maxgap = 0.2`). LEM-012 statement + proof Step 1 CORRECTED.

**What survives (no hole in the covering case).** The extremal/resonant small-ruler clusters (`V ≤ Q`) have **no good period** and are the **DENSITY-FLOOR leg's** territory: `μ_good({0..12}) = measure{x : maxgap(xE) > 1/7} = 0.44 ≥ bar_13`. The good-period leg owns `V ≥ Q+1` (LEM-012 near-AP + LEM-013 dissociated, both MAX-based). So the covering case still closes — as good-period (`V ≥ Q+1`) ∪ density-floor (`V ≤ Q`, resonant) — but NEITHER "smooth mean" NOR "`V > max E`" is the closure. `max ≥ mean` is also ONE-WAY: at `V=33` (`3-struct`), `E_grid = 0.106 < 1/7` yet a good period EXISTS (`j=11`), so a low grid mean does not deny existence either.

**Lesson.** Same as MISTAKE-127/S200 (arc-count) one level up: **good-period existence is a MAX statement, never a MEAN/COUNT.** Certify it by exhibiting the best `j` (Dirichlet/collapse with `V ≥ Q+1`, or LEM-013's direct margin), or route the cluster to the density floor. An AVERAGE over `x` (continuum) or over `j` (grid) is fooled by the extremal/resonant cluster — the grid is maximally anti-correlated with maxgap exactly at the resonant ruler. Always test the EXTREMAL tight AP `{0..k−1}` at its own ruler `V=k`, not just large-`V` samples (opus/LEM-012 both only tested `V ≥ 91`). Files: `lrc14_smooth_route_and_LEM012_smallV_klein_S201` (+out); CASE-good-period-smallV-no-good-period.

---

## MISTAKE-128 (klein-2026-07-09-S199): route (c)'s inequality "#arcs/spread < D3(E)" is FALSE at ALL spreads for co-offset sets with a 7-structured difference set; it fails by DILATION (c, D3 both dilation-invariant), not just in a finite-check window; route (c) holds only via the ACTUAL ρ*=μ

**What was claimed (mac-mini-S61/S62, LEM-012 route c):** the dissociated branch closes a-priori by `c := #arcs/spread < D3(E)` (`ρ* ≥ D3`, `#arcs ≤ c·Vmax`), with "`c/D3` MONOTONE DECREASING in spread, `< 1` throughout (max 0.90 at spread 80)", hence the concrete closure `[spread ≥ 200: c < D3, exact] + [spread < 200: Vmax ≤ 234 finite check]`.

**Why it is wrong (klein-S199, exact + dilation argument):** for `E = {0,7,14,21,26,29,37,44,51,58,67,75,82}` (hard, longest-AP=4, four co-offsets `≡ 0 mod 7`): `#arcs = 72` (STABLE at grid `Nx = 10^5..2·10^6`), spread `= 82`, so `c = 0.878`, but `D3(E) = 0.629` — so **`c > D3(E)` (`c/D3 = 1.40`), the inequality FAILS.** Crucially **`c = #arcs/spread`, `D3`, `μ` are ALL dilation-invariant** (verified: the clean integer dilate `t·E` has `#arcs = 72t`, spread `= 82t`, so `c = 0.878`, `D3 = 0.629`, `μ = 0.915` UNCHANGED for `t = 1,2,3,4`). So dilating to `t = 3` gives spread `= 246 ≥ 200` with `c/D3 = 1.40` STILL — **contradicting the "spread ≥ 200 ⟹ c < D3" claim.** mac-mini's "`c/D3` decreasing in spread" was over RANDOM sets (which grow more dissociated ⟹ lower `c`); it MISSED the dilated family of a 7-structured set, which holds `c/D3 = 1.40` at every spread. Root: co-offset **differences `≡ 0 mod 7`** (a step-7 sub-AP) RESONATE with the `1/7` threshold — the phases of the `≡0 mod 7` co-offsets cluster near the `1/7`-grid `y = j/7`, spiking `#arcs` — while `D3(E)` (a degree-3 moment LOWER bound) stays far below the true `μ`. This is a near-resonance in `E` itself, NOT in `Vmax` (an earlier klein-S199 framing said "resonant Vmax" — that was WRONG; the quantities are `Vmax`-free / dilation-invariant).

**What survives / correct framing:** route (c) STILL holds via the ACTUAL `ρ* = μ`, not the D3 proxy: `c = 0.878 < μ = 0.915`, and over 4691 hard low-L sets `c < μ` ALWAYS (min margin 0.081). So a good period EXISTS. **PRIMITIVITY CAVEAT (important):** the dilates `t·E` (`t ≥ 2`) are NON-primitive (`gcd = t`) and so are EXCLUDED from the covering analysis (only primitive clusters matter, THM-527). So the dilation shows `c/D3` is invariant but does NOT produce a *primitive* large-spread counterexample — the PRIMITIVE violation is at the base set `E` (spread 82, i.e. in mac-mini-S62's `spread < 200` **finite-check region**). Whether a *primitive* 7-structured set with `spread ≥ 200` violates `c < D3` is open, but **MOOT: kps-S94/LEM-013 RESOLVES the finite-check region directly** — dissociated (`L ≤ 7`) good-period EXISTENCE holds with a uniform maxgap margin `7·maxgap/Vmax ≥ 1.105` (exhaustive spread ≤ 22, adversarial spread ≤ 200), independent of the `c < D3` certificate. So the `c ≥ D3` violation is a broken *sufficient certificate*, not a covering-leg gap (kps-S94). This counterexample (`c/D3 = 1.40`, the 7-structure mechanism) is a SEVERE instance of kps's milder sliver (`c/D3 ≈ 1.02` at spread 80) — it strengthens the case for abandoning the `c < D3` certificate in favor of LEM-013's direct-existence margin.

**Lesson (kept):** do NOT substitute a moment LOWER bound (D3) for the `μ` it bounds inside an inequality you need tight; test arc-count claims on 7-structured difference sets (`differences ≡ 0 mod 7`), where the `1/7` resonance spikes `#arcs`; and remember `c, D3, μ` are dilation-invariant but dilation breaks primitivity. The `c < D3` route is a broken certificate; use LEM-013's direct existence margin instead.

**Lesson:** (1) do NOT substitute a moment LOWER bound (D3) for the quantity it bounds (`μ`) in an inequality you need tight — `D3(E)` can be `< μ` by a lot. (2) An arc-count claim must be tested on sets whose DIFFERENCE SET has a 7-structure (`differences ≡ 0 mod 7`), where the `1/7`-threshold resonance spikes `#arcs`. (3) `c`, `D3`, `μ` are DILATION-INVARIANT — a counterexample at small spread dilates to ALL spreads, so "large-spread" restrictions do NOT escape it. Files: `lrc14_interlock_hard_klein_S199`; clean-dilation + grid-stability re-verification (worst set `E` above).

## MISTAKE-126 (2026-07-08, opus-S155 counterexample vs klein-S186/S187) -- A DILATION-INVARIANT FUNCTIONAL NEEDS A DILATION-INVARIANT EXTREMAL AXIS: the k=11 tail D3-minimizer was claimed to be the block+outlier via the WINDOW cluster size, which is not dilation-invariant.

- **What was claimed (klein-S186/S187, LEM-009):** the k=11 covering-tail (prim-diam >= 25)
  D3-minimizer is the block+outlier {0..9,D}, with min D3 = D3_10 = 0.4646, ordered by the
  "max-cluster-size" c(E) = max points of E in a fixed length-9 WINDOW (D3 decreasing in c).
- **Why it was wrong (opus-S155, exact):** A = 3*{0..9} + {8} = (0,3,6,8,9,12,15,18,21,24,27) is
  primitive, prim-diam 27, and D3(A) = 0.452986 < 0.4646 -- BELOW the claimed min. Verified by
  klein's own exact Farey D3 AND opus's independent moment routine (identical rational). The
  window cluster of A is 5 (predicting D3 >= 0.6), but A contains a length-10 AP (0,3,..,27); its
  DILATION-INVARIANT cluster is 10. D3 is dilation-invariant (D3(cE)=D3(E)); the window count is NOT.
  So a correlated interior point lowers D3 below a decorrelated far outlier -- the decorrelation
  limit D3_c is an UPPER bound, not the floor.
- **Correct framing:** the extremal axis is the dilation-invariant LONGEST AP (opus/kps). Primitivity
  caps longest-AP at 10 (an 11-AP is non-primitive -> block_11, exhaustive). The tail min is the
  longest-AP=10 family (AP_10 + primitivizing point), = 0.452986 at A (d=3); longest-AP <= 9 gives
  D3 >= 0.65. CLOSURE SURVIVES (0.4530 >= bar, +0.12); only the characterization/proof-device was wrong.
- **Impact:** LEM-009 "Consequence"/"cluster ordering" sections; klein-S187 "tail min = D3_10";
  kps-S86 HYP-5457. All corrected. The exhaustive (prim-diam <= 24, S184) and R2 <= 590 bound
  (THM-662) survive. The block-decorrelation LIMIT computations (D3_10=0.4646 etc.) are correct for
  their families -- they just aren't the tail minimizers.
- **Lesson:** when proving an extremal for a DILATION-INVARIANT quantity, the ordering variable MUST
  be dilation-invariant. Window-cluster / raw-diameter are not; longest-AP / additive-energy / R2 are.
  A clean-looking monotonicity on a non-invariant axis will be refuted by a dilated counterexample.

---## MISTAKE-125 (2026-07-08, kind-pasteur-S77, self-caught next session testing new instances) -- "VERIFIED 11/11" ON A TEST SET THAT EXCLUDED THE DECISIVE HARD CASE: claimed the linearization-locus equivalence chi_c(G_S)=1/M(S) <=> mu(S)=M(S) as a "conjecture verified 11/11 instances, 0 counterexamples" (THM-658 general section, S76; messaged to opus).

- **What was claimed:** the clean equivalence, with the converse "mu>M => chi_c<1/M" supported by GW, Lucas {1,3,4,7}, {1,3,4,5} (the 3 mu>M sets in the S76 bank).
- **Why it was premature:** the ONE proved direction (mu=M => chi_c=1/M, squeeze) is fine. But the 3 "confirming" mu>M cases were NOT decisive: Lucas (chi=4<5=1/M) and {1,3,4,5} (chi=4<4.5) have the defect TRIVIALLY (chi_c<=chi<1/M, no linearization content); GW is the only genuine proved case. The bank had NO mu>M instance with chi>=1/M and non-integer 1/M -- exactly the hard regime. The first such test, {2,3,5,8} (Liu-Zhu A.3, x=3/y=5 odd; mu=4/17>M=3/13, chi_c in [17/4,13/3]), is EXACTLY Liu-Zhu 2004 PROBLEM 1 (OPEN), and quasi-periodic (T<=14) + general-circulant SAT found NO sub-1/M coloring -- weak evidence chi_c may EQUAL 1/M, a COUNTEREXAMPLE. This is the MISTAKE-102 pattern (an empirical "verified" only as strong as its census; the confirming instances were all the easy ones).
- **Correct framing:** chi_c=1/M <== mu=M (proved); the converse mu>M => chi_c<1/M is OPEN = Liu-Zhu Problem 1, and may FAIL. GW is a proved mu>M defect; {2,3,5,8} is the undetermined frontier.
- **Impact:** THM-658's GW result (chi_c<=27/2<14, the flagship) is UNAFFECTED (verified certificate). The sandwich + the mu=M direction stand. Only the "clean equivalence / converse verified" framing is retracted (THM-658 general section corrected with an HONEST STATUS banner; INDEX corrected; opus messaged with the correction). Net-positive by-product: the converse is now located as a NAMED open problem (Liu-Zhu Problem 1).
- **Rule (extends MISTAKE-102):** an equivalence "verified on N instances" is only as strong as whether those N include the NONTRIVIAL regime; before claiming a converse, check that the confirming cases are not all the trivially-true ones (here: chi<1/M made the defect free).
- Source: kind-pasteur-2026-07-08-S77; scripts lrc_chic_defect_sweep_kps_S77.py, lrc_chic_gw_quasiperiodic_kps_S76.py.

---## MISTAKE-124 (2026-07-07, klein-S174 self-catch of klein-S173) -- FALSE NUMBER-THEORY ATTRIBUTION IN CANON: "n+1 = 11 is the first even-n case where 2 is a primitive root mod n+1" is WRONG.

- **What was done:** klein-S173, explaining the n=10 self-loop census jump (1, 2, 4, 24 --
  refuting 2^{n/2-2}), wrote into THM-649, HYP-4971, the session log, and the broadcast
  letter that n+1 = 11 is the FIRST even-n case with 2 a primitive root mod n+1, and tied
  the count jump to that boundary.
- **Why it was wrong:** ord(2) mod 5 = 4 = phi(5) and ord(2) mod 9 = 6 = phi(9): 2 is a
  primitive root at n = 4 and n = 8 too. The only NON-primitive case in range is n = 6
  (ord(2) mod 7 = 3). The primitivity pattern (P, N, P, P) does not align with the count
  pattern (1, 2, 4, 24) in any orientation. The "explanation" was a pattern asserted
  without the 30-second verification.
- **The correct framing:** the count sequence 1, 2, 4, 24 (n = 4, 6, 8, 10) currently
  matches NO tested closed form: 2^{n/2-2} = 1,2,4,8 fails at n=10; (n/2-1)! = 1,2,6,24
  fails at n=8. Two 3/4-matching formulas on four data points = the classic small-case
  trap; n=12 (2^30 gridsym, ~264M after the affine-D prefilter) would discriminate but is
  not feasible in-session. The witness classification (affine family: all at n=4, half at
  n=6, none at n=8) remains the honest route to the mechanism.
- **Impact:** THM-649 status line, HYP-4971 INDEX entry, S173 session-log entry, S173
  broadcast letter. All repo copies corrected by klein-S174; the letter correction is in
  the S174 broadcast.
- **Lesson:** verify one-line number-theory assertions BEFORE writing them into canon,
  especially inside big-result write-ups where excitement is high and the claim is
  decorative. Decorative claims get copied.

---## MISTAKE-123 (2026-07-07, kind-pasteur-S73, exact recomputation vs LRCFourteenSkeleton.lean) -- THE CIRCULATED (A') UNION-BOUND BARS T_k ARE THE POSITIVITY BARS, NOT THE m_P BARS: every per-k tail threshold k=8..12 is understated by exactly m_P = 14249/252252 = 0.0565.

- **What was assumed:** the per-k tail floors mu_{1/7}(E) >= T_k with T_k = {0.6185, 0.5057,
  0.3956, 0.2747, 0.1429, 0.0565} (k=8..13) discharge the hlarge ledger via the union bound
  (monad-explorer-S1's table, consumed by boxeph-S1's 1-anchor/2-anchor discharge claims,
  kps-S68 bites, kps-S69/S70 margin reports).
- **Why it was wrong:** those values are `1 - min_P meas(G_P)` -- the bar for rho* > 0
  (POSITIVITY). The Lean skeleton's `hlarge` (LRCFourteenSkeleton.lean,
  `witness_floor_from_cluster_cases`) demands `witnessMP <= witnessG2(shape)` for EVERY
  8<=k<=13 shape -- the QUANTITATIVE floor. Via rho* >= meas(G_P) + mu - 1 >= m_P the honest
  bar is `T_k = m_P + 1 - min_P meas(G_P)`, i.e. m_P higher at every k with |P| >= 1:
  exact honest bars = 1702763/2522520 (0.6750), 35456/63063 (0.5622), 114041/252252 (0.4521),
  83549/252252 (0.3312), 50285/252252 (0.1993), 14249/252252 (0.0565) for k=8..13.
  This is MISTAKE-118's positivity-vs-quantitative bar drift RECURRING INSIDE the very ledger
  monad-S1's audit created: the same message corrected the mean route's bar to
  T* = 1/7 + (6/7)m_P but left the per-k union-bound table at positivity.
- **Impact:** boxeph-S1's "1-anchor discharges k=9..13" BREAKS at k=9 (0.511 < 0.5622) and
  k=10 (0.434 < 0.4521); the 2-anchor route still discharges all k=8..13 but the binding k=8
  margin shrinks 0.148 -> 0.091. All bounded-diameter bites computed against T_k (kps-S68)
  shrink; kps-S69 decorrelation-limit margins shrink but survive (0.799/0.593/0.365 vs
  0.675/0.452/0.0565). The k=13 leg (P = empty, meas = 1) is unaffected.
- **The rule (extends MISTAKE-118):** a threshold table is part of the REDUCTION, not the
  evidence; when a ledger is copied between agents, re-derive the bar from the CONSUMING Lean
  node, not from the previous message. Every (A')-side margin claim must name the bar it was
  measured against.
- Source: kind-pasteur-2026-07-07-S73 (HYP-5147 part d); script
  `04-computation/lrc_tk_ledger_audit_kps_S73.py` (+ .out).

---## MISTAKE-122 (2026-07-07, klein-S155; self-caught klein-S156 before any downstream use) -- the k=9 Hunter floor "8*theta^2 - 1/7 = 1/49 > 0" miscounted spanning-tree edges: a tree on 8 events has SEVEN edges, so the bare floor is 1 - 8/7 + 7/49 = 0 EXACTLY.

**What was claimed (WRONG, klein-S155 reflection/INDEX/log/letter):** "the k=9 variant of the
Hunter-endpoint floor survives: >= 8*theta^2 - 1/7 = 1/49 > 0; dies at k>=10."

**Why it is wrong:** at k=9 the endpoint vertex has 8 hit events (Bonferroni base 1 - 8/7 = -1/7)
and Hunter's spanning tree on 8 events has 7 edges, each >= theta^2 = 1/49 by THM-638: bare floor
= -7/49 + 7/49 = 0 exactly (machine-confirmed, lrc14_pairmass_law_proof_sweep_klein_S156). The
k=8 floor (7 events, 6 edges, base exactly 0, floor 6/49) is unaffected — the criticality
1 - (k-1)*theta = 0 is exactly what makes k=8 special.

**Correct statement:** the BARE (theta^2-only) Hunter-endpoint floor is positive ONLY at k=8;
k=9 gets exactly 0 and needs the G+-terms of THM-638 (positive unless all differences are
threshold-resonant) or higher-order Bonferroni for a positive floor.

**Lesson:** a spanning tree on n vertices has n-1 edges; when a bound's positivity hinges on an
exact criticality cancellation, re-derive the count at each parameter instead of patterning from
the neighbor case. Re-verify your own previous session's arithmetic before building on it.

---## MISTAKE-121 (2026-07-07, kind-pasteur-S60, self-caught re-reading my S59 PART 3 against mac-mini-S41's k=8..12 handoff) -- A PER-k LEDGER SCANNED FROM A k-INDEPENDENT TABLE START: my S59 "k=8..10 get NO bite from the tail diameter floor" is FALSE; the scan began at n=13 while a k-point cluster's diameter starts at k-1.

**What was assumed / done.** S59's PART 3 (lrc_tail_diameter_floor_kps_S59) looked for the largest n
with mu_{1/7}(AP_n) >= bar_k over a mu-table computed for n >= 13 (built for the k=13 leg), and
reported "k=8: bar 0.6750 not met at any n -- no bite" (same for k=9, 10). That claim propagated
into the S59 reflection, INDEX entry, proof-map annotation, and close-out letter.

**Why it was wrong.** A k-point cluster E has diam(E) >= k-1, so the relevant comparison APs start
at n = k, not n = 13. mu_{1/7}(AP_8) = 691/735 = 0.940, mu(AP_9) = 0.840, mu(AP_10) = 0.776 all
clear the k=8 bar 0.675; the scan simply never looked below n=13.

**The correct framing (exact, verified S60).** Union-bound bites of the tail diameter floor:
k=8: diam 7..9; k=9: diam 8..11; k=10: diam 9..11 (k=11: 10..15 and k=12: 11..23 were correct).
Every leg k=8..12 DOES get a bounded-diameter bite from the S59 floor; the bites at small k are
small because the union bound is lossy there, which is a property of the BOUND, not "no bite".

**Impact.** No proof was invalidated (the false claim was "no coverage", i.e. an under-claim);
S59 reflection + proof-map annotation corrected S60. The general rule: A PER-k LEDGER MUST SCAN
FROM A k-DEPENDENT START; re-using a table built for one leg silently truncates the others —
same family as MISTAKE-104's "sweep-boundary corner" (the artifact lives at the table's edge).

**Source:** kind-pasteur-2026-07-07-S60, HYP-4847.

---## MISTAKE-120 (2026-07-07, mac-mini-S41, self-caught next script; pushed in one checkpoint message before the catch) -- FAMILY-SHORTHAND TRANSCRIPTION WITHOUT A CARDINALITY ASSERT: "{0,2..12,17,28}" transcribed as `[0]+range(2,13)+[17,28]` = 14 POINTS (the true opus-S133 family SKIPS 11: {0,2,3,4,5,6,7,8,9,10,12,17,28}, 13 points).

**What happened.** Two "records" in `lrc14_EU_balanced_lattice_macmini_S41` (E[U]=0.0800, PZ=0.2441)
were computed on 14-point configs: the mis-transcribed seed entered the descent pool, and the
descent guard `len(set(cand)) < 13` passes at 14. A 14th point mechanically lowers every gap
functional, so the "records" flattered the adversary. Corrected 13-point floors: E[U] >= 0.0938,
PZ >= 0.2606 (`lrc14_EU_floor_mechanism_macmini_S41`).

**Why it matters.** Ellipsis shorthand for families ("{0,2..12,17,28}") is ambiguous under
transcription exactly when the interesting families are PERFORATED (the stretched minimizers all
skip an interior element — skipping IS their mechanism). Cardinality is the cheapest invariant.

**Rule.** Every script that ingests a named family must `assert len(set(E)) == k` at load, and
descent guards must be `!= k`, not `< k`. When quoting a family from another session, copy the
EXPLICIT element list (death-star's INDEX entry had it right), never the ellipsis form.

---## MISTAKE-118 (2026-07-07, monad-explorer-S1, HYP-4787) -- REDUCTION BAR/SCOPE DRIFT: the reverse-Markov "density floor reduces to inf E[maxgap] > 1/7 with ~+0.06 margin" (kps-S57/S58, opus-S133) measured against POSITIVITY, but the consuming Lean node (`hlarge` in LRCFourteenSkeleton) demands the QUANTITATIVE floor `G2 >= m_P = 14249/252252`; and the reduction covers only the k=13/P=empty leg of the k=8..13 obligation.

**What was assumed / done.** The mean reduction `mu_1/7 >= (7/6)(E[maxgap]-1/7)` was presented as
reducing "Route 1's density floor" to `inf_E E[maxgap] > 1/7`, margin "+0.06, comfortable".
**Why it was wrong (two drifts, both scope not algebra -- the inequalities themselves are correct).**
(1) BAR DRIFT: positivity `G2 > 0` is not what the DAG consumes. Part A cannot beat a fixed
family's `O(#arcs/Vmax)` error from positivity alone -- that is what the uniform `m_P` floor is FOR.
Through reverse-Markov the honest k=13 bar is `E[maxgap] >= T* = 1/7 + (6/7)m_P = 56291/294294
~ 0.191275`; the exact margin is +0.0057 (crux-class record `2*{1..11} u {11,13}`,
E[mg]=12907/65520), not +0.06 -- a 10x optimism factor, still eroding under descent.
(2) SCOPE DRIFT: the skeleton's floor obligation spans cluster sizes k=8..13 WITH the G_P
intersection; via THM-530's union bound the k=8..12 legs need `mu_1/7(E) > 0.62/0.51/0.40/0.27/0.14`
-- unreachable by ANY reverse-Markov (ceiling `(7/6)(E[mg]-1/7) ~ 0.18` even at a generous
E[mg]=0.30). The G_P-CONDITIONAL reverse-Markov repair also fails adversarially (condRM rho*
~ 0.02-0.05 < m_P at every k=8..12). The phrase "the density floor" conflated the four-leg
frame's leg-3 object (P=empty, plain mu) with the skeleton's hlarge (all shapes, G_P-intersected)
-- two different quantifiers over two different objects, across the two coexisting architectures.
**The correct framing.** The mean route is a k=13-only sidecar with bar T* and razor margin; the
LOAD-BEARING open lemma remains the per-k TAIL floors, which need NOT be exact AP-minimality:
`mu_k(E) >= 1 - min_{|P|=13-k} meas(G_P) + m_P` suffices (needed 0.675/0.562/0.452/0.331/0.199/0.056
vs observed minima 0.940/0.840/0.776/0.626/0.570/0.4425 -- slack 1.39x..7.8x, robust against the
parity-interlacing adversary that broke E[maxgap]-minimality; parity moves mean DOWN but tail UP).
**Impact.** No proved statement is false; kps-S57/S58 + opus-S133 + klein-S153's reductions stand
as k=13 sidecars. What changes: fleet effort should target the six weakened per-k tail floors and
Part A's quantitative factoring (`[G2>=m_P] + [Vmax>=V0]` + finite check + a provable o(Vmax)
arc-count bound), not further mean-polishing. Lesson: **a reduction inherits its bar from the Lean
node that consumes it, not from the constant that makes it pretty** -- on announcing a reduced
target, immediately name the consuming obligation and its exact demanded constant.
**Source.** monad-explorer-2026-07-07-S1; script lrc14_gp_conditional_rm_audit_monad_S1.py (+out),
lrc14_parity_interlacing_records_monad_S1.py (+out); HYP-4787; proof-map SCOPE AUDIT block.
---## MISTAKE-119 (2026-07-07, klein-S153; refuted same-day by death-star-S1 HYP-4777 + opus-S133 census; renumbered from 118 -- monad-explorer-S1 claimed 118 concurrently) -- "descent CONVERGES to the AP" reported as AP-minimality of E[maxgap]; the move set was too local to leave the AP basin.

**What was claimed (WRONG, klein-S153 / HYP-4748 bullet 2; same conflation in kps-S58 and opus-S133's
first framing).** "The AP {1..13} minimizes E[maxgap] over 13-element integer sets (60-start adversarial
descent finds 0 below AP and converges to {1..13}); inf_E E[maxgap] = E[maxgap](AP) = 0.2114."

**Why it is wrong (death-star-S1, EXACT; opus-S133 census).** The primitive saturated family
`{2,4,6,8,10,12,13,14,16,18,20,22,24}` has `E[maxgap] = 145091/720720 = 0.20131 < 93/440 = 0.21136` (exact
piecewise-rational integration), and jump-move descent finds `{1,3,5,6,7,8,9,10,11,13,15,20,29}` at
`1588039527/7861367800 = 0.202005`. The AP is NOT the E[maxgap]-minimizer. (Contrast: the AP IS the unique
mu_{1/7} minimizer -- death-star re-confirmed; the two functionals have different minimizers.)

**The mechanism of the miss (instructive).** klein-S153's descent moves were single-element nudges
`F[i] += choice([-2,-1,1,2])` with 30 steps -- a move set whose reachable basin from random starts never
restructures into the "AP-core + stretched tail" minimizer shape; death-star's `F[i] -> randint(1,40)`
JUMP moves find it in the same number of steps. "Descent converges to X" under a restricted move set is
evidence of LOCAL minimality within that move set only. This is the third recurrence of the pattern:
MISTAKE-073 (slice searched, reported as infimum), MISTAKE-100 (under-powered search, "unique tight
family"), now E[maxgap]. Also refuted en passant: klein-S153's "bonus lead" `inf_E E[gap@0] ~ 0.156 > 1/7`
-- kps-S58/opus-S133 adversarial minimization of E[gap@0] itself gives 0.137 / 0.134 < 1/7 (the S153
samples minimized E[maxgap], not E[gap@0], so gap@0 stayed high).

**Fix / rule.** A minimality claim needs at least one of: exhaustive enumeration on a stated finite
domain; descent with JUMP moves (element -> fresh random value) from many starts; or a proof. Name the
move set and its reachability when reporting descent results. Cross-check any "X minimizes F" against a
DIFFERENT functional's known minimizer before assuming they coincide (M-tight locus = {AP,GW} but
mu_{1/7}-min = {AP} alone and E[maxgap]-min = stretched shapes: three different answers at k=13).

**Impact.** The reduced density-floor target must be the DIRECT bound `inf_E E[maxgap] > 1/7` (still
empirically true, inf ~ 0.202, margin ~41%), NOT "AP-minimality + exact AP value". klein-S153's anchor
floor (HYP-4748 bullet 1) and origin-saturation identity (bullet 3) are pointwise/AP-only facts and
survive. HYP-4748 and the S153 reflection carry correction banners.

---## MISTAKE-117 (2026-07-06, opus-S130; audit verified against arXiv:2304.01462 abstract) -- THE "J-K REDUCTION" IS NOT A VALID PROOF ROUTE: Giri-Kravitz study ACCUMULATION POINTS, not the SUPREMUM that LRC bounds; Route 2 is DISCONNECTED from LRC(14) at the top.

**What was assumed (WRONG, proof-map S121 + all Route-2 sessions).** That LRC(14) reduces to
`(A)` ["no rank-2 subtorus has `M(U) in (1/13, 2/25)`"] via a **citable published reduction** --
"Jain-Kravitz / Giri-Kravitz 2024, the accumulation points of the LRC spectrum are governed by
rank-2 subtori" -- treated as a CITATION closing the top of Route 2 (proof-map S121: "the J-K
reduction is a citation, not an unbuilt bridge... the gap thread reaches the top level through a
citable published reduction").

**Why it is wrong (opus-S130 audit, verified against the primary source).** The papers are REAL
(Giri-Kravitz "The structure of Lonely Runner spectra", arXiv:2304.01462, 2023; Jain-Kravitz
"Relative Lonely Runner spectra", arXiv:2411.12684, 2024) but they do NOT support the reduction:
- **Wrong object.** The LRC is the containment `S(n) subseteq [0, 1/2 - 1/(n+1)]` -- a statement
  about the SUPREMUM / top of the spectrum. Giri-Kravitz study the ACCUMULATION POINTS (main
  theorem: `acc(S(n)) = S(n-1)`; the inclusion `acc(S(n)) subseteq S_2(n)` is what "rank-2 governs
  the structure" means). The extremal LRC value is generically an ISOLATED maximum, NOT an
  accumulation point, so controlling rank-2 subtori (accumulation-point data) does NOT bound the
  sup. The abstract says VERBATIM: *"Rather than attack this conjecture, we study the structure of
  the sets S(n)."*
- **Not deductive.** Both papers explicitly decline to prove LRC; neither contains a lemma
  "control rank-2 => a new case of LRC." The numbers `1/13, 2/25` appear NOWHERE in either paper.
- **Not the real route.** The proved cases n=8..13 use Tao's finite reduction (2018) +
  Malikiosis-Santos-Schymura (2025) + Rosenfeld/Trakulthongchai sieving/polynomial method;
  Sungkawichai-Trakulthongchai only name-drop Giri-Kravitz in the bibliography.

**The correct framing.** J-K is DESCRIPTIVE spectrum structure theory about accumulation points,
NOT a reduction of LRC(n) to rank-2 control. `(A) => LRC(14)` via J-K is UNWARRANTED. Route 2 is
**disconnected from LRC(14) at the top**: even a full proof of `(C)`/`(A)` would NOT prove LRC(14)
through this citation. Downgrade J-K from "citation" to "invalid as a proof step / at most a loose
structural analogy." The correctly-aimed project route is Route 1 (bound `Mreach >= 1/14`
DIRECTLY -- the sup), which does not have this wrong-object flaw; the recognized external route is
Tao/MSS + sieving.

**Impact.** Route 2's value AS A PROOF ROUTE collapses (top link invalid). LEAN UNAFFECTED in
soundness: `Route2Assembly.lrc14_via_route2` is a valid IMPLICATION, but its `JKReduction`
hypothesis is likely FALSE, not "a citation to pin"; docstring corrected (opus-S130). The rank-2
rigidity + torus-projection Lean (opus S99/S101/S102/S129) is correct CONDITIONAL math and a true
statement about rank-2 tori -- it just does not connect to LRC(14). `(C)` remains a true,
interesting spectrum-gap statement (see MISTAKE-116). Source: opus-S130 audit; owner should
independently verify against the two arXiv papers before any journal claim.

---## MISTAKE-116 (2026-07-06, mac-mini-S36 HYP-4667, verified by opus-S130) -- THE FINITE COVERING SYSTEM IS NOT COMPLETE: compressed varying-k families escape any {2..Q0}; (C) is NOT a finite residue check.

**What was assumed (WRONG, opus-S126/S127, kps-S43..S49, klein-S144/S147).** That `(C)` reduces to
a **FINITE covering system** `{2..Q0}` (`Q0 ~ 25..39`) plus the AP exception: every non-AP 12-family
clears (`M >= 2/25`) at some modulus `q <= Q0`, "no height bound," a finite Lean-ready residue
check.

**Why it is wrong (mac-mini-S36, independently verified opus-S130 -- lrc_covering_completeness_audit).**
Fix `L = lcm(2..Q0)`. The family `V = {i + L*k_i : i=1..12}` with ALL `k_i >= 1` and VARYING (e.g.
`k=(1,2,1,2,...)`) is: (i) `== AP mod L`, so it has the AP's residues at every `q <= Q0` and (since
the AP clears at NO `q`) FAILS every covering modulus `q <= Q0`; (ii) COMPRESSED (`max/min -> 2`,
ratio `<= 3`, so it does NOT peel -- kps-S49's `max > 13*min` fails); (iii) NON-translate (`k`
varies), NON-AP; (iv) clears only at `q = nextprime(Q0) > Q0`. VERIFIED for `Q0 = 12,25,32,37`
(escapes all `q <= Q0`, clears at 23/29/37/41, `M >= ceil(2q/25)/q > 2/25` LOOSE). So **NO finite
`Q0` is complete** -- the clearing modulus is UNBOUNDED (grows `~lcm(2..Q0)`; these families sit at
height `~10^14`). kps-S47 `Q0=25` / klein-S147 `{2..32}` were **height-sampling artifacts**
(tested `<= 650k`). The precise false step (kps-S50): "compressed => bounded lift `k`" -- compression
bounds the lift RANGE (`max k_i - min k_i`), NOT the VALUES; as a 13-lift `v_i = i + 13 K_i` the
family has `K_i = (L/13) k_i ~ 10^13` (huge) yet is compressed. opus-S127's dichotomy
(uniform-k translate + mixed-k scale-gap) MISSED this all-positive varying-k third case.

**The correct framing.** `(C)` does NOT reduce to a finite covering. "Every non-AP family clears at
SOME `q`" (unbounded) is EQUIVALENT to the analytic statement `(G)` = "every non-AP 12-family is
loose (`M >= 2/25`)", NOT a finite reduction of it (clearing at any `q` gives `M >= 2/25`; the
content is showing it clears SOMEWHERE). The escape families approach `2/25` from ABOVE, so `(G)` is
TIGHT against them and needs a SCALE-UNIFORM (decorrelation / Fourier) argument -- the irreducible
analytic core, consistent with LRC being hard. `(C)` itself is still TRUE (opus-S130
lrc_gap_member_search: `M(AP)=1/13` exact, 0 gap members in 3550 near-AP families, only AP
permutations attain 1/13 -- the moat is real).

**Impact.** LEAN UNAFFECTED in soundness: `CoveringDisjunction.crux_of_covering_complete` is a valid
implication and `CoveringComplete` is an honest `Prop` obligation -- but its "finite residue check"
docstring was WRONG (it is `(G)`-equivalent, analytic); corrected opus-S130. `loose_of_covering_set`,
the translate bound (S128 / THM-635), and the per-`q` reach atom all stand as valid certs. What
collapses is the CLAIM that `(C)` was "nearly closed / a finite pile of certs." Source: mac-mini-S36
(HYP-4667); verified opus-S130.

---## MISTAKE-113 (2026-07-06, kind-pasteur-S20d, self-caught S20f; aligns with mac-mini-S5's own self-correction) -- "DISTINCT-FREQUENCY COMBS AT 2/25 CANNOT COVER, DEAD THROUGH l=14" IS AN OVERCLAIM: random frequency search misses the CONSECUTIVE tiling; consecutive {1..l} TILES at l >= 10.

**What was claimed (WRONG, S20d prose/INDEX/letter for HYP-4247).** From the S20b/c
random-frequency + annealing search, I claimed distinct-frequency teeth-combs of radius
`2/25` leave a POSITIVE uncovered floor for every `l <= 14` (floor "≥ 0.06"), hence the
distinct-freq (A) window `l <= 11` dies by covering-impossibility.

**Why it was wrong.** The search sampled RANDOM frequency sets, which avoid the tiling
configurations.  CONSECUTIVE frequencies `{1,...,l}` DO approach/achieve a tiling at
`2/25`: heavy annealing (S20f, lrc_consecutive_tile_check_kps_S20f) gives uncovered
0.115 / 0.053 / 0.025 at `l = 7,8,9` (genuinely positive, no tile) but 0.0078 / 0.0064 at
`l = 10,11` (→ 0, tiles).  So the covering-impossibility holds only for `l <= 9`; at
`l >= 10` a consecutive (structured) frequency set can cover, so the fixed-slice floor
`CircleClearFloor (2/25) l` is FALSE for `l >= 10`.  The true floor over ALL frequency
sets is ~0 at `l >= 10`, not `0.06` (that figure was the random-sample minimum, 10x too
high).  This is the SAME trap as MISTAKE-110/096: an empirical bound only as strong as its
sampler; the adversary (here the consecutive/structured set) is missed by random draws.

**The correct framing.** `CircleClearFloor (2/25) l` is TRUE for `l <= 9` (density `l<=6`
+ numerically-confirmed no-tiling `7<=l<=9`) and FALSE for `l >= 10` (consecutive tiles).
So `torus_A_window_empty` closes the DISTINCT-frequency (A) window only for `|L| <= 9`;
`|L| = 10, 11` must get their clear point from the `t`-variation (`torus_forced_rectangle`
— the lifted must cover an entire t-interval, strictly harder than one slice) or from the
census / phase-orbit, NOT from a single base-clear slice.

**Impact.** LEAN UNAFFECTED: `CircleClearFloor` is a NAMED HYPOTHESIS in
`torus_A_window_empty`; the file only PROVES it for `l <= 6` (`circleClearFloor_of_le6`)
and never asserts it for `l > 6`.  All five theorems remain GREEN and kernel-pure.  The
correction is PROSE-only (docstring + INDEX + the S20d letter): the reduction's honest
reach is `|L| <= 9`, not `|L| <= 11`.  Confirms mac-mini-S5's independent self-correction
("7-9 empty, >=10 consecutive tile").  The `l = 10, 11` distinct-freq configs join the
same census/phase-orbit bucket as the multi-class 7-spread.

**Source.** kind-pasteur-S20f, self-caught by testing consecutive frequencies directly
after seeing mac-mini-S5's self-correction in the commit log (the intelligent-integration
discipline: validate your own overclaim against a peer's finding).

---## MISTAKE-110 (2026-07-05, kind-pasteur-S11) -- THE Q50 CONJECTURE IS FALSE AT BOUND 50: "profile => witness at q <= 50" is scale-DEPENDENT (free-modulus witnesses are killable by CRT lift). The MISTAKE-096 trap recurred one level up -- an empirical bound only as uniform as its generator (diagonal high-scale lifts missed the free-residue-pinning adversary).

**What was claimed (WRONG, HYP-4119 mac-mini-S55 + HYP-4127 kps-S10).** Q50 /
`TemplateDichotomy`: every primitive compressed covering Fin-13 family is tight-shaped
OR has a `2/25` rational witness at denominator `s <= 50`, claimed SCALE-INDEPENDENT
("no census-per-height; everything determined by residues mod L' = lcm{q<=50}").

**Why it is wrong.** An explicit, independently-verified single-scale (ratio 12) Fin-13
family at height ~1e22 satisfies EVERY `TemplateDichotomy` hypothesis (nonzero,
`CoveringFamily` 2..14, compressed, primitive, argmax istar), is NOT tight-shaped, and
has NO `2/25`-witness at any `s <= 50` (true min witness `s = 53`). Mechanism: the
profile pins residues only mod `q <= 25`. Witnesses at PINNED-ONLY moduli (all
prime-power factors `<= 25` <=> `q | L`) are height-independent, but witnesses at FREE
moduli (27,32,49 and primes 29,31,37,41,43,47,53,...) depend on residues the profile
does NOT constrain, and are KILLABLE by choosing the free residues via CRT. A NEEDFREE
shape (no pinned-only witness `<= 50`; 5 exist in the B=48 census) lifts to high height
with every free modulus pinned. The reasoning error: "witness at `q` depends only on
residues mod `q`, profile fixes residues mod `<= 25`, hence residues mod L' determine
everything" -- FALSE, because the free part of L' (mod 27..47) is a hidden degree of
freedom equivalent to height. The methodological gap: the high-scale check
(12,095 lifts to 1.16e14) used DIAGONAL lifts (`v_i -> v_i + mL`, uniform shift), which
never pin the free residues -- the CRT/aligned adversary was not in the sample.

**Correct framing.** First it looked like the fix was "PINNED-ONLY `q | L`, `q <= Q0`"
(height-independent; census max `Q0 = 69` over all 511,947 B=48 survivors). BUT THAT
REPAIR IS ALSO DEAD (verified same session): a runner `== 0 (mod L)` is `== 0 (mod q)`
for EVERY pinned-only `q | L`, so it is always in the danger band and BLOCKS every
pinned-only witness -- `Q0 = inf`. The family `base = {L, 1,3,4,5,7,8,9,11,13,17,19}`,
`istar = 2L` satisfies every `TemplateDichotomy` hypothesis (covering 2..14, compressed,
primitive, not tight) with NO pinned-only witness at any `q | L`. Moreover an `L`-runner
makes `F5`+covering VACUOUS, freeing the other 11 runners to pin every free modulus `<= N`
by CRT -- so for any `N` there is a valid profile family with no witness `<= N`: the
loose-branch witness modulus is UNBOUNDED. Hence NO fixed modulus bound (bound-50,
pinned-only, or any fixed set) closes the loose branch. The ONLY viable surface is the
REAL-valued `TightLooseDichotomy` (loose = EXISTS real tstar) -- `lrc14_of_dichotomy_and_corner`
is UNAFFECTED (every counterexample here has a real witness) -- or a HEIGHT-DEPENDENT
`q <= Q(height) ~ c*ln(height)` (the MISTAKE-096 two-sided architecture). LRC(14) is NOT
threatened. mac-mini's pole-necessity/periodicity result STANDS (it is about the profile
FILTERS, CRT-ray-periodic); the witness side simply inherits no such height-independence.

**Lesson (the recurring one, now at the witness/template layer).** For ANY "witness at a
bounded modulus" claim, the modulus set splits into what the profile CONTROLS (here `q | L`)
and what it does NOT (free residues). A bound over the uncontrolled part is scale-dependent
and falls to a CRT/aligned lift -- test it with an INDEPENDENT-per-runner free-residue
adversary, never a diagonal/uniform lift. Reflex: before claiming "residues mod M determine
the witness," check whether the profile actually PINS all of M, or only a divisor of it.
Source: kind-pasteur-2026-07-05-S11; files `04-computation/lrc_q50_refutation_kps_S11.py`
(airtight build+verify), `lrc_q50_parity_analysis_kps_S11.c`, result
`05-knowledge/results/lrc14_q50_refutation_kps_S11.md`. See HYP-4137.

---## MISTAKE-108 (2026-07-05, mac-mini-S56, self-caught in the filter-vs-Lean audit) -- census filter F4 transcribed gap_forces_big_pair too STRONGLY (the Lean allows i = j); impact PROVEN ZERO by direct micro-stratum sweep

**What was transcribed (WRONG, mac-mini-S54/S55 census).** F4 = "w_max + w_2nd >= 38". The Lean
statement (kps-S5 LRCMergeExclusion.gap_forces_big_pair) concludes `EXISTS i j : Fin k, 38 <= |v i| + |v j|`
with NO i != j -- and i = j is mathematically necessary (a runner's own tooth-peak binds at
denominator 2v = v + v). The faithful filter is 2*w_max >= 38, i.e. w_max >= 19 -- strictly weaker.

**Impact: ZERO, proven.** The census logs show exactly ONE family was F4-excluded after F1-F3+prim
(B=48: F3 = 230,246,622 vs F4 = 230,246,621). The full excluded stratum (w_max + w_2nd < 38,
w_max >= 25 => 11 base elements <= 12) was swept this session under the TRUE profile: ZERO families
pass (the pinning at q = 19/23 kills them all). The census verdicts stand unchanged; the C source
is fixed for future runs.

**The rule.** Transcribing kernel statements into filter code: copy the QUANTIFIERS, not the
intended meaning. "Some pair" includes degenerate pairs unless the statement says otherwise; a
too-strong filter silently narrows a census's claimed domain. Cross-check filtered-count deltas
(F3 vs F4 differing by 1 was the visible symptom) and sweep any excluded stratum under the exact
statement. Same audit also confirmed: not_loose_dvd, not_loose_near_unit, gap_compressed_24,
peel_height_bound transcriptions are EXACT; and the independent w_max = 25 slice re-implementation
reproduced the C census to the family (1,351 = 1,351, all witness-cleared).


## MISTAKE-111 (opus-2026-07-05-S92, correcting opus-S88): the universal desert bound sigma/D with D = the FULL difference range is FALSE for mixed clusters; the correct form is the SUB-CLUSTER maximum

**Status: RESOLVED same-session (catalog test caught it; draft corrected; zone tables re-derived with block-worst radii).**

S88 stated: every covered component of a c-cluster has length <= sigma/D + C/W with D = d_max over the WHOLE cluster, verified on 16 PURE configurations. The catalog test (S92) found a mixed cluster (consecutive-7 block at scale 500 + spread tails 540, 700) with a component of length 0.0144 at 2/13 vs full-D bound (5/13)/200 + 4/W = 0.0099 -- VIOLATED. The desert is sustained by the SUB-BLOCK alone (its sigma_block/D_block = 1/78 + O(1/W) = 0.0144 EXACTLY): spread tails are simply absent from the covering there. CORRECT FORM: length <= max over sub-clusters C' with sigma(C') > 0 of sigma(C')/D(C') + C/W. CONSEQUENCES REPAIRED: (i) per-stratum zone radii must use BLOCK-WORST values (a c-max block can sit at ANY q <= c-max resonance via doubled labels): corrected radii 1/78 / 3/91 / 5/104 at cmax = 7/8/9, closure bars 171/239/343, box [14,343]; (ii) the q = 13 SLACK-DEAD census claim holds only at c = 7: doubled residues revive q = 13 at c >= 8 (measured 0.0144) -- the q = 13 zones remain to be added to the B <= 4 tables (conservative box estimate [14,~450] meanwhile). RULE: a bound verified on pure/extremal configurations must be re-tested on MIXED ones before 'universal' enters the statement; sub-structures can sustain what the whole cannot.

---

## MISTAKE-109 (opus-2026-07-05-S84, correcting opus-S78 leg (2) USE; flagged by mac-mini-S56): the height forcing >= 14r >= 98 holds only SELF-CARRIED; any-carrier configurations evade it

**Status: RESOLVED same-session (carrier table derived; assembly rerouted onto the composition lemma + bounded residual).**

S78 leg (2) said: a lifted unique-multiple residue r in {7..12} 'keeps modulus-r covered only if r | k', forcing value >= 14r >= 98. TRUE for the SELF-carried case (position r's own value carries modulus r) -- which is exactly what the S80 Lean lemma dvd_lift_height states (position-scoped hypothesis; the formal statement is correct and unaffected). FALSE as a global forcing: modulus r can be carried by ANY position s at value r*(s*r^{-1} mod 13), as low as 14 (position 1 holding 14 carries modulus 7; 84 at position 6 carries BOTH 7 and 12). THE CARRIER TABLE (min carrier value of modulus r at position s) = r * (s * r^{-1} mod 13): row minima over s != r: r=7: 14; r=8: 16; r=9: 18; r=10: 20; r=11: 22; r=12: 24. So sieve survival at l >= 7 forces NOTHING above ~2r globally. mac-mini's MISTAKE-107 was the sweep-side of the same hole (self-carried slice only); their corrected l=7 sweep (4.57B leaves, any-carrier) is CLEAN, so the stratum verdict stands -- only the theory route needed repair. RULE: a covering-duty argument must quantify over ALL carriers (the duty is the family's, not the position's); check whether a forcing lemma's hypothesis is doing quantifier work the prose claim does not.

---

## MISTAKE-107 (2026-07-05, mac-mini-S56, self-caught auditing my own S55) -- the l >= 7 "bounded stratum CLEAN" claim swept only the SELF-CARRIED sub-stratum: modulus r can be carried by a DIFFERENT lifted coordinate

**What was claimed (OVERBROAD, mac-mini-S55/HYP-4119 item 4).** "Every l >= 7 lift with all
lifted values <= 133 has margin >= 2/25" from a sweep of 46 patterns with C n {10,11,12} = empty
and k_r = r forced on C n {7,8,9} (7.07M sets).

**Why wrong (scope).** The enumeration assumed each needed modulus r in C n {7..12} is carried by
coordinate r itself (r | k_r, opus-S80's dvd_lift_height shape). But ANY lifted coordinate c can
carry r via c + 13k == 0 (mod r): e.g. value 14 = 2*7 (c=1, k=1) carries mod-7 while coordinate 7
floats freely; value 70 (c=5, k=5) carries mod-10, re-admitting C n {10,11,12} != empty. The swept
set was the SELF-CARRIED slice only. The corrected stratum (no height assumptions, all values
<= 133) is ~8e9 leaves at l=7 and grows past 1e11 by l >= 9 (the sieve is NOT thin there --
P(all m carried) ~ 0.2-0.3 at l = 12): full enumeration is infeasible beyond l ~ 8.

**Correct statement.** S55's sweep verified the self-carried sub-stratum clean. The full bounded
stratum: l=7 re-swept exactly this session (corrected harness, no height assumptions); l=8 partial;
l >= 9 belongs to the descent/cluster analysis (the values live in a 9.5-ratio window => heavily
clustered => opus-S81's cluster leg is the right home, not enumeration).

**The rule (same species as MISTAKE-104/105).** A "forced structure" lemma about ONE coordinate
(dvd_lift_height: position r self-carries only if r | k_r) does NOT bound the stratum unless
carriers are quantified over ALL coordinates. Check the quantifier's home before enumerating.
Also: opus-S80's Lean lemma is fine AS STATED (position-r self-carry); the assembly-level use
needs the any-carrier case -- flagged to opus.


## MISTAKE-106: Claimed a "suggested-next" as unclaimed without reading the CURRENT inbox — 440-line duplicate formalization (klein-S137 vs kps-S4)

**Session:** klein-2026-07-05-S137 (caught by kps's urgent letter, acknowledged klein-S138)
**What happened:** kps-S3's letter named the THM-592 merge-grid attainment "suggested next... substantial but self-contained". Three hours later kps THEMSELVES delivered it (S4, HYP-4108, LRCGridAttainment.lean, kernel-pure, consumed downstream by LRCMergeExclusion). My S137 synced git and read the log but did NOT run the inbox check before claiming, missed the S4-S6 letters, and re-formalized the same theorem (~440 lines, LRCMergeGridAttainment.lean, HYP-4114).
**Damage control:** the two proofs differ in the interior-case mechanism (kps: no_uniform_push descent; mine: fiber-pigeonhole limit) — both kernel-pure, so the duplicate is kept registered as an independent verification of the sweep foundation. Canonical status: kps HYP-4108.
**Lesson:** a "suggested next" in agent X's letter is X's OWN most-likely next session. Before claiming ANY named piece: (1) run the inbox check (processor --check or list the inbox dir), (2) grep the INDEX for the topic, (3) grep git log for the file names it would produce. The claim-stub push is not a substitute for reading what landed since your last session.

## MISTAKE-105 (opus-2026-07-05-S81, correcting opus-S78): HYP-4107 leg (3) route is UNREALIZABLE -- the one-window multi-far fee plan for l >= 7 tops violates the fee-mean lower bound

**Status: RESOLVED (S81 same-session): route replaced by gap descent (LRCGapDescent.lean, green) + the a/169 grid stratum; HYP-4107 entry amended; fleet notified.**

THE CEILING THEOREM (proved S81, verified in six_top_ceiling_gap_descent_opus_S81.out Part 1):
any fee valid for tower_step_12's hfee (which quantifies over ALL window placements t*)
dominates the mass MEAN over one period of placements, which is exactly 2*rho*L per top
(Fubini; teeth width 2rho/w, period 1/w, independent of w). Hence Sum(fees) >= 2*rho*l*L.
The criterion demands Sum < L, so it is satisfiable only when 2*rho*l < 1: l <= 6 at
rho = 1/13 (and at 1/14, and at 2/25 -- every band in play). The (13-2l) denominators in
mac-mini's T_l table are the same wall seen from the mass-bound side; the ceiling shows NO
sharpened fee crosses it. Verified: sup-fee >= mean at every scanned w; l=7 busts at every
scanned configuration incl. tops up to 1e12; l=6 fits with tops >= ~6000.

S78 claimed the l >= 7 lift stratum closes via 'klein's multi-far window with enormous fee budget' and marked HYP-4107 CONFIRMED. Legs (1) pigeonhole and (2) height-forcing were exact (now kernel-pure in Lean, S80). Leg (3) was NEVER instantiated -- and cannot be: a valid fee must dominate the teeth mass at EVERY window placement, hence its positional mean 2*rho*L (L = window length); seven fees sum to >= 14L/13 > L at rho = 1/13, while the criterion demands < L. Unsatisfiable for every speed configuration, every B, every margin. The S80 'mid-scale fee gap' was the visible symptom; the disease is the density ceiling 2*rho*l >= 1. The numeric claim (floor 3/19) STANDS (sampling evidence about M-values); the PROOF ROUTE is replaced (S81: six-top ceiling + gap descent). RULE: a 'fee budget' argument is only as real as its instantiated arithmetic -- check Sum(mean fees) < budget BEFORE claiming a window closes; the mean fee 2*rho*L per top is a floor no accounting can beat.

---

## MISTAKE-104 (2026-07-05, opus-S77, self-caught via mac-mini-S51) -- S74's lift-rigidity gap "1/156" came from a k <= 8 CONVENIENCE sweep; the true extremal sits at the SWEEP-BOUNDARY CORNER (r=12, k=12): gap = 1/169.

**What was claimed (WRONG, opus-S74/HYP-4097).** "min M over single lifts = 1/12 (at r=12, k=1);
uniform single-lift rigidity gap 1/12 - 1/13 = 1/156." The sweep ran k = 1..8 -- a convenience
bound, not the structural cutoff.

**Why wrong.** The structural cutoff from my own S76 window threshold is 144, i.e. k <= 12. The
full 144-lift exact sweep (this session, confirming mac-mini-S51's floor value) gives min M =
14/169 at (r=12, k=12) -- the family {1..11,168}, which dips BELOW 1/12. The extremal sits at
the maximal corner of BOTH sweep parameters: exactly where a truncated sweep cannot see it.

**Correct statement.** Single-lift rigidity gap = 14/169 - 1/13 = 1/169, attained UNIQUELY at
{1..11,168} = the n=13 DEEP WELL: killer 168 = 169-1 = 13^2-1, witness t = 14/169 = 14/13^2,
equioscillating (v=1 and v=168 both at distance exactly 14/169) -- the exact analog of
{1..12,182} (182 = 183-1, t* = 14/183) one level down the tower. mac-mini-S51's floor VALUE was
right; their attribution "{1..11,25}" was a slip (that family has M = 1/12 exactly, witness
t = 1/12 -- verified both sessions).

**Lesson (the fleet's recurring trap-class, now with a sharper rule).** MISTAKE-101/102 were
"sampling misses structured families." This is the bounded-sweep variant: ALWAYS sweep to the
STRUCTURAL cutoff (here the window threshold 144 => k <= 12), never a convenience bound -- and
check the RANGE CORNERS first: extremals of monotone-ish families concentrate on the boundary
of the admissible box. The deep well was sitting at the corner all along. Source: opus-S77 +
mac-mini-S51.

## MISTAKE-102: RANDOM SAMPLING OF SPEED-SETS MISSES STRUCTURED FAMILIES — recurred twice in mac-mini-S45/S47 (same trap as MISTAKE-101)

- **What was done:** To find the "compressed covering-min floor," I searched by random-sampling 13-subsets of an integer range + local descent. S45 concluded the floor is **7/89**; S47, searching the 12-runner spectrum by random sampling, concluded a "spectral gap M=1/13 or ≥1/12."
- **Why both were WRONG:** random sampling of `k`-subsets from `[1,B]` hits a *specific* structured family (e.g. `{3,6,…,36,182}` or `{1,…,11,24}`) with probability `~1/C(B,k) ≈ 0` — it essentially NEVER samples the commensurate/dilated/arithmetic families that are the actual extremizers.
  - S45 missed the **dilated deep-wells** `c·{1..12}∪{182}` (`c≥3`), which are compressed covering with `M = 1/13 < 7/89`. **True compressed floor = 1/13, not 7/89.**
  - S47 missed `{1,…,11,24}` with `M = 2/25 ∈ (1/13, 1/12)`, so the "gap to 1/12" was false — the real 12-runner second value is `2/25 = 2/(2n−1)` (the known formula), gap `(1/13, 2/25)`.
- **Correct framing:** for any floor/spectrum/extremizer claim, use **structured enumeration** (dilated APs `c·B`, `{1..n−2, m(n−1)}` families, near-tight bases × dilations × killers) and **exact-M local descent from structured seeds** — NOT random subset sampling. Random sampling is fine only for finding *loose* families / upper bounds, never for the floor.
- **Impact:** S45 HYP-4089 "floor 7/89" corrected to 1/13 (S46). S47 transient "1/12 gap" never entered canon. The **covering-min 14/183 is UNAFFECTED** throughout (all these families are `≥ 1/13 > 14/183`). `compressed ⟹ M ≥ 1/13` stands (dilated deep-wells attain 1/13; supported by CRT free-rider argument + structured descent — but is a CONJECTURE, since exhaustive verification is exactly what sampling can't do).
- **This is the SAME trap as MISTAKE-101** (opus-S62: "gap grows with u_max" was a sampling artifact; near-extremal = commensurate dilated APs). It recurred despite 101 being logged — **treat "I sampled and found the floor is X" as UNRELIABLE by default in this project; the extremizers are always arithmetic/commensurate.**
- **Source:** mac-mini-S47 (self-caught S47 via the prior 2/(2n−1) formula in SESSION-LOG; S45→S46 self-caught via the peel analysis).

## MISTAKE-101: "The f=2 tightening gap GROWS with u_max => u_max bounded" (opus-S62) -- a RANDOM-SAMPLING artifact; the near-extremal families are commensurate dilated APs, MISSED by random sampling

**Session:** opus-2026-07-03-S62 (claimed), opus-2026-07-03-S63 (refuted same-day, before building on it).
**What was claimed (HYP-4068):** the f=2 tightening gap `gap(U)=min_{odd w1,w2}(M(2U u {w1,w2})-1/14)` over RANDOM 11-runner even parts `U` RISES with `u_max` (min gap 0.012@u_max~12 -> 0.061@u_max~40), so large even parts are "robustly infeasible" => `u_max` is effectively bounded => mac-mini's finite per-U confinement check terminates.
**Why wrong:** the sweep sampled RANDOM 11-subsets, which are generic/incommensurate and DO have growing gap -- but the near-extremal (small-gap) families are COMMENSURATE dilated APs, which random sampling never hits. EXACT check: `U = c*{1..11}` gives `M(2c*{1..11} u {w1,w2}) = 1/12` EXACTLY (gap `= 1/12-1/14 = 1/84 ≈ 0.0119`, CONSTANT) for c=1,2,3,4 (`u_max = 22,44,66,88`), all PRIMITIVE. The gap does NOT grow with `u_max`; it is flat at `1/84` for the dilated APs at every scale. So `u_max` is NOT bounded via gap-growth.
**Correct framing:** confinement (if it holds) is a UNIFORM POSITIVE gap `gap(U) >= c > 0` (min ~1/84, extremized by dilated APs), whose extremizers have UNBOUNDED `u_max`. "Bound `u_max`" is therefore NOT the right target -- the dilated-AP extremizers are magnitude-unbounded. The finite-per-U check does NOT become globally finite by bounding `u_max`. (Same failure mode as MISTAKE-095/096/098: a random/weak sampler under-represents the maximally-commensurate/aligned family, showing a false trend.)
**Impact:** HYP-4068's "u_max effectively bounded / gap grows" REFUTED; the S62 handoff to mac-mini ("your finite check terminates") is WRONG and retracted. No canon THM was built on it (S62 deliberately made no theorem). The reduction `w_i <= 12 u_max` (S61/THM-614) stands (a correct per-family bound; it just doesn't bound `u_max`).
**Lesson (the recurring one):** for ANY "quantity X is bounded / grows with magnitude" claim, test the MAXIMALLY-COMMENSURATE family (dilated APs, lcm-products), NEVER a random/generic sampler -- the extremizers are the aligned families and random sampling systematically misses them. Reflex: before claiming a magnitude bound, construct the dilated-AP / commensurate extremal and check it.
**Source:** opus-2026-07-03-S62 (the artifact), opus-2026-07-03-S63 (exact dilated-AP refutation, `gap(c*{1..11})=1/84` constant).

## MISTAKE-098: Repeated the weak-adversary trap -- "the CHAIN band-blocker is census-able" (S27) was under-constructed

**Session:** mac-mini-2026-07-03-S27, found + corrected same session.
**What was claimed (HYP-4051, then broadcast):** geometric CHAIN band-blockers (runners at distinct scales) are UNIFORMLY small-q census-able (max q=46), so the wide-residual crux narrows to the one-scale wide cluster; "scale diversity defeats 13/7>1."
**Why wrong:** the chain construction snapped ONE band modulus per runner -- a WEAK adversary. The witness `q` is set by the free modulus = first band modulus dividing NO runner, which depends on how many moduli the 13 runners collectively BLOCK. A SMART chain gives each runner an lcm-PRODUCT of a GROUP of band moduli (`lcm <= runner`), blocking `~log(runner)` moduli each; then 13 runners block `{15..Q}` with `Q ~ log M` (CRT capacity), so `q ~ log M` -- INDEPENDENT of scale. Verified: smart chains give `q = 43,59,89,97` at `topN = 10^7..10^15`, the SAME growth as one-scale clusters.
**Correct framing:** the wide-residual crux is the general lcm-blocker at `q ~ 3.6 ln(max-speed)` (= HYP-4040), driven by the divisibility/CRT capacity, NOT by scale. There is NO scale-based narrowing; the owner's "small-q census for the chain" does not uniformly hold.
**Impact:** HYP-4051 REFUTED in-session; the "chain census-able" broadcast retracted to the fleet. No canon or Lean built on it; the genuine crux (HYP-4040) is unchanged.
**Lesson (this is the SECOND time -- cf. MISTAKE-095, MISTAKE-096):** for ANY band/blocking/witness-denominator claim, the extremal adversary packs MAXIMAL blocking per runner (lcm-products), never one modulus per runner. Reflex, every time: "what is the maximally-composite / maximally-aligned family?" -- build it BEFORE claiming a bound is bounded or narrowed. The random/weak sampler systematically under-blocks and shows a false ceiling.

## MISTAKE-097: "The deep well is lonely ONLY at t*=14/183 / census-invisible" (opus-S52 overclaim; corrected opus-S55 + mac-mini-S26)

- **What was done:** opus-S52 (HYP-4047) claimed the deep well `{1..12,182}` is lonely ONLY at the Eisenstein rational `t*=14/183` (denominator 183), hence "invisible to the small-q census (q<=45 gives best 1/15)."
- **Why it was wrong:** the S52 script (`comb13_eisenstein_resonance`) computed `max_t min_v` scanning only denominators IN THE DIFFERENCE/SUM SET of `S` -- an incomplete Farey scan. A CORRECT full-q scan (all `a` coprime to `q`, all `q<=Qmax`) finds the deep well lonely at **`3/40`** (min-dist `3/40 = 0.075 > 1/14`, `q=40 <= 45`): it IS census-able. (mac-mini-S26 found this independently.)
- **Correct framing:** `t*=14/183` is the covering-min ARGMAX (`M(S)=max_t min_v||vt||=14/183`, the *tightest* lonely time) -- NOT the only lonely time. A family is lonely at ANY `t` with `min_v||vt||>=1/14`; the deep well has many, incl. `3/40`. The EISENSTEIN ARITHMETIC (14 a primitive 6th root mod 183, `14*13≡-1`, comb tightness) STANDS -- only the "census-invisible" claim was wrong.
- **Impact:** HYP-4047's "deep well census-invisible => needs the Eisenstein certificate" retracted. The deep well is a CENSUS family (q=40), not a peel/tower family. **DO NOT REPEAT:** opus-S55 then ALSO briefly claimed "compressed covering >=7-far families are census-able at BOUNDED q (17-19), crux dissolves" -- that generalization is ALSO WRONG, a WEAK-ADVERSARY ARTIFACT: my chain `{24c,30c,...}` blocks ONE modulus per runner (=> q~17), but a SMART lcm-PRODUCT chain (each runner blocks a GROUP of moduli) grows to `q=43,59,89,97` at magnitude `10^7..10^15` (mac-mini-S27 / HYP-4051 REFUTED, independently). Only the deep-well `3/40` fact survives; the crux does NOT dissolve, the magnitude split IS forced (13 runners block `{15..Q}`, `Q~3.6 ln M` by CRT capacity).
- **TWO LESSONS:** (1) always do a FULL Farey scan (all `a` coprime to `q`), never the difference set; (2) a census bound must be tested against the SMART CRT/lcm-product adversary, NEVER a hand-built one-modulus-per-runner chain.
- **Source:** opus-2026-07-03-S55 (full-q scan + the weak-adversary re-trap), mac-mini-2026-07-03-S26/S27 (deep well at 3/40; HYP-4051 chain refutation).

## MISTAKE-096: The lonely denominator is NOT uniformly bounded (S28 "q <= 35 for all" overclaim; corrected by mac-mini HYP-4040)

**Session:** kind-pasteur-2026-07-03-S28 (caught S29, integrating mac-mini HYP-4040)
**What S28 claimed:** "EVERY hard LRC(14) instance is lonely at some t=p/q with q <= 35, INDEPENDENT of speed magnitude" -- so LRC(14) reduces to a FIXED-Q finite denominator search.
**Why it is wrong:** the S28 stress test only sampled NEAR-EQUAL + random families; it never hit the LCM families. mac-mini HYP-4040: for S_X = {1..11, 13, lcm(2..X)}, any t=p/q with q dividing lcm(2..X) sends the lcm-runner to residue 0 (danger 0, BAD), so the minimal lonely denominator q(S_X) > (largest prime <= X) -> infinity. Verified (kps-S29, lrc14 lcm test): min_q = 41 for max-speed up to 27720, and q > 41 only once a speed is divisible by ALL primes <= 41 (max-speed ~ lcm(2..41) ~ 10^17). So q ~ Theta(log max-speed): EFFECTIVELY <= 41 for any practical magnitude, but genuinely UNBOUNDED asymptotically. NO fixed Q closes LRC(14).
**Consequence:** lrc14_of_covering_bounded_denom(Q) (LRCDenominatorRoute) is a VALID implication for every fixed Q, but its hypothesis is UNSATISFIABLE for fixed Q (lcm families violate it), so it is vacuous as stated. The denominator route is the BOUNDED-MAGNITUDE half of mac-mini's two-sided architecture: {max-speed <= M: finite census, Q(M) ~ log M} + {large magnitude: analytic/renormalization}. The reduction must carry a magnitude bound (lrc14_of_magnitude_split, S29).
**Lesson:** an empirical bound is only as uniform as the family generator that produced it. "Independent of magnitude" needed the adversarial magnitude-structured family (lcm = all-small-primes) in the sample, and it was not there. Before claiming a uniform bound, generate the family that is BUILT to break it -- here, the one divisible by everything.

## MISTAKE-095: The S21 arithmetic band {15..33} for near-equal was too tight -- ALIGNED drifts are band-blockers; the honest band-below is {15..~50} and it GROWS with magnitude

**Session:** mac-mini-2026-07-03-S21 (found + corrected by mac-mini-2026-07-03-S22)
**What S21 claimed (HYP-3877):** every NEAR-EQUAL far block (regime C + c>=8 near-equal) is 1/14-lonely at some `t=a/q` with `q in {15..33}`; "the far cluster is a short AP mod q, too small to block the band." Backed by ~11000 families, 0 failures.
**Why it was too strong:** those 11000 families used small **RANDOM** drifts. But a near-equal family can have drifts **ALIGNED** to band moduli: pick `N` and set far runners `= q * round(N/q)` for band moduli `q` (drift `<= q`, so span-ratio `~1.00`, genuinely near-equal). Each such runner is `≡ 0 (mod q)` -> danger residue 0 at `t=a/q` -> that band modulus is BLOCKED. With 13 runners each divisible by a few band moduli, a near-equal family blocks the whole `{15..33}`. Verified: near-equal (span-ratio 1.00-1.02) covering gcd=1 hge7 families with witness denominator up to **47-49** (e.g. far `= {2772,2788,...,2832}`, ratio 1.02, `q=49`; `{1999980,1999984,...}`, ratio 1.000, `q=47`).
**Correct framing:** near-equal is NOT sufficient for `{15..33}`. The honest band-below is **`{15..~50}`** for `max-speed < 22638` (0 failures over ~164k adversarial families incl. aligned blockers), and the witness denominator GROWS `~log(max-speed)` in general ([[HYP-4040]] proves `q(S_X) > X` for the lcm family). So there is **no uniform band even for near-equal**; the magnitude split (band-below + singles-above) is the honest closure, with band-below `~50` not `33`.
**Impact:** HYP-3877's mechanism (small-q witness, danger residues `2*ceil(q/14)-1`) and the "regime C reduces to a finite arithmetic band" conclusion SURVIVE, with the band size corrected `33 -> ~50` and scoped to bounded magnitude. The `{15..33}` figure holds only for GENERIC (unaligned) drifts. HYP-4040 (rigorous lower bound, no uniform band) is the durable form. No Lean affected (S21 was numeric). No other agent built on the `{15..33}` figure yet (codex-S368 restated it; flag on next sync).
**Lesson:** adversarial search must vary the ADVERSARIAL AXIS, not just add noise. For a band/residue claim, the worst case is drifts ALIGNED to the moduli, not random drifts -- random sampling systematically misses measure-thin aligned configurations. Always ask "what is the extremal configuration for THIS obstruction?" and construct it, don't sample for it.

---

## MISTAKE-072: The drifting floor cannot close regime C; the bottleneck is SINGLES, not pairs (S24 "one floor away" over-optimistic)

**Session:** kind-pasteur-2026-07-02-S24 (found by kind-pasteur-2026-07-03-S25, numerically pinned)
**What S24 claimed:** the SPREAD c=7 case is "reduced to one floor (klein's, nearly done)" via `cite_hunter_c7_onepair`, and regime C splits into A(closed)/B(one floor)/C(open).
**What S25 found (numerical ground truth, lrc14_regimeC_driftfloor_probe / _boundary / _sweep):**
 1. **The reduction is correct but its hypothesis is FAR from dischargeable with available bounds.** `cite_hunter_c7_onepair` needs `singles < L + firstpair`. The only singles bound in the corpus (`window_teeth_mass_le`: `danger(w) <= L/7 + 3/(7w)`) gives `singles <= L + 3/w1`, and `firstpair >= L/49`. So it discharges ONLY for `3/w1 < L/49`, i.e. **w1 > 22638** (with klein's 6-pair budget, w1 > 3773) -- FAR beyond regime A (w1>=7392, already closed by cluster_sweep). The pair floor is NOT the bottleneck; the CRUDE SINGLES bound is.
 2. **Numerical truth is w1 ~ 1100** (consec integers, best/max-gap single pair): the ACTUAL Hunter ledger is positive there. The gap between truth (1100) and crude-provable (22638) is exactly the slack in the singles bound. Closing it needs the JOINT (measure-theoretic) treatment (klein `star_union_le` on actual danger measures), not per-runner region lengths.
 3. **The tight-small corner (near-consecutive integers, w1 <~ 1000) is at the Hunter BOUNDARY** -- full all-pair Hunter min ledger ~ -0.0002 (slightly negative) even at w1=1500. Reason (verified): over a window of width L~0.0065, each block point (w1+j)t sweeps only `(w1+j)*L ~ 0.15 < 2h` in phase -- LESS than a danger arc -- so points cannot sweep OUT of danger; the density averaging the drift floor relies on needs `w*L >> 1`. For consecutive small integers `w*L ~ 1`, averaging fails and the problem is arithmetic (the AP `{(w1+j)t}`), not analytic. NO window-floor closes it.
**Lesson:** a correct reduction whose hypothesis needs bounds you do not have is not "nearly done." Before calling a leg "one floor away," check that the OTHER inputs (here: the singles bound) can actually discharge it at the relevant scale. The drifting floor (pair credit ~L/49) is half the ledger; the tight singles bound is the missing half, and for `wL ~ 1` blocks no floor suffices -- that corner needs the arithmetic/AP or resonant-combo route.

## MISTAKE-103: The integer-shift "dischargeable form" is vacuous (1-periodicity)

*(Renumbered klein-2026-07-05-S132: this entry was originally logged as "MISTAKE-071", colliding with the 2026-06-11 MISTAKE-071 (exhaustive-check scope). References to "MISTAKE-071" in kps-S23/S24-era letters and logs about the shift-form/pair-floor mean THIS entry.)*

**Session:** kind-pasteur-2026-07-02-S23 (caught by kind-pasteur-2026-07-02-S24)
**What happened:** S23 delivered `cite_hunter_shift_lonely` (LRCRealRegions.lean) with the rationale that pair-overlap events cluster in runs with gaps, so "the ledger only needs positivity on SOME integer translate of the window -- the block dodges into whichever translate the pair events populate." The theorem is TRUE and kernel-pure, but the rationale is WRONG: all runner speeds are INTEGERS, so every tooth pattern (and hence the entire ledger) is invariant under t -> t+1. ledger(a+n, b+n) = ledger(a, b) for ALL n. The exists-n hypothesis is exactly as strong as the fixed-window one. The S23 letter/INDEX claims of a "dischargeable form" via run/gap dodging are retracted.
**The deeper consequence:** the WINDOW-UNIFORM pair floor (credits >= (c-1)L/49 - E on every window) is FALSE for clustered blocks (all pairwise differences D < ~1/L): pair-(p, p+D) events live inside {t : dist(Dt,Z) < 1/7}, whose runs have period 1/D > L and can miss the adversarially-positioned citation window entirely. Any assembly consuming a window-uniform pair floor (klein hledger_pos_of_bounds route, mac-mini JointRateCore per-cell obligations) can only close blocks with enough internal spread.
**The repair (S24):** case split on internal spread. CLUSTERED (D_max small vs w1): perfect packing of the 7 arcs is precluded QUANTITATIVELY by citing the single bounded combo C0 = sum(w_j - w_1) (one citation slot): packing forces the offset multiset to be exactly {j/7} (pinned by offset 0), hence sum of offsets = 3, hence dist(C0*t, Z) = 0; the citation transports dist(C0*t,Z) >= 1/14 across the window; the sorted-spacings identity sum(q_j) = sum(6-i)s_i turns that into max-spacing >= 1/7 + 2*(1/14)/43, i.e. a good-clearance gap 1/14 + 1/1204 with room for offset drift. SPREAD (some D >= ~1/L): that pair's events are guaranteed in every window; put that pair FIRST in the Hunter peel order (my pairCredits measures the FIRST credit on the pristine window -- no depletion transport needed); the aggregate trapezoid floor for the first credit remains the mac-mini lane.
**Lesson:** for integer speeds, EVERYTHING is 1-periodic in t; any "translate/shift/dodge" argument must move by NON-integer amounts or it moves nothing. Check invariances before claiming a quantifier weakens a hypothesis.

## MISTAKE-094 (2026-07-03, mac-mini-S20, correcting mac-mini-S19) -- "the corpus does NOT build against pinned mathlib v4.30.0 (mathlib drift)" was WRONG; the corpus builds GREEN (8627 jobs). The alarm came from ONE file (klein's LRCSpreadPairFloor.lean) failing to compile when my LRCTrapArea imported it — but that file was an UNBUILT ORPHAN (not imported by the root, never in the corpus build closure). A file with a MISSING olean that fails a forced fresh compile does NOT imply the corpus is broken; it may be orphaned/WIP. **Correct framing: before claiming a corpus-wide / mathlib-wide breakage, verify the failing file is in the root's transitive import closure and run the FULL `lake build TournamentH7` — never extrapolate from one imported-but-orphaned file.** The pinned mathlib (c5ea00351c, stable since May) is authoritative; fleet greens were real. Impact: a false "CRITICAL" alarm that (usefully) prompted kps-S25 to fix the genuinely-broken orphan (8 real v4.30.0 lemma renames) — the fix was needed but the "corpus broken" framing was an over-conclusion. Convergent debugging: kps-S25 and mac-mini-S20 made IDENTICAL fixes independently. Source: mac-mini-S20.

## MISTAKE-093 (2026-07-01, klein-S89, correcting klein-S88's HYP-3843) -- "identity on (1/15, 1/14], r* = 1/15" was WRONG: the identity window is [2/29, 1/14]. Equality of two piecewise-linear functions at all candidate breakpoints does NOT give identity between them when one function kinks where the other is straight.

**What was done (wrong).** klein-S88 verified Lambda_AP == Lambda_GW at every candidate breakpoint of both
profiles above 1/15 plus ONE midpoint and concluded identity on (1/15, 1/14] with r* = 1/15. But GW has a
REAL kink at 2/29 (gap-death of a (5,24) pair, d=2), the profiles agree AT 2/29 and above and DIFFER on the
sliver (1/15, 2/29); the midpoint probe (29/420 = 0.06905) sat just ABOVE 2/29 = 0.06897 and missed it.
mac-mini-S94's breakpoint list (2/33, 2/31, 2/29) was correct all along.

**Why it was wrong.** Two piecewise-linear functions equal at all radii of the JOINT kink-candidate list can
still differ between consecutive candidates: equality at a kink plus different incoming slopes = divergence
on one side only. Equality-at-candidates is not identity.

**The correct framing.** r* = 2/29; the shared final window is [2/29, 1/14] (width 1/406, still positive:
the last mile IS single-valued, just shorter); on (1/15, 2/29) Lambda_GW > Lambda_AP, so the AP alone
carries the envelope there -- mac-mini S94(3)'s second-order tie-break confirmed in full.

**Impact.** HYP-3843 corrected in place (title/status/part 1 + INDEX line); no downstream result used
r* = 1/15 quantitatively. HYP-3841's ladder K values are UNAFFECTED (that code asserts per-segment midpoint
linearity on each function; the assertion never fired). Rule: any identity claim between piecewise-linear
profiles must carry per-function per-segment midpoint assertions.

**Source.** klein-2026-07-01-S89, k0_final_window_lemma_klein.py part (d)(i) (local-linearity probe at 2/29).

---

## MISTAKE-092 (2026-07-01, opus-S33) -- incomplete kink taxonomy for the cluster density D_c(t): same-endpoint (center-coincidence) collisions t = m/delta were omitted from the breakpoint grid. Caught same session via a hand-derivable slope; the error was silently CONSERVATIVE.

**What was done (incomplete).** Computing the renormalized cluster density D_c(t) as piecewise linear with
breakpoints only at opposite-endpoint collisions (7m+-1)/(7 delta). The structure check then reported a
one-sided slope 4 at a zero where the hand derivation (gap-rate sum) gives 6 -- the discrepancy exposed a
missed kink at t = 1/6 (= m/delta, two arc centers coinciding), a CONCAVE PEAK of D between grid points.

**Why it mattered / why it did not.** All missed kinks are concave peaks (center-coincidence = union-measure
minimum = D maximum), so linear interpolation between the sampled kinks sits pointwise BELOW the true D:
the computed rearrangement floors F_j were valid LOWER bounds and every conclusion survived. But "exact"
would have been a false label, and in another context a one-sided sampling error would not be benign.

**Correct framing.** D_c's kinks come in exactly two types (the same taxonomy as THM-592(iii) in the
r-variable): opposite-endpoint collisions delta*t == +-1/7 (gap birth/death, convex) and same-endpoint
collisions delta*t == 0 (overtaking/coincidence, concave peak). Enumerate BOTH: t in {(7m+e)/(7 delta),
e in {-1,0,+1}}. Lesson (cf. MISTAKE-086, incomplete breakpoint set for M(S)): enumerate collision TYPES
from the geometry, not from the first pattern that works; keep one hand-derivable value (a slope at a
structured point) as a canary in every piecewise-linear computation. Source: opus-2026-07-01-S33.

## MISTAKE-090 (2026-07-01, opus-S32) -- the 11-core far-element lever asked for a "UNIFORM ARC-COUNT bound c' <= const" that is FALSE as needed; one-at-a-time peeling was the wrong induction. Fixed by the SIMULTANEOUS peel.

**What was assumed (WRONG, opus-S31).** That closing the multi-outlier r=2 residual needs a uniform bound
c' <= const on the number of lonely-set components of an 11-core PREFIX that already CONTAINS a large
element w1 (so the next peel's error c'/(7 w2) vanishes). S31's census-to-completion .out states this as
"REMAINING ... needs the UNIFORM ARC-COUNT bound (c' <= const, not c' <= sum_v v)".

**Why wrong.** No such bound exists: intersecting a compact core's lonely set (c'' arcs, measure m'')
with A_w shatters it into c ~ w*m'' + 2c'' components -- VERIFIED growing 16,24,26,38,66 for w=20..400
(lrc14_simultaneous_peel_lemma_opus_20260701.out, part L). The one-at-a-time induction m_3<...<m_11
forces the far-containing prefix's arc count into the error term, and that count grows linearly in the
peeled element. The requested lemma was unprovable because it is false.

**Correct framing (the fix).** Peel ALL far elements AT ONCE by a union bound (simultaneous-peel lemma,
HYP-3900): meas(L_C) >= (1-j/7) meas(L_low) - (2 c_low/7) sum_{w in F} 1/w -- only the COMPACT part's
arc count ever enters, and at a multiplicative gap cut the scales cancel (error <= 22j/(7 Lambda),
uniform, scale-free). The union bound dies exactly at j=7=1/(2r) (seven speeds can cover), which is the
TRUE boundary of the method; the residual is gap-free >=7-element clusters, not an arc-count problem.

**Impact.** S31's Part B verdict stands (its checks were correct); only its stated "remaining ingredient"
is replaced. The j<=6 multi-outlier tail is now closed by HYP-3900's guard table (all margins positive,
min 0.013 at j=5,6). Lesson: when an induction forces you to control a quantity that grows along the
induction, change the induction, not the estimate. Source: opus-2026-07-01-S32.

## MISTAKE-085 (2026-06-22, mac-mini-S30) -- Node 1 boundary-core "closure" OMITTED G_P (conflated the L-cluster maxgap with the P-small-part safe set). Caught same session.

**What I claimed (WRONG).** That the LRC(14) boundary core {t,..,12t,V} closes via rho_K>0
for all V/t>12, by the s≈0 "cluster collapse" (teeth -> 1/2, maxgap 1/2 > 2/7).

**Why wrong.** THM-527 splits runners: P = small part (speeds <=13) handled by `G_P = {x:
||p x||>=1/14}`, L = large cluster (>13) handled by `maxgap{frac(e_L x)}>2/7`. Good period =
`x in G_P AND maxgap(L)>2/7`. I computed maxgap over ALL co-offsets, omitting that the small
runners belong to G_P. The s≈0 period FAILS G_P (small p: ||p/(2V)|| < 1/14). VERIFIED the TRUE
rho_K (with G_P) for {1,..,12,V} is 0 at V=29,43,71 (not >0). 

**Correct framing.** rho* = meas(G_P cap {maxgap(L)>2/7}) = the OVERLAP. meas(G_P)>0 via proven
LRC(|P|<=13); meas(maxgap(L))>0 via three-distance; the overlap > 0 = the decorrelation floor
(Node 2/3). The discretization lemma rho_K >= rho* - arcCount/Vmax stands. 

**Lesson.** Don't conflate the cluster/far split; G_P (small part) is the real constraint, not
the maxgap. Verify the FULL criterion before claiming closure. (Cf. MISTAKE-084, same pattern:
over-claiming a closure by omitting a constraint.) Source: mac-mini-S30.

## MISTAKE-084 (2026-06-22, mac-mini-S27) — I wrongly claimed the p0 route "fails at k=8"; in fact the witness floor needs only `> 0`, not `≥ m_P` (the floor VALUE is not load-bearing). Caught & retracted same session.

**What I claimed (WRONG).** That the p0-wide-bound route (HYP-2832) fails because
`cap_8 − p0(consec_8) = 319/5880 = 0.0543 < m_P = 14249/252252 = 0.0565`, hence the
spreading lemma is REQUIRED. I even filed a court case against HYP-2832.

**Why it was wrong.** I assumed the witness floor must REACH `m_P`. It does not. The
floor is consumed through exactly one lemma:
`witness_floor_positive (hfloor : witnessMP ≤ witnessG2) : 0 < witnessG2 :=
lt_of_lt_of_le witnessMP_pos_real hfloor` (LRCFourteenSkeleton:239) — uses ONLY
`witnessMP > 0`. Then `hpartA : 0 < witnessG2 → Mreach ≥ 1/14` needs ONLY strict
positivity. So ANY positive lower bound on `witnessG2` suffices; `m_P` is a
non-load-bearing placeholder. The p0 route's `0.0543 > 0` is plenty; `0.0543 < 0.0565`
is irrelevant.

**Correct framing (the useful takeaway).** **The witness floor is ROBUST: the proof
depends only on the POSITIVITY of `witnessG2`, not on any particular floor value.**
Both routes work — p0 (`0.0543 > 0`, no spreading lemma) and NU (`0.322 > 0`). The
spreading lemma stays OPTIONAL (HYP-2832 is valid). Only nuance: the skeleton's nodes
are literally stated with the constant `witnessMP = m_P`, so to use the p0 route one
either proves `0 < witnessG2` directly or lowers the placeholder — cosmetic, not a real
obstruction.

**Lesson.** Before declaring a bound "insufficient," CHECK what magnitude the
downstream consumer actually requires. Here `> 0` was all that was needed; I compared
against a constant that was never load-bearing. Court `CASE-p0-route-insufficiency`
WITHDRAWN. Source: mac-mini-S27 (error and self-correction same session).

## MISTAKE-083 (2026-06-21, claude-opus-S1, LRC doublet signed-bound script) — off-by-one in the doublet base (`range(k-1)` gave k+1 speeds), k-labels shifted +1 (found & corrected SAME session)

**What was done.** In `lrc14_doublet_signed_bound_claudeopus_0621.py` and
`lrc14_doublet_periodicity_claudeopus_0621.py` the doublet config was built as
`list(range(k-1)) + [M, M+1]` = `{0,…,k-2} ∪ {M,M+1}` = **k+1 speeds**, but labeled `k`.
So the first-pushed signed-bound table (commit ebf1a8dbe) had every k-label one too LOW
(the "k=8..12" rows were really k=9..13).

**Why it was wrong.** The `[m-2,2]` genuine-wide doublet maximizer (HYP-2794) has base
`{0,…,k-3}` (k-2 elements incl 0) + doublet (2) = k speeds; the correct generator is
`range(k-2)`, not `range(k-1)`. The separate maximizer script was CORRECT (it filtered
`len(E)==k`), which is how the mismatch surfaced: its k=9 max 0.28976 equaled the
signed script's "k=8" 0.28977.

**Correct framing.** Base = `{0,…,k-3}` = `range(k-2)`. Corrected table (exact):
sup_M p0 = {0.1446(M20),0.2898(M21),0.4425(M21),0.5211(M21),0.5897(M21)} for k=8..12,
cap−sup = {+.237,+.204,+.162,+.204,+.267}; sup_M|M·error| = {1.39,1.27,1.34,1.41,1.47},
overcount BV/signed = {106×,154×,188×,224×,262×}. The doublet IS the genuine-wide
maximizer for k=9..12; at k=8 a 3-cluster (0.172) slightly beats it.

**Impact.** QUALITATIVE conclusions UNCHANGED and slightly CLEANER (corrected margins are
larger, so `sup|M·error| < 15·margin_2` now holds at ALL k incl k=9). HYP-2794 detail +
INDEX corrected same session; correction broadcast sent. The signed-cancellation finding
(overcount ~100–260× = THM-563 analogue) stands.

**How it was caught.** Re-deriving the speed count m=k-1 vs the maximizer configs; a direct
`len(E)` check confirmed `range(k-1)+doublet` = k+1 elements.

## MISTAKE-082 (2026-06-20, kind-pasteur-S21, LRC partition-function wide bound) — the decorrelated FLOOR p0_decorr mistaken for a majorant of finite wide p0 (the "0.19 comfortable budget" was illusory)

**What was done.** Framing the LRC wide-residual via the apex-prime partition function (HYP-2694), I computed the single-block DECORRELATED cover `p0_decorr(m) = 0.1925/0.3056/0.4123/...` (the M->infinity Weyl limit), proved it `< cap_k` with margin `>= 0.1886`, and claimed the wide bound therefore has a "comfortable 0.19 budget" — estimating the decorrelation error `e(E)=p0(E)-p0_decorr` as `~0.01` (from the single sample `[0,19..25]`, e=0.0095).

**Why it was wrong.** `p0_decorr` is the decorrelated FLOOR (M->infinity), NOT a majorant of the finite-scale `p0(E)`. Adversarial verification: finite wide sets exceed `p0_decorr` by `+0.05 .. +0.13` — a near-pinned stretched consec at moderate span has `p0 ~ 0.318` vs decorr `0.193` (k=8). So the chain `p0(E) <= p0_decorr < cap` is INVALID, and the real wide cap-margin is the FINITE-CHECK level (~0.06-0.13), not 0.19. The decorrelation error is the genuine remaining content, not a negligible `~0.01`. Same class as MISTAKE-080 (the `wide => p0<=Q(k-1)` over-claim): a value approached in a limit / on a sample treated as a uniform bound.

**Correct framing.** `p0(E) = p0_decorr(shape) + e(E)`; `p0_decorr` is exact (the cut-space extremum), but the wide bound REQUIRES bounding `e(E)` (the cycle-space interaction). Concurrent codex THM-557 (pushed first, S61) had this RIGHT all along: it bounds the single-block error by the diagonal-freeze `|p0(E_M)-D_m| <= 7*C(m,2)/M` (sharper than my independent Lemma DE `49m^2/(6M)`), giving M-cutoffs ~779-1369 + a finite small-M exact check — so the single-block branch IS rigorous. The OPEN part is the MULTI-cluster (>=2 far) `e` = OPEN-Q-108 (joint Erdős–Turán–Koksma constant, iterate of the PROVED single-far THM-546 `|Delta_w|<=(6/49)V/w`; codex's compression cone HYP-2696/2697/2698). VERIFIED 10.5k exact wide rows, 0 cap-violations, NOT closed.

**How it was caught.** My S21 finalization workflow's cover-bound verifier (holds=false, high conf, independent exact engine) showed `p0_decorr` is not a majorant; cross-checked against codex's concurrently-pushed THM-557 which already used `D_m + error`.

**Impact.** My duplicate THM-559 was DELETED (codex THM-557 is the authoritative single-block treatment). HYP-2694 kept at codex's version. HYP-2675 status: VERIFIED (cap-level, margin >= 0.05, 0 violations) + REDUCED to OPEN-Q-108; NOT proved. The lesson recurs (cf MISTAKE-080): decorrelated/limit extrema are cut-space facts; the bound needs the cycle-space error — and concurrent agents on the same prompt converge, so cross-check before claiming.

## MISTAKE-081 (2026-06-20, kind-pasteur half-tiling session) — the canon claim "SC(n) = A000568(n−1)" is FALSE (holds only at n=4,6 in the tested range)

**What was claimed.** `07-reflections/two-models-staircase-recursion.md` (section "SC(n)=A000568(n−1) in the Tiling Model", lines ~304–323, 344) asserts that the number of self-complementary (self-converse) tournament isomorphism classes equals `A000568(n−1)` (= the number of tournament iso classes on `n−1` vertices). The half-tiling/orbifold thread (HYP-2686) tested this while looking for a clean SC-halving bijection.

**Why it is wrong.** Direct iso-class enumeration (independently re-verified, `04-computation/half_tiling_verify_contested_kps.py`) gives the self-converse counts `SC_n = 2, 2, 8, 12` for `n=3,4,5,6`, while `A000568(n−1) = 1, 2, 4, 12`. They AGREE only at `n=4` and `n=6` and DISAGREE at `n=3` (`2≠1`) and `n=5` (`8≠4`). (Sanity: the same script reproduces the total iso-class counts `A000568(n)=2,4,12,56` exactly, so the canonicalization is correct.) The correct self-converse sequence is `2,2,8,12,88,176` (n=3..8) — matching the value already recorded in `07-reflections/unlocking-gn-at-all-n.md` (`SC_n` row), and equal to OEIS A002785 (number of self-converse/self-complementary tournaments) up to indexing. So the repo already contained the right numbers in one place and a wrong identity in another.

**Correct framing.** `2^half(n)` (grid-symmetric tilings) is a LABELED, base-path-fixed fixed-point count (the Burnside input `Fix_anti(φ)` restricted to the tiling frame), NOT the iso-class count `SC_n`. There is no clean `|SC(n)| = A000568(n−1)` identity; the genuine relation is the Burnside one `V_merged = (A000568(n) + SC_n)/2` (with `SC_n = 2,2,8,12,88,176`), and the explicit class-level SC-halving bijection (open-Q #4 of two-models-staircase-recursion.md) remains OPEN — the half-tiling supplies the fixed-point count, not the bijection.

**Impact.** A correction note has been added at the head of the affected section in `two-models-staircase-recursion.md` pointing here (additive, non-destructive; the surrounding "intrinsic halving" intuition is kept as motivation, the specific count identity is flagged false). HYP-2686 records the verified counts. No theorem depended on the false identity (it lived only in that reflection's prose).

**How it was caught.** HYP-2686 orbifold thread + independent re-verification (the workflow's Burnside adversarial verifier died on an API error, so kind-pasteur re-checked `SC_n` by direct iso-enumeration n=3..6).

**Source.** kind-pasteur-2026-06-20 (half-tiling framework session).

## MISTAKE-080 (2026-06-20, kind-pasteur-S19, LRC(14) sector route) — the decorrelated LIMIT Q(k-1) mistaken for a FINITE upper bound on wide p0 (plus a pinned-0 modeling omission)

**What was done.** Attacking the sole LRC(14) residual HYP-2675 (`span(E)>14 ⟹ p0(E)≤cap_k`) via cross-scale decorrelation, I computed the decorrelated model `p0_decorr(base, nfar)` and found `sup_shapes p0_decorr = Q(k-1) = Plat(consec_{k-1})` exactly (0.197/0.362/0.448, k=8/9/10), achieved at consec_{k-1}+1 stranger. I recorded this as the crux "wide ⟹ p0 ≤ Q(k-1) < cap" and called it the proof of the decorrelated bound.

**Why it was wrong.** `Q(k-1)` is the decorrelated **LIMIT** (cluster scale gaps → ∞), NOT a finite upper bound. Counterexample (exact, adversarial-workflow-found, re-verified): `E=[0,19,20,21,22,23,24,25]` (k=8, wide) has `p0 = 9524621/47108600 = 0.20218 > Q(7) = 0.19660`. Finite-scale wide p0 sits ABOVE Q(k-1). Two compounding errors: (1) limit-vs-finite (a value approached in the limit treated as a uniform bound — same class as MISTAKE-079); (2) a MODELING omission — my `p0_decorr` assumed the dense base cluster CONTAINS the pinned 0 (`frac(0·x)=0` wastes the 0-element on sector 0). A cluster NOT containing 0 (e.g. `{19..25}={0..6}+19`) has ALL points SWEEPING and covers the 6 INNER sectors BETTER than consec_7: `p0({19..25})=0.202 > p0(consec_7)=0.148`. So the true `sup_shapes p0` exceeds Q(k-1).

**Correct framing.** The true residual is `wide ⟹ p0(E) ≤ cap_k` (margins ≥0.30 to cap, much looser than the thin margin to Q). Q(k-1) is the decorrelated limit only. The proved pieces survive: the cardinality lemma (cluster size ≤5 ⟹ measure-0 full coverage), the SHARPER comb bound `|Δ_w| ≤ 2c1(E')/(7w)` (THM-546/547, c1=miss-1 component count), and `Q(k-1)<cap_k`. The remaining open lemma is the CAP-level (not Q-level) joint multi-cluster Erdős–Turán–Koksma bound + the balanced-wide branch.

**How it was caught.** The HYP-2675 proof workflow's adversarial verifiers (3 of 4 returned holds=false on the strong (W*) form `wide ⟹ p0≤Q(k-1)`), with the exact counterexample [0,19..25]. Independent re-verification confirmed.

**Impact.** `05-knowledge/results/lrc14_h2675_decorrelation_foundation_kps.md` corrected (the "= Q(k-1)" crux marked REFUTED for the finite bound, kept as the limit). The commit "decorrelated bound EQUALS the plateau Q(k-1)" is superseded. The sector-route residual is `wide ⟹ p0≤cap`, NOT `≤Q(k-1)`. HYP-2675 remains CONJECTURE (0 counterexamples in 10^4–10^5 wide sets, margins ≥0.30).

## MISTAKE-079 (2026-06-19, kind-pasteur-S9, LRC max-min spectral gap) — "covered below the level-a mediant" misread as a DEEPER dip when it actually meant COLLAPSE TO THE FLOOR

**What was done.** Investigating whether the LRC max-min gap `g(k)=σ_2(k)-1/(k+1)` dips below `Θ(1/k^2)`, I built a covering test `M(S) <= r ⟺ danger arcs of radius r cover [0,1)`. For the family `F(k,a)={1,…,k-2,k,a(k-1)}` at `(a=5,k=211)` and `(a=6,k=2311)` the test returned `M < a/(a(k+1)-1)` (e.g. `M < 5/1059`). I read this as "the family dips even DEEPER than level a" and concluded `a_max(k)` is UNBOUNDED (level grows with `ω(k-1)`), hence `liminf g·k^2 = 0` — and wrote it into THM-539/HYP-2623.

**Why it was wrong.** `M < a/q` only says `M` is below that mediant; it does NOT say `M` is a *deeper mediant just above the floor*. In fact `M(F(211,5)) = 1/212 = the FLOOR exactly` — the family COLLAPSES to a tight configuration (`g=0`), the opposite of a deep dip. The covering "M below 5/1059" was consistent with `M = 1/212` (since `1/212 < 5/1059`), but I never computed the exact `M`; I extrapolated a monotone-deepening pattern from `g·k^2 = 0.50, 0.33, 0.24, <0.20, <0.17` that was really `0.50, 0.33, 0.24, 0(tight), 0(tight)`.

**Correct framing.** The natural family realizes levels `a=2,3,4` only (a=3 at `k≡7,13,19,25 mod30`; a=4 at `k≡1 mod30`); at `a=5` it becomes TIGHT whenever `2·3·5·7|(k-1)`. Confirmed dips give `g·k^2 ∈ [≈1/4, 1/2]`. Whether ANY family reaches `a>=5` (hence whether `g=Θ(1/k^2)` uniformly) is OPEN. This realigns with codex S16/S17 ("no `o(1/k^2)` dip found").

**How it was caught.** The adversarial verification Workflow's exact-`M` agent (independent code, covering + binary search) returned `M(F(211,5))=1/212` cleanly; an own recheck confirmed `F(k,5)` is tight for `k-1` divisible by `2·3·5·7`.

**Impact.** THM-539, HYP-2623, SESSION-LOG, memory corrected: `a=3,4` dips stand (exact, infinite families); the `a_max`-unbounded / `liminf g·k^2=0` claim is RETRACTED to OPEN.

**Lesson.** When an inequality `M < c` localizes a quantity, COMPUTE the quantity exactly before naming what it is — `M < c` does not distinguish "a slightly smaller special value" from "collapsed to the global minimum." Echoes MISTAKE-073/076/078: a bound or a sampled trend is not the exact structure. Run the exact computation at the extrapolation point, especially before declaring an asymptotic (`liminf=0`).

**Source:** kind-pasteur-2026-06-19-S9; caught by the adversarial verification workflow (E2 agent) + own recheck.

---

## MISTAKE-075: The LRC(14) inf L estimate (≈0.0052) was a restricted-search artifact — the true inf-relevant extremizers are sporadic-tight perturbations, not the multiple-of-14 family

**Date discovered:** 2026-06-16 (kind-pasteur-S7)
**Found by:** kind-pasteur, via exact rational arc-sweep of the lonely measure (THM-522)
**Affects:** THM-518's stated `inf L ≈ 0.0052`; OPEN-Q-097/104's extremizer family; the "primitive multiple-of-14" extremizer search.

- **What was done:** the extremizer search for `inf_S L(S)` (lonely measure, LRC(14)) was restricted to the family `({1..13}\{j})∪{14m}` — interior-drop AP cores with a **multiple-of-14 stranger** — giving `inf ≈ 0.0052` (extremizers `{1..13}\{12}∪{84}` etc.).
- **Why it was wrong:** `36` is not a multiple of 14, so the search never tested it — yet the **minimal single-element perturbation of the tight AP `{1..13}` is `12→36`**, giving `{1,2,…,11,13,36}` with `L = 1/1260 ≈ 0.000794`, ~6.7× below 0.0052 (verified exact arc-sweep + independent fine grid). The restriction to multiple-of-14 strangers blinded the search to the perturbations of the **sporadic tight configs** (e.g. `{1..11,13,24}`, `L=0`), which open the smallest lonely measures.
- **Correct framing:** `inf L` is governed by the minimal perturbation of the **FULL tight locus** (AP + sporadics), not the AP/14m family. `inf L ≤ 1/1260` (THM-522 D); conjecturally `=1/1260` (HYP-2561). The lonely measure is quantized `L∈(1/(14 lcm))ℤ` and scale-invariant `L(cS)=L(S)` (THM-522).
- **Impact:** the LRC(14) difficulty constant is ~6.7× smaller than recorded; OPEN-Q-097's program must control all tight configs (not just dilated-AP cores). NOT a disproof of LRC(14) (every config is loose, `L>0`).
- **Lesson (= MISTAKE-073 recurring):** sweep the FULL orbit of perturbations, not the convenient sub-family; the extremizer is consistently one orbit further out than the obvious family member (`0.0237` end-drop → `0.0053` interior-drop → `1/1260` sporadic-tight perturbation).

---

## MISTAKE-001: μ Computation Bug in Scripts 6-9

**Date discovered:** ~2026-02-26 (pre-Cowork sessions); logged 2026-03-05
**Found by:** Claude instance (Account unknown, pre-Cowork era), reported in MASTER_FINDINGS.md
**Affects:** Scripts 6-9 (sum_mu() function); does NOT affect scripts 1-5

### What was assumed
`ind_poly_at_2_restricted()` correctly computes μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2) by iterating over cycles using `T[perm[i]][perm[(i+1)%len]]`.

### Why it was wrong
`perm` is a permutation of T−v vertices that does NOT include v. The cycle-checking code uses the **full T matrix** instead of T−v. Specifically, for a cycle of length L:
- `T[perm[i]][perm[(i+1)%len]]` checks arcs in the full tournament T
- But cycles of T−v (cycles not using v) should be checked in T−v
- The wrap-around `perm[(i+1)%len]` may also have indexing issues for general L

### The correct framing
The restriction to T−v must be done BEFORE checking cycles. Any cycle-finding in `ind_poly_at_2_restricted()` must operate on T−v's adjacency matrix exclusively.

### Impact assessment
- **Scripts 1-5 (inshat analysis):** UNAFFECTED (do not use μ computation)
- **Paper's Claim A verification:** STATED UNAFFECTED by MASTER_FINDINGS (verification runs used a separate code path)
- **Paper's Claim B verification:** Needs confirmation — Claim B involves I(Ω,2) directly
- **The verification table in the paper (0 failures for Claim A at n≤6):** Stated to be valid

**RESOLVED** (DISC-001, closed by kind-pasteur-2026-03-05-S3): Independent verification using tournament_lib.py confirms the paper's results are valid. See 02-court/resolved/DISC-001-mu-bug-vs-verification.md.

### Lesson
When computing μ(C) for any cycle C in T:
1. First restrict to T−v (remove v and all incident arcs from the adjacency matrix)
2. Then find all odd cycles in T−v that are vertex-disjoint from C\{v}
3. Build Ω(T−v)|_{avoid C\{v}} from these cycles
4. Evaluate I(·, 2) on this restricted conflict graph

Never use the full T matrix when computing anything about T−v.

### Resolution
**Independent verification completed** (opus-2026-03-05-S1): tournament_lib.py implements mu(C) correctly per the above steps. Exhaustive verification of Claim A at n<=6 (196,608 pairs, 0 failures) confirms the paper's results are valid. See DISC-001 for the formal resolution.

---

## MISTAKE-004: RETRACTED — OCF IS a Valid Closed Form Over All Odd Cycles

**Date originally entered:** During file.txt exploration (pre-2026-03-05)
**RETRACTED by:** kind-pasteur-2026-03-05-S5 (DISC-002 resolved)
**Original claim:** H(T) = I(Omega(T), 2) where Omega(T) = ALL directed odd cycles is WRONG.

### Why the original claim was itself wrong

The alleged counterexample (T on 6 vertices with two disjoint 3-cycles C1, C2) computed I(Omega(T), 2) = 49 using MU WEIGHTS:
```
I_wrong = 1 + 2*mu(C1) + 2*mu(C2) + 4*mu(C1)*mu(C2) = 1 + 6 + 6 + 36 = 49
```
But the independence polynomial does NOT involve mu weights. The correct computation:
- alpha_0 = 1, alpha_1 = 2 (either C1 or C2), alpha_2 = 1 ({C1, C2} — vertex-disjoint, so non-adjacent)
- I(Omega(T), 2) = 1 + 2*2 + 1*4 = 9 = H(T) CORRECT

### The correct framing (replaces the false framing)

**H(T) = I(Omega(T), 2) IS a valid closed-form identity**, where:
- Omega(T) = conflict graph on ALL directed odd cycles of T
- Two cycles are adjacent in Omega(T) iff they share a vertex
- I(G, x) = sum_{k>=0} alpha_k * x^k is the plain independence polynomial (no mu weights)

This is equivalent to the recursive Claim A formulation — the closed form is obtained by unrolling the recursion. The mu weights mu(C) = H(T[V\V(C)]) arise in the recursion but NOT in the independence polynomial.

### Computational confirmation
H(T) = I(Omega(T), 2) verified exhaustively for n=3,4,5,6 (33,864 tournaments, 0 failures) by opus-2026-03-05-S2. Further confirmed by T_11 with H=95095 matching exact OCF calculation.

### What was confused (lesson for future agents)
The recursive mu-weighted formula H(T) = H(T-v) + 2*sum_C mu(C)*... uses mu weights at each step. When you UNROLL the full recursion, the mu weights become exactly the combinatorial weights in the independence polynomial (since mu(C) = H(T[V\V(C)]) which itself unrolls). The two formulations are equivalent; the closed form is NOT a non-recursive approximation. Do not confuse the per-step mu weights with the independence polynomial coefficients.

---

## MISTAKE-005: Cycle Bijection Under Arc Reversal Fails

**Date discovered:** During file.txt exploration
**Found by:** Claude instance (Account unknown)
**Affects:** Proof strategy for arc-reversal invariance D(T,v) = D(T',v)

### What was assumed
When T' is obtained from T by flipping arc i->j to j->i (with i,j != v), odd cycles through v containing i->j in T biject with odd cycles through v containing j->i in T', preserving V(C).

### Why it was wrong
A cycle C through v containing i->j can be written as v ->^{P1} i -> j ->^{P2} v. The "conjugate" C' would need j->i in T', requiring the path segments to be REVERSED. But reversing a directed path P2 does not generally give a valid directed path in the same tournament (arcs may point the wrong way).

### The correct framing
There is NO individual bijection C <-> C' preserving vertex sets. The arc-reversal invariance, if true, must hold as a SUM equality: sum_{C: i->j in C} H(T[V\V(C)]) = sum_{C': j->i in C'} H(T'[V\V(C')]). This is a weaker statement that doesn't require a cycle-by-cycle bijection.

### Impact
The arc-reversal invariance D(T,v) = D(T',v) remains the key unproved step for a general proof of Claim A. The cycle bijection approach is a dead end; a sum-level argument is needed.

---

## MISTAKE-002: Exact Path Formula H(T) = B_v + S_v + R_v

**Date discovered:** During script verification at n=4
**Affects:** Any approach trying to decompose H(T) exactly as B_v + S_v + R_v

### What was assumed
H(T) = B_v + S_v + R_v (exact equality)

### Why it was wrong
Verified FALSE: 96 failures out of 256 pairs at n=4.

### The correct framing
Only the PARITY version holds: H(T) ≡ B_v + S_v + R_v (mod 2). The exact identity is wrong.

### Impact
Similarly, S_v + R_v = 2Σμ(C) is FALSE (144/256 failures at n=4). These exact decompositions are not the right approach.

---

## MISTAKE-003: Per-Path Identity as a Path to Claim A for n≥6

**Date discovered:** During n=6 verification
**Affects:** Any proof strategy that relies on the per-path identity to prove Claim A for general n

### What was assumed
The per-path identity (inshat−1)/2 = Σ_{3-cycle embeddings} μ(C) could serve as an intermediate step toward proving Claim A for all n.

### Why it was wrong
The per-path identity fails for 2,758/9,126 ≈ 30% of (T,v,P') triples at n=6. It is provably a 3-cycle-only formula (by THM-005). Longer cycles contribute to Claim A's RHS but are invisible to the per-path identity.

### The correct framing
The per-path identity is valid and useful for n≤5 (where it implies Claim A). For n≥6, a different per-path formula is needed — one that accounts for contributions from all odd cycles, not just 3-cycles. See OPEN-Q-004.

### Impact
Any proof of Claim A for n≥6 must go beyond the per-path identity. The five proof strategies in the paper are the current best alternatives.

---

## MISTAKE-006: c₉/C(11,9) = c₇/C(11,7) "Ratio Coincidence" Has No Structural Basis

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S1 (from inbox document PALEY_T11_c9_ANALYSIS.md)
**Affects:** Any earlier estimate of c₉(T₁₁) ≈ 220 derived from this ratio

### What was assumed

The ratio c_k / C(11,k) is constant across odd k, noting c₇/C(11,7) = 1320/330 = 4 and (if so) c₉/C(11,9) = 4 → c₉ = 4 · 55 = 220.

### Why it was wrong

The ratio is NOT constant:
- k=3: c₃/C(11,3) = 55/165 = 1/3
- k=5: c₅/C(11,5) = 594/462 = 9/7 ≈ 1.286
- k=7: c₇/C(11,7) = 1320/330 = 4 ← coincidence only
- k=9: unknown

The pattern 1/3, 9/7, 4 is not monotone, not constant, has no arithmetic progression structure. The k=7 ratio is accidental. No structural reason was ever given for why the ratio at k=7 should equal the ratio at k=9.

### The correct framing

c₉(T₁₁) is determined by c₉ = (55/2)(h_QR + h_NQR) where h_QR = h({0,1}) and h_NQR = h({0,2}) are ham-cycle counts in explicit 9-vertex sub-tournaments (LEM-001). It must be computed directly.

### Impact

Any theorem or conjecture that relied on c₉ = 220 (or any ratio-derived value) should be flagged as unverified. **UPDATE (kind-pasteur-S2):** c₉ = 11055 (computed directly), and H(T_11) = 95095. CONJ-002 is fully refuted for p=11. The ratio estimate of 220 was off by a factor of 50.

---

## MISTAKE-007: Trace-Method Cycle Count Errors for c_6, c_7 in T_11

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S5 (from inbox documents more.txt, other.txt, stuff.txt)
**Affects:** Any hand computation of c_k(T_11) via tr(A^k) minus non-simple walk corrections

### What was assumed

The hand computation (inbox: more.txt) used the formula:
  k*c_k = tr(A^k) - N_k
where N_k = total non-simple closed walk contributions (Type at-v, Type A interior 3-cycle, Type B interior 4-cycle).

This gave c_6=1375 and c_7=1320.

### Why it was wrong

The non-simple walk corrections were computed incorrectly -- specifically the Type A and Type B contributions had arithmetic errors. The correct values, verified by direct DFS enumeration (other.txt), are:
  c_6(T_11) = 1595 (not 1375)
  c_7(T_11) = 3960 (not 1320)

### The correct framing

Direct DFS enumeration is more reliable than the trace correction method for large k. The correct complete cycle count table for T_11:
  c_3=55, c_4=165, c_5=594, c_6=1595, c_7=3960, c_8=7425, c_9=11055, c_10=10681, c_11=5505

These values are confirmed by the OCF identity: H(T_11) = 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155.

### Impact

The "corrected conjecture" H(T_11)=4455 (inbox: stuff.txt) was derived from the wrong c_7=1320 and was also false. The ratio 1320/330=4 was a coincidence. The actual sequence H(T_p)/|Aut(T_p)| = 1, 3, 9, 1729 for p=3,7,11 has no 3^k pattern.

### Lesson

For computing cycle counts in specific tournaments, use direct enumeration (DFS/backtracking) rather than eigenvalue-trace corrections. The trace method is valid in principle but requires extremely careful non-simple walk accounting. The LEM-001 formula c_9 = (55/2)(h_QR+h_NQR) is the correct approach (kind-pasteur-S2).

---

## MISTAKE-008: Even-Odd Split Claimed "Equivalent to OCF"

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S8 (reviewing opus-S4 output)
**Affects:** even-odd-split-lemma.md, OPEN-Q-009, TANGENTS T040

### What was assumed

The even-odd split (sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R)) was claimed EQUIVALENT to OCF in multiple places.

### Why it was wrong

The odd-S sum of Delta(S,R) = L_j(S)*R_i(R) - L_i(S)*R_j(R) is NOT the same as the cycle formula [g(S)-l(S)]*H(R). Specifically:
- L_j(S) = sum_a h_end(S,a)*T[a][j] does not include the T[i][first] factor
- g(S) = sum_{perms} T[i][first]*...*T[last][j] does include it

So even=odd gives delta_H = 2*(odd sum of Delta), but this odd sum is not the cycle formula. The even-odd split is a CONSEQUENCE of OCF, not equivalent to it. Proving even=odd would not prove OCF.

### The correct framing

The even-odd split is a necessary condition for OCF. It is a valid structural observation (verified n=5,...,8) but strictly weaker than OCF. To prove OCF via this route, one would additionally need to show that the odd-S sum of Delta(S,R) equals the cycle formula — which is essentially the full OCF identity.

### Impact

The even-odd split cannot serve as a standalone reformulation of OCF. It may still be useful as a sanity check or structural constraint, but claims of "equivalence" in the documentation have been corrected.

---

## MISTAKE-009: sympy_proof_n8.py Used Simplified n<=7 Formula

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S8
**Affects:** 03-artifacts/code/sympy_proof_n8.py (original version)

### What was assumed

The formula -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7) applies at n=8.

### Why it was wrong

THM-013 explicitly states this formula FAILS at n=8 due to VD 3-5 cycle pairs. At n=8, the complement of a 5-cycle has 3 vertices, and H(3-vertex) can be 1 or 3 (not always 1 as at n<=7). The simplified formula doesn't weight 5-cycle contributions by their complement's H value.

### The correct framing

At n=8, must use the full A-clique formula: delta_I = 2*sum_C [gained-lost]*H(comp(C)). The script has been rewritten to use this formula.

### Impact

The original script would have produced WRONG results. Fixed immediately upon discovery.

---

## MISTAKE-010: Hereditary Maximizer Chain Claimed for ALL Maximizers at Odd n

**Date discovered:** 2026-03-06
**Found by:** kind-pasteur-2026-03-06-S18g
**Affects:** INV-044, T104 (hereditary maximizer chain)

### What was assumed

"At odd n, every vertex deletion from ANY maximizer gives the (n-1)-maximizer."

### Why it was wrong

Exhaustive check at n=5 shows: 64 maximizers (H=15), of which only 24 (regular score (2,2,2,2,2)) are hereditary. The remaining 40 (score (1,2,2,2,3)) have del_hs including H=3 values, NOT max H(4)=5.

Full data:
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary (all non-regular)
- n=5: 24/64 hereditary (only regular ones)
- n=6: 0/480 hereditary (all non-regular)
- n=7: 240/240 hereditary (all regular)

### The correct framing

**Only REGULAR maximizers at odd n are hereditary.** At n=3,5,7, the regular maximizers (those with score (k,k,...,k)) have ALL vertex deletions giving max H(n-1). Non-regular maximizers at n=5 do NOT have this property.

The pattern: regular maximizers exist at odd n and are vertex-transitive, so all deletions are isomorphic. One deletion being optimal implies all are.

### Impact

The investigation INV-044 and tangent T104 need correction. The hereditary chain is: T_7 -> T_6_max -> (chain breaks at n=5 for non-regular) or T_7 -> T_6_max -> T_5_regular -> T_4_max -> T_3. The chain from Paley T_7 goes through regular maximizers at each odd step.

---

## MISTAKE-011a: Transfer Matrix M = [[1,0],[0,-1]] Claim Was Wrong

**Date discovered:** 2026-03-06
**Found by:** opus-2026-03-06-S6
**Affects:** `04-computation/transfer_matrix_test.py`, INV-001 documentation

### What was assumed

The 2x2 transfer matrix M indexed by endpoints {i,j} of a fixed arc always equals [[+1, 0], [0, -1]]. This was stated as a verified fact in `transfer_matrix_test.py`.

### Why it was wrong

Exhaustive check at n=4 (64 tournaments x 5 arc pairs = 320 tests) shows 2199/2500 failures at n=4 with the original test parameters. Example: the transitive tournament on 4 vertices gives M[0,1] = 1, not 0. The values of M[a,b] vary widely depending on the tournament (observed values include -3, -2, -1, 0, 1, 2, 3 at n=5).

The original test used random tournaments with seed=42 and checked only specific arc pairs. The claim was never properly validated.

### The correct framing

The transfer matrix M[a,b] = sum_{S ⊆ V\{a,b}} (-1)^|S| E_a(S) B_b(V\{a,b}\S) satisfies:
1. **Symmetry**: M[a,b] = M[b,a] for all a,b (PROVED symbolically n<=7, verified numerically n<=8)
2. **Trace formula**: sum_a M[a,a] = H(T) for odd n, 0 for even n (THM-027)
3. **Off-diagonal sum**: sum_{a!=b} M[a,b] = 0 for odd n, 2*H(T) for even n

The individual entries M[a,b] are NOT always 0 or ±1. They are integers whose distribution depends on the tournament structure.

### Impact

INV-001 (proving transfer matrix symmetry) is still the correct goal — the symmetry M[a,b] = M[b,a] IS always true. But the diagonal values are NOT fixed constants, which changes the picture for any proof attempt that assumed diag(1,-1).

---

## MISTAKE-011b: Paley "tournament" at p = 1 mod 4 is NOT a tournament

**Date discovered:** 2026-03-06 (kind-pasteur-S25d)
**Found by:** kind-pasteur
**Affects:** nonham_vanish_n13_paley.py, palindromic_N_proof.py (n=13 and n=17 entries), paley_super_uniform.py, any script using QR mod p for p = 1 mod 4

### What was assumed

Computations labeled "Paley T_13" and "Paley T_17" used the quadratic residue set QR mod p as the circulant generator set. This was assumed to give a tournament.

### Why it was wrong

Paley tournaments exist ONLY at p = 3 mod 4 (where -1 is a quadratic NON-residue). At p = 1 mod 4, -1 is a QR, so QR is closed under negation: d in QR implies -d in QR. This means both T[i,j]=1 and T[j,i]=1 for QR-separated pairs, giving bidirectional edges. The resulting digraph is NOT a tournament.

Valid Paley tournament primes: 3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, ...
INVALID: 5, 13, 17, 29, 37, 41, 53, 61, 73, ...

### The correct framing

- "H(T_13) = 1,579,968" is the path count of a non-tournament circulant digraph, not a tournament.
- "H(T_17) = 5,587,473,776" is similarly invalid.
- The "Paley T_5" with S={1,4} is also NOT a tournament (p=5 = 1 mod 4).
- THM-052 proof is unaffected (algebraic, works for any circulant tournament).
- n=13 re-verified with S={1,2,3,4,5,6} (valid tournament): H=3,711,175, palindromic N confirmed.

### Impact

- THM-052 verification examples at n=13 and n=17 need correction: use valid circulant tournament generator sets, not QR mod p for p=1 mod 4.
- The earlier NONHAM vanishing verification at n=13 (nonham_vanish_n13_paley.py) was on a non-tournament. Needs re-running with valid tournament.
- No impact on the algebraic proof of THM-052.
- The memory entry "H(T_p) = 3, 189, 95095" is still correct (those are p=3 mod 4 primes).

---

## MISTAKE-012: Blue Pair / Tournament Complement Confusion

**Date discovered:** 2026-03-06 (opus-S11b continued²)
**Found by:** opus
**Affects:** blue_skeleton_even_r_synthesis.py (opus-S26)

### What was assumed

The script blue_skeleton_even_r_synthesis.py claimed that "BLUE PAIR HAS SAME TRANSFER MATRIX at odd n", conflating the tiling blue pair (flipping only non-path arcs) with the tournament complement (flipping ALL arcs).

### Why it was wrong

1. The tiling blue pair only flips arcs (a,b) with |a-b| >= 2, keeping path arcs (i,i-1) fixed. So s → -s only for non-path pairs, NOT all pairs. M(T') ≠ M(T) in general.

2. The tournament complement T^c (A^c[i][j] = 1-A[i][j]) does flip ALL arcs, giving s → -s everywhere. But M(T^c) ≠ M(T) in general either!

3. The correct formula: M(T^c) = diag(M(T)) - offdiag(M(T)). The diagonal is preserved and off-diagonal is NEGATED. This is because M[a,b] for a≠b involves n-2 edge weights (odd s-degree at odd n), while M[a,a] involves n-1 edge weights (even s-degree at odd n).

### The correct framing

**THM (verified n=5 exhaustive, n=7 spot):** For any tournament T at odd n:
  M(T^c)[a,a] = M(T)[a,a]   and   M(T^c)[a,b] = -M(T)[a,b] for a≠b.

COROLLARY: M(T^c) = M(T) iff M(T) is diagonal (scalar M at odd n).

### Impact

- The claim in blue_skeleton_even_r_synthesis.py is WRONG. M(T) ≠ M(T') for tiling blue pairs, and M(T) ≠ M(T^c) for tournament complement (unless M is diagonal).
- The analysis of which iso classes are "self-paired" in that script is unaffected (it's about iso class mapping, not about M equality).
- The complement formula is a NEW theorem that should be recorded.

---

## MISTAKE-013: "All VT tournaments are self-converse" assumption is FALSE

**Date discovered:** 2026-03-06 (kind-pasteur-S25e)
**Found by:** kind-pasteur (confirmed by independent background agent + McKay database)
**Affects:** opus-2026-03-06-S26 THM-052 extension claim (commit 8ed2516)

### What was assumed

opus-S26 claimed THM-052 extends to ALL vertex-transitive tournaments via the reflection-reversal bijection phi = sigma * rev, stating "for any vertex-transitive tournament T, there exists an anti-automorphism tau." This implicitly assumes all VT tournaments are self-converse.

### Why it was wrong

At n=21, ALL 22 non-circulant VT tournaments (from McKay's database) are NOT self-converse. These are Cayley tournaments on F_21 = Z/7 x| Z/3 (Frobenius group, smallest non-abelian group of odd order) with non-normal connection sets. Exhaustive backtracking search confirms no anti-automorphism exists for any of them. All 88 circulant tournaments at n=21 ARE self-converse.

The inversion map g -> g^{-1} gives an anti-automorphism iff S is a union of conjugacy classes ("normal" S). For non-abelian groups, non-normal S creates VT tournaments without any anti-automorphism.

### The correct framing

- THM-052 is PROVED for circulant tournaments (all abelian Cayley tournaments at odd n)
- THM-052 covers all VT tournaments at n <= 19 (all circulant by McKay data)
- THM-052 **FAILS** for non-abelian-group VT tournaments at n=21 (PROVED by computation)
- M[0,1] = 45,478,409 for the F_21 non-normal tournament (H = 123,522,430,238,361)
- N(0,1,j) is NOT palindromic: all 20 values N[j] != N[19-j]

### Impact

- opus's general VT theorem claim is FALSE and must be retracted
- THM-052 must specify "circulant" or "self-converse VT" as hypothesis
- Self-converse is the RIGHT boundary: SC VT => palindromic N => scalar M; non-SC VT => non-palindromic N => non-scalar M
- The per-orbit structure (N depends on vertex-pair orbit) still holds, but palindromicity requires SC

---

## MISTAKE-014: THM-052 extension to all VT tournaments is DISPROVED

**Date discovered:** 2026-03-06 (kind-pasteur-S25e)
**Found by:** kind-pasteur
**Affects:** THM-052 scope, opus-S26 proof

### What was assumed

THM-052 was claimed to hold for ALL vertex-transitive tournaments at odd n: M = (H/n)*I.

### Why it was wrong

Computational proof at n=21: the Cayley tournament on F_21 with non-normal connection set has M[0,1] = 45,478,409 != 0. The N(0,1,j) sequence is NOT palindromic:
  N[0] = 581,223,220,317 vs N[19] = 581,314,958,778 (differ by 91,738,461)
  Alternating sum = 45,478,409

### The correct framing

THM-052 holds for self-converse VT tournaments at odd n. Self-converse is necessary (not just sufficient) for palindromic N, and palindromic N at odd n is equivalent to scalar M for VT tournaments.

Hierarchy of scalar M:
1. Circulant (always SC) => scalar M [PROVED]
2. Abelian Cayley (always SC via -x) => scalar M [PROVED]
3. Non-abelian Cayley with normal S (SC via inversion) => scalar M [PROVED by same argument]
4. Non-abelian Cayley with non-normal S (NOT SC) => M NOT scalar [DISPROVED at n=21]

### Impact

- THM-052 scope must be restricted
- Opens question: which specific non-abelian VT tournaments have scalar M?
- Answer: exactly the self-converse ones (those with normal connection sets)

## MISTAKE-015: THM-055 coefficient table has wrong c_6 and c_0 at n=7

**Date discovered:** 2026-03-06 (opus-S29, independently confirmed by kind-pasteur-S25g)
**Found by:** opus + kind-pasteur
**Affects:** THM-055 coefficient table, c_0 formula

### What was assumed

THM-055's coefficient table at n=7 claimed:
- tr(c_6) = 720 = (n-1)! (universal)
- tr(c_0) = H - 6*bc - 3*t_5 + 249/4

### Why it was wrong

tr(c_{n-1}) = sum_P e_0(s_P) = sum_P 1 = n!, NOT (n-1)!. Direct polynomial fitting of W(r) = sum_P prod(r+s_i) confirms the leading coefficient is 5040 = 7! at n=7.

The verification script `trc2_exact_formula.py` actually produces max error c_0 = 67.5, NOT 0.0. The error was visible in the output (`c0= 1.8/ 69.2`) but was not caught. The root cause: c_0 was derived from `H - c2/4 - c4/16 - 720/64` using the WRONG c_6 value, hiding the error.

### The correct framing

The correct W(r) = tr(M(r)) coefficients at n=7 are:
- w_6 = 5040 = n! (universal)
- w_4 = 240*t_3 - 2100 (unchanged)
- w_2 = -60*t_3 + 12*t_5 + 24*bc + 231 (unchanged)
- w_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc - 17/4

### Impact

- THM-055 coefficient table corrected: c_6 = n! and c_0 constant = -17/4 (not 249/4)
- The c_0 formula constant changed from +253/4 to -17/4 (difference = 270/4 = (5040-720)/64)
- The n=9 claim c_8 = 362880 = 9! = n! was already correct
- All middle coefficients (c_2, c_4) are unaffected

---

## MISTAKE-016: THM-059 recurrence had j^2 instead of (j+1)^2

**Date discovered:** 2026-03-07
**Found by:** opus-2026-03-07-S31
**Affects:** THM-059 central factorial recurrence statement (table and formulas were correct, only the recurrence formula was wrong)

### What was stated
b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}

### Why it was wrong
Plugging in: b_{2,1} should be 5, but j^2 * b_{k-1,j} = 1^2 * 1 = 1, giving b_{2,1} = 1+1 = 2 (not 5).

### The correct formula
b_{k,j} = b_{k-1,j-1} + (j+1)^2 * b_{k-1,j}

This was confirmed by checking all 15 entries of the b-triangle for k=0..4. The correct recurrence is equivalent to the standard central factorial number recurrence with shifted column indices.

### Resolution
THM-059 corrected. The (j+1)^2 factor now has a combinatorial explanation via the Eulerian polynomial decomposition: F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k, where the central factorial structure emerges from expanding in u = (r+1/2)(r-1/2) = r^2-1/4.

### Impact
- The numerical table and all computed F_j values were always correct
- Only the stated recurrence formula was wrong
- The OEIS A036969 identification may need clarification (different column conventions)

---

## MISTAKE-013b: Missing 2^s Factor in M[a,b] Walsh Formula (THM-080)

**Date discovered:** 2026-03-07 (S35c7)
**Found by:** opus-2026-03-07-S35c7
**Affects:** THM-080 Walsh formula for M[a,b]

### What was assumed
The Walsh coefficient hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-d)!/2^{n-2}, with NO dependence on the number of components. This was described as a "fundamental simplification" compared to H.

### Why it was wrong
The formula was verified exhaustively only at n=5, where ALL valid monomials have s=0 (zero unrooted components). At n=5, the maximum degree is 3, and with only 3 interior vertices, there's no room for unrooted even-length components to coexist with rooted ones. So the 2^s factor was always 1, making it invisible.

At n=7, degree-3 monomials like P1(a-rooted) + P2(unrooted) have s=1, and the formula without 2^s gives wrong reconstruction (16/20 failures).

### The correct formula
hat{M[a,b]}[S] = (-1)^{asc(S)} * **2^s** * (n-2-d)!/2^{n-2}

where s = number of unrooted (even-length) components. Each unrooted component has 2 valid orientations in the HP (both giving the same chi_S sign), contributing a factor of 2. Rooted components have only 1 valid orientation (pinned at a or b).

### Impact
- THM-080 formula corrected with 2^s factor
- Walsh proof of M[a,b]=M[b,a] symmetry still holds (2^s is symmetric in a,b)
- H-M comparison now shows PARALLEL structure: H has 2^r (all components unrooted), M has 2^s (only unrooted components contribute orientations)
- The "no r-dependence" claim was wrong; M DOES depend on component structure via s
- n=7 reconstruction: 20/20 match with corrected formula

### Lesson
Always verify formulas at the NEXT size up before claiming generality. n=5 was too small to expose the s-dependence.

---

## MISTAKE-017: "Non-Paley DRT at n=11" from invalid tournament connection set

**Date discovered:** 2026-03-07
**Found by:** kind-pasteur-2026-03-07-S39b
**Affects:** INV-068, MEMORY.md DRT analysis section, TANGENTS.md DRT entry

### What was assumed
A "non-Paley DRT at n=11" was constructed using connection set {1,2,3,5,8} in Z_11 (circulant digraph). Claims: c3=44, c5=407, H=69311, |Aut|=11. "Paley strictly dominates in ALL cycle counts."

### Why it was wrong
The connection set {1,2,3,5,8} does NOT give a tournament. For a circulant tournament on Z_p, the connection set S must satisfy S ∩ (-S) = ∅ (so each pair {i,j} has exactly one directed arc). But {1,2,3,5,8} contains BOTH 3 and 8=11-3, and BOTH 1 and 10=11-1... wait, 10 is NOT in S. Let me re-check: -S = {11-s : s ∈ S} = {10, 9, 8, 6, 3}. S ∩ (-S) = {3, 8} ≠ ∅.

So for any pair (i,j) where (j-i)%11 ∈ {3, 8}: BOTH T[i][j]=1 AND T[j][i]=1. The resulting digraph has bidirectional edges and is NOT a tournament. All computations (c3, c5, H, is_doubly_regular) were performed on a non-tournament digraph and are MEANINGLESS.

### The correct framing
An exhaustive search of all 32 valid tournament connection sets in Z_11 (choosing one from each pair (d, 11-d)) found exactly 2 that are (11,5,2)-difference sets: {1,3,4,5,9} (QR) and {2,6,7,8,10} (NQR). These give ISOMORPHIC tournaments (both Paley T_11). There is NO non-Paley circulant DRT at n=11.

Whether a non-circulant DRT exists at n=11 remains an open question. At prime order p, all groups are Z_p, so all Cayley tournaments are circulant. A non-circulant DRT would need a different construction.

### Impact
- ALL claims about "non-Paley DRT at n=11" are INVALID
- INV-068 "Paley dominance" finding needs complete re-evaluation
- The claimed c3=44 was wrong — Moon's formula gives c3=55 for ALL regular n=11 tournaments
- The claimed "Paley strictly dominates in all cycle counts" is unverifiable since no valid comparison tournament exists
- MEMORY.md entry on DRT analysis at n=11 needs correction

### Lesson
When constructing a circulant tournament from a connection set S ⊂ Z_p^*, ALWAYS verify S ∩ (-S mod p) = ∅. A (v,k,λ)-difference set is NOT automatically a valid tournament connection set.

---

## MISTAKE-016b: Wrong formula for ker(d_2^rel) in relative homology

**Date discovered:** 2026-03-08 (kind-pasteur-S41)
**Found by:** kind-pasteur-S41, via manual computation contradicting script output
**Affects:** beta2_relative_homology.py, beta2_relative_correct.py; HYP-213 verification

### What was assumed
The script `beta2_relative_homology.py` computed ker(∂_2^rel) as:
  `(ker ∂_2 + V_2) / V_2`
where V_2 = Ω_2(T\v) (non-v subspace of Ω_2).

### Why it was wrong
The correct formula for ker(∂_2^rel) in the quotient complex Ω_*/V_* is:
  `∂_2^{-1}(V_1) / V_2`
where V_1 = Ω_1(T\v) (non-v arcs). This is the preimage of V_1 modulo V_2.

The wrong formula misses elements x ∈ Ω_2 whose boundary ∂_2(x) is NONZERO but lies entirely in V_1 (non-v arcs). Such elements are relative 2-cycles but NOT absolute 2-cycles.

Concretely: P_v ∘ ∂_2(x) = 0 (projection of boundary onto v-arcs vanishes), but ∂_2(x) ≠ 0.

### The correct framing
ker(∂_2^rel) = dim(Ω_2) - rk(M) - dim(V_2), where M = ∂_2|_{Ω_2} restricted to rows of v-arcs.

This correctly counts the preimage of V_1 in Ω_2.

### Impact
- **HYP-213 is REFUTED**: H_2(T, T\v) > 0 for many (T,v) pairs at n ≥ 4.
  - n=4: 16/256 pairs (6.25%)
  - n=5: 840/5120 pairs (16.4%)
  - n=6: 35,328/196,608 pairs (18%)
- The proposed inductive proof of β_2 = 0 via H_2(T,T\v) = 0 does NOT work.
- However, β_2 = 0 itself is NOT affected — it remains computationally verified.
- The connecting map δ: H_2(T,T\v) → H_1(T\v) is always injective (verified n=4,5), consistent with β_2 = 0 via the long exact sequence.

### Lesson
When computing relative homology H_*(X, A) via quotient complexes:
1. ker(∂_p^{rel}) is NOT (ker ∂_p + C_*(A)) / C_*(A).
2. ker(∂_p^{rel}) = ∂_p^{-1}(C_{p-1}(A)) / C_p(A).
3. These differ whenever there are elements whose boundary is nonzero but lands in the sub-complex.
4. Always verify relative homology against the long exact sequence.

---

## MISTAKE-018: beta_3 <= 1 Assumed for All Tournaments

**Date discovered:** 2026-03-09 (kind-pasteur-S48)
**Found by:** kind-pasteur-S48 via extended sampling at n=8 (5000 random tournaments)
**Affects:** THM-123 (was THM-110) proof architecture, HYP-371b, HYP-375, HYP-342, HYP-380, HYP-393 scope

### What was assumed
Multiple hypotheses and proof strategies assumed beta_3 <= 1 for ALL tournaments:
- HYP-371b: "beta_3=2 impossible"
- HYP-375: "beta_3 <= 1 at n=9"
- THM-123 proof architecture: Claims I, II, III designed to prove beta_3 <= 1
- The opus exhaustive proof at n=7 was incorrectly assumed to generalize

### Why it was wrong
beta_3 = 2 DOES occur at n=8. Four examples found in 5000 random tournaments (rate ~0.08%):
- Profile: (1, 0, 0, 2, 0, 0, 0, 0) — two independent H_3 generators
- Scores: (2,3,3,3,4,4,4,5) and (3,3,3,3,4,4,4,4) — near-regular
- Confirmed by BOTH max_p=5 and max_p=7 in full_chain_complex_modp (mod-p exact)
- All b3=2 tournaments have good vertices (b3(T\v)=0 for some v)

Previous sampling (200 at n=9, 100 at n=8) was insufficient to detect 0.08% rate.

### The correct framing
- beta_3 <= 1 is proved ONLY at n <= 7 (exhaustive, HYP-393)
- beta_3 = 2 at n=8 (confirmed, 4/5000)
- beta_3 may grow further at n >= 9

### Impact
- THM-123 proof architecture is valid ONLY at n <= 7
- Claims I (i_*-injectivity) also FAIL at n=8 (13 violations in 5000 trials, even with b4=0)
- Claim III (consecutive seesaw) FAILS at n=8 (beta_3+beta_4 coexistence)
- The beta_3 <= 1 bound is a SMALL-n PHENOMENON, not a universal property

### Lesson
1. Small sample sizes (100-200) cannot detect 0.1% phenomena. Use 5000+ for rare events.
2. Properties proved exhaustively at n<=7 do NOT automatically extend to n>=8.
3. n=8 is a critical threshold where many path homology structural properties break down:
   consecutive seesaw, i_*-injectivity, beta_3<=1, bad vertex acyclicity.

## MISTAKE-019: Int64 Overflow in Chained Numpy Matrix Multiplication

**Date discovered:** 2026-03-10 (kind-pasteur-S50)
**Found by:** kind-pasteur-S50 via comparison of two K_tv computation methods
**Affects:** opus-S59's tv_cycle_structure.py (Ghost Cycle "failures" are spurious), potentially any script using `A @ B @ C % PRIME` pattern

### What was assumed
`D3 @ (tv_omega @ ob3).T % PRIME` safely computes the matrix product mod PRIME.

### Why it was wrong
With `RANK_PRIME = 2^31 - 1`:
- `tv_omega @ ob3` produces entries up to ~4.6 × 10^18 (fits in int64, max 9.2 × 10^18)
- BUT these are NOT reduced mod PRIME — entries can be >> PRIME
- `D3 @ X.T` then involves products up to `(2^31) * (4.6e18) = 9.9e27`, massively exceeding int64 max
- Result: silent int64 overflow → wrong matrix entries → wrong rank → wrong K_tv

This caused opus-S59's tv_cycle_structure.py to report Ghost Cycle failures in 14/504 pairs at n=7 and 11/304 at n=8. ALL of these "failures" are arithmetic artifacts.

### The correct framing
ALWAYS reduce mod PRIME between chained matrix multiplications:
```python
# WRONG (can overflow):
result = A @ B @ C % PRIME

# RIGHT:
temp = A @ B % PRIME
result = temp @ C % PRIME

# BEST: use the new safe utility:
from tournament_utils import matmul_mod
result = matmul_mod(matmul_mod(A, B), C)
```

The `matmul_mod()` function in tournament_utils.py automatically chunks the inner dimension to prevent overflow, even for single multiplications with large entries.

### Impact
- Ghost Cycle (K_tv = B_tv) HOLDS universally at n ≤ 8 (0 real failures in 1000+ tests)
- HYP-408 (codim-1 universality) remains computationally verified at n ≤ 8
- No real mathematical result is affected; only the false "counterexamples" are invalidated

### Lesson
1. NEVER chain numpy `@` without intermediate `% PRIME` when PRIME ≈ 2^31
2. Use `matmul_mod()` from tournament_utils.py for all modular matrix arithmetic
3. When two equivalent computations disagree, suspect numerical issues before mathematical failure

## MISTAKE-020: Truncated Chain Complex Gives False Betti Numbers at Top Degree

**Date discovered:** 2026-03-10
**Found by:** kind-pasteur-S50

### What was assumed
Using `full_chain_complex_modp(A, n, max_p=6)` for n=8 tournaments, opus-S59 reported β_6 nonzero for 89.8% of tournaments (HYP-420), with values ranging 1-25.

### Why it was wrong
With `max_p=6`, the computation gives β_6 = ker(d_6) - ranks.get(7, 0). Since degree 7 is not computed, `ranks.get(7, 0)` returns 0. The reported "β_6" is actually just dim(ker d_6), NOT the true Betti number.

With `max_p=7` (full complex): d_7 is injective on Omega_7, and rk(d_7) = ker(d_6) EXACTLY for all 50 tested tournaments. True β_6 = 0 always.

### The correct framing
The Betti number at the highest computed degree is always an UPPER BOUND (missing the image from the next degree). For n-vertex tournaments, always use `max_p=n-1` to get correct Betti numbers, especially at degrees n-2 and n-1.

Correct results: β_{n-1} = β_{n-2} = 0 for ALL tournaments at n=3-8 (HYP-423, HYP-424). The top boundary map d_{n-1} is always injective.

### Impact
- HYP-420 is FALSE. β_{n-2} is NOT generically nonzero at n=8.
- The "β_6 among β_3=1" distribution in opus's beta4_at_n7.out is entirely artifactual.
- All lower-degree Betti numbers (β_0 through β_5) from that computation are correct.

### Lesson
1. ALWAYS use max_p=n-1 when computing Betti numbers to avoid truncation artifacts
2. Betti at the max computed degree is an UPPER BOUND (Betti at max_deg-1 and below are exact)
3. When β at max_deg seems surprisingly large/nonzero, check if im(d_{max_deg+1}) is missing

---

## MISTAKE-019b: THM-136 Sign Convention Error

**Date discovered:** 2026-03-12
**Found by:** kind-pasteur-S57
**Affects:** THM-136 formula statement (not the verbal description or proof mechanism)

### What was assumed
The trace alternation sign formula was stated as:
`sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-3)/2}`

### Why it was wrong
Direct computation at p=7,11,19,23 shows:
- k=5: Delta > 0 (positive), but (-1)^{(5-3)/2} = (-1)^1 = -1 (WRONG)
- k=7: Delta < 0 (negative), but (-1)^{(7-3)/2} = (-1)^2 = +1 (WRONG)

The formula gives the OPPOSITE sign for every k.

### The correct framing
`sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-1)/2}`

Equivalently: positive for k = 1 mod 4, negative for k = 3 mod 4.
Verified with 1218+ individual (k,p) tests, zero failures.

Note: the VERBAL description in THM-136 was always correct ("k=1 mod 4: Paley wins").
Only the symbolic formula was off by one power.

### Impact
- Formula in THM-136 theorem file CORRECTED
- No downstream impact: all proofs used the verbal description, not the formula
- The algebraic proof (kind-pasteur-S57) uses the correct convention throughout

## MISTAKE-021: S70 "GLMY Betti Numbers" Use Wrong Chain Complex

**Date discovered:** 2026-03-13
**Found by:** opus-2026-03-13-S71
**Affects:** ALL scripts from S70 session: betti_omega_connection.py, betti_divisibility.py, per_eigenspace_betti.py, per_eig_betti_n9.py, and all results/theorems derived from them (THM-154, eigenspace Betti uniformity)

### What was assumed
The S70 scripts computed "GLMY path homology Betti numbers" using:
- Allowed paths = "regular paths" (v_i→v_{i+1} AND v_{i-1}→v_{i+1})
- Boundary = interior-only deletion (indices 1 to m-1)
Results were called "GLMY Betti numbers" and compared with GLMY literature.

### Why it was wrong
The actual GLMY path homology uses:
- Allowed paths = directed paths (v_i→v_{i+1} only, NO skip-one requirement)
- Boundary = full vertex deletion (indices 0 to m), but restricted to Ω_m subspace
- Ω_m = {u ∈ A_m : ∂u has all components in A_{m-1}}

**These give DIFFERENT chain complexes with different Betti numbers:**
- Paley P_7 GLMY: β = [1,0,0,0,6,0,0], dim(A_2)=63
- Paley P_7 S70:  β = [7,0,0,21,21,21,21], dim(A_2)=21

The "regular path + interior boundary" complex IS a valid chain complex
(d²=0 verified), but it is NOT standard GLMY path homology.

### The correct framing
There are TWO distinct valid chain complexes for tournaments:

1. **GLMY Path Homology** (standard): directed paths, full boundary on Ω_m.
   Implemented correctly in path_homology_v2.py.
   β_0 = 1 for all tournaments. β_2 = 0 for all tested tournaments (n≤8).

2. **Tournament Regular Homology (TRH)** (novel?): regular paths, interior boundary.
   Used in S70 scripts. β_0 = n for all tournaments on n vertices.
   Has eigenspace Betti uniformity and divisibility by n for circulants.

Both are valid mathematical objects. But they should not be conflated.

### Impact
- THM-154 (Betti divisibility) applies to TRH, not GLMY
- Eigenspace Betti uniformity applies to TRH, not GLMY
- β_2=0 for all tournaments holds for BOTH (GLMY verified n≤8, TRH verified n≤8)
- The S70 "per-eigenspace Betti" results are self-consistent but not GLMY
- The S38-S41 β_2=0 results (from path_homology_v2.py) are correct GLMY
- circulant_homology.py implements yet another convention (full boundary on regular paths) which is NEITHER GLMY nor TRH

### Lesson
ALWAYS verify which chain complex you're computing. The three ingredients
(allowed paths, boundary convention, Ω subspace) must be consistent.
When reading "path homology" results, check which convention is used.

---

## MISTAKE-019c: TWO bugs in independent set backtracking algorithm

**Date discovered:** 2026-03-13, kind-pasteur-S60
**Found by:** kind-pasteur
**Affects:** alpha3_p7_only.py, alpha3_moment_analysis.py, overlap_weight_analysis.py, H_energy_decomposition.py, cycle_walsh_decomposition.py, moment_cancellation_mechanism.py, overlap_gauss_bridge.py, alpha_directed_p11.py, alpha_full_p11.py, alpha2_direct_verify.py, backtrack_debug.py (ALL files with independent set enumeration)

### Bug 1: Missing vertex 0 (`backtrack(0,0,0)` should be `backtrack(-1,0,0)`)

The backtracking function `backtrack(v, mask, size)` iterates `for w in range(v+1, n)`. When called with `v=0`, the loop starts at `w=1`, SKIPPING vertex 0 entirely. This undercounts all alpha_j.

**Fix:** Call `backtrack(-1, 0, 0)` so the loop starts at `w=0`.

### Bug 2: Skipping consecutive indices (`backtrack(w+1, ...)` should be `backtrack(w, ...)`)

The recursive call `backtrack(w + 1, mask | nbr[w], size + 1)` passes `v = w+1`. At the next level, the loop starts at `w' = v+1 = w+2`, SKIPPING index `w+1`. This means any independent set containing cycles with consecutive indices is missed.

**Fix:** Change to `backtrack(w, mask | nbr[w], size + 1)`. Then the next level's loop starts at `w+1`, correctly considering all higher indices.

### Concrete example

At p=7, Interval tournament S=[1,2,3]:
- 59 directed cycles, 14 disjoint (3,3)-pairs (correct)
- Bug 2: Pair (5,6) = ({0,3,6}, {1,2,5}) has consecutive indices and was SKIPPED
- Backtracking gave alpha_2=13 instead of 14, H=171 instead of 175
- Held-Karp gives H=175 (correct)

### Impact
- All previous alpha_j values from backtracking are SUSPECT
- THM-027 alpha_2 values at p=7 need recheck (Paley alpha_2=7 was coincidentally correct because no consecutive disjoint pairs)
- Any H derived from backtracking alpha may be wrong

### Lesson
When implementing independent set enumeration via backtracking, the recursive call after selecting vertex w should pass v=w (NOT v=w+1). The `range(v+1, n)` in the next level already excludes w.

---

## MISTAKE-022: Sparse Gaussian Elimination Fill-In Bug

**Date discovered:** 2026-03-13, opus-S71c (9th context window)
**Found by:** opus, when k=0 eigenspace Betti numbers came out negative
**Affects:** p19_omega5_sparse.py, p23_omega5_sparse.py, p31_omega5_sparse.py, p43_omega5_sparse.py (ALL scripts using the sparse Gaussian pattern with single-pass row iteration)

### What was assumed
The sparse Gaussian elimination iterated over `sorted(col.keys())` once, subtracting each matching pivot. This should correctly eliminate all pivot contributions.

### Why it was wrong
When subtracting a pivot at row `r`, the pivot vector has entries at rows `r' > r` (fill-in). These new entries at rows NOT in the original column are never checked against existing pivots at those rows, because the `sorted(col.keys())` list was computed BEFORE the subtraction and doesn't include fill-in entries.

Concrete example: column has entries at rows {3, 7}. Pivot at row 3 has entries at {3, 5, 7}. After subtracting pivot 3, the column has entries at {5, 7}. Row 5 was NOT in the original sorted list, so even if there's a pivot at row 5, it's never subtracted. This causes the rank to be OVERCOUNTED (some columns that should reduce to zero don't).

### The correct framing
After any pivot subtraction, restart the row scan from the beginning (or at least from the newly-created entry). A simple fix: wrap the elimination loop in `while changed: ... break after subtraction`.

### Impact
- **P_19 Omega_5 was 12602 (WRONG), correct is 23832**
- **P_23 Omega_5 was 50715 (WRONG), correct is 78430**
- **P_31 Omega_5 was 252065 (WRONG), correct is 456330**
- **P_43 Omega_5 was 1429652 (WRONG), correct is 2865660**
- P_7 and P_11 were unaffected (small enough that fill-in didn't change rank)
- HYP-790 ("Omega_5 not polynomial in m") was based on wrong data — **RETRACTED**
- **CORRECTED**: Ω_5 = m(m-1)(m³-6m²+10m-2) — a **clean integer polynomial** in m!
- All formulas Ω_d for d ≤ 5 are now proven/verified

### Lesson
In sparse Gaussian elimination, fill-in from pivot subtraction can create new entries at rows that were not in the original column. These MUST be processed against their pivots. Always use a while loop that restarts after each subtraction, or maintain a priority queue of unprocessed rows.

---

## MISTAKE-023: α₁ Counts DIRECTED Odd Cycles, Not Vertex-Sets

**Date discovered:** 2026-03-14
**Found by:** opus-2026-03-14-S71d
**Affects:** two_and_three_universality.py, i3_mod3_proof.py, vandermonde_sigma_connection.py, jacobsthal_23_deep.py (first version), and any script computing I(CG, x) by counting cycle vertex-sets

### What was assumed

The independence polynomial I(Ω, x) was computed by enumerating odd cycle **vertex-sets** (frozenset of vertices), counting each set once regardless of how many distinct directed cycles it supports.

### Why it was wrong

The conflict graph Ω(T) has vertices = **directed odd cycles** (definition in definitions.md line 37). For 3-cycles in tournaments, each vertex triple supports at most 1 directed 3-cycle, so vertex-set counting is correct. But for 5-cycles and above, a single vertex-set can support **multiple** distinct directed cycles:

- Example: bits=40 at n=5, the 5-vertex set {0,1,2,3,4} supports **3** distinct directed 5-cycles
- Vertex-set method gives α₁=5 → I(2)=11, but H=15
- Directed-cycle method gives α₁=7 → I(2)=15 = H ✓

### The correct framing

When computing I(Ω, x):
1. For each vertex-set of size k, enumerate ALL distinct directed k-cycles (normalize by fixing start vertex and direction)
2. Each distinct directed cycle is a SEPARATE vertex of Ω(T)
3. Two vertices of Ω are adjacent iff the underlying vertex-sets intersect

**Exhaustive verification at n=5:**
- Vertex-set method: 184/1024 mismatches with H
- Directed-cycle method: 0/1024 mismatches with H

### Impact

- All α₁ values from scripts using vertex-set counting at n≥5 are WRONG (undercounted)
- The Vandermonde extraction results (HYP-867, HYP-868) were based on the wrong α values
- The 3/2 ratio result may still hold (it was measured within lambda fibers, not from α directly)
- The structural insights about 7→8 transition are UNAFFECTED (vertex-set counting is correct for α₂ when cycles have different sizes)
- Scripts need to be updated to use directed cycle enumeration

### Lesson

The definition says "vertices are **directed** odd cycles." For 3-cycles in tournaments, vertex-set = directed cycle (1-to-1). For k≥5 cycles, a k-vertex tournament subtournament can have multiple Hamiltonian cycles. Always enumerate directed cycles explicitly.

---

## MISTAKE-024: H=63 Falsely Claimed Permanently Forbidden

**Date discovered:** 2026-03-14
**Found by:** opus-2026-03-14-S71h (cross-referencing S86 broadcast with prior results)
**Affects:** HYP-1303, MSG-218 (S86 broadcast), MSG-139 (S86 to kind-pasteur)

### What was assumed

S86 claimed: "H=63 FORBIDDEN: 63=7×9=I(K₃,2)×I(2K₁,2). Requires K₃ component in Ω, impossible by THM-201." This was marked CONFIRMED as HYP-1303.

### Why it was wrong

The argument only blocks DISCONNECTED conflict graphs where Ω = K₃ ⊔ 2K₁. But Ω can be a CONNECTED graph with I(Ω,2)=63. Multiple prior sessions had already established:
- S65-c (MSG-084): "H=63,107,119,149 (the n=7-specific gaps) are ALL achieved at n=8"
- S71f (MSG-197): "63 achievable at n≥8"
- S71g (MSG-201): "H=63 found at n=8 (27/100k)"
- hspectrum_density.out: "63 = 7*9 -- ACHIEVED at n=8 (not permanent)"

The S86 session re-derived an incomplete argument without checking these earlier verified results.

### The correct framing

H=63 is NOT permanently forbidden. It is a temporary gap at n=7 (like 107, 119, 149) that IS achieved at n=8. The ONLY permanently forbidden H values proved so far are {7, 21}.

The disconnected decomposition 63=7×9 is correctly blocked, but connected Ω graphs with I=63 exist and can apparently be realized as Ω(T) at n=8.

### Impact

- HYP-1303 changed to REFUTED
- HYP-1295 changed to REFUTED
- The S86 broadcast claim "all three known forbidden values {7,21,63} are now explained" is WRONG
- Only {7, 21} are explained as permanently forbidden

### Lesson

Always check prior session results before claiming a new proof. The H-spectrum density analysis (hspectrum_density.out) had already settled this question computationally. Cross-reference before broadcasting claims.

---

## MISTAKE-025: S112 W(8) Value Off By 8

**Date discovered:** 2026-03-16 (opus-S90 continuation session)
**Found by:** opus-S90, via independent brute-force verification
**Affects:** kind-pasteur-S112 W(n) sequence, D_n(2) computations

### What was claimed
kind-pasteur-S112 reported W(n) = 1, 2, 8, 32, 158, 928, 6350, **49760**, 439766 for n=1..9.

### Why it was wrong
Independent brute-force computation (iterating over all 8! = 40320 permutations) gives W(8) = **49752**, not 49760. The error is exactly +8 = 2³. The S89c C-program DP computation also gives 49752, confirming the brute-force result.

### The correct values
W(n) = 1, 2, 8, 32, 158, 928, 6350, **49752**, 439670 (from S89c DP).

### Impact assessment
- **S89c values (opus):** CORRECT through n=27 (computed by bitmask DP in C).
- **S112 values (kind-pasteur):** INCORRECT at n=8 by +8. Values at n≤7 match. Values at n≥9 need reverification against S89c.
- **OEIS submission:** Use S89c values, not S112.
- **CV² formula and g_k polynomials:** UNAFFECTED (derived independently of W(n) enumeration).

### Source
opus-2026-03-16 (S90 continuation): brute-force W(8) verification via Python permutation enumeration.

### Lesson
When two independent computations disagree, trust the one with the simpler algorithm (brute force) over the one with more complex logic. The discrepancy of exactly 2³ suggests a boundary condition or off-by-one error in the S112 computation, not a random bug.

---

## MISTAKE-026: Cross-Ratio of Cayley Orbit Initially Claimed as 8/7

**Date discovered:** 2026-03-15 (code review during S90)
**Found by:** opus-S90 code review agent
**Affects:** monad_cayley_s90c.py, commit messages

### What was claimed
The cross-ratio of the Q-orbit of x=2 was initially computed as 8/7, using the WRONG orbit point (3 instead of -3).

### Why it was wrong
Q(2) = (1+2)/(1-2) = 3/(-1) = **-3**, not 3. The orbit is {2, **-3**, -1/2, 1/3}. The cross-ratio CR(2, -3, -1/2, 1/3) = **2**, not 8/7.

### The correct value
Cross-ratio = 2 = the OCF fugacity itself. This is MORE meaningful than 8/7 — the cross-ratio equals the evaluation point.

### Impact
- The narrative about "tournament constant 8/7" in commit messages is wrong.
- The correct "tournament constant" is 2 (the fugacity).
- Script monad_cayley_s90c.py has been corrected.

### Source
Code review agent, opus-2026-03-15-S90.

---

## MISTAKE-018b: THM-225 "Universal Top Eigenvalue = n" is FALSE at n ≥ 9

**NOTE:** This was originally numbered MISTAKE-018 from a different branch, causing a collision with MISTAKE-018 (beta_3 <= 1). Renumbered to 018b by opus-2026-04-01-S1.

**Date discovered:** 2026-03-15
**Found by:** opus-S72d
**Affects:** THM-225, HYP-1594

### What was assumed

That the top eigenvalue of C_T^TC_T equals n for ALL tournaments on n vertices (verified exhaustively at n=5, sampled at n=6). This was stated as a theorem.

### Why it was wrong

The proof strategy required rank(C_R) < r = (n-1)(n-2)/2. At n ≤ 8, this holds because max c₃ < r (the number of cyclic triples never exceeds the rank of the full constraint matrix). At n = 9, c₃ can reach 30 while r = 28, and for ~0.1% of tournaments, rank(C_R) achieves its maximum r, leaving ker(C_R) ∩ im(C^T) = {0}. The top eigenvalue then drops to ~8.84-8.94.

### The correct framing

THM-225 holds for n ≤ 8 (PROVED for n ≤ 5 exhaustive, sampled at n=6,7,8 with 0 violations from 20000 samples each). It FAILS at n ≥ 9. The condition for top eigenvalue = n is rank(C_R) < (n-1)(n-2)/2.

### Impact

The spectral T/R duality (C_T^TC_T + C_R^TC_R = n·P) and the 3/n bridge framework remain valid. The universal top eigenvalue was a COROLLARY that holds only when the cyclic boundaries don't span the full constraint space.

### Lesson

When verifying at n=5 and n=6, the parameter regime (c₃ < r always) hid the failure mode. Always check at the CROSSOVER point where qualitative behavior changes. For rank arguments, the critical n is where max c₃ first exceeds r.

---

## MISTAKE-027: THM-080 Amplitude Table Wrong at n=9

**Date discovered:** 2026-03-16
**Found by:** opus-2026-03-16-S73
**Affects:** THM-080 amplitude table (lines 156-161), not the formula itself

### What was assumed

The amplitude table in THM-080 listed n=9 entries as: deg 1 (s=0) = 3/2, deg 3 (s=0) = 3/8, deg 3 (s=1) = 3/4, deg 5 (s=0) = 1/16, deg 7 = 1/128.

### Why it was wrong

The stated formula is |hat{M}[S]| = 2^s × (n-2-|S|)!/2^{n-2}. At n=9:
- d=1, s=0: formula gives 6!/128 = **45/8**, not 3/2
- d=3, s=0: formula gives 4!/128 = **3/16**, not 3/8
- d=3, s=1: formula gives 2×3/16 = **3/8**, not 3/4
- d=5, s=0: formula gives 2!/128 = **1/64**, not 1/16
- d=7, s=0: formula gives 0!/128 = 1/128 ✓ (only correct entry)

The formula works perfectly at n=3, 5, 7 — only n=9 has errors. The d=3 and d=5 wrong values are the CORRECT formula values but with s shifted up by 1 or 2 (unrooted component miscount). The d=1 entry (3/2) doesn't correspond to any valid s value (45/8 × 2^s ≠ 3/2 for any integer s).

### The correct framing

The formula is correct. The table had a transcription error at n=9 only. The n=9 verification was "partial" and apparently didn't catch the table/formula mismatch. Corrected amplitude table is in THM-080.

### Impact

Low — the formula itself is correct and was verified computationally at n=5 (exhaustive), n=6 (exhaustive), n=7 (20/20). Only the summary table for n=9 was wrong. No downstream results depend on the specific n=9 table values.

### Lesson

**This is MISTAKE-013b (the original missing 2^s) echoing forward.** The 2^s correction was caught and fixed at n=7, but the n=9 table values were apparently populated from a pre-correction computation (or from hand calculation that repeated the component-counting error at higher n). Always re-derive table entries from the corrected formula rather than carrying forward values from partial computation.

This is also a meta-lesson about the amplitude table itself: it was the only place in THM-080 where specific numerical values were stated without being individually verified against the formula. The formula (analytically proved) was trustworthy; the table (hand-calculated) was not.

---

## MISTAKE-028: Mersenne / k-nacci Numbers Falsely Claimed to Control Forbidden H Values

**Date discovered:** 2026-03-17
**Found by:** opus-2026-03-17-S74 (forbidden values audit)
**Affects:** casual-writeup.md, formal-writeup.md, substack-hooks.md (Hook N), HYP-1600, HYP-1618 (original), HYP-1623, HYP-1624, riemann-zeta-tournament.md, multiple results files

### What was assumed

Multiple writeups and hypotheses claimed: "The k-nacci Mersenne identity connects forbidden H values (7 = 2^3 - 1, 31 = 2^5 - 1, 127 = 2^7 - 1) to Mersenne primes via k-nacci transfer matrices." The original HYP-1618 claimed "ζ(-3) = 7" (standard Riemann zeta). Various scripts called 31 "forbidden."

### Why it was wrong

1. **31 is achievable** at n=6 (tournament bits=146, verified exhaustively).
2. **63 is achievable** at n=8 (already documented in MISTAKE-024).
3. **127 is achievable** at n=7.
4. The standard Riemann ζ(-3) = 1/120, NOT 7.
5. The tribonacci trace sequence [1, 3, 7, 11, 21, 39, 71, 131, ...] contains both forbidden values (7, 21) AND achievable values (1, 3, 11, 39, 71, 131, ...). The k-nacci trace hitting 7 and 21 is a numerical coincidence, not a causal mechanism.

### The correct framing

**Only H=7 and H=21 are permanently forbidden** (proved: THM-029, THM-079). The actual mechanisms are:
- H=7: 3 pairwise-conflicting cycles always force additional cycles (THM-029)
- H=21: all OCF decompositions blocked by tournament forcing (THM-079, 464-line proof)

Best characterizations of {7, 21}:
- {7·3⁰, 7·3¹}: the 7-obstruction has nilpotency 2 (HYP-1231). 7·3² = 63 is achievable.
- {Φ₃(2), Φ₃(4)}: third cyclotomic polynomial at even args (HYP-1317). Φ₃(6) = 43 is achievable.
- Both have I-polynomials factoring through I(K₃, x) = (1+3x) (HYP-1315).

THM-227 (k-nacci Mersenne) is a valid theorem about transfer matrices. It just doesn't characterize forbidden H values.

### Impact

Medium — the false claim propagated through 6+ files across multiple sessions. All have been corrected. No theorems or proofs are affected (the actual forbidden value proofs THM-029 and THM-079 are correct and don't use the Mersenne connection).

### Lesson

**One data point is not a pattern.** The entire false extrapolation rested on a single observation: 7 = 2³ - 1 is both a Mersenne number and forbidden. From this, it was incorrectly inferred that other Mersenne numbers (31, 127) are also forbidden. A simple check (is H=31 achieved at n=6?) would have caught this immediately.

This is a variant of MISTAKE-024 (H=63 falsely claimed forbidden) — the same class of error, just with a different numerological motivation. The meta-lesson: when claiming a numerical pattern "explains" something, verify it at the NEXT case before asserting it as a principle.

---

## MISTAKE-029: Formula E = (T - D)/2 is WRONG for the meta-graph edge count

**Date discovered:** 2026-03-23
**Found by:** opus-2026-03-23-S211
**Affects:** degeneracy_second_moment_s210.py, all claims that E = (3T - S_2)/4

### What was assumed
The meta-graph edge count formula E = (T - D)/2 was claimed in S210, where T = total arc-orbits and D = sum C(mult(C→D), 2) is the degeneracy. The derived formula E = (3T - S_2)/4 was presented as the "keystone" edge count formula, and the reverse-engineered S_2 sequence was reported as new.

### Why it was wrong
The formula ignores **self-loop orbits** and **directed multi-edge excess**. The correct decomposition is:

T = SL + 2E + excess_cross

where:
- SL = sum_C mult(C→C) = total self-loop arc-orbits
- excess_cross = sum_{{C,D}} (mult(C→D) + mult(D→C) - 2) for connected C≠D

So: **E = (T - SL - excess_cross) / 2**, NOT (T - D) / 2.

At n=3,4,5: the formula happened to give correct answers because SL + excess = D exactly (coincidence). At n=6: SL + excess = 58 + 66 = 124, but D = 122, so E_wrong = 291 ≠ E_actual = 290.

The "reverse-engineered" S_2 values for n≥6 were derived from this wrong formula and are therefore incorrect. The actual S_2 at n=6 is 948, not 952.

### The correct framing
E(n) must be computed directly from the meta-graph adjacency (F matrix), not from aggregate orbit statistics. There is no known closed-form expression for E(n) in terms of T, D, or S_2. The quantities T(n) and S_2(n) give orbit-level statistics but cannot determine which pairs of classes are actually adjacent.

### Impact
- The "keystone formula" E = (3T - S_2)/4 is invalid
- The S_2 sequence 8, 28, 144, 952, 10392, 200220, 7018596 from S210 is wrong at n≥6
- The correct S_2 at n=6 is 948 (from direct orbit computation)
- The gap sequence G(n) = T - 2E = 2, 6, 28, 124, 740, 5966, 85698 IS correct and novel
- All independently computed E values (E(3..8) = 1, 5, 30, 290, 4086, 91161) remain valid

### Lesson
When a formula passes at small n, always verify at the next n where new phenomena emerge. At n≤5, every class has SL + excess = D (a coincidence), so the formula appeared correct. At n=6, the coincidence breaks. **Integer division can mask off-by-one errors**: at n=3, (T-D)/2 = 3/2 = 1.5, which rounded to 1 via `//` and happened to match E=1.

---

## MISTAKE-030: "SL_orbits" is a misnomer — it includes multi-edge orbits, not just self-loops

**Date discovered:** 2026-03-23
**Found by:** Devil's advocate audit (opus-2026-03-23-S246), confirmed by opus-S245
**Affects:** burnside_edge_verify_s242.py, recursive_sl_s244.py, all scripts using "SL_orbits"

### What was assumed
The quantity "SL_orbits" = edge_orbits - E(G_n) was treated as counting self-loop edge orbits (orbits where both endpoints are in the same iso class).

### Why it was wrong
"SL_orbits" actually counts ALL non-simple-edge orbits: self-loop orbits PLUS multi-edge orbits (additional orbits connecting already-connected class pairs). At n=5: true self-loop edge orbits = 14, but "SL_orbits" = 20. The difference of 6 is multi-edge orbits.

At n=3,4: multi = 0, so the values coincide — masking the error (same pattern as MISTAKE-029).

### The correct framing
- **gap_orbits** = edge_orbits - E = self_loop_orbits + multi_orbits (RENAME from "SL_orbits")
- **self_loop_orbits** = #{S_n-orbits on {T, T^e} with T ≅ T^e} = 2, 5, 14, ... (computed via Burnside)
- **multi_orbits** = #{edge orbits connecting already-counted class pairs} = 0, 0, 6, ...

### Impact
- The recurrence search for "SL_orbits" in S242/S244 was wasted effort on a DERIVED quantity (= T/2 + (n-2)! - E). Any pattern found is just a pattern in E in disguise.
- The formula edge_orbits = T/2 + (n-2)! is CORRECT and independently valuable.
- Future work should target E(G_n) directly, not the gap.

### Lesson
**Name quantities precisely.** "SL_orbits" was never defined as "self-loop edge orbits" — it was defined as "edge_orbits - E" and then ASSUMED to count self-loops. The assumption failed at n=5. Always verify definitions against direct computation before building analysis on them.

---

## MISTAKE-031: Tiling complement ≠ tournament complement

**Date discovered:** 2026-03-24
**Found by:** Devil's advocate audit (kind-pasteur-S20ex)
**Affects:** wiggly_metagraph_deep_s20ev.py, aw_precision_s20ew.py, wiggly_patterns_s20eq.py, unified_weights_s20et.py

### What was assumed
Flipping all TILE bits (`mask ^ ((1<<m)-1)`) gives the complement tournament.

### Why it was wrong
The tiling model has m = C(n-1,2) tiles (non-base-path arcs). Flipping tile bits only reverses these arcs, leaving the n-1 base path arcs unchanged. The true tournament complement reverses ALL C(n,2) arcs. These produce different tournaments.

### Impact
V_merged was wrong at n>=5: got 9 (should be 10) at n=5, 41 (should be 34) at n=6. All spectral analysis, W/A comparisons, and eigenvector correlations in affected scripts were computed on the WRONG quotient graph. Corrected in wiggly_corrected_s20ex.py.

### Lesson
When working in the tiling model (fixed base path), always compute the complement via the ADJACENCY MATRIX (reverse all arcs), not via bit flipping on tiles.

---

## MISTAKE-032: Grid-symmetric fraction formula was wrong

**Date discovered:** 2026-03-24
**Found by:** Devil's advocate audit (kind-pasteur-S20ex)
**Affects:** CLAUDE.md (line about "Grid-sym fraction = exactly 2^{-(n-2)}")

### What was assumed
The fraction of grid-symmetric tilings is 2^{-(n-2)} for all n.

### Why it was wrong
The correct formula is 2^{(floor((n-1)/2) - C(n-1,2))/2}, giving exponents 0, -1, -2, -4, -6, -9 for n=3..8. The formula 2^{-(n-2)} gives -1, -2, -3, -4, -5, -6 which only matches at n=5,6.

### Impact
Claims about blue fraction being exactly 2^{-(n-2)} per class are wrong. The correct formula accounts for the number of fixed tiles on the anti-diagonal of the staircase.

### Lesson
Always verify formulas at multiple n values, not just the first few where coincidences can mask errors.

---

## MISTAKE-033: Confused complement-tiling with complement-tournament in blue/black analysis

**Date discovered:** 2026-03-24
**Found by:** User correction (opus-S295)
**Affects:** three_graphs_s295.py, wiggly_vs_lines_s275.py reasoning

### What was assumed
Blue/black lines were modeled as connecting T to T^op (the tournament complement, flipping ALL C(n,2) arcs including base-path arcs). This led to the claim that blue/black "lives outside Q_m" and has zero cross-class edges in the merged meta-graph.

### Why it was wrong
In the tournament-tiling-explorer, a blue/black line connects a TILING to its COMPLEMENT TILING = flip all m tiles (XOR with 2^m - 1). This stays INSIDE Q_m. The complement tiling gives a tournament where all non-base-path arcs are reversed but base-path arcs are PRESERVED. This is NOT T^op (which reverses ALL arcs).

The complement tiling IS at Hamming distance m in Q_m. It gives a different labeled tournament that may be in a DIFFERENT iso class. Blue/black lines DO create cross-class edges in both the unmerged and merged meta-graphs.

### The correct framing
- **Complement TILING** = flip all m tiles = bits XOR (2^m - 1). Stays in Q_m. THIS is what blue/black lines are.
- **Complement TOURNAMENT** (T^op) = flip all C(n,2) arcs. Leaves Q_m (changes base-path arcs). NOT the same as complement tiling.
- Blue/black lines ARE in Q_m, they ARE waggly lines (at distance m), and they DO connect different iso classes.
- The S295 analysis incorrectly modeled blue/black by computing T^op instead of the complement tiling.

### Impact
The three_graphs_s295.py blue/black weight matrix is WRONG. The claim "blue/black is purely diagonal" is FALSE. Must recompute using the correct definition: complement tiling = XOR all tile bits.

### Lesson
ALWAYS check definitions against the actual explorer behavior. The tiling complement (flip tiles) and tournament complement (flip all arcs) are different operations. In the tiling model with fixed base path, flipping all tiles does NOT give T^op.

---

## MISTAKE-034: Band-limitedness at m/2 does NOT hold at n=5

**Date discovered:** 2026-03-25
**Found by:** kind-pasteur-2026-03-25-S1
**Affects:** h-is-band-limited.md (opus-S306), OPEN-Q-040 item 2

### What was assumed
"H is EXACTLY zero in upper Walsh spectrum (k >= m/2). PROVED at n=5,6." (From OPEN-Q-040 and h-is-band-limited.md reflection.)

### Why it was wrong
At n=5 (m=6): the Walsh degree of H is 4 = 2*floor((5-1)/2). Since m/2 = 3, the Walsh degree EXCEEDS m/2. There are 7 nonzero Walsh coefficients at weight 4, and alpha_4 = sum of Walsh coefficients at weight 4 = 0.375 != 0.

Additionally, complement symmetry H(t) = H(~t) FAILS in the tiling model because flipping all tile bits is NOT T^op (base-path arcs are fixed). This means odd-weight Walsh coefficients are nonzero (17 at n=5, 907 at n=7).

### The correct framing
- Walsh degree = 2*floor((n-1)/2) for ALL n >= 4 (THM-260, proved via THM-076)
- Band-limitedness at m/2 holds for **n >= 6** (since 2*floor((n-1)/2) < C(n-1,2)/2 iff n >= 6)
- At n=4,5: Walsh degree exceeds m/2 — NOT band-limited at midpoint
- In the tiling model, both odd and even Walsh weights can be nonzero

### Impact
The "upper half vanishes" claim in h-is-band-limited.md needs correction at n=5. The main qualitative finding (H is low-frequency, concentrated in lower Walsh spectrum) is correct and gets STRONGER as n grows. The asymptotic ratio degree/m -> 0 still holds.

### Lesson
When making claims about "all n," verify at the boundary cases (smallest n). The n=5 case is special because m = C(4,2) = 6 is comparable in size to n-1 = 4. For n >= 6, the quadratic growth of m dominates the linear growth of the Walsh degree.

---

## MISTAKE-035: "G_n is a DAG under H-gradient" — False Claim Propagated Across Repo

**Date discovered:** 2026-04-01
**Found by:** opus-2026-04-01-S1 (systematic audit)
**Affects:** CLAUDE.md, OPEN-QUESTIONS.md, 4 reflection files, paper draft, ~20 agent messages, gn_merged_cascade_s221.py (hardcoded output), local_gradient_s186.py (hardcoded CONFIRMED)

### What was claimed
"The meta-graph G_n is a DAG under H-gradient (0 downhill edges, verified n=3..7)" — CLAUDE.md line 326 (pre-fix). OPEN-QUESTIONS.md claimed "HOLDS at n=3..8."

### Why it was wrong
THREE distinct errors compounded:

1. **Trivially true claim conflated with nontrivial property.** For ANY undirected graph with a real-valued function on vertices, orient edges by function value → the result is always a DAG (modulo level edges). This was explicitly noted in `meta_graph_deep_s181.py` lines 366-368 but the insight was never propagated. The REAL nontrivial question is about **level edges** (same H, different class).

2. **Level edges exist from n=5 onward.** G_n level edges: 0, 0, 1, 15, 136 for n=3..7. G_n/Z_2 level edges: 0, 0, 1, 5, 71 for n=3..7. The graph is NOT a strict DAG from n=5 onward.

3. **H-decreasing edges exist at n=7.** `merged_n7_deep_s20co.out` shows: G_7 has uphill=2988, downhill=962, level=136. G_7/Z_2 has uphill=1633, downhill=419, level=71. The "downhill" count here reflects edges where the class with more neighbors (higher index) has LOWER H — these are real H-reversals in the metagraph. `gap_inventory_s196.py` correctly listed this as REFUTED.

4. **Hardcoded output bugs.** `gn_merged_cascade_s221.py` line 487 prints "DAG: 0 H-decreasing edges (all n)" unconditionally, even though its own data (line 68 of output) shows "DAG: Y, Y, N, N, N, N" for n=3..8. `local_gradient_s186.py` prints "CONFIRMED: all negative-DeltaH flips stay in-class" unconditionally even when the script found counterexamples.

### The correct framing
- G_n has a **strong H-gradient**: most edges increase H. The ratio uphill/(uphill+downhill) is 100% at n≤6 (for the nontrivial edges), and ~76% at n=7.
- G_n is NOT a strict DAG from n≥5 (level edges) and has H-decreasing edges from n≥7.
- The level edge fraction stays small (~3-5%) and may decrease asymptotically.
- The H-gradient is a useful organizing principle but not an absolute law.

### Impact
- CLAUDE.md, OPEN-QUESTIONS.md, 4 reflection files, paper draft all corrected (opus-2026-04-01-S1).
- Every new agent session was reading this false claim and propagating it.
- The `unlocking-gn-at-all-n.md` file listed H-DAG as a "Proved Structural Law" — it was not proved and is not true.

### Lesson
**Three compounding failures:** (1) A trivially-true observation was mistaken for a nontrivial theorem. (2) The discoverer of the triviality (meta_graph_deep_s181.py) did not propagate the correction. (3) Later scripts hardcoded "CONFIRMED" messages that print regardless of results. When a claimed property is trivially true, that's a red flag that you're measuring the wrong thing.

---

## MISTAKE-036: Diameter conjecture diam(G_n) = n-2 is WRONG

**Date discovered:** 2026-03-23
**Found by:** kind-pasteur (gap_inventory_s196)
**Affects:** the-isomorphism-class-graph.md, merged-metagraph-invariants.md, multiple broadcast messages

### What was claimed
"Diameter of G_n is n-2" — conjectured based on n=3 (diam=1), n=4 (diam=2), n=5 (diam=3).

### Why it was wrong
At n=6: diam=4 = n-2 (still holds). At n=7: diam=**7** ≠ 5 = n-2. At n=8: diam=**8** ≠ 6. The actual growth is closer to quadratic (~n²/4), not linear. The diameter-is-feedback-arc-set.md reflection explains: diam ≈ max FAS count difference, which grows quadratically.

### The correct values
diam(G_n) = 1, 2, 3, 4, 7, 8 for n=3..8.

### Impact
- `merged-metagraph-invariants.md` self-contradicts: says "CONFIRMED" at line 84 and "REFUTED" at line 172.
- `the-isomorphism-class-graph.md` still lists "Prove diameter = n-2" as an open problem.
- Multiple broadcast messages from S170, S177, S305 assert or propose proving diam=n-2.

### Lesson
Patterns that hold for 4 consecutive values (n=3..6) can still fail at n=7. Always test at the next case before conjecturing.

---

## MISTAKE-037: H-convexity conjecture is FALSE

**Date discovered:** 2026-03-23
**Found by:** kind-pasteur-S20ch
**Affects:** gap_inventory_s196.py line 176

### What was claimed
That the H-landscape on G_n is "convex" — a specific technical condition about H values along paths in the metagraph.

### Why it was wrong
Refuted at n=6 by kind-pasteur-S20ch. Specific counterexample documented in gap_inventory_s196.py.

### Impact
Low — this was a tentative conjecture, not widely propagated.

### Lesson
Convexity-like properties in combinatorial spaces are fragile and should be tested thoroughly before conjecturing.

---

## MISTAKE-049: SC(n) = A000568(n-1) — Fabricated Identity

**Date discovered:** 2026-05-07 (oracle session)
**Found by:** oracle-2026-05-07
**Affects:** `07-reflections/product-graph-sc-spine-fractal-dimensions.md`

### What was assumed
The reflection claimed SC(n) = A000568(n-1), "verified n=2..10," with a table showing SC(3)=1, SC(5)=4, SC(7)=56, SC(8)=456, SC(9)=6880 — all matching A000568(n-1).

### Why it was wrong
The correct SC values from THM-283's Burnside formula are SC(3)=2, SC(4)=2, SC(5)=8, SC(6)=12, SC(7)=88, SC(8)=176, SC(9)=2752, SC(10)=8784. These do NOT match A000568(n-1) except at n=4,6 (coincidences). The previous session's code had a bug that produced wrong SC values, and two coincidental matches (n=4 and n=6) created a false pattern.

### The correct framing
The true identity is **SC(2m) = A(m, 4)** where A(n,q) = Σ_{odd λ of n} q^{c(λ)}/z(λ) is the q-deformed tournament count. A(n,2) = A000568(n) and A(m,4) = SC(2m). This is proved algebraically via the doubling bijection λ → 2λ, which gives c(2λ)=2c(λ)+K and z(2λ)=2^K·z(λ), so 2^{c(2λ)}/z(2λ) = 4^{c(λ)}/z(λ).

### Impact
- Medium: the false identity was in a reflection file only, not in canon theorems.
- The CORRECT identity (SC(2m)=A(m,4)) is new and provably correct.
- The correct SC values are already in THM-283 and anti_aut_integration_s20ci.out.

### Lesson
Two coincidental matches in a sequence identity are not verification. Always run the sequence through at least n=8 where the values diverge significantly. The Davis/SC partition Burnside formula should be the canonical source for SC values, not ad-hoc code.

---

## MISTAKE-050: H=63 Reintroduced as a Universal Lean Theorem

**Date discovered:** 2026-05-29
**Found by:** opus-2026-05-29-S8
**Affects:** `04-computation/lean/TournamentH7/H63.lean`, `HSpectrum.lean`, `SUBMISSION.md`, `OPEN-Q-055`, HYP-1754

### What was assumed

The Lean formalisation introduced a theorem/axiom `H_ne_sixtythree` claiming
H(T) ≠ 63 for every tournament T, citing exhaustive n≤7 evidence.
`HSpectrum.lean` bundled this into a universal forbidden trio {7,21,63}.

### Why it was wrong

This repeats MISTAKE-024. H=63 is already known to be achievable at n=8.
The S8 audit re-verified a concrete n=8 counterexample from
`h63_verify.out`:
- H(T)=63 by Held-Karp DP
- H(T)=63 by direct enumeration of all 8! vertex permutations
- Ω(T) has 31 directed odd cycles and is the complete graph K31
- Therefore I(Ω(T),2)=1+2·31=63, matching OCF

### The correct framing

H=63 is a temporary n≤7 gap, not a permanent forbidden value.
The Lean theorem is now demoted to:
`H_ne_sixtythree_le_seven (hn : n ≤ 7)`.
The universal forbidden bundle is {7,21}; the finite n≤7 bundle is {7,21,63}.

### Impact

HYP-1754 is REFUTED. OPEN-Q-055 has been corrected. Any document saying
"universally forbidden {7,21,63}" should be treated as stale unless it explicitly
means n≤7.

### Lesson

Finite exhaustive evidence must carry its finite quantifier into Lean. A theorem
with no `n≤7` hypothesis turns computational evidence into a false universal
axiom. Also: H=63 unlocks in the simplest possible OCF shape, Ω=K31, so the
old disconnected-factor obstruction was measuring the wrong graph shape.

---

## MISTAKE-051: Universal TRRT Revived Despite THM-025 Counterexample

**Date discovered:** 2026-05-29
**Found by:** opus-2026-05-29-S8 during repo scour
**Affects:** OPEN-Q-047, OPEN-Q-051/052/053 priority labels, INV-189/INV-186, HYP-1729

### What was assumed

Newer notes revived the Tournament Real-Rootedness Theorem (TRRT): for every
tournament T, I(Ω(T),x) has all real negative roots. The revived entries cited
small samples at n=9,10 with zero failures and treated the Hermite-Biehler
program as a route to a universal theorem.

### Why it was wrong

Canon THM-025 already disproves universal real-rootedness at n=9. The explicit
counterexample has score sequence [1,1,3,4,4,4,6,6,7] and
I(Ω,x)=1+94x+10x²+x³. Newton's inequality fails at k=2:
10² < (3/2)·94·1, so the polynomial has non-real roots.

### The correct framing

Real-rootedness is proved for n≤8 via claw-freeness and is common in samples,
but it is not universal. The right open problem is to characterize the
real-rooted subclass and locate the THM-025 failure inside any
Hermite-Biehler/interlacing framework.

### Impact

OPEN-Q-047 is retitled as a characterization problem. The HB lemmas are
downgraded from "critical to prove universal TRRT" to "important for the
real-rooted subclass program." HYP-1729 is marked REFUTED as a universal
theorem.

### Lesson

Sampling cannot override a canon counterexample. Before reviving a conjecture,
search `01-canon/theorems/` and `MISTAKES.md` for explicit disproofs.

---

## MISTAKE-052: THM-390 claimed twice in one day (codex-S547 vs codex-S548)

**Date discovered:** 2026-06-01
**Found by:** monad-reviewer-2026-06-01 (QC startup audit)
**Affects:** `01-canon/theorems/THM-390-*`, HYP-2036, HYP-2038, definitions.md,
TANGENTS.md, results/INDEX.md, hypotheses/INDEX.md, reflections, SESSION-LOG

### What happened

Two **distinct, both-PROVED** LRC theorems were independently filed under the same
id THM-390 on the same day:
- codex-2026-06-01-**S547** — `lrc-padic-zero-branch-cover-core` (committed fa44a9d):
  the denominator-sieve semantics (`z_q=0 ⇒ t=1/q` witness) and the minimum AP
  open cover `U_n={u: 2u≥n}` of size `floor(n/2)`.
- codex-2026-06-01-**S548** — `lrc-zero-branch-star-core-peeling` (committed 2264cf3):
  a single q-grid zero-branch star has empty strict endpoint-protection core, with
  explicit peel layers `|C|·m_s`.

S548 did not notice S547 had already taken THM-390 (concurrent sessions, both under
the `codex` line). The collision made every `THM-390` reference ambiguous — HYP-2036
in particular cited both theorems under the one number.

### Why it matters

Ambiguous canon ids break `depends_on` graphs and citations: a reader cannot tell
which theorem a reference means. This is the same class of issue as the
THM-013/THM-082 filename collisions (resolved as THM-012b / THM-084) and the
MISTAKE-018/018b renumber.

### Resolution

First claimant keeps the number. **S547 cover-core stays THM-390; S548 star-peeling
renumbered to THM-391.** File renamed, all star-pointing references updated
(definitions.md, TANGENTS, results/INDEX, hypotheses/INDEX, HYP-2036 [now depends on
both], HYP-2038, two reflections, the verifier script, SESSION-LOG entry). Both
proofs were independently re-derived and are correct (see verification notes in each
theorem file). Historical inbox/broadcast messages left as-is.

### Lesson

Before filing a new `THM-N`, run `ls 01-canon/theorems/ | grep THM-N` to confirm the
id is free — especially in concurrent multi-agent sessions where two agents may pick
the same "next" number on the same day. The repo still carries older unresolved id
collisions (THM-260×3, THM-338×2, THM-336/337 dups); those are latent debt that
should likewise be renumbered when next touched.

---

## MISTAKE-053: Systemic HYP-number collisions — five `HYP-N` reused in one 30-hour LRC burst

**Date discovered:** 2026-06-02
**Found by:** monad-reviewer-2026-06-02 (QC startup audit)
**Affects:** HYP-2050, HYP-2052, HYP-2058, HYP-2061, HYP-2063 (and their INDEX
entries, files, reflections). This is MISTAKE-052 (the THM-390 collision)
repeating at scale for the `HYP-*` namespace.

### What happened

Between 2026-06-01 and 2026-06-02, three concurrent agent lines (`opus`,
`oracle`, `codex`) ran the LRC@14/n=17 frontier in parallel and each picked the
same "next" HYP number within **3–12 minutes** of one another. Five collisions:

| HYP | First claimant (UTC) | Second claimant (UTC) | Both have a file? |
|-----|----------------------|------------------------|-------------------|
| 2050 | codex-S551 tetration 20:53 | oracle-S549o Lean 20:56 | only codex |
| 2052 | opus-S551 sieve-no-completeness 21:11 | oracle-S552 loneliness-spectral-gap 21:21 | **BOTH** |
| 2058 | oracle-S553o almost-lonely 15:03 | opus-S556 proof-lite-and-tension 15:21 | only opus |
| 2061 | oracle-S555o pinch-time-pigeonhole 17:41 | codex-S558 small-pinch-shield 17:54 | only codex |
| 2063 | opus-S559 2q-tight-tuple-apex 18:03 | codex-S559 n17-prime-gate 18:15 | **BOTH** |

### Why it matters

Same as MISTAKE-052: an ambiguous id breaks `depends_on`/citation graphs — a
reader cannot tell which hypothesis "HYP-2061" means. THM-396 already
`depends_on: HYP-2059, HYP-2060`, and HYP-2059's INDEX entry chains into HYP-2061,
so the ambiguity reaches a canon theorem's dependency closure.

### Resolution (this session)

- **HYP-2063 (both-file collision, newest):** fully renumbered. First claimant
  opus keeps `HYP-2063` (2q-apex); codex's n17-prime-gate → **HYP-2069**. File
  renamed, INDEX/SESSION-LOG/TANGENTS updated, 0 stray refs remain.
  **Caution — the frontier is a live race:** my first reassignment to `HYP-2064`
  immediately collided *four ways* — a rebase mid-session pulled in oracle-S557o
  (gap-bound), codex-S560 (gate-skip-transfer, has file), and monad-researcher-S560
  (A000568-Burnside), all independently filed under `HYP-2064` within hours. I moved
  codex-n17 clear of the contested 2050–2068 band to **HYP-2069**. The residual
  three-way `HYP-2064` (oracle-S557o / codex-S560 / monad-researcher-S560 — not my
  artifacts) is left to its owning sessions + the cleanup session, banner-flagged in
  the INDEX (suggest HYP-2070/2071 by first-commit timestamp). monad-researcher-S560
  already self-flagged it as a known 3-way collision in its SESSION-LOG entry.
- **HYP-2052 (both-file collision, older, 16 refs):** documented but **NOT yet
  renumbered** — the reference web is dense and a botched mass-rename would create
  more inconsistency than it removes. Canonical assignment: opus-S551
  `lrc-sieve-no-finite-completeness` is first claimant and keeps `HYP-2052`;
  oracle-S552 `lrc-loneliness-spectral-gap` is the duplicate and should be
  renumbered (suggested **HYP-2065**) in a focused future cleanup. Until then,
  always disambiguate by the file slug, not the bare number.
- **HYP-2050 / 2058 / 2061 (single-file collisions):** the idea that owns the file
  keeps the number (minimizes churn); the file-less duplicate (always an `oracle`
  index/reflection entry) is latent debt — suggested reassignments HYP-2066
  (oracle almost-lonely, ex-2058), HYP-2067 (oracle pinch-pigeonhole, ex-2061),
  HYP-2068 (oracle Lean-formalization, ex-2050). Disambiguate by slug meanwhile.

### Lesson

The MISTAKE-052 lesson ("`ls | grep` before filing") was logged for `THM-*` but
not adopted for `HYP-*`, and the failure rate is far higher because HYP numbers
advance many times per day across ≥3 concurrent lines. **Reserve the id first
(Step 5c checkpoint) before doing the work**, and `grep "HYP-N" 05-knowledge/hypotheses/INDEX.md`
+ `ls 05-knowledge/hypotheses/ | grep HYP-N` immediately before `finish_session`.
A sub-300-second reservation push at session start would have prevented all five.
Latent renumber debt remaining: HYP-2052 (both-file), and the three single-file
oracle duplicates above.

**Additional pre-existing two-file HYP collisions found in the same audit** (older
than this 24h window — full latent debt list for the future cleanup session):
- HYP-1969: `lrc-h-phase-plateau` vs `lrc-proof-route-currencies`
- HYP-1992: `lrc-n18-observer-source-gate-battlefield` vs `lrc-rapidity-formal-group-bridge`
- HYP-1995: `lrc-exact-gap-race-wall-ledger` vs `lrc-twin-roots-of-unity-bridge`
- HYP-2009: `lrc-polygon-outside-inside-arcs` vs `resonance-debt-conjecture`
- HYP-2040: `lrc-conditional-clearance-wedge-transitivity` vs `lrc-n4-measure-gap-unique-tight`

These confirm the collision rate has been chronic across the whole LRC era, not a
one-off. The cleanup session should resolve all of them by first-commit-timestamp
and rebuild a contiguous HYP index.

## MISTAKE (oracle-2026-06-03-S576o): pinch-M with a gcd(m,C)=1 filter gives SPURIOUS LRC counterexamples
When computing the loneliness radius M(S)=max_t min_i ||v_i t|| as a max over PINCH times
t=m/(v_a+v_b) (HYP-2075: the optimum is a pair-sum pinch), you MUST range over ALL
m=1..C-1, NOT only the coprime m (gcd(m,C)=1). The optimal pinch need not be in lowest
terms: e.g. S=(1,4,5) has M=1/3 attained at t=2/6 (pair-sum C=1+5=6, m=2, gcd=2), which a
coprimality filter drops, yielding a false M=2/9 < 1/4 -- a spurious "counterexample" to the
PROVEN LRC(4). Symptom: bounded-speed censuses report min M < 1/n at small n where LRC is a
theorem. Fix: drop the gcd filter (evaluate every m). Caught in lrc_even_ladder_selfconverse_proof_s576.py.

## MISTAKE (monad-compute-2026-06-03-S4): minH_strong(m)=m²−5m+9 is a 4-point coincidence; true value at m=7 is 25 not 23
HYP-2180 (opus-S599s) fit the strong-tournament Hamiltonian-path minimum minH_strong(m)=3,5,9,15 (m=3..6, exhaustive) to the quadratic m²−5m+9 and used a *near-transitive scan* to assert minH_strong(7)=23. EXHAUSTIVE enumeration of all 2^21 tournaments on 7 vertices (reversal-halving, `strong_H_spectrum_m7_exhaustive_monad_s4.py`) gives **minH_strong(7)=25, not 23** — and 23 is not a strong-tournament H-value at m=7 at all. The quadratic matched m=3..6 only by coincidence (same trap as MISTAKE-028/036: a pattern holding for 4 values then failing). ~~The CORRECT law is the known **Busch (2006) recurrence p(n)=p(n−1)+p(n−2)+1** for the minimum number of Hamiltonian paths in a strong tournament, giving 3,5,9,15,25,41,67,….~~ **[SUPERSEDED by MISTAKE-055, monad-compute-2026-06-06-S5/S6:** this recurrence is a MIS-CITATION of Busch. Exhaustive iso-class enumeration gives minH_strong(8)=**45** (not 41) and minH_strong(9)=**75** (not 67). The recurrence p(n)=p(n−1)+p(n−2)+1 itself fits only m≤7 then breaks at m=8 — the very same trap it was logged to correct. Busch's TRUE published sequence is **3,5,9,15,25,45,75,125,225,…**, which the exhaustive computation reproduces exactly.]** Everything else in HYP-2180 survived the exhaustive check: 7,21,63 are NOT strong values at m=7; 35=7·5 and 49=7² do fill in; only {7,21} are permanent H-gaps (63 achievable at n=8). Lesson (again): fit a candidate closed form only after it is verified at the FIRST genuinely new case, and trust a near-transitive scan for nothing more than a lower bound.

## MISTAKE-054: Incremental 3-cycle counter swapped i-beats-j / j-beats-i (under-pruning)

**Date discovered:** 2026-06-04 (monad-compute-2026-06-04-S4)
**Found by:** monad-compute, via ground-truth disagreement with the direct-count engine
**Affects:** the FIRST version of `h21_finite_check_v2_monad_s4.py` (the DFS-pruned
extension `extend()`); FIXED before any result was reported.

### What happened
The fast engine v2 builds each new vertex's orientation by DFS, accumulating the
new 3-cycles `{j, i, new}` incrementally and pruning when partial `c_3 > CAP`.
The triple's out-degrees were coded as
  `dj = ij + (1-nj)`,  `di = ji + (1-ni)`
i.e. vertex `j`'s out-degree used `ij` (i beats j) instead of `ji` (j beats i),
and symmetrically for `i`. Because the cycle test requires BOTH `di==1` and
`dj==1`, this is **not** a harmless relabel — with `nj`/`ni` attached to the
wrong term it tests a different condition, so some true 3-cycles were not counted.

### Symptom
v2 reported MORE iso-classes with `c_3<=10` than the direct-count engine v1
(m=7: 453 vs 339; m=9: 17,667 vs 2,575). Both engines still reproduced A000568
with the cap removed (the bug only affects the *capped* count), which is why the
no-cap self-validation did not catch it.

### The fix / correct framing
For triple `{j, i, new}` with `j<i`:
  `dj = ji + (1-nj)`  (j beats i? + j beats new?),
  `di = ij + (1-ni)`  (i beats j? + i beats new?),
  `dn = ni + nj`.
3-cycle iff `dj==di==dn==1`. After the fix, v2 matches v1 EXACTLY for m<=10 and
runs ~10x faster.

### Lesson
A no-cap / total-count self-check does NOT validate threshold/pruning logic.
Always cross-check a fast pruned engine against a slow direct-count engine on the
ACTUAL filtered quantity (here `#{iso classes with c_3<=10}`), not just the total.

## MISTAKE-055: Busch (2006) strong-min recurrence mis-cited as p(n)=p(n−1)+p(n−2)+1 (gives 41,67); true minH_strong is 3,5,9,15,25,45,75

**Date discovered:** 2026-06-06 (monad-compute-2026-06-06-S5/S6)
**Found by:** monad-compute, via exhaustive iso-class enumeration at m=8 and m=9
**Affects:** the MISTAKE-(2026-06-03-S4) entry above; HYP-2180; HYP-2271's "Busch-type" reduction; opus-S699j/k's strong-min(8)≤45 search bound; any downstream use of the 41/67 values.

### What was assumed
The prior monad-compute session (2026-06-03-S4), while correcting an *earlier* bad fit (the quadratic m²−5m+9), asserted that the minimum number of Hamiltonian paths in a strong tournament obeys the recurrence p(n)=p(n−1)+p(n−2)+1, giving 3,5,9,15,25,**41**,**67**,… and attributed this to Busch (2006).

### Why it was wrong
That recurrence matches the EXHAUSTIVE values 3,5,9,15,25 (m=3..7) but BREAKS at m=8 — the identical "holds for several values then fails" trap (cf. MISTAKE-028/036/054) the entry was written to warn against. EXHAUSTIVE enumeration of all non-isomorphic strong tournaments (generated by canonical augmentation, validated against A000568 = …,456,6880,191536 for n=7,8,9) gives:

  minH_strong(m) = 3, 5, 9, 15, 25, **45**, **75**   for m = 3..9   (NOT …25,41,67)

opus-S699j/k's non-exhaustive reversal-search bound strong-min(8) ≤ 45 was therefore TIGHT (=45), not loose; and strong-min(9)=75.

### The correct framing
Busch, "A Note on the Number of Hamiltonian Paths in Strong Tournaments" (Electron. J. Combin. 13 (2006), #N3) proves the minimum equals Moon's (1972) upper bound, with sequence **3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625, …** (n≥3). The exhaustive computation reproduces this EXACTLY through m=9. Empirically the data satisfies p(n)=3·p(n−2) for every step except n=7 (25 vs 27); the asymptotic growth is ~(√3)^n. (Do NOT re-fit a clean recurrence here without checking against Busch's closed form.)

### Impact — POSITIVE for the program
- HYP-2271 (opus-S699j/k) reduced the delta-field polarization / "7,21 never H" theorem to the lower bound **strong-min(m) ≥ 22 for all m≥7**. Busch (2006) proves the minimum is 25,45,75,… (strictly increasing, ≥25 for m≥7) FOR ALL n ⟹ the reduction is CLOSED by a published theorem, not just "Busch-type, to be proven".
- {7,21} are confirmed absent from the strong H-spectrum exhaustively for m≤9 (7,21,35 below strong-min; 49,63 ARE strong values at m=8). Combined with strong-component multiplicativity H=∏H(Cᵢ), this verifies the phantom-volume theorem (only {7,21} are durable forbidden H, genus-2 multiplicative semigroup) for all tournaments whose strong components have ≤9 vertices, and reduces the all-n statement to the published Busch bound.

### Lesson
When citing a literature recurrence, verify its VALUES against the first genuinely new exhaustive case before propagating it. The "41/67" recurrence was adopted as the fix for a bad fit and itself silently inherited the same coincidence failure mode. Exhaustive iso-class enumeration (via gentourng/nauty-style canonical augmentation — here a pure-Python canonical-augmentation generator validated by A000568) makes m=8,9 cheap (6880 / 191536 classes) where labeled enumeration (2^28 / 2^36) is not.

---

## MISTAKE-056: Signed-LRC worry-set "split" claimed first at n=14 — it is first at n=8

**Date discovered:** 2026-06-06
**Found by:** monad-explorer-2026-06-06-S708b
**Affects:** opus-S699 reflection `signed-lrc-theory-sign-is-a-cut-and-the-worryset-splits-s699.md`, HYP-2262 (the "MAIN RESULT" narrative), and the broadcast MSG-001 ("n=14 is the FIRST n whose C=2n−1 admits a doubled-speed shell-partner"). Does NOT affect the theorems T1–T4 (all correct).

### What was claimed
"Through n=7 every tight (M=1/n) config is shell-partner-free; it FAILS at n=14 (V*=(1..11,13,24), shell-partner 3+24=27). n=14 is the FIRST n whose C=2n−1 admits a doubled-speed shell-partner (24=2·12)."

### Why it is wrong
S699 verified n=4,5,6,7 (shell-partner-free) and then jumped straight to the *known* n=14 frontier, never checking n=8,10,12. But **n=8 already carries shell-partner tight configs.** Exhaustively (exact M, and independently the S592 floor test), the n=8 worry-set has 3 floor-tight primitive configs and **two carry a shell-partner**:
- `(1,2,3,4,5,7,12)` = AP_8 with 6→12, where 12=2·6≡−3 (mod 15), shell-partner (3,12), 3+12=15=2·8−1. M=1/8.
- `(1,4,5,6,7,11,13)`, shell-partner (4,11). M=1/8.

The first is the SAME "double the (n−2) speed" mechanism as n=14's V* (double 12→24). So n=8 (C=15=3·5) is the first n whose C admits a doubled-speed shell-partner tight config — not n=14.

### The correct framing
"tight ⟹ no shell-partner" holds for n≤7 and FAILS first at **n=8**. The V*-type (shell-partner-carrying tight) configs form the infinite **Family II** = AP_n with (n−2)↦2(n−2), floor-tight ⟺ **n≡2 (mod 6)** = every even n with 3∣(2n−1) = {8,14,20,26,…} (verified exact n≤29). n=14 is special only as the first such n whose C is a pure prime power (3³). The shell-partner is always (3, 2(n−2)). See HYP-2281 / reflection `the-worryset-split-is-at-n8-shell-transversality-as-the-gauge-invariant-s708.md`.

### Impact
- The "split exists / is finer than M" conclusion STANDS — only the "first n" is corrected (8, not 14).
- POSITIVE: gives a minimal, SOLVED (LRC(8) is true) laboratory for the prime-2×prime-3 doubling mechanism that recurs unsolved at n=14; the (3,24) carry attack should be prototyped on n=8's (3,12).
- Also reframes the carrier as a purely UNSIGNED, gauge-invariant property: "carries a shell-partner" ⟺ "S mod 2n−1 is not a shell-transversal" (HYP-2281 L1–L2).

### Lesson
When a property is verified up to n=k and then claimed to "first fail" at some larger known-frontier n=N, CHECK every n in (k,N). The interesting frontier (here n=14, C=3³) is rarely the *first* instance of a phenomenon; the first instance (n=8, C=3·5) is usually smaller, more tractable, and already solved.

## MISTAKE-057: THM-427 + HYP-2294 + T765 triple-claimed by two concurrent monad-explorer-S3 instances

**Date discovered:** 2026-06-06 (monad-explorer-2026-06-06-S3, the gcd-torsion lane, at close-out)
**Found by:** monad-explorer (self, on post-push `ls` of theorem dir)
**Affects:** `01-canon/theorems/THM-427-*`, HYP-2294 (INDEX), T765 (TANGENTS), and the two-tower reflection/script. Same class as MISTAKE-052 (THM-390 dup) / MISTAKE-053 (HYP-* dups).

### What happened
Two DISTINCT, both-good LRC results — both responding to the same dispatched seed's "find the unifying statement" — were filed by two concurrent `monad-explorer-2026-06-06-S3` instances under the SAME ids THM-427 / HYP-2294 / T765:
- **gcd-torsion lane** (commit 63ed166, 2026-06-07 01:38:09 UTC): composite-LRC cell-leak `= N_i·n − g·W_i(g)`, a function of `gcd(r,n)=n/ord(r)`.
- **two-tower lane** (commit dba3832, 2026-06-07 01:46:44 UTC): the clock ℤ/n × shell ℤ/(2n−1) coprime-CRT witness group.

The two-tower commit landed ~8.5 min later, when the gcd-torsion THM-427 was already on origin — it did not rebase-detect the taken id (the live-race failure mode of MISTAKE-053).

### Resolution
First claimant keeps the number (gcd-torsion, earlier commit + already on origin). The two-tower lane renumbered: **THM-427→THM-428, HYP-2294→HYP-2295, T765→T766**. Theorem file `git mv`'d; self-references flipped in the two-tower theorem file, its reflection, its script+`.out`, and the shared INDEX table-row / TANGENTS entry. 0 stray refs remain (the two results are complementary: gcd-torsion = mod-n leak face, two-tower = the mod-n ⟂ mod-2n−1 CRT product — they reinforce, not conflict).

### Lesson
The MISTAKE-053 fix ("reserve the id at Step 5c BEFORE the work; `ls 01-canon/theorems | grep THM-N` immediately before finish") still was not adopted. Sub-300s reservation pushes at session start would have prevented this. When two agents share a machine-name line (`monad-explorer`) and a date, the `[machine]-[date]-S[N]` id does NOT disambiguate concurrent instances — both became "S3". Consider a per-instance random suffix when a line is run in parallel.
---

## MISTAKE-058: a THIRD concurrent monad-explorer-S3 lane (signed-pairwise) also hit THM-427/HYP — the collision was 3-way, not 2-way

**Date discovered:** 2026-06-06 (monad-explorer-2026-06-06-S3, the **signed-pairwise** lane)
**Found by:** monad-explorer (this session), at session-end rebase — after MISTAKE-057 (the two-tower
lane) had already documented the *gcd-torsion vs two-tower* pair.
**Affects:** completes MISTAKE-057. The same window saw a THIRD distinct LRC result claim `THM-427`.

### What happened
MISTAKE-057 recorded TWO concurrent `monad-explorer-2026-06-06-S3` lanes colliding on `THM-427`. There
was a **THIRD**: this signed-pairwise lane (`THM-427-signed-pairwise-floor-is-a-maxcut-LRC`,
`Gstar ≥ 1/(2 r_min)`), committed 20:51:32 -0500 — after gcd-torsion (20:37) and two-tower (20:46).
Three distinct, all-good LRC theorems under one id `THM-427`, all from the same instance name.

### Resolution (first-come keeps the number; consistent with MISTAKE-057)
- `THM-427` → **gcd-torsion** (first). `HYP-2294`, `T765` → gcd-torsion.
- `THM-428`, `HYP-2295`, `T766` → **two-tower** (second; self-renumbered, MISTAKE-057).
- `THM-429`, **`HYP-2296`**, T764-update → **signed-pairwise** (third, this lane): file `git mv`'d,
  id + `signed_lrc_rmin_bound_monad_s3.py` docstring updated; HYP renumbered 2295→**2296** (2295 is
  the two-tower's), references flipped in THM-429, the reflection, INDEX, TANGENTS, SESSION-LOG. My
  already-pushed *commit messages* still say THM-427/HYP-2295 (immutable history); the canon files are
  **THM-429 / HYP-2296**.

### Lesson
Even after a collision is "resolved," re-check before finishing: a 2-way resolution can be incomplete
if a third concurrent instance is in flight. And renumber by first-commit author-date end-to-end
(here gcd-torsion < two-tower < signed-pairwise ⟹ 427/428/429, 2294/2295/2296). The deeper fix
remains MISTAKE-053's: reserve ids at Step 5c before doing the work; when a `[machine]-[date]` line is
run ≥3-way in parallel, the `S[N]` suffix does not disambiguate — use a per-instance random tag.
The Step-5c "reserve the id first, `ls | grep` before filing" rule (MISTAKE-053) must run **even
against your own instance id** — concurrency can duplicate the *session name*, not just the number.
A one-line reservation push at session start (claiming THM-N/HYP-N as honest stubs) would have
prevented all three. When three files share `THM-N`, resolve by first-commit author-date, not by
who notices last.

**ADDENDUM (same session, on rebase):** it was a THREE-way race, not two. A *third* concurrent
`monad-explorer-2026-06-06-S3` filed THM-427 = "signed pairwise floor is a max-cut LRC"
(commit 20:58 -0500, 20 min after the gcd-torsion claim; it also forward-referenced HYP-2294 for
its asymptotic question). First-claimant rule again: gcd-torsion keeps THM-427; the max-cut lane →
**THM-429**, its HYP-2294 forward-ref → **HYP-2296** (free). Three independent S3 instances,
three THM-427 claims, all within 20 minutes — the strongest evidence yet for per-instance id
suffixes when one agent line is fanned out in parallel.

---

## MISTAKE-059: "Exactly 3-to-1" inferred from a count ratio without checking the map (caught + corrected same session)

**Date discovered:** 2026-06-07 (monad-explorer-S6, self-caught)
**Affects:** THM-436 ADDENDUM (2″) as first checkpointed; HYP-2305; reflection `the-icosahedral-fifteen-s6.md` (all corrected before any agent built on them)

### What was assumed
The commutator map {60 oriented overlapping cyclic-triangle pairs on a 5-set} → {20 three-cycles of A_5}
was stated as "**onto and exactly 3-to-1**" — inferred purely from `60 / 20 = 3`, and dressed up as the
icosahedral **face-vertex flag** incidence (`20` faces × `3` vertices = `60` flags, flag→face uniformly
3-to-1).

### Why it was wrong
The fibers are **not uniform**. Direct enumeration (`04-computation/icosahedral_flag_fibers_monad_s6.py`)
gives fiber sizes `{2 (×3), 3 (×14), 4 (×3)}` (sum 60 over 20 three-cycles). The `3`-to-`1` holds only
on average. So `60 = |A₅|` is the group order, NOT a flag count, and the commutator covering is NOT the
icosahedral flag map.

### The correct framing
What is actually true and robust: **every one of the 60 oriented overlapping pairs has a commutator of
cycle-type a 3-cycle** (conjugation/inversion-invariant ⇒ order-convention-independent), and the 60
commutators are **onto all 20** three-cycles. That type-uniformity — not any multiplicity-uniformity — is
the real content of "A_n perfect realized by overlapping triangle pairs."

### Lesson
A matching TOTAL (`60 = 20·3`) does not certify a uniform MAP. When a count coincides "too cleanly" with
a known structure (here, icosahedral flags), verify the **fibers / the map**, not just the cardinality.
This is the project's own "too clean ⇒ test it" rule applied to itself; the test refuted the clean story
and left the honest one (type-uniformity + a 15-fold canonical bijection) standing.

---

## MISTAKE-060: THM-438 "bigon-trees ARE the Catalan count" — top order is a SIGNED cactus cancellation, not a +1-per-tree count

**Date discovered:** 2026-06-07 (monad-explorer-2026-06-07, deep-research / analytic lane, 3rd session)
**Found by:** monad-explorer, while attempting the "small remaining write-up" (the +C_k sign) flagged OPEN in THM-438's Honest-status section
**Affects:** THM-438 Part B proof MECHANISM + error term; the reflection `the-paley-cluster-integrals-are-catalan-numbers-tree-walks-and-the-moment-method.md` ("Patterns with any non-bigon cycle ... are O(p^{k+1/2})"; "the top order is an all-bigon graph ... a tree of bigons ... counted by C_k"); the reflection's stated O(1/sqrt p) convergence. **Does NOT affect** the STATEMENTS A_{2k}=C_k p^{k+1} or R(p)->e — both stand (verified).

### What was assumed
The leading order p^{k+1} of the cluster integral `A_{2k} = sum_{distinct x_0..x_{2k}} prod chi(x_{i+1}-x_i)`
is reached ONLY by all-bigon coincidence patterns; a **tree** of k bigons maximizes V=k+1; each such bigon-tree
(= Euler tour of a plane tree) contributes **+1**, so the leading coefficient is literally the Catalan count
C_k. "Patterns with any non-bigon cycle ... are O(p^{k+1/2})." Error term O(p^{k+1/2}).

### Why it was wrong (verified exactly, `04-computation/paley_cluster_cactus_census_monad.py`)
Three things are false:
1. **Bigon-trees do NOT each contribute +1.** Via the partition-lattice Moebius inversion
   `A_{2k} = sum_sigma mu(0,sigma) M_sigma`, a bigon-tree pattern `sigma` carries Moebius weight
   `mu(0,sigma) = prod_blocks (-1)^{|B|-1}(|B|-1)!`, which is NOT 1 when a vertex is visited >=3 times.
   The bigon-tree leading coefficient (sum over non-crossing edge-pairings of `prod_v (b_v-1)!`) is
   **1, 3, 13, 69, 421, 2867** (k=1..6) = **OEIS A088368** (g.f. `A=sum n! x^n A^n`, `a(n)~e*n!`) —
   NEITHER C_k NOR (2k-1)!!. At k=2 bigon-trees give **3**, at k=3 they give **13** (census confirms).
2. **Even cycles DO reach the top order p^{k+1}.** The single 2k-cycle pattern (`x_0=x_{2k}`) equals
   `tr(M^{2k}) = (-p)^k(p-1) ~ (-1)^k p^{k+1}` — the SAME order as bigon-trees, not O(p^{k+1/2}).
   It enters with `mu=-1` and SUBTRACTS. More generally every **even cactus** (connected graph whose
   biconnected blocks are all even cycles, total half-edges k) contributes at p^{k+1}.
3. **The Catalan number is a signed cancellation.** Census:
   `k=2: bigons(+3) + 4-cycle(-1) = 2 = C_2`;  `k=3: bigons(+13) + {bigon+4cyc} + {6cyc} = 5 = C_3`.

### The correct framing
**Closed form (PROVED via Gauss-sum inversion `chi(w)=g^{-1} sum_t chi(t) omega^{tw}`, verified exactly):**
```
M_sigma = (-1)^k * p^{V-k} * F(sigma),    F(sigma) = sum over F_p-flows t on G_sigma of prod_e chi(t_e),
```
V = #blocks, flow space = cycle space (dim m = 2k-V+1). A pattern reaches p^{k+1} iff F reaches full order
p^m; those are exactly the **even cacti**. The leading coefficient of A_{2k} is the **signed sum over even
cacti** `sum mu(0,sigma) * lead(M_sigma) = C_k` — an inclusion-exclusion that converts the all-pairings
overcount (A088368, ~e*n!) into the **non-crossing** count C_k. This is the genuine free-probability /
moment-method content: the two-point Gauss spectrum's even-cycle terms are PRECISELY the corrections that
take Gaussian-style pairings to the semicircle's non-crossing pairings.

**Error term:** `A_{2k} = C_k p^{k+1} + O(p^k)` (NOT O(p^{k+1/2})). Verified: `(A_4-2p^3)/p^2` is STABLE
(~ -7.1..-7.8 -> ~-8), while `/p^{2.5}` drifts to 0. Hence R(p)-e has relative correction **O(1/p)**,
resolving the reflection's stated O(1/sqrt p) vs the close-out's "favors 1/p" IN FAVOR OF 1/p.

**Part C simplifies:** R(p)->e needs **NO Weil bound**. The only V=2k no-leaf pattern is the single
2k-cycle = `tr(M^{2k}) = (-p)^k(p-1)` (elementary); V<2k is trivially `O(p^{2k-1})=o(p^{2k})`.

### Impact
- THM-438 Part B mechanism CORRECTED (addendum added). Statements A_{2k}=C_k p^{k+1}, R(p)->e UNCHANGED.
- Part C upgraded: fully elementary (no Weil).
- Error term corrected p^{k+1/2} -> p^k; convergence rate of R->e pinned to 1/p (feeds HYP-2307 #2).

### Lesson
The project's own rule again: a clean final count (C_k) reached by a clean-sounding mechanism
("bigon-trees = Catalan") does not certify the mechanism. The Moebius weights and the
equal-order even-cycle patterns were invisible at the level of "count the bigon-trees." Always
decompose the inclusion-exclusion and check which patterns share the leading order — here the
cancellation `A088368 -> C_k` is the actual phenomenon, and it is the free-probability fingerprint
the moment-method slogan was pointing at.

---

## MISTAKE-061: THM-438 — the top-order patterns are NOT "even cacti"; they are the larger class of EVEN-SERIES patterns (even theta graphs included)

**Date discovered:** 2026-06-07
**Found by:** monad-explorer-2026-06-07 (deep-research / analytic lane, 4th session)
**Affects:** THM-438 ADDENDUM and MISTAKE-060 (the *characterization* of which coincidence
patterns reach the leading order `p^{k+1}`). Does NOT affect the Catalan law `A_{2k}=C_k p^{k+1}`
itself (re-confirmed here, rigorously) nor `R(p)->e`.

### What was assumed (MISTAKE-060 / THM-438 ADDENDUM)
"`M_sigma` reaches the top order `p^{k+1}` **iff** `F(sigma)` reaches full order `p^m` — exactly
the **even cacti** (connected, all biconnected blocks even cycles)." The census then grouped the
leading coefficient as bigon-trees (+A088368) corrected by even-cycle **cacti** down to `C_k`.

### Why it was wrong
`F(sigma) = sum_{flows} prod_e chi(t_e)` reaches full order `p^m` iff the flow-form product
`P(s) = prod_e ell_e(s)` is a **perfect square** (then `chi(P)=chi(Q^2)=+1` off the zero locus,
so `F ~ +p^m`). `P` is a perfect square iff **every series-class of edges has even size** (each
distinct flow-line occurs an even number of times). The even cacti satisfy this — but so do
**even theta graphs** (two vertices joined by three even paths; biconnected block is NOT a single
cycle) and, generally, all "even series-parallel" 2-connected patterns. These are NON-cacti yet
reach `p^{k+1}` and MUST be counted. Verified (`04-computation/paley_cluster_theta_check_monad.py`):
at `k=3` the `V=5, m=2` top-order patterns are **6 even-cacti{2,4} + 1 even-theta(2,2,2)** — the
even theta (mu=+1) was invisible to the "even cacti" census (it sat in the `(6,)` biconnected
bucket, silently cancelling the single 6-cycle, so the *total* still came out right).

### The correct framing (VERIFIED k<=4; the `g` step PROVED)
Let `c0 = lim A_{2k}/p^{k+1}`. Then
```
c0 = (-1)^k * sum_{rho : connected, EVERY series-class even}  mu(0,rho) * g(rho),
```
and `g(rho) := lim F(rho)/p^m = +1` for EVERY such pattern. **`g==+1` is PROVABLE:** within each
series-class the closed Euler walk passes straight through the degree-2 internal vertices, so all
edges of the class get the SAME orientation sign `s in {+1,-1}`; the class is even, so
`prod_{e in class} sign_e = s^{even} = +1`; hence `P = (prod sign_e) Q^2 = +Q^2` and
`g=chi(P)=+1`. Therefore the entire character/Gauss-sum content collapses and
```
($$)   sum_{rho : even-series pattern}  mu(0,rho)  =  (-1)^k C_k        (number-theory-FREE).
```
RIGOROUSLY CONFIRMED `c0 = 2, 5, 14 = C_2, C_3, C_4` by clean Richardson (`1/p`) extrapolation of
the exact flow-Moebius value (`04-computation/paley_cluster_topterm_monad.py`) — this also REPLACES
the prior slowly-converging census (which read `1.56, 2.77, 3.11` at `p<=23` and only *looked* like
it might reach `5`). The breakdown:
```
k=3:  bigon-trees(m=3) +13,  (m=2: cacti+theta) -9,  (m=1: 6-cycle) +1   = 5 = C_3
k=4:  bigon-trees(m=4) +69,  (m=3) -72,  (m=2) +18,  (m=1: 8-cycle) -1   = 14 = C_4
```
(bigon-tree sub-sums `+13, +69` = OEIS A088368, the all-pairings overcount, as before).

### Impact
- THM-438 ADDENDUM-2 added: Catalan law `A_{2k}=C_k p^{k+1}` RE-CONFIRMED (rigorous, k<=4), error
  `O(p^k)` unchanged, `R(p)->e` unchanged.
- The MECHANISM is corrected a SECOND time: top-order class = **even-series patterns** (perfect-square
  flow product), strictly larger than even cacti. `g==+1` is proved, reducing handoff #1 to the
  clean number-theory-free Moebius identity `($$)`.
- Free-probability reading SHARPENED: the random skew-Rademacher matrix gives `C_k` *directly* from
  non-crossing pairings (each `+1`, no factorials); the deterministic Paley Moebius expansion
  over-counts to A088368 in the bigon sector and the even cacti + even thetas + ... cancel it back to
  `C_k`. The equality is Wigner quasirandomness; `($$)` is its exact combinatorial fingerprint.

### Lesson
MISTAKE-060 corrected the *value* mechanism (bigon-trees -> A088368 -> C_k) but inherited a wrong
*support*: "even cacti." A pattern can saturate the flow character-sum without being a cactus — any
even series-parallel skeleton does. When a leading-order census gives the right TOTAL, that does not
certify the per-class STRUCTURE: a missing pattern (the even theta) can hide inside a coarse bucket,
cancelling against another, leaving the total correct and the story wrong. Characterize the support
by the actual saturation condition (perfect-square flow product = even series-classes), not by the
most familiar sub-family.

---

## MISTAKE-062: even-series pattern count is NOT OEIS A215257 — a 5-term coincidence that breaks at k=6

**Date discovered:** 2026-06-07 (monad-explorer, 8th session)
**Found by:** monad-explorer-2026-06-07 (deep-research, 8th session)
**Affects:** THM-438 ADDENDUM-3 point (2); HYP-2308; the reflection `the-drt-engine-is-S-squared-equals-J-minus-nI-the-catalan-is-genus-zero.md`; INDEX/SESSION-LOG entries asserting "even-series count = A215257"

### What was assumed
THM-438 ADDENDUM-3 (5th session) identified the number of EVEN-SERIES patterns of the
path `[0..2k]` (= the unsigned support of `(**)`) as **OEIS A215257**: the values for
`k=1..5` are `1, 3, 13, 67, 383 = A215257(k+1)` (indecomposable deque-sortable
permutations). The recursion script hardcoded the *predicted* next value `A215257(7)=2345`
for `k=6` but NEVER actually computed `k=6` (its `KMAX` default was 5).

### Why it was wrong
A direct exhaustive count at `k=6` (fast integer enumerator
`04-computation/paley_starstar_triangle_fast_monad.py`, cross-validated against the original
SVD test `04-computation/paley_starstar_crosscheck_monad.py` with **0 disagreements** over
all `Bell(13)=27.6M` partitions exhaustively at `k<=5` and a 300k sample at `k=6`) gives
```
   even-series count, k=1..6  =  1, 3, 13, 67, 383, 2351.
```
The OEIS b-file gives `A215257(7) = 2345 != 2351`. An OEIS search for
`1,3,13,67,383,2351` returns **no results** — the unsigned even-series count is not (yet)
a catalogued sequence. The A215257 match was a **5-term small-number coincidence**.

### The correct framing
- The unsigned even-series pattern count is `1, 3, 13, 67, 383, 2351, ...` (computed,
  rigorous through k=6), NOT A215257, and presently matches no OEIS sequence.
- This does NOT touch any headline result. The Moebius-SIGNED sum
  `(**) S_k = sum_{even-series} mu(0,sigma) = (-1)^k C_k` is independently re-verified
  exhaustively at `k=6` (`S_6 = 132 = C_6`), as is the cycle-rank triangle row
  `t(6,m) = 1, 45, 560, 2626, 4845, 2867` and the loop equation
  `S_k = -sum_{i+j=k-1} S_i S_j`.
- If anything the refutation SHARPENS the thread's thesis: the *unsigned* count is so
  unstructured it is not even a known sequence, while the *signed* sum is the cleanest
  possible (Catalan). "The Catalan is a cancellation, not a count" is now literal.

### Impact
- THM-438 ADDENDUM-3 (2) corrected (see ADDENDUM-6). HYP-2308 / INDEX A215257 cells updated.
- The "indecomposable deque-sortable permutations" bijection program (ADDENDUM-3/4 handoff)
  is moot — there is no A215257 bijection to find because the counts differ.

### Lesson
A 5-term OEIS hit is weak evidence — A215257 and the even-series count share five terms by
chance. NEVER hardcode an "expected next" OEIS value as if computed; compute it. Generic
divergence of two integer sequences after a short common prefix is the default, not the
exception (cf. MISTAKE-006 ratio coincidence, MISTAKE-010 small-n pattern break).

---

## MISTAKE-063: THM-438 ADD-9 wrongly "refuted" `A088368(m) ~ e*m!` — the original asymptotic is CORRECT (Kotesovec/OEIS); ADD-9 sampled only the pre-peak rising side

**Date discovered:** 2026-06-07 (monad-explorer, 12th session)
**Found by:** monad-explorer-2026-06-07 (deep-research, 12th session), via direct OEIS lookup of A088368
**Affects:** THM-438 ADDENDUM-9 point (6) ("CORRECTION FLAG"); any reflection/INDEX/SESSION-LOG
line asserting "A088368(m) ~ m!(m+2)/2, NOT e*m!"

### What was assumed (the erroneous "correction")
ADD-9 point (6) flagged the long-standing claim "`A088368(m) ~ e*m!`" (the diagonal
`t(k,k)` of the cycle-rank triangle = the all-pairings overcount) as **NOT supported by the
data**, observing `A088368(m)/m! = 1, 1.5, 2.17, 2.875, 3.51, 3.98, 4.45` (m=1..7) is
"monotonically increasing PAST e ≈ 2.718", and proposed instead the empirical
`A088368(m) ≈ m!(m+2)/2` (ratio `(m+2)/2`).

### Why it was wrong
The asymptotic `a(n) ~ e * n!` is an **established OEIS result** for A088368 (Vaclav
Kotesovec, Apr 10 2019; verbatim on the A088368 page). Computing the ratio with the
**true** b-file values (ADD-9 also had a transcription slip: `A088368(7) = 21477`, not
"22417") shows `a(n)/n!` does NOT increase monotonically — it OVERSHOOTS e, **peaks at n=8
(≈ 4.359)**, then **strictly DECREASES** back toward e:
```
   n:        2     3     4     5     6     7     8(peak) 9    10    12    16    20
   a(n)/n!: 1.50  2.17  2.88  3.51  3.98  4.26  4.36   4.32  4.19  3.85  3.36  3.14  ... -> e
```
ADD-9 sampled only `m<=7` — entirely on the **rising side, before the peak** — and mistook
the slow, overshooting approach for divergence. ADD-9's `(m+2)/2` fits the rising side only
and diverges (it predicts 11.0 at n=20, where the true ratio is 3.14 and falling).
Verified: `04-computation/paley_starstar_diagonal_noncrossing_monad.py`.

### The correct framing
- `A088368(m) ~ e * m!` is CORRECT. The diagonal of the cycle-rank triangle grows like
  `e * m!` (this is the "wild end" of the bridge polynomial; equivalently `h_m(m) =
  A088368(m)/m! -> e`).
- A088368 = **"number of partitions of [n] into sets of NONCROSSING LISTS"** (Callan,
  arXiv:0711.4841), G.f. `A(x/F(x)) = F(x)` with `F(x) = sum n! x^n`. It is a named,
  closed-form object — the diagonal is NOT "uncatalogued" (only the OFF-diagonal columns are).
- The original `~e*m!` slogan (ADD, ADD-8) should be RESTORED; ADD-9 point (6) is retracted.

### Impact
- THM-438 ADD-9 point (6) retracted; ADDENDUM-10 records the restoration + the noncrossing-
  lists identity. No headline result was ever affected (`A_{2k}=C_k p^{k+1}`, `R(p)->e`,
  `(**)`, column rationality all stand).
- The "tame<->wild bridge" of ADD-8/9 now has explicit asymptotic endpoints: `h_m(m) -> e`
  (wild/A088368) and `h_m(-1) -> 0` (tame/Mersenne, super-exponentially).

### Lesson
The MIRROR of MISTAKE-062. There, a sequence MATCH was over-trusted from 5 small terms;
here, an asymptotic was REFUTED from 6 small terms. **A factorial-scale ratio that is still
changing at n<=8 has converged to nothing** — slow asymptotics routinely overshoot and turn
around (here the turn is at n=8). Before declaring an `~ c*n!` claim false from a finite
ratio table, (i) check OEIS for an established asymptotic, and (ii) extend the ratio far
enough to see whether it is still rising — a monotone prefix is not a limit.

---

## MISTAKE-064: misread Erdős Problem 64 as "even cycle" (2k) when it is "power-of-2 cycle" (2^k); proved the trivial settled statement and wrongly framed it as the open problem

**Author:** opus-2026-06-08-S708 (caught by the user same day)
**Severity:** framing/attribution error (the math proved is correct but answers the WRONG question)

### What happened
The user asked to "Work On: Does every finite graph with minimum degree at least 3 contain a
cycle of length 2𝑘 for some 𝑘 ≥2?" — **Erdős Problem 64 = the Erdős–Gyárfás conjecture (1995):
every graph with min degree ≥3 contains a cycle whose length is a POWER OF TWO** (`2^k`: 4, 8,
16, …). This is **OPEN and falsifiable.** I read "2𝑘" as "2·k" (an even number ≥4), proved the
classical longest-path result *"min degree ≥3 ⟹ an even cycle of length ≥4"* (TRUE and settled),
and **wrongly presented it as the user's problem** in THM-443/HYP-2313/S708.

### Why it is wrong
"Min degree ≥3 ⟹ even cycle" is settled-true, so it CANNOT be an open falsifiable problem — that
alone should have flagged the misreading. The real condition is **multiplicative/2-adic** (length
exactly a power of 2), enormously stronger than "even." Petersen (girth 5) has no 4-cycle but an
8-cycle; a counterexample would be a min-degree-3 graph avoiding 4, 8, 16, 32, … *all at once*.

### Fix
- THM-443 **rescoped**: it correctly proves the EVEN-cycle lemma (classical) — relabeled to state
  explicitly that this is NOT Erdős 64. The power-of-2 problem is treated honestly as OPEN in
  THM-444/S709 (verified computationally on small/structured graphs; no counterexample; proven for
  cubic planar by Heckman–Krakovski; Markström computational searches).
- HYP-2313's parity-covering lens stands as a lens but its even-cycle leg is downgraded to the
  weaker statement; the power-of-2 ("dyadic") version is the open one.
- The S707 Pfaffian/even-dicycle bridge applies to EVEN cycles (Pólya), NOT to power-of-2 cycles.

### Lesson
**When a problem is stated to be OPEN, first check that your proof does not settle it — if it
does, you have misread the problem.** A one-character misread (`2^k` vs `2k`) flips a trivial
exercise into a famous open conjecture. Parse the *difficulty claim* as a constraint on the
interpretation.

---

## MISTAKE-065: Tile-Bit Negation Under Path Reversal — T^op Along a Reversed Path is a Grid Transpose WITHOUT Negation

**Date discovered:** 2026-06-09 (caught in-session by branch-C verification before propagating)
**Found by:** kind-pasteur-2026-06-09-S1 (hand derivation wrong; computational branch corrected it)
**Affects:** THM-447(5) original wording, HYP-2335 original wording, T767 original wording (all corrected in place)

### What was assumed
In the canonical frame of the skew-Sylvester double D(T) (path, twin, reversed path), the copy-2
tile block was claimed to be the "grid-transposed NEGATED" copy of T's tiling x — reasoning:
copy 2 is T^op (all arcs reversed), so its tile bits must be complemented.

### Why it was wrong
Copy 2 is traversed along the REVERSED base path. Reversing the path also reverses the
upper/lower convention of every tile, which complements each bit a second time. The two
negations cancel exactly: t(X,Y) = x(σ_n(X,Y)) — grid transpose with NO negation. Verified in
100% of 1098 framed + 1096 labeled cases (n=3..6).

### The correct framing
The single negated Sylvester copy in the tile schema [[H,H],[H,−H]] lives in the CROSS block:
σ-partner cross tiles (i,j) ↔ (j,i) carry complementary bits A[i][j] vs A[j][i]. See THM-452.

### Lesson (same family as MISTAKE-033)
Arc reversal (T^op), path reversal, and tile-bit complementation are THREE involutions that
compose in non-obvious ways. T^op + reversed path = grid transpose, NOT complement-transpose.
Whenever a claim involves "negated copy" at the tiling level, track ALL involutions explicitly
— two of them silently composing to the identity is the recurring trap. (See also the broader
NUMBER-conflation hygiene reference [[polysemous-constants-bridges-traps-and-homonyms]] (klein-S7):
the "2" homonym above is one of several constants — 2, n, 7, 14, 28, 6, 3 — that wear arithmetic vs
dimensional hats; run the PERSISTENCE TEST before treating any numeric coincidence as structure.)

### Impact
None propagated: corrected in THM-447(5-CORRECTED), THM-452(1), HYP-2335 status, T767 note,
all within the same session.
## MISTAKE-066: Erdős 592 finite-bridge direction stated BACKWARDS in first tree-grid script

**Date:** 2026-06-09
**Found by:** mac-mini-2026-06-09-S1 (same session, caught while writing the pysat version)
**Affects:** `04-computation/erdos592_treegrid_dichotomy_macmini_s1.py` original docstring
(corrected in place); draft doc §3.2 v1 (corrected)

### What was assumed
"An infinite witness for ω^n ↛ (ω^n,3) restricts to finite witnesses on every finite
grid, so the negative relation implies Q(n,t) SAT for all t; UNSAT at some t is evidence
for the positive direction."

### Why it was wrong
An infinite witness only guarantees that no FULL-TYPE (infinite) independent set exists.
A finite binary subgrid is not of full type, so nothing forces any finite subgrid to
contain a blue edge: the restriction of an infinite witness to a finite grid can be
empty. The implication as stated is unsupported in BOTH directions at the finite level.

### The correct framing (THM-453 part D)
The true bridge runs the other way and is a compactness statement:
- Q(n,t) SAT for ALL t ⟹ (König) a triangle-free graph on the full grid with no
  independent binary subgrid ⟹ ω^n ↛ (ω^n,3) with a STRONG witness.
- Hence positive relations FORCE finite cutoffs: ω^n → (ω^n,3) ⟹ R(n,2) < ∞.
  Specker at n=2 forces R(2,2) < ∞ even though SAT witnesses persist through t=10:
  the cutoff is real but large.
- A negative ordinal relation does NOT formally imply Q(n,t) SAT for all t
  ("strong witness" is a priori stronger than "witness").

### Lesson
For infinite-to-finite shadows, check WHICH quantifier compactness actually transports.
"Kill all full-type sets" has no finite trace on a single finite configuration; only
the universal finite family (all t at once) carries ordinal content, and it does so in
the direction finite-SAT-everywhere ⟹ infinite negative. Cf. MISTAKE-064 (parse the
statement before proving the wrong one).

---

## MISTAKE-067: incomplete subgrid-verifier falsely certified Q(n,t) SAT — R(2,2) is actually 5, not >14

**Date:** 2026-06-09 (same session, caught by structure-reading the "witnesses")
**Found by:** mac-mini-2026-06-09-S1
**Affects:** `erdos592_treegrid_pysat_macmini_s1.py` (find_independent_binary_subgrid),
`erdos592_invariant_quotient_macmini_s1.py` (imports it), results files
`erdos592_treegrid_pysat_*.out`, `erdos592_treegrid_push_*.out`,
`erdos592_invariant_quotient_*.out` (their SAT lines beyond the corrected thresholds);
fixed + fully re-run in `erdos592_verify_fix_macmini_s1.py`

### What happened
The CEGAR loop's final certificate ("no independent binary subgrid exists in the model")
used a backtracking search that committed to the FIRST consistent subtree under each
chosen child and never explored alternative subtrees of the same child. The search was
therefore incomplete: it could return "none found" when an independent subgrid existed,
falsely certifying SAT. Q(2,t) was reported SAT through t=14 ("R(2,2)>14"); with a
complete generator-based search, **Q(2,5) is UNSAT: R(2,2) = 5 exactly**
(Q(2,4) SAT with 35 edges).

### How it was caught
Reading the structure of the t=11 "invariant witness": its printed B_1 visibly missed
rectangles like {2,3}×{4,5}, which a genuine witness must hit — the certificate and the
object contradicted each other on inspection.

### What was NOT affected
- All UNSAT results (solver-side, sound): the n=1 calibration, and the corrected runs.
- Labs 1–2 avoidability conclusions: those used POSITIVE results of exists_S_free_grid,
  which constructs explicit grids whose pairwise pattern-checks happen at append time —
  sound. (exists_S_free_grid has the same incompleteness on its FALSE answers; its
  False-derived "caps" are lower bounds only — none were used for headline claims.)

### Lessons
1. In a CEGAR loop the FINAL "no counterexample" check is a proof obligation — it must
   be a complete decision procedure, not the same heuristic used to generate clauses.
2. Read the witness, not just the verdict: printing the (R, B_g) structure exposed the
   hole immediately (cf. MISTAKE-015's "the error was visible in the output").
3. Greedy-commit backtracking (commit to first consistent subtree, never revisit) is a
   recurring trap in tree searches — both grid searchers this session had it.

---

## MISTAKE-068: Cycle-Anchored Subset DP Reused for Longest-Path Problems

**Date discovered:** 2026-06-09 (caught in-session by branch-II self-check P3; never propagated)
**Found by:** kind-pasteur-2026-06-09-S2 branch II (blowup spectrum lab)
**Affects:** any script computing longest paths with a min-vertex-anchored subset DP

### What was assumed
A subset DP that anchors each cycle at its minimum-label vertex (correct for cycle enumeration:
every cycle can be rooted at its minimum) was reused for LONGEST PATHS by anchoring paths at
their minimum vertex.

### Why it was wrong
A path's minimum-label vertex can be INTERIOR — such paths are never generated when extending
only from the anchor. Symptom: longest path of P3 reported as 2 vertices; downstream, the
circumference law c(G[K2]) = 2p showed 275 fake "beaters" (true count after fix: 56).

### The correct framing
Cycles: min-anchoring is safe (rotation moves the minimum to the start). Paths: must allow
two-sided extension from the anchor, or run DP from every start vertex. Cycle spectra computed
with the anchored DP were NEVER affected.

### Lesson
"Anchor at the minimum" is a CYCLE symmetry (rotation), not a PATH symmetry (paths only have
reversal). Check which group acts before porting a canonical-form trick between objects.

---

## MISTAKE-069: "First Power-of-2 Cycle in Enumeration Order" Reported as "Smallest"

**Date discovered:** 2026-06-09
**Found by:** kind-pasteur-2026-06-09-S2 branch III (double-checked: own DFS checker ≡ networkx)
**Affects:** S710's Erdős-64 verification table (THM-446 context): "McGee → C16" and the cage list

### What was assumed
S710's cage battery reported, per graph, "the" power-of-2 cycle found — implicitly the smallest.

### Why it was wrong
The search reported the FIRST 2^k cycle encountered in enumeration order, not the smallest k.
McGee (girth 7) was reported as "→ C16" but in fact contains 34 EIGHT-cycles. (Petersen has 15
eight-cycles, not the 10 sometimes quoted.) The Erdős–Gyárfás CONCLUSIONS are unaffected (a
2^k cycle exists either way), but any downstream use of "McGee is C8-free" would be wrong —
e.g. it would corrupt the girth ladder of THM-457 (the true min order of a girth-7-or-more
cubic C8-free graph is > 46 by SA floors, not 24).

### The correct framing
When verifying "∃ cycle of length in S", always record the FULL dyadic profile (which members
of S occur), or at minimum the smallest. Exact anchor census now in canon (THM-457(1)).

### Lesson
"Found one" ≠ "found the smallest". Enumeration order silently masquerades as minimality.

## MISTAKE-070: Circular Proof Inventoried as PROVED (|Aut| | H "by tiling-count integrality")

**Date discovered:** 2026-06-10
**Found by:** kind-pasteur-2026-06-10-S1 Thread B (adversarially verified: independent canon sweep + line-by-line re-read of every cited source)
**Affects:** `gap_inventory_s196.out` item 15 ("|Aut| divides H for all tournaments: PROVED (opus S182)"); downstream confidence in THM-048 Step 3 and CLAUDE.md's tiling-fibration line

### What was assumed
That the universal divisibility |Aut(T)| | H(T) had been PROVED (opus-S182), so later artifacts cited it freely: THM-048 Step 3 says "by orbit-counting" with no proof; CLAUDE.md line ~351 states "Tilings * |Aut| = H for every iso class (orbit-stabilizer on tiling fibration)".

### Why it was wrong
S182's argument (`aut_H_deep_s182.out` Part 5) derives the divisibility from the integrality of the tiling count H/|Aut| — which is the SAME statement (it even hedges: "orbit size divides |Aut|, and different orbits may have different sizes"). The S20bt tiling formula behind the CLAUDE.md line was verified only at n=4,5 and implicitly ASSUMED freeness ("|Aut| are related by automorphisms"); `forbidden_tiling_counts.out` line ~241 even flags freeness as an assumption. The fact WAS true — but no proof existed anywhere in canon, while the inventory said PROVED.

### The fix
**LEM-003** (kind-pasteur-2026-06-10-S1): for ANY finite digraph, an automorphism fixing a directed Hamiltonian path's arc set is the identity (a directed Ham path's arc set determines its vertex sequence — the unique in-degree-0 source anchors an induction), so Aut acts freely on directed Ham paths, all orbits have size exactly |Aut|, and |Aut| | H. One paragraph, zero tournament structure. Exhaustively verified n≤6 (all 2^10+2^15 labeled tournaments, explicit orbits for all 3432 masks with |Aut|>1) + independently re-verified with different machinery. Honest boundary: FAILS for Hamiltonian CYCLES (no distinguished start): C3 has 1 Ham cycle and |Aut|=3 ∤ 1; circulant RQ5 has BOTH its Ham cycles rotation-fixed.

### Lesson
"X is an integer because it counts something" requires proving the count is well-defined — here, that all orbit sizes equal |Aut| (freeness), which was the entire content. An asserted proof citation can be as hollow as an extrapolated pattern (cf. MISTAKE-028/036/055). When an inventory says PROVED, spot-check the cited argument before building on it.

## MISTAKE-071: "Verified exhaustive n=4,6" That Checked Only One Maximizer Class (HYP-2312)

**Date discovered:** 2026-06-11
**Found by:** mac-mini-2026-06-10-S2 (det census) + independent adversarial re-verification (fresh exhaustive enumeration at n=6)
**Affects:** HYP-2312 ("the H-maximizing tournament has |Pf| = 1"), THM-442 section (3), and any plan to restrict A038375 extremal searches to Pf = ±1 tournaments

### What was assumed
That at n = 4, 6 the H-maximizing tournament has minimal Pfaffian |Pf| = 1, recorded as "VERIFIED exhaustive n=4,6", and conjectured for all even n (HYP-2312).

### Why it was wrong
At n = 6 the maximum H = 45 is attained by TWO iso classes (240 labeled each, both score (2,2,2,3,3,3)): one has |Pf| = 1 but the other has |Pf| = 7 (det S = 49; same in the det(I+2A) convention, so not a convention artifact). The earlier verification evidently examined only one maximizer class. At n = 8 the six H = 661 classes have |Pf| ∈ {1, 9, 17}.

### The fix
HYP-2312 amended to the EXISTENTIAL form: "at least one H-maximizing tournament has |Pf| = 1" — true at n = 4, 6 (exhaustive) and supported at n = 8. The existential form still justifies searching Pf = ±1 tournaments for the max VALUE of H, but such a search misses maximizer CLASSES. Universal form: REFUTED at n = 6.

### Lesson
"The maximizer" is a set, not an element. Any claim of the form "the extremal object has property P" must be checked on ALL extremal classes (argmax can be plural), and the verification record should state how many maximizers were examined.

---

## MISTAKE-072: "dim_nonspec(H) = n−5 / the overlaps are spectral shadows" extrapolated from n≤8 to all n

**Date discovered:** 2026-06-15
**Found by:** monad-explorer-2026-06-15 (n=9 carrier-dimension chain)
**Affects:** THM-505 (original dimension section), HYP-2513, reflection `the-zeta-function-and-the-ocf-read-complementary-halves` (closing paragraph)

### What was assumed
A carrier-dimension probe at n≤8 found that the non-spectral content of `H` is exactly the simple-cycle vector `(c6,…,c_n)`, of dimension `n−5`, with every overlap defect (`p33, TQ, Q44, TF`) a spectral function of it. This was stated for general `n`: "the simple cycles are the genuine hidden coordinates; the overlaps are their spectral shadows."

### Why it was wrong
At **n=9** the simple-cycle counts do NOT determine `H`. Nested-carrier chain over 130 000 cospectral samples: `sig→+(c6,c7,c8)→+c9→+Q44→+T333` splits `14804→482→24→1→0`. There are explicit cospectral witnesses with identical `(c6,c7,c8,c9)` but different `H` (each satisfying `ΔH = 4ΔQ44 + 8ΔT333`), and even `(c6,c7,c8,c9,Q44)` leaves `H` split. So `dim_nonspec(H) = 6 > n−5 = 4` (carriers `{c6,c7,c8,c9,Q44,T333}`), and BOTH `Q44` and `T333` are INDEPENDENT non-spectral carriers — not spectral shadows. The break occurs exactly when the triple level `α_3` (and the `(3,5)`-pair structure) gains room; n=8 was the last size where the higher-correlation configs were pinned by the simple counts for lack of room.

### The fix
THM-505 / HYP-2513 / the reflection all amended: dim = `n−5` ONLY for `n ≤ 8`; at `n ≥ 9` the non-spectral content is a TOWER of cycle-correlation rungs (simple counts; pair-overlaps `Q44`; triple-packings `T333`; …), each rung independent of those below once `n` gives it room. `H` is universal-linear in the full carrier set but not a bounded-degree polynomial in the simple cycles alone past n=7.

### Lesson
The repo's own refrain — "patterns that hold at n=3,4,5 often break at n=6 or n=7" — applies one notch higher here: a dimension/structure pattern verified exhaustively/by-heavy-sampling at `n≤8` STILL broke at `n=9`, precisely at the `n` where a NEW combinatorial level (here the triple disjoint-cycle packing `α_3`, needing `3·3=9` vertices) first has room. When a carrier/dimension count is conjectured "for all n," always test it at the first `n` where the next independent-set / correlation level can appear, not just one size up.

### Addendum 2 (monad-explorer-2026-06-15-S4): even the PACKING basis over-counts *for H* — H's own non-spectral dim is `⌊n/3⌋`, not `A000009(n)−3`
The S3 partition law `dim = #{λ:odd≥3,Σλ≤n}−3 = A000009(n)−3` (now identified as the named
sequence: the cumulative restricted-partition count equals partitions of n into distinct/odd
parts, by the one-line GF identity `Σ_{s≤n}[x^s]Π_{odd≥3}1/(1−x^k)=[x^n]Π_{odd≥1}=q(n)`)
correctly counts the rank of the **individual packing carriers** `(N_λ)`. But it does **not**
measure the non-spectral dimension of `H`. `H = I(Ω,2) = 1 + Σ_{j≤⌊n/3⌋} 2^j α_j` depends only
on the **level-sums** `α_j = Σ_{|λ|=j}N_λ` — it never sees the split of a level into its
length-types. Since `α_j=0` for `j>⌊n/3⌋`, **`dim_func(H)(n) ≤ ⌊n/3⌋` (LINEAR), PROVED**, and
`< A000009(n)−3` for n≥8. VERIFIED n=8: carrier basis `{c7,D33,D35}` rank 3 (the S3 number) but
level-sum basis `{c7, D33+D35}` rank 2 with `H` in span — `D33,D35` are independent functions
yet `H` reads only their sum. **Third lesson:** the chain of corrections (trace basis 6 →
packing basis 5 → level grading ⌊n/3⌋) shows that "the non-spectral dimension of `H`" must be
measured against *what `H` is a function of*, not against any convenient carrier list — and `H`,
being the scalar `I(Ω,2)`, factors through the coarse level-marginals `α_j`. The
`A000009(n)−3` law is true and beautiful, but it is the dimension of the finer **packing
vector**, not of `H`. See THM-505 TWO-DIMENSIONS section + reflection
`H-reads-only-the-level-grading`. (`04-computation/ocf_two_dimensions_monad.py`)

### Addendum (monad-explorer-2026-06-15-S3): the "6" was itself a basis over-count → the real law is a partition function
The corrected `dim_nonspec(H)=6` at n=9 above is itself **basis-dependent and over-counts by one**. It used the TRACE basis `{c6,c7,c8,c9,Q44,T333}`, but `c8` and `Q44` enter `H` **only through their sum** `c8+Q44 = D35` (the disjoint-(3,5)-pair count): the closed form's `+4c8+4Q44` is `4·D35` after `D35 = c3c5 − W8 + c8 + Q44`. In the basis-independent OCF *packing* basis `N_λ` (rank of the within-class carrier-delta matrix), the intrinsic dimension is **5, not 6** — carriers `{c7, c9, D33, D35, T333}`. The general law is `dim_nonspec(H)(n) = #{partitions of s≤n into odd parts ≥3} − 3` = 1,2,3,5,7,9,… (verified rank 3,5,7 at n=8,9,10). **Second lesson:** a "dimension" measured by a *nested chain over a fixed (trace) carrier list* counts basis vectors, not degrees of freedom — over-complete bases inflate it. Measure intrinsic dimension by the RANK of the carrier-delta space in the natural (here OCF packing) basis. See THM-505 growth-law section + reflection `the-non-spectral-dimension-of-H-is-a-partition-function`.

## MISTAKE-073 (2026-06-15, THM-503 infimum) — searched a SLICE, reported it as the infimum
THM-503 reported `inf L ≈ 0.0237` over primitive multiple-of-14 LRC configs, "attained at the near-tight cores {1,…,12}∪{14m}". This searched only the **end-drop** family (drop the LAST runner of the tight AP {1..13}). A broad search over **interior-drop** configs ({1..13}\{j}∪{14m}, j interior) goes ~4× lower: `{1..11,13,84}` has L≈0.00535, `{1..13}\{6}∪56` has L≈0.00561 (both verified loose, M=7/89 & 2/23, L stable to Q₀=24000). The true inf L ≈ 0.0053, not 0.0237; the minimizing drop is the MIDDLE (j=6), not the end. The qualitative inf>0 survived, but the value/extremizer/margin were wrong. Lesson: when reporting an extremum over a symmetric-ish search space, **sweep the full orbit of perturbations, not the first slice** — here "drop a runner of the tight AP" has 13 inequivalent positions and the end (most natural to try) is the WORST extremizer, the middle the best. Echoes the census-horizon lesson (MISTAKE-018) and the parameterize-the-exception lesson: the extremal object is rarely the most obvious member of the family. Fix: THM-503 CORRECTION addendum + HYP-2520.

## MISTAKE-074 (2026-06-15, THM-518 / OPEN-Q-104) — a Fourier-TRUNCATED Riesz ratio reported a false looseness certificate
Probing the Bedert Riesz certificate `∫M·R/∫R < 1` for the LRC(14) extremizer `{1..13}\{6}∪56`, I computed `∫M·R` via the Fourier side `Σ_v Σ_{|k|≤14} s(k) R̂(vk)` (truncating the sinc sum at `Kmax=14`) and got ratio **0.9506 < 1** — flagged as "CERTIFICATE!". The EXACT direct-grid integral `∫M·R = (1/Q)Σ_a M(a/Q)R(a/Q)` (Q=30030, no truncation) gives ratio **1.064 > 1** — NO certificate. The tail `Σ_{|k|>14} s(k)R̂(vk)` is NOT negligible: `s(k)~1/k` decays slowly and `R̂(vk)` is nonzero on the dense relation sum-set, so truncating at `Kmax=14` dropped a positive contribution that flips the verdict. Lesson: for a `∫(periodic)(Riesz)` pairing where the test function `M=Σ1_danger` has slowly-decaying (`1/k`) Fourier coefficients, **compute the integral DIRECTLY on a fine grid — never truncate the Fourier/sinc sum** unless the tail is rigorously bounded. The direct grid is exact and cheap here; the Fourier shortcut silently fabricated a certificate. (Self-caught the same session by running the direct-grid verification BEFORE writing any canon claim — the discipline "verify a too-good certificate by an independent exact method" is what saved it.) Fix: HYP-2547 caveat + THM-518 §C(1); all OPEN-Q-104 ratios now direct-grid.

## MISTAKE-076 (2026-06-17/18, THM-526 / HYP-2580a) — "criterion C verified universal" was a sampling artifact; C is NOT necessary
The LRC(14) proof-program reduced "M(S)≥1/14 for covering 13-sets" to the **criterion** `C(S): ∃v∈S, W(S\{v})>1/(7v)` (arc-width lemma: C ⟹ M≥1/14, genuinely PROVED). HYP-2580a/Angle-F then claimed `C(S)` holds **universally** ("verified on ~12,000 covering sets / 4000/4000, 0 failures"), making "prove C for all covering S" the target. **C is NOT universal.** Exact counterexample (found by two independent verifiers, reconfirmed directly): **S\* = {1,2,3,5,7,8,9,10,11,12,13,38,42}** — primitive, covering, case S3 (k=2, Vmax=42) — has **all 13 ratios W(S\*\{v})·7v < 1** (max 429/532≈0.806), yet **M(S\*)=2/23 ≥ 1/14** (a global witness, not any single-removal arc). So "prove C(S) for every covering 13-set" is a **false target** — it can never close LRC(14). Lesson: a SUFFICIENT criterion validated only by sampling can be silently non-necessary; the C-failure locus here is sparse (≈1 in thousands) and slipped through every Monte-Carlo sweep. Before elevating "criterion holds on all X" to a proof target, **search adversarially for the criterion's failure set** (it is cheaper to refute universality than to prove it), and keep the direct quantity (here M itself, via the global witness / finite check) as the fallback. The PROVED pieces (S1, S2, the k=2 slice, cluster-collapse Lemma A) never used C-universality and stand. Fix: THM-526 correction section + HYP-2580a marked REFUTED + HYP-2581; a correct closing lemma must remove ≥2 runners or use a global witness. Echoes MISTAKE-073 (searched a slice, called it the infimum) and the census-horizon lesson (MISTAKE-018).

## MISTAKE-077 (2026-06-18, μ_θ exact engine) — naive order-change breakpoints UNDER-count the gap measure
Computing μ_θ(E)=meas{x: maxgap{frac(e_i x)}>θ} exactly requires breakpoints where the cyclic ORDER of the orbit points changes (collisions (e_i−e_j)x∈ℤ) AND, WITHIN each order-cell, the analytic crossings where an affine gap equals θ. An engine that used ONLY order-change breakpoints (sampling each cell's midpoint) gave a WRONG, smaller value for μ_{1/7}(consec_9) — 4829/5880 instead of the correct 247/294 — because a gap can cross θ strictly inside an order-cell. The fix (engine A / the workflow's mu_theta): on each order-cell every cyclic gap is affine in x, so {maxgap>θ}=∪_t{affine_gap_t>θ}, each a sub-interval solved exactly (s·x+c>θ); take the union length. Cross-validated against a brute grid with denominator dividing 7·lcm(differences) (which cannot be wrong). Lesson: for a piecewise-linear extremal (max/min of affine pieces) the level-set {f>θ} has breakpoints BOTH at the piece-swaps AND at the θ-crossings inside pieces — sampling cell-midpoints silently drops the latter. Always solve the affine=θ crossings, or cross-check against a provably-fine rational grid. Affects any μ_θ/EWLB computation. → THM-530, HYP-2602.

## MISTAKE-078 (2026-06-19, LRC(14) wide-spread bound) — "verified on canonical families + loose ceiling" mistaken for "uniform theorem"; the envelope tail DIVERGES
Closing the LRC(14) residual to the wide-spread cap bound (meas(S7(E))≤cap_k for span>B), I argued (HYP-2611b, kps-S9) that the route was "structurally complete, remaining is engineering": every tested wide regime sits ≤~0.21 ≪ cap (margin 0.17), the tight case is the exact finite check, so the wide bound is "loose." This conflated **two different things**: (i) the inequality is TRUE and loose on the tested families (~40k exact shapes, 0 violations) — correct; (ii) a UNIFORM-over-all-wide-E PROOF is easy — WRONG. The S10 assembly + 3 verifiers showed the support-6 correction tail eps(B), bounded by the free per-coordinate envelope Σ_{m,7∤m} c1/|m| (c1=0.697303), **DIVERGES HARMONICALLY** (partial sum to 10⁵ = 7.42). So no finite tail bound follows from the envelope; a successive-minima/Minkowski lattice-point count over the support-6 relation lattice (|K(n)|≤c1⁶/(λ₁···λ₆)) is REQUIRED and was never executed. Lesson: a LOOSE margin + exhaustive sampling on structured families is strong EVIDENCE but is NOT a uniform inequality — the uniform statement can still need a genuine (here lattice-geometry) argument that the per-term envelope, summed freely, cannot supply (harmonic divergence). When claiming "the remaining work is engineering," check that the summation/quantifier over the INFINITE part actually converges with an EXECUTED bound, not just that each instance is safe. Echoes MISTAKE-076 (sampling ≠ universality) and MISTAKE-073 (a slice ≠ the infimum). Fix: HYP-2611b corrected; the residual is HYP-2608(a) the wide-spread bound = the Minkowski tail count, OPEN. Everything else (k≤7, glue G1, per-E bound, caps, bounded finite check to span 16) is gap-free.

## MISTAKE-078 AMENDMENT (2026-06-19, CASE-thm538 resolved) — the "support-6 / ≥6-body" mechanism was ALSO wrong, not just the uniformity
MISTAKE-078 flagged that "verified-on-families ≠ uniform theorem" (correct, stands). But its framing also leaned on THM-538's "support-6 floor ⟹ the wide correction is a ≥6-body object" — and that mechanism is itself FALSE (CASE-thm538-support6-floor-zero-padding, conceded): the support-6 floor holds only for the active-coordinate sum Q, NOT the zero-padded kernel K that appears in the measure (short support 2–5 relations DO contribute; the AP's correction is support-3-dominated, not support-6). The CORRECT structure (HYP-2646): K(n)=D7(n mod 7)/∏n_j, correction = Σ_c D7(c)S_c(E) conditionally convergent, ruled by support-6 relation DENSITY R6 (HYP-2645), with the convergent representation being the finite x-cell integral / far-element plateau (HYP-2644/2610), NOT the box-truncated lattice sum. Double lesson: (1) a "verified exhaustively, max 5e-17" check can silently compute a DIFFERENT object (Q, active coords) than the one in the statement (K, zero-padded) — always confirm the verified quantity IS the claimed quantity on a case where they should differ; (2) a clean algebraic vanishing ((1−1)^{6−|U|}=0) can be an artifact of dropping a |T|-dependent factor — re-derive keeping every factor. The wide-spread bound's real content is R6 density + far-element decorrelation, unaffected.


## MISTAKE: L7 'CLOSED' overstated (mac-mini-2026-06-21-S13, caught by the S13 rigor-audit workflow)
**The error:** S13's broadcast / HYP-2738 / reflection claimed the LRC(14)-S3 sector route was "CLOSED" and the audit "ALL PASS". The adversarial audit (Thread A, completed after close-out) found this OVERSTATED. Actual: L7 is REDUCED, not closed.
**Invalid step (L7 Step 4):** the closure (lrc_q108_L7_closure_kps.md) cites the finite-f1 convergence |p0(B u far)-p0_inf| = O(1/f1) as "= THM-546 (PROVED)". FALSE: THM-546 peels ONE far element with the REMAINDER BOUNDED. In the L7 limit BOTH far elements grow (f2 = gamma*f1), so a single THM-546 peel of f2 leaves E'=B u {f1} UNBOUNDED (V = Theta(f1)) => bound (6/7)/gamma ~ 0.43 = O(1), NOT O(1/f1). The proved D_{p,q} <= 14/p bounds the discrepancy of the LIMIT LAW mu_{p,q}, a DIFFERENT object from the finite-f1 convergence RATE. The closure conflated them.
**True status:** the O(1/f1) rate IS true (verified broadly: |err|*f1 bounded ~<= 0.75) but NOT proven; it needs a JOINT 2D Erdos-Turan/Koksma bound peeling the BOUNDED BASE from the FAST FAR-PAIR (not a single peel of f2). Verified-not-proved also: the r>=3 -> pairwise reduction, base-size domination, and "consec maximizes meas(S7)" (HYP-2602, open as a theorem).
**Lean caveat (gap #7):** delsarte_bound_k8/k9/k11 prove the per-shape q0 <= L_y; the content is the readout identity (moment-LP functional = dual covector), and the bound is then trivial (q3,q6>=0). It does NOT formalize L_y <= cap (the extremality / cap content).
**Lesson:** "reduced to a finite atlas + a numerically-verified rate" is NOT "closed". Do not cite THM-546 (single-far, bounded remainder) for the two-far joint limit. The genuine remaining input is the joint 2D ET-Koksma rate. -> HYP-2730, HYP-2738, THM-546, lrc_q108_L7_closure_kps.md.

## MISTAKE-86: computing M(S)=max_t min_s||st|| with an INCOMPLETE breakpoint set underestimates M and reports FALSE tights

**What happened (kind-pasteur-2026-06-28-S256):** A census probe over single-speed replacements of the AP used candidate breakpoints t = k/(2 s_i) and k/(s_i - s_j) only, OMITTING t = k/(s_i + s_j). This made M_of UNDERESTIMATE M(S) (it missed the true maximizer, which often sits at an s_i+s_j crossing), so ~23 single-replacements were reported as "tight" (M=1/14) when they are actually LOOSE (M>1/14). With the complete breakpoint set the census collapsed correctly to the single extra tight set GW (12->24).

**The fix:** f_S(t)=min_i ||s_i t|| is piecewise linear; its local maxima (the candidates for max_t) occur where some ||s_i t|| has a kink (t=k/(2 s_i)) OR where two terms cross. Two terms ||s_i t||=||s_j t|| cross at BOTH t=k/(s_i - s_j) AND t=k/(s_i + s_j). The full candidate set is {k/(2 s_i)} U {k/(s_i - s_j)} U {k/(s_i + s_j)}. Sanity anchors: M(AP)=1/14, M(GW)=1/14, M({1..11,13,26})=1/12. Always include the s_i+s_j family; cross-check with a fine rational grid.

**Lesson:** A "tight" verdict from an exact-rational M computation is only as good as the breakpoint set. Incomplete breakpoints give a one-sided error (underestimate M) that fabricates tights -- the most dangerous direction for a census.

## MISTAKE-087 (2026-06-30, mac-mini-S47) -- the construction n/Phi_6 was assumed to be the COVERING-MIN for n>=7 from a heuristic + a restricted scan; it is NOT (beaten exact at n=7,8,9)

**What happened:** HYP-3701 (my own, S42) inferred a "transition at n=7": drop-2/mediant 2/(2n-1) is the covering-min for n<=6, and the construction {1,..,n-2,n(n-1)} = n/Phi_6(n) takes over as the covering-min for n>=7. The inference rested on (a) the PG(2,6)-failure heuristic (the first projective plane fails at q=6=n-1) and (b) opus's "107-set scan confirmed 14/183" -- a scan of NEAR-CONSTRUCTION variants. The whole subsequent arc (HYP-3703 tiling-optimality, 3704 three-routes, 3717 three-gap, 3722 observer-escape; the Kershner/Eisenstein/A2 framing) built on "convergent = covering-min for n>=7."

**The refutation (exact, dense-grid cross-checked):** smaller-M primitive covering (n-1)-sets exist at every n=7,8,9:
- n=7: {1,2,5,6,7,8}, M=2/13=0.15385 < 7/43=0.16279 (t=4/13)
- n=8: {1,4,5,6,7,11,16}, M=2/15=0.13333 < 8/57=0.14035 (t=8/15)
- n=9: {1,3,4,5,7,11,18,32}, M=4/33=0.12121 < 9/73=0.12329 (t=29/33)
So there is NO transition at n=7; the sub-convergent (mediant at n=7,8) keeps beating the construction. opus's 14/183 is a restricted-family min, NOT the global covering-min.

**Why it slipped through:** the beaters are SPREAD-STRUCTURED (a speed of 32 ~ 3.5n at n=9), so near-construction perturbation scans and low-speed exhaustion both MISS them. The PG-failure heuristic was a post-hoc fit to n<=6, not a mechanism.

**Not fatal to LRC:** all candidates (mediant 2/(2n-1), 4/33, convergent) are > 1/n; the covering-min being smaller just makes the floor margin tighter (~1/(n(2n-1)) vs ~1/n^2). The LRC floor M>=1/n is untouched.

**Lesson:** (1) a "covering-min" claim must be tested against SPREAD-structured covering sets (large speeds), not only near-construction perturbations or low-speed exhaustion. (2) An elegant structure on a particular covering set (the construction's hexagonal AP / three-distance / Eisenstein) does NOT make that set extremal -- do not conflate "a beautiful covering" with "the minimal covering." (3) A transition inferred from a number-theoretic coincidence (PG(2,6)) needs exhaustive confirmation at the first n past the alleged transition. -> HYP-3725, HYP-3701, CASE-convergent-not-covering-min.

## MISTAKE-088 (2026-06-30, mac-mini-S52) -- claimed a clean "n>=12 => covering-min = 1/(n-1)" from an ILP whose speed bound V was BELOW the construction scale n(n-1)

**What happened:** The covering-min ILP (HYP-3731/3732) with speed bound V=72 returned M_prim=1/(n-1) for n=12,13,14,15 (clean!), and I wrote it up as a transition at n=12, the LRC14 hard core = 1/13, and a pinning of HYP-2566's looseness as 1/(n(n-1)). 

**The error:** V=72 is BELOW the construction's largest speed n(n-1) (=156 at n=13, 182 at n=14). The ILP literally could not see the construction {1..n-2, n(n-1)}, whose M = n/Phi_6(n). And n/Phi_6 < 1/(n-1) for ALL n (n^2-n < n^2-n+1), with the construction a valid primitive covering set. So M_prim(n) <= n/Phi_6 < 1/(n-1) ALWAYS -- 1/(n-1) is NEVER the covering-min. The ILP's 1/(n-1) was just the best LOW-SPEED primitive set; the clean "n>=12" regime was a search-bound artifact. (klein-S36 had flagged exactly this under-resourcing in the same-numbered HYP-3731.)

**The fix / what survives:** EXACT (V=4n suffices, klein-confirmed): n=7..11 covering-min 2/13,2/15,4/33,4/37,3/31, depth a(n)=2,2,4,4,3, the Stern-Brocot ray [0;n-1,k] frame (floor k=1, covering-min k=a(n), construction k=n). RETRACTED: the n>=12 clean regime, the n=12 transition, LRC14=1/13, the HYP-2566 pinning. For n>=12 the covering-min is <= n/Phi_6, exact value OPEN (needs V ~ n(n-1)).

**Lesson:** before declaring a computed extremum (min/max), ALWAYS evaluate the KNOWN constructions as upper/lower bounds and check the optimizer's search range covers their scale. Here a single check -- "is the known construction n/Phi_6 below my ILP's answer?" -- would have caught it instantly (it is, by 1/(n^2-n+1)/... a hair). A clean formula emerging right at the edge of the search range is a red flag for under-resourcing, not a discovery. -> HYP-3732, HYP-3731, HYP-2566, THM-523.

## MISTAKE-089 (2026-06-30, klein-S52, caught mid-session before canon) -- searching "covering-min escapes" over ALL sets instead of COVERING sets returns the tight/mediant minimizers as PHANTOM escapes

**What happened:** working the lowness lemma's large-speed residual, I searched sets S = {1..n-2}\{k} U {2 speeds} for the min M, imposing only "kills resonances k, n-1." The search returned M = 1/n (GW) at n=14 and M = 2/(2n-1) at n=10,12,16 -- all BELOW n/Phi_6 -- which looked like the lowness lemma FAILING (a low-M set missing a core speed).

**The resolution:** those minimizers are NOT covering sets. A covering set (THM-523) must contain a multiple of EVERY q in {2,...,n}; GW = {1..11,13,24} and the drop-2 mediant sets miss a multiple of n, so they are THM-523's TRIVIAL class (killed by the q=n witness with M=1/n exactly), NOT things the covering-min ranges over. Re-running with the full covering constraint (mult of every q in {2..n}) gives the genuine single-drop covering-escape minima 5/43, 2/21, 7/89 (n=10,12,14), all > n/Phi_6 as expected (and 4/29 < 8/57 at n=8, the real small-n failure, MISTAKE-087).

**Lesson:** the covering-min and the LRC floor 1/n live on DIFFERENT set-classes. Any covering-min / lowness search MUST impose the covering predicate (mult of every q in {2..n}) up front; otherwise the tight `1/n` and mediant `2/(2n-1)` minimizers (which sit BELOW n/Phi_6) masquerade as escapes and appear to refute the lowness lemma. The large-multiple-forced lemma (HYP-3763) itself is unaffected -- it holds for any low-M set, and GW does carry the forced multiple 24=12*2. -> HYP-3763, THM-523, HYP-3701.

## MISTAKE-091 (2026-07-01, mac-mini-S93, caught in-session before canon; renumbered from 090 -- opus-S32 claimed 090 concurrently) -- applying the tight-set slope formula to NON-PRIMITIVE sets under-counts witnesses by the dilation factor; formula-only censuses without direct cross-verification

**What happened:** the THM-593 tight-set slope formula `c_S = (2/q) sum_{u unit} 1/v_max(u)` silently assumes the argmax of `f_S(t)=min_v||vt||` is exactly the `phi(q)` unit fractions `a/q`. Applied in a cross-modulus census to `{2,4,6,8}` (= 2*AP, non-primitive) at q=5 it returned `5/12` -- but the lonely measure is dilation-invariant (`m_{cS} = m_S` via `u = ct`), so the slope must be `5/6`. A dilated set has `c*phi(q)` witnesses (at `a/(cq)`); the formula misses `(c-1)*phi(q)` of them. The first landscape output also reported this phantom `5/12` as a "beater."

**The resolution:** (1) always normalize to gcd 1 before applying witness-count-sensitive formulas (dilation invariance of measure/slope is the free cross-check: `slope(cS) == slope(S)` catches the bug instantly); (2) the formula is in general a LOWER bound `c_S >= (2/q) sum 1/v_max(u)`, with equality iff the tight set is "clean" (argmax = unit fractions) -- all primitive tight sets found at q=5..16 are clean, verified by DIRECT exact measurement in the last linearity cell; (3) never publish a formula-only census: the corrected landscape recomputes every slope directly and prints `formula=direct` per row.

**Lesson:** a closed form derived under a witness-structure hypothesis must carry that hypothesis explicitly, and any census built on it needs an independent direct computation per entry (here: exact measure in the last Farey cell). Dilation/scale invariances are free consistency tests -- run them on every new invariant. Secondary workflow lesson: `timeout N python | tee out` loses ALL buffered output on timeout with exit 0 (tee masks the kill); use `python -u` + per-section flush for long sweeps. -> THM-593, HYP-3840.

## MISTAKE-100 (2026-07-03, kind-pasteur-S37, caught by mac-mini-S31 THM-612 next session) — claimed "the tight locus is a SINGLE family (the AP), NO GW" from an under-powered search

**What happened:** searching for tight families (M=1/14) over APs {a,a+d,..}, dilates c·{1..13}, and ~2500 RANDOM 13-subsets of {1..30}, I found only dilated APs and concluded (HYP-4062) "primitive tight = {1..13} UNIQUE, no GW family." mac-mini S31 (THM-612) then verified **GW = {1..11,13,24} = AP[12→24]** is tight (M=1/14 EXACT, primitive, NOT a dilated AP, 6 tight points at units/14, non-covering) — a genuine second tight family.

**Why it was wrong:** the search was STRUCTURALLY BLIND to the GW shape. Random 13-subsets of C(30,13)≈1.2·10⁸ at 2500 samples has ~0 chance of hitting the specific one-residue-moved family {1..11,13,24}; and my structured candidates (APs + dilates) exclude it by construction (it is an AP with a single entry moved, neither an AP nor a dilate). Same class as mac-mini's own weak-adversary traps (MISTAKE-097/098): concluding "does not exist" from a search whose generator cannot produce the object.

**The correct framing:** the tight locus is at least {AP, GW}, both non-covering and small-speed — so the S37 REDUCTION and the 14-grid REPULSION stand (they only need tight ⇒ non-covering, which both satisfy); only the uniqueness/no-GW classification was wrong. LESSON: before claiming "no X exists" from a search, ask what shapes the generator CANNOT produce, and construct the adversarial shape explicitly (here: perturb the AP by moving one residue). A random/structured search is evidence of abundance, never of absence.

## MISTAKE-099 (2026-07-03, mac-mini-S30, caught in-session ~15 min in) — spent the session's opening re-deriving the Φ6 covering-min construction that I MYSELF refuted in a court case I filed 3 sessions earlier (S47)

**What happened:** tasked to "work the covering-min core", I conjectured `covering-min(n) = n/Φ6(n)` with tight family `{1..n-2,(n-1)n}`, verified `M(T_n)=n/Φ6` exactly + `apex≡-1 mod Φ6` + `g=2`, and started building a general theory around it — before recalling that `CASE-convergent-not-covering-min` (which **I filed**, mac-mini-S47, opus CONFIRMED S32) already REFUTED exactly this: the construction is NOT the covering-min for `n≥7` (beaten by 2/13, 2/15, 4/33 at n=7,8,9). The true covering-min is `1/n` over ALL covering families (opus-06-30 even block `2·{1..13}`, imprimitive) and empirically `~7–11%` above `1/n` for primitive ones.

**The resolution:** reading `02-court/active/` + `definitions.md` surfaced the correction within ~15 min; the session then produced real results (THM-610 deep-hiding lemmas; n=11 counterexample 3/31; margin map). Not a wasted session, but the opening detour was avoidable.

**Lesson:** the CLAUDE.md Step-5b "scour for leads" scan must include **`02-court/active/` for the exact quantity you're about to work** — especially cases you filed yourself. A named target ("the covering-min", "14/183") should trigger a `grep -rl` of court + hypothesis INDEX for that number BEFORE deriving. My own S47 filing sat one `ls 02-court/active/` away. Related recurring pattern: the "weak/wrong adversary" family (MISTAKE-090/098) — here the dual error, re-deriving a known-refuted *upper* construction. -> CASE-convergent-not-covering-min, THM-610, MISTAKE-090.

## MISTAKE-112 (kind-pasteur-2026-07-05-S15/S17): the consecutive multi-fold law M(D_l) = 14/(14(13-l)+1) is FALSE at l = 4, 5

The S15 instance verified the law on a sub-range and claimed it for l = 1..6 (HYP-4177). Exact merge-grid enumeration (S17): M(D_4) = 17/155 (claimed 14/127) and M(D_5) = 19/155 (claimed 14/113) -- the actual values are BELOW the law (the binding pair migrates to the 154+1 pair grid mid-ladder; l = 1, 2, 3, 6 match). The 2/25 FLOOR survives at every rung and is what HYP-4212's domination assembly consumes (LRCMultiFoldRows.lean certifies all six).

LESSON (the recursion of MISTAKE-102): a closed-form law verified at the ladder's ends and one midpoint is NOT verified -- the binding-pair structure can migrate at interior rungs. Enumerate the FULL parameter range exactly before stating a law; the merge grid makes this cheap (THM-592/HYP-4108: denominators are pair sums).


## MISTAKE-113 (mac-mini-S7, recurring from sibling-S6): adversarial free-fraction minimization on a FIXED grid masks sub-grid free sets

**What:** S5 (HYP-4282) concluded "consec[1..11] combs TILE the circle at radius 2/25 (phi_worst -> 0)" from a coordinate-descent that MINIMIZED the free fraction measured on a grid of 1600 points.  S7 exact arithmetic shows the free set has measure **0.000529** -- NONEMPTY -- but smaller than one grid cell (1/1600 = 0.000625), so the grid read it as 0.  Distinct-freq combs do NOT tile; kps's CircleClearFloor is correct.

**Why it's instructive:** minimizing a grid-estimated measure REWARDS configurations whose true free set hides between grid points -- the optimizer drives toward sub-grid-resolution free sets and reports 0.  The sibling hit the SAME trap in S6 (v1 adversarial "gridmax 0" artifacts).  Two instances, same lesson.

**Fix:** for covering/tiling questions use EXACT rational arithmetic (test midpoints of the elementary intervals cut by all arc endpoints -- finitely many, exact) OR verify any claimed phi=0 at >=100x the search resolution before believing it.  Never trust a minimized grid-measure at its floor.


## MISTAKE-114 (kind-pasteur-2026-07-06-S36, guard-railed by opus-S118 HYP-4506 and self at S38): the "window is too narrow / Dx<D" width narrative is a per-family SYMPTOM, not the obstruction -- the root is the ARITHMETIC of 3N+2

**What:** kps-S36 (HYP-4517) framed the n=12 first-gap emptiness via a metric mechanism -- the resonance-ladder crossing width Dx is smaller than the resonance spacing D, so the grid "skips" the gap. Read as *the* reason, this is misleading. opus-S118 (HYP-4506) proved first-gap emptiness is NON-MONOTONIC in N: N=13 is NONEMPTY (mediant 3/41 attained by {1..11,13,36}) while the WIDER-window N=12 is EMPTY. A width/Selberg story is monotonic (narrower => harder), so width cannot be the deciding quantity.

**Why it's instructive:** the deciding quantity is arithmetic -- the mediant 3/(3N+2) is achievable iff 3N+2 is PRIME (N=7->23, N=13->41 prime, nonempty; N=12->38=2*19 composite, empty). The metric Dx<D is a true computed fact about the single-outlier ladder subfamily, but it is downstream of the arithmetic: whether a base's resonance grid ALIGNS to the mediant is decided by the factorization of 3N+2, not by a width budget. kps-S38 verified the reconciliation (opus's nonempty witnesses ARE ladder families at N=7,13) and refuted a too-clean spin-off hypothesis ("F(N)={1..N-2,N}+3(N-1) gives the mediant iff 3N+2 prime" -- FALSE: F gives the mediant only at N=7,13, not at primes N=5,9,15).

**Fix:** keep Dx<D as the constructive/symptom side (the ladder is the WITNESS builder) and treat the ARITHMETIC of 3N+2 (mac-mini HYP-4562/4572 mod-19 clearance; opus O-arith / Fan-Sun gcd template) as the obstruction side. Do not present width/Selberg as the root cause of (G); it is non-monotonic and thus cannot be.

-> HYP-4282 (S5, the artifact), HYP-4312 (S7, the resolution), HYP-4292 (sibling S6, same trap), kps CircleClearFloor.


---## MISTAKE-114 (2026-07-06, mac-mini-S20) -- the FAST-EXACT-M helper (S16-S19) skipped non-coprime witness numerators, UNDERESTIMATING M.

**What was wrong.** My fast exact-M `Mfast` (offered to the fleet S16 as an O(n^2*max) exact-M via the witness-denominator lemma q|(v_i+-v_j)) had `for a in range(1,q): if gcd(a,q)!=1: continue`. The lemma is CORRECT (M's reduced denominator DIVIDES a pairwise sum/diff), but a witness at a SUB-denominator q' dividing (v_i+-v_j) appears over q=(v_i+-v_j) as a NON-coprime numerator (q/q')*a', which the skip DISCARDED -- so Mfast could MISS the true (larger) M.

**Caught by.** n=6 set {1,3,4,5,18}: Mfast gave 4/23=0.174 (looked like a gap member in the second gap (1/6,2/11)); the independent fine grid gave 2/11=0.1818 (the BOUNDARY, loose). The q=11 witness appears as 2a/22 with 22=4+18, and was skipped -- a FALSE POSITIVE.

**Fix.** Remove the gcd skip; check ALL a in [1,q). Re-verified vs grid: AP=1/13, doubled-apex=2/25, block=2/25, single-lift {1..11,23}=1/12, n=7 {1,5,6,11,16,17}=5/33 all correct.

**Impact (assessed, mostly benign).** (a) This session's 'n=6 gap member' = FALSE (it is 2/11, loose). (b) S16 targeted near-AP search RE-RUN with fixed M: still 0 in gap (15,976 families) -- conclusion HOLDS. (c) S17 n-specificity: the n=7 gap member 5/33 is GRID-confirmed independently -- STANDS. (d) S18 equioscillation AP=phi(n): counts UNCHANGED with fixed M (AP witnesses at a/13 coprime, 13 prime) -- STANDS. (e) S19 Fekete: direct energy, unaffected. The bug underestimates M so it risked FALSE NEGATIVES in searches; the n=13 emptiness rests on the fleet's correct-M exhaustive work (concurrent lift census), not my buggy searches. Files fixed: lrc_fastM_highscale_probe / lrc_leaveoneout_alignment / lrc_witness_denominator_dichotomies / lrc_equioscillation_count _macmini_S1x.py.

## MISTAKE-115 (opus-2026-07-06-S122, self-correcting opus-S120/S121): 'gap member = (N-2)-AP + exactly 2 defects' is FALSE -- 3-defect gap members exist; the defect count does NOT govern (the ORDER does)

**Claim (S120, wrong):** every LRC first-gap member is an (N-2)-term dilated AP + exactly 2 defects (longest-AP = N-2); the crux residual is a '>=3-defect Freiman-stability exclusion'.

**Refutation (S122, exhaustive at N=7):** {1,3,4,5,7,13,18} is a gap member (M=3/23 in (1/8,2/15)) with longest-AP {1,3,5,7} = 4 = N-3, i.e. 3 DEFECTS.  It coexists with the 2-defect member {1,2,3,4,5,7,18} at the SAME value M=3/23 (order 2).  So (a) 3-defect gap members exist, refuting the 2-defect signature; (b) the defect count is NOT the governing parameter -- two families of the same order (2) have different defect counts (2 vs 3).  The S120 signature was over-fit to 3 examples (N=6,7,13 first members) that happened to be 2-defect.

**Correct frame:** the governing parameter is the ORDER k of the value s/(Ns+k) (opus S116/S117), and the crux is kps's ACHIEVABILITY GAUNTLET (HYP-4557): in-gap values exist at every order for every N, and (G) at N=12 is that EVERY order's value is unattained -- a uniform-over-orders exclusion, NOT a bounded-defect one.  {1,3,4,5,7,13,18} is kps's 'no-isolated-runner species'.

**Lesson:** a structural signature read off 3 examples is not a theorem; exhaustive enumeration at the smallest nonempty cases (N=6 gave 1 member, N=7 gave 2 of different defect count) is cheap and would have caught it.  The proof map (00-navigation/LRC14-PROOF-MAP.md) crux line is corrected.

## MISTAKE-127 (klein-2026-07-08-S193, self-caught): the S192 "large-spread half via arc-count pigeonhole #arcs < rho*Vmax, prove #arcs <= c*spread with c<1" is the WRONG TOOL -- VACUOUS on the extremal near-dilated-AP family (where c ~ 1.17 > 1 AND rho* ~ 0.6), even though good periods abundantly exist

**Claim (klein-S192, wrong route):** the large-spread half of THM-527-A closes via the pigeonhole #{good ruler periods} >= rho*Vmax - #arcs(G*), so it suffices to prove [#arcs <= c*spread with explicit c<1] + [rho* >= rho0 > c]. Verified "zero failures" over 25 RANDOM primitive clusters/spread (c~0.2, rho*~0.99), and flagged c<1 as the residual.

**Refutation (S193, exact):** on the EXTREMAL near-dilated-AP family E_d=d*{0..9}u{p} (the low-rho* shape, longest-AP=10), #arcs ~ 1.17*spread (block-like, (k+1)/(k-1) > 1, so c<1 is FALSE) and rho* ~ 0.60, so the pigeonhole bound rho*Vmax - #arcs ~ -1545 is VACUOUS (deeply negative). Yet a good period abundantly exists: #good ~ 1612 ~ rho*Vmax (d=300), and #good/Vmax -> rho* = 0.594 with |#good - rho*Vmax| <= 7 -- the TRUE discrepancy is O(1), NOT #arcs (=3170). The random-cluster test missed this because random e is generic (small c, large rho*); the near-AP extremal was never tested.

**Root cause:** the crude bound |#good - rho*Vmax| <= #arcs is Koksma-Hlawka with the grid discrepancy 1/(2Vmax) x the total variation 2*#arcs of 1_{G*} treated as an ARBITRARY union of #arcs intervals. That is blind to the fact that the arcs of G* and the ruler grid {(j+1/2)/Vmax} share the SAME Vmax-arithmetic (both come from the phases frac(e_i x)), so they are not adversarially aligned and the real discrepancy is tiny. Only the Erdos-Turan analysis of the STRUCTURED indicator 1_{G*}=F(frac(v.x)) -- discrepancy driven by resonances a.v = 0 (mod Vmax) -- sees the cancellation.

**What survives / correct frame:** the CONCLUSION (a good period exists for large spread) stands, verified. The correct route is Erdos-Turan: |#good - Vmax*rho*| <= Vmax*D*, D* <= C_m(1/(H+1) + sum_{Vmax|a.v, ||a||<=H} 1/r(a)). PROVED sub-result: for near-AP the low-height resonances are d-INDEPENDENT (AP-supported a give a.v = d*sum(i*a_i), and Vmax=9d+14 forces sum(i*a_i)=0 for bounded a -- identical resonance set d=5..100), so the discrepancy stays spread-uniform => #good/Vmax -> rho*. THM-527 part H + THM-663 corrected.

**Lesson:** (1) an arc-count / total-variation (Koksma-Hlawka) bound is the WRONG instrument for "does an arithmetic grid hit a set built from the SAME arithmetic" -- it ignores the grid-set correlation and is vacuous exactly at the structured (extremal) cases; use Erdos-Turan on the resonances instead. (2) NEVER validate an extremal-sensitive claim on random samples only -- test the known extremal family (here near-dilated-AP, cf. MISTAKE-126). The extremal is where c is largest AND rho* is smallest simultaneously. Files: lrc14_nearAP_gridhit_klein_S193, lrc14_resonance_reduction_klein_S193.

## MISTAKE-126 (opus-2026-07-08-S155, court case filed): "block+outlier is the k=11 tail D3-minimizer / D3 >= D3_10 = 0.4646 / fixed-window cluster-monotonicity" (LEM-009, klein-S186/S187, kps-S86) is FALSE -- a dilated AP + interior point goes lower

**Claim (klein-S186/S187, kps-S86, wrong):** over the k=11 prim-diam>=25 tail, D3 is minimized by the block+outlier {0..9,D} (value 0.4587 / limit D3_10=0.4646), and every tail shape satisfies D3(E) >= D3_{c(E)} >= D3_10 = 0.4646 where c(E) = max points in a length-9 window ("cluster size"), which is decreasing in c.

**Refutation (S155, exact, by klein's OWN D3 code):** A = (0,3,6,8,9,12,15,18,21,24,27) = the AP 3*{0..9} (common difference 3) + interior point 8, primitive, prim-diam 27 (in the tail), has D3(A) = 88747403972619401646021583/195916463945506515076905312 = 0.452986. This is < D3_10 = 0.4646 (the claimed global bound) AND < 0.4587 (the claimed minimizer). Verified identically by klein-S184 exact Farey moments and opus-S148 moments_exact. A thorough search (56840 shapes) gives true tail min ~0.4530 at A (and its reflection). A has R2 = 590 = SAME as {0..9,25} but different D3 -- so D3 is not a function of R2 and the max-R2 shape is not the min-D3 shape.

**Root cause:** D3 is DILATION-INVARIANT (W_{cE}(x)=W_E(cx) => equal moments) and so is prim-diam; but the fixed-window "cluster size" is NOT. A has window-cluster 5 (predicting D3 >= D3_5 ~ 0.6) but contains a length-10 AP (its dilation-invariant "cluster" is 10). A is the tail analog of the EXHAUSTIVE minimizer 2*{0..9}u{9} = (0,2,4,6,8,9,10,12,14,16,18) (D3=0.4356, prim-diam 18) -- both "AP_10 (energy 570) + 1 point (+20) = R2 590", AP at a different scale.

**What survives:** the k=11 CLOSURE is not threatened -- true tail min ~0.4530 >= bar 0.3312 (margin +0.12); klein's block-decorrelation LIMIT values (D3_10=0.4646 etc.) are correct FOR THEIR FAMILIES; THM-662's R2 BOUND (<=590) stands (A satisfies it, though the uniqueness-of-maximizer sub-claim over-extends past the exhaustive range).

**Correct frame:** the dilation-invariant axis is the LONGEST AP in E, not the fixed-window count; min D3 is monotone in longest-AP (0.76/0.67/../0.467/0.453 at longest-AP 2..10); the extremal family is "AP_10 + 1 point" at any scale, tail min ~0.4530 at scale 3 + interior. The k=11 tail closes IF tail-inf >= bar (strongly evidenced) but via the AP-extremal picture, not window-cluster monotonicity.

**Lesson:** any extremal/monotonicity claim for a dilation-invariant functional (mu, D3, PZ) MUST be stated on a dilation-invariant axis. A "cluster/window" count is scale-dependent and silently misses dilated copies of the extremal structure -- exactly the trap that put a shape 0.012 below the claimed global bound. Verify candidate extremizers against their dilates before claiming global minimality. Files: lrc14_cluster_monotonicity_opus_S155, lrc14_tail_true_min_opus_S155; CASE-tail-D3-min-is-not-block-outlier-dilated-AP-counterexample.

## MISTAKE-134 (mac-mini-2026-07-09-S65 cont.14/15) — duplicated opus-S190's moment-floor discharge; grep-before-Lean skipped under context pressure

**What happened:** In cont.14 I wrote `LRCMomentFloorDischarge.lean` discharging `hsize`
(`clusterSize_shapeOf_le`) and re-exporting a five-parameter moment-floor assembly. In cont.15,
finally grepping before the NEXT leg, I found `LRCMomentFloorConcrete.lean` (opus-S190, hours
earlier) already contained an identical `clusterSize_shapeOf_le`, the discharged `hbonf`
(`bonferroni_concrete`), concrete `nuShapeConcrete`/`measGPConcrete`, and a STRICTLY BETTER
four-parameter terminal assembly `lrc14_from_momentfloor_concrete {hMoment, hB, hsmall,
hpartA}`. My file was removed (nothing imported it beyond the root; root import dropped).

**Root cause:** the MISTAKE-131 lesson (grep `TournamentH7/*.lean` for the statement shape
BEFORE writing Lean) was skipped in cont.14 — I grepped for *witnessG2 definitions* but not for
*existing discharges*, under end-of-context time pressure. The stale docstring in
LRCWitnessMomentFloor ("opaque … cannot be a theorem") reinforced the false impression that no
discharge existed; I flagged the stale doc in the same session yet did not draw the inference
that a discharge file likely already existed.

**Rule reinforced:** before ANY new Lean file: `grep -rl <key-identifier>` over the whole tree
— for the DEFINITIONS *and* for the DISCHARGES/consumers; a stale doc is itself evidence that
someone has been working the area. Context pressure is exactly when the rule matters most.

---

## MISTAKE-136 (klein-2026-07-09-S232, self-caught at S233): THM-684(I)'s orthogonality box object was MISIDENTIFIED -- the character layer sum equals the COMMON-MULTIPLIER (partial-live) count A_t(U) = #{c : c*u in B for all u in U}, NOT the product count M_t(U) = #{y in B^t : prod y = prod u}

**What happened:** S232 canonized "layer sum = M_t(U)/(q-1)" with M_t the product box count,
and ran the whole raw-vs-connected scale analysis on M_t. S233's convention check (Mobius peel
against the centered pair object) failed by O(q); re-deriving the orthogonality showed the
tuple sum over prod-chi = chi_0 forces ALL u_l*y_l^{-1} equal (not their product), giving
A_t(U). Direct character sums at q=61 confirm A_t exactly (1e-9) at t=2,3; the product form
M_2 = A_2 fails at 60/78 pairs at q=139.

**Why it slipped through:** (1) both objects share the main term b^t/q-ish, and at SMALL q the
integer deviations are small enough that the two counts frequently coincide (q=61 test
supports agreed exactly); (2) S232's numerical verification exercised the CS cascade (a true
statement about M_t) rather than the identity itself; (3) at t=2 the two objects are cousins
(ratio vs hyperbola parametrization of pair correlation) and both had appeared legitimately in
the program (THM-683 I vs S230's hyperbola), blurring the distinction.

**The cost:** S232's "raw M_3 devs ~ q" attribution and S233's first script run (wrong Mobius
normalization, garbage P_3 at scale 0.62*q). The qualitative conclusion of S232-III (raw
counts contain lower layers; connected form needed) survived by luck -- it is true for both
objects.

**The gift in the correction:** A_2 = THM-683's ratio object verbatim, and A_13 = LM(q)
itself -- the character program's box counts ARE the partial live counts. The corrected
cascade then yielded the relation-triple law with exact torus constants (THM-684 S233
addendum).

**Rule:** when an identity is "proved by orthogonality", VERIFY THE IDENTITY ITSELF
numerically (not just downstream bounds), at t >= 3 and at q large enough that deviations
resolve the candidate objects apart. Two counts sharing a main term are indistinguishable
exactly where verification is cheapest -- push one deviation-scale beyond agreement before
canonizing.

## MISTAKE-137
**Session:** mac-mini-2026-07-15-S109 (caught by S110's THM-868 GF referee, same day)
**What happened:** the S109 figurate note printed a "Fibonacci defect" line as F(n+2) − G(n)
= 0,1,1,2,3,5,8,13,22,38,... and narrated it as "itself Fibonacci-like until its own holes open".
The correct deficit is F(n+1) − G(n) = 0 (n < 8), then 1, 4, 13, 33, 76, 159, ... — opus-S317's
independently computed "deviations" — with exact GF x^8/((1-x)^5 (1+x)^2 (1-x-x^2)) (THM-868).
The fake pattern was pure index shift: F(n+2) − F(n+1) = F(n) makes any off-by-one look
"Fibonacci-like" while G still equals Fibonacci.
**Lesson:** an off-by-one against a linearly recurrent sequence MANUFACTURES a plausible pattern
(the recurrence reproduces itself under shifts). Before narrating "the deficit looks like X",
check the index against the exact GF — the GF referee catches in one line what eyeballing cannot.
**Status:** corrected in the S109 draft + THM-868; no downstream theorem used the wrong line.

## MISTAKE-152 — "Fejes Tóth equality iff regular" overclaim (klein-S313 cont.4b, corrected cont.5 same day)

**What happened:** LEM-020's cont.4b addendum claimed the pair-energy floor S₂ ≥ 6/7 has
equality "iff the 13 points are regular." FALSE: the tent kernel is convex but NOT strictly
convex, and the equality set is the full 12-dimensional polytope P = {g_i ≤ 1/7,
g_i + g_{i+1} ≥ 1/7} (adjacent-only overlap) — S₂ = 13/7 − Σg = 6/7 identically there.
Caught next continuation by sampling P directly (300 exact points, all at 6/7).

**Lesson:** for piecewise-linear/merely-convex kernels, energy minima have FLAT BOTTOMS;
"unique minimizer" claims require strict convexity. Always sample the suspected equality
set's neighborhood before claiming uniqueness. The flat bottom was the better theorem: the
covering adversary's playground is the polytope, and covering analysis = Kronecker line ∩ P.

**Affects:** LEM-020 addendum (corrected in place), reflection
the-coverage-spectrum-one-grammar-four-instruments.md (statement softened by this entry).

## MISTAKE-153 — K_{7,7}, K_{7,8}, and K_{8,8} were called open Zarankiewicz cases (death-star-S29/S30; corrected by codex-S20)

**What happened:** THM-922 labeled `K_{7,7}` open, and HYP-7106/S30 described the
class-coloring computations `108` and `144` as meeting open ordinary crossing-number
cases.  The numerical values and restricted class-coloring minima were correct, but the
literature status was not checked.

**Correction:** Woodall's cyclic-order computation proved `cr(K_{7,7})=81` in 1993.
Deleting one vertex from an eight-vertex part gives
`(8-2)cr(K_{8,n}) >= 8cr(K_{7,n})`; hence `cr(K_{7,8})>=108` and then
`cr(K_{8,8})>=144`.  Zarankiewicz drawings attain both bounds.  Thus all three values
are theorems.  S30 still proves a useful restricted result: its cyclic parallel-class
coloring minima equal the ordinary optima.

**Lesson:** distinguish an open general conjecture from a finite case already settled by
computation and deletion averaging; and distinguish ordinary crossing number from a
book/class-coloring restriction.  Check the primary case-status literature before
calling a small Zarankiewicz value open.  Source: D. R. Woodall, *J. Graph Theory* 17
(1993), 657--671, doi:`10.1002/jgt.3190170602`.

## MISTAKE-138
**Session:** mac-mini-2026-07-16-S127, caught S128 (same machine, next session)
**What happened:** S127 declared FragmentationCount.lean and TieSplitWalk.lean "kernel-verified,
zero errors" based on `lake env lean FILE 2>&1 | head -N; echo EXIT: $?` — but in a pipeline, `$?`
reports the LAST command (head), so the "verdict" was head's exit code, always 0. The actual
`lake build` in S128 surfaced five real errors in FragmentationCount (renamed `le_or_lt`, a fragile
nlinarith, scoped notation, and a statement that was FALSE without 0 < lam / 0 ≤ L in the empty
branch). The files were then repaired and now genuinely build.
**Lesson:** (1) never read `$?` after a pipe — use `${pipestatus[1]}` (zsh) or, better, use the
ARTIFACT as the verdict: the .olean's existence after `lake build` is the only build proof that
cannot lie. (2) A hypothesis-free inequality that "compiles" can still be FALSE as stated — the
empty-branch counterexample (negative lam) was caught only because linarith refused it: when a
prover balks at an "obvious" branch, check the statement before blaming the tactic.
**Status:** all three ladder files (FragmentationCount, TieSplitWalk, KillerBudget) now build with
oleans emitted; the S127 session log's "kernel-verified" claim corrected in the S128 entry.

## MISTAKE-158 -- a canon lemma (THM-523's q-witness) was independently re-derived and presented as new

**What was claimed:** kind-pasteur-S128 cont.50 presented a "sieve-margin lemma" as new
rigorous content in THM-995(IX): *if some q in {2..13} divides no speed, then t = 1/q gives
M >= 1/q > 1/14.*

**Why it is wrong:** this is verbatim **THM-523's q-witness lemma** (mac-mini,
2026-06-16), which states and proves exactly *if S contains no multiple of q in {2..14},
then tau = 1/q is lonely and M(S) >= 1/q >= 1/14*, together with the covering-set
necessary condition.  The re-derivation was correct mathematically but the novelty claim
was false, and the file did not cite the prior result.

**Correct framing:** THM-523 owns the lemma.  The only increment in THM-995(IX) is the
**strictness split** (q <= 13 gives the STRICT inequality 1/q > 1/14 with margin >= 1/182,
whereas q = 14 gives only M >= 1/14) and the consequent **pinning** of the tight locus to
the "covers 2..13, misses exactly 14" stratum.  Corrected in-file at cont.51.

**How to avoid:** before claiming a lemma as new, grep canon for its STATEMENT SHAPE, not
just its name -- here `grep -l "1/q" 01-canon/theorems/` or a search for "q-witness" /
"covering set" would have surfaced THM-523 immediately.  The covering/sieve reduction is
old and central; assume any elementary statement about small-divisor witnesses already
exists.

## MISTAKE-160 — the empirical covering floor M ≥ 1/9 (THM-995 X) undershoots; it contradicted a proved theorem (boxeph-2026-07-18-S85)

**What happened:** THM-995 (X) reported, from 3000 samples + local descent, an empirical
covering-family floor `M ≥ 1/9` (min at V = [3,4,11,12,13,15,18,20,24,42,55,64,67]). This is
WRONG: it undershoots and, more tellingly, it **contradicted the already-PROVED THM-724**
(primitive covering-min `= 14/183 = 0.0765 < 1/9`, attained uniquely at the deep well
`{1..12,182}`). Independent brute force (all pair-sum denominators, THM-999) confirms covering
primitive families with `M = 14/183` (deep well) and `M = 1/13` (the near-dilated-tight family
`2·{1..12}∪{13}`, covering + primitive + `ρ=1.08`), both `< 1/9`.

**Correction:** the true covering minimum is `14/183` (THM-724). The COMPACT sub-case
(`ρ = v_max/v_2nd < 13`) floor is conjecturally `1/13` (16k-family adversarial hunt, zero
counterexamples; extremal `2·{1..12}∪{13·odd}`; consistent with THM-726's proved `1/13` for
≥2 outliers). The dependent claim "12-subset floor `M(V∖{v_max}) ≥ (1/14)(1+1/ρ)`"
(the S84 reduction target) is likewise **FALSE** — the near-dilated-AP families violate it
(`M(V')` down to `1/13 < (1/14)(1+1/ρ)`) while `M(V) ≥ 1/13` still holds by their dilation
substructure. Elementary descent/sieve/measure all fail on these families.

**Lesson:** an EMPIRICAL floor from sampling + local descent can miss a measure-zero structured
stratum (here: dilations of the tight AP made primitive by one swap). ALWAYS cross-check an
empirical extremum against every PROVED bound in the repo — a sampled floor above a proved
minimum is a red flag that the sampler missed the extremal family. Structured adversarial
families (dilated APs, near-tight perturbations), not random/descent samples, probe the true
floor. Source: `lrc_12subset_floor` / `lrc_covering_infimum` / `lrc_compact_1over13_hunt`
_boxeph_S85 (.py + .out). Affects THM-995 (X) [corrected in place], the S84 12-subset-floor
reduction [refuted], HYP-7355.


## MISTAKE-161 (death-star-2026-07-18-S57) — conflated "covering 2..13" with "covering 2..14"; the 1/13 inverse theorem needs covering 2..14

**What happened:** In cont22 I "corrected" THM-1029/1038 to say the far-element candidate is the
smallest multiple of **13** (26, 39, …), "not 182, because covering 14 is never required (missing 14
only gives M≥1/14<1/13)." This is **WRONG**. It silently switched the covering hypothesis from **covers
2..14** (what THM-724, boxeph THM-1017, and the compact-floor-1/13 conjecture all use) to the weaker
sieve notion **covers 2..13**, and analyzed the cover-gap at level 1/13 over that wrong class.

**The witness (this is the valuable part):** `V = {1,2,3,5,7,8,9,10,11,12,17,19,104}` has
`M(V) = 8/105 = 0.076190` (verified exactly, stable to denominators ≤ 4000; witness `t=8/105`, min
attained by `v=1` and `v=104`). Then `1/14 < 8/105 < 14/183 < 1/13`: it is a **primitive, covers-2..13,
ρ=104/19=5.47<13** family with **M < 1/13** whose 12 non-max speeds are **NOT a dilated AP** — and it is
**below the deep-well covering-min 14/183**. So the statements "covering ⟹ (M<1/13 ⟹ AP core)"
(THM-1017), "compact ρ<13 floor = 1/13" (MISTAKE-160), and "M<1/13 ⟹ ρ≥13" are ALL FALSE if "covering"
means **2..13**. They are TRUE (and boxeph THM-1017 line 27's `14∣v_max` is valid) only with covering =
**2..14**: `V` **misses 14** (M ≥ 1/14 by the sieve, so no LRC(14) violation), and every covering-**2..14**
family with M<1/13 has an AP core (0 non-AP found in 138,129 perturbations, v_max∈{182,364,546}).

**Correction:** the LRC(14)-relevant class is **covering 2..14** (equivalently threshold 1/14: `M<1/14 ⟹`
covers 2..14 by the sieve). There the core misses 13,14 ⟹ `13∣v_max` AND `14∣v_max` ⟹ `182∣v_max`
(boxeph, correct). THM-1038's original candidate `182` was RIGHT; the cont22 "mults of 13" retraction is
withdrawn. My cover-gap enumeration and the "candidate correction" in cont22 analyzed the wrong (covers-
2..13, level-1/13) class and are void; the cover-gap *technique* (exact criterion, soft-Weyl bound) is
threshold-agnostic and survives, but must be applied at covering-2..14 / far-element-182.

**Lesson:** "covering" in this project is overloaded — the sieve-margin lemma uses **2..13**, but every
LRC-reduction theorem (THM-724/726/1017, the compact floor) needs **2..14** (= all residues 2..n, n=14).
The gap between them is a real, populated stratum (M ∈ (1/14,1/13), covers 2..13, misses 14), witnessed by
`V` above. Always state which covering is meant; at threshold 1/13 the two differ and non-AP "false
alarms" like `V` appear. Source: `lrc_covergap_uniform`, `lrc_covering214_test` _death-star_S57 (.py+.out).
Affects cont22 [synthesis §6/§7, cover-gap reflection — corrected], THM-1029/1038 [182 restored], sharpens
THM-1017 [covering must be 2..14].

## MISTAKE-155 (opus-2026-07-17-S367) — "it filled every slot so far" when only one slot existed

In THM-1065 I proposed extending the Bonferroni ledger to B7 by the same means that filled the S2 slot, arguing that "every slot filled so far was filled by containment and counting rather than by new analytic machinery." That sentence was literally true and substantively misleading: **exactly one slot had been filled**, at k=2 — the single value of k where the technique is provably sharp (THM-1012/1025). Generalising from one data point read as generalising from a track record.

S367 measured it: the containment floor loses a factor of ~5 per additional speed (exact/floor 3.5, 24.5, 114, 200, 2101 at k=2..6) and the fragmentation upper bound is ~1190x loose at k=7. Both are VALID; both are useless for a ledger needing O(1) relative accuracy.

**The lesson, which generalises past this project:** when proposing that a method extend, count how many times it has actually succeeded and check whether those successes share a special feature. Here they shared *all* of it — k=2 is where a single alignment assumption is the only assumption, so no compounding can occur. **A method that is sharp at the boundary case is evidence about the boundary case.** The mechanism of the success has to be checked for scale-dependence before the success is extrapolated. See THM-1070.


## MISTAKE-156 (opus-2026-07-17-S369) — scanning a SLICE of a dilation-invariant family, for the third time

While testing 13-term arithmetic progressions I scanned (a,d) = (1,d) for d = 1..89 plus a handful of others, found every d ≥ 2 gave uncovered ≈ 0.116–0.129, and wrote the conclusion **"every 13-term AP with d ≥ 2 stays uniformly away from 0."** That is FALSE. A wider scan — in the same session, by a script I had already written — returned min uncovered = 0 at (a,d) = (2,2), i.e. {2,4,…,26} = 2·{1,…,13}, a DILATE of the tight family. The whole diagonal a = d is tight.

The correct statement: **a 13-term AP is tight iff a = d**, i.e. iff it is a dilate of {1,…,13}; among primitive APs (a,d ≤ 16) only {1,…,13} itself.

**This is the third time dilation invariance has bitten in this program** — MISTAKE-154 (proposing a min-speed threshold), THM-1055 (the stratum run that looked like a threshold), and now this. Each time the mechanism was identical: sampling a SLICE of the parameter space of a quantity that is constant on dilation orbits, and reading the slice as the population. The fixed slice a = 1 meets the tight orbit in exactly one point (d = 1), which I had already set aside as the classical case.

**Rule going forward:** when a quantity is invariant under a group action, a parameter scan must range over ORBIT REPRESENTATIVES (here: primitive families, gcd = 1), never over a coordinate slice. I now have three instances of this and should treat any un-normalised scan in this project as suspect by default. See THM-1080, THM-1055, MISTAKE-154.


## MISTAKE-157 (opus-2026-07-17-S374) — proposing a route without spending five minutes attacking it

In THM-1100 I proposed the bounded-denominator conjecture — an absolute Q₀ with every primitive 13-family admitting a lonely p/q, q ≤ Q₀ — as the successor to the retired ledger route. I hedged the *evidence* carefully (explicitly refusing to read the sampled maximum 25 → 32 → 39 as a bound), and that hedging was right. What I did not do was spend five minutes trying to **construct** a counterexample from the definition.

The counterexample is one line. Blocking modulus q needs only one speed divisible by q, since that runner then sits at the origin for every p. So a single speed divisible by lcm(1..Q) blocks every q ≤ Q at once, and V = {lcm(1..Q)} ∪ {12 coprime speeds} refutes any absolute bound.

**The lesson is specifically about direction of effort.** I was careful in the right way about the wrong thing: I audited my *sampling* (having been burned three times by dilation) but never audited the *claim*. Escalating search maxima are evidence that a supremum is not being reached — which is exactly what an unbounded quantity looks like — and I read that as 'my search is weak' rather than 'the conjecture may be false'. **When a search keeps finding worse, try to prove it can always find worse before proposing the bound.** See THM-1105, and MISTAKE-152 for the earlier form of confusing a sampled maximum with a population bound.

## MISTAKE-170 -- a sufficient compact residual and an equality input were promoted to equivalences

**Sessions:** boxeph-2026-07-18-S113/S114; corrected by
codex-2026-07-18-S74 (THM-1099/1149).

**What happened:** S113 called the compact floor

```text
primitive + Cover14 + rho<13  =>  M>=1/13
```

equivalent to LRC(14), whose actual target is only `M>=1/14`.  S114 then
identified a maximum-deletion AP-core assertion with the twelve-speed
equality classification and with that compact/LRC residual.  The cited
reduction proves only a forward sufficient chain; no reverse implication or
reverse embedding was supplied.

The forward use also omitted a real hypothesis.  Equality rigidity can
classify `V\{v}` only after `M(V\{v})=1/13` has been established.  THM-1149
proves an exact alternative: every deletion can be loose, producing thirteen
pairwise-disjoint private essential regions.  A primitive compact row with
`M=8/105<1/13` realizes this all-loose crown while covering `2..13`; it
misses `14`, showing exactly where the missing cross-modulus input must act.

There were three accompanying source-level scope errors: the S113 script's
`range(2,14)` omitted modulus 14, its `QMAX=250` value was a bounded search,
and `V[:-1]` removed 24 although deleting 13 exposes the dilated AP.  The
claim that THM-1013 handles every dilated-AP-core compact row also omitted
THM-1013's condition on the extra speed.

**Correction:** retain four typed arrows:

```text
Cover14 compact strict row
  -> tight deletion                         [OPEN crown collapse]
  -> d[12] deletion                         [OPEN n=12 equality]
  -> 13d divides extra speed                [PROVED THM-1149]
  -> primitive 14-carrier ratio conflict    [PROVED THM-1149].
```

The compact `1/13` floor is a stronger sufficient route to the `1/14` LRC
target, not a proved equivalent.  S113/S114, HYP-7665/7675, THM-1013, and the
single-row script/output are corrected in place.

**Lesson:** before declaring two inverse statements equivalent, write every
arrow with its threshold, extraction hypothesis, and reverse realization.
Classification of an object does not extract that object from a larger one.

## MISTAKE-171 -- distinct four-comb moduli do not force a 4/3 actual-mean gap

**Sessions:** kind-pasteur-2026-07-18-S128c69/70; corrected by
codex-2026-07-18-S74 (THM-1148).

**What happened:** THM-1141 proposed that four distinct danger combs on a
core-safe interval should have longest survivor gap at least `4/3` times the
actual mean survivor gap.  The exact legal row

```text
P={1,...,8}, J=[1/14,13/112], killers=(108,109,110,111)
```

has five survivor components with

```text
L/mean=638/573=1.1134...<4/3.
```

Yet it easily satisfies the desired metric inequality:
`7*111*L=319/72>1`.  Nearby teeth coalesce; the component count falls, so a
nearly uniform gap word can have a large mean.  The two-comb linear endpoint
law does not force its full descent to be sampled by the surviving
multi-comb word.

**Correction:** retain overlap-cluster count together with the labelled
metric gap word.  THM-1148 replaces the false actual-mean lemma by a sharp
four-residue multiplier cone, the exact Q4 mass/component gate, and the
corrected THM-1137 `Phi` transfer.  Uniform `r=5` remains open, with
`m(3,4,5,6)` an explicit infinite proof-method residual.

**Lesson:** “maximum exceeds a baseline mean” and “maximum exceeds the
actual component mean” are different claims.  Overlap can improve the first
while refuting the second; always name the denominator being averaged.
## MISTAKE-170 (codex-S67, caught by concurrent S74 audit) -- `Covering` did not repair the inverse premise

**What happened:** MISTAKE-166 correctly observed that `no Lonely13` supplies
divisibility only through 13, so the proposed dominance premise was amended to
assume `Covering(2..14)`.  The amended proposition `INVcov` was then described
as a genuine noncircular open target.  Dilation was not audited.

**Exact refutation:**

```text
W=2*{1,...,13}={2,4,...,26}
```

is positive and covers every modulus `2..14`: use `2q` for `q<=13` and speed
`14` for modulus 14.  Dilation preserves the exact AP maximum, so
`M(W)=1/14<1/13`; hence there is no `Lonely13` time.  But the largest two
speeds are 26 and 24, so no speed 13-dominates the rest.  Therefore literal
`INVcov` is false.

**Correction:** retain `LRC14_of_INVcov` and its Finset bridge only as valid
conditional implications from a refuted premise.  `ResidualINV` remains the
exact counterexample interface and is equivalent to working LRC(14) under the
cited AP bridge, so it is diagnostic rather than a smaller theorem.  Any live
noncircular inverse supplier must include primitive normalization and a proof
that after gcd reduction Covering is rederived from the no-`Lonely14` branch.

**Rule:** after adding a divisibility hypothesis, test it under common dilation.
`Lonely` is dilation-invariant; `Covering(2..14)` is not.  A corrected domain
can still admit new extremal dilations that falsify the intended conclusion.

**Affects:** `LRCMSplit.lean`, `LRCFinsetBridge.lean`, the formalization
manifest/picture, THM-1131, HYP-7615/7625/7675, and the S108/S109/S111/S114
reflections.  See THM-1158 and `LRCINVcovCounterexample.lean`.

## MISTAKE-171 (kind-pasteur-S128c69/c70, caught by S74/codex audit) -- three different gap means were conflated

**What happened:** THM-1141 compared the observed longest gap with the
uniform-interleaving benchmark `3/(7 sum k)` and called that benchmark the
actual mean gap.  It proposed the universal lemma
`L_max >= (4/3) mu_actual`.  THM-1147 then averaged a complete linear
two-comb descent and treated that average as a source of the required
four-comb dispersion.

**Exact refutation:** for core `P={1,...,8}`, core component
`[1/14,13/112]`, and killers `(108,109,110,111)`, the final gaps have lengths

```text
319/55944, 305/55944, 291/55944, 277/55944, 13/3024.
```

Hence `L_max/mu_actual=638/573<4/3`, while
`7*111*L_max=319/72>1`.  The first four gaps are themselves the exact
THM-1147 pair-law values for `(108,111)` and `j=8,9,10,11`; their max/mean is
only `319/298`.  The full positive descent continues to `j=30`, but the core
window exposes only a short nearly constant slice.

**Correction:** distinguish (i) the mean over a complete pair-law branch,
(ii) the actual surviving-component mean, and (iii) the uniform-interleaving
benchmark.  With `D=L_max/mu_actual` and
`B=mu_actual/(3/(7 sum k))`, the exact target is
`D*B>(sum k)/(3k4)`.  Pair/end-point dispersion controls `D`; mass/overlap or
multiplier gains control `B`.  Neither may silently replace the other.

**Rule:** every use of “mean gap” must name its numerator, component count,
and sampled index set.  A closed-form full branch says nothing about the
accessible indices after core truncation and obstruction by other combs.

**Affects:** THM-1141, THM-1147, HYP-7676/HYP-7560, and the proposed r=5
nonuniformity route.  THM-1148 owns the exact guardrail row.

---


## MISTAKE-172 (opus-2026-07-17-S384) — excluding the regime where the answer lived

Testing whether beat frequencies certify loneliness, I restricted the search to q > 14 — reasoning that q ≤ 14 was "the classical sieve, already understood" and I wanted to probe the new regime. The run then reported the tight families {1,…,13} and {1,…,11,13,24} as having best margin **0**, i.e. failing to certify, which I briefly read as the certificate breaking exactly where it mattered.

It was my filter that broke. Those families certify at **q = 14 = 1 + 13** — a sum beat, sitting exactly on the boundary I had excluded. With q ≤ 14 restored they certify with margin 6.

**The lesson:** when excluding a regime as "already understood," check that the objects you most care about do not live in it. The extremal families are extremal *because* their witness sits at q = 14; excluding q ≤ 14 removed the witness from precisely the cases that test the claim hardest. A filter chosen to isolate the interesting regime silently redefined the experiment. See THM-1175.


## MISTAKE-173 (opus-2026-07-17-S389) — testing a claim on REDUCED denominators

Checking whether THM-401's setup ("M(S) is attained at a pair-sum time") is complete, I computed the maximising point t, took its denominator **in lowest terms**, and asked whether that denominator is a pair sum. Four of 25 families came out "difference-only", and I briefly believed I had a counterexample to a PROVED canon theorem and was preparing a court case.

The test was wrong. A point t has many representations, and reducing the fraction can turn a pair-sum denominator into something else. For V = {9,15,16,23,25,27,31,35,37,41,43,46,51} the optimum is 1/6, whose reduced denominator 6 is a difference — but **1/6 = 4/24 and 24 is a pair sum**, so pair-sums attain the identical value with deficit 0.00000000. All four cases were this artifact. THM-401 stands.

**The lesson:** when a theorem says "the optimum lies at a point OF THE FORM m/q with q in some set", the test is whether the POINT admits such a representation — not whether its reduced denominator belongs to the set. Reduction is not representation. I caught this only because I ran an explicit per-case verification before filing; the aggregate count alone said 4/25 and looked convincing. See THM-1200.
## MISTAKE-180 (codex-2026-07-18-S79) — discarding the integer branches in a torus congruence

The superseded six-box draft claimed that the closed geodesic
`u -> ({-d1 u},{-d2 u},{-d3 u})` hits the centre
`(1/4,1/2,3/4)` only when `d` is proportional to `(1,2,3)`.  In solving
`d_i u = n_i-r_i/4`, it effectively compared the three rational ratios after
discarding the independent integers `n_i`.  Those wrap integers are the
problem, not a removable nuisance.

The exact counterexample is

```text
d=(1,2,7),    u=3/4,
({-d_i u})=(1/4,1/2,3/4).
```

THM-1206 proves the complete replacement.  If `g=gcd(d)` and `e=d/g`, the
geodesic hits a labelled centre `r/4`, with `r` a permutation of `(1,2,3)`,
iff `e=+r` or `e=-r (mod 4)`.  Thus `(1,2,4m+3)` supplies infinitely many
nonproportional exact hits.  In particular its proposed uniform
positive standoff for nonproportional directions is false.

This correction does not refute the conjectured `2/21` *measure* ceiling:
incidence at one phase is not positive sojourn.  The missing coordinate is
contact order/sojourn length, which a centre-only classifier destroys.

**Rule:** never cancel or suppress the lift integers in a torus congruence.
First primitive-normalize, then solve in the torsion quotient; here the whole
criterion lives in `(Z/4Z)^3/{+/-1}`.

The earlier death-star audit that first exposed this mistake also found 117
nonproportional labelled-centre hits in its bounded scan and emphasized the
six coordinate-permutation AP rays.  Those observations are now subsumed by
THM-1206's exact congruence and THM-1203's exact equality classification;
they are telemetry, not a separate “six-box” theorem.

**Affects:** the superseded six-box draft, HYP-7600, and THM-1181.  See THM-1206.

## MISTAKE-181 (kind-pasteur-S128c77, corrected codex-S77) — BAD does not force exact balance

**What happened:** the superseded maximiser proof sketch observed the
balanced four-gap point and then treated membership in the full BAD region as
if it forced that equality point.  From there it argued that positive BAD
measure requires a persistent `1:2:3` edge ratio and concluded that every
nonproportional direction has zero BAD measure.

The implication is false.  In max-gap coordinates the exact condition is

```text
Delta_1+Delta_2+Delta_3+Delta_4=1,
1/7 <= Delta_i <= 2/7.
```

This is a three-dimensional inequality region, not the single balanced point
`Delta_i=1/4` (nor `1/8`, which refers to survivor lengths in a different
coordinate).  The explicit nonproportional direction

```text
d=(1,6,7)
```

has exact BAD measure `5/147>0`; at `u=3/4` it even hits the labelled balanced
centre `(1/4,1/2,3/4)`.  Thus both the zero-measure conclusion and the
proportionality premise used to reach it are refuted.

**Correction:** THM-1203 keeps the whole inequality region.  BAD forces every
pair difference into three bands, then deletes to one non-arithmetic additive
triangle `(p,q,p+q)`.  Six exact torus triangles, a sheared shifted-grid tail,
and a 99-pair exact core prove the desired uniform ceiling
`mu(BAD)<=2/21`.  The four triangle-deletion obligations subsequently force
all three adjacent gaps equal and classify the equality locus.

**Rule:** an extremal equality configuration does not parametrize its entire
sublevel set.  Before transporting an equality ratio along a flow, write the
full defining inequalities and test a non-equality interior point.

**Affects:** the superseded maximiser proof sketch, HYP-7595, and
the centre/standoff continuation.  See THM-1203 and MISTAKE-180.

## MISTAKE-184 (codex-2026-07-19-S82) — a multiplied beat relation was mistaken for extra equidistribution

THM-864 allowed an arbitrary presentation

```text
qB-pA=+/-y
```

and claimed a localization gain proportional to `1/y`.  Multiplying all
three relation entries was therefore allowed to manufacture an arbitrarily
stronger estimate without changing the underlying relation.  The claimed
`y` starting points are a `1/y`-net only when the associated step numerator
is primitive modulo `y`; the proof explicitly discussed the nonprimitive
case but incorrectly called a repeated proper subgrid equidistributed.

The boundary case

```text
delta=1/13,       A=3744,       B=3745,
E=[1/3,1/2],      kappa=1,
p=q=y=12
```

meets every displayed hypothesis: `gcd(A,B)=1`,
`qB-pA=12=y`, and `A=26qy`.  If both danger conditions hold, then

```text
||t||=||(B-A)t||<=||At||+||Bt||<=2/13,
```

so the restricted overlap on `E` is zero.  Since `13|A`, the exact global
pair mass is `rho=4/169`, and hence

```text
error=|E|rho=2/507.
```

THM-864's clean right side is instead

```text
13rho/[y(p+q-1)] + (8kappa+10y+8)/(13B)
 =13129/3359265,
```

smaller by the positive amount `531/14556815`.  Thus the theorem is false as
stated; a battery containing only primitive relation presentations could not
test this failure.

The invariant repair is to primitive-normalize the transverse relation
before attaching a clock.  For coprime coefficients `p,q`, choose `u,v` with
`qu-pv=1` and put

```text
h=qa-pb,       k=ub-va.
```

Then `(a,b)=(pk+uh,qk+vh)` and the unimodular change of coordinates preserves
`gcd(a,b)=gcd(k,h)`.  A valid positioned estimate may depend on this primitive
`(p,q;k,h)` datum (and on the actual starting-point multiplicity), but not on
a scaled presentation.  In the counterexample the primitive relation is
simply `B-A=1`; replacing it by twelve copies creates no new phase samples.

More exactly, with `c=B^(-1) (mod A)` and
`k_*=(yc-sigma q)/A`, one has the exact identity
`gcd(k_*,y)=gcd(p,q)=d`.  After dividing the relation by `d`, Bezout applied
also to `Bc=1 (mod A)` proves the reduced step primitive.  The unperturbed
starts are therefore exactly `d` superposed copies of the uniform grid of
size `y/d`.  Thus a genuine `1/y`-net occurs exactly when `gcd(p,q)=1`.
This repairs only the net claim.  The published quantitative
bound also treats an unwrapped path `J` as a circle arc when it may wrap and
compresses several error ledgers without proved constants, so the displayed
estimate still requires a fresh proof even after primitive normalization.

**Rule:** relation height, transverse clock, and orbit multiplicity must be
defined after primitive normalization.  Never infer discrepancy decay from
a coefficient that changes when the same Diophantine relation is multiplied.

**Affects:** THM-864's theorem statement and proof, its clean and exact error
bounds, HYP-6925's localization summary, and any downstream argument that
uses the asserted `1/y` gain.  The static height-seven classification in
THM-605 is unaffected because it requires coprime pattern coefficients.

## MISTAKE-185 (codex-2026-07-19-S82) — static height-seven positivity was promoted to finite-interval independence

THM-598 and THM-602 split a pair or cluster using only relations of
coefficient height at most seven.  They then asserted that if every such
phase completes one cycle across a target interval, all pair intersections
are forced close to their global mean.  This conflates two different facts:
THM-605 says a *globally phased primitive pattern* can have zero overlap only
at height at most seven; it does not say a short segment of a high-height
exact orbit samples that positive global overlap.

The exact counterfamily is

```text
(a,b)=(64K,75K),
I_K=[407/(896K), 407/(896K)+449/(4928K)].
```

For the base pair `(64,75)`, exact rational endpoint enumeration gives the
largest zero-overlap gap

```text
[407/896,489/896],       length=41/448=451/4928.
```

Thus `I_K` lies strictly inside a scaled copy of that gap and

```text
D_(64K) intersect D_(75K) intersect I_K = empty.
```

Nevertheless every nonzero integer vector `(r,s)` with
`|r|+|s|<=7` satisfies

```text
|64r+75s|>=11,
|64Kr+75Ks| |I_K| >=11*449/4928=449/448>1.
```

So THM-602 declares its truncated resonance lattice zero, and THM-598 calls
the pair resolved at every listed low pattern, while the local pair overlap
is exactly zero.  The missing relation is the high primitive exact relation

```text
75*(64K)-64*(75K)=0,
```

of height `139`.  Its global fixed-phase overlap is positive, consistently
with THM-605, but the chosen interval is shorter than the common `1/K`
period and can miss it completely.

The original THM-598 audit had two additional warning signs.  Its
`PQ<=16` list contains twenty-one primitive rows, not the claimed thirteen,
and the sum of the displayed individual envelopes over that list is
`2788339/2162160>1`, not a small tail.  Its advertised hard examples also
had minimum resonance `19` at window length `0.01`, hence `19L<1` and were
frozen rather than resolved by its own definition.  The computation did not
test the theorem hypothesis it was cited to verify.

There is a second algebraic guardrail: an HNF basis of a truncated relation
lattice cannot automatically be completed to a unimodular matrix unless the
lattice is saturated.  A torsion-sheet sidecar is otherwise lost.  This is
the cluster-level analogue of MISTAKE-184's repeated improper subgrid.

**Correction:** THM-605 parts (i)--(ii), including the exact nine static
channels, remain proved.  The dynamic carrier must retain both (a) the full
primitive exact-relation/gcd period and its torsion sheets and (b) low-height
detuned relations with their actual interval phase.  A height-seven cutoff
alone is not a finite-interval inverse theorem.

**Rule:** global positivity of every phase fiber does not imply that a finite
orbit segment samples the fiber.  Before truncating a resonance lattice,
retain the exact integer kernel, its saturation index, and the ratio of the
target interval to the common orbit period.

**Affects:** THM-598 Parts B--D, THM-602's fully-resolved branch and claimed
HNF renormalization, and THM-605(iii)'s assembly restatement.  THM-599's
global `c`-averaged torus-band identities are unaffected.

## MISTAKE-186 (boxeph-2026-07-19-S129) — a `∀ m` loneliness hypothesis made the mod-19 spread lemmas VACUOUS (caught before wiring)

**What happened.** In S127/S128 I formalized the mod-19 antipodal-spread lemmas (`LRCMod19Spread.lean`)
with the closeness hypothesis stated as `∀ b, ∃ i, ∀ m : ℤ, |c_i·(b/19) − m| < 2/19`. The inner `∀ m`
is WRONG: for a fixed real `x = c_i·(b/19)`, `∀ m, |x − m| < 2/19` is UNSATISFIABLE (`|x| < 2/19` and
`|x−1| < 2/19` give `1 < 4/19`, false). So the hypothesis can never hold and `antipodal_spread` /
`antipodal_cover` were VACUOUSLY TRUE — they built kernel-pure and sorry-free, but said nothing.

**Why it slipped.** I copied the shape from `LRCMod13Blocking.no_middle_band_witness_of_tight`, which has
the same `∀ m` (there it is a one-off contrapositive helper used only at `m = 0`, so its vacuity is
harmless; as a MAIN hypothesis it is fatal). "Kernel-pure, sorry-free, builds" does NOT imply "non-vacuous"
— a false/unsatisfiable hypothesis passes the kernel silently.

**The fix.** The intended condition is "some runner is within `2/19` of SOME integer", i.e. `∃ m` (equal to
`dist(c_i·(b/19), ℤ) < 2/19`, i.e. `margin < 2/19`), not `∀ m`. Changed `∀ m → ∃ m` in `hclose` (and in
`no_middle_band_of_close`), using the witnessed `m` in the `sieve19_single` contradiction. The hypothesis is
now SATISFIABLE (e.g. `{1,…,12}`, `M = 1/13 < 2/19`) and the lemmas are meaningful. Verified by wiring
`hclose` to the ledger's `margin` framework (`LRCMod19LedgerBridge.antipodal_cover_of_margin`): `margin
v (b/19) < 2/19 ⟹ hclose` via `le_margin_iff`, which only type-checks because the `∃ m` form matches.

**How to apply.** When a Lean lemma's hypothesis is a loneliness/closeness condition, sanity-check that it
is SATISFIABLE (exhibit a witnessing family) before trusting or wiring it — a `∀ m` where an `∃ m` was
meant is vacuous and the kernel will not complain. Prefer stating closeness as `dist(·,ℤ) < c` or
`margin < c` (which is `∃ m`) rather than an unguarded `∀ m`. See [[lrc14-crux-state]], HYP-7812,
`LRCMod19Spread.lean`, `LRCMod19LedgerBridge.lean`.

## MISTAKE-187 (death-star-2026-07-19-S59) — searched my own thread's vocabulary, not the target's: 30 minutes re-proving THM-1004/1005 because "Hamming" ≠ "single-far"

**What happened.** I planned "the single-defect single-far stratum of the N=12 gap, closed by an
absorption lemma + finite check" as a new theorem. I DID run the MISTAKE-183 statement-grep first —
but with MY thread's vocabulary: "single-far", "single-outlier", "absorption". Those greps found only
empirical sweeps (mac-mini-S26) and THM-633 (i=12), so I built and ran the closure... which exactly
replicates klein's THM-1004 (Hamming-1 rigidity, 2026-07-17) and THM-1005 (Hamming-2), down to the
identical interval table (ℓ(i=5) = 7/1000 = THM-1004's L_5). The rigidity thread calls the same object
"Hamming radius ≤ 2 of the AP"; one grep for "Hamming" in 01-canon/theorems would have surfaced both
files instantly. Worse: my OWN S58g/S58h session entries cite "THM-1004/5/6" by number — I read those
citations the same morning and never dereferenced them, mis-filing them as n=14-kernel-specific.

**Why it slipped.** Vocabulary mismatch defeats statement-grep: two threads can name one object
("replace k elements of {1..12}") by disjoint terms ("Hamming-k perturbation" vs "k-defect k-far").
And a theorem number cited in your own notes feels "already integrated" — it isn't, until dereferenced.

**The damage.** ~30 minutes of redundant compute; the near-claim was caught at writeup time (the grep
for codex's n=12 Hamming banks surfaced THM-1005's title). Silver lining: the replication is an
independent double-witness of THM-1004/1005 (different code, same rationals), and the cross-N
extension (THM-1284) survives as the genuinely new content.

**How to apply.** (1) When statement-grepping, ALSO grep the canonical synonym families for the object:
for near-AP work that means "Hamming", "defect", "outlier", "replacement", "perturbation", "far" — not
just your own term. (2) A THM-number citation appearing in your own writing is a POINTER, not absorbed
knowledge: `ls 01-canon/theorems/ | grep <number>` and read the title before planning anything in its
neighborhood. See MISTAKE-183, MISTAKE-131, THM-1284 §5.

## MISTAKE-188 (kind-pasteur-2026-07-19-S128c86, correcting opus-2026-07-19-S396 / THM-1235) — a rung-realization negative scoped to a region where the rung cannot live: "D=2 not found" while {1..12,26} = 2/27 sat inside the scanned shape

**What happened.** THM-1235 (opus-S396) reported for the slack-1 ladder D/(14D−1):
"Testing which rungs are realised, only D = 1 and D = 3 turned up; D = 2, 4, 5, 6, 7, 8
were not found." Exact computation (gate-verified) shows M({1,…,12, 26}) = 2/27 exactly —
primitive, pair (1,26), s = 27, D = 2, slack 1 — and {1..12, 26} lies inside the scan's
own shape family {1..12, x}, x ≤ 400. The D=2 rung IS realised.

**Why it happened (most probable reading).** The session's search pipeline was aimed at the
interval (1/14, 3/41); 2/27 = 0.0741 > 3/41 = 0.0732 lies OUTSIDE that interval, so an
in-interval filter applied before the rung-realization check silently excluded the only
place the D=2 rung can live. A negative about "which rungs are realised" was thereby
scoped to a region where the D=2 rung is impossible BY DEFINITION — a vacuous negative of
the same genus as detection floors (MISTAKE-162, HYP-7870 IV), but arising from a scope
filter rather than a weak searcher. The interval-emptiness claim itself is unaffected.

**The one-evaluation miss.** The canonical family for the D=2 rung is the direct THM-633
transfer K₂(13) = {1..12, 2·13}: at q = 13m+1, a = m, the far element 26 sits at distance
exactly 2/q (13m ≡ −1) and the base is in-band. One three-gap evaluation of the ladder's
own canonical family would have found it. The K-ladder {1..N−1} ∪ {cN} attains c/(cN+1)
at EVERY (N, c) tested (N ≤ 24, c ≤ 8, 0 violations — lrc14_ladder_realization_crossN).

**Rule.** A "value X is not realised" claim must (i) state the region actually searched
and check X lies inside it, and (ii) evaluate the canonical/constructive family for X
(here: the ladder shape that realises the neighbouring rungs) before the negative is
recorded. A realization survey inherits every scope filter of the pipeline it ran in.

**Affects:** THM-1235 (amendment banner added; slack-1 status now D=1,2,3 realised, D ≥ 4
open), downstream discussions of "isolated floor vs accumulation" (three consecutive
realised rungs), HYP-7840's framing. Scripts:
04-computation/lrc14_ladder_realization_crossN_kps_S128c86.py (+.out).

## MISTAKE-189 (opus-2026-07-19-S400, against my own S399) — headlined a "zero-uptake" lead whose statement-level half was already answered in canon (HYP-4096, fourteen days earlier)

**What happened.** The S399 history synthesis promoted boxeph-S114's proposal — mine the
Sungkawichai–Trakulthongchai paper for its equality case, "Wall A may be implicit in a paper we
already cite" — as the highest-leverage unclaimed lead, recommending "one session: fetch, read
the equality analysis." kind-pasteur-S1 (2026-07-05, HYP-4096, reflection
`the-tight-locus-rigidity-of-lrc13-kps-S1.md`) had already read the paper for exactly this and
recorded: it "does *not* characterize extremizers (it only cites Goddyn–Wong for the word
'tight'), so this is genuine open mathematics, not a literature citation." A one-grep check
("Sungkawichai" in 07-reflections/) would have surfaced it; the S399 session — whose own §3/§8
preached MISTAKE-183's grep-the-statement rule — ran six subagent sweeps and still missed it,
because none of the sweeps targeted the lead's own statement.

**Why it is instructive (beyond MISTAKE-183 again).** (1) Synthesis sessions are not immune:
aggregating 175 mistakes into a taxonomy does not apply the taxonomy to the synthesis's own
outputs. Rule: before RECOMMENDING a lead, grep its statement exactly as if about to work it.
(2) The recovery matters: the S400 session proceeded anyway because the PROOF-level half
(what the paper's internal structure pins, vs what it states) was genuinely unmined — and that
half yielded the S-T tightness cage + height-258,276 rigidity (HYP-7920). A lead can be
half-dead and still be the right thing to work; the mistake is mis-stating WHICH half is open.
(3) Corollary for backlog hygiene: a lead entry should cite the canon items that ALREADY bear
on it (HYP-4096 should have been in the S399 backlog entry).

**Cost:** none in session time (the deeper mining was the correct move regardless); the cost is
a mis-calibrated recommendation broadcast to the fleet for one day.

## MISTAKE-190 (opus-2026-07-19-S402, against my own S401/HYP-7930) — propagated a mis-extracted inclusion chain into a HYP headline: the "1/13 floor" route rode a garbled star-index; the proven chain ends in gridmax-land, not lines-one-dimension-down

**What happened.** S401's HYP-7930 headline asserted "acc(13-speed M-spectrum) ⊆ [1/13, 1/2]
via their acc(S(n)) ⊆ S*(n−1) + settled LRC(≤13)". The pinning session (S402) shows the
proven chain in the PUBLISHED v4 is acc(S_k(n)) ⊆ S_{k+1}(n) ⊆ S*_{k+1}(n) = S*₀(n−k−1):
for lines the terminal set is S*₀(n−2) — D-values of FINITE proper subgroups = grid-loneliness
values of the (n−2)-torus — which has NO 1/13 floor (gridmax((1,…,11); 14) = 1/14 exactly,
verified; interior window values realizable). The "= S*(n−1)" k=1 line in my S401 source
fetch was a mis-rendering. Consequences of the error: corollaries (C2)/(C3) (window
finiteness; 12-speed gap = finite list) are CONDITIONAL on their Conjecture 1.5, not
derived; only (C1) floor-isolation survives unconditionally — and by a DIFFERENT (stronger)
route: the proven "only upper accumulation points" phrase (now THM-1289). Also the
abstract-vs-body contradiction (the /abs fetch served the v1/v2 abstract whose equality
claim was WITHDRAWN in v3 and demoted to Conjecture 1.5 in v4) went unnoticed for one day.

**Instructive because:** (1) two extractions of the SAME source can disagree; the tie-breaker
is internal cross-validation (here: Theorem 1.3's mirror statement validated the 1.4 phrase,
and the version history + footnote-2 error mechanism explained the discrepancy) — one fetch
is never a pin. (2) A DERIVED-status headline still propagates: agents read headlines, not
status fields; keep conditional clauses IN the headline. (3) The garbled star-index was
precisely the paper's own historical error (their footnote 2: subtorus∩subspace is a
subGROUP — disconnected); when an extraction smooths over exactly the subtlety that broke
the source's earlier version, suspect the extraction. Cost: one day of a mis-scoped
headline; the S402 pinning was already scheduled, and the designed process caught it.

## MISTAKE-192 (kind-pasteur-2026-07-19-S128c91, against my own HYP-7955 H-G2) — hunted the "extinction prime" c(p) > 13 while my own gate, two sessions earlier, had already proved the sea never empties: the AP inhabits I(13,p,1) at EVERY p coprime to 14

**What happened.** HYP-7955 seam G proposed the extinction principle "c(p) > 13 ⟹
I(13,p,1) = ∅" and HYP-8000 (ex-7975) hunted the crossing with compiled code (gated 9/9, greedy
triage over 130 primes to 1200, parallel exhaustive sweeps). The triage found every
prime ALIVE — because the answer was a one-line theorem I already possessed: the
acceptance-test gate (HYP-7930, two sessions earlier) verified "(1..13) improper mod
every tested p", and the general proof is immediate from the pinch/located-maximizer
theory (THM-401/1002): every maximizer of the AP {1..13} sits at a sum-14 pair, hence
at denominator 14, so at any t = a/p with gcd(p,14) = 1 the min-distance is STRICTLY
below 1/14 and the AP is improper. **c(p) ≤ 13 identically; the extinction prime does
not exist; H-G2 was vacuous as posed.** Same for k=12: {1..12}'s maximizers live at
denominator 13, so I(12,p,1) ∋ AP12 always — the "S-T retrodiction control" had a
predetermined answer too.

**Why it happened.** The extinction idea was born inside the covering-number frame
(seam G's c(p) growth data) and never re-checked against the improper-tuple frame's
own gates — the two frames are one object, and the gate FACT ("AP improper at every
tested p") was recorded as an instrument check, not read back as a structural theorem
quantified over all p. A vocabulary-internal blind spot: MISTAKE-183's grep rule
would not have caught it (same session family, same files); what catches it is asking
"does a known object inhabit the set I want to empty, for structural reasons?"

**Rule (the eternal-inhabitant check).** Before hunting for the emptiness/extinction
of a set, list the known structured objects and check each for membership AT THE
GENERAL PARAMETER — a gate verified "at every tested p" whose proof is p-uniform is a
theorem, not a check. (Kinship: THM-1289's only-upper-accumulation logic; the
tight families are eternal inhabitants of every level-1 sieve.)

**What survives (real yields, not salvage-spin).** (1) c(p) ≤ 13 for all p coprime to
14 — now PROVED, with the AP cover as the uniform witness; c(p) data pins the band
{≤12 → 12 → ≤13} with the 12-onset at p = 181. (2) The ragged greedy-frontier (15
greedy-hard primes among 130) survives as a SIEVE-COST predictor: greedy-difficulty
ranks how thin I(13,p,1) is per prime. (3) The instrument lesson: exhaustive-MRV
order is terrible at FINDING covers (minutes) where randomized greedy finds them in
0.1s — finding and refuting need different search orders. (4) The REFRAMED question
(the one the hunt was groping toward): not "does I(13,p,1) empty" (never) but "does
it COLLAPSE TO THE TIGHT-FAMILY SHELL as p grows" — i.e. does c(p) reach and stick at
13 with few, classifiable minimal covers — the large-p ansatz-completeness question,
now sharply posed with instruments in hand.

**Affects:** HYP-7955 (H-G2 struck), HYP-8000 (ex-7975) (resolution recorded in place), the
S128c90 reflection (seam-G consequence paragraph amended), CONSTANTS-INDEX (no
extinction-prime entry to add — the non-existence noted at the c(p) line).

## MISTAKE-195 — Claim-check and claim-write executed in one compound command (death-star-S59m)
Ran `grep -c HYP-8070` and the stub-write in the SAME Bash call; the grep revealed the collision (kind-pasteur had pushed 8070 minutes earlier) but the write had already executed before I read the output. Cost: one renumber cycle (8070 -> 8075) and an INDEX edit race. Rule sharpened: the fetch-grep-claim sequence must be SEQUENTIAL TOOL CALLS — inspect the check's output BEFORE the write ever runs; never batch a namespace check with the write that depends on it.

## MISTAKE-196 — Emptiness claims from un-validated search instruments (death-star-S59n, caught in-session twice)
(a) The general-k ansatz carried y^(k+1) where the weight arithmetic demands y^k; the row-solver then "proved" absurd rigidity (dim 0 spaces EXCLUDING the known witness). Caught because the known solution must lie in every claimed solution space — an impossibility check that costs one substitution. (b) The first numeric Newton hunt returned 0/800 hits INCLUDING at k=2 where the witness exists: the search was drowning in the degenerate attractor (C = delta = 0 gives det = 0 = "constant"), and the filter silently rejected everything. The k=2 validation gate exposed it; the fix (Keller-forcing equation c0(0)*nu = 1 deleting the degenerate branch) is the general pattern. RULE: a search that cannot rediscover the known witness has ZERO evidential weight for emptiness; run the witness-rediscovery gate BEFORE believing any negative result, and prefer surgically deleting degenerate branches (reciprocal-variable equations) over post-hoc filtering.

## MISTAKE-197 — Presenting a classically-trivial identity as a novel theorem (death-star-S59o/p)
THM-1320's headline factorization det = −E₀(0)A(0)C(0) is det JF(0) in disguise: any Keller map's constant det equals the det of its linear part, and the six-function form's linear part is the antidiagonal of the unit constants. The identity was derived through the c₀-formula and felt like discovered structure; one normalization (G = L⁻¹F, next session) exposed it as linear-part invertibility. Caught in-session-chain and amended. RULE: before filing a "factorization/necessity" theorem, evaluate it at the simplest distinguished point (the origin, the identity, the trivial case) and check whether the claim reduces to a textbook fact there; the FRAME can still be canon-worthy — file it as a frame, not a theorem.

## MISTAKE-198 — Three-instance recurrence fit claimed as a verified law (S419 diagonal sums)
- **Who/when:** opus-2026-07-19-S419, caught by opus-2026-07-20-S420 (same agent, next day).
- **What:** S419 stated the owner-triangle diagonal sums 1,2,4,7,12,21,37 "satisfy a(n)=a(n-1)+a(n-2)+a(n-4) (verified on all three available instances)". The statement was literally true and still a mirage: the law breaks at the very next term (m=8: predicts 65, actual 68 — both on the pure Faulhaber grid and on the owner's corrected triangle). The true sequence (1,2,4,7,12,21,37,68,129,...) has NO constant-coefficient recurrence up to order 8 (kind-pasteur-S128c103 lists it OEIS-new).
- **Genus:** category E (small-case extrapolation; Moser's circle / width-of-G_n). A 3-term-deep fit of an order-3-support recurrence has zero spare verification instances.
- **Rule:** a recurrence of order r fit on a sequence needs >= 2r+2 terms with the fit verified on the surplus; for families with super-exponential drivers (Faulhaber columns), ALWAYS compute at least one continuation term beyond the data used, before claiming any law.
- **Repair:** S420 reflection section 4 (museum of impersonations); banner added to the S419 reflection and HYP-8155.

## MISTAKE-199 (Nth recurrence, death-star-S61b/c/d) — three ID collisions in one session cluster on fleet-wide owner prompts
The owner's arborescence and odd/even/Babai-Cameron prompts went fleet-wide, and I collided THREE times by filing before re-pulling: (1) THM-1445 (switching H-sum) vs opus/kp THM-1445, pushed 11 min earlier — renumbered to THM-1460 then again; (2) THM-1460 (arborescence det-shadow) vs mac-mini THM-1460, pushed 23 SECONDS earlier and carried further (two poles, ordinal-sum log) — renumbered mine to THM-1467; (3) THM-1465 (canonical member / Babai-Cameron 7.4 = 0 at every odd n) vs kp THM-1465 (5 min earlier) AND opus THM-1460 (10 min earlier), BOTH identical (all-even anchor n≡1, all-ODD anchor n≡3, via klein's score-parity law) — CEDED entirely, my file deleted. Net: of my three S61b-d "theorems," two were independent rediscoveries of same-day fleet work and one (THM-1467 switching-sum) is the only distinct survivor, plus the 3/8-mass confirmation of boxeph HYP-8295. HARDENED RULE (again): on any owner prompt, `git fetch && rebase` IMMEDIATELY BEFORE the checkpoint that claims an ID, not just at session start — the fleet moves in minutes, not sessions. And when a prompt is visibly fleet-wide (Babai-Cameron had klein+opus+kp already), default to CONFIRMATION/synthesis and do not file a competing theorem number at all. The distinct-contribution test must be applied BEFORE filing, not after the collision.

## MISTAKE-204 (boxeph-2026-07-20-S182) — THM-1635's ladder was verified by injecting its own ansatz: a numeric that models β_eff(m) = β(1 + c/m) tests the ladder's ARITHMETIC, not its {m^{-k}} scale PREMISE

**What happened:** THM-1635 §3 built a 1/m-Vandermonde ladder on the premise
that tied jump coefficients vary on the scale set {m^{-k}}, and the S182
machine check "confirmed" it — by SYNTHESIZING data of exactly that form
(β(1 + c/m)) and watching the tied sum come out (c₁−c₂)/m. Circular: the
check could not have failed on the premise, only on the algebra.

**The truth (referee, machine-verified in
`04-computation/tie_ladder_scale_referee_boxeph_S182r.py` + frozen `.out`):**
the scale premise is FALSE for generic tie arcs. Inverting a generic branch
r(t) = C²/t² + A/t + ... produces HALF-STEP Puiseux corrections
t = C r^{-1/2}(1 + u r^{-1/2} + ...), and Laplace over the germ turns those
into β_eff(m) = β · e^{-u√(2m)}(1 + O(m^{-1/2})): the honest scale set is
{e^{a√m} · m^{-k/2}} (measured: log β_eff/√(2m) → −u, −0.2957 at m = 512 for
u = 0.3). Even INTEGER-step corrections v/r, which do give a true 1/m ladder,
dress the base constant to e^{-2v}β (0.5488 measured) — the ladder's "β" was
never the raw fold coefficient. Re a ≠ 0 breaks "slowly varying" entirely;
imaginary a puts √m-drifting phases under the Vandermonde's nose.

**Genus:** the numerical cousin of MISTAKE-203. That one: representation-level
ARGUMENTS dodge (only function-level invariants survive). This one:
ansatz-level NUMERICS dodge — a check whose test data is generated FROM the
model's own form verifies closure of the model under its own assumptions,
nothing more. Also kin to MISTAKE-186 (vacuous hypothesis: the check was
true, but of a vacuously self-referential statement).

**Rule:** before building (or "verifying") an asymptotic ladder, DERIVE the
scale set from the geometry (here: invert the branch, Laplace the germ, read
off which powers of m appear in the exponent AND the prefactor). A
verification numeric must be fed data from the GEOMETRY (integrate the actual
germ), never from the model's ansatz. If the only available test data is
synthetic, say "arithmetic check" in the writeup, never "machine-verified."

**Repair:** THM-1635 §6 (verdict + the √m-graded repair route), THM-1630 §6
(dressed-constant amendment; domination survives since e^{O(√m)} never
crosses (C_j/C_k)^m), HYP-8505.

## MISTAKE-203 (boxeph-2026-07-20-S178) — THM-1615's pinch dichotomy was false: fold-curve crossings are Stokes crossings of the representation, not singularities of the function

**What happened:** the pinch bridge (S177) inferred "the sweep of |t*(r)| forces
a pinch or endpoint singularity of A(t) at finite t." FALSE: as t moves, the
branch point r_j(t) crosses the REAL r-axis and the contour simply deforms
around it (a Stokes crossing of the representation). Genuine singularities need
a TRAPPED pinch (two branch points colliding — finitely many t) or endpoint
contact (excluded near finite t by the sweep's own r -> 0 blow-up). A generic
ray misses the finite trapped set — the theorem's own ray-genericity deleted
the singularities it needed. Also: the sweep itself fails when f_0(0) != 0
(t*(r) -> 1/f_0(0) finite), and the Gevrey/rotation citations from THM-1565 do
not transfer (integrand not entire in r).

**How caught:** owner-ordered adversarial pass (S178), referee 1, before any
completion claim was made (the theorem was filed as a skeleton with the attack
surface listed — the process worked).

**Lesson (two bridges dead by the same class of error):** both failed bridges
compared or tracked REPRESENTATION-level structure (term domination; fold
crossings) instead of FUNCTION-level invariants (actual singularities, actual
values). The Radial Lemma and the orbit-product proof survived review because
their contradictions are function-level (A -> 0 vs A == 1; (ct)^r vs constant).
Rule: a bridge argument must exhibit a function-level invariant that
distinguishes mixed from one-sided — anything contour- or term-level will
dodge.

## MISTAKE-202 (boxeph-2026-07-20-S176) — THM-1605's original local lemma overclaimed: equal products over disjoint subsets do NOT force equal subsets

**What happened:** THM-1605's first proof (S175) concluded from the Puiseux-DFT
expansion that disjoint index sets I, J with equal branch-products must be
EQUAL ("Fourier inversion => I = J"). FALSE as stated: the mu^j coefficient of
the log-difference is C_j * (Sigma_I zeta^{ij} - Sigma_J zeta^{ij}) with a
UNIVERSAL factor C_j = [mu^j] log(W(mu)/w*) that can VANISH; when it does, the
j-th character sum is unconstrained. Example: I = {0, 2} in Z_4 has S_1 = S_3
= 0 automatically, so with C_2 = 0 all constraints hold with I != J. The
step "no monodromy element moves the cluster" was therefore NOT established
(and is not even true in general — the identity only forces equal PRODUCTS
along the orbit).

**How it was caught:** in self-review while PREPARING the fleet adversarial
review (owner-ordered), before the referees ran. The repaired proof (orbit-
product: permanence + transitivity + Vieta + the c != 0 patch) is simpler and
was then confirmed by two hostile referees, who found one further pinhole
(c = 0 exclusion — patched in place) and two exposition debts (paid).

**Lesson:** when a lemma concludes SET equality from SYMMETRIC-FUNCTION
equality, check the kernel: products/power sums see subsets only through
characters, and universal prefactors (here C_j) can kill exactly the
characters you need. Also: the strongest review step is trying to write the
attack yourself before delegating it.

## MISTAKE-201 — Over-estimating the VC-witness reduction dimension 4× by a per-monomial heuristic, then reading noise as corroboration (death-star-S61)
- **Who/when:** death-star-2026-07-20-S61, caught same session by the reduction agent's concrete computation.
- **What:** I estimated the de Bondt–van den Essen reduction dimension of F at **≈76** using a crude "~⌈(deg−1)/2⌉ auxiliary variables per nonlinear monomial × 13 monomials" heuristic, then treated a concurrent agent's partial-output "76/77" hits as an independent cross-check. Both were wrong: (a) F's 13 nonlinear monomials are **not independent** — they all share the single quadratic u=1+xy, so **6** helper coordinates reduce everything to cubic (N≈10, not 35–38; M=2N≈20, not 70–77); (b) the "76/77" in the partial output was the ballpark of Zhao's *a-priori VC bound* (3/2)(3^{M−2}−1), NOT the reduction dimension — I pattern-matched a number to the conclusion I'd already written.
- **Genus:** category C (heuristic upper bound mistaken for an estimate) + confirmation bias (misreading unrelated partial output as agreement). The heuristic ignored shared sub-expressions; the "corroboration" was manufactured from noise before the agent finished.
- **Rule:** (1) an auxiliary-variable count for degree reduction must count **distinct building blocks after common-subexpression elimination**, never monomials independently — sparse/structured maps (here ℂ*-equivariant, all powers of one quadratic) collapse the count. (2) NEVER cite a running agent's *partial* output as corroboration; wait for its final result, and match numbers to their actual definitions (dimension vs a-priori bound) before claiming agreement.
- **Repair:** corrected §2 of the S61 reflection (≈76→≈20), the HYP-8265 sub-entry, the SESSION-LOG, and the PROBLEM-LEDGER. Correction broadcast to the fleet.
- **ADDENDUM (same session, on executing the transport):** the corrected "≈20, feasible, gate clears, no new math" was ITSELF overstated. Exact ℚ(i) computation (`polylib_exact_deathstar_S61.py`) showed: (a) Yagzhev G's Jacobian is NOT nilpotent, so the lift needs a genuine cubic-**homogeneous KELLER** reduction; (b) the "6-helper stacking" reduces degree and transports the collision but is Keller only ON the section {W=X^β} — off-section det varies (Schur complement = JG only there), so the reduced map is NOT globally Keller and NOT nilpotent-Jacobian; (c) naive homogenization breaks Keller too (det(I+x₀JH₂+JH₃) not constant). So the agent's stacking (which I relayed) conflated the easy degree-reduction with the hard Keller-preserving homogeneous reduction — the witness EXISTS (BCW) but its construction is real math (global-determinant control), not "engineering." **Deeper rule:** when relaying a subagent's feasibility claim, re-derive the load-bearing step yourself before writing "feasible/no new math" into canon — a dimension count is not a construction. The exact lift+rotation machinery WAS validated (produces Hessian-nilpotent quartics from nilpotent cubic-homog maps), so the honest status is "machinery ready, nilpotent input still to be built," not "done modulo typing."

## MISTAKE-199 — Re-deriving concurrent fleet work across two sessions on a fleet-wide owner prompt (death-star-S59s/S59t)
The owner's n·2^x+1 / figurate-triangle prompt went to the WHOLE fleet. I worked it in S59s (HYP-8165) and S59t (HYP-8175, THM-1355/1360) without first grepping for concurrent same-prompt HYP claims — and kind-pasteur (S128c102 Rosetta triangle, S128c103 shear catalog / HYP-8170), opus (S420 shear duality), and mac-mini (S1 Pisot tower) had gone DEEPER on the same objects: the Pisot ladder, the Proth spectrum, the products=ordered-tournaments identity, and the EXACT triangle (Faulhaber + 3 deviations) were all already theirs. My results were correct but ~90% duplicative; I even collided on HYP-8165 (both mine-S59s and kp-S128c102) and HYP-8170 (kp's; I took 8175). My repo-mining agent surfaced the overlap only AFTER I closed both sessions. RULE: on any prompt likely broadcast to the fleet (owner prompts usually are), the FIRST action is `grep -i <topic-keywords> 05-knowledge/hypotheses/INDEX.md` and a scan of same-day SESSION-LOG / inbox for concurrent claims — BEFORE claiming IDs or doing deep work. Convergent independent verification has value, but it must be FRAMED as convergent (credit + priority banner), never presented as primary. Amended THM-1360 + the reflection with priority banners; credited kind-pasteur/opus/mac-mini. **RECURRENCE (S59u, same day):** the NEXT owner prompt (the problem-ledger) was ALSO fleet-wide — klein-S332 wrote a peer ledger, kp-S128c104 + mac-mini-S140 stubbed. This time a mining agent surfaced the collision DURING the session and I consolidated into one canonical PROBLEM-LEDGER.md crediting all rather than competing — the discipline working. Rule: on any owner prompt assume fleet-wide, grep INDEX + same-day SESSION-LOG FIRST; on collision, CONSOLIDATE with credit, never produce a parallel artifact.

**MISTAKE-199 4th recurrence (S59v):** the weight-sign reduced JC was proven by mac-mini-S123 (THM-1370) and pushed WHILE I computed it, though the INDEX was clean when I claimed. Escalated rule: on a hot fleet-wide prompt re-pull IMMEDIATELY BEFORE filing any theorem, and default to synthesis/confirmation over primary claims once 2+ agents are visibly on the prompt. Credited mac-mini; kept only a verification cross-check + the three-reduced-JC synthesis.

## MISTAKE-200 — naming an invariant before grepping canon for it (d_sat = the diameter)

**What happened (mac-mini-S126, caught by kind-pasteur-S128c108 / THM-1400 SS I).** I
introduced the "saturation depth" `d_sat(n)` — the least `d` at which the waggly truncation
`G^(<=d)` becomes complete — called it a **new metagraph invariant**, computed `2,3,4,7` for
`n=4..7`, refuted `d_sat = n-2` at `n=7`, and handed off "compute n=8 before conjecturing."

**Why it was wrong.** `G^(<=d)` is complete exactly when every pair of classes is within `d`
flips, so `d_sat = diam(G_n)` — trivially, once stated. And the diameter was already canon:
`07-reflections/diameter-is-feedback-arc-set.md` (opus-2026-03-24-S306) records
`diam(G_n) = max_T min-FAS(T) = A003141(n)`, growth `~n^2/4`, with `OPEN-QUESTIONS.md`
listing it **RESOLVED** and `README.md` carrying it as the Waggly Completeness Theorem.
So: the invariant was a rediscovery, the handoff needed no computation (`A003141(8) = 8`),
"no linear formula" was a known **quadratic**, and my `n=7` refutation restated opus-S306's.

**The rule.** *Before naming an invariant, grep the reflections, OPEN-QUESTIONS, and README
for it — including under other names.* A quantity defined by "the least `d` at which X
becomes complete" is a **diameter**; ask what metric, then search for that metric. I searched
for "saturation" and "d_sat" (my own coinage, so of course no hits) and not for "diameter."

**What survived.** The map-graph *framing* of THM-1390 (point-contact vs edge-contact, the
clique-explosion reading) is untouched and THM-1400 explicitly grants it; the point-adjacency
thread continues correctly in THM-1405. The cost was one invented name and one bad handoff,
both retracted in place. See `01-canon/theorems/THM-1390-...md` (correction banner).

## MISTAKE-202 — Trusting an asymptotic domination claim measured only to m=8 (death-star-S61g/S61h)

**What happened.** For GMC(2) I declared "GMC(2) is complete" (S61g) on top of klein-S351's
Gamma-Bridge domination step, and separately claimed (S61h) that for two-sided `P` the top
charge-0 term `γ_{a_max} a_max!` carries `> 50%` of the mass of `E[P^m] = Σ_a γ_a a!`, so
`|top| > |rest|` and the triangle inequality forces `E[P^m] ≠ 0`. I measured the share only
across **m = 2..8**, where it sits at ≈ 0.60–0.67, and read that as domination.

**Why it was wrong.** kind-pasteur-S128c120 (THM-1585) refuted klein's version exactly — the
top-term share falls to **0.04%**, the consecutive-term ratio grows to **45×** — and opened
`CASE-gamma-bridge-domination-step`. Prompted by their addendum I reran **my own** statistic
to **m = 24**: the share collapses `0.67 → 0.068` (`Z²+W+ZW²`) and `→ 0.0002` on a b-sweep;
top-dominance is False for every `m` past ≈ 8. **A domination claim measured to m = 8 cannot
tell "share → 1" from "share → 0" — both look like ⅔ there.** `E[P^m] ≠ 0` still holds, so
the *conclusion* (NC2) survives, but the *mechanism* was false: "domination was an analytic
strategy for an algebraic fact" (kp). The sound route is kp's THM-1605 Hermite no-common-root.

**The rule.** An asymptotic-domination / "top term wins" claim is only tested by pushing the
statistic until the ratio's *trend* is unmistakable (here m ≈ 20+), and by a parameter sweep
that can *amplify* the competing term (the b-sweep). Never certify "X dominates as m → ∞" from
a window where the sub-dominant term is still the same order. When the honest tool is a
factorial-weighted sum `Σ c_k k!`, suspect that nonvanishing is *algebraic* (orthogonality /
no-common-root of a classical sequence), not a size comparison. Conceded in full in the court
case; S61g headline withdrawn, S61h §1 retracted (both banners in place); S61h §2 (the Lean
NC2 ⇒ GMC(2) reduction, pure charge arithmetic, assumes NC2) is unaffected. See also
MISTAKE-199 (the same over-eager pattern, there in fleet-coordination form).


## MISTAKE-207 (boxeph-2026-07-20-S186r) — THM-1765's pair-sum constant used the symmetric-collision model (missing the Λ''' midpoint shift), and its evidence verified the RATIO, never the CONSEQUENCE object

**What happened:** THM-1765 §1 claimed the far-end pair-sum limit
(1/t)·2/(u_c²Λ'') and hence the universal constant −2/(d₁d₂). Two errors:
(a) the merging roots are NOT symmetric about u_c at second order — the
Λ''' midpoint-shift enters at the SAME order; the correct law is
PS = (1/t)[2/(u_c²Λ'') + (2/3)Λ'''/(u_cΛ''²)], giving −(2/3)(d₁+d₂)/(d₁d₂)
on a two-term edge (referee-measured: (1,−1) → 0 — forced exactly by the
global residue identity Σ_{all roots} 1/(uΛ') = 0; (2,−1) → 1/3; (3,−1) →
4/9; claimed: 2, 1, 2/3). (b) The frozen evidence (F1–F4) measured only
the RATIO Λ''u²/v — the identity's hypothesis — and never once measured
the pair-sum, the object the Consequence is about. The identity was true;
the Consequence drawn from it was wrong; the check could not have caught
it because it tested the wrong object.

**Also refuted downstream (same session, same verdict):** "O(1)
universally" — the charge-0 coefficient p₀ enters v but never Λ'/Λ'':
VALUE-HIJACKED ends (witness P₄ = ZW + Z⁹W⁷ + W) realize the S183r threat
with non-integrable s^{−2} pair-sums; and the (L2) two-mechanism
classification + O(t^{−1/2}) rate (zero-drift mechanism; Θ(T^{−2/5})
example).

**Genus:** the leading-term trap, third instance (MISTAKE-202: leading
products; MISTAKE-206: leading odd-coefficient; 207: leading collision
model), COMPOUNDED with MISTAKE-204's cousin: the numeric verified an
intermediate quantity, not the concluded one.

**Rule:** (1) at a colliding pair, never use the symmetric local model
without checking the next Puiseux coefficient (Λ''' here) — or better,
compute pair-sums by the GLOBAL RESIDUE IDENTITY (sum over all roots is
zero: pair-sum = −spectator residues), which is exact and model-free.
(2) Every verification numeric must measure the CONSEQUENCE'S object, not
the hypothesis' — if the theorem says "hence X is bounded", the script
must print X.

**Repair:** THM-1765 amended in place (corrected law; value-hijack
sub-case named; ledger reverted; residue identity canonized; §6 verdict
archive); THM-1680 §4 pointer corrected; referee checks frozen at
04-computation/thm1765_referee_{hijack,momenttest}_S186r.py + .outs.

## MISTAKE-206 (boxeph-2026-07-20-S183r) — THM-1680's deletion biconditional read only the LEADING odd coefficient: "B ≡ 0 ⟺ removable" is false one rung down, and the dichotomy missed the boundary-truncated class

**What happened:** THM-1680 §1/§2 filed "DEFECT ≡ 0 ⟺ B ≡ 0 ⟺ removable" and
"per germ exactly two cases (delete / graded-visible)". The hostile referee
(S183r) refuted both AS STATED: (a) the local odd sector at a fold is a full
Puiseux vector Σ c_k τ^{k/2}; B = c₋₁ is only its leading entry. A germ with
B ≡ 0 but c₁ ≢ 0 is NOT removable and stays Γ-graded exactly one rung down
(T3: C(1/2,m)/C(−1/2,m) = −1/(2m−1) exact; T4: its monodromy defect 2|c₁|√r
sits twelve orders above the instrument floor — THE INSTRUMENT WAS FINER
THAN THE PROSE). (b) The ε-signs can hop at sub-germ boundaries, so
"vanishing on an interval ⟹ vanishing on the germ" propagates only within a
sub-germ; a sign flip mid-germ creates BOUNDARY-TRUNCATED arcs with pure
exponential moment grade I_m ~ e^{−s₂}v(s₂)^m/m (T2) — a third class.

**Genus:** the leading-term trap, one floor up from MISTAKE-204. 204: a
check fed by the model's own ansatz verifies arithmetic, not premises. 206:
a FUNCTION-LEVEL invariant (good, per MISTAKE-203) that is only the LEADING
coefficient of the honest invariant (the full odd-sector vector) inherits
exactly the failure it was built to avoid, one rung down. Cousin of
MISTAKE-202 (leading products vs full set data).

**Rule:** when the invariant is a Puiseux/graded expansion, the deletion or
vanishing criterion must quantify over the WHOLE graded vector, and every
"identically zero" must name the domain it propagates over (sub-germ, germ,
family). If an instrument (here the monodromy defect) is sensitive to the
full vector, let the instrument's zero-set DEFINE the criterion rather than
the leading coefficient's.

**Repair:** THM-1680 §1/§2/§4 amended in place (odd-sector vector; per-
sub-germ trichotomy with the truncated exponential class; §4 hypothesis =
defect ≡ 0 everywhere; finiteness lemma added); referee verdict archived as
§8; checks frozen at 04-computation/thm1680_referee_hostile_S183r.py + .out.
Referee: "no scenario survives the repaired statements."

## MISTAKE-205 — the "Alpoge-Mathew" attribution of the JC counterexample (THM-1300) is a HALLUCINATION (owner-corrected 2026-07-20)

**What happened.** THM-1300's attribution blocks (mac-mini-S127 and S129) confidently attributed the
dim-3 Jacobian-Conjecture counterexample to **Levent Alpoge** with co-credit to **Akhil Mathew**,
"obtained with Claude Fable 5, announced 2026-07-19," citing web searches (X posts, a Wikipedia edit,
an arXiv absence check). The OWNER has flagged this as a **hallucination**: those searches were run on a
result then roughly ONE DAY OLD, and the specific-name attribution is not reliable — it should not be
restated as established fact.

**The corrected provenance (per the owner).** The counterexample was **discovered by Claude**. It became
public on **2026-07-19 via a tweet from an Anthropic employee**, which is where the "shared publicly
yesterday" signal came from — that is a SHARING event, not a discovery attribution. So the map is a
Claude discovery that an Anthropic employee surfaced; it is NOT an Alpoge-Mathew result.

**The lesson (this is the reusable part).** Do NOT attribute a very recent result from web searches with
confidence. A result that is ~1 day old has essentially NO reliable public record: no peer review, no
stable arXiv, and search engines index tweets and rumors that conflate "who shared it" with "who found
it." When mac-mini's own search notes flagged "no arXiv preprint, no peer review, no journal," that was
the tell that NO reliable attribution was available — and the correct move was to record the map as
**provenance-uncertain**, not to fill the vacuum with plausible names. Confabulated attribution is exactly
the failure mode a fresh, un-refereed result invites.

**Standing correction for the fleet.** Treat the JC counterexample as **Claude-discovered, provenance
otherwise uncertain**. Do NOT write "Alpoge-Mathew" (or "Alpoge") as the discoverers in any file, letter,
or write-up. THM-1300's external-attribution blocks are kept as history but are now marked contested by
the banner added this session. What the repo legitimately holds is unchanged: the INDEPENDENT exact
verification (det J ≡ -2, triple collision, all identities) and the equivariant/Dixmier/elliptic anatomy.

---

## MISTAKE-208 — a "poly-time invariant" that is not isomorphism-invariant (arborescences rooted at a fixed vertex)

**kind-pasteur-2026-07-21-S128c143.** In the invariant-lattice work I used `arb(T)` = number of
spanning arborescences **rooted at vertex 0** (one reduced-Laplacian minor). This is NOT an
isomorphism invariant: the root "vertex 0" depends on the labeling, so relabeling changes it. On the
exhaustive census it *looked* fine because I computed it on ONE canonical representative per class
(deterministic, but an artifact of the canonicalization). The bug surfaced when a random-labeling
SAMPLER reported spurious "same-invariant, different-H" collisions at n=6 that the exact exhaustive
census (poly-tuple determines iso at n=6) said were impossible.

**Why it's instructive:** (1) a per-root / per-vertex quantity is only an invariant after
symmetrization (sum over roots, or the *sorted tuple* of per-root values). The fix is
`arb_inv = tuple(sorted(arb_root(T,r) for all r))`. (2) Two code paths (canonical-rep vs
random-label) disagreeing is a RELIABLE detector of a non-invariant — always cross-check an
"invariant" on random relabelings. (3) The scare (a suspected floating-point artifact) was itself a
false lead — the numbers matched numpy exactly; the real bug was mathematical (root-dependence), not
numerical. Corrected: `|arb_inv|=55` at n=6 (vs 32 rooted-at-0); `arb_inv` refines score. Headline
findings (score ⟂ cyc, the 2-adic edge THM-1980) use only genuine invariants and are unaffected.
Banner added to THM-1965.
## MISTAKE-216 -- confusing Poisson rank two with a DC(2) or planar-JC counterexample

**Error.** From a nonautomorphic Poisson endomorphism of
`C[x,q,p,z]`, infer either a counterexample to the planar Jacobian conjecture or
an endomorphism counterexample of the second Weyl algebra merely by replacing
Poisson brackets with commutators.

**Correction.** The four-variable symplectic map is a Keller counterexample in
affine dimension four, not two. Quantizing its four output polynomials directly
in `A_2` requires exact control of ordering corrections; equality of principal
Poisson symbols does not give exact Weyl commutators. The safe cotangent or
Hamiltonian-dual construction doubles four classical coordinates to four Weyl
pairs and therefore lands in conventional `A_4`. To reach `DC(2)`, verify an
actual endomorphism of `A_2`; to reach planar JC, construct or exclude a
two-coordinate Keller pair. Keep Poisson rank, affine dimension, Weyl rank, and
number of algebra generators separate.

## MISTAKE-282 (2026-07-27, mac-mini-S144, self-report on mac-mini-S131's THM-1315 §3-§4) — the surjectivity claim was FALSE and the on-caustic fiber count was wrong (2 is impossible): a MULTIPLICITY error in the shared-root case analysis

- **What was claimed (THM-1315 §4, S131):** F is surjective, with fiber count 3 off K and 2 on K; "φ shares at most ONE root with the denominator conic ⟹ at least two honest preimages" for a ≠ 0, c ≠ 0.
- **Why it is wrong:** the shared root is the DOUBLE root of the fiber cubic (that is exactly the Δ ∝ Res ∝ a²K coherence I had myself proved in §3), so one shared VALUE removes TWO root-copies: the fold pair escapes jointly and the count on K is 1, never 2. Deeper, on the rational curve E = {K = 0, 3bc = 4} ≅ ℂ* the remaining root also degenerates and the fiber is EMPTY: **F(ℂ³) = ℂ³ ∖ E — F is not surjective.** Both facts are in opus's THM-2473 (the depressed x-eliminant Lx³ + (4−3bc)x − 2c with x²-coefficient identically zero), whose companion cube disc_a(L) = −4(3bc−4)³ has E's wall as its root. Independent verification this session: the fiber over (4/27, 4/3, 1) is empty (exact Gröbner); (2/27, 1, 1) has exactly one preimage (2, 5/6, −7/8).
- **The compounding process error:** klein-S324's master-quartic note ("fiber spectrum {3,1} never 2") was ALREADY ON DISK when I wrote §4 — a grep for the statement would have caught the clash before canonization. This is the MISTAKE-183 pattern (derive-instead-of-grep) recurring on my own output; and the mechanical lesson is new: **count fiber drops in root MULTIPLICITY, not in shared root values** — a shared value at a double root costs two.
- **What survives:** THM-1315 §1 (syzygies, fiber cubic), §2 (S₃ pin — generic targets), §3's caustic identity and ramification-at-infinity; the corrected Euler ledger now balances (1 = 3·0 + 1·(1−0) + 0 with χ(K) = 1, χ(E) = 0). Banner added to THM-1315; the S142 reflection updated.
- **Source:** opus-2026-07-27 repair letter (MSG-1605) + THM-2473; verified independently in-session.

## MISTAKE-318 (2026-07-28, root audit of klein-S691) — a changed-diagram unknot certificate proves `u<=1`, not `u=1`

**What was claimed.**  The first version of
`04-computation/unknot1_decider.py` returned `TRUE_CERTIFIED` whenever one
crossing change produced a diagram which greedily reduced to the empty
diagram by Reidemeister R1/R2 moves.  The implementation attempted to reject
the input unknot first, but used the same incomplete greedy R1/R2 reducer.

**Why it is wrong.**  The crossing-change certificate proves only

```text
u(K)<=1.
```

To conclude `u(K)=1`, one must separately prove `K` is nontrivial.  Failure of
one simplification heuristic is not such a proof.  The exact hostile PD

```text
[[1,11,2,10],[6,10,7,9],[3,8,4,9],
 [11,5,12,4],[7,2,8,3],[5,1,6,12]]
```

is an unknot obtained from the one-crossing unknot by legal reverse moves
`R1,R2,R2,R3,R3`.  Greedy R1/R2 stalls on this six-crossing input.  Changing
crossing four then reduces by `R2,R1,R2,R1`, so the old engine returned
`TRUE_CERTIFIED` although the represented knot has `u=0`.

**Repair.**  A true certificate now has two independent halves:

1. an input-nontriviality certificate (currently `det(K)!=1` or
   `sigma(K)!=0`); and
2. a one-crossing change followed by an explicit R1/R2 unknot certificate.

If (2) holds without (1), the engine returns `UNKNOWN` and reports only the
proved upper bound `u<=1`.  The hostile PD is frozen in the 16-check
ordinary/optimized regression suite.  The Murasugi and gated Lickorish
`u>=2` false-certificates are unaffected.

**Reusable rule.**  An exact witness for the result of a move does not certify
the starting object's distance is exactly one.  Every exact-distance
certificate needs both an upper-bound path and an independent lower-bound
obstruction; a failed search is never the latter.

## MISTAKE-319 (2026-07-28, root audit of klein-S691) — the `3+3+1=7` fragment is Keller escape, not an Arithmetic-Kakeya `mu_3` identification

**What was claimed.**  A truncated owner fragment about a depth-two tree with
seven rather than nine points was provisionally decoded as a `mu_3`-fixed
branch on which two points coincide.  The AK workbench then treated that
putative collision as evidence for identification-gluing.

**Exact resolution.**  The fragment is THM-2473's inverse tree for the
sporadic Keller map.  Its three level-one points have `3,1,3` finite
preimages.  On the middle branch the cubic eliminant loses its leading
coefficient and becomes linear; two sheets escape on the Jelonek
nonproperness surface.  Because the Jacobian determinant is the nonzero
constant `-2`, no finite collision or ramification occurs.  The symmetry
fixing the middle point is the order-two involution
`diag(-1,-1,1)` and the generic monodromy is `S3`, not `mu_3`.

**Repair and survivor.**  The fragment note and AK workbench now withdraw the
transfer.  Identification-gluing remains a legitimate independently
testable AK mechanism: a merged vertex may inherit several incident species.
It receives no support from this Keller fibre.  A transfer based only on the
shared count `9 -> 7` loses the decisive coordinate—finite vertex
identification versus loss of properness at infinity.

**Reusable rule.**  Reconstruct a fragment's generating object before using
its numerical pattern elsewhere.  Equal cardinality defects do not identify
mechanisms; record whether mass collides, cancels, is quotiented, or escapes
the ambient space.
