# Planar Jacobian extended session: audited results and next frontier

**Status: OPEN research program; individual results carry separate statuses.**
Started 2026-09-06, based on `2d3c53942d51`, with the user-requested addition
of Yinjie Li, *The S-matrix conjecture*, arXiv:2608.29750v1.
The planar Jacobian conjecture remains open.

The assigned anchor is planar JC. The inherited mechanisms are
THM-4438 / `jc2-row-fifteen-relative-response-on-boundary-gm`,
THM-4364 / `sharp-binomial-diagonal-annihilator-hierarchy`, and the
[ninth-source collision period](overnight8_20260906_jc_collision_period.md).
Canonical hostiles are omitted weight-eighteen memory (THM-4426), the
collision-breaking ninth jet `t^9*x`, and a complex isotropic vector `(1,i)`
against an unqualified sum-of-squares argument. The corrected near miss is
to call a frozen-prefix obstruction a source obstruction. The least-used
sidecars are complete diagonal polynomial traces, simultaneous moving
collision sections, and the entire raw response cokernel.

| Lane / concept | Object and question | Operation / invariant | Lost information and cheapest test |
|---|---|---|---|
| Anchor: earlier source memory | THM-4438 partial source; can lower-weight changes replace the weight-24 payment? | Re-solve earlier rows and retain joint bracket/depth response | A same-row quotient forgets earlier source coefficients; test complete omitted faces |
| Anchor instrument: full depth image | All depth modules `P_d` in source-normal coordinates | Split by `ell=2n-r`; test divisibility of signed diagonal polynomials by `(1-z)^ceil(ell/3)` | Individual annihilators need not be complete; compare exact full matrix ranks and source lifts |
| Niche: moving collision period | Exceptional Russell compiler and compensated branch sections | Differentiate the full collision identity; retain common motion and tangent relation | One chosen form hides incompatible branch motion; test every two-form slot and earlier clocks |
| Wildcard: global defect/projection | Li's real matrix defect identity and JC response equations | Replace scalar diagnostics with full measurement maps and dual certificates | Positivity, coefficient bounds, and integrality do not transfer automatically; test `(1,i)` and arbitrary rescaling |
| Termination versus projection | Finite polynomials versus compatible finite jets | Seek an effective finite recognition bound once support is fixed | Unbounded support allows new payments forever; construct late-row hostiles |

## Connection contract for the new paper

Source: the exact defect budget and Gram projection in Sections 2, 4, and 6
of [Li v1](https://arxiv.org/html/2608.29750v1). Target: finite JC bracket,
depth, and collision measurement modules. The proposed map retains the
complete linear measurement vector before eliminating repair variables.
It preserves exact compatibility, provided the actual repair image is used.
It discards the real norm and binary rounding structure; these need a
separate positive form, bounded region, and integral separation if used.
The cheapest algebraic test is equality of the proposed annihilator kernel
with the actual source image, including every raw compatibility.

This is an inspiration and a research contract, not an S-matrix-to-JC
implication. No external Lean build or whole-paper proof audit is claimed.

## Audited first checkpoint

1. [Sharp finite recognition](planar_jc_long_20260906_depth.md) is
   **PROVED + INDEPENDENTLY AUDITED**. For `T>=1`, a polynomial of t-degree at most T,
   with the inherited depth-d caps, belongs to P_d iff it passes the
   projected test through `floor(4T/3)+d+1`. The cutoff is optimal, with
   actual-source hostile `t^T`. The all-height surface version is
   `floor(h^2*T/(h+1))+d+1` under its explicitly defined weighted caps.
   Producer45,048 and independent828,684 exact gates pass. The principal
   diagonal image was recovered from THM-4369/3973; the sharp recognition
   consumer is the new result. Support must be fixed, including zero tails.
2. [Complete source response](planar_jc_long_20260906_memory.md) is
   **PROVED FINITE-ROW + INDEPENDENTLY AUDITED**. In valuation>=13 and
   weight<=23 the full eight-channel response has three conditions. The
   minimum replacement weight is22, and the source family at that weight
   is `G_m x A^2` with terminal `A^10`. The actual transport is credited
   to incoming `d0208173b2 / planar_jc48_sep06_weight14.md`; its full
   filtered classification and minimality are new. Producer191 and
   independent248 exact gates pass, using different row/depth coordinates.
3. [Moving collision periods](planar_jc_long_20260906_collision.md) is
   **PROVED + INDEPENDENTLY AUDITED**. The exact identity
   `Pi+chi' Xi=0` holds at every order. It excludes constant fixed-source
   densities throughout the inherited three-parameter compensated family
   and every nonzero formal clock. An explicit formal unit density and
   source change show the precise failure boundary of coordinate-free
   overstatements. The independent audit has133 exact gates.
4. [Algebraic dual testers and Student torsion](planar_jc_long_20260906_smatrix.md)
   is **PROVED + INDEPENDENTLY AUDITED**. Generic dual residual coefficients
   generate exactly the one-right-hand-column mixed maximal-minor ideal.
   The actual unrestricted Student operator has an all-index Smith formula
   from two weighted paths. Its row14 integral response exponent720720
   is attained by an explicit seven-column bank. This is an arithmetic
   sidecar, not a characteristic-zero Keller obstruction. Producer2,455
   and independent13,868 exact gates pass.

The [source reading](../reference/CORE-PAPERS-SMATRIX-2026-09-06.md) records
the paper's hypotheses and the exact limit of the transfer. No external
S-matrix theorem is used as a mathematical dependency of these results.
The [session manifest](planar_jc_long_20260906_manifest.json) pins the frozen
artifacts and their reproduction commands. The original scouts have now
been proved and independently audited in the second round below.

## Audited second round

5. [Earlier source memory](planar_jc_long_20260906_memory_earlier.md) is
   **PROVED FINITE-ROW + INDEPENDENTLY AUDITED**. Allowing coordinate
   changes from row ten and source changes of valuation at least eleven
   gives the complete weight-22 source family `G_m x A^9`, with separate
   `A^10` coordinate fibres. Seven exact equations classify the response.
   The background genuinely varies, but a polynomial inverse with a
   square-zero correction proves a constant response quotient at every
   specialization. Weight 21 still fails. The odd `py^5` channel is neutral
   through row thirteen and forced to vanish at fourteen, refuting a parity
   shortcut. The intermediate [valuation-twelve family](planar_jc_long_20260906_memory_valuation12.md)
   is `G_m x A^5`; valuation thirteen gave `G_m x A^2`. These are nested,
   fully classified finite horizons, not an infinite continuation theorem.
   The valuation-eleven producer/background have329/276 gates and the
   independent smaller-coordinate path5,868; valuation twelve has1,646
   producer and5,146 independent gates. All normal/optimized pairs agree.
6. [Completed Hamiltonian repairs](planar_jc_long_20260906_hamiltonian.md)
   are **PROVED + INDEPENDENTLY AUDITED**. Every inherited universal
   carrier `S=c+p^2(p^3-y^2)R` lowers positive actual depth and satisfies
   `delta(D^N P_d) subset D^(N+1)P_d`. The scalar-time exponential is an
   inverse automorphism of the full actual `D`-adic completion, preserving
   the completed source family and volume. Literal source and Laurent
   comparisons have separate proofs; their complete rings are not equated.
   The [maximal depth carrier](planar_jc_long_20260906_depth_carrier.md)
   is exactly `c+p^2R`: preserving depth zero already preserves every depth.
   Its intersection with preservation of even one affine source is the
   smaller universal carrier. Every nonconstant depth carrier is non-LND.
   The completed carrier has528 exact controls; the maximal carrier has
   1,589 producer and56 independent generic-formula controls.
7. [Every nonzero time is non-rational](planar_jc_long_20260906_nonrational.md)
   is **PROVED + CITED CURVE INPUTS + INDEPENDENTLY AUDITED** for the
   explicit infinite family `S=f(p^a(p^3-y^2))`, `a>=2`, nonconstant `f`.
   These operations retain source depth and source form at all orders,
   but cannot specialize to rational source maps at any nonzero time.
   Their generic invariant curves have genus `ceil(a/2)+1`; finite curve
   automorphism groups contradict the flow's nonzero displacement at every
   positive iterate. A separate injective source-series comparison closes
   the topology obligation. The non-LND rational flow generated by `x^2t`
   is the hostile to a weaker argument. Producer226 and independent2,268
   exact controls pass; the universal conclusion rests on the analytic
   proof and its precise curve inputs, not a finite rational-map search.

## Incoming connections and remaining obligations

Fetched incoming `d0208173b2` and read its current board, relevant proofs,
and correction additions before this checkpoint. Its weight22 transport
is the same compensated line, and its Hamiltonian carrier/non-nilpotence
result prevents repeating a failed additive-flow search. Its collision
quadrics describe a different, three-dimensional spanning-tangent case;
our moving-clock theorem uses only a necessary relation in a tangent plane.
The new cusp passports and explicit higher-cusp target remain separate
geometric questions and are not imported as a consequence of our jets.

- **Earlier memory:** valuation ten would start coordinate corrections at
  row nine and expose the seventh actual background row. The recovered
  affine coordinate `xi10(s)` already parametrizes the boundary, but the
  next full response is uncomputed. No lower-weight universal claim is made.
- **Formal repair versus the finite response:** the completed carrier
  preserves depth and source form together. It remains to determine whether
  a needed finite response lies in its image; the inherited five-term finite
  Hamiltonian does not automatically lie in that carrier. Even a positive
  response test would leave polynomial specialization and termination.
- **A sharper termination test:** the explicit high-genus families are now
  excluded at every nonzero rational time. Search for a useful carrier with
  genus-zero or genus-one generic invariant, or a composition whose
  termination can be proved. No composition exclusion was established.
- **Fixed support:** the sharp recognition theorem decides depth once a
  finite polynomial candidate and its zero tail are given. An a priori
  support bound remains absent. Collision-clock exclusion retains its fixed
  source volume and specified compensated family; chart entry remains open.

Incoming `ebc62e7d7d` repaired a false all-order ballot-ratio claim and a
confusion between recurrence roots and polynomial root parameters. Its
specific arithmetic is orthogonal here. The shared methodological control
is concrete: finite True flags cannot certify an untested prose implication.
The all-index recognition, Smith, and flow claims above each have a proof
separate from their finite control universes.

Every resulting object was compared with fixed support, full raw
compatibility, moving common data, integral saturation, and the actual
formal-versus-polynomial boundary on the concept board above. The source
response and flow lanes now meet at a precise remaining image test, rather
than an assumed implication from finite solvability to termination.
