# Euler mechanisms integrated with the research portfolio

**Status: PROVED scoped results and FINITE-EXACT controls; major-problem
closures remain unclaimed.** Session date: 2026-09-08. The mathematical base
was `687644309fd3`; work proceeded in an isolated worktree and was integrated
through small mainline checkpoints. The external PDE claim is not a dependency
of either new theorem.

## What this session established

1. [THM-4457: sharp transverse-shear distance budget](../../01-canon/theorems/THM-4457-euler-sharp-transverse-shear-distance-budget.md)
   proves, for every trace-free real matrix `M` with nonzero vorticity axis
   `xi`,

   `|xi^T M xi| <= sqrt(7/6) dist_F(M, rank-one trace-free shears)`.

   The constant is sharp. Along a classical Euler trajectory, amplification
   by `A>1` requires an integrated departure from pure shear of at least
   `sqrt(6/7) log A`. Pure shear cannot stretch its own vorticity, even when
   its frame moves. The proof includes an exact compact quartic optimization
   for the distance and classical continuation consequences.

2. [THM-4458: one-sided adverse-leak budget](../../01-canon/theorems/THM-4458-lrc-one-sided-adverse-leak-budget.md)
   improves the [THM-3682 tariff](../../01-canon/theorems/THM-3682-lrc-pure-role-to-lawful-target-leak-tariff.md).
   For real profiles `C=R+L`, centered `r`, and `E=mean(r^2)>0`, put

   `B=min_b mean(-r_s(L_s-b))_+`.

   A weighted median computes `B` exactly, and
   `Var(C)>=(B-E)^2/E`. This is sharp at every budget. It dominates the old
   reverse-triangle bound. An exact mass-preserving indicator-density
   example has `Var(L)>Var(R)` but `B/Var(R)=1625/2508<1`: the old gate is
   silent while the new gate certifies nonconstancy.

3. The [source audit and Poisson refinement](../reference/OPENAI-EULER-AUDIT-2026-09-08.md)
   identifies the explicit Lean arctangent profile with a Poisson kernel.
   Its sharp derivative lower bound is `-1/(2+delta)`, improving the
   source's used bound `-1` for that particular profile. This sharpens a
   favorable-sign pressure constant; it does not validate the entire PDE
   construction or apply to every manuscript-admissible profile.

The Navier--Stokes shear-distance continuation condition is already implied
by [Miller's middle-strain-eigenvalue criterion](https://arxiv.org/abs/1710.05569).
THM-4457 explicitly proves this inclusion. The pointwise material estimate
is a different statement, with its own proof and sharpness example. No
worldwide priority claim is made for the elementary results here.

## Inheritance and concept board

The anchor was a usable transfer from the supplied PDE work to a live
major-problem obligation. The niche was the geometry of shear matrices and
their material vorticity axis. The wildcard was an optimally shifted
positive-part cost. The wildcard produced the LRC improvement.

Closest mechanisms: THM-3682's corrected profile tariff and the classical
vorticity equation. Canonical hostiles: complete shifted-danger cancellation,
AP13, and a pure shear carrying a passive vector. Corrected near misses:
present-Q is not delayed-Q(Rx), and an ambient shear product is not an Euler
trajectory. Least-used sidecars: signed role/leak alignment and the difference
between the actual vorticity axis and the shear's own axis.

| Live concept | Operation and result | Information still needed |
|---|---|---|
| One-sided adverse cost | Replace a full norm by a convex, constant-invariant hinge budget; THM-4458 succeeds | A structural bound for the same canonical observer and actual leak |
| Rank-one shear | Minimize distance over orientation and amplitude; THM-4457 gives a sharp growth cost | Full gradient on the same material trajectory |
| Independent phase versus physical restriction | Keep the evaluation map before averaging | Resonances, owner, word, and source/target action |
| Harmonic smoothing | Identify the explicit profile derivative with half a Poisson defect | Nonnegative flux and the chosen profile |
| Repeated amplification | Separate passive, stagewise, and same-solution growth | A supplier that carries the target through every stage |

Every positive result was compared back against these five concepts. The
adverse budget keeps the sign discarded by a norm. The sharpness family shows
that vorticity zeros matter. The source's profile choice makes a constant
improvement available but supplies no physical LRC map. The final tests keep
the original target fixed.

## Genuine maps and their limits

**PDE sign budget to finite profiles.** The transferred operation charges
only the adverse contribution to a fixed scalar response. It becomes the
functional `L -> B_R(L)` above. Sublinearity and a weighted-median dual prove
the variance consequence. The PDE, its pressure, and its moving phase do not
map to LRC by this construction. The auxiliary density example establishes
that the new inequality is useful without establishing canonical entry.

**Shear decomposition to material amplification.** The map takes the full
gradient to its nearest pure shear while retaining its own curl axis. A
matrix norm or passive-vector propagator would lose that axis. For
`M=S+E`, the identity `omega^T S omega=omega(E)^T S omega(E)` exposes the
mechanism: the shear annihilates its own vorticity on both sides. The actual
material logarithmic growth, not a proxy norm, is the consequence object.

**Shear maps to LRC.** If integer `a,m,s` satisfy `m.a=m.s=0`, then
`I+k a m^T` fixes every point of the orbit `ts`, despite arbitrarily large
operator norm. Thus even exact unimodular amplification can leave all safe
times unchanged. Covariantly transported inequality normals preserve the
target; deleting them changes it. This is a proved stopping obstruction,
not an Euler-to-LRC identification.

**Whole-stage preservation.** The repo's
[THM-2022 balanced-face mechanism](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
and the external packet construction both motivate retaining the complete
object needed by the next operation. Their algebraic and PDE predicates are
different. No Frobenius-to-fluid, Jacobian-to-flow, or tournament equivalence
is asserted. There is no intrinsic tournament needed for these inequalities.

**Incoming Jacobian result, commit `4329ebb75`.** The newly integrated
[(7,1) formal-inverse theorem](planar_jc48_sep08_seven_one_formal.md) gives
an entire inverse hierarchy with a formal Jacobian mate and a literal
example having no rational mate. This is a stronger hostile than failure
at some bounded jet order. It reinforces the fifth board item: preservation
at every formal stage does not supply the required target object. Here the
missing condition is rational algebraization; in a PDE expansion it is
analytic realization. This is a comparison of proved boundaries, not a map
between those predicates or a claimed defect in the external construction.

## Exact verification and explicit stopping points

- Shear audit: 275,562 matrix/shear presentations from 6,561 trace-free
  matrices, plus rational rotations, a sharp limiting family, and pressure
  sign controls. [Code](../../04-computation/euler_shear_distance_budget_20260908.py),
  [output](euler_shear_distance_budget_20260908.out).
- Leak audit: all 7,020 nonconstant-observer pairs over `{-1,0,1}` in
  lengths 2, 3, and 4. Hinge minimization, weighted median, and a separate
  dual-polytope enumeration agree. The 13-phase control is reconstructed
  from rational circle intervals. [Code](../../04-computation/lrc14_euler_adverse_leak_20260908.py),
  [output](lrc14_euler_adverse_leak_20260908.out).
- Normal and optimized Python runs agree with checks enabled. The finite
  universes supplement written universal proofs; neither computation is a
  PDE solver or an LRC counterexample search.
- The external audit pins a source commit and the actual 57-page PDF hash.
  A successful dependency-cache fetch is not a completed Euler build. Its
  exact build boundary and static import counts are recorded separately.

The strongest unresolved LRC step is now concrete: for the same present-Q
observer, establish `B_R(L)<=b0<121 rho^2/2028` by source/owner structure.
That would yield the already typed variance and drift tariffs. Supplying
the delayed-Q(Rx) intertwiner remains a separate obligation. Changing the
observer to lower its budget would repeat MISTAKE-547.

For Euler, compute or bound the integrated shear distance along a **single
fixed solution** before using THM-4457. Comparing two different stage
solutions is not its hypothesis. In the packet decomposition the background
belongs in `E`; it cannot be discarded because the child shear is large.
The criterion does not currently contradict or independently establish the
external singularity claim. Obtaining an a priori bound on this budget is
OPEN, and existing strain criteria already subsume the stated viscous
continuation specialization.

META-PATTERNS used: search the statement before the method; correct the
object before sharpening the technique; type every analogy and implication.
No new general method card was promoted from a single imported analogy.
