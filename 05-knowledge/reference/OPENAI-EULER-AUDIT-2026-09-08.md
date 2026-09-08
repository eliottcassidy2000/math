# Bounded external audit, 2026-09-08

## Scope and immutable inputs

- Git source: `8937a8f4cbc7abaab5e9e97d1cc7f5d2319d9538`, cloned from <https://github.com/openai/NavierStokesAndEuler>.
- Directly downloaded Euler PDF: 57 pages; SHA256 `a0c234518e6c489e16996805023eb2e75c00b7c03455f7a3a5be2c124954bfdd`.
- The browsing tool served a different, cached 45-page extraction at the identical PDF URL. All page references below use the 57-page download. The abstract/main conclusion agrees, but section and equation numbers changed; the URL alone is insufficient provenance.
- This is a bounded source-level audit, not independent certification of the full PDE result. No demonstrable mathematical error or top-level vacuity was established in the inspected paths.

## Main theorem and reach

The external README's Navier--Stokes statements permit forcing. Its
whole-space conclusion also specifies a uniform kinetic-energy class.
Those statements are not imported as theorems here. THM-4457's viscous
continuation corollary explicitly treats unforced finite-energy solutions;
forcing would add a curl-of-force term to its enstrophy calculation.

`Euler/Solution.lean:43`, `Euler.exists_compact_smooth_euler_singularity`, has no construction premises: it existentially supplies an ordinary nonzero compact smooth datum, positive lifespan at most one, ordinary velocity and scalar pressure, Sobolev Euler on the lifespan, maximality in its all-order Sobolev class, local norm finiteness, and endpoint norm divergence. The separate global nonexistence clause is in the broader smooth finite-energy class.

The definitions in `Euler/SolutionDefinitions.lean` use actual Frechet derivatives, within-time derivatives, pointwise curl and L2 spatial derivatives. The fallback `toL2=0` outside MemLp is protected by explicit MemLp assumptions wherever used. Thus it is not an evident vacuity mechanism.

Static import traversal finds (all 2486 tracked Lean files were checked: 6868 import lines, zero multi-module/public/private/meta import syntax exceptions):

| Root | Tracked local modules including root | Comparator challenge imports | Lexical sorry/admit/axiom declarations |
|---|---:|---:|---:|
| Euler | 1771 | 0 | 0/0/0 |
| Euler.Solution | 1829 | 0 | 0/0/0 |
| NavierStokes | 581 | 0 | 0/0/0 |
| NavierStokes.ComparatorSolution | 580 | 0 | 0/0/0 |

`Euler` imports the internal singularity result but does not reach `Euler.Solution`. The lake target's `Euler.+` glob includes that wrapper, and Comparator explicitly targets `Euler.Solution`. The root-import distinction is therefore not evidence of an unbuilt theorem under the default build. The Comparator references intentionally contain four `sorry` placeholders across Euler and Navier--Stokes; neither proof dependency graph imports those modules.

The repository supplies `#print axioms` in `Euler/Solution.lean` and permits only `propext`, `Classical.choice`, and `Quot.sound` in `ComparatorChallenges/Euler.json`. These are author-supplied checks until an independent build prints their results; lexical scanning is not a substitute.

## Mechanisms checked

1. `Euler/CylinderDirichletData.lean:30` specifies a positive time, a pointwise lower frame bound, actual first/second frame derivative laws, the Jacobi law, an upper Hessian quadratic bound, and `potential*T^2/2 <= 1/2`. No inverse or solution is inserted into this input structure. `:94` constructs the coordinate solver through `Euler/TransverseFixedSpaceInverse.lean:74` and Lax--Milgram in `Euler/EulerProof.lean:599`. Zero Hessian and constant identity frame give an elementary satisfiable test instance.
2. `Euler/PacketInfiniteConstruction.lean:68` chooses scales through `exists_scales` and `:72` recursively constructs the actual stages from a first stage and successor maps. The final theorem does not merely assume an infinite successful family; fully auditing the long supplier graph remains outstanding.
3. `Euler/PacketPhysicalLowBounds.lean:97` preserves the sign of the flux `c=<m,Mv>`, and `:139` obtains a positive-Hessian cost reduced by the profile parameter delta. The corresponding manuscript mechanism is Section 5.7, pages 52--53. Its early/history intervals retain absolute costs; their separate decay estimates are essential.

## Independent exact refinement: Poisson profile and sharp sign cost

The source's specific periodic profile is

`f_delta(theta) = atan(sin(theta)/(1+delta-cos(theta)))`, delta > 0.

Its derivative is

`g_delta(c) = ((1+delta)c-1)/((1+delta)^2-2(1+delta)c+1)`, c=cos(theta).

Set r=1/(1+delta). Direct algebra gives

`g_delta(c) = (P_r(theta)-1)/2`,

where `P_r(theta)=(1-r^2)/(1-2r cos(theta)+r^2)` is the Poisson kernel. Hence the absolutely convergent Fourier series is `f_delta'(theta)=sum_{n>=1} r^n cos(n theta)`.

The derivative of g_delta with respect to c is

`(1+delta)*delta*(delta+2)/((1+delta)^2-2(1+delta)c+1)^2 > 0`.

Thus its exact range is `[-1/(delta+2), 1/delta]`, with equality at theta=pi and theta=0 modulo 2pi respectively. The source only uses the weaker lower bound `-1` at `Euler/EulerProof.lean:11825`.

For a>=0, c>=0, m nonzero, let `Q=-2*a*c*f_delta'(theta)*(m tensor m)/|m|^2`. Its sharp uniform quadratic upper bound is

`Q <= (2*a*c/(delta+2))*I`.

Equality for a*c>0 occurs at theta=pi in direction m. At theta=0, Q is negative semidefinite of norm `2*a*c/delta`. This supplies an explicit elementary refinement of the source's `2*a*c` positive cost; it changes a constant, not the PDE proof status or an open-problem verdict.

The sign sidecar is necessary. Choose a=1, m=e1, v=e2 and M(v)=-e1, so c=-1. For delta=1/100 at theta=0 the positive eigenvalue is 200; a proposed sign-free bound 2*a*abs(c)=2 fails.

Scope: the PDF Section 3.1 uses an arbitrary bump-primitive profile meeting its abstract bounds. The Poisson identity and sharp constant apply to the explicit Lean instantiation, not every admissible profile.

Reproduce static traversal and 20 exact-rational controls with:

```powershell
python3 04-computation/openai_euler_static_audit_20260908.py --repo EXTERNAL_CHECKOUT --pdf DOWNLOADED_PDF
```

## Build record

`lake build Euler.Solution` was started in this external clone. It installed the required Lean 4.34.0-rc2 toolchain and fetched pinned mathlib `85e3a25e006c35636f0e53b0e9296caca2685bc0` and Comparator `19e111e2141cf333c7daff0f64c5f24acc91dd2e`. The command began compiling uncached Mathlib; it was voluntarily interrupted with Ctrl-C after approximately 500 preliminary tasks, before compiling an Euler theorem. This was a bounded-scope stop, not a discovered source/build failure.

The recommended prerequisite `lake exe cache get` then completed successfully with exit 0: 8747 of 8747 cache files downloaded and decompressed. The full 1829-module Euler dependency build was not restarted within this session. There is **no independent Euler theorem axiom output and no Comparator/nanoda result** from this audit. Downloaded toolchain, dependencies, and cache remain available for a resumed audit.

The separate small scalar certificate [standalone scalar root](../../04-computation/lean/standalone/EulerSharpProfileAudit.lean) imports only Mathlib, not the Euler project. The final command

```powershell
lake env lean -DwarningAsError=true ABSOLUTE_MATH_REPO/04-computation/lean/standalone/EulerSharpProfileAudit.lean
```

completed with exit 0. It proves the exact rational-expression lower/upper bounds, both endpoint equalities, and the improved scalar pressure cost. Its printed axiom dependencies for `slope_lower`, `slope_upper`, and `sharp_scalar_pressure_cost` are exactly `[propext, Classical.choice, Quot.sound]`. The output is retained in [frozen scalar build output](../results/euler_sharp_profile_lean_20260908.out). This independently validates the elementary refinement; it does not validate the external PDE construction or formalize its arctangent derivative identification.

All owned build/cache/check sessions have terminated; no background audit work remains running.

## Reproduction scope in this math repository

[Static checker](../../04-computation/openai_euler_static_audit_20260908.py)
and [frozen output](../results/openai_euler_static_audit_20260908.out) are
preserved together. The checker counts literal import/admission syntax in
Git-tracked source only; it is not a semantic scanner or a kernel audit.
All 6868 import lines in the 2486 tracked Lean files use the supported
single-module `import` syntax. The code was replayed normally and under
`-O` with checks enabled; outputs match. Fetch the exact external commit
before replaying; a PDF with a different hash is a different source version.

The scalar Lean file itself is the standalone root passed to Lean: all its
lemmas are reachable there. It is deliberately not advertised as imported
by the repo's separate Lean-4.30 TournamentH7 root. The scalar hypotheses
are satisfiable, for example `delta=1`, `c=0`, `a=flux=1`. The certificate
proves a rational-expression inequality only; the trigonometric derivative
identification and Poisson equality have separate elementary proofs above.
The complete Euler formalization has not passed this session's validity gate.

Primary links: [pinned external source](https://github.com/openai/NavierStokesAndEuler/tree/8937a8f4cbc7abaab5e9e97d1cc7f5d2319d9538),
[manuscript URL](https://cdn.openai.com/pdf/315b36cd-ec98-4023-8342-93345194ece1/euler.pdf).

## Classical continuation inputs and prior-art boundary

- **CITED:** Beale--Kato--Majda, *CMP* 94 (1984), 61--66,
  [DOI](https://doi.org/10.1007/BF01212349), gives Euler continuation under
  finite integrated maximum vorticity in the classical Sobolev setting.
- **CITED:** [Tao, APDE 6 (2013), Corollary 5.8](https://msp.org/apde/2013/6-1/apde-v6-n1-p02-s.pdf),
  supplies whole-space Navier--Stokes H1 continuation.
- **CITED:** [Miller, arXiv:1710.05569v4, Theorem 1.1](https://arxiv.org/abs/1710.05569v4),
  gives scale-critical control by the positive middle strain eigenvalue.
  [THM-4457](../../01-canon/theorems/THM-4457-euler-sharp-transverse-shear-distance-budget.md)
  proves its NS sufficient condition is already implied by Miller. Its
  sharp pointwise trajectory estimate is separately proved; no general
  regularity theorem or literature-priority claim is inferred.
