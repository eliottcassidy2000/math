# Plane-coloring blueprint: exact extensions and a finite spectral hostile

**Session date:** 2026-10-20 UTC.  
**Input:** owner-supplied *Blueprint: monochromatic triangles in two-colorings of the plane*, 670 lines, SHA-256 `e71c29ce6fba832425140a604d66189d21754ebc2d4525a9b873529125553ea1`. Entire input read. No Lean build was requested or performed. Its large global theorem, generalized carrier rigidity and formal-check claims remain **SOURCE CLAIMS**.

## Results

1. **A stronger circle obstruction, PROVED.** A measurable circle sign function is antiperiodic under all shifts `α_j` exactly when one integer `n` satisfies `nα_j≡π mod2π` for every j. In particular, simultaneous alternation under `π/3` and `γ` permits precisely reduced odd/odd ratios `γ/π=p/q`. An independently proved density-boundary argument transports incompatibility to measurable weak avoidance of equal-leg isosceles plane templates. [Full circle proof](circle/circle_antiperiodicity.md).
2. **Exact quantitative circle optimization, PROVED.** Rational angle constraints reduce to a finite weighted maximum-cut problem. With `π/3` alternation imposed, the minimum γ-defect is zero in the odd/odd case and `gcd(3,q)/q` otherwise. Irrational incompatibility has defect infimum zero, so qualitative exclusion supplies no uniform positive circle gap. [Arithmetic checks](circle/validation.md).
3. **Loss-one robust binary transfer, PROVED abstractly.** Fine independent cell thresholds preserve `Q/4` in expectation for separated triples; some measurable coloring achieves any prescribed weighted defect objective no worse than the bounded field. The coefficient one is asymptotically sharp for the stated general placement-measure formulation. Application to the blueprint's invariant carrier process is **CONDITIONAL** on its nonzero pair-correlation/projection input. [Transfer proof](proofs/binary_transfer_and_degree_boundary.md).
4. **Exact arity boundary, PROVED.** The monochromatic k-point Boolean function has degree `k` for even k and `k-1` for odd k. No nonnegative quadratic certificate can vanish at every nonmonochromatic k-tuple for `k≥4`. Four-bit parity gives an exact hostile with identical marginals through order three and different monochromatic probabilities. [Degree proof](proofs/binary_transfer_and_degree_boundary.md#4-boolean-degree-where-pair-transfer-stops).
5. **A finite-field threshold, FINITE-EXACT + elementary proof.** The equilateral norm hypergraphs over `F_(5²)` and `F_(11²)` admit the supplied explicit two-colorings, while `F_(17²)` does not. The 289-point obstruction has a complete 3,451-node refutation independently verified from its norm graph. The eight-point gadget then proves that every distinct three-point template is forced in every two-coloring of `F_(17²)`. [Finite theorem](proofs/finite_norm_triangle_bridge.md#3-exact-equilateral-threshold-at-q17).
6. **A precise spectral hostile, PROVED using the finite certificate.** At the same q17, a rotationally invariant positive-definite correlation has `ρ(1)=−1/3`, exactly the local avoidance value. It violates the explicitly proved cut inequality `Cut≤1732` by two and cannot come from an actual Boolean field. Thus radial PSD plus local triangle consistency loses global Boolean realizability. [Construction](proofs/finite_norm_triangle_bridge.md#4-the-radial-positive-semidefinite-hostile-at-the-same-q).
7. **A small all-shapes spectral bound, certified rationally.** At q23 every distinct template has at least 16 monochromatic ordered placements (at least six unordered equilateral triples). All 22 Kloosterman sums are enclosed by exact rational arithmetic, with maximum strictly below 8. This is a finite norm-plane bound, not an archimedean-plane result. [Certificate](results/kloosterman_bounds.json).

## Inheritance and concept board

| Lane | Inherited mechanism | Hostile / missing coordinate | New operation and outcome |
|---|---|---|---|
| Anchor: circle boundary | Blueprint §7 null boundary, apex forcing, self-support | Rationality alone omits parity; irrationality supplies no uniform gap | Common Fourier support; exact odd/odd criterion, finite max-cut, sharp constrained defects |
| Niche: HN finite norms | Corrected THM-4628 Kloosterman multiplier; THM-4663 polygon/walk ladder | Reduction creates incidences; PSD is not a Boolean law | Change observable to Q/4; q17 checked refutation and explicit PSD false positive |
| Wildcard: robust rounding | Blueprint §6 bounded projected sign | Common threshold synchronizes unrelated values; degree four retains parity | Fine cell addresses restore independent rounding and preserve weighted defect at factor one |
| Algebraic carriers | Blueprint §§3,8 one-parameter algebraic closures; THM-4700 unit-path ring | Higher transcendence degree and prime-ideal data are not automatic | Keep carrier assumption explicit; transport only the algebraic ζ6 gadget unconditionally |
| Finite automata / proof carrier | NAE constraint propagation and circle shift action | A manufactured tournament loses the symmetric ternary predicate | Native finite cyclic Cayley graph or checked binary proof tree, with actual maps |

The canonical hostile in the real plane remains the alternating strip coloring, which prevents lifting the finite q17 equilateral obstruction as a faithful real unit-equilateral hypergraph. The corrected prior scopes were read before using the residue and polygon results. No Hadwiger--Nelson bound, LRC statement, or global source theorem is promoted by this packet.

## Complete source-to-extension map

| Blueprint part | Equations or mechanism read | Use in this session / boundary |
|---|---|---|
| §1 reduction / eight points | `ζ+barζ=1`, six congruent copies, equilateral forcing | Six formal identities and all256 labelings checked independently; gives q17 all-template theorem and `D_T≥D_eq/6`. Source global assembly not imported. |
| §2 two scales | x=s/r, far/near split | Preserved as source architecture; circle test is strengthened but no claim replaces all carrier/cyclotomic inputs. |
| §3 carriers | Q-span template selection, wild-character rigidity | Explicit conditional input to robust transfer; higher-transcendence extension not assumed. |
| §4 label / radial spectrum | `1+ab+ac+bc=4 mono`; `ρ(r)=−1/3`; radial Fourier measure | Exact B1/F1; shared-frequency Kloosterman analogue F2; q17 PSD hostile blocks the reverse implication. |
| §5 Bessel bounds | angular J0 multiplier, pair sum, energy and Taylor bounds | The operation maps to norm-circle Fourier averaging; rational cosine certificates at17/23 are independently checked. The source's `−83/125` Bessel bound is not revalidated. |
| §6 binary transfer | clamp, projected coefficient identity, Q≥0, threshold/Tonelli | Loss-one cell rounding; general three-point conditional transfer; exact degree-four parity obstruction. |
| §7 measurable exclusion | alternation, apex, density boundary, self-support, rational angles | Independent proofs of all measure steps used; common-frequency classification and finite rational defect optimization. No quantitative plane-density bound inferred. |
| §8 transcendental template | Cayley escape, Cauchy frames, branch rationality, annihilator quotient | Fully read, retained as source-claim prerequisite; no false higher-transcendence generalization. The genuinely transported algebra is instead the finite ζ6 forcing gadget. |
| §9 comparison with original | alternative proof implementations | Status/provenance retained. No Lean theorem names are cited as independently checked certificates. |

## Reproduction

Run from this packet's directory:

    python3 scripts/check_binary_and_controls.py
    python3 scripts/check_eight_point_forcing.py
    python3 circle/check_finite_cases.py
    python3 scripts/verify_q17.py
    python3 scripts/certify_kloosterman_bounds.py

To regenerate the finite search and full refutation:

    python3 scripts/probe_finite_triangles.py 5 11 17
    python3 scripts/certify_q17.py
    python3 scripts/verify_q17.py

All scripts use the Python standard library. Main runs were on Python3. The search is not an unstructured UNSAT assertion: the emitted tree cites each forced assignment and each monochromatic conflict. The verifier uses a different geometric construction, norm adjacency rather than rotated templates. Finite checks do not certify the analytic proofs or any source Lean build.

## Remaining work

- Reduce the support or simplify the proof of the explicit q17 separating inequality `Cut≤1732`.
- Minimize the q17 obstructing hypergraph, and determine its exact monochromatic minimum.
- Classify avoiding non-equilateral shapes at q11; its equilateral control is known.
- Extend robust rounding to exact attainment without a common separation scale; an arbitrary ε-loss is already proved.
- Supply a genuine fourth-moment preservation condition for four-point transfer.

The independently proved circle exclusion and finite norm-plane theorem are scoped completed results. The higher-order, quantitative-plane and general-carrier targets remain open.
