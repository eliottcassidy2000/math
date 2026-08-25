# `p-adic-zeta-irrationality` source audit (2026-08-25)

**Status: AUTHOR-CLAIMED / UNREFEREED RESEARCH DRAFT; SPECIALIST AUDIT OPEN.
Numerical margin certificate VERIFIED-EXACT relative to the manuscript
formulas. Do not cite the 22-value theorem as proved canon.**

## Frozen source and provenance

- Christopher D. Long,
  [*Hybrid arithmetic holonomy and twenty-two individual p-adic zeta values*](https://github.com/octonion/p-adic-zeta-irrationality),
  frozen at [commit `b46a1770901551961710e155d775aae7c5ea39e7`](https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7),
  commit time 2026-08-24 21:05:53 -0400.
- The repository is MIT-licensed and contains TeX/PDF plus one standard-library
  Python interval certificate. At the frozen commit it has no tagged release,
  CI proof harness, or executable geometric worked example.
- Frozen SHA-256:
  - TeX: `f4d8e0582d3ce3ca88d0553b9ae5a212007fc01da169188d22ed3e0311974fd3`;
  - PDF: `1e10c2eb9d9f30ec7f357ffc0ab276236f807e91d15c5f32d8a1d5a3cf6ac75f`;
  - certificate: `1408c9092d0c34917253c8f520db56853112c58640d15c98d58b76a61a0478f5`;
  - stored output: `9e1d8f74eb198c7e7b0c7420a1e8021d7bad4d2bc8014cc939b6b353111c0da0`.

## Claimed theorem and mechanism

The manuscript claims irrationality of 22 Kubota--Leopoldt values:

```text
zeta_2(s),  s=3,5,...,29       (14 values),
zeta_3(s),  s=3,5,...,11       (5 values),
zeta_5(3), zeta_5(5), zeta_7(3).                          (1)
```

Its proposed mechanism is hybrid arithmetic holonomy for the rank-`s+1`
Eisenstein--Eichler connection on the genus-zero curves `X_0(p)`,
`p in {2,3,5,7}`:

1. at small auxiliary primes, Hasse algebraization plus a genuine-source
   rank-nullity kernel saves one determinant power on many directions;
2. at large auxiliary primes, a de Rham frame-torsor calculation, divided
   Frobenius, and coefficient-ring Taylor no-backflow recover a prime window;
3. Bost/CDT slopes and modular Jensen energy assemble the local costs.

The displayed conditional cost is

```text
tau_(p,s)(xi)=2[S_(p,s)(xi)+s I_s(xi)]/(s+1)^2,           (2)
```

and rationality would force

```text
M_(p,s)(xi,Y)
 =s Lambda_p(Y)-(s+1)tau_(p,s)(xi)-C_p(Y) <= 0.           (3)
```

The code certifies positive values of `(3)` at 22 rational pairs `(xi,Y)`.

## What was reproduced

In a clean temporary clone:

```bash
sha256sum -c hybrid_rational_interval_certificate_sha256.txt
python3 -B hybrid_rational_interval_certificate.py
python3 -B -O hybrid_rational_interval_certificate.py
```

The manifest passed; both Python modes byte-matched the frozen output. A
separate 110-digit direct reconstruction put every one of the 22 values inside
the printed intervals. The smallest was

```text
M_(5,5)(29/27,1/16)
 =0.131799356827016832557664457131... >0.               (4)
```

A direct eta-product quadrature for the `(5,1/16)` Jensen term converged toward
the displayed value. The certificate's `Fraction`, integer fixed-point,
Machin/atanh, `isqrt`, and outward-rounding arithmetic appears internally
sound for this numerical task.

This verifies **only `(3)--(4)` conditional on the paper's formulas and
theorems**. The program hardcodes the 22 cases. It neither constructs the
geometric source nor checks any new local theorem.

## Published-input spot checks

The following cited interfaces exist and broadly match their stated roles:

- [Calegari, *Irrationality of certain p-adic periods for small p*](https://arxiv.org/abs/math/0408214):
  the overconvergent Eisenstein family and constant term;
- Buzzard's finite-slope continuation theorem: arbitrary integer weight is
  allowed;
- Graham--Pilloni--Rodrigues Jacinto, Propositions A.14--A.15
  ([DOI](https://doi.org/10.1112/S0010437X25102479)): neat-level `B`-torsor,
  filtration, torus weights, and Gauss--Manin structure;
- [Calegari--Dimitrov--Tang 2024](https://arxiv.org/abs/2408.15403): scalar
  zero estimates, pivot heights, and prime-window formulas;
- [Calegari--Dimitrov--Tang 2025](https://arxiv.org/abs/2510.04156): adelic
  assembly used by the manuscript.

These checks validate citation existence and broad normalization, not the
new internal bridges.

## Highest-risk proof obligations

No short contradiction was found, but the following gates need specialist
line-by-line review.

1. **Fixed-bundle comparison.** The internal `lem:fixed-bundle` transfers
   scalar CDT source and pivot estimates to `H^0(P^1,V(DP+B))`, including
   integral finite-place behavior, by Birkhoff--Grothendieck splitting and
   bounded frame changes. The leading archimedean `O(D)` comparison is
   plausible; the integral pivot/window transfer is not a direct cited
   consequence.
2. **Descent and one-action bridge.** Uniform descent/spreading, one flat
   integral torsor chart, `f`-support `[-C,D+C]`, and the claimed single
   degree-`s-1` unipotent action are compressed into short internal arguments.
3. **Small-prime Hasse kernel.** Rank-nullity is correct once the exact coarse
   degree, fixed pole clearer, base change, saturation, and identification of
   the genuine deepest multiplier all hold uniformly in `ell,D`.
4. **Large-prime Cartier strictness.** The coefficient-ring no-backflow
   algebra is coherent on paper, but no finite worked case constructs the
   actual torsor algebra, full-product digits, scalar pivot, and pole grade.
5. **Global assembly.** The functional/BGG lift, modular Jensen
   multiplicities, and admissible adelic template passage remain internal
   proof obligations.

The cheapest decisive test is one explicit `(p,s)=(5,5)` case naming the
frame, torsor chart, divisor `B`, exceptional set, normalizer, source bases,
one small-good-prime Hasse matrix with its kernel, and one large-prime
pole-grade pivot with `n=j+ell r`. Re-running `(4)` is not that test.

## New exact formula frontier

[THM-4089](../../01-canon/theorems/THM-4089-hybrid-padic-zeta-margin-optimization-and-next-case-obstruction.md)
proves, without importing the arithmetic-holonomy theorem, that the displayed
two-parameter margin has an explicit global cutoff minimizer and a unique
analytic-radius maximizer. It independently certifies the 22 formula signs
and proves negative global maxima for the four immediate next cases:

```text
(2,31), (3,13), (5,7), (7,5).                            (5)
```

Thus extending `(1)` one odd weight in any row requires a stronger cost or
energy theorem, not finer search for `xi,Y`.

The same theorem also identifies a distinct source-level wall. In the first
active arithmetic chamber the stationary point satisfies

```text
xi*-1=(s-1)(11-p)/(12s^2+(s-1)(p+1)),                   (6)
```

so it is interior exactly for `p<11`. A complete five-chamber calculation at
`p=13,s=3` then proves that even the idealized bound `J_13(xi)<=xi` has
`tau>=613/288`; since `log(13)<8/3` and the remaining analytic terms are
nonpositive, every continuation radius satisfies

```text
M_(13,3)<-37/72.                                        (7)
```

For the actual formal integrand, `J_13(xi)=3/7` on `xi>1`,
`inf tau=515/224`, and the bound sharpens to `M_(13,3)<-67/56`. This is a
no-go for the displayed one-power architecture, not evidence that
`zeta_13(3)` is rational or that another method cannot reach it.

## Lawful repo interfaces and firewalls

| source object | repo target/map | preserved | destroyed / needed sidecar |
|---|---|---|---|
| `v_ell(lcm(1,...,N))` | THM-4056 denominator clock / THM-390 prime-power obligation | valuation depth and cutoff | LRC carrier, phase, owner, safe time |
| binomial atom mod `ell^a` | Sun role-labelled residue tree | congruence and local multiplicity | top-index height, exact lift equality, cancellation |
| raw depth `n^e` / `L_N^e` | Apéry denominator condition | per-prime denominator bound | explicit recurrence, approximant, nonvanishing, decay roots |
| pivots versus prime windows | weighted bipartite incidence graph | which place can charge which pivot | magnitudes if reduced to orientation |

- THM-4027 proves universal modular solubility for Sun's 2-4-6-8 form, so no
  fixed p-adic congruence obstruction explains the counterexample. A lawful
  p-adic search must retain bounded lift heights.
- `L_5=60` is a genuine common LCM object with valuation vector `(2,1,1)`,
  but it does not identify the Fibonacci Pisano period, Sturmian AP phase, or
  their owner clocks.
- Tournament Analysis is not native here. Prime/pivot incidence has ties and
  hyperedges; a forced total orientation would discard the additive adelic
  budget which makes `(3)` meaningful.

The 22-value theorem remains **AUTHOR-CLAIMED / UNREFEREED** until the
geometric and adelic gates are independently verified.
