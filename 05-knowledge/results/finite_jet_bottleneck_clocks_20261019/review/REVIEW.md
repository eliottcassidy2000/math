# Independent review: finite-jet bottleneck clocks

19 October 2026. **PASS for the stated local, fixed-endpoint, constant-additive family.** No unresolved mathematical defect was found. This review makes no novelty claim. Exact finite algebra and numerical diagnostics remain distinct from the analytic proof and the cited divergence theorem.

## Reviewed sources and replay

The reviewed mathematical note is `FINITE_JET_BOTTLENECK_CLOCKS.md`, SHA-256:

    383adcb99f6bf2c516afee94600502d972df3c075d96f39d8b36ef976d9fe318

The reviewed author verifier is `verify_finite_jet_clocks.py`, SHA-256:

    528bfb29876f3fbdbb7c29e5a2f400838bdf023ad6572cf53bc2abb28ec98527

I read the full note and verifier, independently checked the analytic estimates, inspected the cited sources, and replayed a copied author verifier outside the repository checkout. Normal and optimized Python both passed all 157 explicit checks. Their stdout and JSON matched byte-for-byte, including the author's captured JSON. The independent checker passed 107 exact checks in both modes, with byte-identical stdout. Versions and reviewed hashes are recorded in `reviewed-source-sha256.json`.

## Analytic audit

1. **Quantifiers and stopping.** The theorem chooses one sufficiently small neighborhood for each fixed analytic germ, then permits each fixed negative/positive endpoint pair in it. Constants may depend on the germ and endpoints. Positivity makes the orbit monotone, the lower increment bound epsilon proves finite exit, and the stated derivative/epsilon bounds keep the full final step inside the analytic domain. No arbitrary entry point, moving endpoint, epsilon-dependent germ, or uniformity over a jet class is smuggled in.

2. **Corrected clock.** For `q(t)=f(y+t f(y))/f(y)`, the identity `q'(t)=f'(y+t f(y))` is correct. Taylor and quadratic derivative bounds give the stated remainder estimates. The `-A/2` integral defect and `+A/2` logarithmic defect cancel. The coefficients `5/12` and `7/12` in the error constant are correct. Summing `D f(y_n)` telescopes to a bounded displacement. The fixed nonzero endpoints bound the endpoint logarithm; the final-step comparison bounds the overshoot integral by `1/(1-eta)`. Hence `N=integral+O(1)` is genuinely uniform in epsilon.

3. **Complete expansion and absence of a logarithm.** Exact flattening `u=y(g(y)/y^r)^(1/r)` is valid locally with the positive real-analytic branch. The transformed interval contains a fixed symmetric inner interval. Odd monomials cancel there; the two outer pieces are uniformly bounded. The remainder `u^r R(u)` is integrably bounded, and each even degree through `r-2` contributes the stated convergent beta integral. The borderline `u^(r-1)` term is odd, with the exact bounded endpoint log ratio shown in the note. Asymmetric endpoints and odd original Taylor coefficients therefore cause no omitted divergent logarithm.

4. **Coefficient formula and completeness.** Lagrange inversion yields precisely `c_k=[y^k]Q(y)^(-(k+1)/r)`. It uses the germ through degree `r+k`; the last divergent coefficient therefore uses degree `2r-2`. Distinct strictly divergent powers prove the bounded-difference equivalence for germs of the same even order. The last-jet perturbation has the displayed nonzero epsilon^(-1/r) coefficient. The separate quadratic sharpness statement is necessary and correct.

5. **Examples and excluded extensions.** All quartic and sextic constants and exponents agree with independent exact derivations. The example `y^4+y^5+(7/8)y^6` cancels the second divergent coefficient while changing the degree-six jet. The one-sided pure quadratic retains `+(1/2)log(1/epsilon)` through the corrected-clock endpoint term, as stated. No inference about integer rounding, residues, primes, or global Collatz entry follows from bounded absolute error.

## Independent exact controls

`independent_controls.py` imports no author code. It derives generic coefficients for orders 2, 4, 6, 8, 10, and 12 in two separate ways:

- the triangular coefficient recurrence obtained from `Q B'=alpha Q' B`;
- a scaled geometric expansion of the original denominator, with beta moments indexed by integer partitions.

These agree for every possible divergent coefficient in the tested orders. Further exact checks cover the printed quartic/sextic formulas, above-cutoff independence, last-jet sharpness, embedded-flow jets, borderline primitives, and tail integrability.

A useful additional hostile test verifies the coordinate-family restriction. For `phi(y)=y+c y^2`, conjugating `y+y^4` gives the increment

    u^4 - 2c u^5 + 4c^2 u^6 + O(u^7).

Reimposing a fresh constant-additive epsilon in this coordinate gives `c_2=-3c^2/8`, hence `A_1=-3*pi*c^2/(8*sqrt(2))`; the original pure quartic has `A_1=0`. Thus the coefficient vector is not an invariant of the unperturbed conjugacy class. This agrees with, rather than weakens, the note's stated scope.

For the positive quartic control on `J=[-1/4,1/4]`, the checker also gives exact coefficient-majorant certificates: `Q>=89/128`, `eta=357/4096<1/4`, `B=1193/1024`, and `C=127449/45568`. These imply the complete-step domain and corrected-clock bounds without numerical grid sampling.

## Embedding contrast and source verification

The explicit embedded germ is the time-one map of `a y^r partial_y`; its first increment difference from `a y^r` occurs in degree `2r-1`, beyond the divergent-clock cutoff. The companion cotangent theorem was inspected: it applies to a nonidentity holomorphic parabolic germ and obtains convergence of the normalized iterative logarithm from a nonconstant meromorphic invariant using its first nonzero Laurent coefficient. These hypotheses hold for both examples.

The polynomial nonembedding claim is supported by the classical result explicitly recorded on printed page 3 of [Aschenbrenner–Bergweiler](https://www.mat.univie.ac.at/~maschenbrenner/pdf/itlog-final.pdf): the convergence classification applies to functions meromorphic on the whole complex plane. The polynomial qualifies; for higher order the embedded algebraic local germ need not be globally meromorphic. The note correctly cites this classification rather than inferring divergence from differential transcendence alone. The comparison concerns epsilon zero, not flow embedding of the perturbed family.

Also checked: [DLMF Lagrange inversion](https://dlmf.nist.gov/1.10#vii), [DLMF beta integrals](https://dlmf.nist.gov/5.12), [Kuehn Proposition 2.1](https://arxiv.org/html/0807.1546), and the broader context in [Brockmoeller–Scherz–Srkalovic](https://arxiv.org/html/2505.07138v1). The note uses these within their actual roles. Both repository dependency hashes match `SOURCE-DEPENDENCIES.json`.

One initial companion-directory link was incorrect; the author corrected it to `parabolic_cotangent_integrability_20261019` before the hashes above were captured. No unresolved correction remains.

## Reproduce

From this review directory:

    python3 independent_controls.py
    python3 -O independent_controls.py

The script needs Python 3 and SymPy and writes `independent-results.json` beside itself. It has no dependency on a checkout location. In a publication bundle with this directory named `review`, replay the author verifier from its parent directory using the commands in the parent's README. Numerical orbit counts, including Decimal cross-checks, are supporting controls only; no empirical bound or finite coefficient list is offered as a proof of the analytic claims.
