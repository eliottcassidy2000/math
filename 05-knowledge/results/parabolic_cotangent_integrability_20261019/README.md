# Parabolic clocks, cotangent Hamiltonians, and the missing analytic invariant

2026-10-19. **PROVED, scoped:** formal invariant-ring and meromorphic embedding criteria, the full rational-invariant obstruction for polynomial cotangent lifts, a finite-jet hostile family, and the exact period-three residue 51/49. **CITED:** divergence/differential transcendence of polynomial iterative logarithms and the existence of sectorial Fatou coordinates. No historical-priority claim is made.

The quadratic cusp recurrence becomes

    q_(n+1)=q_n+q_n²+epsilon,  q=x−1/2.

At epsilon=0 its canonical cotangent lift is

    C_f(q,p)=(q+q², p/(1+2q)).

Its formal Hamiltonian logarithm H=p·itlog(f) exists, but it has no nonconstant meromorphic first-integral germ. This is a sharper category distinction than the absence of a polynomial Hamiltonian. The same cotangent mechanism applies to the parabolic period-three return at c=−7/4; its exact time-form residue is 51/49, whereas the basic cusp residue is 1.

Read [parabolic_cotangent_invariant_barrier.md](parabolic_cotangent_invariant_barrier.md) for the proofs, scope and literature. The separate extension [formal_quadratic_unfolding_period.md](formal_quadratic_unfolding_period.md) identifies the exact formal singular time-form coefficients A(delta), B(delta) from the multipliers 1±2i delta; it predicts the pi delta/3 correction while leaving analytic remainder estimates to the real escape proof. Run `python3 verify_parabolic_geometry.py` and `python3 verify_formal_unfolding.py` for exact controls. Finite coefficient checks support identities, not divergence or all-degree nonintegrability.

## Verification and reproduction

The [independent review](review/REVIEW.md) approves the proofs within their stated scope. The main checkers pass 186 + 15 exact checks; the independent standard-library checker passes 129. Normal and optimized Python runs give identical output and receipts. The main checkers require Python 3 and SymPy (verified with SymPy 1.14.0); the independent checker uses only the Python standard library.

Run from the repository root:

```bash
python3 05-knowledge/results/parabolic_cotangent_integrability_20261019/verify_parabolic_geometry.py
python3 -O 05-knowledge/results/parabolic_cotangent_integrability_20261019/verify_parabolic_geometry.py
python3 05-knowledge/results/parabolic_cotangent_integrability_20261019/verify_formal_unfolding.py
python3 -O 05-knowledge/results/parabolic_cotangent_integrability_20261019/verify_formal_unfolding.py
python3 05-knowledge/results/parabolic_cotangent_integrability_20261019/review/independent_controls.py
python3 -O 05-knowledge/results/parabolic_cotangent_integrability_20261019/review/independent_controls.py
(cd 05-knowledge/results/parabolic_cotangent_integrability_20261019 && sha256sum -c SHA256SUMS)
```

The [initial real escape theorem](../parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md) and [next-order real escape theorem](../parabolic_escape_next_order_20261019/parabolic_escape_next_order_20261019.md) provide the separately proved transit estimates. This packet supplies their formal-geometric interpretation.

[Publication overlays](publication-overlays.json) record relative-link, reproduction-path and presentation-only changes from the independently reviewed sources. The theorem statements and symbolic computations are unchanged.
