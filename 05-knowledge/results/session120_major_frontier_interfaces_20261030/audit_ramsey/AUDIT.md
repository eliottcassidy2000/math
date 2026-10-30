# Independent audit of the valuation-preserving spectral compiler

Status: **PASS**, 2026-10-30. This audits `../ramsey/SPECTRAL_COMPILER.md`; it does not establish a Ramsey record or a circuit lower bound.

## Analytic audit

1. With `S[a,b]=sigma(b-a)` and positive-exponent Fourier vectors, the eigenvalue uses the positive exponential, as printed. No unnoticed conjugation occurs.
2. CRT at modulus `m` requires the local inverse `(m/p^a)^(-1) modulo p^a`. The displayed local argument contains this inverse. Its prime-adic valuation is zero, so it preserves the stated valuation tests.
3. Summing lifts of a nonzero residue modulo `p` gives zero unless `p^(a-1)` divides the frequency. If divisible, the twisted factor is zero at exact frequency zero modulo `p^a`, and otherwise is `p^(a-1) chi_p(h/p^(a-1)) tau_p`. Thus nonunit/valuation strata are retained, not discarded.
4. The Gauss-sum norm proof is correct, and conjugation is `-tau_p` since `p=3 mod 4`. Consequently `tau_p^2=-p`. The calculations need no choice of a numerical square-root sign.
5. The square classes `-p` are independent in `Q*/Q*²`: in any nonempty product, each selected prime has odd valuation. The multiquadratic monomials therefore constitute a basis; constant-coefficient extraction equals the normalized field trace.
6. At frequency `N/d`, a twisted local term can survive only when `A_p-v_p(d)=a-1` with `a<=A_p`; hence its prime divides `d`. The eigenvalue is in the multiquadratic subfield with generators indexed by `p|d`. Units modulo `d` realize all independent Legendre signs equally often by CRT. This justifies the weight `phi(d)` multiplying a constant coefficient, rather than an unnormalized field trace. The zero frequency vanishes. Odd powers have no constant coefficient because the eigenvalue lies in the odd monomial subspace.
7. Squaring `a tau3+b tau7+c tau11+d tau3 tau7 tau11` yields constant `-(3a²+7b²+11c²+231d²)` and pair coefficients `2ab-22cd`, `2ac-14bd`, `2bc-6ad`. Their squared basis elements are `21,33,77`, verifying the quartic formula and every cross-term sign.

## Independent finite check

`check_independent.py` imports no producer modules. It enumerates odd truth tables using a binary cube and antipodal XOR pairing, rather than using the producer's response indices or character/majority formulas. It classifies signed characters by normalized multiplicativity on the cube. It builds actual positive arc sets separately on every additive-order stratum using literal Legendre values.

For every one of the 8192 N=231 symbols it computes cyclic support intersections. For nonzero shift `h`, expanding `sigma(x)=2*1_B(x)-1+1_{x=0}` gives

`C(h)=4|B intersect (B+h)|-N+2`, and `C(0)=N-1`.

Since `S` is skew and `SS^T` is the symmetric correlation circulant, `tr(S^4)=N sum_h C(h)^2`. The independent result agrees with **each complete histogram** in the producer receipt, not only the extrema:

| Class | Count | Minimum | Minimum multiplicity | Maximum | Distinct moments |
|---|---:|---:|---:|---:|---:|
| Locally signed characters | 4096 | 17225670 | 16 | 645693510 | 468 |
| Nonlinear unit response | 4096 | 29082438 | 8 | 175030086 | 429 |

Both independently found minimizing masks also pass literal double-loop scalar autocorrelation. All 52 symbols in the seven controls `N=3,7,9,21,27,33,49` pass literal scalar autocorrelation, and their complete histograms agree with the producer receipt. In particular the N=21 values are `15876,19908,23940,44100`, each with multiplicity four.

Normal and optimized executions are byte-identical:

```
python session120_research/audit_ramsey/check_independent.py > session120_research/audit_ramsey/normal.json
python -O session120_research/audit_ramsey/check_independent.py > session120_research/audit_ramsey/optimized.json
cmp session120_research/audit_ramsey/normal.json session120_research/audit_ramsey/optimized.json
```

Both receipts have SHA256 `cea2dc9e4e8c88ec0a7c82b30ec789686ae9b00ebbf736635e8d7ec44de13522`. Explicit exceptions, not removable assertions, enforce checks.

## Scope and remaining boundary

The all-parameter theorem is supported by the analytic audit; the finite census is independently verified. Agreement of complete histograms does not establish per-index correspondence between the two deliberately different enumerations. It does establish the stated family extrema and multiplicities. This is a comparison of the fourth spectral moment only: no implication for optimal transitive-subtournament size or a Ramsey-number improvement has been proved. No novelty or priority claim beyond the producer's stated scope was audited.
