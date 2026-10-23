# Independent Gaussian/phase compiler audit

Audit date: 2026-10-23. Source packet: `session43-gaussian-jacobian`. Final reviewed manifest SHA-256: `c417b8986ce5eadecb731fb6ce8a54a2e318cf99ffbd304f4867a09bfc867c93`.

## Verdict

PASS for the elementary scoped theorems, with one minor scope repair made by the producer during review: the cyclic shear is nonlinear only for L >= 3. L=2 is a linear, non-torus-equivariant shear. No source or Git mutations were made by this auditor. Both independent test modes pass 529 checks with byte-identical output. General theorems rest on the arguments below, not those finite checks.

## Proof review

1. Finite Reynolds null powers. The exact hypothesis is a reduced commutative Q-algebra, not merely any ring of characteristic zero. Newton identities are valid with all relevant integer denominators invertible. Vanishing of the first N orbit power sums implies the orbit characteristic polynomial is T^N. Substituting f, one orbit root, gives f^N=0; reducedness gives f=0. Zero divisors do not invalidate this argument. The C_N character on C[z] demonstrates uniform sharpness, including N=1. A bare characteristic-zero reduced ring can have positive-characteristic components and is not an acceptable replacement for a Q-algebra.

2. Pointwise finite orbit. Over C, grouping distinct nonzero values gives an invertible Vandermonde system. Positive weights are sufficient; nonzero grouped weights actually suffice. Signed weights can cancel at the same value, so the distinction from signed cancellation is meaningful. The nonreduced C2 action on dual numbers is a valid failure control.

3. Angular compilation. The phase operator is diagonal on integral charges. Intersection of divisibility conditions gives the least common multiple, not product in general; for 8 and 9 it is 72. Exact equality with continuous averaging requires no surviving nonzero charge multiple of L. The support-sum bandwidth criterion is sufficient, not necessary where coefficients cancel. The source correctly states this distinction.

4. Radial normalization. For Z=(X+iY)/sqrt(2) with independent standard real Gaussians, |Z|^2 is exponential with mean 1 and its kth moment is k!. Finite 72-angle sampling of Z^72 therefore produces E[r^72]=36!, whereas the continuous angular expectation is zero. This is not a claim that an algebraic radial functional is canonically defined on all aliased charged polynomials. The source correctly separates radial sampling from the invariant-ring contraction.

5. Shared phase. Z1 W2 has total charge zero, but individual vector charges +1 and -1. Simultaneous averaging retains it; independent averaging kills it. Setting the vectors equal is an explicit correlated specialization and gives expectation 1. One must not claim expectation 1 for two independent centered Gaussian vectors. The source does not make that error.

6. Exact torus rigidity. For positive integer w, weight-(1,-w) equivariance forces (x A(t),y B(t)), t=x^w y. The determinant equals AB+t AB'+w t A'B. For nonzero A,B of degrees a,b, its highest coefficient is (1+b+wa) times the product of their leading coefficients. In characteristic zero this cannot vanish. A nonzero constant determinant forces a+b=0. This gives an especially direct version of the source's quotient-degree proof; no negative-power ambiguity is needed. The conclusion is diagonal linearity for this exact equivariance, not JC2.

7. Finite symmetry control. The corrected L>=3 shear has determinant one, a polynomial inverse, and the right congruence weight. It violates exact integer weight, proving that finite residue guards alone do not furnish the rigidity argument. For L=2 it remains non-torus-equivariant but is linear.

## Inherited claims and primary-source scope

The supplied THM-2022 snapshot and the repository GMC2Main.lean front door were inspected. The front door delegates its gmc2 theorem, with only the moment-null premise, to GMC2DvdKOmegaWiring.gmc2_unconditional. This is an inspected repository claim. Neither the producer nor this auditor independently rebuilt Lean, inspected all transitive proof obligations, or reproduced the axiom report. No GMC2 resolution is certified by this audit.

Primary source checked: Christopher D. Long, https://arxiv.org/html/2607.18186v1, Theorem 5.1, Corollary 5.2, Section 7 (accessed 2026-10-23). Its three-variable polynomial has vanishing pure moments and mixed moments m!; unused-variable padding extends failure to larger dimensions. Section 7 explicitly leaves dimension two unsettled. This verifies the producer's citation scope, not the paper's unrelated historical claims or Jacobian discussion. The independent checker directly contracts the displayed three-variable polynomial for orders 1 through 8. Its all-order identity remains attributed to the cited theorem, not newly proved by finite checking.

## Reproduction

Run `python3 session43-audit-gaussian/check.py` and `python3 -O session43-audit-gaussian/check.py`. Both produce the committed receipt exactly. SymPy is required. The script first checks all producer manifest file hashes, uses explicit RuntimeError guards, checks actual cyclotomic sums, discrete convolution multiplicities, a separate real-Gaussian radial expansion, direct Jacobians, and direct Wick contractions. The audit is tied to the producer manifest above. Files in this audit packet are covered by MANIFEST.json and SHA256SUMS; those two inventory files are intentionally not recursively self-hashed.

## Root bounded-degree repair extension

Also independently reviewed `session43-root/FINITE_GUARDS_AND_DEGREE.md` and its `check_degree_guards.py`; exact source hashes are in `degree.receipt.json`. Verdict PASS. For total degree D, each input charge lies in [-D,D], so difference from target +1 or -1 is at most D+1 in absolute value. L>D+1 prohibits any nonzero congruence alias. Exact equivariance and constant nonzero determinant give diagonal linearity as proved above. The least degree of a wrong-charge monomial in either output congruence class is exactly L-1 for every L>=2. Thus this is a sharp universal threshold, not merely a loose sufficient bound. At L=2 the boundary shear is linear, as the root explicitly notes. Combined exact C8/C9 actions give C72 and diagonal linearity through degree70; the degree71 shear shows the boundary. No arbitrary mod8/mod9 tournament-statistic guard is interchangeable with those actions.

The extension checker independently searches degree layers for the first wrong-charge monomial for both targets and all L=2,...,120, and checks the combined action through degree74. Normal and optimized outputs agree. This extends the audit but does not alter the producer's packet or the root's files.

## Publication portability

Both checkers support the original workspace layout and the published layout, where `audit-gaussian/` and `gaussian-jacobian/` are sibling directories and root degree files sit in their parent. A temporary mock publication tree was tested explicitly in both normal and optimized modes; all four outputs matched the workspace receipts byte for byte. This repair changes path resolution only, not mathematical claims or checks.

Final source repin: the producer normalized one copied relative Markdown link to a commit-pinned GitHub URL and documented the original snapshot digest. No theorem changed. All producer file hashes and independent checks were rerun against this final source manifest.
