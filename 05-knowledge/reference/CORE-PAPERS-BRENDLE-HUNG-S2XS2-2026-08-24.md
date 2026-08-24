# Brendle--Hung `S2 x S2` positive-curvature preprint: source and referee ledger -- 2026-08-24

## Truth boundary

Simon Brendle and Pei-Ken Hung, [*A metric on `S2 x S2` with positive
sectional curvature*](https://arxiv.org/abs/2608.19068), arXiv:2608.19068v1,
submitted 2026-08-19, is **PREPRINT CLAIM / UNDER AUDIT**.

The paper claims a smooth metric of positive sectional curvature on
`S2 x S2`.  If correct, this settles the positive-curvature question often
called the Hopf product problem.  It is distinct from Hopf's problem asking
whether `S6` admits an integrable complex structure.

This repository does **not** promote the headline theorem.  The paper and its
companion notebook contain several concrete transcription/state errors.  Most
are locally repairable, but one identity used in the load-bearing quadratic
positivity calculation is not certified by the attached notebook as saved.
The correct status is therefore:

```text
headline existence theorem                 PREPRINT CLAIM / UNDER AUDIT
abstract minimization framework, Section 2 READ / NO IMMEDIATE CONTRADICTION
displayed perturbation mechanism            READ / NO IMMEDIATE CONTRADICTION
V_bc^(2)=0                                  OPEN; exact point + conditional symmetry route
z(h_a^(1))=0                                VERIFIED-EXACT on both displayed Sigma pieces
other notebook identities                   PARTIALLY AUDITED, not imported
Poisson repair principle                    PROVED independently in THM-3990
```

The source PDFs and archive used for this intake had hashes

```text
PDF SHA-256     46bb66824fa239b95dbfa6a661de675a2d15d0a1e38332433827bf9f31f4e4f7
source SHA-256  81cf36b7fa6f05fcb7f1dbf825cd391b2b95fae71fa016d73f21e871c848d2b1
notebook SHA-256 c607493b837ff3a1a73b74dfdc2337e1a45bb7290e5ac3fc59d98a3f12c82c58
```

The source archive contains `counterexample.tex` and
`Mathematica_Code_Companion.nb`.  All 34 PDF pages were text-extracted; pages
1, 11, 24, 26, and 28 were also rendered and inspected.

## 1. Claimed mechanism

Let `M=S2 x S2`.  The construction begins with a Cheeger--Mueter metric `g`
of nonnegative sectional curvature.  On

```text
M_generic=M minus (Delta_+ union Delta_-),
```

there is exactly one zero-curvature two-plane at each point.  Those planes
form a smooth four-dimensional submanifold `Z` of the eight-dimensional
Grassmann bundle.

The paper perturbs

```text
g_s=g+s h^(1)+s^2 h^(2)+s^3 h^(3).
```

After minimizing over transverse plane coordinates, the minimum curvature has
the stated expansion

```text
s^2 V^(2)+s^3 V^(3)+O(s^4).
```

The intended proof spine is:

1. the first coefficient vanishes on `Z`;
2. `V^(2)>=delta(1-<q_3,e_z>^2)` and its residual zero set is the completed
   two-torus `Sigma`;
3. the integral of `V^(3)` over `Sigma` is nonzero for a generic parameter
   choice;
4. solve `V^(3)-Delta_Sigma chi=mu`, where `mu` is its average, and take the
   conformal third-order correction `h^(3)=6 chi g`;
5. choose the sign of `s` from the sign of `mu`, then use the compact abstract
   framework and density of `M_generic` to obtain strict positivity globally.

The last Poisson step is conceptually complete once the preceding curvature
coefficients are known.  [THM-3990](../../01-canon/theorems/THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md)
proves the exact general statement: on a connected closed residual locus the
mean is the complete Laplacian cokernel, while on a disconnected locus every
component mean must have one common nonzero sign.

The two open pieces of `Sigma intersect M_generic` are not separate closed
components: they join through `Delta_+ union Delta_-` in the completed torus.
Their individual open-piece integrals have boundary flux and are not Poisson
invariants.  The relevant invariant is the mean over the completed `Sigma`.

## 2. Paper transcription errors

The following are visible in arXiv v1's TeX source.

| source location | printed text | required repair | role |
|---|---|---|---|
| TeX 658 | an unmatched `\rangle` in `h_a^(1)` | delete the stray bracket | display transcription |
| TeX 881 | `V_cd^(2)` uses `L h_bd^(2)` | use `L h_cd^(2)` | unique index repair |
| TeX 885--886 | the expansion ends in `2 lambda_d V_dd^(2)` | use `2 lambda_c lambda_d V_cd^(2)` | load-bearing displayed expansion |
| TeX 917--919 | equation labelled `V_bd^(2)=0` prints `V_bb^(2)=0` | print `V_bd^(2)=0` | label transcription |
| TeX 928--930 | from `|V_dd^(2)|<=C rho`, choose `lambda_0=sqrt(C/delta)` | require `lambda_0^2 C<=delta`, e.g. `lambda_0=sqrt(delta/C)` | load-bearing inequality, immediate local repair |

The last correction follows directly from the intended estimate

```text
V^(2) >= 2 delta rho-lambda_d^2 C rho.
```

Thus the printed choice has the ratio inverted unless `C=delta`.  This does
not expose a conceptual obstruction, but the proposition is not correct as
printed.

## 3. Companion-notebook state errors

### 3.1 The missing `V_bc^(2)` verification

In the raw notebook, lines 8240--8244 compute

```text
rb=P1[L[hb]];
rc=P1[L[hc]];
zc=-Inverse[H].rc;
Vbc=P0[L[hbc]]+P0[Q[hb,hc]]+1/2 rb.zc;
Vbc=FullSimplify[Vac,0<theta<Pi/2];
```

The fifth line overwrites the just-computed expression with the already
computed `Vac`.  The later cell at notebook lines 8269--8274 reports
`Vbc==0` as `True`, but it is therefore a second check of `Vac==0`, not a
check of the intended `Vbc` expression.

This is load-bearing.  The paper's quadratic coefficient contains
`2 lambda_c V_bc^(2)`, and the claimed lower bound in Proposition 5.2 uses
`V_bc^(2)=0`.  No later cell consumes the corrupted variable, so the error is
localized, but localization does not supply the missing identity.

The independent exact companion
[`brendle_hung_vbc_exact_point_audit_20260824.py`](../../04-computation/brendle_hung_vbc_exact_point_audit_20260824.py)
rebuilds the three terms at one generic rational point. All three summands are
nonzero and cancel exactly. This is a sensitive **FINITE-EXACT point control**,
not a proof of the identity.

There is also a conditional symmetry route. Let `R` be the quarter-turn about
`e_z` and put `Phi(p1,p2)=(R p1,-R p2)`. Direct tensor transformation gives

```text
Phi^*h_b=-h_a,       Phi^*h_c=h_c,       Phi^*h_bc=-h_ac.
```

Naturality of the curvature/minimization operators then gives
`Phi^*V_bc=-V_ac`. Thus a clean independent global verification of the
source's `V_ac=0` would also prove `V_bc=0`. The saved notebook's `Vac` check
is source evidence, not that independent verification, so the obligation
remains open.

### 3.2 The missing `z(h_a^(1))=0` verification

At notebook lines 10384--10389, the Lemma 5.4 block begins

```text
ra=P1[L[hc]];
rb=P1[L[hb]];
rc=P1[L[hc]];
```

Thus `ra` repeats `hc` instead of using `ha`, and the subsequent `za` check
does not certify `z(h_a^(1))=0`.  Earlier notebook blocks did compute
`ra=P1[L[ha]]`, but the saved Lemma 5.4 block does not provide a clean
fresh-state certificate for the stated lemma.

This error is also localized.  The later `A1`--`A4` cells defining the
coefficient of `lambda_c lambda_d^2` independently rebuild their `c,d`
objects and do not consume the corrupted `ra` or `za`.  Moreover the
`lambda_c lambda_d^2` monomial only needs the checked `c,d` pieces.  Therefore
the saved-state error does not by itself refute the later displayed cubic
coefficient. The independent exact companion
[`brendle_hung_lemma54_independent_audit_20260824.py`](../../04-computation/brendle_hung_lemma54_independent_audit_20260824.py)
now verifies `P1(L(h_a))=0` symbolically on both displayed open pieces of
`Sigma`, with four exact point controls and a nonzero hostile tensor. This
closes the omitted `h_a` branch only; the cubic coefficient still awaits an
independent clean-kernel audit.

## 4. Exact remaining referee obligations

The shortest path to a trustworthy v1 certificate is:

1. **Fresh-kernel quadratic audit.** Recompute `Vbc` globally without
   assignment reuse, or independently rebuild `Vac=0` and the displayed
   symmetry. The one-point exact cancellation is only a positive control.
2. **Complete Lemma 5.4 audit.** The omitted `z(ha)=0` branch is now
   independently exact; retain or rebuild the source's `hb,hc` branches when
   certifying the full lemma on completed `Sigma`.
3. **Rebuild all ten `V^(2)` summands.** Use the corrected `h_cd`, corrected
   cross term, and a valid smallness constant
   `lambda_0<=sqrt(delta/C)`.
4. **Rebuild the cubic witness.** Evaluate the four `A1`--`A4` contributions
   in a clean state and independently recover the stated integral
   `pi^2/(18 sqrt(3))`.
5. **Global smoothness audit.** Check that every infinite-series coefficient
   used in `h^(2)` extends smoothly across `Delta_+ union Delta_-`, with the
   convergence and differentiation needed by the curvature calculation.
6. **Abstract-to-geometric interface.** Verify all hypotheses of the compact
   minimization theorem for the corrected global tensors, including the
   uniform remainder and density statements.

Running the notebook top-to-bottom without fixing the two assignment errors
is not an independent path.  A repaired notebook should clear its namespace
or use immutable names for every asserted identity.

## 5. Transferable themes and new frontiers

The useful mathematical object is not “third order” by itself.  It is the
pair

```text
(residual zero stratum, cokernel of the allowed correction operator).
```

This suggests four bounded extensions.

1. **Nested residual strata.** If the second coefficient vanishes on `S_2`
   and the repaired third coefficient vanishes on `S_3 subset S_2`, repeat
   the construction with the operator induced on `S_3`.  The next obstruction
   belongs to a new cokernel and cannot be inferred from the previous mean.
2. **Weighted and equivariant repairs.** For a correction operator `A`, replace
   component averages by the pairing with `ker(A^*)`.  Symmetry can enlarge
   that adjoint kernel, so averaging only constants may miss an obstruction.
3. **Singular residual loci.** When components meet, branchwise means require
   interface/flux conditions.  The proper analogue is a graph or conductor
   Laplacian, retaining integral Smith torsion when gluing is discrete.
4. **Proof-carrying symbolic calculations.** Export each load-bearing identity
   as a separately named exact expression with dependency hashes and hostile
   nonzero controls.  This is a better target than one mutable notebook state.

The third extension is the precise bridge to the `S6` clutch and the planar
Jacobian conductor work.  It is a bridge of repair quotients, not evidence
that either preprint's headline claim is true.

Reproduce the two independent controls from the repository root:

```text
python3 04-computation/brendle_hung_lemma54_independent_audit_20260824.py
python3 -O 04-computation/brendle_hung_lemma54_independent_audit_20260824.py
python3 04-computation/brendle_hung_vbc_exact_point_audit_20260824.py
python3 -O 04-computation/brendle_hung_vbc_exact_point_audit_20260824.py
```

Their script/output SHA-256 pairs are respectively
`2b8295e8...ef4d / c504027f...d655` and
`45f7c0b6...1102 / 3c6dcf70...12dc`.
