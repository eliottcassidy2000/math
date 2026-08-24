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
are locally repairable.  Two identities are not certified by the attached
notebook as saved; the targeted independent symbolic repairs below now close
those two identities relative to a reconstructed implementation, without
replaying the remaining quadratic or headline calculation.  The correct
status is therefore:

```text
headline existence theorem                 PREPRINT CLAIM / UNDER AUDIT
abstract minimization framework, Section 2 READ / NO IMMEDIATE CONTRADICTION
displayed perturbation mechanism            READ / NO IMMEDIATE CONTRADICTION
V_ac^(2)=V_bc^(2)=0                         VERIFIED EXACT GLOBAL SYMBOLIC IDENTITIES
V_bc^(2) at one generic exact point         FINITE-EXACT positive control
P1(L(h_a^(1)))=0 on Sigma intersect M_generic VERIFIED EXACT SYMBOLIC IDENTITY
generic normal Hessian positive definiteness VERIFIED EXACT SYMBOLIC IDENTITY
z(h_a^(1))=0 on Sigma intersect M_generic   VERIFIED EXACT SYMBOLIC IDENTITY
other notebook identities                   PARTIALLY AUDITED, not imported
Poisson repair principle                    PROVED independently in THM-3990
```

The source PDFs and archive used for this intake had hashes

```text
PDF SHA-256     46bb66824fa239b95dbfa6a661de675a2d15d0a1e38332433827bf9f31f4e4f7
source SHA-256  81cf36b7fa6f05fcb7f1dbf825cd391b2b95fae71fa016d73f21e871c848d2b1
notebook SHA-256 c607493b837ff3a1a73b74dfdc2337e1a45bb7290e5ac3fc59d98a3f12c82c58
TeX SHA-256     9d123f843555ed49b38c17d651867fd7bd94d0e1da851025c673f97bde62c800
global audit script SHA-256 9939f4fdb2ca7a7a6c639747de753d3bfe0b7827c7a1579bc66dd862ae13e772
global audit output SHA-256 9dc962b4a2de1f8b64706a83826e6bf9030d23f2e9b3353c26513d68ee036b2b
Sylvester audit script SHA-256 69bfb466ffd26ad88465e4fc9f6adb71adfc0bf621e3a9fcd649d45f7035677d
Sylvester audit output SHA-256 82cae029b85efddcbd6404585cb8d52483944c2f9a3e2de87fd990f937f6d999
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
| TeX 320, 324 | the set `S` is used but never defined | replace `S` by the defined subset `Lambda` | notation repair |
| TeX 476--477 | the intermediate `g_std,t` display has `+d theta tensor d theta` | use `+2 d theta tensor d theta`; the final scaled metric at 482--484 is already consistent | non-load-bearing intermediate inconsistency |
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

As saved, this is load-bearing.  The paper's quadratic coefficient contains
`2 lambda_c V_bc^(2)`, and the claimed lower bound in Proposition 5.2 uses
`V_bc^(2)=0`.  No later cell consumes the corrupted variable, so the error is
localized.  The symmetry reduction, exact point control, and independent
global reconstruction below now close this identity relative to the displayed
source formulas.

First, there is a global symmetry argument.  Let `R` be the oriented
quarter-turn

```text
R e_x=e_y,       R e_y=-e_x,       R e_z=e_z,
```

and put

```text
F(p_1,p_2)=(R p_1,-R p_2).                               (1a)
```

The map is an isometry of the background Cheeger metric.  If `alpha_A` is the
diagonal `SO(3)` action, then

```text
F o alpha_A=alpha_(R A R^(-1)) o F,
F_*K_omega=K_(R omega),       F_*J_omega=J_(R omega),
<F(p_1),F(p_2)>=-<p_1,p_2>.                               (1b)
```

Writing the target point with primes, one also has

```text
R(p_1+p_2)=-(p'_2-p'_1),
R(p_2-p_1)=-(p'_1+p'_2).                                  (1c)
```

Substitution in the displayed tensor definitions gives

```text
F_*h_a^(1)=-h_b^(1),
F_*h_c^(1)= h_c^(1),
F_*h_ac^(2)=-h_bc^(2).                                    (1d)
```

In the paper's zero-plane charts this also gives

```text
F_*m_1=-m'_2,   F_*m_2=-m'_1,   F_*m_3=-m'_3,   F_*m_4=-m'_4,
z'=(z_3,z_4,z_1,z_2).                                    (1e)
```

The bi-invariant `SO(3)` metric used in the Cheeger construction is
`Ad`-invariant, while the curvature operators `L` and `Q(h,k)` are natural
under isometries.  Together with transport of the zero-plane bundle, its
normal Hessian, and its minimizing normal vector, this makes the curvature
expansion and normal minimization natural under `F`.  Therefore

```text
V_ac o F^(-1)=-V_bc.                                      (1f)
```

This is the **pushforward** convention.  Direct pullback by `F` instead gives

```text
F^*h_a^(1)=h_b^(1),   F^*h_c^(1)=h_c^(1),
F^*h_ac^(2)=h_bc^(2),   hence   V_ac o F=V_bc.             (1f')
```

There is no sign discrepancy: `F_*=(F^(-1))^*`.  The inverse quarter-turn
has the two minus signs in (1d), whereas the forward pullback has the plus
signs in (1f').

The notebook's preceding `V_ac=0` block is not overwritten.  Conditional on
that saved exact calculation, (1f) proves `V_bc=0` throughout `M_generic`.
This is a conceptual repair of the assignment error, not an independent
replay of the `V_ac` block.

Second, the
[exact point audit](../../04-computation/brendle_hung_vbc_exact_point_audit_20260824.py)
is independent of the saved Mathematica state.  It rebuilds the paper's
metric, connection, Riemann tensor, `L`, `Q`, Hessian, and normal minimizer at
`theta=pi/6` and the rational `SO(3)` frame obtained from quaternion
`(1,2,3,4)`.  It shares the same reconstructed background-geometry
implementation as the `h_a` audit below, so it is not a second independent
implementation of that geometry.  It finds the nontrivial cancellation

```text
P0 L(h_bc)             =  581/20250,
P0 Q(h_b,h_c)          = -1064/50625,
(1/2) r(h_b).z(h_c)    = -259/33750,
V_bc                   = 0.                               (1g)
```

All three summands are nonzero, so the control detects a missing or duplicated
term.  Its [frozen output](../results/brendle_hung_vbc_exact_point_audit_20260824.out)
LF-normalized byte-matches normal and optimized Python. This is
**FINITE-EXACT at one
generic point**, not by itself an identity on `M_generic`; the conditional
generic-locus repair is (1f).

Finally, the independent
[global symbolic audit](../../04-computation/brendle_hung_vac_vbc_global_identity_audit_20260824.py)
reconstructs both mixed coefficients and then sets `x=tan(theta)>0`.  It
is independent of the saved Mathematica state but shares the repository's
reconstructed moving-frame implementation with the point and `h_a` audits;
it is not a second implementation of that geometry.  It first obtains

```text
det(H)=256 x^2 (x^2+1)^6 (3x^2+5)(5x^2+3)
       / (81 (x^2+3)^5 (3x^2+1)^5) > 0,                 (1h)

z(h_c)=(-q_1z q_3z x/(3 sqrt(x^2+1)), 0,
        -q_2z q_3z/(3 sqrt(x^2+1)), 0).                 (1i)
```

Writing `B=q_2z q_3y q_3z` and `D=(1+x^2)^(5/2)`, the three independently
reduced `V_bc` summands are

```text
-B x(x^2-1)(2x^2+27)/(15D),
 B x(6x^6+159x^4+150x^2-203)/(45D(x^2+3)),
-2B x(33x^4+3x^2+20)/(45D(x^2+3)).                     (1j)
```

Their common-denominator numerator is the literal polynomial identity

```text
-3(x^2-1)(2x^2+27)(x^2+3)
 +(6x^6+159x^4+150x^2-203)
 -2(33x^4+3x^2+20)=0.                                  (1k)
```

The audit independently gives the analogous three-term identity for
`V_ac=0`. No algebraic `SO(3)` constraint on the final `q`-components is used:
the sums vanish with those components treated as formal variables. The
background reconstruction itself still uses the `SO(3)` Lie-frame brackets.
Dropping or flipping the
minimizer term leaves a recorded nonzero polynomial, and the exact symmetry
check verifies both (1f) and (1f').  The
[frozen output](../results/brendle_hung_vac_vbc_global_identity_audit_20260824.out)
LF-normalized byte-matches normal and optimized Python and records script and
dependency hashes. Thus `V_ac=V_bc=0` on `M_generic` is now a **VERIFIED EXACT GLOBAL
SYMBOLIC IDENTITY relative to the reconstructed source formulas**.  This
closes the overwritten assignment, not the other quadratic summands or the
headline theorem.

A genuinely separate
[Sylvester audit](../../04-computation/brendle_hung_normal_hessian_sylvester_independent_audit_20260824.py)
rebuilds the moving-frame geometry directly in `x`, without importing the
repository geometry helper, and uses the source `P2` coordinate order
`(z1,z2,z3,z4)`. Its four leading principal minors are

```text
Delta1=2(x^2+1)^2(3x^2+7)/(3x^2+1)^3,
Delta2=16(x^2+1)^3(3x^2+5)/(9(x^2+3)^2(3x^2+1)^3),
Delta3=32(x^2+1)^5(3x^2+5)(7x^2+3)
       /(9(x^2+3)^5(3x^2+1)^3),
Delta4=det(H) from (1h).                                (1l)
```

All are strictly positive for `x>0`, so Sylvester's criterion proves that the
reconstructed generic normal Hessian is positive definite. At `x=1` the
minors are `(5/4,1/9,5/36,1/81)`. The determinant tends to zero at both
deleted endpoints, so this gives neither endpoint extension nor uniform
coercivity. A determinant-one indefinite hostile with `Delta2=-1` confirms
that the earlier determinant check alone did not prove positivity. The
[frozen output](../results/brendle_hung_normal_hessian_sylvester_independent_audit_20260824.out)
LF-normalized byte-matches normal and optimized Python.

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
the saved-state error invalidates the notebook cell as a certificate but does
not by itself refute the later displayed cubic coefficient.

The omitted `P1(L(h_a))=0` identity has now been independently repaired.  The
[symbolic audit](../../04-computation/brendle_hung_lemma54_independent_audit_20260824.py)
reconstructs the moving-frame geometry from the paper and proves on both
parametrized components of `Sigma intersect M_generic` that

```text
P1(L(h_a))=(0,0,0,0).                                    (2a)
```

This is a **VERIFIED EXACT SYMBOLIC IDENTITY relative to the reconstructed
formulas**, not a finite census.  The script parametrically simplifies all
four components on both generic pieces.  Equation (1h), independently
reconstructed in the global mixed-coefficient audit, proves that the same
normal Hessian is invertible everywhere on `M_generic`.  Consequently
`z(h_a)=-H^(-1)P1(L(h_a))=0` on both pieces is now also **VERIFIED EXACT
relative to the reconstructed formulas**.  Four exact point controls agree.
A deliberately non-special constant
symmetric tensor instead gives

```text
P1(L(h_hostile))=(0,0,-1/3,2/9)                           (2b)
```

at the recorded torus point, ruling out a vacuous all-zero implementation.
The [frozen output](../results/brendle_hung_lemma54_independent_audit_20260824.out)
LF-normalized byte-matches normal and optimized Python. This promotes only
the omitted
`P1(L(h_a))=0` clause to **VERIFIED / EXACT SYMBOLIC relative to the
reconstructed formulas**; it does not reach the nongeneric endpoints, replay
the other Lemma 5.4 clauses, audit the generic Hessian invertibility, or
recompute the cubic coefficient.

## 4. Audit matrix after the targeted repairs

| dependency | status | evidence / remaining boundary |
|---|---|---|
| headline metric with `sec>0` | **PREPRINT CLAIM / UNDER AUDIT** | consumes every row below |
| Section 2 minimization framework | **READ / NO IMMEDIATE CONTRADICTION** | no full formal replay |
| arXiv v1 TeX defects | **VERIFIED DEFECTS / locally repairable** | notation, metric, index, cross-term, label, and cutoff corrections above |
| `V_ac=V_bc=0` | **VERIFIED EXACT GLOBAL SYMBOLIC IDENTITIES** | direct reduction (1h)--(1k), symmetry (1a)--(1f'), and FINITE-EXACT point hostile (1g) |
| generic normal Hessian positive definiteness | **VERIFIED EXACT SYMBOLIC IDENTITY** | independent four-minor Sylvester factorization (1l) for every `x=tan(theta)>0`; no endpoint/uniform claim |
| Lemma 5.4 `P1(L(h_a))=0` clause | **VERIFIED EXACT SYMBOLIC IDENTITY** | relative to reconstructed formulas; hostile (2b) |
| Lemma 5.4 `z(h_a)=0` consequence | **VERIFIED EXACT SYMBOLIC IDENTITY** | combine (1h) with the exact zero response (2a) |
| remaining quadratic identities | **NOTEBOOK CLAIM / PARTIALLY AUDITED** | clean immutable replay still required |
| cubic coefficient `pi^2/(18sqrt(3))` | **LOAD-BEARING NOTEBOOK CLAIM / OPEN** | highest-value independent computation |
| Poisson repair from a nonzero mean | **PROVED ABSTRACTLY** | THM-3990; paper-specific mean remains claimed |

## 5. Exact remaining referee obligations

The shortest path to a trustworthy v1 certificate is:

1. **Fresh-kernel quadratic audit.** Rebuild the remaining eight `V^(2)`
   summands with immutable names and assemble all ten with the corrected
   `h_cd`, cross term, and cutoff.  The globally verified `V_ac,V_bc` pair and
   the exact generic point are positive and hostile controls for that replay.
2. **Complete Lemma 5.4 replay.** Retain the independently verified
   `P1(L(h_a))=0` identity and recompute `z(h_b),z(h_c)` from a cleared
   namespace on both components of `Sigma intersect M_generic`, followed by a
   separate boundary/extension audit on the connected completed `Sigma`.
3. **Rebuild the full lower bound.** Use the corrected `h_cd`, corrected
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

## 6. Reproduction of the targeted exact repairs

From the repository root:

```bash
python3 04-computation/brendle_hung_lemma54_independent_audit_20260824.py
python3 -O 04-computation/brendle_hung_lemma54_independent_audit_20260824.py
python3 04-computation/brendle_hung_vbc_exact_point_audit_20260824.py
python3 -O 04-computation/brendle_hung_vbc_exact_point_audit_20260824.py
python3 04-computation/brendle_hung_vac_vbc_global_identity_audit_20260824.py
python3 -O 04-computation/brendle_hung_vac_vbc_global_identity_audit_20260824.py
python3 04-computation/brendle_hung_normal_hessian_sylvester_independent_audit_20260824.py
python3 -O 04-computation/brendle_hung_normal_hessian_sylvester_independent_audit_20260824.py
```

The scripts use exact SymPy arithmetic.  They declare their universes,
positive controls, hostile controls, conclusions, and nonconsequences.  They
do not call or translate cached Mathematica output.

## 7. Transferable themes and new frontiers

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
