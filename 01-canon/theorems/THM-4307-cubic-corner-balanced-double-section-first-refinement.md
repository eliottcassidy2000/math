---
id: THM-4307
title: "Cubic-corner balanced double-section one-coordinate refinement extinction"
status: >
  PROVED RELATIVE TO THM-4301 AND THM-4304 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. The actual balanced repeated first face is an A1
  double root with one intrinsic deformation coordinate, not the provisional
  two-coordinate A2 model. Its discriminant is one smooth formal graph. The
  unique splitter beta=5s has a smooth rational normalized chart, and all
  completed-local exceptional divisorial refinement carriers above either
  balanced double section are rational and hence Keller-constant. This does
  not close the
  three Regime-A double sections: an exact k=3 hostile produces a genus-one
  j=1728 refinement, although that hostile is independently Keller-constant.
source: root / planar-Jacobian continuation session, 2026-09-01
depends_on:
  - THM-4301-cubic-corner-first-face-keller-extinction
  - THM-4304-cubic-corner-repeated-first-face-rationality
related:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4299-d-zero-square-face-elliptic-splitting-and-off-corner-extinction
primary_script: 04-computation/jc23_cubic_corner_balanced_double_section_first_refinement_thm4307.py
primary_output: 05-knowledge/results/jc23_cubic_corner_balanced_double_section_first_refinement_thm4307.out
primary_script_sha256: f3c5111c6f2cc5af8a2db5427795ad45a7c2a2ddf0c21fb8efcf24c2332a59a9
primary_output_sha256: 9e64179fc25929e3d0f0fa9ab3f0b4ea7be7a358df6f20f33154cf72ae613228
independent_audit_script: 04-computation/jc23_cubic_corner_balanced_double_section_first_refinement_independent_audit_thm4307.py
independent_audit_output: 05-knowledge/results/jc23_cubic_corner_balanced_double_section_first_refinement_independent_audit_thm4307.out
independent_audit_script_sha256: aa9a8c1673e227b5dd8b27ed09f2665d7aecb25842a219edc867f55bdecf8034
independent_audit_output_sha256: bfddfd95d34e1d7862fd57cb5898fbe427e5a4076d98b7c35c8d642e4e3c6e16
hash_basis: raw LF bytes
audit: >
  PASS. A SymPy reconstruction and a dependency-free arithmetic path agree on
  the complete strict transform, local A1 algebra, both quadratic-field
  coefficient values, smooth discriminant graphs, splitter at beta=5s,
  rational balanced normalization, and the Regime-A elliptic hostile. Normal
  and optimized runs have byte-identical output on both paths.
---

# THM-4307 -- Cubic-corner balanced double-section one-coordinate refinement extinction

**PROVED RELATIVE TO THM-4301 AND THM-4304 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE BALANCED COMPLETED-LOCAL REFINEMENT TOWER IS
KELLER-CONSTANT. THE THREE REGIME-A DOUBLE SECTIONS, `U=0`, `Z=0`, SEAM
ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Use THM-4304's exact-`M=12` cubic corner

```text
D=Lambda=0,                 UZ!=0,                 beta>2s,           (1)
t=sigma*z,                  q=t^4 y,
Delta=2048/45,              Theta=-Delta,
alpha_11=beta_11=upsilon_5=xi_10=eta=zeta_3=Phi=0.                  (2)
```

The balanced first-face cubic has a double root precisely at the two values

```text
U_epsilon=-315392/3645
          +epsilon*(217088 sqrt(265))/18225,
rho_epsilon=-(15+epsilon*3sqrt(265))/256,   epsilon in {+1,-1}.      (3)
```

> **Theorem.** In the completed local ring at either section
> `y=rho_epsilon`, the repeated locus is one smooth formal graph. Its only
> valuation splitter is `beta=5s`. The normalized tie chart is smooth with
> rational exceptional residue, and every later exceptional divisorial
> refinement carrier centered at its closed point is rational. Consequently
> every such balanced refinement carrier maps constantly to the good elliptic
> target.

The assertion is formal-local above the two balanced sections. It is not a
seam-entry theorem and does not classify the three Regime-A sections.

The inheritance pass and live board were:

- closest proved mechanism: THM-4304 reduces every repeated first face to an
  exact double section with rational reduced carrier;
- canonical hostile: THM-4301's separable genus-four cubic shows that positive
  genus alone is not excluded by the source;
- corrected near miss: the two-dimensional algebra retained in THM-4304
  belongs to a triple root, which THM-4304 itself excludes;
- least-used sidecar: the critical-root graph, rather than the raw pair of
  translated Taylor coefficients;
- live concepts: double section, critical recentering, discriminant contact,
  normalized tie chart, good-form order, and seam-entry firewall.

## 2. Exact balanced strict transform

Put

```text
x=t^2,                    m=z^12/t^8.                     (4)
```

Direct substitution in THM-4301's literal source gives

```text
F/t^12=H_U(y,x)-m y,                                       (5)
H_U=p_U+xR_2+x^2R_4+x^3R_6+x^4R_8+x^5R_10+x^6R_12+x^8R_16,
p_U=Uy^3+(2048/45)y^2+(8/3)y-1/2,                         (6)
R_2 =-(1376/135)y^2-3y,
R_4 =4Uy^4+(2048/15)y^3+(16/3)y^2,
R_6 =-(2752/135)y^3-3y^2,
R_8 =6Uy^5+(2048/15)y^4+(8/3)y^3,
R_10=-(1376/135)y^4,
R_12=4Uy^6+(2048/45)y^5,
R_16=Uy^7.                                                 (7)
```

The exact scripts reconstruct `(5)--(7)` from the unexpanded source rather
than taking these rows as input.

## 3. The intrinsic coordinate is one-dimensional

At `U=U_epsilon`, the root `rho_epsilon` is double and the third root is
distinct. After translating `Q=y-rho_epsilon`, the simple factor is a unit,
so the local hypersurface singularity has leading term `Q^2`. Hence

```text
T^1=C[[Q]]/(partial_Q Q^2)=C[[Q]]/(Q)=span{1}.             (8)
```

The formerly displayed quotient `C[[Q]]/(Q^2)=span{1,Q}` is the algebra of
the excluded triple root `Q^3`, not of the actual double section. The raw
constant and linear Taylor rows remain legitimate ambient data, but critical
recentering combines them into one critical-value coordinate.

For `P(y,x,m)=H_U(y,x)-my`, the common-root equations are

```text
E:=H_U-y partial_y H_U=0,              m=partial_y H_U.    (9)
```

At `(y,x)=(rho_epsilon,0)`, define

```text
a_epsilon=(1/2)p_U''(rho_epsilon)
          =-6784/135+epsilon*(128sqrt(265))/135,
g_epsilon=R_2(rho_epsilon)
          =(-707+epsilon*65sqrt(265))/3072.               (10)
```

All of `rho_epsilon`, `a_epsilon`, and `g_epsilon` are nonzero, and

```text
partial_y E(rho_epsilon,0)=-2rho_epsilon a_epsilon!=0.    (11)
```

Thus `(9)` defines a unique smooth formal discriminant graph

```text
m_epsilon(x)=chi_epsilon x+lambda_epsilon x^2+O(x^3),     (12)
chi_epsilon=-173/72+epsilon*(43sqrt(265))/360,
lambda_epsilon=11975/27648
               -epsilon*(53621sqrt(265))/7326720.
```

Here `chi_epsilon=g_epsilon/rho_epsilon` and
`Norm(chi_epsilon)=269/135!=0`. Equivalently, the critical value begins
`g_epsilon x-rho_epsilon m`, and the local discriminant is a unit times
`m-m_epsilon(x)`.

## 4. Splitter, normalized chart, and extinction

For a divisorial valuation with `v(sigma)=s` and `v(z)=beta`,

```text
v(m)=4(beta-2s),             v(x)=2(s+beta),
v(m)-v(x)=2(beta-5s).                                  (13)
```

There are exactly three zones:

```text
2s<beta<5s:     m is first and is a square monomial;
beta=5s:        the two terms tie;
beta>5s:        x=t^2 is first and is a square monomial. (14)
```

On the tie put `z=sigma^5 w`. Then

```text
x=sigma^12 w^2,              m/x=w^2.                    (15)
```

Write `psi_epsilon(x)=m_epsilon(x)/x`. The one-variable formal Morse lemma
and a square root of its unit coefficient put the normalized double cover in
the form

```text
V^2=w^2-psi_epsilon(sigma^12 w^2).                       (16)
```

Its exceptional residue is

```text
V^2=w^2-chi_epsilon,                                     (17)
```

a smooth projective conic. The total chart `(16)` is smooth: where `V` is a
unit use `partial/partial V`; where `V=0`, equation `(16)` makes `w` a unit
and `partial/partial w` is a unit. Its two discriminant arcs are

```text
z=sigma^5[+/-sqrt(chi_epsilon)
          *(1+(lambda_epsilon/2)sigma^12+O(sigma^24))].  (18)
```

In the two strict zones the leading radicand is already a square monomial;
on the tie the only new residue is the rational conic `(17)`. Since `(16)`
is a smooth surface germ, every exceptional divisor obtained by refining its
closed point comes from a finite sequence of point blowups, whose exceptional
curves are projective lines. This proves the asserted completed-local
rationality for exceptional refinement carriers. It makes no claim that an
arbitrary height-one strict transform on the smooth chart is rational.

There is also a differential check away from higher contact with `(12)`. If
`d=v(m-m_epsilon(x))`, then

```text
v(F_q)=8(s+beta)+d/2,
v(phi^*eta_0)>=s+3beta-d/2.                              (19)
```

In the two strict zones the resulting lower bounds are respectively
`5s+beta` and `2beta`, both positive. On the discriminant arcs the
smooth-chart rationality argument, not an unjustified generic value of `d`,
supplies constancy.

## 5. Exact Regime-A hostile

Balanced rationality is not a universal property of double-section
refinements. Set

```text
U=135/28672,                   eta=1,        zeta_3=-1,
alpha_11=beta_11=upsilon_5=xi_10=Delta=Theta=Phi=0.       (20)
```

The `k=3` Regime-A face at

```text
(s,beta,gamma)=(1,2,9),          rho=-14336/135           (21)
```

is a literal double section. Its discriminant graph begins
`m_3^*=(8/3)t^2+...`; on the tie chart `w=z/sigma^2`, the first normalized
exceptional curve is

```text
V^2=w^4-8/3.                                             (22)
```

The quartic is squarefree, so `(22)` has genus one. Its binary-quartic
invariants are `I=-32`, `J=0`, hence `j=1728`. The good elliptic target in
THM-4230 has `j=0`; the two CM endomorphism algebras differ, so they are not
isogenous and every map between them is constant. Independently, the literal
ledger gives

```text
v(F_q)=21,                    v(phi^*eta_0)=31-21=10>0.   (23)
```

Thus positive genus can genuinely be born above a cubic-corner double
section, but this exact hostile still does not carry a nonconstant Keller map.
The remaining Regime-A problem is to classify all three contact graphs and
their possible higher refinements, not to assume rationality.

## 6. Reproduction and scope

Run

```bash
python3 -B 04-computation/jc23_cubic_corner_balanced_double_section_first_refinement_thm4307.py
python3 -B -O 04-computation/jc23_cubic_corner_balanced_double_section_first_refinement_thm4307.py
python3 -B 04-computation/jc23_cubic_corner_balanced_double_section_first_refinement_independent_audit_thm4307.py
python3 -B -O 04-computation/jc23_cubic_corner_balanced_double_section_first_refinement_independent_audit_thm4307.py
```

The universe is the literal exact-`M=12`, `D=Lambda=0`, `UZ!=0` source at
the two balanced repeated coefficient values `(3)`. The theorem neither
crosses `U=0` or `Z=0`, classifies the three Regime-A double sections, proves
seam entry, nor proves `JC(2)` or `DC(2)`.

**QED.**
