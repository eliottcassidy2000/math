---
id: THM-4251
title: "W=0 hidden degree-twenty-four attachment exclusion"
status: >
  PROVED RELATIVE TO THM-4230/4241/4247 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. The saturated W=0 hidden degree-24 Hom shell has 24 maps, is
  exactly twice the degree-6 shell, and has two free target-unit/source-T
  orbits. Each representative has reduced coefficient denominator
  t(t^2-1)L(t^2)^2 with reciprocal gcd t^2-1. Its only common roots t=+-1
  force Z/U=0 outside the gate. Thus no hidden degree-24 map collapses the
  twelve attachments. This independently reproves a row already contained
  in the stronger THM-4249 sieve; it does not further shrink the current
  residual or close W=0, M=12, or JC(2).
source: root/jc-planar/2026-08-26
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4247-w0-involution-degree-twelve-attachment-exclusion
related:
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
mistake_firewall:
  - MISTAKE-521
  - MISTAKE-522
  - MISTAKE-525
scripts:
  - 04-computation/jc23_w0_hidden_degree24_attachment_audit_thm4251.py
  - 04-computation/jc23_w0_hidden_degree24_attachment_independent_audit_thm4251.py
outputs:
  - 05-knowledge/results/jc23_w0_hidden_degree24_attachment_audit_thm4251.out
  - 05-knowledge/results/jc23_w0_hidden_degree24_attachment_independent_audit_thm4251.out
script_sha256:
  - 7adc5b697f7007ca19ec08766f60c378d9987e6203dd705adb4d8068a029bd50
  - b8c571d4eeef84e1ef5bb30ba5fbd1b8012fcac3c0ea5eca6f46b11b65e6e0ff
output_sha256:
  - 643e01337a5512157f05fa42f654293085326d03514978e3bc8924e8d505adbb
  - 0293b9aa8ceca98da25556aa68e3fc8e7f9ad454f294f658cf1b3f20f29b36d5
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The primary imports only THM-4247's maintained coefficient
  group law. The clean-room referee imports no repository computation and
  independently reconstructs both Grams, a widened shell census, generated
  symmetry orbits, characteristic-zero noncancellation norms, direct group
  laws in four good reductions, reciprocal gcds, the full glued projection
  histograms, and the overlap ledger. Normal and optimized outputs byte-match.
---

# THM-4251 -- hidden degree-twenty-four attachment exclusion

**PROVED RELATIVE TO THM-4230/4241/4247 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THM-4249 ALREADY SUBSUMES THE ROW CONSEQUENCE.**

## 1. Statement and inheritance

Use the `W=0` normalization

```text
C0: x^6+y^4=1,                 E0: Y^2=X^3+1,
O=Z[omega],                    omega^2+omega+1=0,
L=Hom(J(C0),E0)^(iota=-1),     iota:(x,y)->(x,-y).
```

THM-4241 proves that, for `g=Tf`,

```text
L=O f direct-sum O g,          T^2=-omega,
H_L=[6,-4-2omega; -4-2omega^2,6].                    (1)
```

> **Theorem.** The degree-`24` shell of `L` contains exactly `24` maps.
> It equals twice the degree-`6` shell and is the disjoint union of the two
> free size-twelve orbits represented by
>
> ```text
> [2]f,                         [2](f+Tf)              (2)
> ```
>
> under target `mu_6` and source `T`. No map in this shell sends all twelve
> gate-interior attachments to the origin of `E0`.

The closest mechanism is THM-4247's odd-denominator reciprocal test. The
hostile is THM-4241's index-four full-lattice glue; MISTAKE-521 forbids using
rational isogeny data as integral exhaustion. The new sidecar is the exact
duplication of the degree-six maps.

## 2. Complete shell

For `a,b in O`, the degree of `af+bg` has integer matrix

```text
[ 6,-3,-3, 0;
 -3, 6, 3,-3;
 -3, 3, 6,-3;
  0,-3,-3, 6].                                      (3)
```

The relative Hermitian eigenvalues are `6+-sqrt(12)`, so degree `24` gives
`N(a)+N(b)<12`. Since `N(m+n omega)>=(m^2+n^2)/2`, `[-4,4]^4` is complete.
Exact enumeration, independently widened to `[-6,6]^4`, gives

```text
#{q=6}=24,                  #{q=24}=24.               (4)
```

Every degree-`24` coordinate is even and halving bijects the two shells.
Generating the actions

```text
epsilon:(a,b)->(epsilon a,epsilon b),
T:(a,b)->(-omega b,a)                                (5)
```

produces exactly the two free orbits in `(2)`. The ordinary real matrix in
`(3)` has a different spectrum; only the relative Hermitian eigenvalue is
used in the completeness bound.

## 3. Exact doubled denominators

Put `u=t^2` and work over

```text
z^4-z^2+1=0,
p^2-(1+2z-z^3)p+1=0.                                (6)
```

Up to nonzero target-weight scalars, write a degree-six representative as

```text
X_r/x=A_r(u)/t,                 Y_r/y=L_r(u)/t.       (7)
```

For `r=f`, take `A_f=u-p^2`, `L_f=u+p^3`. For `r=f+Tf`, coefficient-form
addition has constant slope

```text
lambda=(-z^3p^3-1)/(z^2p^2-1),                       (8)
```

and again produces linear `A_r,L_r`. Duplication on `E0` gives, after
clearing nonzero field scalars,

```text
N_r(u)=9A_r(u)^4-8A_r(u)(u-1)L_r(u)^2,
d_raw(t)=8t(u-1)L_r(u)^2.                            (9)
```

The clean-room characteristic-zero path normalizes both representatives and
computes

```text
Norm(L(0))=Norm(lc(L))=1,
Norm(L(1))=2916,                 Norm(L(-1))=4,
Norm(A(1))=144,                  Norm(Res_u(A,L))=36,
Norm(N_r(0))=1.                                      (10)
```

Thus none of `t`, `u-1`, or `L(u)` cancels. The exact reduced denominator is

```text
d_r(t)=constant*t(t^2-1)L_r(t^2)^2,       deg d_r=7. (11)
```

The primary symbolic path agrees. Direct exact group-law reductions at

```text
(q,z,p,s)=(313,29,135,21),(349,24,246,28),
          (373,69,297,33),(397,157,161,27)            (12)
```

retain denominator degree seven in all eight representative/prime cases.

## 4. Reciprocal obstruction

Write `L(u)=b0+b1u`. Apart from `u-1`, its reciprocal line is `b1+b0u`, and

```text
Res_u(b0+b1u,b1+b0u)=b1^2-b0^2=-L(1)L(-1).           (13)
```

Equation `(10)` makes this nonzero. Hence

```text
gcd(d_r(t),t^7d_r(-1/t))=t^2-1.                      (14)
```

All four reductions in `(12)` independently reproduce the monic gcd. The
only finite common roots are `t=+-1`, but the marked ratio satisfies

```text
Z/U=((t^2-1)/(2t))^2.                                (15)
```

They force `Z=0`, outside `U*Z*(U+Z)!=0`. Target units fix the origin and `T`
permutes the attachments while sending `t` to `-1/t`; the orbit audit makes
the two cases complete. This proves the theorem.

## 5. Projection consequence and subsumption

THM-4247 proves that full-map collapse implies origin collapse of
`ell=m-m composed_with iota`. The independent full-gluing reconstruction gives

```text
q(m)=34, q(ell)=24: 2304 vectors,
q(m)=42, q(ell)=24:  288 vectors.                    (16)
```

They are disjoint from THM-4247's prior rows. Compared only with that earlier
ledger, the raw remainders would be `31,872/15,840`. Current canon is stronger:
THM-4249 independently excludes this same low-hidden shell, then leaves
`176/132` symmetry classes in `55/34` ratio envelopes. Equation `(16)` is an
independent hostile control, not a second subtraction from that frontier.

## 6. Scope and reproduction

Later hidden shells, the mixed evaluations in THM-4249's ratio envelopes,
the hidden-Hom locus away from `W=0`, exact `M=12`, entry, `JC(2)`, and `DC(2)`
remain open.

```bash
python3 -B 04-computation/jc23_w0_hidden_degree24_attachment_audit_thm4251.py
python3 -B -O 04-computation/jc23_w0_hidden_degree24_attachment_audit_thm4251.py
python3 -B 04-computation/jc23_w0_hidden_degree24_attachment_independent_audit_thm4251.py
python3 -B -O 04-computation/jc23_w0_hidden_degree24_attachment_independent_audit_thm4251.py
```

Both normal/optimized pairs byte-match the frozen outputs. **QED.**
