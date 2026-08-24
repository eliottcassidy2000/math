---
id: THM-4016
title: "Sharp five-by-five formal weight-six point non-torsion"
status: >
  PROVED ARITHMETIC + VERIFIED-EXACT + TWICE INDEPENDENTLY AUDITED; FORMAL
  TRUNCATION POINT; DIRECT SHARP-FACE APPLICATION REFUTED. The weight-six
  formal truncation of the exact THM-4007 sharp 5x5 coefficients gives six
  normalized points on E:Y^2=X^3+1, represented by
  alpha^3=43/84 and beta^2=127/84. The representative reduces to points of
  exact orders 12 and 6 at good places above 11 and 17, so it is non-torsion;
  all six automorphic images are non-torsion. THM-4017 proves that the same
  sharp data force a nonzero p^4 term with rho-exponent -2 on the proposed
  weight-six scale. Hence these formal points are not proved attachment
  points of the sharp stable reduction, and the old conditional face
  exclusion is REFUTED. The full reduced 2:3 cell and JC(2) remain OPEN.
source: root + two no-import audits / reduced 2:3 cell continuation, 2026-08-24
audit: >
  PASS. The primary verifier reconstructs the sharp residual, normalized
  point, degree-six coordinate field, two good places, exact group law, and
  reduction-kernel contradiction. A separate search enumerates 38 degree-one
  good reductions below 500 and selects 11/17 as the first incompatible pair
  with unique cube roots. A second no-import referee checks every constant,
  the cyclic-kernel quotient, ramification/residue-extension boundary,
  six-point orbit, Eisenstein unit, and arithmetic scope boundary. All six
  normal/optimized streams match their frozen outputs. Two later independent
  Newton audits preserve every arithmetic gate but refute the geometric
  identification at the forced p^4 coefficient.
depends_on:
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go
related:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4017-sharp-weight-eight-specialization-obstruction-and-newton-ledger
script: 04-computation/jc2_sharp_attachment_nontorsion_thm4016.py
output: 05-knowledge/results/jc2_sharp_attachment_nontorsion_thm4016.out
script_sha256: 25540456322fa01541efd6fa5657927f01f14c9b3a21b26e00b4a96f32898240
output_sha256: 224919176779dbe19c226466d01104e57c05337baaf64d3b4edf229e5b4531be
semantic_sha256: 1484315bd0c7f31e9f27e1a0cdbe88fba0541af06e53bcb20c22fe57cb042574
search_script: 04-computation/jc2_sharp_attachment_reduction_search_thm4016.py
search_output: 05-knowledge/results/jc2_sharp_attachment_reduction_search_thm4016.out
search_script_sha256: 5f23e9e1aa963a94ce31a77db1b387f2d2053f3f78841acc65aec369b227fdc0
search_output_sha256: 0123bd1344d890187ec35aa4874594324a308509850c2351f9429ad683ba5b49
independent_audit_script: 04-computation/jc2_sharp_attachment_nontorsion_thm4016_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_sharp_attachment_nontorsion_thm4016_independent_audit.out
independent_audit_script_sha256: cc23b5b7bf8588b6b41de58ae7b3648215c77d2c46edfd3cd4e0879f20860513
independent_audit_output_sha256: a9f4303e6563310f24eca0dd7524953f122f6793df222e9b9a69651fe4c5bfc0
independent_audit_semantic_sha256: 6d788e109e52f7a3c191ea755892dc19becd5c89d3ba6fe28c9c6e677a4e4d5d
hash_basis: raw LF bytes
---

# THM-4016 -- the sharp formal weight-six point is non-torsion

**PROVED ARITHMETIC + VERIFIED-EXACT + TWICE INDEPENDENTLY AUDITED; FORMAL
TRUNCATION POINT; DIRECT SHARP-FACE APPLICATION REFUTED.** Work over
characteristic zero in THM-4007's sharp formal `5x5` hostile. The theorem
proves an unconditional arithmetic fact about the point obtained by retaining
the weight-six coefficients. THM-4017 shows that this truncation is not the
leading stable model of the exact sharp residual.

## 1. Exact sharp coefficients

Write `A5=a^5` and `Rtilde=R/gamma`. THM-4007 gives

```text
[p^3]Rtilde=2752/(135 A5^3),
105 A5^4 c40+90 A5^3 c02+11392=0,                    (1)
```

where `c40=[p^4]Rtilde` and `c02=[y^2]Rtilde`. On the sharp `5x5`
face, `c21=[p^2y]Rtilde=0` and cancellation of the new source-normal weight
gives

```text
A5 c40/2+256/(9 A5^3)=0.
```

Therefore

```text
c40=-512/(9 A5^4),             c02=-8128/(135 A5^3). (2)
```

Since `gamma=-a^3/2`, the raw weight-six coefficients are

```text
epsilon=[p^3]R=-1376/(135 a^12),
kappa=[y^2]R=4064/(135 a^12),
epsilon+kappa=2688/(135 a^12)!=0.                     (3)
```

In the formal truncation that discards the simultaneously forced `p^4` term,
put `P_source=S^2`, `T=S^3`, and `(epsilon+kappa)S^6=1`. After normalizing
its elliptic component to `E:Y^2=X^3+1`, the six formal points satisfy

```text
X^3=-epsilon/(epsilon+kappa)=43/84,
Y^2=kappa/(epsilon+kappa)=127/84.                     (4)
```

Choose one representative

```text
P=(alpha,beta),              alpha^3=43/84,
                             beta^2=127/84.            (5)
```

## 2. Two explicit good places

Put `K=Q(alpha,beta)`. The `43`-valuation shows that
`T^3-43/84` is irreducible, and the `127`-valuation does the same for
`T^2-127/84`. The degree-three and degree-two fields have trivial
intersection, so `[K:Q]=6`. In `R=Z[1/84,alpha,beta]`, the evaluations

```text
mod 11: alpha -> 9,       beta -> 2,
mod 17: alpha -> 2,       beta -> 3                  (6)
```

respect both defining equations and define maximal ideals; choose places of
`K` above them in the normalization. The reduced coordinates already lie
in the prime fields, and residue-field extension does not change a fixed
point's order. Since

```text
Delta(E)=-432=-2^4*3^3,                                (7)
```

both are good-reduction places and `P` is integral there.

## 3. Exact reduced orders

At the place above `11`, `P_11=(9,2)` and direct group-law arithmetic gives

```text
2P_11=(2,8),     3P_11=(5,4),
4P_11=(0,10),    6P_11=(10,0).                         (8)
```

Thus `6P_11` is nonzero two-torsion and `4P_11` is nonzero. Every proper
divisor of `12` divides `4` or `6`, so

```text
ord(P_11)=12.                                          (9)
```

At the place above `17`, `P_17=(2,3)` and

```text
2P_17=(0,1),              3P_17=(16,0).               (10)
```

The second point is nonzero two-torsion and the first is nonzero, whence

```text
ord(P_17)=6.                                           (11)
```

The independent search checks all `38` degree-one good reductions below
`500` and finds `(9)--(11)` as its first incompatible unique-cube-root
pair.

## 4. All-order non-torsion

At a good-reduction place of residue characteristic `p`, the reduction
kernel has only `p`-primary torsion: for every prime `ell!=p`, the
formal-group series `[ell](T)=ell T+O(T^2)` has unit linear coefficient.

Suppose `P` had finite order `N`. On the cyclic group generated by `P`,
the kernel of reduction has cardinality `N/ord(Pbar)`, so `(9)` forces

```text
N=12*11^r                                             (12)
```

for some `r>=0`, while `(11)` forces

```text
N=6*17^s                                              (13)
```

for some `s>=0`. These give simultaneously `v_2(N)=2` and `v_2(N)=1`,
a contradiction. Therefore

```text
P is non-torsion in E(Qbar).                           (14)
```

This is an all-order argument; no bounded torsion classification or finite
division-polynomial ceiling is used.

## 5. All six formal points

Let `zeta_3` be a primitive cube root of unity. The curve automorphism

```text
sigma(X,Y)=(zeta_3 X,-Y)                              (15)
```

fixes the origin, and its six iterates on `P` exhaust the three cube-root
and two square-root choices in `(4)`. As a group automorphism it preserves
torsion, so all six attachments are non-torsion.

Equivalently, the corresponding Eisenstein unit is a primitive sixth root
`u`; `u-1` is again an Eisenstein unit. This is the exact unit used in
the conditional invoice below.

## 6. The former sharp-face join is refuted

The abstract geometric implication remains sound with the owners stated in
THM-4017: if an elliptic component owns positive degree and all six branches
meet one connected contracted rational clutch, their images coincide. After
translation the elliptic restriction is an isogeny `psi`, so equality of
adjacent orbit points gives

```text
psi((sigma-1)P)=0.                                      (16)
```

Since `sigma-1` is an Eisenstein unit, `(16)` would make `P` torsion.

But the exact sharp residual does not have the assumed weight-six model.
Equation `(2)` forces `[p^4]R=256/(9a^17)!=0`. Under the weight-six scaling
this term has rho-exponent `-2`, so the family is nonintegral and clearing
the denominator changes the special fibre. Consequently, the points in
`(4)` are not proved attachment points of the actual sharp stable reduction.
The old direct sharp-face application is **REFUTED**, not merely conditional.

## 7. Scope firewall and reproduction

This theorem does **not** prove:

1. that the formal THM-4007 hostile is an atlas point or a `B_2` pair;
2. that its formal weight-six points are actual stable attachments;
3. any exclusion of the formal `5x5` survivor;
4. exclusion of other values of `[y^2]R`, resonant or higher weighted faces,
   or faces with additional leading monomials;
5. emptiness of the reduced `(2,3)` cell or `JC(2)`.

Reproduction:

```bash
python3 -B 04-computation/jc2_sharp_attachment_nontorsion_thm4016.py
python3 -B -O 04-computation/jc2_sharp_attachment_nontorsion_thm4016.py
python3 -B 04-computation/jc2_sharp_attachment_reduction_search_thm4016.py
python3 -B -O 04-computation/jc2_sharp_attachment_reduction_search_thm4016.py
python3 -B 04-computation/jc2_sharp_attachment_nontorsion_thm4016_independent_audit.py
python3 -B -O 04-computation/jc2_sharp_attachment_nontorsion_thm4016_independent_audit.py
```

All streams match their frozen outputs. The remaining honest target is the
complete residual and its actual highest-face/stable reduction, including the
forced weight-eight support recorded in THM-4017. **QED.**
