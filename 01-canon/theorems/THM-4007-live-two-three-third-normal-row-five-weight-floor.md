---
id: THM-4007
title: "The live reduced 2:3 third normal row forces five A-weights"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the oriented reduced pole-depth 2:3 cell, in the fixed THM-3992 target
  gauge and on the continued THM-3997 seam, the third source-normal Jacobian
  row has degree caps 4 and 5. On the only old minimal-support stratum b=d=0,
  its constant equation forces [x^0*t^3](A/a)=2176/(135*A5^3), hence a new
  nonzero A-weight -6. Together with THM-4005 this gives at least five
  retained A-weights on every live atlas stratum, rejects the fixed-gauge
  4x5 invoice, and moves the first unrejected fixed-gauge floor to 5x5. The
  A-floor transfers across exactly the THM-3992 scalings, translations, and
  C-by-A shear. The same row also forces
  [p^4]Rtilde+(6/(7*A5))*[y^2]Rtilde=-11392/(105*A5^4), while
  [p^2*y]Rtilde remains free at this order. No atlas-point existence,
  all-row lift, reversed orientation, other reduced cell, or JC(2) conclusion
  follows.
source: root + next_hasse_row / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (primary + independent no-import audit, 2026-08-24). The exact
  certificate reconstructs the third
  Laurent/source degree caps, solves the complete bounded t^2 Jacobian row,
  checks the load-bearing constant independently at x=0, derives the full
  invariant third rows and t^4 residual coupling, and passes a direct
  full-series fourth-jet hostile. It also records the distinct-row
  determinant and wrong-elimination-sign traps. Normal and optimized streams
  agree exactly across 54 gates. A SymPy-free verifier independently derives
  the axis coefficient, tests the residual row at three exact hostile atlas
  points, and detects the determinant-factor swap across 21 gates. Its normal
  and optimized streams match the frozen transcript and semantic hash.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4005-reduced-two-three-live-seam-invariant-support-atlas
related:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-3999-companion-divisor-boundary-endpoint-and-class-ledger
  - THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go
  - THM-4016-sharp-five-by-five-elliptic-attachment-nontorsion
script: 04-computation/jc2_live_23_third_normal_row_thm4007.py
output: 05-knowledge/results/jc2_live_23_third_normal_row_thm4007.out
script_sha256: 15f7b4d45175ddc69e8550748a919cc512abe9536c9b16b917827a6580cf0f1b
output_sha256: 1d31f05af4e536db702c73458e5034c710611fcbaa21bc017ad3583a98cf94b6
independent_audit_script: 04-computation/jc2_live_23_third_normal_row_thm4007_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_live_23_third_normal_row_thm4007_independent_audit.out
independent_audit_script_sha256: bf6ff9ca628c3caa22c57624d23e9e03129ce8d4fa5dbe83907156cfe7e54fa4
independent_audit_output_sha256: 488d42fadaa15c3dfe2023628c8ea923f77ba3bd049de30e45909c672e84fbe6
independent_audit_semantic_sha256: 148a234ce88f33ca899a1b0e6b234d19a9e866837cb1e3eb6968c69647bbfe27
hash_basis: raw LF bytes
---

# THM-4007 -- the third normal row forces five A-weights

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. Let
`(A,C) in B_2^2` be a hypothetical polynomial Darboux pair with

```text
J_(x,t)(A,C)=1.
```

Assume the oriented reduced pole-depth `(2,3)` cell, the specific target
normalization of THM-3992, and the live seam and invariant atlas of THM-3997:

```text
gamma=-a^3/2,       a!=0,       A5=a^5,
b=[y](R/gamma),     d=[py](R/gamma).                    (1)
```

THM-4005 left one minimal fixed-gauge stratum:

```text
b=d=0,
{2,0,-2,-4} subset supp(A),
{3,1,-1,-3} subset supp(C),                             (2)
```

with an additional even negative `C`-tail. The theorem below continues
exactly that stratum. It proves necessary identities; it does not assert that
any formal atlas point lifts to a pair in `B_2`.

## 1. Laurent support gives the complete third-diagonal caps

Use THM-3997's source-normal notation

```text
A=A_0+tN+t^2M+t^3U+O(t^4),
C=C_0+tK+t^2L+t^3V+O(t^4).                              (3)
```

For a finite Laurent polynomial `F=sum_i f_i(s) tau^i`, with `s=xtau`
and `tau=t`, the coefficient identity is

```text
[t^r]F(x,t)=sum_i [s^(r-i)]f_i(s) x^(r-i).              (4)
```

In the reduced cell, the extreme rows are exactly

```text
a_-2=gamma^2*s^2,             c_-3=gamma^3*s^3.         (5)
```

For `r=3`, the would-be `x^5` coefficient of `A` comes only from
`[s^5]a_-2`, which vanishes by `(5)`; every remaining Laurent index is at
least `-1` and hence contributes degree at most four. Similarly, the
would-be `x^6` coefficient of `C` is `[s^6]c_-3=0`, and every remaining
index is at least `-2` and contributes degree at most five. Therefore

```text
deg_x U<=4,                    deg_x V<=5.               (6)
```

The caps are sharp as support statements: `[s^4]a_-1` and `[s^5]c_-2`
can occupy their top slots. No converse from `(6)` to Laurent or `B_2`
membership is used.

## 2. The constant equation forces a new A-weight

On `b=d=0`, the first three THM-4005 rows, before division by their residual
characters, are

```text
A_0=a*(a^5*x^2+4)/4,
C_0=-a^4*x*(a^5*x^2+6)/8,

N=2*(3a^5*x^2+2)/(3a^4),
K=-x*(3a^5*x^2+8)/(2a),

M=-4*(9a^5*x^2+40)/(45a^9),
L=-4x*(9a^5*x^2-22)/(15a^6).                           (7)
```

Write

```text
U=sum_(j=0)^4 u_j*x^j,       V=sum_(j=0)^5 v_j*x^j.     (8)
```

The coefficient of `t^2` in the Jacobian is

```text
3(A_0'V-U C_0')
 +2(N'L-K'M)+(M'K-L'N)=0.                               (9)
```

At `x=0`, the exact lower-row values are

```text
A_0'=N'=K=M'=L=0,
C_0'=-3a^4/4,       N=4/(3a^4),
K'=-4/a,            M=-32/(9a^9),       L'=88/(15a^6). (10)
```

Thus the lower part of `(9)` is

```text
-2M(0)K'(0)-N(0)L'(0)=-544/(15a^10),                   (11)
```

and `(9)` gives

```text
(9a^4/4)u_0-544/(15a^10)=0,
u_0=2176/(135a^14).                                    (12)
```

After dividing `A` by its residual character `a`, this is the invariant
identity

```text
[x^0*t^3](A/a)=2176/(135A5^3) !=0.                     (13)
```

The monomial `x^j t^n` has source weight `j-2n`; hence `(13)` has weight
`-6`. It is distinct from all four weights in `(2)`, and a later diagonal
cannot cancel its coefficient at `t^3`. Consequently

```text
b=d=0  ==>  {2,0,-2,-4,-6} subset supp(A).              (14)
```

Combining `(14)` with the two nonminimal strata already proved in THM-4005
gives the repaired support atlas

```text
stratum             |supp(A)|     |supp(C)|
b!=0                    >=5           >=6
b=0,d!=0                >=5           >=6
b=d=0                   >=5           >=5.              (15)
```

The last `C` bound still uses THM-4005's imported even negative tail. Thus
the first support packet not rejected by the displayed rows in this fixed
gauge is `5x5`, not `4x5`. This is a necessary floor, not a construction.

## 3. The full third row and the first new residual coupling

Let

```text
e=[y^2]Rtilde,       f=[p^2*y]Rtilde,
g=[p^4]Rtilde,       Rtilde=R/gamma.                    (16)
```

Solving all coefficients of `(9)`, then comparing the coefficient of `t^4`
in the nodal defect with its `k[p,y]` expansion, gives the complete invariant
third rows

```text
[t^3](A/a)=
  2176/(135A5^3)
 +(A5*f/4)x
 +(1088/(315A5^2)+(2A5*e)/7)x^2
 -(32/(15A5))x^4,                                      (17)

[t^3](C/a^4)=
 -3f/8
 +(-8128/(315A5^3)-3e/7)x
 -(3A5*f/16)x^2
 +(736/(105A5^2)-(3A5*e)/14)x^3
 +(8/(5A5))x^5.                                        (18)
```

In particular, the first three missing residual coefficients are not all
independent:

```text
g+(6/(7A5))*e=-11392/(105A5^4).                        (19)
```

The source-odd scalar `f=[p^2y]Rtilde` remains free at this row. Equation
`(19)` implies that `e` and `g` cannot both vanish, but `(13)`, not that
weaker residual observation, is the support obstruction.

For completeness, before residual comparison the bounded solution of `(9)`
has free coordinates `v_2,v_3,v_4,v_5` and

```text
u_0=2176/(135a^14),
u_1=-4(a^5v_2-2v_4)/(3a^8),
u_2=-4(5a^6v_3-10av_5-32)/(15a^9),
u_3=-4v_4/(3a^3),              u_4=-4v_5/(3a^3),

v_0=2(a^5v_2-2v_4)/a^10,
v_1=2(45a^6v_3-90av_5-752)/(45a^11).                   (20)
```

The `t^4` comparison fixes

```text
v_2=-3a^9f/16,
v_3=736/(105a^6)-3a^9e/14,
v_4=0,                         v_5=8/(5a),              (21)
```

and its constant term is exactly `(19)`.

## 4. The determinant-sign hostile

The elimination used in deriving `(17)--(19)` has a sharp sign trap. Put

```text
q=gamma*x^2+3a/(2gamma),
grad P(A_0,C_0)=q*(-C_0',A_0').                         (22)
```

If `(R_4,S_4)` is the fourth source-normal vector and `J_low` is the part of
`[t^3]J` made from the first four displayed vectors, then

```text
4 det((A_0,C_0)',(R_4,S_4))+J_low=0.                   (23)
```

Therefore its contribution to `[t^4]P` is

```text
q det((A_0,C_0)',(R_4,S_4))=-q*J_low/4.                (24)
```

The minus sign in `(24)` cancels every power above `x^4`. Reversing it
creates a spurious degree-eight polynomial whose `x^8` coefficient is
`2a^4v_5`. Likewise, for distinct vectors `P,Q`, the lawful expression is

```text
det(P',Q)=P_A'Q_C-P_C'Q_A,                              (25)
```

not the same-row Wronskian `P_A'Q_C-P_AQ_C'`. Those two expressions agree
when `P=Q`, so a same-row test is a deliberately blind control.

The primary certificate solves `(23)` directly in the caps
`deg R_4<=5`, `deg S_4<=6` for the hostile specialization

```text
a=1,       v_2=v_3=v_4=0,       v_5=z.
```

Full expansion then gives `J=1+O(t^4)` and

```text
[t^4]P=(25x^4z-160x^2z-448x^2+60z+928)/15,             (26)
```

which has degree four for arbitrary `z`. Thus the apparent `x^8,x^7,x^6`
contradiction is a determinant/sign artifact, not a no-go for the entire
`b=d=0` seam.

## 5. Exact gauge transfer and scope

The target operations actually used in THM-3992 are nonzero diagonal
scalings, translations, and one triangular linear operation

```text
C_old=c*C_new+kappa*A_new+constant,       c!=0.         (27)
```

They do not mix anything into `A`; its nonconstant source-weight support is
only scaled. Hence the five-`A` floor in `(15)` transfers across exactly this
recorded normalization. In particular, the oriented live seam has no
pre-normalization packet with only four retained `A`-weights, including the
old `4x5` invoice.

The fixed-gauge `C>=5` floor need not transfer across `(27)`, because a copy
of an `A`-component in `C_new` can cancel. No invariance under arbitrary
nonlinear target automorphisms is asserted.

This theorem does **not** prove:

1. existence or all-row consistency of an `(A5,e,f)` atlas point;
2. that the displayed finite jets lift to elements of `B_2`;
3. anything about the reversed `(3,2)` orientation or another reduced cell;
4. a pre-gauge five-weight floor for `C`;
5. companion irreducibility, component ownership, or address completeness;
6. emptiness of the reduced `(2,3)` problem or `JC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_live_23_third_normal_row_thm4007.py
python3 -B -O 04-computation/jc2_live_23_third_normal_row_thm4007.py
python3 -B 04-computation/jc2_live_23_third_normal_row_thm4007_independent_audit.py
python3 -B -O 04-computation/jc2_live_23_third_normal_row_thm4007_independent_audit.py
```

The two primary streams agree exactly after `54` gates. The independent
streams match their frozen `21`-gate transcript and have semantic SHA-256
`148a234ce88f33ca899a1b0e6b234d19a9e866837cb1e3eb6968c69647bbfe27`.
