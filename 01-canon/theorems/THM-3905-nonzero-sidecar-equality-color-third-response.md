---
id: THM-3905
title: "Nonzero-sidecar equality colors have an exact third response"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every positive
  equal y-degree seam of the THM-3881 residual, the equianharmonic color law
  extends exactly through epsilon^3.  Seven residual monomials comprise the
  complete deficit-at-most-three support, with marked-source boundaries at
  n=1,2,3 and an extra n=1 KTf^2-KT^3 source.  The named address-compatible
  n=1 two-jet control from THM-3902 has no polynomial third response.  This is
  a necessary filtration, not equality-seam emptiness; JC(2) remains OPEN.
source: root / post-THM-3902 third-response audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  The 37-gate primary expands all sixteen
  residual terms, proves the support ledger for n=1,...,8, verifies the compact
  law symbolically and on hostile substitutions, denominator-clears the named
  n=1 control, and checks n=3/n>=4 controls.  A clean-room 39-gate path checks
  the nine excluded rows, the exact seven-row list, 28 raw coefficients for
  n=1,...,7, and the third polynomial sidecar.  A separate 22-gate quadratic-
  field audit proves the named n=1 denominator obstruction.  All normal,
  optimized, and frozen streams match.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3902-nonzero-sidecar-equality-color-two-jet-response
related:
  - THM-3899-nonzero-sidecar-y-degree-tariff-and-equianharmonic-equality-colors
  - THM-3904-nonzero-sidecar-constant-y-seam-emptiness-and-primitive-equality-colors
script: 04-computation/jc2_nonzero_sidecar_equality_color_third_response_thm3905.py
output: 05-knowledge/results/jc2_nonzero_sidecar_equality_color_third_response_thm3905.out
script_sha256: bbe12c8754d78f788f7dc0949bd2ea4519448f279f06e8312d58bdcb3cc353a6
output_sha256: 0f148b04bf73385320aaddadcd4804807edecd00e48383f3ed1a8e4c5f867f23
semantic_sha256: e47f4c9dedd0d95c39a119b9f85630b71965df4d915a2b96747e3311081a18ff
independent_audit_script: 04-computation/jc2_nonzero_sidecar_equality_color_third_response_independent_audit_thm3905.py
independent_audit_output: 05-knowledge/results/jc2_nonzero_sidecar_equality_color_third_response_independent_audit_thm3905.out
independent_audit_script_sha256: 063c67b6dfe6a51a474d78d867eab103f58e8bd71b03c9c1e9cda979d4422624
independent_audit_output_sha256: ae276e12f81a420b89a8856a8ddc3574b6f56bb473d8a06e641e36152b352aed
independent_audit_semantic_sha256: 15b2c85974162cd4ae967b3a74be5031a6f931f85159784a57b474a942c100da
hostile_audit_script: 04-computation/jc2_nonzero_sidecar_equality_color_third_response_n1_hostile_audit_thm3905.py
hostile_audit_output: 05-knowledge/results/jc2_nonzero_sidecar_equality_color_third_response_n1_hostile_audit_thm3905.out
hostile_audit_script_sha256: a9f93ee00852eec540d04463422576705c6251090a3f07bd2523520daf7449c0
hostile_audit_output_sha256: 40d982b574c4f22a68155924d43750559ca1ec17427e296267c563258b9caaac
hostile_audit_semantic_sha256: 4cd89863a111096eeda3fa5304939c87e1edefa3b4f113bb3a5bbadd4b2ca0b3
hash_basis: raw LF bytes
---

# THM-3905 -- the equality colors have an exact third response

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Retain the THM-3881
coordinates

```text
D=k[x,y],                  a=x+1,
L=9x+4,                    K=y^2+kappa,
kappa=-15x^2-15x-4,        P=aL^2,
r=aT+Kf,                   A=KT+aPf.                       (1)
```

Suppose

```text
f!=0,                      S(T,f)=G^2,
deg_y(T)=deg_y(f)=n>=1.                                   (2)
```

The response identity below does not use the THM-3881 source address, so it
applies in particular to the addressed residual lane.

## 1. Coordinates and colors at infinity

Put `epsilon=y^(-1)` and write

```text
T=y^n U,       U=u+u_1 epsilon+u_2 epsilon^2+u_3 epsilon^3+O(epsilon^4),
f=y^n V,       V=v+v_1 epsilon+v_2 epsilon^2+v_3 epsilon^3+O(epsilon^4),
G=y^(2n+2) Gamma,
Gamma=g+g_1 epsilon+g_2 epsilon^2+g_3 epsilon^3+O(epsilon^4).             (3)
```

Set `u_i=v_i=0` whenever `i>n`.  As in THM-3902, `v|g`; write `g=vh` and
choose `d in k` with `d^2=-3`.  Then

```text
h^2=3(aL^2v^2-u^2),
H=Gamma/V=h+j_1 epsilon+j_2 epsilon^2+j_3 epsilon^3+O(epsilon^4),
C_-=H-dU,                  C_+=H+dU.                       (4)
```

## 2. Complete deficit-three support

For a residual monomial `T^p f^q K^r`, its deficit from degree `4n+4` is

```text
Delta(p,q,r)=(4-p-q)n+4-2r.                               (5)
```

The complete union of rows with deficit at most three is

| residual term | `(p,q,r)` | deficit | active through this window |
|---|---:|---:|---|
| `3aL^2K^2f^4` | `(0,4,2)` | `0` | every `n` |
| `-3K^2T^2f^2` | `(2,2,2)` | `0` | every `n` |
| `2L^2K^2f^3` | `(0,3,2)` | `n` | `n=1,2,3` |
| `6a^2L^2KTf^3` | `(1,3,1)` | `2` | every `n` |
| `-6aKT^3f` | `(3,1,1)` | `2` | every `n` |
| `12aL^2KTf^2` | `(1,2,1)` | `n+2` | `n=1` |
| `-8KT^3` | `(3,0,1)` | `n+2` | `n=1` |

Each of the other nine residual monomials has nonnegative deficit slope and
deficit at least four already at `n=1`.  Thus no omitted row can enter the
third response.

## 3. Exact compact law modulo `epsilon^4`

Put

```text
E=aL^2V^2-U^2.                                             (6)
```

Then the color product is

```text
C_-C_+
 =3aL^2V^2
  +2L^2 epsilon^n(1+2kappa epsilon^2)V
  +6epsilon^2(kappa+aU/V)E
  +[n=1]epsilon^3(12aL^2U-8U^3/V^2)
  +O(epsilon^4).                                           (7)
```

Here `[n=1]` is one in degree one and zero otherwise.  The four boundary
forms are therefore

```text
n=1:   3aL^2V^2+2L^2 eps V+6eps^2 Q E
       +eps^3(4kappa L^2V+12aL^2U-8U^3/V^2),
n=2:   3aL^2V^2+2L^2 eps^2 V+6eps^2 Q E,
n=3:   3aL^2V^2+2L^2 eps^3 V+6eps^2 Q E,
n>=4:  3aL^2V^2+6eps^2 Q E,              Q=kappa+aU/V,     (8)
```

all modulo `epsilon^4`.  Thus the marked `K^2f^3` source arrives exactly at
response depth `n`; the third window sees it precisely for `n<=3`.

## 4. Denominator-cleared third response

Let

```text
E_0=aL^2v^2-u^2,
E_1=2aL^2vv_1-2uu_1,
rho_1=(u_1v-uv_1)/v^2,
delta_i=[n=i].                                             (9)
```

The coefficient of `epsilon^3` on the right of `(7)` is

```text
B_3=
 6aL^2(vv_3+v_1v_2)
 +2L^2(delta_1v_2+delta_2v_1+delta_3v)
 +4kappa L^2 delta_1v
 +6{(kappa+au/v)E_1+a rho_1 E_0}
 +delta_1(12aL^2u-8u^3/v^2).                            (10)
```

Retain THM-3902's polynomial sidecars `J_1=vj_1`, `J_2=v^2j_2` and
`R_(+/- ,i)=J_i +/- d v^i u_i` for `i=1,2`.  Define

```text
J_3=v^2(g_3-hv_3)-v_1J_2-vv_2J_1=v^3j_3,
R_(+/- ,3)=J_3 +/- d v^3u_3.                              (11)
```

Coefficient extraction from `(7)` gives the exact polynomial identity

```text
c_-R_(+,3)+c_+R_(-,3)
 +R_(-,1)R_(+,2)+R_(-,2)R_(+,1)=v^3B_3,                 (12)
```

where `c_-=h-du`, `c_+=h+du`.  Although `(10)` is written in the fraction
field, multiplication by `v^3` clears every displayed denominator in `(12)`.

## 5. Sharp controls and the killed named lift

Use THM-3902's canonical leading payment

```text
v=1,
h_*=(a+3L^2)/2,          u_*=(3L^2-a)/(2d),
h_*-du_*=a,              h_*+du_*=3L^2.                 (13)
```

With zero lower sidecars, the `n=3` marked source forces

```text
j_3=L^2/h_*,                                                (14)
```

which is not polynomial because `h_*` is quadratic and coprime to `L`.  For
every `n>=4`, the same zero-sidecar choice has `j_3=0`, giving an exact
positive third-jet control.  Neither assertion concerns a complete square.

More sharply, take the address-compatible `n=1` positive two-jet control in
THM-3902, Section 6.3.  Direct extraction of `[y^5]S` and the color law `(7)`
agree on its forced `g_3`.  In `k[x,d]/(d^2+3)` its reduced denominator is

```text
249143169618(243x^2+217x+49),                              (15)
```

and its numerator modulo the quadratic factor is

```text
(46118408/3)(x+1)(-1+4d),                                 (16)
```

which is nonzero.  Hence this named two-jet control has no polynomial third
response.  This kills one control, not the general `n=1` seam.

## 6. Scope

The theorem is a necessary third-jet law on `deg_y(T)=deg_y(f)>=1`.
THM-3904 closes only the x-only boundary.  Positive equality, all later
responses, the strict-degree lanes after THM-3901's tariffs, a Keller atlas,
and `JC(2)` remain **OPEN**.

The operation-depth analogy with the tournament response calculus is only a
research heuristic: no tournament theorem enters the proof.  Here the native
observable is the residual-monomial deficit `(5)`, and the marked source is a
literal polynomial term whose arrival depth is `n`.

## 7. Exact replay

Run from the repository root:

```bash
python3 -B 04-computation/jc2_nonzero_sidecar_equality_color_third_response_thm3905.py
python3 -B -O 04-computation/jc2_nonzero_sidecar_equality_color_third_response_thm3905.py
python3 -B 04-computation/jc2_nonzero_sidecar_equality_color_third_response_independent_audit_thm3905.py
python3 -B -O 04-computation/jc2_nonzero_sidecar_equality_color_third_response_independent_audit_thm3905.py
python3 -B 04-computation/jc2_nonzero_sidecar_equality_color_third_response_n1_hostile_audit_thm3905.py
python3 -B -O 04-computation/jc2_nonzero_sidecar_equality_color_third_response_n1_hostile_audit_thm3905.py
```

Each normal/optimized pair must byte-match its raw-LF frozen output.  The
finite degree samples audit coefficient extraction; the affine deficit proof
establishes the support cutoff for every `n>=1`.  **QED.**
