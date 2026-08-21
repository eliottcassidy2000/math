---
id: THM-3571
title: "Quadratic target-graph Euler no-go"
status: >
  PROVED + VERIFIED-EXACT.  Every collision-compatible genuinely quadratic
  target graph c+phi(a,b)=0 for the fixed THM-1300 Keller map has irreducible
  complete pullback X with compactly supported Euler characteristic at least
  two.  Hence X is not A2 and the pullback polynomial has no source-coordinate
  factor.  The proof is an exact discriminant-support classification: the
  strict Jelonek section has chi=2-rho(delta)-nu, where rho(delta)>=3 and
  nu>=1, while the omitted support has at most five points.  All sixteen
  possible two-root multiplicity patterns are eliminated over QQ; the only
  two nonempty ideals force the quadratic part of phi to vanish.
source: kps-s188 + subagent-zeno
depends_on:
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
  - THM-3564-target-graph-trisection-degree-resonance-irreducibility
  - THM-3565-resonant-linear-a-target-graph-factor-classification
related:
  - THM-3570-universal-pell-conic-target-graph-factor-compiler
companion: 04-computation/jacobian_quadratic_target_graph_euler_no_go_kps_s188.py
output: 05-knowledge/results/jacobian_quadratic_target_graph_euler_no_go_kps_s188.out
script_sha256: feb93461bb523404b8720385758ecb66df22a511c68bf5cd2fbd7a52ee5fe129
output_sha256: 895e98e2c9a427ad7416121d0a3bddb8c778e096700d3a4cbc5ca3d795c75b9d
hash_basis: LF-normalized bytes
---

# THM-3571 -- quadratic target-graph Euler no-go

**PROVED + VERIFIED-EXACT.**  The complete collision-compatible quadratic
cell of triangular target coordinates is closed: none pulls back to a source
coordinate for the fixed THM-1300 map.

All varieties are over `C`, and `chi` denotes compactly supported topological
Euler characteristic.

## 1. Statement and collision cell

Every polynomial of total degree at most two whose graph contains the unique
collision value `(-1/4,0,0)` has the form

```text
phi(a,b)=A a^2+B ab+C b^2+D a+E b+D/4-A/16.             (1)
```

Call it genuinely quadratic when

```text
(A,B,C)!=(0,0,0).                                      (2)
```

Let

```text
T_phi=V(c+phi(a,b)) ~= A2,
X_phi=F^(-1)(T_phi)=V(F3+phi(F1,F2)).                   (3)
```

Then

```text
chi(X_phi)>=2.                                         (4)
```

Moreover, `X_phi` is irreducible.  Therefore it is not `A2`, its defining
polynomial is not a source coordinate, and it has no source-coordinate
factor.

## 2. A polynomial strict transform of the Jelonek section

THM-3560 writes the Jelonek hypersurface in target coordinates as

```text
L=27a^2c^2-18abc+b^3c+16a-b^2.                         (5)
```

Resolve its quadratic discriminant by

```text
a=s(2b-s)/12.                                          (6)
```

After `(6)`, equation `(5)` factors as

```text
48L=(4b+3cs^2-8s)
     (12b^2c-12bcs-12b+3cs^2+8s).                     (7)
```

The two factors are exchanged by `s -> 2b-s`.  On the first sheet,

```text
c=4(-b+2s)/(3s^2).                                     (8)
```

Substitute `(1)` and clear the harmless sheet denominator.  The strict graph
section is

```text
G(s,b)=q2(s)b^2+q1(s)b+q0(s)=0,                        (9)

q2=s^2(As^2+6Bs+36C)/48,

q1=-(As^5+3Bs^4-6Ds^3-36Es^2+48)/48,

q0=s(As^5-9As-12Ds^3+36Ds+384)/192.                  (10)
```

The finite map

```text
(s,b) |-> (s(2b-s)/12,b)                               (11)
```

restricts to a constructible bijection from `V(G)` to the target Jelonek
section

```text
D_phi=(T_phi intersect V(L))_red.                      (12)
```

Indeed, the two roots of `s^2-2bs+12a` correspond to the two factors in
`(7)`.  Exactly one factor vanishes unless the roots coincide, in which case
they give one point.  On the apparent denominator boundary `a=0`, the root
`s=0` is excluded when `b!=0` because `G(0,b)=-b`; at `(0,0)` it gives one
point and `partial G/partial b=-1`.  Thus no affine point is inserted or lost.
Consequently

```text
chi(D_phi)=chi(V(G)).                                  (13)
```

## 3. The vertical-fibre-corrected Euler formula

Let

```text
Delta=q1^2-4q2q0=delta_tilde/256,                      (14)
```

where

```text
delta_tilde =
 256-384E s^2-64(24C+D)s^3
 +(36AC-224B-144CD+144E^2)s^4
 +(6AB-32A-24BD+48DE)s^5
 +(A^2-4AD-24BE+48CD+4D^2)s^6
 +4(BD-2AE)s^7+(B^2-4AC)s^8.                         (15)
```

For a nonzero one-variable polynomial `P`, write `rho(P)` for the number of
its distinct complex roots.  Define

```text
nu = #{alpha:q2(alpha)=0 and
              (q1(alpha),q0(alpha))!=(0,0)}.           (16)
```

Euler integration of the projection `V(G)->A1_s` gives the exact identity

```text
chi(D_phi)=2-rho(delta_tilde)-nu.                       (17)
```

The correction `nu` cannot be dropped.  The complete fibre table is:

```text
q2!=0, Delta!=0:                 2 points;
q2!=0, Delta=0:                  1 point;
q2=0, q1!=0:                     1 point;
q2=q1=0, q0!=0:                  0 points;
q2=q1=q0=0:                      one whole A1, chi=1.
```

Relative to the generic contribution two, a discriminant root costs one;
exactly the non-entire vertical fibres in `(16)` cost one more.  This proves
`(17)` even when the three coefficients have a common root.

There is always a vertical correction:

```text
q2(0)=0,                 q1(0)=-1,                     (18)
```

so

```text
nu>=1.                                                    (19)
```

## 4. At least three discriminant roots

The discriminant polynomial has

```text
delta_tilde(0)=256,       delta_tilde'(0)=0,            (20)
```

and degree at most eight.  We prove

```text
rho(delta_tilde)>=3                                      (21)
```

under the genuinely quadratic hypothesis `(2)`.

First suppose `delta_tilde=256` is constant.  Its coefficients successively
give

```text
E=0,
D=-24C,
B=9C(A+96C)/56.                                        (22)
```

After these substitutions, the remaining exact coefficient ideal in
`Q[A,C]` has Groebner basis

```text
[A,C^2].                                                (23)
```

Thus every complex point has `A=C=0`, and `(22)` then gives `B=0`, contrary
to `(2)`.

A nonconstant polynomial satisfying `(20)` cannot have one distinct root:
after normalizing its constant term it would be `256(1+ps)^m` with `p!=0`,
whose derivative at zero is nonzero.

It remains to exclude exactly two roots.  Order their positive
multiplicities `m<=n`.  Since the degree is at most eight, there are exactly
sixteen possibilities, and `(20)` forces the normal form

```text
delta_tilde=256(1+ps)^m(1-(m/n)ps)^n,   p!=0.           (24)
```

For every pair, equate coefficients through `s^8`.  The coefficients of
`s^2,s^3,s^4` solve uniquely for `E,D,B`; saturate the remaining ideal by
`p!=0` using `zp-1`.  Exact Groebner reduction over `Q` gives

```text
m=1, n=1,...,7:             unit ideal;
m=2, n=3,5,6:               unit ideal;
m=3, n=3,4,5:               unit ideal;
m=4, n=4:                   unit ideal.                 (25)
```

The only nonunit cells are `(2,2)` and `(2,4)`.  For `(2,2)`,

```text
E=4p^2/3,  D=-24C,  B=9C(A+96C)/56,
GB=[pz-1,A,C].                                          (26)
```

For `(2,4)`,

```text
E=p^2,  D=-24C-2p^3,
B=9C(A+96C+8p^3)/56,

GB=[-Cp+4Cz,pz-1,7A+144C,C^2,Cp^2-4C].                (27)
```

At every complex point, `(26)` or `(27)` forces `A=B=C=0`.  Both survivors
are therefore affine boundary rows, not genuinely quadratic rows.  This
exhausts the one-, two-, and zero-root alternatives and proves `(21)`.

## 5. Omitted support and the final inequality

Parameterize the omitted target curve using its nonzero `b` coordinate:

```text
a=b^2/12,                   c=4/(3b).                   (28)
```

On `(28)`, the graph equation is

```text
c+phi(a,b)=f_tilde(b)/(144b),                           (29)

f_tilde=
 A b^5+12B b^4+(144C+12D)b^3+144E b^2
 +(-9A+36D)b+192.                                      (30)
```

Because `f_tilde(0)=192` and `deg(f_tilde)<=5`, its reduced omitted support
`e_phi=T_phi intersect E` satisfies

```text
chi(e_phi)=rho(f_tilde)<=5.                             (31)
```

Apply THM-3560's exact `3/1/0` fibre Euler formula to `T_phi~=A2`:

```text
chi(X_phi)
 =3-2chi(D_phi)-chi(e_phi)
 =-1+2rho(delta_tilde)+2nu-rho(f_tilde).               (32)
```

Equations `(19)`, `(21)`, and `(31)` now give

```text
chi(X_phi)>=-1+2*3+2*1-5=2.                            (33)
```

Since `chi(A2)=1`, the complete pullback is not `A2`.

Finally, it is irreducible.  If `A!=0`, then `deg_a(phi)=2`, and THM-3564's
degree-resonance gate applies.  If `A=0,B!=0`, then `deg_a(phi)=1`; the only
reducible rows are THM-3565's `h`-family, whose nonconstant members have total
degree at least four and whose constant members are affine.  If `A=B=0`,
then `C!=0` and `deg_a(phi)=0`, again excluded by THM-3564.  Thus an alleged
source-coordinate factor would have to be the whole irreducible pullback,
contradicting `(33)`.  This proves the theorem.

## 6. Exact verification

Run

```bash
python3 04-computation/jacobian_quadratic_target_graph_euler_no_go_kps_s188.py
python3 -O 04-computation/jacobian_quadratic_target_graph_euler_no_go_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion verifies the
sheet factorization `(7)`, strict transform `(9)--(10)`, discriminant `(15)`,
omitted polynomial `(30)`, constant-discriminant ideal `(23)`, all sixteen
saturated ideals `(25)--(27)`, and three independent reduced-root/Euler
controls.  In particular, a pure `Cb^2` row gives
`(rho(delta),nu,rho(f),chi(D),chi(X))=(3,1,3,-2,4)`, while two generic rows
give `(8,3,5,-9,16)`.

**QED.**
