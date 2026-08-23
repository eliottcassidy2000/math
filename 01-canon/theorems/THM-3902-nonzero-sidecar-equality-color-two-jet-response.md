---
id: THM-3902
title: "Nonzero-sidecar equality-color two-jet response"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  On the positive equal y-degree seam of the
  THM-3881 residual, the equianharmonic leading colors carry an exact
  two-jet marked-response law.  The term 2K^2L^2f^3 enters as a source at
  response depth n.  For n=1 its first response has exactly one less total
  a=x+1 valuation than the leading color product.  This is a necessary
  response filtration only: explicit controls lift both displayed jets,
  and neither the equality seam nor JC(2) is closed.
source: tournament-jc-breakthrough / post-THM-3899 marked-response scout, 2026-08-23
audit: >
  SELF-AUDITED PROVISIONAL CANDIDATE.  The focused exact companion expands
  the full THM-3881 residual, proves the all-degree support cutoff through
  two y-jets, verifies the epsilon and denominator-cleared color laws, and
  checks the a-valuation mechanism.  It includes zero-sidecar hostiles at
  n=1,2, zero-sidecar positive lifts for n>=3, and an address-compatible
  n=1 positive two-jet control.  Normal and optimized runs must byte-match
  the frozen transcript.  Independent audit must recheck the support
  exhaustion, formal division by V, polynomial sidecar clearing, and the
  exact distinction between a two-jet lift and a square residual.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
related:
  - THM-3899-nonzero-sidecar-y-degree-tariff-and-equianharmonic-equality-colors
  - THM-3377-path-colour-deletion-compiler-and-skew-current
  - THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries
script: 04-computation/jc2_nonzero_sidecar_equality_color_two_jet_response_thm3902.py
output: 05-knowledge/results/jc2_nonzero_sidecar_equality_color_two_jet_response_thm3902.out
script_sha256: e37e39b49a757efac8b8d7cfec9eda02eb944dbe3033e8222c27128c41cca526
output_sha256: 210392b7b81a2b04d0c5252fea4852bb3410be02fd6cfca35f31a96126d7eb5b
semantic_sha256: a92f4bcb59ddfe7488b0383916d4a11dc44a6f2fe6634e7b86cb821d9c8d8a5d
hash_basis: raw LF bytes
---

# THM-3902 -- the equality colors have a mandatory response jet

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Use the THM-3881 coordinates

```text
D=k[x,y],                  a=x+1,
L=9x+4,                    K=y^2+kappa,
kappa=-15x^2-15x-4,        P=aL^2,
r=aT+Kf,                   A=KT+aPf.                       (1)
```

The complete residual is

```text
S(T,f)=L^4
 +2(3A+3P+r^2)L^2f
 +(8A+6P+3r^2)(Pf^2-T^2).                                (2)
```

Suppose

```text
f!=0,                      S(T,f)=G^2,
T(0,0)=4f(0,0),            deg_y(T)=deg_y(f)=n>=1.        (3)
```

This candidate determines the next two `y`-coefficients after the leading
equianharmonic norm.  It does not assert that the resulting jets extend to a
square.

## 1. Coordinates at `y=infinity`

Put `epsilon=y^(-1)` and write

```text
T=y^n U,       U=u+u_1 epsilon+u_2 epsilon^2+O(epsilon^3),
f=y^n V,       V=v+v_1 epsilon+v_2 epsilon^2+O(epsilon^3),
G=y^(2n+2) Gamma,
Gamma=g+g_1 epsilon+g_2 epsilon^2+O(epsilon^3),           (4)
```

where the displayed coefficients lie in `k[x]` and `u,v` are nonzero.  When
`n=1`, set `u_2=v_2=0`; there are no such coefficients in `T,f`.

The leading coefficient equation obtained directly from `(2)` is

```text
g^2=3v^2(aL^2v^2-u^2).                                   (5)
```

As in the leading calculation behind THM-3899, the factor in parentheses
cannot vanish: `u^2=aL^2v^2` has odd `a`-valuation.  UFD valuations give
`v|g`.  After changing the sign of the square root if necessary, write

```text
g=vh,                    h^2=3(aL^2v^2-u^2).             (6)
```

Formal division is now legal in `k(x)[[epsilon]]`.  Define

```text
H=Gamma/V=h+j_1 epsilon+j_2 epsilon^2+O(epsilon^3),       (7)
j_1=(g_1-hv_1)/v,
j_2=(v(g_2-hv_2)-v_1(g_1-hv_1))/v^2.                    (8)
```

The `j_i` need not lie in `k[x]`; their denominator-cleared numerators will
be retained below.

Fix `d in k` with `d^2=-3` and define the two color series

```text
C_-=H-dU,                 C_+=H+dU.                       (9)
```

Their leading coefficients are

```text
c_-=h-du,                 c_+=h+du,
c_-c_+=3aL^2v^2.                                           (10)
```

## 2. Exact epsilon response

Modulo `epsilon^3`, the square equation has the compact form

```text
C_-C_+
 =3aL^2V^2
  +2L^2 epsilon^n V
  +6epsilon^2(kappa+aU/V)(aL^2V^2-U^2)
  +O(epsilon^3).                                          (11)
```

The middle term is the marked response.  It appears at the first jet for
`n=1`, at the second jet for `n=2`, and below the displayed window for every
`n>=3`.

To prove `(11)`, assign to a residual monomial `T^p f^q K^r` its deficit from
the leading degree `4n+4`:

```text
delta(p,q,r)=(4-p-q)n+4-2r.                               (12)
```

The only monomials with deficit at most two for some `n>=1` are

```text
3aL^2K^2f^4-3K^2T^2f^2,         deficit 0,
2L^2K^2f^3,                      deficit n,
6a^2L^2KTf^3-6aKT^3f,           deficit 2.               (13)
```

Every other monomial in the exact expansion of `(2)` has nonnegative slope
in `n` and deficit at least three already at `n=1`.  Now use

```text
K=y^2(1+kappa epsilon^2).                                 (14)
```

After division by `y^(4n+4)V^2`, the first line of `(13)` contributes

```text
3(1+kappa epsilon^2)^2(aL^2V^2-U^2),                     (15)
```

the second contributes

```text
2L^2epsilon^n(1+kappa epsilon^2)^2V,                     (16)
```

and the third contributes

```text
6a epsilon^2(1+kappa epsilon^2)(U/V)(aL^2V^2-U^2).       (17)
```

The left side is `H^2`.  Add `3U^2`, use
`C_-C_+=H^2+3U^2`, and truncate after `epsilon^2`; this gives `(11)`.

## 3. The two raw square-root equations

Let

```text
E=aL^2v^2-u^2,          delta_i=[n=i].                   (18)
```

Write `S_i=[y^(4n+4-i)]S`.  Direct coefficient extraction gives

```text
S_0=3v^2E,                                                (19)

S_1=6v((2aL^2v^2-u^2)v_1-uv u_1)
     +2delta_1 L^2v^3,                                   (20)

S_2=3{aL^2(4v^3v_2+6v^2v_1^2)
       -u^2(v_1^2+2vv_2)
       -4uu_1vv_1-(u_1^2+2uu_2)v^2}
     +6(kappa v^2+auv)E
     +6delta_1 L^2v^2v_1+2delta_2 L^2v^3.                (21)
```

Thus `(3)` requires

```text
S_1=2vhg_1,             S_2=g_1^2+2vhg_2.                (22)
```

Equations `(20)--(22)` are the unnormalized form of the response law.  They
make polynomial divisibility visible but combine the two colors.

## 4. Marked color responses with polynomial sidecars

Write

```text
c_(+/-),1=j_1 +/- d u_1,
c_(+/-),2=j_2 +/- d u_2.                                 (23)
```

Taking coefficients in `(11)` gives the normalized laws

```text
c_- c_(+,1)+c_+ c_(-,1)
 =6aL^2vv_1+2delta_1L^2v,                                (24)

c_- c_(+,2)+c_+ c_(-,2)+c_(-,1)c_(+,1)
 =3aL^2(v_1^2+2vv_2)
  +2L^2(delta_1v_1+delta_2v)
  +2(kappa+au/v)h^2.                                     (25)
```

To preserve all denominator and overlap data, put

```text
J_1=g_1-hv_1=vj_1,
J_2=v(g_2-hv_2)-v_1J_1=v^2j_2,

R_(+/-),1=J_1 +/- dv u_1,
R_(+/-),2=J_2 +/- dv^2u_2.                               (26)
```

These are polynomials in `k[x]`.  Clearing denominators in `(24),(25)` gives
the exact polynomial response laws

```text
c_-R_(+,1)+c_+R_(-,1)
 =2v^2(3aL^2v_1+delta_1L^2),                             (27)

c_-R_(+,2)+c_+R_(-,2)+R_(-,1)R_(+,1)
 =3aL^2v^2(v_1^2+2vv_2)
  +2L^2v^2(delta_1v_1+delta_2v)
  +2v(kappa v+au)h^2.                                    (28)
```

In particular, modulo either color, `(27)` prescribes the response attached
to that color after multiplication by the other.  If `c_-` and `c_+` share
factors, their gcd is indispensable sidecar data; one cannot recover the
individual congruences from the product `(10)` alone.  This is the precise
content imported from the marked-response idea.  No tournament statement is
used as a proof dependency.

## 5. Exact deletion of one `a`-order when `n=1`

Let

```text
r_-=ord_a(c_-),           r_+=ord_a(c_+),
e=ord_a(v).                                                (29)
```

Since `a` and `L` are coprime, `(10)` gives

```text
r_-+r_+=1+2e.                                             (30)
```

For `n=1`, the right side of `(27)` is

```text
2L^2v^2(1+3av_1).                                        (31)
```

Both `L` and `1+3av_1` are units at `a`; therefore

```text
ord_a(c_-R_(+,1)+c_+R_(-,1))=2e.                         (32)
```

The first marked response has exactly one less total `a`-order than the
leading product.  In the clean case `e=0`, exactly one color is divisible by
`a`.  If, for example, `r_-=1` and `r_+=0`, reduction of `(27)` modulo `a`
forces `R_(-,1)` to be an `a`-unit.  Thus the response belonging to the
`a`-marked color deletes its simple `a`-factor.  When `e>0`, equation `(32)`
retains the overlap rather than assigning it cosmetically to one color.

## 6. Exact hostile and positive controls

The canonical leading payment from THM-3899 is

```text
v=1,
h_*=(a+3L^2)/2,          u_*=(3L^2-a)/(2d),              (33)
h_*-du_*=a,              h_*+du_*=3L^2.                  (34)
```

Here `h_*` is a nonconstant quadratic and `gcd(h_*,L)=1`.

### 6.1. Zero-sidecar hostiles

Take `u_1=v_1=u_2=v_2=0`.

* If `n=1`, equation `(20)` requires

  ```text
  h_*g_1=L^2,                                               (35)
  ```

  impossible because `h_*` is a nonunit coprime to `L`.

* If `n=2`, the first response lifts with `g_1=0`, but `(21)` requires

  ```text
  g_2=(kappa+au_*)h_*+L^2/h_*,                             (36)
  ```

  which is again not polynomial.

These pairs satisfy the address because their constant terms vanish.  They
show that the leading colors do not determine the response.

### 6.2. A simple positive lift for every `n>=3`

For `n>=3`, use the same zero sidecars and take

```text
f=y^n,                    T=u_*y^n,
g=h_*,                    g_1=0,
g_2=(kappa+au_*)h_*.                                      (37)
```

Equations `(19)--(22)` hold through the second response.  This is a positive
two-jet lift, not a square residual and not an existence statement.

### 6.3. An address-compatible `n=1` positive two-jet lift

The marked sidecar can also repair both displayed equations at `n=1`.  Put

```text
p=-2/(147(1+4d)),
v_1=ph_*,
j_1=av_1+1/3,             u_1=dj_1/3,                    (38)

j_2=(kappa+au_*)h_*
     +(L^2p/2)(3aph_*+2),                                 (39)

g_1=h_*v_1+j_1,           g_2=j_2+v_1j_1.                (40)
```

The scalar `1+4d` is nonzero because `d^2=-3`.  Since
`h_*(0)=49/2`, direct reduction in `k[d]/(d^2+3)` gives

```text
u_1(0)=4v_1(0).                                           (41)
```

Hence

```text
f=y+v_1,                  T=u_*y+u_1                     (42)
```

satisfies the THM-3881 address.  Substitution verifies

```text
S_0=h_*^2,
S_1=2h_*g_1,
S_2=g_1^2+2h_*g_2.                                       (43)
```

Thus even the exceptional `n=1` source is a response requirement, not an
automatic contradiction.  Nothing here asserts that coefficients below
`g_2` can be chosen.

## 7. Scope

The result is a necessary two-jet law on `m=n>=1`.  It does not prove that
any positive control extends to `S(T,f)=G^2`, close the equality seam, treat
`m>n` or `m=n=0`, produce a Keller atlas, or settle JC(2).  The theorem also
does not identify colors through their product alone: the two valuations,
their gcd, and the marked responses in `(26)` are load-bearing coordinates.

The next exact task is the response at depth `n`, where the source
`2L^2epsilon^nV` first appears for every common degree and where additional
residual monomials must be reintroduced according to `(12)`.

## 8. Reproduction

```bash
python3 04-computation/jc2_nonzero_sidecar_equality_color_two_jet_response_thm3902.py
python3 -O 04-computation/jc2_nonzero_sidecar_equality_color_two_jet_response_thm3902.py
```

Both streams must byte-match
`05-knowledge/results/jc2_nonzero_sidecar_equality_color_two_jet_response_thm3902.out`.
The companion has `83` active gates.  Its finite cases `n=1,2,3,4` freeze
the coefficient formulas; the monomial-deficit gate proves the support
cutoff for all `n>=1`.
