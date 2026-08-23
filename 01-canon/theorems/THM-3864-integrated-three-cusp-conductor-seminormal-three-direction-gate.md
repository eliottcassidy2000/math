---
id: THM-3864
title: "The integrated three-cusp front has a three-dimensional seminormal defect"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED CORE; PROVISIONAL
  VERIFIED-EXACT SECOND-LINE STRENGTHENING AWAITING INDEPENDENT AUDIT.  The
  conductor of the THM-3854 integrated
  quintic has double zeros at its three cusp addresses and simple zeros at
  all six node addresses.  Its seminormalization retains exactly the three
  node-pair gluing conditions, and the quotient by the branch ring is
  three-dimensional with cusp derivatives as coordinates.  A previously
  missing third square/cube-descent direction is explicit.  Its canonical
  depressed-cubic residual and every minimum-degree parity-preserving lift
  remain nonsquare.  The audited core closes the canonical `<h_1,h_3>` line;
  the provisional strengthening closes `<h_1,h_2>` by an exact specialized
  quartic square ideal.  General mixed/higher lifts, `<h_2,h_3>`, and the
  projective-plane interior remain OPEN.
source: jc_zero_debt_lift / integrated-cusp alternative cubic-order lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit checked that
  the node-value and cusp-derivative functionals are independent and exhaust
  k[t]/R, then rederived the conductor multiplicities from arbitrary Hermite
  data and the local-to-global seminormal condition.  It verified the
  derivative-coordinate isomorphism, square/cube membership of all three
  directions, the UFD nonsquare tests, and completeness of the bounded
  parity-preserving representatives, including the a=84 boundary.  It also
  checked the canonical projective-line factorization and both endpoint
  arguments.  The companion verifies the conductor
  generator and address multiplicities, the ranks of all node/cusp
  conditions modulo the conductor, the three node-equality packets, the
  derivative-coordinate determinant, the explicit third square/cube
  descent, its residual factor, and the complete bounded parity-preserving
  lift coefficients.  It also verifies the symbolic factorization and both
  boundary specializations for the complete canonical projective line
  spanned by the first and third directions.  The new provisional section
  additionally verifies every cross descent on `<h_1,h_2>`, its specialized
  quartic, the complete square-coefficient Groebner basis, and both
  endpoints.  Normal and optimized runs must byte-match the frozen
  transcript.  Independent hostile audit is required only for this new
  second-line strengthening.
depends_on:
  - THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction
related:
  - THM-3849-russell-arm-conductor-polynomial-and-residual-contact-graph
  - THM-3855-formal-inverse-discriminant-lift-and-algebraization-gate
script: 04-computation/jc2_integrated_three_cusp_seminormal_defect_thm3864.py
output: 05-knowledge/results/jc2_integrated_three_cusp_seminormal_defect_thm3864.out
script_sha256: b0a3d3fb6804db5445e2414bfd1ee0d0b53a6da86b8fe468d8d7fd0d45581cb2
output_sha256: 71a35831d408b123e6dbbec7a1c899056c8d584bbc217a0f3cea178c938323bd
semantic_sha256: 03a59be34df363c1a050bcc5a0ff65a72c67651d9d9e76d36726ce6b10b828fc
hash_basis: raw LF bytes
---

# THM-3864 -- all three seminormal cusp directions

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED CORE, WITH A
PROVISIONAL SECOND-LINE STRENGTHENING.**  Work over an algebraically closed
field `k` of characteristic zero.  Retain the
normalization of THM-3854,

```text
x=t^4-2t^2,                         y=3t^5-5t^3,               (1)
R=k[x(t),y(t)] subset k[t].                                      (2)
```

The conductor of `R` in its normalization `k[t]` is exactly

```text
c=t^2(t^2-1)^2(3t^2-5)(9t^4-18t^2+4) k[t].                  (3)
```

Let `S` be the seminormalization of `R` in `k[t]`.  Then

```text
dim_k(S/R)=3,                                                  (4)
```

and cusp derivatives give a canonical coordinate isomorphism

```text
S/R -> k^3,
[f] |-> (f'(0),f'(1),f'(-1)).                                 (5)
```

Two directions were used as controls in THM-3854.  A convenient complete
basis is

```text
h_1=t^2(t^2-1)(9t^2-14),
h_2=t(t^2-1)(9t^2-4)(3t^2-5),
h_3=t^3(t^2-1)(3t^2-5)(9t^4+24t^2-38).                       (6)
```

Every `h_i` has square and cube in `R`.  The previously missing third
direction has the explicit descents `(23)-(24)` below.  Its simplest
depressed-cubic residual is not a square.  More strongly, no
minimum-total-degree parity-preserving change of its two polynomial lifts
makes the residual a square, and no canonical element on the entire
projective line spanned by `h_1,h_3` has square residual.

The last statements are deliberately bounded.  They do not exclude
higher-degree changes, parity-breaking changes, or combinations involving
`h_2`.

## 1. Exact global coordinate conditions

THM-3854 proves that `(1)` is the normalization of an affine rational
quintic with three `A2` cusps at

```text
t=0,1,-1,                                                     (7)
```

and three transverse nodes.  One node pairs the two roots of

```text
3t^2-5=0                         by t |-> -t.                  (8)
```

The other two pair the four roots of

```text
9t^4-18t^2+4=0                   by t |-> -2/(3t).            (9)
```

For `f in k[t]`, membership in `R` is therefore exactly the following six
linear conditions:

```text
f has equal values on each of the three node pairs,            (10)
f'(0)=f'(1)=f'(-1)=0.                                          (11)
```

Here is a dimension-safe proof that there is no hidden condition.  Every
element of `R` satisfies `(10)` because the paired addresses have the same
image, and satisfies `(11)` because both coordinate derivatives in `(1)`
vanish at the cusp addresses.  The node conditions have rank three and the
three Hermite derivative conditions add rank three: arbitrary values and
first derivatives at these distinct normalization addresses can be
interpolated by a polynomial.  Thus the subspace cut out by `(10)-(11)` has
codimension six in `k[t]`.

On the other hand,

```text
dim_k(k[t]/R)=sum_p delta_p=3 delta(A2)+3 delta(A1)=6.         (12)
```

The inclusion of `R` in the condition space and equality of codimensions
prove `(10)-(11)` in both directions.

## 2. The conductor is exactly `(3)`

The conductor consists of those `f in k[t]` for which `fh` satisfies
`(10)-(11)` for every `h in k[t]`.  At a cusp address `a`,

```text
(fh)'(a)=f'(a)h(a)+f(a)h'(a).                                 (13)
```

Since `h(a)` and `h'(a)` can be prescribed independently, `(13)` vanishes
for every `h` exactly when `f` has a double zero at `a`.  At a node pair
`a!=b`, the equality

```text
f(a)h(a)=f(b)h(b)                                             (14)
```

for independently prescribed `h(a),h(b)` holds exactly when `f` vanishes
at both addresses.  Hence the conductor has double zeros at `(7)` and
simple zeros at the six roots in `(8)-(9)`, with no other zeros.  Since
`k[t]` is a PID, this proves `(3)` including every multiplicity.

The generator is visibly not just an abstract annihilator.  It belongs to
the branch ring because

```text
27x^3+33x^2+10x+y^2
 =t^2(t^2-1)^2(3t^2-5)(9t^4-18t^2+4).                        (15)
```

Thus `(3)` is an equality of conductor ideals in the actual normalization.

## 3. Seminormalization and derivative coordinates

An `A2` cusp has local ring `k[[u^2,u^3]]`, whose seminormalization is
`k[[u]]`.  An ordinary node is already seminormal and retains the equality
of the two branch residues.  Consequently

```text
S={f in k[t]: f satisfies the three node equalities (10)}.    (16)
```

Equations `(11)` show that the derivative map `(5)` has kernel exactly `R`.
It is surjective because the three columns `(6)` have derivative matrix

```text
          h_1'  h_2'  h_3'
t=0       0     -20    0
t=1      -10    -20   20
t=-1      10    -20   20,                                   (17)
```

whose determinant is `-8000`.  This proves `(4)-(5)` and shows that the
earlier even/odd pair left one genuine cusp direction untested.

For completeness, every polynomial in `(6)` obeys both node equalities and
vanishes at all three cusp addresses.  Hence `h_i in S`, while

```text
(h_i^m)'(a)=m h_i(a)^(m-1)h_i'(a)=0,          m=2,3,          (18)
```

at every cusp.  Characterization `(10)-(11)` proves

```text
h_i^2,h_i^3 in R                                                (19)
```

without a quotient or completion loss.

## 4. The third square/cube descent

Put

```text
P_3=81x^3y^2+657x^2y^2+2356xy^2+84y^4+1444y^2,              (20)

Q_3=8991x^4y^3+21897x^3y^3+81x^2y^5+92226x^2y^3
    +5364xy^5+134292xy^3+5308y^5+54872y^3.                  (21)
```

Exact substitution into `(1)` gives

```text
P_3(t)=h_3(t)^2,                       Q_3(t)=h_3(t)^3.       (22)
```

For the quintic equation `Delta` of THM-3854, exact division yields

```text
P_3^3-Q_3^2=Delta y^6 C_3,                                   (23)

C_3=6561x^4-845640x^3-2056104x^2-10708704x
    -592704y^2-7133360.                                      (24)
```

The factor `y^6` is a square, but `C_3` is not.  As a polynomial in `y`, it
has degree two, nonzero `y^2` coefficient, no linear term, and a nonzero
constant part.  A square would be `(a(x)y+b(x))^2`; the missing cross term
forces `a(x)b(x)=0`, contradicting the other two coefficients.  Since
`k[x,y]` is a UFD and `k` is algebraically closed, the same conclusion holds
in `k(x,y)`.

Thus the canonical depressed cubic attached to `h_3` has discriminant
`108 Delta y^6 C_3`, not `Delta` times a square.

## 5. Complete minimum-degree parity-preserving lift no-go

The degree bounds in this section are part of the theorem.  The displayed
`P_3` has total degree five.  Since the branch ideal is the principal ideal
`(Delta)` with `deg Delta=5`, every other representative of `h_3^2` of
total degree at most five is

```text
P=P_3+a Delta,                          a in k.                (25)
```

The displayed `Q_3` has total degree seven.  Every representative of
`h_3^3` of total degree at most seven differs by `Delta B` with
`deg B<=2`.  Requiring the same parities `P(x,-y)=P(x,y)` and
`Q(x,-y)=-Q(x,y)` makes the complete family

```text
Q=Q_3+Delta y(b_0+b_1x),                 b_0,b_1 in k.        (26)
```

Define

```text
R_(a,b_0,b_1)=(P^3-Q^2)/Delta.                                (27)
```

It is a polynomial and is even in `y`.  Put `q=y^2`.  Three exact
coefficients control every parameter:

```text
[q^4]R=(a-84)^3,
[x^4q^3]R=6561,
[q^0]R=a^3x^6(9x+5)^4.                                      (28)
```

Suppose `R=H^2`.  Invariance under `y |-> -y` gives
`H(x,-y)=H(x,y)` or `H(x,-y)=-H(x,y)`, because `k[x,y]` is a domain.

If `a!=84`, equation `(28)` gives `q`-degree four.  An odd `H` has the form
`yK(x,q)`, so its square has odd `q`-degree, impossible.  If `H` is even,
write

```text
H=s_0(x)+s_1(x)q+s_2(x)q^2.                                  (29)
```

Since `deg R<=10`, one has `deg H<=5`.  The nonzero constant coefficient of
`q^4` in `(28)` forces `s_2` to be a nonzero constant.  Therefore the
coefficient `2s_1s_2` of `q^3` has `x`-degree at most three, contradicting
`[x^4q^3]R=6561`.

If `a=84`, the `q^4` term vanishes but the second line of `(28)` makes the
`q`-degree exactly three.  This excludes an even square.  An odd square
vanishes at `q=0`, whereas the last line of `(28)` is then nonzero.  This
excludes the final case.  Hence no member of `(25)-(26)` has square
residual.

## 6. The complete canonical line spanned by `h_1,h_3`

The third direction also closes one full projective family without changing
representatives.  In addition to `P_1,Q_1` from THM-3854 and `P_3,Q_3`
above, put

```text
P_13=81x^3y-369x^2y-266xy+46y^3,                              (30)

Q_113=2835x^4y+4797x^3y+81x^2y^3+1862x^2y
      +424xy^3+368y^3,

Q_133=5913x^4y^2-6147x^3y^2+81x^2y^4-22268x^2y^2
      +2172xy^4-10108xy^2+2116y^4.                            (31)
```

These pull back to `h_1h_3`, `h_1^2h_3`, and `h_1h_3^2`, respectively.
For `[L:N] in P1(k)`, define the canonical descents

```text
P_(L,N)=L^2P_1+2LNP_13+N^2P_3,
Q_(L,N)=L^3Q_1+3L^2NQ_113+3LN^2Q_133+N^3Q_3.                 (32)
```

Thus they pull back to `(Lh_1+Nh_3)^2` and `(Lh_1+Nh_3)^3`.  Direct exact
division and factorization give

```text
(P_(L,N)^3-Q_(L,N)^2)/Delta=(L+Ny)^2 K_(L,N),                (33)
K_(L,N)(x,0)=243L^4x^3(27x+16),                              (34)
K_(0,N)=N^4y^4 C_3.                                          (35)
```

No member of this projective line has square residual.  If `L!=0` and the
right side of `(33)` were a square in `k(x,y)`, then `K_(L,N)` would be a
square there.  The UFD property reduces this to a polynomial square up to a
constant, and every constant is already a square in the algebraically
closed field `k`.  Specializing at `y=0` would make `(34)` a square in
`k[x]`, impossible because the distinct factors `x` and `27x+16` have odd
valuations three and one.

If `L=0`, then `N!=0`, and `(33),(35)` reduce to

```text
N^6y^6 C_3.                                                   (36)
```

The scalar and `y^6` are squares, whereas `C_3` is nonsquare by Section 4.
This handles the projective endpoint and completes the line uniformly.

## 7. Provisional strengthening: the canonical `h_1,h_2` line

**PROVISIONAL + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**  The
second coordinate line also closes.  Besides the endpoint descents from
THM-3854, the cross terms are

```text
P_12=81x^2y+137xy+56y,                                       (37)

Q_112=-891x^3y-1183x^2y+81xy^3-392xy+56y^3,

Q_122=4536x^4+5040x^3-1539x^2y^2+1400x^2
      -1247xy^2+81y^4-256y^2.                                (38)
```

They pull back to `h_1h_2`, `h_1^2h_2`, and `h_1h_2^2`.  Hence the
canonical descents for `[U:V] in P1(k)` are

```text
P_(U,V)=U^2P_1+2UVP_12+V^2P_2,
Q_(U,V)=U^3Q_1+3U^2VQ_112+3UV^2Q_122+V^3Q_2.                 (39)
```

Let `R_(U,V)=(P_(U,V)^3-Q_(U,V)^2)/Delta`.  On `V!=0`, scale to
`V=1` and put `zeta=U^2`.  The specialization at `y=0` is the quartic

```text
G_zeta(x)=
 (6561zeta^3-157464zeta^2+1259712zeta-3359232)x^4
+(3888zeta^3-108864zeta^2-124416zeta-7464960)x^3
+(-9576zeta^2-1304640zeta-6220800)x^2
+(-470400zeta-2304000)x-320000.                              (40)
```

If `R_(U,V)` were a square in `k(x,y)`, the UFD argument used in Section 6
would make it a polynomial square up to a constant, and the constant is a
square in `k`.  Since `(40)` is never zero, it would have the form

```text
G_zeta(x)=(A_2x^2+B_1x+C_0)^2.                               (41)
```

Coefficient comparison in `(41)` has the following reduced lexicographic
Groebner basis:

```text
400A_2+243C_0zeta-1296C_0,
200B_1-147C_0zeta-720C_0,
C_0^2+320000,
zeta^2.                                                       (42)
```

Thus `zeta=0`, so `U=0`.  This endpoint is the pure `h_2` residual

```text
-(9x+5)^2
 (41472x^2+6561xy^2+46080x+3888y^2+12800),                  (43)
```

which is nonsquare: after removing the displayed square, the remaining
factor has `y`-degree two, no linear `y` term, and nonzero `y^2` and
constant parts.  If `V=0`, one is at the pure `h_1` residual

```text
6561x^4+3888x^3-512y^2,                                     (44)
```

nonsquare by the same argument.  Therefore every point of the canonical
`<h_1,h_2>` projective line has nonsquare residual, subject only to the
pending independent audit of this section.

## 8. Scope and next exact search

The result converts the seminormal boundary problem from an incomplete
two-direction scout to a complete three-dimensional defect space.  It does
**not** classify binary cubic orders with discriminant `Delta`.  The first
live searches are:

1. the remaining `<h_2,h_3>` coordinate line and the interior of the full
   projective plane;
2. parity-breaking representatives of the same boundary elements;
3. representatives above total degrees five and seven; and
4. the all-orders algebraization mechanism of THM-3855 applied to this
   one-place quintic rather than its four-ray front.

Any positive result must still pass the constant-unit and Keller-atlas gates
of the current JC frontier.  No cubic completion or Jacobian counterexample
is claimed.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_integrated_three_cusp_seminormal_defect_thm3864.py
python3 -O 04-computation/jc2_integrated_three_cusp_seminormal_defect_thm3864.py
```

Both runs must byte-match
`05-knowledge/results/jc2_integrated_three_cusp_seminormal_defect_thm3864.out`.
