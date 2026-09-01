---
id: THM-4312
title: "Source-normal cubic-corner repeated-face and first-splitter collapse"
status: >
  PROVED RELATIVE TO THM-4301, THM-4304, AND THM-4308 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Inside THM-4308's exact row-eight projected
  bracket/depth gate, the exact-M=12 D=Lambda=0, UZ!=0 corner has an explicit
  one-parameter response. Its balanced, Regime-A k=2, and Regime-A k=3
  repeated first faces are impossible; only k=1 survives. The surviving
  discriminant graph has first coefficient L1. Off the two equality rays the
  first carriers have rational normalization. On equality, L1 nonzero gives
  an elliptic j=0 carrier; on the genuine L1=0 locus, exact coprimality forces
  L2 nonzero and gives an elliptic j=1728 carrier. The formal good differential
  has positive order on both elliptic carriers, so any actual exact-M=12
  Keller lift realizing the finite datum maps every first carrier constantly.
  This is a finite row-eight and first-exceptional theorem, not an all-row
  lift, completed-local extinction, seam-entry, JC(2), or DC(2) theorem.
source: root / planar-Jacobian continuation session, 2026-09-01
depends_on:
  - THM-4301-cubic-corner-first-face-keller-extinction
  - THM-4304-cubic-corner-repeated-first-face-rationality
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
related:
  - THM-4307-cubic-corner-balanced-double-section-first-refinement
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4298-weighted-face-source-normal-unimodular-visibility-transform
primary_script: 04-computation/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312.py
primary_output: 05-knowledge/results/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312.out
primary_script_sha256: f5af8c0bd43d138aedbae0ed917b3925ab3b12b1581440446cc95c4a320ca5f4
primary_output_sha256: b0a960487d9afde5eba0d959f4d322c9fd20c671f14984dea43677e8fea77807
independent_audit_script: 04-computation/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312_independent_audit.out
independent_audit_script_sha256: ad30ded2c6728be70c2d401db6dd82d4392c7765176dcfebae34a89a998fcc93
independent_audit_output_sha256: 95b0dc8e62d8e036f1611b6955b9d7228179385a471bb56f80c13c556b8d0bf1
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy path derives the row-eight response, corner
  intersection, regime exclusions, literal k=1 critical graph, both weighted
  quotients, differential ledgers, and exact cancellation coprimality. A
  dependency-free Fraction/sparse-polynomial implementation reconstructs the
  source rows and repeats the same calculation, including a quadratic-field
  positive control and a firewall against the ramified-cover genus error.
  Normal and optimized streams byte-match both frozen outputs.
---

# THM-4312 -- Source-normal cubic-corner repeated-face and first-splitter collapse

**PROVED RELATIVE TO THM-4301, THM-4304, AND THM-4308 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. WITHIN THE EXACT ROW-EIGHT PROJECTED GATE, ONLY THE
REGIME-A `k=1` DOUBLE SECTIONS SURVIVE AT THE CUBIC CORNER. EVERY FIRST
REFINEMENT CARRIER IS RATIONAL OFF-TIE OR HAS POSITIVE FORMAL
GOOD-DIFFERENTIAL ORDER ON AN EQUALITY RAY; ANY ACTUAL EXACT-`M=12` KELLER
LIFT REALIZING IT MAPS IT CONSTANTLY. ALL-ROW LIFT, COMPLETED-LOCAL
EXTINCTION, SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Work in THM-4304's exact-`M=12` corner

```text
D=W^2-4UZ=0,             Lambda=U+W+Z=0,             UZ!=0. (1)
```

Also impose exactly THM-4308's finite conditions: the bracket equations
through source row eight and the projected memberships

```text
pi_8(A) in pi_8(P_2),                 pi_8(C) in pi_8(P_3). (2)
```

> **Theorem.** Among THM-4304's complete five repeated first-face regimes,
> conditions `(1)--(2)` exclude both balanced sections and the Regime-A
> `k=2,3` sections. The surviving locus is exactly the Regime-A `k=1`
> locus in Section 4. Its formal discriminant graph has one of exactly two
> first nonzero rows:
>
> ```text
> L_1!=0:       beta=s/3,       elliptic residue j=0;
> L_1=0:        L_2!=0, beta=s/2, elliptic residue j=1728. (3)
> ```
>
> In both cases the formal pullback of the good target differential has
> positive order. Consequently, any actual exact-`M=12` Keller lift realizing
> the finite datum maps the first refinement carrier constantly.

This is conditional on the finite row-eight projection `(2)`. It neither
asserts that a point extends to an all-row `B_2` lift nor that an arbitrary
hypothetical counterexample enters this exact corner.

The inheritance pass and concept board were:

- closest proved mechanisms: THM-4304's exhaustive repeated-face census and
  THM-4308's exact finite gate;
- canonical hostile: THM-4307's genuine positive-genus Regime-A refinement;
- corrected near miss: MISTAKE-540 records that a ramified root-coordinate
  cover is not the weighted exceptional function field;
- least-used sidecar: the invariant coordinate on `P(s,beta)` together with
  the finite stabilizer lost by the cover;
- live concepts: corner response, repeated-regime ideal, A1 critical graph,
  weighted quotient, elliptic `j`, good-form order, and lift/seam firewall.

## 2. Exact row-eight corner response

THM-4308 gives

```text
Delta=896/15,                  Theta=512/75,
K=-32/5,                       zeta_3=-3Phi/2,
upsilon_5=-731648/2025,                                    (4)

U=(475515904-109350xi_10)/200475,
W=-(4343625Phi^2-17172000xi_10+143826305024)/4009500,
Z=(12506118074368-173745000Phi^2-195463125Phi eta
   -926883000xi_10)/108256500.                              (5)
```

The identity

```text
W^2-4UZ=(W+2U)^2-4U(U+W+Z)                                (6)
```

shows that `(1)` is equivalent to `W=-2U`, `Z=U`. Solving these two linear
response equations gives

```text
xi_10=(4343625Phi^2+124805668864)/12798000,

Phi eta=2091705253888/107983125-(2839/1185)Phi^2,

U=Z=-13(820125Phi^2+13056802816)/57591000,
W=-2U.                                                     (7)
```

Necessarily

```text
Phi!=0,                820125Phi^2+13056802816!=0.         (8)
```

Indeed `Phi=0` contradicts the middle equation in `(7)`, while the second
inequality is exactly `U!=0`.

Put

```text
c_2=upsilon_5+xi_10
   =11(14625Phi^2+404652032)/474000,

c_3=eta+zeta_3
   =(4183410507776-841357125Phi^2)/(215966250Phi),

c_4=Delta+Theta=1664/25.                                  (9)
```

Equations `(7)--(9)` are an exact one-parameter **finite response**, not an
existence theorem for a polynomial Keller pair.

## 3. Four of the five repeated regimes disappear

THM-4304 proves that the complete repeated locus consists of two balanced
sections and Regime-A `k=1,2,3`.

For either balanced section, its forced coefficient is
`Delta=2048/45`, whereas `(4)` gives

```text
Delta-2048/45=128/9!=0.                                  (10)
```

For `k=3`, the necessary cancellation `c_4=0` contradicts `(9)`.

For `k=2`, the necessary cancellation `c_2=0` gives
`xi_10=-upsilon_5`; equation `(5)` then gives

```text
U=39636992/18225.
```

But its repeated-square condition fails by the exact residual

```text
upsilon_5^2-4U c_4=-1839105572864/4100625!=0.             (11)
```

Thus only Regime-A `k=1` remains. This conclusion uses the completeness of
THM-4304's census; it is not a generic sampling statement.

## 4. Exact surviving `k=1` locus and critical graph

The surviving locus is

```text
beta_11=-alpha_11,               alpha_11!=0,
alpha_11^2=4U c_2
 =-143(14625Phi^2+404652032)
      (820125Phi^2+13056802816)/6824533500000.             (12)
```

Hence `14625Phi^2+404652032!=0` as well. Equivalently, choose `rho` with

```text
rho^2=c_2/U,
alpha_11=-2U rho,                  beta_11=2U rho.          (13)
```

The repeated section is `q=rho t`. Since the base field is algebraically
closed, `(12)` defines a genuine two-sheeted one-parameter finite-gate locus
away from the three forbidden factors in `(8)` and `(12)`.

Put `q=ty` and `m=z^12/t^2`. Direct expansion of the literal source gives

```text
F/t^3=P_0(y)+tP_1(y)+O(t^2)-ym,

P_0=y(Uy^2+(5alpha_11+4beta_11)y+c_2)
   =Uy(y-rho)^2,

P_1=y(4Uy^3+(10alpha_11+6beta_11)y^2
       +(5upsilon_5+4xi_10)y+c_3).                        (14)
```

At `y=rho`, the simple factor `y` is a unit and the local equation is A1.
Critical recentering therefore produces one formal graph

```text
m_*(t)=L_1t+L_2t^2+O(t^3),                               (15)

L_1=c_3+rho upsilon_5
   =c_3-alpha_11 upsilon_5/(2U),                          (16)

L_2=[4U^2c_4-2Ualpha_11(4eta+3zeta_3)-Uupsilon_5^2
     +4alpha_11^2upsilon_5]/(4U^2).                       (17)
```

The primary script obtains `(14)--(17)` by expanding the source and solving
the critical-root equation; they are not inserted as unchecked graph data.

For an exact nondegenerate control, take `Phi=1`. Then

```text
xi_10=124810012489/12798000,
eta=2091446550013/107983125,
U=-169749098233/57591000,
rho^2=-98333997651/30863472406.                           (18)
```

In `Q(rho)`, the norm of `L_1=c_3+upsilon_5rho` is

```text
270260378011253985379632330934787603
------------------------------------------------ !=0.     (19)
        719758107151040278729687500
```

Thus the `L_1!=0` case is nonempty.

## 5. The high-contact locus stops at `L_2`

Let `X=Phi^2`. Eliminating `rho` from `L_1=0` and `(13)` gives

```text
R(X)=7547170421607067494140625X^3
    +164114458618573873612800000000X^2
    +2284603892441775363795663716352000X
    +2970579390109346274816679296272171008=0.             (20)
```

The exact gcd of `R` with

```text
X(820125X+13056802816)(14625X+404652032)                  (21)
```

is one. Therefore all roots of `(20)` avoid the forbidden factors, and the
high-contact locus is genuinely nonempty over the algebraic closure.

On `L_1=0`, equation `(17)` reduces to

```text
L_2=c_4+c_3zeta_3/upsilon_5-upsilon_5^2/(4U)

   =-S(X)/[676262246400(820125X+13056802816)],             (22)

S(X)=8970234157828125X^2
    +61293210070929408000X
    -1395571970793868500140032.
```

Now `gcd(R,S)=1`. A compact exact witness is the identity

```text
(-5X-8)R+(-X^2+X-5)S=1                    in F_17[X].      (23)
```

Consequently `L_1` and `L_2` never vanish together on the allowed locus.
The first nonzero discriminant row is exhausted by the dichotomy `(3)`.

## 6. Actual weighted quotients and conditional Keller extinction

The critical root itself moves. Through the order needed on both splitter
rays it is

```text
y_c(t)=rho-[upsilon_5/(2U)]t
            -[(4eta+4rho upsilon_5+3zeta_3)/(2U)]t^2+O(t^3),
Q=y-y_c(t).                                               (24)
```

The linear translation cannot be dropped: on the high-contact ray below it
has the same weight as `Q`, while the displayed `t^2` coefficient is needed
to identify the critical value through `L_2`. For `j=1` when `L_1!=0`, and
`j=2` when `L_1=0`, the local first form is

```text
aQ^2=m-L_jt^j,                           a!=0.            (25)
```

If `v(m)<jv(t)` or `v(m)>jv(t)`, this reduces respectively to `aQ^2=m` or
`aQ^2=-L_jt^j`. Each irreducible normalization is rational: indeed
`m=(z^5/sigma)^2`, while `t` for `j=1` and `t^2` for `j=2` give binomial
parametrizations. Thus only the equality rays can carry positive genus.
For an actual Keller lift, each rational carrier maps constantly to the good
elliptic target.

For `L_1!=0`, compare

```text
v(m)=10beta-2s,                  v(t)=s+beta.
```

Their difference is `3(3beta-s)`, so the splitter is `beta=s/3`, with
primitive weights

```text
(s,beta,gamma)=(3,1,4).                                  (26)
```

The invariant coordinate on the actual weighted exceptional divisor
`P(3,1)` is

```text
w=sigma/z^3.                                              (27)
```

The critical coordinate `Q` in `(24)` has excess two, so put `V=Q/z^2`.
After absorbing a nonzero local unit and setting `Y=wV`, the residue is

```text
aY^2=1-L_1w^3,                       a!=0.                (28)
```

It is a squarefree cubic, hence a genus-one curve with `j=0`. The tempting
substitution `sigma=lambda^3`, `z=lambda w_0` instead gives a degree-ten
genus-four curve on a ramified `mu_3` cover. MISTAKE-540 records why that
cover is not the carrier in `(28)`.

For `L_1=0`, equation `(22)` gives `L_2!=0`. Now

```text
v(m)-v(t^2)=4(2beta-s),
```

so the splitter is `beta=s/2`, with primitive weights

```text
(s,beta,gamma)=(2,1,3).                                  (29)
```

On `P(2,1)` use `w=sigma/z^2` and `V=Q/z^3`. After `Y=wV`,

```text
aY^2=1-L_2w^4,                       a!=0.                (30)
```

This squarefree quartic has genus one. Its binary-quartic invariants have
`J=0`, so `j=1728`. Both equality rays satisfy the inherited `k=1` strict
range `s<5beta` from THM-4304.

The inherited THM-4301 good-differential estimate is positive on both formal
elliptic carriers. For a `k=1` double root, let `d` be the common tie weight
`d=v(m)=jv(t)`. Then

```text
v(F_q)=2(s+beta)+d/2,
v(phi^*eta_0)>=9s+11beta-v(F_q).                          (31)
```

At `(26)`, `d=4`, so these orders are `10` and at least `28`. At `(29)`,
`d=6`, so they are `9` and at least `20`. A nonconstant map of smooth curves
in characteristic zero has nonzero differential. Hence any actual
exact-`M=12` Keller lift realizing either finite datum maps its first
refinement carrier constantly. In particular, the `j=0` carrier in `(28)` is
killed by the differential even though its `j` matches the good target; a
`j`-mismatch argument would miss the sharper mechanism.

## 7. Reproduction and scope

Run

```bash
python3 -B 04-computation/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312.py
python3 -B -O 04-computation/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312.py
python3 -B 04-computation/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_cubic_corner_repeated_face_collapse_thm4312_independent_audit.py
```

The exact universe is the intersection of THM-4304's exact-`M=12` corner
with THM-4308's row-eight projected gate. The theorem classifies the first
refinement carrier above every repeated face surviving that intersection.
Its differential conclusion is conditional on an actual exact-`M=12` Keller
lift realizing the finite datum. It does **not** prove membership in the
infinite depth modules, an all-row `B_2` lift, polynomial termination, later
completed-local refinements above the elliptic carriers, seam entry, `JC(2)`,
or `DC(2)`.

**QED.**
