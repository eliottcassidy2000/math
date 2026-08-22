# Nonlinear integral Hamiltonian response: one-root classification and the split-root first window

**Status:** derivation companion for
[THM-3422](../01-canon/theorems/THM-3422-one-root-nonlinear-integral-hamiltonian-response.md),
now **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**.  The one-root module
classification and the repeated split-root first-window formula below have
complete proofs and exact replays.  A full multiple-root module decomposition
is deliberately **not** claimed.

## 1. Inheritance, portfolio, and connection contract

The closest proved mechanism is the exact sector colimit in
[THM-3418](../01-canon/theorems/THM-3418-one-monomial-nonlinear-fiber-keller-classification.md).
[THM-3419](../01-canon/theorems/THM-3419-generic-kummer-response-regular-sector-rank.md)
then proves that every sector has generic rank `N=deg(rad(g))`.  The canonical
hostile is `P=x+x^2z^2`: all generic ranks are one, but its two integral
sectors have different module types.  The corrected near miss is
MISTAKE-374: localization can kill the distinguished class `[1]` without
producing a polynomial mate.  The least-used relevant sidecar is the wrap
sector's quotient tower `K[x]/(g^k)`.

The applicable method card is **Separate descent, ambient scale, and
regularity debt** from `META-PATTERNS.md`: the generic regular representation
is the descended rank statement, while the integral intersection owes a
separate torsion/divisibility calculation.  The split-root hostile below is a
new instance of the card's counterwarning that equal rank hides integral
regularity debt; it is not yet evidence for a new meta-pattern card.

The live concept board was:

| object | predicate | invariant / operation | missing coordinate | cheapest test |
|---|---|---|---|---|
| THM-3418 sector colimit | integral module type | monomial transition | orbit label across stages | `g=x^e` |
| THM-3419 regular packet | equal generic ranks | localization | boundary thickness | compare `e=2,3` at `d=2` |
| THM-3412 Prüfer arms | persistent torsion | denominator depth | nonlinear character shift | regress to `d=1` |
| split Kummer fiber | first boundary window | character local system | all root monodromies | `(d;e_1,e_2)=(4;3,1)` |

The connection is exact:

| item | content |
|---|---|
| source | the integral sector `C_(s-1)` of `C_P=K[x,z]/D_P(K[x,z])` |
| target | its generic THM-3419 line and its boundary/Prüfer lattice |
| map | localization `K[P] -> K(P)`, and reduction modulo `P-beta_i` |
| preserved | sector character and generic rank |
| destroyed by localization | `(P-beta_i)`-primary torsion, integral divisibility, and filtration depth |
| destroyed by one-root localization | monodromy around every other root |
| required sidecars | the orbit label `ell=m-k(e-1)` and, for split roots, the full Kummer exponent vector |
| decisive hostiles | `d=2,e=2` versus `d=2,e=3`; then `(4;3,1)` at character `s=2` |

This is the generic-regular-packet/boundary-nonuniformity pattern in exact
form: multiplicity is invisible to every generic rank but controls which
characters carry torsion and how quickly that torsion thickens.

## 2. Exact one-root theorem

Let `K` have characteristic zero, let `d>=2`, and take

```text
P=ax+b+c(x-alpha)^e z^d,       a,c in K*,       e>=1.       (1)
```

Put

```text
y=x-alpha,       beta=a alpha+b,
lambda=(P-beta)/a=y+gamma y^e z^d,       gamma=c/a,
h=e-1.                                                     (2)
```

Then `K[P]=K[lambda]`.  Let `C_(s-1)` denote the Hamiltonian target sector of
fiber weight `s-1 mod d`, where `1<=s<=d` and `s=d` names the wrap sector.

If `e=1`, then every sector is free:

```text
C_(s-1) ~= K[lambda].                                     (3)
```

Suppose `e>1`.  Then

```text
C_(s-1) ~=
  K[lambda] direct-sum K[lambda,lambda^(-1)]/K[lambda]
       if d divides s(e-1),

  K[lambda,lambda^(-1)]
       if d does not divide s(e-1).                        (4)
```

The isomorphisms in `(3)--(4)` require scalar choices along monomial chains;
the sector, its torsion submodule, and the filtration below are canonical.
In particular, if

```text
g_0=gcd(d,e-1),       Pr_lambda=K[lambda,lambda^(-1)]/K[lambda],
```

then

```text
C_P ~= (K[lambda] direct-sum Pr_lambda)^g_0
       direct-sum K[lambda,lambda^(-1)]^(d-g_0),            (5)
tors_(K[P])(C_P) ~= Pr_lambda^g_0.                          (6)
```

Every summand in `(4)` has `K(lambda)`-rank one, so `(5)` recovers THM-3419's
`d`-dimensional regular generic packet while exposing a nonuniform integral
character distribution.

### Finite torsion thickness

In a selected sector put

```text
q_s=s(e-1)/d.
```

Let `F_(s,k)` be the image in the colimit of its depth-`k` polynomial stage,
whose target fiber exponent is `s-1+kd`.  Then

```text
dim_K(F_(s,k) intersect tors(C_(s-1)))=k(e-1)+q_s.          (7)
```

Thus all selected characters have the same multiplicity slope `e-1`, but
their intercepts `q_s` differ.  Multiplicity does **not** create equal Prüfer
arms in all `d` characters: it selects the shifted subgroup

```text
{s-1 mod d : s(e-1)=0 mod d},                              (8)
```

which has size `gcd(d,e-1)`.

### The distinguished unit observer

Let `theta=[1] in C_0`.  It is always nonzero when `e>=1`, in agreement with
THM-3418's absence of a polynomial mate for nonconstant `g`.  Its exact
annihilator is

```text
Ann_(K[lambda])(theta)=
  (lambda^((e-1)/d))       if d divides e-1 and e>1,
  (0)                      otherwise.                     (9)
```

Consequently

```text
theta tensor K(lambda)=0  iff  d divides e-1 and e>1.     (10)
```

This is a sharp nonlinear instance of MISTAKE-374.  Generic vanishing in
`(10)` is only necessary for a mate.  The integral class remains nonzero, so
no polynomial `Q` with `D_P(Q)=1` follows.

When `e-1=qd`, the upper annihilator in `(9)` has an explicit primitive.  If
`delta=D_P/a`, then

```text
Q_q=sum_(j=0)^(q-1)
      [binom(q-1,j)/(1+jd)] gamma^j
      y^(q+j(e-1)) z^(1+jd)                               (11)

delta(Q_q)=lambda^q.                                      (12)
```

The orbit proof below shows that no smaller power kills `[1]`.

## 3. Proof of the one-root theorem

Divide the Hamiltonian differential by the unit `a`.  Direct differentiation
gives

```text
delta=(1+gamma e y^(e-1)z^d) partial_z
      -gamma d y^e z^(d-1) partial_y.                      (13)
```

For `n=s+kd`, define the target monomial class

```text
E_(k,m)=[y^m z^(s-1+kd)].                                  (14)
```

Applying `(13)` to `y^m z^n` gives the exact relation

```text
n E_(k,m)+gamma(ne-dm)E_(k+1,m+e-1)=0.                    (15)
```

For `s<d`, one has `m>=0`.  In the wrap sector `s=d`, the derivative of the
fiber-exponent-zero stage gives exactly the THM-3418 quotient, so

```text
0<=m<e(k+1).                                               (16)
```

Assume first that `h=e-1>0`.  Relation `(15)` preserves

```text
ell=m-kh.                                                  (17)
```

For each fixed integer `ell`, admissible states occur for every sufficiently
large `k`, both in the ordinary sectors and under the wrap bound `(16)`.  The
transition coefficient is, up to a nonzero scalar,

```text
dm-en=d ell-es-dk.                                         (18)
```

It vanishes at most once.  A one-dimensional chain with one zero transition
has a killed initial segment and a one-dimensional surviving tail; a chain
with no zero has the same one-dimensional colimit.  Hence the sector colimit
is the algebraic direct sum of one line for every `ell in Z`.

Multiplication by `lambda` sends the `ell` line to the `ell+1` line.  Indeed,
write

```text
lambda E_(k,m)=A+gamma B,
A=E_(k,m+1),       B=E_(k+1,m+e).
```

Relation `(15)` applied to the input monomial `y^(m+1)z^n` gives, at any
sufficiently late stage,

```text
lambda E_(k,m)
 =[d(m+1)-n(e-1)]/[d(m+1)-ne] A.                          (19)
```

The numerator loses all dependence on `k`:

```text
d(m+1)-n(e-1)=d(ell+1)-s(e-1).                            (20)
```

Therefore the `lambda` arrow vanishes exactly when

```text
d divides s(e-1),       ell=ell_0=s(e-1)/d-1.             (21)
```

If there is no break, the bilateral chain of lines is
`K[lambda,lambda^(-1)]`.  If `(21)` holds, the half-chain
`ell<=ell_0` is the Prüfer module and the half-chain `ell>=ell_0+1` is a free
polynomial line.  Rescaling the nonzero arrows to one proves `(4)`.

When `e=1`, the shift `h` is zero.  The orbit labels are instead
`ell=m>=0`.  No ordinary-sector transition can vanish because `1<=s<d`; the
wrap quotient admits the `m` line from depth `m` onward.  Formula `(20)` is
now `d(m+1)`, which never vanishes.  Every sector is the unilateral free chain
`K[lambda]`, proving `(3)`.

The congruence in `(21)` has `gcd(d,e-1)` solutions for `s mod d`, proving
`(5)--(6)` and `(8)`.  Its endpoint has a direct polynomial witness.  With
`q_s=s(e-1)/d`,

```text
tau_s=[y^(q_s-1)z^(s-1)],
lambda tau_s=delta(y^q_s z^s/s)=0.                         (22)
```

In a selected nonwrap sector, `(18)` cannot vanish: divisibility of both
`s(e-1)` and `se` by `d` would imply `d|s`.  In the wrap sector the formal
zero occurs exactly at the excluded monomial `m=e(k+1)`.  Thus the depth-`k`
stage meets the Prüfer half-chain in the uninterrupted interval

```text
-k(e-1)<=ell<=ell_0,
```

whose length is `(7)`.

Finally `[1]` occupies `s=1,ell=0`.  If `d` does not divide `e-1`, it lies on
a Laurent line and is nontorsion.  If `e-1=qd`, the unique zero arrow is at
`ell=q-1`, so exactly `q` applications of `lambda` kill it.  This proves the
minimal annihilator `(9)`.  In `(11)`, the low coefficient of the `j`th term
under `delta` is `binom(q-1,j)`, while its high coefficient is the same number;
adjacent terms add by Pascal's identity and give `(y+gamma y^e z^d)^q`.
This independently proves `(12)`.  QED.

## 4. Exact split-root first-window theorem

The one-root answer does **not** tensor together root by root.  There is,
however, a complete first boundary window.

Let `K'/K` be a faithful splitting extension containing the `d`th roots of
unity, put `C'=C_P tensor_K K'`, and write

```text
g(x)=c product_(j=1)^N (x-alpha_j)^e_j                  (23)
```

with distinct roots.  Put `beta_i=a alpha_i+b` and fix a repeated root
`e_i>1`.  For `1<=s<=d`, define

```text
v_i=e_i-1,              v_j=e_j for j!=i,
c_i=gcd(d,v_1,...,v_N).                                  (24)
```

Then the first `lambda_i=P-beta_i` window of target sector `s-1` is

```text
dim_(K') C'_(s-1)/(P-beta_i)C'_(s-1)
 = N       if d divides s v_j for every j,
 = N-1     otherwise.                                    (25)
```

Exactly `c_i` of the `d` characters take the value `N`; the other `d-c_i`
take the value `N-1`.

### Proof

Set `y=x-alpha_i`.  The special fiber factors as

```text
P-beta_i
 =y [a+c y^(e_i-1) product_(j!=i)(x-alpha_j)^e_j z^d].   (26)
```

Because `e_i>1`, the two factors in `(26)` are comaximal.  On the vertical
component `y=0`, the differential is `a partial_z`, whose cokernel vanishes in
every sector.  On the other component, `y`, `z`, and every
`x-alpha_j` are units.  The component is the finite etale Kummer cover of

```text
U=A1 minus {alpha_1,...,alpha_N}                           (27)
```

given by

```text
z^d=-a/[c y^(e_i-1) product_(j!=i)(x-alpha_j)^e_j].        (28)
```

Here `P_z` is a unit.  The THM-3419 map from relative one-forms to the
Hamiltonian cokernel is therefore an isomorphism on this component, with its
same one-character shift: target weight `s-1` corresponds to de Rham
character `s`.

The character-`s` rank-one local system on `U` has monodromy exponents
`s v_j mod d`.  It is trivial exactly under the first line of `(25)`.  Since
`U` is a wedge of `N` circles, its trivial rank-one system has `H^1` dimension
`N`; every nontrivial rank-one system has `H^0=0` and `H^1` dimension `N-1`.
This proves `(25)`.  Equivalently, the covering graph has `d` vertices, `N`
translation edges out of each vertex, `c_i` connected components, and first
Betti number

```text
dN-d+c_i=d(N-1)+c_i,                                      (29)
```

the sum of the character dimensions in `(25)`.  QED.

## 5. The sharp gluing hostile and collision ledger

Take two split roots with

```text
d=4,             (e_i,e_j)=(3,1).                        (30)
```

At the repeated root, the one-root congruence sees `e_i-1=2`, so character
`s=2` is locally resonant.  A false rootwise direct-sum rule would assign one
free special-fiber dimension to this root and one to the other root, predicting
dimension two.  The actual exponent vector is `(2,1)`.  Character `s=2` still
has nontrivial monodromy around the other root, and `(25)` gives dimension
`N-1=1`.  Only `s=4` has dimension two.

Thus:

- local one-root arms do not determine the global integral module;
- the missing coordinate is the complete puncture-monodromy vector;
- zero local residue does not make a global Kummer period exact; and
- `(25)` is a quotient window, not a torsion-kernel or Prüfer-persistence
  theorem.

There are no literal boundary-value collisions among distinct split roots:
`beta_i=a alpha_i+b` is injective because `a!=0`.  Over a nonsplit base field,
conjugate roots instead descend as an irreducible support packet.  The genuine
obstruction is therefore global period coupling, not coincident `P`-values.

## 6. Boundaries and non-consequences

- If `g` is a nonzero constant, or if `g=0`, then `K[x,z]=K[P,z]` and
  `D_P=a partial_z`; the complete response is zero.  These are not `e=0`
  instances of `(3)--(6)`.
- The proof works at `d=1` as a regression: for `e>1` it gives exactly one
  free line plus one Prüfer arm, agreeing with the one-root case of THM-3412.
- The split-root theorem requires the chosen root to be repeated.  At a simple
  root the two factors in `(26)` intersect at critical points, so the product
  and etale-cover proof is invalid.
- Formula `(25)` computes only reduction modulo one boundary prime.  It does
  not determine higher powers, the torsion kernel, extensions between generic
  periods, or a global decomposition of `C_P`.
- The class `[1]` is the mate obstruction, but its generic vanishing in `(10)`
  is not an integral primitive.  THM-3418 remains the Keller classification.
- Nothing here treats intermediate `z` coefficients or proves `JC(2)`.

## 7. Exact companion

The standard-library companion uses only `Fraction` and integer arithmetic.
It checks `34,504` direct bivariate monomial relations, `280` direct torsion
endpoint identities, `50` explicit unit-observer primitives, `30,576` orbit
transition identities, `30,576` lambda-arrow/break identities, `2,268` finite
torsion-thickness windows, and `1,260` simple-root controls.  Independently,
it checks `21,060` repeated split-root packets and `126,360` special-fiber
character dimensions against the covering-graph component and Betti counts.
The normal and optimized outputs are byte-identical.

Reproduce with

```text
python3 04-computation/jc_one_root_nonlinear_integral_response_thm3422.py
python3 -O 04-computation/jc_one_root_nonlinear_integral_response_thm3422.py
```

Artifacts:

- `04-computation/jc_one_root_nonlinear_integral_response_thm3422.py`
  (`sha256 1c862dd5dfaee00a3ee5827a9d004c160a574deb103cd3b878a6615aeac7a766`)
- `05-knowledge/results/jc_one_root_nonlinear_integral_response_thm3422.out`
  (`sha256 1e0abbff4bbd2ab601cbf2aee4c840789f4ebf36a1444bf2aa8c8506579e485d`)
