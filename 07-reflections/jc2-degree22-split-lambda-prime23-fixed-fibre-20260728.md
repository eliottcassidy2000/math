# Degree twenty two after the nonsplit closure: the split even-Faber first-flux constant creates a prime-23 curve

**Status:** VERIFIED-EXACT STRUCTURAL SCOUT ON THE SPLIT EVEN-FABER SUBCHART.
This note identifies and checks the continuation of THM-2692 obtained by
retaining the split first-flux constant while setting every odd Faber seed to
zero.  The full split deck permits an additional odd-seed bank.  This note
does **not** prove that bank absent, prove the split branch empty, close degree
twenty two, or prove `JC(2)`.

Companion:

```text
04-computation/jc2_degree22_split_lambda_prime23_scout.py
05-knowledge/results/jc2_degree22_split_lambda_prime23_scout.out
script SHA-256 1bbadb900e27112f57f600b0b5b73a3b85fc5b02f23e3f8f687dac4fa1c41fc3
output SHA-256 2c22072ac4e0419a83fe72c2529aa8f05922ac20c4582a443b0d8f96cd7cdfa4
```

## 1. Inheritance pass and the missing coordinate

The closest proved mechanism is THM-2692: when nonsplit parity forces the
anti-invariant first flux to vanish, the two normalized flux equations are
`N1=N2=0`; their full-support eliminant is absolutely irreducible, and the
fixed-place Kummer argument closes the inherited branch.

The canonical hostile is THM-2214 Section 4 / THM-2206: on the split etale
deck `K x K`, the constant field has an anti-invariant line

```text
{(lambda,-lambda):lambda in C}.
```

Within the inherited even-Faber bank, the object forgotten by the nonsplit
calculation is not a sixth coefficient-support parameter.  It is the
first-flux constant `lambda`.

This is not the whole split deck.  THM-2230 Section 3.1, equations (11f)--
(11g), permits every Faber degree not divisible by four, while target shears
remove exactly `E_(4j)=P^j`.  THM-2202 Section 3, equation (34), then displays
the odd bank explicitly.  At degree twenty two it is

```text
E1,E3,E5,E7,E9,E11,E13,E15,E17,E19,E21.              (0)
```

THM-2189 equations (10)--(15) delete these terms only on the nonsplit field
deck: `sigma(E_j)=(-1)^j E_j` and all constants there are fixed.  THM-2206
equation (16) explains the failure on the split etale deck, whose constants
contain the anti-invariant line `{(a,-a)}`.  Each odd `E_j` can therefore have
an anti-invariant coefficient and still contribute a deck-fixed product.  On
one chosen sheet its coefficient is an ordinary complex constant.

The prime-23 scout is the exact coordinate section where all eleven odd
coefficients vanish.  This restriction is stable under `P -> P+c`, because
translation changes a Faber degree only by multiples of four.  Any full split
argument must restore the bank (0) or prove an order-raising/descent theorem
that kills it.
The live concept board is

```text
lambda; the odd-seed bank; the q -> T=q^2 -> Z=T^2 quotient ladder;
N1/N2; weight 23; the old L5 and new G3 fixed branches;
the Z_2[C2] filtration.
```

## 2. Exact split-sheet equation

After restricting to THM-2411's even-Faber bank, the exact Faber rows give

```text
Phi_Q=q*(-N1/7496192),
Psi_Q-Psi0=N2/1319329792.                                (1)
```

Here `q` is the deck-character coordinate, `q^2=T`, and `T^2=Z`.
Choose one split sheet and write `Phi_Q=lambda`.  On this subchart, equation
(1) is exactly

```text
q N1=-7496192 lambda,                 N2=0.             (2)
```

No squaring is used in (2).  The companion reconstructs all six Faber rows
in degrees `2,6,10,14,18,22` and verifies (1) against THM-2411's primitive
integer fluxes.

Forgetting the chosen `q`-sheet gives the necessary invariant equations

```text
T N1^2=7496192^2 lambda^2,                              (3)
Z N1^4=7496192^4 lambda^4.                              (4)
```

The script supplies literal ideal-membership identities deriving (3) from
`(2,q^2-T)` and (4) from `(3,T^2-Z)`.

At `lambda=0`, the inherited open branch has `q!=0`, so (2) is equivalent to
`N1=0`.  Together with `N2=0`, this is exactly the flux fibre used by
THM-2411/2636/2692.  This statement concerns the eliminant fibre only:
THM-2692's final Kummer contradiction still uses a genuine nonsplit deck and
must not be transferred to the split deck.

## 3. Why 23 appears

Use THM-2692's now-typed scaling

```text
B_0=rho^2, C_0=c rho^3, D_0=d rho^4,
E_0=e rho^5, W_0=w rho^6,
t=rho/y, v=u/y^2, zeta=Z/y^3.                           (5)
```

The primitive fluxes have weighted degrees

```text
wt(N1)=5,                 wt(N2)=6,                  (6)
```

while `wt(Z)=3`.  Consequently the left side of (4) has weight

```text
3+4*5=23.                                                (7)
```

Writing `N1=y^5 f1` and `N2=y^6 f2`, equation (4) becomes

```text
f2(t,v,zeta)=0,
zeta f1(t,v,zeta)^4=eta t^23,
eta=7496192^4 lambda^4/rho^23 !=0.                      (8)
```

The exact scaled polynomials are

```text
f1=
 9370240c t^3v-232320c t^3+2044416d t^4-14992384e t^5
 -2981440t^2v+819896t^2zeta+24640t^2
 +3689532v^2-1449459v zeta-101640v+83853zeta+252,

f2=
 449771520c t^3v-206145280c t^3zeta-1239040c t^3
 -1978994688d t^4v+16355328d t^4-239878144e t^5
 -1319329792w t^6+1443016960t^2v^2-71554560t^2v
 +65591680t^2zeta+98560t^2-1190488992v^3
 +147581280v^2-162339408v zeta-1219680v
 +15944049zeta^2+2236080zeta+672.
                                                               (9)
```

Their frozen canonical digests begin `54017f31` and `9aecbd91`.

The prime `23` is therefore not numerology.  It is the first-flux weight
after the two successive sign quotients.  Equivalently, the unsquared form
of (2) is the Bring--Jerrard-shaped relation

```text
1331 mathcal A q^5+4 mathcal K q+7496192 lambda=0,     (10)
```

because `N1=1331 mathcal A Z+4 mathcal K` and `Z=q^4`.

## 4. The exact fixed fibre is three new plus five old

At `t=0`, eliminate `zeta` from `f2=0` and
`zeta f1^4=0`.  Resultant multiplicativity gives

```text
Res_zeta(f2,zeta f1^4)
 =Res_zeta(f2,zeta) Res_zeta(f2,f1)^4
 .=G3(v)L5(v)^4.                                         (11)
```

The old factor is THM-2636's squarefree irreducible quintic `L5`.  The new
factor is the squarefree cubic

```text
G3(v)=(121v-1)(14641v^2-1694v+1),                      (12)
```

with roots

```text
1/121,            (7+4sqrt(3))/121, (7-4sqrt(3))/121. (13)
```

The script verifies `gcd(G3,L5)=1`.  Thus the fixed eliminant has degree

```text
3+4*5=23.                                               (14)
```

The multiplicity four on `L5` is geometric information, not a nonreduced
claim about the normalization.

## 5. Exact local normalization indices

At every `L5` point, `zeta` is a unit and the Jacobian

```text
det partial(f1,f2)/partial(v,zeta)
```

is a unit.  Hence `f1` and `t` are local coordinates on the `f2=0`
surface, and (8) has the local form

```text
unit*f1^4=unit*t^23.
```

Since `gcd(4,23)=1`, its normalization has one branch with

```text
ord(t)=4,                    ord(f1)=23.               (15)
```

At every `G3` point, `f1` and `partial_v(f2|_(zeta=0))` are units.  Equation
(8) has local form

```text
zeta=unit*t^23,
```

so there is one branch with `ord(t)=1`.  Therefore the normalized fibre over
`t=0` has exact order word

```text
4,4,4,4,4,1,1,1,
deg div_0(t)=5*4+3*1=23.                                (16)
```

This is the split analogue of the nonsplit five-place picture: five old
branches acquire index four, while three previously invisible `zeta=0`
branches complete the prime-23 fibre.

## 6. Loss ledger

| map | preserved | destroyed | required sidecar |
|---|---|---|---|
| `(q,N1,lambda)->(T=q^2,N1,lambda^2)` | necessary split first-flux equation | chosen split-sheet sign | `q` or the anti-invariant line |
| `(T,N1)->(Z=T^2,N1)` | fourth-power invariant equation and weight 23 | square root `T`; four roots replace the physical two-step tower | `T` plus `q^2=T` |
| `lambda->0` | the exact THM-2692 flux fibre | split constant obstruction | proof that `lambda` vanishes, not an assumption |
| restriction to the even bank | THM-2411's five coefficient parameters and `lambda` | every odd Faber seed allowed on the full split deck | the odd-seed coefficient/local-system bank |
| fixed resultant | branch count and multiplicities | global components, infinity, and the Keller one-form | normalization and pole divisor |

Equation (4) is necessary but not equivalent to the chosen-sheet equation
(2) without the `q,T` sidecars.  It must not be used to manufacture a split
trajectory.

### The inherited Kummer gate is exactly neutralized

There is nevertheless no missing fourth-root extension on the nonzero-
`lambda` curve.  The physical spectral coordinate is

```text
Z=rho^3 zeta/t^3=q^4.
```

Equations (2), (5), and `N1=rho^5 f1/t^5` reconstruct the chosen sheet
rationally:

```text
q=-7496192 lambda t^5/(rho^5 f1).                    (17)
```

Indeed (8) gives the global identity

```text
rho^3 zeta/t^3
 =(7496192 lambda t^5/(rho^5 f1))^4.                 (18)
```

Thus the Kummer cover used in THM-2636/2692 is not merely unramified at the
fixed fibre; it is identically split on this even-Faber `lambda!=0` curve.
The local numbers make the mechanism visible.  At an old branch,
`ord(t)=4`, `ord(zeta)=0`, so `ord(Z)=-12=4(-3)`.  At a new branch,
`ord(t)=1`, `ord(zeta)=23`, so `ord(Z)=20=4(5)`.  Formula (17) gives the
corresponding `q`-orders `-3` and `5`.

This is a sharp stopping certificate: the five-old/three-new prime-23 grammar
does **not** feed the old odd-place Kummer genus gate.  Any genus obstruction
must come from the intrinsic normalized curve (8), its other branch fibres,
or interactions with the omitted odd bank.

## 7. Cheapest order-raising test and next geometry

THM-2206 identifies the proposed uniform way to force both the odd channels
and `lambda` to vanish.  In
`Z_2[C_2]`, with `Delta=sigma-1`,

```text
Delta^j=(-2)^(j-1)Delta.
```

The companion checks this through `j=8`.  The cheapest decisive test is still
THM-2206's degree-six bank: retain `E6,E5,E3,E2,E1`, clear the boundary,
`Phi`, and `Psi` rows integrally, and compute whether eliminating the odd
seeds raises every surviving anti-invariant defect on

```text
I^j/I^(j+1) isomorphic F_2.
```

A failure there stops the integral route.  A success is only a finite
positive control; it does not prove uniform raising.

Independently, (8) gives a concrete geometric route on the even-Faber
subchart.  The next exact probes are:

1. determine the generic components/irreducibility of the prime-23 complete
   intersection;
2. compute the infinity fibre and its local Newton polygons;
3. combine (16) and the infinity indices in Riemann--Hurwitz; and
4. seek intrinsic monodromy or a different cover: (18) proves that the
   physical fourth-root/Kummer reformulation itself is trivial.

The fixed fibre alone does not determine genus.  It supplies exact local
ramification and the degree of `t`, but a positive-genus conclusion still
requires global connectedness and enough total ramification.

## 8. Reproduction

```bash
python3 04-computation/jc2_degree22_split_lambda_prime23_scout.py
python3 -O 04-computation/jc2_degree22_split_lambda_prime23_scout.py
```

Both modes byte-match the frozen output.  Every calculation is over exact
rational polynomial rings; no finite-field or floating-point inference is
used.
