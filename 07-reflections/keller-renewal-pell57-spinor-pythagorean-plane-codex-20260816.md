# The Keller renewal ray is multiplication by a norm-minus-eight Pell element

**Status: PROVED ALGEBRAIC COROLLARY OF THE THM-3522 PACKET TRANSFORM;
VERIFIED-EXACT.**  The abstract renewal orbit

```text
(e_(n+1),m_(n+1))=(7e_n-2m_n,3e_n-2m_n),
(e_0,m_0)=(1,0),
```

is exactly a multiplicative orbit in the quadratic integer ring of
`Q(sqrt(57))`.  The named fixed Keller polynomials `L,H,J,G,R5,R6,R7,R8`
realize rows `n=0,...,7`, and THM-3528 now realizes every later row as a raw
polynomial complete packet.  Later finite-sheet units, prime/image status,
separability, and every general Jacobian-conjecture statement remain open.

## Inheritance pass

THM-3522 proves the packet transform

```text
A(e,m) -> A(7e-2m,3e-2m)
```

whenever the next fixed-chart norm is polynomial.  THM-3523 pays that gate
through `R8` after THM-3527.  THM-3509 identifies `(e,m)` with a reduced fraction,
harmonic `K4` face, primitive Pythagorean current, and an arithmetic self-map
of the even Berggren tree.  Its canonical hostile shows that the renewal
matrix is not a fixed Berggren branch word.

The least-used sidecar is the quadratic field singled out by the renewal
matrix's characteristic polynomial

```text
t^2-5t-8.
```

Its discriminant is `57`, the same radical governing the attracting packet
slope.

## Exact Pell semiconjugacy

Put

```text
M = [[7,-2],[3,-2]],
X=6e-9m,
Y=m.
```

A direct matrix identity gives

```text
2 [X']   [5 57][X]
  [Y'] = [1  5][Y].                                  (1)
```

Since `57=1 mod 4`,

```text
alpha=(5+sqrt(57))/2
```

is a quadratic integer, and (1) is precisely

```text
X'+Y'sqrt(57)=alpha(X+Ysqrt(57)).                     (2)
```

The initial value is `X_0+Y_0sqrt(57)=6`.  Therefore

```text
X_n+Y_nsqrt(57)=6 alpha^n.                            (3)
```

Because `Norm(alpha)=-8`, taking norms gives the all-`n` identity

```text
X_n^2-57Y_n^2=36(-8)^n.                              (4)
```

Returning to packet coordinates,

```text
3e_n^2-9e_nm_n+2m_n^2=3(-8)^n.                       (5)
```

Thus the packet ray lies on a tower of Pell-type conics whose level changes
by the norm grade `-8`.  This is not an approximation or a finite pattern.

The conjugate roots

```text
alpha=(5+sqrt(57))/2,
beta =(5-sqrt(57))/2
```

also give

```text
m_n = 3(alpha^n-beta^n)/sqrt(57),
X_n = 3(alpha^n+beta^n),
u_(n+2)=5u_(n+1)+8u_n
```

for either scalar coordinate `u=e,m`.  This is a genuine Fibonacci-style
linear recurrence, but it is not the ordinary Fibonacci sequence.

## Cassini and primitivity

Since `det M=-8`, consecutive packet columns obey

```text
e_nm_(n+1)-m_ne_(n+1)=3(-8)^n.                       (6)
```

Equation (6) is the generalized Cassini law.  For every `n>=1`, induction
modulo six gives

```text
(e_n,m_n)=(1,3) mod 6.
```

Any common divisor divides a preceding Cassini determinant, hence has no odd
prime outside `{3}`; the congruence excludes `2` and `3`.  Therefore

```text
gcd(e_n,m_n)=1,
```

and both coordinates are odd.  The fractions `m_n/e_n` are reduced and
converge alternately to

```text
(9-sqrt(57))/4.
```

Their consecutive determinants have magnitude `3*8^n`, so they are not
Farey neighbors and do not form a Stern--Brocot edge ray.

## Moving planes of primitive Pythagorean triples

For `n>=1`, define the parity-divided Euclid triple

```text
A_n=(e_n^2-m_n^2)/2,
B_n=e_nm_n,
C_n=(e_n^2+m_n^2)/2.                                 (7)
```

The mod-six and gcd facts make (7) a primitive Pythagorean triple.  Equation
(5) becomes the exact moving-plane law

```text
A_n-9B_n+5C_n=3(-8)^n.                               (8)
```

The first fixed-tower triples are

```text
H : (20,21,29),
J : (812,645,1037),
G : (31820,26829,41621),
R5: (1254188,1044885,1632413),
R6: (49372940,41233821,64326629),
R7: (1944120812,1622829285,2532425837),
R8: (76548298700,63904114029,99716487221).
```

The renewal step acts directly on triples by

```text
T = [[20,-8,20],
     [17,-20,25],
     [25,-20,33]].                                   (9)
```

For the Lorentz form `J=diag(-1,-1,1)`, exact multiplication gives

```text
T^T J T = 64 J,
det T=-512.                                           (10)
```

Thus `T` is an integral Lorentz similitude of multiplier `64`, not a
unimodular Berggren isometry.  It preserves the light cone and carries this
special primitive orbit to itself, but it is not one edge or one fixed word
in Berggren's ternary tree.  This sharpens the earlier determinant-square
class obstruction: the recurrence is a nonunimodular arithmetic jump across
the same fraction universe.

## Harmonic `K4`, tournaments, and the missing ancestry

Write `x=e-m`, `y=m`.  THM-3509's harmonic face is

```text
(u,v,z)=(xy,x(x+y),y(x+y)),
1/u=1/v+1/z.                                          (11)
```

The exact gcd decoder is

```text
(gcd(u,v),gcd(u,z),gcd(v,z))=(x,y,x+y),
```

so the Pell packet, reduced fraction, harmonic face, and primitive triple
remain mutually typed.  Numeric comparison gives strict transitive `T4` and
`T6` carriers on these nonroot rows, but those tournaments record order only;
they do not encode the nonunimodular step (9).  The ternary Berggren ancestry
word remains a separate address.

This also clarifies the harmonic-series analogy.  Equation (11) is a genuine
three-term harmonic subseries identity, while the recurrence lives on its
primitive decoder `(x,y)`.  Replacing the decoder by the scalar reciprocal
sum loses the Pell level, Cassini sign, and tree ancestry.

## Two gradings, not one

The fixed inverse tower has multiplicative fibre degree `3^n`.  The renewal
packet has quadratic norm grade `(-8)^n`.  Equations (2)--(4) make the latter
a literal multiplicative monoid generated by `alpha`, while composition of
maps supplies the former.

```text
map-composition grade: 3^n,
packet/Pell grade:     Norm(alpha^n)=(-8)^n.
```

Their shared index does not identify the gradings.  THM-3528 now proves that
every raw cleared norm is polynomial; it does not prove the stronger finite-
sheet unit needed for `L`-coprimality or a new image prime.

THM-3527 now realizes the formerly predicted row

```text
n=7: (419839,152211) = R8.
```

THM-3528 makes the next two rows actual raw polynomial packets:

```text
n=8: (2634451,955095),
n=9: (16530967,5993163).
```

They are not proved irreducible, `L`-coprime, or image equations.

## Connection contract

| field | exact answer |
|---|---|
| source | THM-3522 packet state `(e,m)` |
| target | Pell coordinate `X+m sqrt(57)` and primitive Pythagorean current |
| map | `(e,m)->(6e-9m)+m sqrt(57)` |
| preserved | recurrence index, primitivity, reduced fraction, harmonic K4 decoder, light-cone nullity |
| transformed | quadratic norm multiplies by `-8`; Lorentz form by `64` |
| destroyed by triple alone | packet sign/level unless the moving-plane value is retained; Berggren ancestry word |
| hostile | wrong `-2m` coefficient breaks (4); `det T=-512` blocks Berggren-isometry promotion |
| fixed realization | raw polynomial packets at all levels by THM-3528; named/L-coprime through `R8` |
| open boundary | next finite sheet, factors, images, irreducibility, separability, general JC |

## Reproduction

```text
python -B 04-computation/keller_renewal_pell_spinor_sidecar_20260816.py
python -B -O 04-computation/keller_renewal_pell_spinor_sidecar_20260816.py
```

The semantic digest is
`c5758a441d60c45f254edbfea3d2ec06c34886fdcd0fc9bc27685dc0aaa4a5af`.
