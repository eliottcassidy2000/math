---
id: THM-3470
title: "Mixed-character triangular quartics have a sharp fifth-moment gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the
  reflection-symmetric three-real-parameter triangular family
  g=aM+b(U^2-1/6)+c(U^4-1/15)+iU on the normalized triangle, every a!=0
  gives an explicit polynomial automorphism with constant Jacobian -2a.
  Exactly four real automorphisms annihilate moments 1 through 4; all have
  nonzero support in every C3 character space and all exit at moment 5.
  Over C, the moment-2-through-4 locus has 24 distinct reduced points, while
  moments 2 through 5 generate the unit ideal.  This is a theorem only for
  the displayed shear family, not HFC(3), FC(3), or JC(2).
audit: >
  Exact grevlex gives the unit ideal.  An independent lex route gives a
  squarefree degree-24 terminal polynomial and a degree-23 fifth-moment
  remainder with gcd one.  Independent Sturm replay proves exactly four real
  roots; rational Horner boxes exclude a=0, c=0, and moment 5=0 at each.
  A direct simplex expansion agrees with the mixed-coordinate moment table
  in 231 cases.  The quadratic boundary has six reduced complex prefix
  points, exactly two real Keller points, and a sharp fourth-moment exit.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
related:
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
  - THM-3448-weighted-keller-cyclic-jelonek-inertia-family
  - THM-3465-nonreal-cyclic-character-keller-rigidity-and-hfc-separation
script: 04-computation/factorial_mixed_character_triangular_quartic_moment_gate_thm3470.py
output: 05-knowledge/results/factorial_mixed_character_triangular_quartic_moment_gate_thm3470.out
script_sha256: 614f473562ba36da297947c1420929d5d0d84ca5744df549123c7540e602720c
output_sha256: 195ce0b8200254697563e55487a3094801efa585b80d6bc66fe98228fe44aa86
semantic_sha256: ce17222a7166129cf39b9f6b10c5bad061d06137a1ccb8297c21ec78a42173d0
hash_basis: raw bytes
---

# THM-3470 -- mixed triangular quartics and the fifth-moment gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Family and statement

Let

```text
Delta={(x,y): x>=0, y>=0, x+y<=1},
<h>=2 integral_Delta h(x,y) dx dy,                          (1)
```

so `<1>=1`.  Put

```text
U=x-y,                         M=x+y-2/3,                   (2)
```

and, for real `a,b,c`, define

```text
g=aM+b(U^2-1/6)+c(U^4-1/15)+iU.                            (3)
```

Write `M_r=<g^r>` and `F_g=(Re g,Im g)`.  Then:

1. `M_1=0`, every `M_r` is real, and

   ```text
   Jac(F_g)=-2a.                                            (4)
   ```

   If `a!=0`, `F_g` is an explicitly triangular polynomial automorphism.

2. In `Q[a,b,c]`, after clearing positive integer denominators from
   `M_2,...,M_5` and calling the primitive numerators `P_2,...,P_5`,

   ```text
   (P_2,P_3,P_4,P_5)=(1).                                  (5)
   ```

   Thus no complex parameter, and in particular no real constant-J point,
   annihilates moments `M_1,...,M_5`.  No member of this family is an HFC
   null-moment candidate.

3. The prefix scheme

   ```text
   V(P_2,P_3,P_4) subset A^3_C                              (6)
   ```

   consists of exactly `24` distinct reduced points.  Exactly four are real.
   At every real point `a c!=0`, so it is a genuine degree-four polynomial
   automorphism, and its degree-four part has nonzero support in all three
   `C3` character spaces.  Every one has `M_5!=0`.  Hence four moments are
   sharply insufficient and the fifth is decisive within this family.

4. On the quadratic boundary `c=0`, the locus `M_2=M_3=0` has exactly two
   real polynomial automorphisms.  Both have `M_4!=0`, and

   ```text
   (P_2|_(c=0),P_3|_(c=0),P_4|_(c=0))=(1).                 (7)
   ```

The theorem concerns this reflection-symmetric one-shear cell.  It does not
cover odd powers of `U` in the real part, rotated source frames, arbitrary
mixed degree-at-most-four coefficients, or two-shear Henon cells.

## 2. Exact triangle moment table

Use

```text
x=s(1+u)/2,             y=s(1-u)/2,
0<=s<=1,                -1<=u<=1.                           (8)
```

The normalized measure is `s ds du`, while `M=s-2/3` and `U=su`.  Therefore

```text
<M^j U^k>=0                                      if k is odd,

<M^j U^k>
 =sum_(h=0)^j C(j,h)(-2/3)^(j-h)
    2/((h+k+2)(k+1))                            if k is even. (9)
```

In particular `<M>=0`, `<U^2>=1/6`, and `<U^4>=1/15`, proving `M_1=0`.
Reflection `x<->y` fixes `M`, sends `U` to `-U`, and sends `g` to its complex
conjugate.  The measure is invariant, so every `M_r` is real.

The companion independently expands standard simplex monomials

```text
<x^i y^j>=2 i! j!/(i+j+2)!                                 (10)
```

and compares them with (9) at all `231` pairs of total degree at most `20`.

The first two primitive moment numerators are

```text
P_2=350a^2+280ab+160ac+245b^2+310bc+112c^2-1050,
M_2=P_2/6300,                                               (11)

P_3=-100100a^3-75075a^2b-10725a^2c+214500ab^2
    +351780abc+149760ac^2-900900a+157300b^3
    +388245b^2c+331110bc^2-1576575b+96448c^3-997425c,
M_3=P_3/13513500.                                           (12)
```

The exact, larger `P_4` and `P_5` are frozen in the declared transcript.

## 3. Constant Jacobian and explicit inverse

Put `P=Re g` and `Q=Im g=U`.  Direct differentiation in `x,y` gives (4).
When `a!=0`, recover

```text
U=Q,
M=(P-b(Q^2-1/6)-c(Q^4-1/15))/a,
x=(M+Q+2/3)/2,
y=(M-Q+2/3)/2.                                             (13)
```

These are polynomial formulas.  After the affine change `(x,y)<->(M,U)`,
the map is a triangular shear followed by a nonzero scaling, so it is a tame
polynomial automorphism.  Thus the prefix witnesses below are benign
automorphisms, not putative Jacobian counterexamples.

## 4. Two exact elimination paths

The first route computes a grevlex Groebner basis over `Q` and obtains

```text
Groebner(P_2,P_3,P_4,P_5)=[1],                              (14)
```

which proves (5).

The second route first triangularizes the prefix ideal.  After normalizing
nonzero rational coefficients, its lex basis has the shape

```text
a-A(c),                 b-B(c),                 C_24(c),   (15)
```

where `deg A=deg B=23`, `C_24` is even and squarefree, and
`deg C_24=24`.  Hence the quotient is isomorphic to
`Q[c]/(C_24)`: (6) is a reduced graph of length `24`, not merely a vector-space
dimension count.

Reduce `P_5` by (15).  Its normal remainder `R_5(c)` has degree `23` and

```text
gcd(C_24,R_5)=1.                                           (16)
```

This independently proves that moment five misses every prefix point.  The
transcript freezes the exact polynomials through the hashes

```text
lex basis:       8c27de52c5378cec1375f297b8ae44e3db5f64640948844765495d824a4bd397,
C_24:            2cafb60b9a45045f7b9c7b36816e3474c4a0ff413980ca0667dcb8776910f545.
                                                                    (17)
```

## 5. Exact real locus and character support

Sturm isolation gives exactly four real roots of `C_24`, one in each rational
box printed by the companion.  Exact Horner interval evaluation of `A(c)`,
`B(c)`, and `R_5(c)` gives the approximate centers

```text
(a,b,c; M_5)
 ( 1.835799943,  0.157867847,-2.1656644225; -0.026453454),
 (-1.882953388,  3.444298723,-3.9376205842;  0.002401252),
 (-1.835799943, -0.157867847, 2.1656644225;  0.026453454),
 ( 1.882953388, -3.444298723, 3.9376205842; -0.002401252).  (18)
```

The rational boxes, not these decimals, are the proof objects.  They exclude
`a=0`, `c=0`, and `M_5=0` at all four points.

Let `omega^3=1`, `omega!=1`, and use THM-3310's conjugate Fourier coordinates
`z,w`.  Then

```text
U=((1-omega^2)z+(1-omega)w)/3.                              (19)
```

Both linear coefficients are nonzero.  All five terms of `U^4` therefore
occur, and their weights modulo three meet `{0,1,2}`.  Since `c!=0`, the
degree-four part has nonzero projection to every `C3` character space.  This
is character support, not `C3` equivariance, monodromy, or inertia.  It is
precisely why THM-3465's pure-character rigidity does not apply.

Coefficient conjugation in the real `(M,U)` coordinates is the Hermitian
dagger here; in Fourier coordinates it is coefficient conjugation together
with `z<->w`.

## 6. Quadratic boundary hostile

At `c=0`, the lex terminal polynomial for `P_2,P_3` is

```text
64095 b^6-999792 b^4+5752320 b^2-6272000.                  (20)
```

It is squarefree and has exactly two real roots.  If `beta` is the unique
root in

```text
(118352913031470/10^14,118352913031480/10^14)
```

and

```text
alpha=-(20959065 beta^5-237198984 beta^3+1147229440 beta)
       /513990400,                                         (21)
```

then

```text
g=alpha M+beta(U^2-1/6)+iU                                (22)
```

has constant nonzero Jacobian and exactly `M_1=M_2=M_3=0`, while

```text
0.060663464<M_4<0.060663465.                               (23)
```

The other real witness is `(-alpha,-beta)` and has the same positive fourth
moment.  Adding `P_4` produces the unit ideal, proving (7).  This is the sharp
lower-degree hostile to any proposed three-moment contradiction.

## 7. Information contract and reproduction

```text
source:      the full Hermitian degree-at-most-four Keller/moment system
target:      the displayed reflection-even one-shear slice
map:         coefficient restriction to a,b,c
preserved:   real dagger, exact Jacobian, degree, moments, explicit inverse
destroyed:   odd shear direction, rotated frames, two-shear coefficients
sidecar:     full C3 coefficient chart and Keller normal form
next test:   adjoin one odd shear or rotate U before generalizing
```

Run

```text
python3 04-computation/factorial_mixed_character_triangular_quartic_moment_gate_thm3470.py
python3 -O 04-computation/factorial_mixed_character_triangular_quartic_moment_gate_thm3470.py
```

and compare raw bytes with the declared output.  Both routes are
byte-identical.  The script uses exact rational arithmetic, RuntimeError truth
gates, two moment-table derivations, two elimination paths, Sturm counts, and
rational interval arithmetic.  It contains no floating-point or random
truth gate.

**QED.**
