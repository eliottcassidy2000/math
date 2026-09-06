# Hamiltonian transport, source-form preservation, and two stopping obstructions

**Status: PROVED symbolic identities and analytic argument;
[independent audit PASS](planar_jc48_sep06_hamiltonian_audit.md).**
September 6, 2026. The compensated weight22 transport has a
Hamiltonian interpretation through row15. Its exact displayed representative
has a generator in `P_1`; a terminal change gives a five-term generator in
`P_0=k[p,y]`. This generator neither defines a polynomial additive-group
action nor preserves the fixed source form to all orders. A separate exact
criterion identifies all polynomial generators that preserve that source
form infinitesimally for a fixed source polynomial or for **every** source
polynomial. Every nonconstant generator in either class fails local
nilpotence.

All assertions are over a characteristic-zero field. No full `B_2` Keller
pair, polynomial termination, or planar Jacobian conjecture conclusion is
asserted.

## 1. Inheritance and the connection being tested

The source is the already independently audited
[weight22 finite transport](planar_jc48_sep06_weight14.md), on the fixed
low jet of [THM-4308](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md).
Its inherited row15 boundary remains
[THM-4438](../../01-canon/theorems/THM-4438-jc2-row-fifteen-relative-response-on-boundary-gm.md).
The conic `Q=0` clears projected depth there; the remaining no-response
obstruction is its bracket quartic. Those obligations are distinct.

The proposed connection sends the finite correction of `(A,C)` to a
Hamiltonian derivation in the original `(x,t)` source chart. It preserves
the bracket formally. It does not automatically preserve polynomial
termination, the fixed affine source form, or depth. These are three separate
side conditions. The five-concept board is: finite symplectic jets; terminal
tangent freedom; exact depth representatives; divisorial denominators; and
local nilpotence.

The cheapest hostile is the tempting leading generator `p³y⁵/15`, which
already misses the next generator row. The corrected near miss is to regard
a Hamiltonian polynomial as giving a polynomial flow. The least-used
sidecar is the terminal tangent kernel, which lets us move the generator
from `P_1` into `P_0` without changing the row15 source response. Targeted
searches in the source-normal and Hamiltonian canon and JC2 result notes
found no prior statement of the exact generator or carrier criterion below;
no external priority claim is made.

Use the fixed source coordinates and sign convention

```text
p=t(1+x²t),       y=xtp,       u=x²t,
{f,g}=f_x g_t-f_t g_x,
G=-u/2+H(p,y),    P(A,C)=C²-A³+(3/4)A+1/4.
```

For depth, use exactly THM-4308's span
`P_d=span{x^a u^b p^c y^e : a+b<=d}`. A first generator term in row13
produces a coordinate correction beginning in row12. Generator rows and
unknown rows must therefore not be identified.

## 2. Exact old generator, then a terminal repair into `P_0`

Define the five-term polynomial

```text
S0 = p³y⁵/15 +16p⁴y⁵/135 -y⁷/15
     -2848p⁵y⁵/6075 +2py⁷/621.                         (1)
```

It lies in `P_0`; its largest residual weight for weights `(2,3)` on `(p,y)`
is 25. This is the generator's weight, not the weight22 bound for the source
packet it pays.

First, the **exact previously displayed** correction `(dA,dC)` has generator

```text
Sstar = S0 -4606216p¹⁴y/83835 +2303108p¹⁰y³/83835
        -64798p⁶y⁵/5589 +1367399p²y⁷/139725
        -52451xy⁸/3105.                                (2)
```

Thus `Sstar in P_1`. Its five added terms first appear in generator row16.
For every background with the declared fixed coefficients in rows zero through three,

```text
{A,Sstar}=dA mod t¹⁶,       {C,Sstar}=dC mod t¹⁶.       (3)
```

Here and below only the four fixed low rows enter. Later background rows
begin at four, whereas `S_t` and `S_x` begin at twelve and thirteen.

For the smaller generator set

```text
dA0=pi15({A,S0}),           dC0=pi15({C,S0}).
```

These agree with `dA,dC` through row14. Put

```text
R(x)=7080885x⁸-4102197x⁶+4859850x⁴-11515540x²+23031080.
```

Their only difference through row15 is

```text
dA0-dA = (8x²R(x)/419175)t¹⁵,
dC0-dC = -(2x(x²+2)R(x)/139725)t¹⁵.                   (4)
```

The leading pair is a tangent multiple of `(A0_x,C0_x)`. Direct substitution
shows that it changes neither the source equation through row15 nor the
bracket through row14.

Depth is checked independently of that tangent observation. The source
contains five literal `P_2` monomials and five literal `P_3` monomials for
the two polynomials in (4), each of valuation exactly fifteen. For example,
the leading `x¹⁰t¹⁵` in the first difference is represented by
`u² p y⁶`; the leading `x¹¹t¹⁵` in the second is represented by `u³p²y⁵`.
Combining these with the previous transport's explicit depth representatives
proves

```text
pi15(dA0) in pi15(P_2),     pi15(dC0) in pi15(P_3).     (5)
```

Both corrections satisfy the inherited degree caps. For the same packet

```text
T = p⁵y⁴ -(7/10)py⁶ -(508/135)p²y⁶
    +(235202/27945)p³y⁶,
```

the complete coefficient identities are

```text
{A,dC0}+{dA0,C}=0 mod t¹⁵,
2C dC0-3A²dA0+(3/4)dA0=T mod t¹⁶.                    (6)
```

Hence `A+w dA0,C+w dC0,H+wT` has exactly the same finite-row and depth
properties as the old transport for every scalar `w`. Nonlinear terms start
at rows23 in the bracket and24 in the source equation. Choosing the same
`w=-j/beta` as in the prior note still removes the old weight24 payer and
leaves the weight22 packet, now with the terminal representative (4).

The formal derivation `L={-,S0}` raises `t`-adic valuation by at least twelve.
Its exponential is therefore a well-defined formal symplectic source change,
and `exp(wL)(A)=A+w dA0 mod t¹⁶`, similarly for `C`. This is a formal
integration statement in the completed source ring, not a polynomial
automorphism assertion.

### Two exact generator hostiles

The bare leading generator fails because

```text
[t¹⁴](S0-p³y⁵/15)=x⁵(16-9x²)/135 != 0.               (7)
```

Also, the **unchanged old** correction cannot have a `P_0` generator through
the needed order. Its Hamiltonian jet is unique modulo `t¹⁷` and a scalar:
the low Jacobian is a unit, so two generators of the same coordinate
correction have both derivatives zero to the requisite orders. Expand a
putative `P_0` polynomial by increasing valuation. The leading term of
`p^c y^e` is `x^e t^(c+2e)`, and different `e` at a fixed valuation are
linearly independent. Matching rows13 through16 recursively leaves exactly

```text
-(52451/3105)x⁹t¹⁶.                                    (8)
```

A monomial first appearing in row16 has `e<=8`, so (8) cannot be supplied
in `P_0`. Formula (2) supplies it in `P_1`. The terminal repair (4), rather
than an unjustified depth identification, is what permits (1).

## 3. Correct Poisson denominators and the source-form obstruction

Put `D=p³-y²`. This symbol is unrelated to the fixed low-jet parameter
`Delta=896/15`. The source rational functions satisfy

```text
t=D/p²,       x=yp/D,       u=y²/D,
{u,p}=2xp=2yp²/D,       {u,y}=3up,
{p,y}=-tp=-D/p.                                        (9)
```

A draft probe mistakenly called `2xt(1+u)` by the name `2y`; it is `2xp`.
The root audit caught this before any frozen statement. The certificate
checks all three identities in (9) by literal `(x,t)` differentiation and
retains the false `2y` identity as a rejected control. The finite-jet
computations always used literal derivatives and were unaffected.

For any `S,H in k[p,y]`, let `E=2p partial_p+3y partial_y`. Then

```text
{ -u/2+H,S }
  = -py E(S)/(2D) -(D/p)(H_p S_y-H_y S_p).             (10)
```

For the particular generator (1), reduction in the cusp ring gives

```text
E(S0) mod D = p¹⁰y(14/5-2848p/243) != 0.              (11)
```

The second term of (10) is regular at the generic point of `D=0`, since
`p` is a unit there. The first term has a nonzero pole by (11). Thus for
**every** polynomial `H`, the first variation fails to lie in `k[p,y]`.
The formal flow of this generator does not retain the fixed source form
`-u/2+polynomial(p,y)` to all orders. Its successful finite source response
is exactly (6), not an all-order closure.

## 4. Exact fixed-source and universal-source repair criteria

For a **fixed** polynomial `H`, the exact criterion is

```text
{-u/2+H,S} in k[p,y]
   iff S=c+D R and p divides J_(p,y)(H,S).             (12a)
```

Indeed, in (10) the second term is regular at the generic point of `D=0`,
so the first forces `D|E(S)`. The cusp grading argument below says that
this is equivalent to `S=c+D R`. After this substitution, the first term
is polynomial. Since `p` is coprime to `D`, the second term is polynomial
exactly when the divisibility in (12a) holds. Thus the condition is both
necessary and sufficient. The positive control `H=0,S=D` gives
`{-u/2,D}=-3py`; it lies outside the universal carrier below.

For an arbitrary polynomial generator `S in k[p,y]`, one has the equivalence

```text
{ -u/2+H,S } in k[p,y] for every H in k[p,y]
        iff
S in k + p²(p³-y²)k[p,y].                              (12)
```

This is a theorem about the entire affine space of source polynomials.
It is stronger than the fixed-`H` criterion (12a), as the positive control
shows.

**Necessity.** Subtract the conditions at `H=0,p,y`. Formula (9) forces
`p|S_p` and `p|S_y`. The second condition makes `S(0,y)` constant; the
first kills the coefficient of `p¹`. Thus `S=c+p²R`.

The condition at `H=0`, or just regularity at the generic point of `D=0`,
forces `D|E(S)`. The derivation `E` descends to
`k[p,y]/(p³-y²)=k[z²,z³]` and acts as multiplication by each positive
weighted degree on its corresponding homogeneous component. In
characteristic zero its kernel is exactly the constants. Consequently
`S=c'+D R'`. Evaluating at `(0,0)` gives `c=c'`. Since `p²` and `D` are
coprime in `k[p,y]`, their principal ideals intersect in their product.
This proves the right side of (12).

**Sufficiency.** If `S=c+p²D R`, then both `S_p,S_y` are divisible by `p`,
and

```text
E(S)=p²D(10R+2pR_p+3yR_y).
```

Every denominator in (10) therefore cancels. More explicitly the image is

```text
-p³y(10R+2pR_p+3yR_y)/2
-D[H_p(-2pyR+pD R_y)-H_y((2D+3p³)R+pD R_p)],           (13)
```

which is a polynomial. This proves (12). Such a generator preserves
`k[p,y]` under its derivation and maps `-u/2` into that ring, so the same
source form is preserved coefficient by coefficient under the formal
exponential in the time parameter. No projected-depth or polynomial-flow
conclusion follows from (12).

These give concrete repaired generator spaces for later finite or formal
investigation. Paying a prescribed response and verifying actual depth
remain separate obligations. The following argument disposes of polynomial
additive-group flows for the entire fixed-source carrier.

### No nonconstant fixed-source carrier is locally nilpotent

Suppose a nonconstant `S=c+D R` belongs to (12a), or more generally just
has this displayed divisibility, and pull it back to `k[x,t]`. If
`L={-,S}` were locally nilpotent, its kernel would be factorially closed.
For completeness, define `deg_L f` as the largest `n` with `L^n f!=0` for
nonzero `f`. The Leibniz formula in a characteristic-zero domain gives
`deg_L(fg)=deg_L f+deg_L g`: at the sum of the two degrees the only surviving
term is a nonzero binomial coefficient times the two highest derivatives.
Thus an invariant nonzero product has every factor invariant.

But `S-c` is invariant and

```text
S-c=t³(1+x²t)² R(p,y).
```

Factorial closure would make both `t` and `1+x²t` invariant. Their difference
from the constant one then gives invariant `x²t`; factorial closure gives
invariant `x` as well. Therefore `L` annihilates both polynomial source
coordinates and is zero. A zero Hamiltonian derivation in characteristic
zero means that the pulled-back `S` is constant, contrary to the hypothesis.
The substitution in `(p,y)` is injective by the rational inverse (9).

Consequently **every nonconstant fixed-source carrier**, and in particular
every nonconstant universal-source carrier, fails local nilpotence. Such
generators cannot give polynomial additive-group actions on the source.
They can still provide finite-jet or completed formal repairs; this argument
does not rule out polynomial maps produced by a different mechanism.

## 5. Nontermination of the displayed generator

The unique highest-total-degree term of (1), after substitution in `(x,t)`,
is

```text
c(xt)²⁵,       c=-2848/6075.
```

All other terms have smaller total degree. For a monomial `x^a t^b`, its
bracket with this leading term is
`25c(a-b)x^(a+24)t^(b+24)`. Starting from `x`, the exponent difference
remains one. Induction on the total-degree leading part gives, for every
integer `n>=0`,

```text
leading term of L^n(x)
   =(-2848/243)^n x^(24n+1)t^(24n).                    (14)
```

It is nonzero. Hence `L` is not locally nilpotent and `exp(wL)(x)` does
not terminate as a polynomial in the time parameter `w`. In particular
this derivation cannot be the infinitesimal generator of an algebraic
additive-group action on the polynomial source ring.

This excludes the polynomial-family interpretation of this specific
generator. It does not exclude other generators with the same finite jet,
does not assert anything about a single exceptional specialization of an
infinite formal flow, and does not obstruct the finite additive translations
already proved in the row15 quotient.

## 6. Reproduction and frozen universe

The source checks all relevant polynomial coefficients, not sampled
parameters. It reads only the literal old `WA,WC` depth tables through
`ast.literal_eval`, with the previous source's full SHA256 pinned; it does
not import or execute that producer. The source independently reconstructs
the low jet, both Hamiltonians, every new terminal representative, the
rational Poisson identities, the weighted Euler residue, and the repaired
universal polynomial image.

```sh
python3 -B 04-computation/planar_jc48_sep06_hamiltonian.py
python3 -B -O 04-computation/planar_jc48_sep06_hamiltonian.py
```

Normal and optimized output agree byte for byte: **241 always-active gates**.
The all-degree leading-term induction, carrier equivalences, and general
factorial-closure obstruction are proved above; the finite controls
corroborate their exact formulas.

```text
source SHA256 ba65f911e64907cd41d7a950263e281a8da370045c9da40668796f1245a8df6b
output SHA256 842d3d0b6c8300900c44cb34317dde4f7f64adac72ac18067cb8eafc736952a3
semantic SHA256 270f7657c66819230a5990639875d1faedb15e453b8ac6b22e22d6220120ea57
prior depth source SHA256 a4a140ab538620e7885d8b77758c02e16a5c3a8de7e53e3daac51f01bda321f7
```

The output is [planar_jc48_sep06_hamiltonian.out](planar_jc48_sep06_hamiltonian.out).
The independent referee also reconstructed the finite responses by Fraction
convolution and checked both carrier equivalences and the factorial-closure
argument. The proof, source, output, and audit are frozen.
Root owns canonical integration. No theorem ID is reserved here.
