---
id: THM-2433
title: "Operation-fibre deletion incidence and the startup scar"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every
  commutative operation, deleting a finite hole set changes its weak
  and strict unordered fibres by one universal Burnside formula. For
  the strict summand closure S=N minus {1,4,6}, the additive loss is
  eventually the constant three, whereas the multiplicative loss
  propagates along the principal ideals 4N and 6N with exact overlap
  and square corrections. The internal addition and multiplication
  shadows have explicitly classified transitive closures, while the
  induced divisibility poset has exactly four composite-labelled cover
  edges. Relative multiplicative atoms are primes plus {8,12} in the
  weak fibre and primes, prime squares, plus {8,12} in the strict
  fibre; 12 is the unique retained A014574 centre made artificially
  atomic by the startup deletion.
source: mac-mini-2026-07-26-operation-fibre-startup-scar
depends_on:
  - THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry
related:
  - THM-362-natural-operation-graph-shadows
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
script: 04-computation/operation_fibre_deletion_startup_scar_thm2433.py
output: 05-knowledge/results/operation_fibre_deletion_startup_scar_thm2433.out
script_sha256: 6ca2ba6bfb7bd613f6e50f5043a95148723a3e048e03eea2c72074271f49ddda
output_sha256: 7c81187d3f08e4e1efb38e70b476cf7cbada922a4f9de1318eeb7b5e1e7eeb74
hash_basis: working-tree bytes (LF)
---

# THM-2433 -- finite startup holes leave an infinite operation scar

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2422 separates an existential operation shadow from the labelled
parent fibre and proves

```text
S=Cl^<({2,3})=N_(>0) minus {1,4,6}.                    (1)
```

The finite missing module has very different descendants under the two
operations:

```text
addition:
  three translated rays overlap only during startup
  -> a constant loss of three parent pairs;

multiplication:
  the holes generate the infinite ideals 4N and 6N
  -> divisor-sensitive losses, one overlap, and two square walls. (2)
```

This theorem makes that propagation exact. It also distinguishes three
relations that are easily conflated:

```text
the internal operation shadow retaining a legal co-parent,
its positive transitive closure,
the induced order obtained after forgetting the co-parent.       (3)
```

## 1. Universal swap-fixed deletion calculus

Let `odot` be a commutative operation on a countable carrier `U`, with
every target fibre locally finite, and let `H subset U` be finite. For
functions on `U`, define the ordered fibre convolution

```text
(f star_odot g)(z)
  =sum_(a odot b=z) f(a)g(b),                              (4)
```

and the swap-fixed operator

```text
psi_(2,odot)f(z)
  =sum_(a odot a=z) f(a).                                  (5)
```

Put `S=U minus H`, and write `u,h,s` for the three indicators. Let
`F_T^<=(z)` and `F_T^<(z)` count, respectively,
unordered weak and unordered distinct pairs from `T` landing at `z`.
Then

```text
F_T^<=(z)
  =[(t star_odot t)(z)+psi_(2,odot)t(z)]/2,

F_T^<(z)
  =[(t star_odot t)(z)-psi_(2,odot)t(z)]/2.             (6)
```

Consequently the exact fibre losses caused by deleting `H` are

```text
F_U^<=-F_S^<=
 =u star_odot h-(h star_odot h)/2+psi_(2,odot)h/2,

F_U^< -F_S^<
 =u star_odot h-(h star_odot h)/2-psi_(2,odot)h/2.      (7)
```

**Proof.** Transposition acts on the ordered parent fibre. Burnside's
lemma gives (6): the fixed points are exactly the diagonal pairs. Since
the operation is commutative and `s=u-h`,

```text
u star u-s star s=2u star h-h star h,

psi_2 u-psi_2 s=psi_2 h.
```

Substitution in (6) proves (7). QED.

The diagonal term is not cosmetic. Under addition it lies on even
targets; under multiplication it lies on perfect squares.

## 2. Additive propagation from `{1,4,6}`

Take positive-integer addition and

```text
H={1,4,6},                 S=N_(>0) minus H.             (8)
```

Put

```text
U_+(X)=X/(1-X),             H_+(X)=X+X^4+X^6.
```

Ordinary generating functions turn (7) into

```text
D_+^<=(X)
 =U_+(X)H_+(X)-H_+(X)^2/2+H_+(X^2)/2,

D_+^<(X)
 =U_+(X)H_+(X)-H_+(X)^2/2-H_+(X^2)/2.                (9)
```

Here `D` is ambient fibre size minus `S`-fibre size. For every
`z>=13`, the lost weak and strict pairs are exactly

```text
{1,z-1},             {4,z-4},             {6,z-6}.      (10)
```

Thus

```text
D_+^<=(z)=D_+^<(z)=3,                  z>=13.            (11)
```

The cutoff is sharp for equality of the two tails: at `z=12`, the
diagonal pair `{6,6}` appears in the weak loss but not in the strict
loss. The finite startup irregularity therefore propagates as a
constant labelled-fibre deficit, not as missing later vertices.

## 3. Multiplicative propagation along principal ideals

For proper factor fibres use

```text
U_x={n>=2},                  H_x={4,6},

S_x=U_x minus H_x.                                      (12)
```

Write

```text
K(s)=4^(-s)+6^(-s).
```

Dirichlet series turn multiplicative fibre convolution into ordinary
multiplication. Equation (7) gives, for `Re(s)>1`,

```text
sum_(z>=2) D_x^<=(z)z^(-s)
 =(\zeta(s)-1)K(s)-24^(-s),                            (13)

sum_(z>=2) D_x^<(z)z^(-s)
 =(\zeta(s)-1)K(s)-16^(-s)-24^(-s)-36^(-s).           (14)
```

Coefficientwise,

```text
D_x^<=(z)
 =[4|z and z>=8]+[6|z and z>=12]-[z=24],              (15)

D_x^<(z)
 =D_x^<=(z)-[z=16]-[z=36].                            (16)
```

Thus:

- `24=4*6` is the overlap correction;
- `16=4^2` and `36=6^2` are the swap-fixed corrections; and
- unlike the additive tail, the scar occupies two infinite
  multiplicative ideals.

This is the precise sense in which multiplication sits on a
nonuniform base: deleting one factor affects every multiple of that
factor, and the divisor lattice records the overlaps.

## 4. Internal operation shadows are not transitive

On the vertex set `S` from (8), define

```text
x R_+ z
  iff x<z and z-x in S,

x R_x z
  iff x<z, x|z, and z/x in S.                           (17)
```

These are internal cospan shadows: the second parent must itself be a
retained vertex.

For a composable two-step additive chain with witnesses `a,b in S`,
the shortcut fails exactly when `a+b` is a hole. The complete witness
list is

```text
(a,b)=(2,2),                  (3,3).                   (18)
```

For multiplication, failure means `ab` is a hole, giving exactly

```text
(a,b)=(2,2), (2,3), (3,2).                             (19)
```

Minimal controls are

```text
3 R_+ 5 R_+ 7,                    3 not R_+ 7,

5 R_x 10 R_x 20,                  5 not R_x 20.        (20)
```

The visual similarity with transitivity is therefore real but
conditional: composition multiplies or adds the *witnesses*, and a
missing composite witness destroys the shortcut.

## 5. Exact positive transitive closures

Let `TC` denote positive transitive closure. Then

```text
TC(R_+)
 ={(x,z) in S^2:z-x>=2}.                                (21)
```

Indeed a difference in `S` gives a direct edge. Difference four is
`2+2`. Difference six is `3+3`, except at the only obstructed
intermediate, where `3->5->7->9` uses `2+2+2`. Difference one cannot
be written as a positive sum of legal witnesses, all of which are at
least two.

For multiplication,

```text
TC(R_x)
 ={(x,z) in S^2:x|z, x<z}
   minus {(2,8),(2,12),(3,12)}.                        (22)
```

A quotient in `S` gives a direct edge. Quotient four factors only as
`2*2`; its intermediate is missing precisely for bases two and three.
Quotient six factors as `2*3` or `3*2`; both possible intermediates
are missing only at base two. These give the three exceptions in
(22), and no other quotient is absent.

## 6. Forgetting the co-parent changes the Hasse diagram

Now retain ordinary divisibility on the induced vertex set `S`,
without requiring the quotient to lie in `S`. Its cover edges have a
generic and an exceptional part:

```text
x prec z
iff
  z/x is prime

or

  (x,z;z/x) is one of
  (2,8;4), (2,12;6), (3,12;4), (2,18;9).              (23)
```

**Proof.** Prime quotients have no intermediate divisor. Conversely,
let `q=z/x` be composite and let `p` be its least prime divisor. The
full divisor poset has the intermediate `xp`; an induced cover forces

```text
xp in {4,6}.
```

Thus `(x,p)` is `(2,2)`, `(2,3)`, or `(3,2)`. Moreover every proper
divisor `d` of `q` must have `xd in {4,6}`. In the first two cases this
leaves `q in {4,6,9}`; in the last it leaves `q=4`. Direct inspection
gives exactly the four displayed edges. QED.

Primes remain the uniform generic cover labels, but the finite startup
deletion inserts four nonprime Hasse labels. This is a feature of the
induced order, not of the internal multiplicative cospan.

## 7. Relative atoms and the twin-centre scar

For targets retained in `S`, the proper weak multiplicative fibre is
empty exactly for

```text
primes union {8,12}.                                      (24)
```

The strict fibre is empty exactly for

```text
primes union {prime squares} union {8,12}.                (25)
```

To prove this, let `p` be the least prime divisor of a composite
target and put `q=z/p`. If `p<q` and `q in S`, then `{p,q}` is a legal
strict pair. A missing `q` must be four or six, which gives only eight
or twelve after the least-prime condition is enforced. The remaining
strict case `p=q` is a prime square. The weak case retains that
diagonal.

By THM-2413, A014574 consists of the centres `c` for which `c-1,c+1`
are prime. The first terms are

```text
4,6,12,18,30,...
```

After deleting the startup holes, `12` is the unique retained centre
made into a composite weak atom:

```text
12=2*6=3*4
```

has both proper pairs killed. Every twin centre `c>=18` is divisible
by six and has the legal strict factorization

```text
c=2*(c/2),                       c/2>=9 and c/2 in S.     (26)
```

The startup defect therefore hits the twin-centre carrier exactly
once; it does not generate an infinite family of false prime-like
centres.

## 8. Chain difference versus divisor Möbius inversion

Return here to the full startup set `H={1,4,6}`. The linear marked-hole
term in (7) is an incidence-zeta transform. For addition put

```text
J_+(z)=#{h in H:h<z}.
```

Then

```text
Delta J_+(z)=J_+(z+1)-J_+(z)=1_H(z).                  (27)
```

This is Möbius inversion on the chain. It uses precisely the forward
difference whose lowering basis is the falling-factorial
Gregory--Newton basis in THM-2412.

For multiplication put

```text
J_x(z)=sum_(d|z)1_H(d)=mathbf_1 star_D 1_H(z),          (28a)
```

where `mathbf_1(n)=1` is the divisor-poset zeta function, not the
Dirichlet-convolution identity. Classical divisor Möbius inversion gives

```text
mu star_D J_x=1_H.                                    (28)
```

Define the finite hole Dirichlet series

```text
mathcal H(s)=sum_(h in H)h^(-s).
```

The two fixed-point transforms in (7) become

```text
H_+(X) -> H_+(X^2),

mathcal H(s) -> mathcal H(2s).                        (29)
```

Equations (27)--(29) are the exact additive/multiplicative incidence
analogy. They are **not** a Stirling identification: the Stirling
transform conjugates derivative and forward difference on polynomial
bases, but it does not turn Cauchy convolution into Dirichlet
convolution. The singleton hole `{2}` is the cheapest hostile—its
additive marked profile is a tail step, while its multiplicative
profile is the even-number divisibility indicator.

## 9. Exact companion

Run

```bash
python3 04-computation/operation_fibre_deletion_startup_scar_thm2433.py
python3 -O 04-computation/operation_fibre_deletion_startup_scar_thm2433.py
```

The dependency-free companion:

- checks (7) for all `56` three-hole subsets of `{1,...,8}` through
  thirty additive targets;
- verifies (9)--(16), including the `16,24,36` diagonal/overlap
  controls;
- builds both internal graphs through vertex `500`, computes
  reachability independently, and recovers (18)--(22);
- enumerates relative atoms and induced divisibility covers through
  `100,000`;
- sieves A014574 through centre `1,000,000` and finds `12` as the only
  retained artificial atom; and
- verifies both incidence inversions, including `10,000` exact divisor
  Möbius cases.

Every executable assertion uses an explicit optimization-safe
`require`. Normal and optimized transcripts must byte-match

```text
05-knowledge/results/operation_fibre_deletion_startup_scar_thm2433.out
```

after LF normalization.

## 10. Independent hostile audit

An independent derivation obtained (7) from the transposition action,
then separately reconstructed the two generating functions, the two
transitive closures, the four composite Hasse labels, and the relative
atom classification. Its finite scouts passed all `56` additive
three-hole sets, both multiplicative coefficient formulas through
`10,000`, both graph closures through `500`, and the atom/cover
classifications through `100,000`.

The audit also checked the main typing boundary: (22) concerns an
internal cospan whose quotient must be retained, whereas (23) concerns
the induced divisibility order and deliberately forgets that quotient.
No claim about prime density, twin-prime infinitude, or a causal
passage from additive holes to prime occurrence is made.
