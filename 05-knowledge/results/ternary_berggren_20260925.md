# Berggren recursion, its ternary ray decoder, and the Collatz clock boundary

**Status: PROVED elementary conjugacy and obstruction statements;
INHERITED exact parameter and ordinal charts; FINITE-EXACT controls.**
Collatz convergence remains **OPEN**. The new strengthening of the prior
edge transport is the arbitrary-distance incomparability theorem in section 5.
The ternary adding machine is an explicit synthesis of two inherited results,
not a new claim about p-adic dynamics or an identification of the full ordinary
Berggren tree with Collatz.

## 1. Inheritance and conventions

The closest mechanisms are
[THM-3756, odd-square ordinal Berggren affine descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md),
[THM-4057, CW18–CW25](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md),
the [Berggren edge-transport note, sections 6–7](collatz_mod6_20260921_berggren_edge_transport.md),
and the [inverse-fibre triadic odometer, B1–B3](arithmetic_braids_20260917_collatz.md).
The first two already give a complete intrinsic decoder of all primitive
Pythagorean triples. The last two already give the Collatz sibling ray and
its ternary completion, separately.

The canonical hostile is the legal edge `7->11`, whose Berggren parent
has roots `(7,3)` and is not an ordinary plus Collatz edge. The corrected
near miss is equating a sparse inverse fibre with its entire Berggren ray,
or forgetting that its k=1 element lies off that ray. The least-used
sidecar here is the clock measured in ordinary Berggren edges.
The [semicircle note](collatz_mod6_20260917_pythagorean_semicircle.md)
also separates arithmetic halving from Gaussian angle halving; no such
identification is made here.

| Live concept | Retained structure | Exact boundary |
|---|---|---|
| Odd square roots | primitive triple and its inverse | both roots are necessary |
| Ternary Berggren tree | unimodular parameter branches | neither branch count nor address is Collatz |
| Calkin–Wilf ordinal | complete guarded radix recursion | arithmetic labels change under the chart |
| Fixed-target ray | actual sibling return points | geometric gaps grow by four |
| Ternary refinement | every finite residue cylinder | completion is not ordinary-ray coverage |
| Consecutive orbit edges | two actual endpoint pairs | plus triangles never share an ancestry chain |

Anchor: an exact recursion with an inverse. Niche: locate the ternary
adding machine intrinsically on a Berggren ray. Wildcard: test whether
allowing arbitrarily long Berggren paths repairs the old edge obstruction.

Use positive coprime odd roots `s>t`, the ordered primitive triple

```text
Xi(s,t)=(st,(s^2-t^2)/2,(s^2+t^2)/2),
s=sqrt(C+B),             t=sqrt(C-B),                    (B1)
```

and the inherited Berggren branches

```text
B1(s,t)=(s+2t,t),
B2(s,t)=(2s+t,s),
B3(s,t)=(2s-t,s).                                       (B2)
```

Each branch preserves positivity, oddness and coprimality. Its determinant
is respectively 1,-1,1. The first coordinate strictly increases; the second
never decreases. The unique parent, except at `(3,1)`, is

```text
(s-2t,t)       if s>3t,
(t,s-2t)       if 2t<s<3t,
(t,2t-s)       if t<s<2t.                               (B3)
```

The missing equalities cause no unlisted nodes: `s=2t` contradicts oddness,
and `s=3t` with coprimality is the root. This is the proved canon mechanism,
not a census-derived claim. Classical triple matrices are recorded in
[Janičková–Csókási, section 2, page 2](https://arxiv.org/pdf/2304.05230);
the exact ordinal maps and the Collatz claims below come from the specified
repo results and the present proofs, not from that external paper.

## 2. The already exact ternary recursion on natural ordinals

Let Stern's sequence satisfy
`s(0)=0,s(1)=1,s(2k)=s(k),s(2k+1)=s(k)+s(k+1)`.
To avoid confusion with the root s, write `S(k)` for this sequence here.
The Calkin–Wilf pair at ordinal k is `(S(k),S(k+1))`. Then

```text
(m,n)=(S(k),S(k+1)) is a primitive Euclid pair m>n>0
with opposite parity       iff       k=3 or 5 mod6.     (B4)
```

THM-4057 proves this using `S(k)` even iff `3|k` and the order of adjacent
Stern values. The maps `m=(s+t)/2,n=(s-t)/2` and their inverse identify
these pairs with (B1). If `ell=floor(log2 k)` and
`k*=3*2^ell-1-k`, the three branches (B2) become

```text
k -> 2k-1,       k -> 4k*+3,       k -> 4k+3.            (B5)
```

Their inverses, on (B4), are

```text
k -> (k+1)/2             if k=1 mod4,
k -> ((k-3)/4)*          if k=3 mod8, k>3,
k -> (k-3)/4             if k=7 mod8.                   (B6)
```

They terminate at k=3. This is a complete conjugacy of labelled rooted
trees, retaining primitivity, parameter parity and a one-edge Berggren
clock. The middle branch needs reflection of the full binary suffix;
it is not an extra affine `3n+1` branch. It supplies a real ternary decoder,
with an explicit inverse, before any Collatz interpretation is imposed.
The present script rechecks the inherited chart through depth 8 and its
ordinal domain through 65535.

## 3. The induced return map on a fixed-target ray

Use the actual odd map

```text
U_sigma(x)=(3x+sigma)/2^v2(3x+sigma),    sigma=+1,-1,
```

on its positive integer domain. Fix an odd target y not divisible by 3.
Choose `k0>=2` with `2^k0 y=sigma mod3` and
`x0=(2^k0 y-sigma)/3>y`. The strict inequality excludes the degenerate
plus fixed edge at y=1,k0=2. All subsequent same-fibre sources are

```text
k_j=k0+2j,
x_j=(2^k_j y-sigma)/3,
x_(j+1)=4x_j+sigma,                    j>=0.            (B7)
```

They are coprime to y, and `(x_j,y)` are primitive root pairs. On the ray
through `(x0,y)`, height h means `B1^h(x0,y)=(x0+2hy,y)`. Put

```text
A=2^(k0-1),
H(j)=A(4^j-1)/3.
Then (x_j,y)=B1^H(j)(x0,y).                              (B8)
```

Equation (B8) is inherited from edge transport. Its exact first-return
interpretation is useful: at a legal source of exponent k, moving h
ordinary B1 edges changes `(3x+sigma)/y` from `2^k` to `2^k+6h`.
The next power of two with the same residue modulo 3 is `2^(k+2)`.
Consequently the next legal same-fibre point occurs **first** at

```text
h=2^(k-1).                                              (B9)
```

Thus the sibling operation is the first return to these marked points,
not one ordinary Berggren edge. One sibling step adds 2 to the exact
halving exponent, while its Berggren clock is `2^(k-1)`. These gaps
quadruple on successive returns. A possible k=1 fibre element uses B2
to enter the ray, as established in the inherited note; it is not included
in the present fixed-ray coordinate.

The ordinary inverse of H is explicit. A nonnegative height h belongs
to the marked fibre iff `A|h` and `1+3h/A` is a power of 4; if so,
the exponent of that power is j. This does not require tracing an orbit.
Most ordinary heights fail this guard. In fact the number through height
M is exactly

```text
1+floor(log_4(1+3M/A)),                                 (B10)
```

where the logarithm can be evaluated by integer power comparisons. The
set has ordinary density zero along its ray.

## 4. The ternary refinement tree and its constructive inverse

There is an exact commuting square

```text
       j  ----------------->  j+1
       | H                    | H
       v                      v
       h  ----------------->  4h+A.                    (B11)
```

Moreover, for distinct nonnegative i,j,

```text
v3(H(j)-H(i))=v3(j-i).                                  (B12)
```

For j>i the difference is `A*4^i(4^(j-i)-1)/3`; A and 4 are units
modulo 3. The elementary identity
`v3(4^d-1)=1+v3(d)` proves (B12). It follows by factoring the geometric
sum for `3` not dividing d, and then cubing `1+3^e u` to increase the
valuation exactly once at each factor of 3. This is the inherited
arithmetic-braids proof expressed in Berggren height.

For every r, H therefore induces a bijection on `Z/3^r Z`, compatible
with reduction to lower levels. It gives an isomorphism of rooted
ternary **residue-refinement trees**: the three children of a phase class
`j mod3^r` are its three lifts modulo `3^(r+1)`, and H maps them bijectively
to the three lifts of the corresponding height class. The child digits
may be permuted; literal digit labels are not asserted to be preserved.

The inverse is constructive. Starting from a solved phase `j_r mod3^r`,
test the three lifts `j_r+d*3^r`, d=0,1,2. Exactly one satisfies

```text
H(j_r+d*3^r)=h mod3^(r+1).                              (B13)
```

This gives each successive ternary digit, with no Collatz convergence
assumption. Passing to compatible sequences of residue classes extends H
uniquely to an isometric bijection `Z_3 -> Z_3`. Under that bijection
`h->4h+A` is addition by one. In particular it is a single cycle on every
finite quotient `Z/3^r Z`.

This is a genuine ternary-recursion isomorphism associated with the
Collatz sibling fibre and its Berggren ray. Its exact type is the
residue-refinement tree and its adding-machine action. It does **not**
identify that tree with the full ordinary PPT ancestry tree. The ordinary
membership guard in section 3 and the ternary inverse (B13) answer
different questions: every ternary height has a completed phase, while
only the power-of-four heights have a nonnegative integer phase. Losing
that distinction would turn residue completeness into a false coverage
claim. Even the complete original sibling fibre has only one fixed
Collatz target y, as the inherited odometer note already emphasizes.
For example, at k0=2 the ordinary heights are `0,2,10,42,170,...`;
their first three residues modulo 3 are `0,2,1`, exhausting that quotient
despite the rapidly increasing gaps on the ordinary ray.

## 5. A stronger obstruction: consecutive plus triangles are incomparable

**PROVED.** Let `x>1`, `y=U_+(x)>1`, and `z=U_+(y)`, all odd positive.
Form the two nondegenerate edge triangles with root pairs
`(max(x,y),min(x,y))` and `(max(y,z),min(y,z))`.
They are distinct, and neither is a Berggren ancestor of the other at
any distance. Only consecutive edges are covered by this theorem;
no claim about all pairs of triangles along an orbit is made.

This strengthens the old edge-transport section 7, which excluded only
immediate parent/child relations. The proof uses no convergence theorem.

**Two elementary path facts.** Along every nonempty Berggren path the
upper root strictly increases and the lower root never decreases.
A path that keeps the lower root fixed consists only of B1 steps.
Finally, a descendant of `(s,t)` whose lower root equals the old upper
root s must be

```text
(2l*s+t,s) or (2l*s-t,s),       l>=1.                   (B14)
```

Indeed its first step must be B2 or B3: any initial B1 would increase
the upper root above s before it could become the lower root. Thereafter
only B1 is allowed. This gives B2/B3 followed by B1^(l-1), proving (B14).

Write `3x+1=2^k y` and `3y+1=2^j z`. For an odd source greater than
one, growth is exactly exponent 1; a falling edge has exponent at least 2.
There are four exhaustive cases.

1. **Two rises, `x<y<z`.** Only `(y,x)` could be ancestor of `(z,y)`.
   Since `z<2y`, (B14) forces `z=2y-x`. The two exponent-one equations
   then give `y=2x+1=(3x+1)/2`, hence x=-1, impossible.

2. **A peak, `x<y>z`.** The two upper roots equal y, so a nonempty path
   is impossible. Equality would require z=x and
   `(2^(j+1)-9)x=5`. A positive coefficient on the left is at least 7,
   so this has no x>=3. Thus the nodes are also unequal.

3. **A valley, `x>y<z`.** Their lower roots equal y, so comparability
   would require `x-z` to be a multiple of `2y`, including zero for
   equality. Since j=1,
   `6(x-z)=(2^(k+1)-9)y-5`; therefore y divides 5. The only possible
   y>1 is 5, but `U_+(5)=1` is not a rising edge.

4. **Two falls, `x>y>z`.** The only possible direction is from `(y,z)`
   to `(x,y)`, so (B14) requires `x=2l*y+epsilon*z` with epsilon=+/-1.
   Reducing the two edge equations modulo y gives
   `y|(2^j+3epsilon)`. Here j>=2 and `2^j=(3y+1)/z<=3y+1`.
   If epsilon=+1, the positive odd quotient `(2^j+3)/y` is 1 or 3.
   Quotient 3 is impossible modulo 3. Quotient 1 gives `y=2^j+3`
   and `z=3+10/2^j`, nonintegral for j>=2.
   If epsilon=-1, the positive odd quotient `(2^j-3)/y` is less than
   3 and must be 1. Then `y=2^j-3,z=3-8/2^j`. At j=2 the middle
   value y is 1; at j=3 the endpoint z is even; at larger j it is
   nonintegral. All contradict the stated domain. This completes the proof.

The sign is decisive. The minus orbit segment `27->5->7` gives

```text
(27,5)=B1^2(7,5),                                       (B15)
```

so its consecutive triangles are comparable at distance two. In addition
the minus cycle `5<->7` gives the same unmarked triangle twice.
These are direct positive and hostile controls for the claimed scope.
An independent read-only proof audit checked the path lemma, all four
order cases, exact exponent placement, and the optimized replay: PASS.

Thus allowing a state-dependent number of forward or backward Berggren
steps does not repair the usual **consecutive plus edge** encoding: its
successive triangles are not on a common ancestry chain at all.
This does not obstruct the different oriented rational-triple lift in
[the companion note](ternary_triples_20260925.md), which defines another
edge relation on the same primitive geometric carrier.

## 6. No fixed projective conjugacy to a Berggren word

For an invertible real 2-by-2 matrix M, the quantity

```text
J(M)=tr(M)^2/det(M)                                     (B16)
```

is invariant under conjugation and under nonzero scalar multiplication.
Every Berggren parameter word is an integral matrix with determinant +/-1,
so its J is an integer. This also holds after any fixed projective change
between root ratios and Euclid ratios.

In contrast, an actual nonempty accelerated Collatz word with r odd
steps and total halving exponent K has affine matrix

```text
M=[[3^r, signed carry],[0,2^K]],
J(M)=(3^r+2^K)^2/(3^r*2^K).                             (B17)
```

Here r,K>=1 and the numerator is coprime to both 2 and 3, so J is a
noninteger rational. The sibling map `x->4x+sigma` has J=25/4, also
nonintegral. Therefore no fixed real projective coordinate change can
identify any such fixed affine word or sibling step with any fixed
Berggren word. The result is independent of carry sign; it does not
distinguish convergence behaviour.

The restriction to projective coordinates and fixed words is essential.
The Calkin–Wilf ordinal decoder is an arithmetic encoding, not a projective
change of slope. The induced sibling map uses a changing return time.
The companion rational lift also includes primitive content cancellation
by three, which is absent from a unimodular Berggren word. All three
survive this obstruction for explicitly different reasons.

## 7. Reproduction and scope

```text
python 04-computation/experiments/ternary_berggren_20260925.py
python -O 04-computation/experiments/ternary_berggren_20260925.py
```

Both match [the frozen output](ternary_berggren_20260925.out). The code uses
only exact standard-library arithmetic and explicit checks active under
optimization. Its bounded universes are:

- All 9841 Berggren words of depth 0..8: root and Euclid charts, ordinal
  conjugacy and both inverse maps. Independent 3-by-3 matrices verify
  every generated child triple.
- Every Calkin–Wilf ordinal 1..65535: exactly 21845 have the primitive
  Euclid parity/order predicate, agreeing with the two residue classes.
- Every odd x from 3 through 100001 with `U_+(x)>1`: 49992 consecutive
  pairs pass arbitrary-distance ancestry checks. All 4994 such starts
  through 10001 are checked independently using inverse 3-by-3 Lorentz
  matrices on triples. The signed controls (B15) and `5<->7` are retained.
- All 9841 parameter words through depth 8 have integral J; 1808 affine
  controls with r=1..16, K=r..64 and both carry signs have the predicted
  nonintegral denominator. The all-word obstruction follows from (B17),
  not the finite table.
- k0=2..9 and all phase residues through ternary depth 6: 8736 total
  residue images; 25920 distinct phase-pair valuation comparisons for
  i,j in 0..80. All 1944 height classes modulo 243 for those k0 values
  pass the constructive inverse and commuting-square checks. Ordinary
  heights 0..2000 test the separate power-of-four membership guard.
- Both signs, all valid fixed targets y<=101 and k=2..8: 237 first
  returns, with every intervening ordinary ray height tested directly.

The preserved predicates differ by object: (B4)–(B6) preserve full PPT
ancestry; (B7)–(B9) preserve a fixed Collatz target and exact halving
exponents; (B11)–(B13) preserve ternary congruence distance. Their required
sidecars are respectively the full ordinal, the marked return clock,
and the distinction between ordinary and completed phase. None alone
gives a root certificate for every positive integer.
