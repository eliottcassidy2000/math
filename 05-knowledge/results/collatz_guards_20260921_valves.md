# The 23 cylinder, first-reset boundary, and an exact Beatty density

**Status: PROVED elementary cylinders, spatial density, and infinite family;
CITED one-cycle exclusion; FINITE-EXACT replay.** Collatz convergence and
complete positive/signed cycle classification remain OPEN. No priority claim.

## Inheritance and the live board

The [affine audit](collatz_blueprint_20260921_affine.md) retains guard and
carry; the [completion theorem](arithmetic_braids2_20260917_inverse_completion.md)
separates arbitrary finite prefixes from one fixed integer's orbit. The
closest hostile is an arbitrarily long all-one halving word. The corrected
near miss is a ternary residue advertised as forcing a binary valuation.
The least-used coordinate is the odd cofactor of n+1.

| Object | Operation / exact coordinate | Decisive boundary |
|---|---|---|
| 23 trajectory | exact word (1,1,5) | one dyadic cylinder, not the whole row 5 mod9 |
| First reset | n+1=2^(r+1)u | actual k depends on u, not on r alone |
| Carry | gap D and odd integrality | growth, strict descent, or return |
| Spatial frequency | disjoint cylinders / CRT | does not prescribe orbit time frequency |
| Mechanical word | Beatty digit positions | slope log_2(3), not a golden-ratio coding |

Anchor: replace the proposed valve by a precise descent certificate. Niche:
derive the complete first-reset counting law. Wildcard: identify what the
23 family does after its first large halving. Results below are rederived;
only the explicitly marked one-cycle corollary imports deeper literature.

Throughout, T is the accelerated map on positive odd integers, and k is
the exact exponent v2(3n+1). A reset means an exponent at least two after
a maximal initial run of exponents one. It does not mean descent below
the block's starting value.

## 1. The exact 23 valve occupies one cylinder

The first three exact exponents are (1,1,5) iff

```text
n=23 mod256.
```

Indeed the ordered identity is `128*T^3(n)=27n+19`. Requiring its last
quotient odd gives `27n+19=128 mod256`, equivalent to n=23 mod256;
backward substitution gives all three exact valuations. Explicitly,

```text
23+256t -> 35+384t -> 53+576t -> 5+54t,       t>=0.
```

The endpoint is below n because the carry threshold is `19/101<1`.
Every member of this cylinder ends at 5 mod9, even when its source is
in a different row. Intersecting the source with n=5 mod9 gives

```text
n=23+2304q,       T^3(n)=5+486q,             q>=0.
```

Thus its relative density is 1/128 among odd inputs in row 5 mod9,
or 1/256 among all integers in that row. These are exact progression
densities. Other row members need not have this word: the smallest one
with three successive growth steps is `95 -> 143 -> 215 -> 323`.
The smaller 41 has third endpoint 71 but descends at its first step;
it is not an all-growth example.

## 2. The exact mod-nine support graph has no privileged valve

For input residue a and exponent k,

```text
T(n)=(3a+1)*2^(-k) mod9.
```

Since 2 generates the six units modulo nine, every input residue has
all six units as possible outputs. Every proposed k>=1 occurs in each
row, by CRT between its exact dyadic cylinder and n=a mod9. No row
divisible by three is a possible output.

| Input a mod3 | Outputs for k=1,2,3,4,5,6 mod6 |
|---|---|
| 0 | 5,7,8,4,2,1 |
| 1 | 2,1,5,7,8,4 |
| 2, including a=5 mod9 | 8,4,2,1,5,7 |

Spatial weights of these six branches are `(32,16,8,4,2,1)/63`, from
`sum_{j>=0}2^(-k-6j)`. This is a counting measure on source integers,
not a claim that any fixed orbit samples a Markov chain. The support
graph forgets k, the binary guard, and ordinary height.

## 3. Every first reset has an exact size boundary

Write uniquely

```text
n=2^(r+1)u-1,       r=v2(n+1)-1>=0,       u odd positive.
```

The first r exponents are one, with intermediate values
`n_j=2^(r+1-j)3^j u-1`. The next exponent and endpoint are

```text
k=1+v2(3^(r+1)u-1)>=2,
y=T^(r+1)(n)=(3^(r+1)u-1)/2^(k-1),
2^(r+k)y=3^(r+1)n+B_r,       B_r=3^(r+1)-2^(r+1).       (1)
```

For a specified n,r the minimum integer exponent satisfying strict
descent in (1) is

```text
k_min=1+floor(log_2((3^(r+1)n+B_r)/(2^r n))).             (2)
```

This is an inequality threshold, not an exponent one may choose: the
actual k is already forced by u. Integer comparisons compute (2) exactly.

An integrality improvement is useful. Put `D=2^(r+k)-3^(r+1)`, nonzero.
Equation (1) becomes

```text
y-n=1-(D*u+1)/2^(k-1).                                 (3)
```

If D<0, the positive carry in (1) proves y>n. If D>0, (3) is less
than one; because y,n are odd integers, y<=n. Equality holds precisely
when `D*u=2^(k-1)-1`. Otherwise y<=n-2. These statements are elementary
and do not assume a cycle classification.

**CITED strengthening.** Equality would be a one-rise/one-fall cycle
of shortcut Collatz. Steiner's exclusion of nontrivial one-cycles,
as recorded and incorporated in [Simons--de Weger, sections 1.4 and 2.1,
Theorem 3(b)](https://math.deweger.net/papers/%5B35a%5DSidW-3n%2B1-v1.44%5B2010%5D.pdf),
therefore gives y<n iff D>0 for every n>1. Their rise/fall counts are
r+1 and k-1. This imports a restricted cycle theorem, not a proof of
all-cycle uniqueness. The following density theorem does not need it.

## 4. General cylinders and an elementary descent density

For fixed r>=0,k>=2, the exact block `(1 repeated r times, k)` has

```text
u=3^(-(r+1))*(1+2^(k-1)) mod2^k,
n=2^(r+1)u-1 mod2^(r+k+1).                              (4)
```

This is one odd source class, hence relative density `2^(-r-k)` among
odd integers. CRT gives the same density inside every fixed odd-modulus
row, including odd n=5 mod9. Consequently r has law `2^(-r-1)`, and
conditional on r the reset exponent has law `2^(1-k)`, k>=2. These two
coordinates are independent for this spatial counting measure.

Let `alpha=log_2(3)` and

```text
kappa(r)=ceil((r+1)alpha)-r.
```

The region D>0 is exactly k>=kappa(r). Its density among odd starts is

```text
P=sum_(r>=0) 2^(1-r-kappa(r))
 =sum_(L>=1) 2^(-floor(L*alpha))
 =0.7137254976758912... .                               (5)
```

This is also the density of **strict descent at the first reset**, with
no external cycle theorem needed. For each fixed (r,k), the equality
equation in (3) has at most one integer source. Truncate to finitely
many r,k; their equality sources are finite. The discarded r and k
tails have arbitrarily small density by (4), so all equality sources
together have density zero. The same proof works in every fixed
odd-modulus row. The tail estimates also establish the existence of
the density of these infinite unions, rather than merely adding
infinitely many densities without justification.

For an exact numerical certificate set
`e_L=(3^L).bit_length()-1=floor(L*alpha)` and `S_N=sum_{L=1}^N 2^(-e_L)`.
Then

```text
S_N < P < S_N+3^(-N) <= S_N+2^(-N).                     (6)
```

The bound follows from `2^(-floor(L*alpha))<2*3^(-L)`.
The replay records rational endpoints for N=64; no floating-point
logarithm chooses any exponent.

**An exact mechanical-word connection.** Since 1<alpha<2 is irrational,
the positions floor(L*alpha) are distinct. Thus (5) is the binary real
whose 1-digits occupy exactly those Beatty positions. For j>=1 its
jth digit is `floor((j+1)/alpha)-floor(j/alpha)`, a mechanical coding
of slope 1/alpha. The threshold increments kappa(r+1)-kappa(r) form
the corresponding 0/1 staircase with slope alpha-1. These codings
describe first-reset thresholds and source density; neither describes
the halving word along a specified Collatz orbit. No golden ratio enters.

**Elementary corollary: P is irrational.** Its binary digit frequency
is `1/alpha`, which is irrational. A rational number has an eventually
periodic binary expansion, whose digit frequency is rational. The
nonperiodic mechanical coding therefore rules out rational P; no claim
about transcendence is needed here.

## 5. The 23 micro-example has an infinite opposite regime

For a>=1 put `n_a=3*2^(2a+1)-1`. Every n_a is 5 mod9 and has exactly
2a initial exponents one. Its first reset is determined explicitly:

```text
k_a=4+v2(a+1),
T^(2a+1)(n_a)=(3^(2a+2)-1)/2^(3+v2(a+1)).              (7)
```

For proof, (1) has u=3; the elementary two-adic identity
`v2(3^(2a+2)-1)=3+v2(a+1)` follows by factoring 9^(a+1)-1.
If a+1 is odd the geometric sum is odd; each doubling of the exponent
adds one valuation because 9^j+1=2 mod8.

| a | Start | Growth steps | Reset exponent | Reset endpoint |
|---|---:|---:|---:|---:|
| 1 | 23 | 2 | 5 | 5 |
| 2 | 95 | 4 | 4 | 91 |
| 3 | 383 | 6 | 6 | 205 |
| 4 | 1535 | 8 | 4 | 7381 |

The first three reset endpoints descend. **For every a>=4 the first
reset endpoint exceeds n_a.** For a=4 compare `2^12<3^9`. For a>=5,

```text
2^(2a+k_a) <= 16(a+1)4^a < 3*9^a=3^(2a+1).
```

The strict inequality holds at a=5 and its right/left ratio increases
by `(9/4)*(a+1)/(a+2)>1`. Thus D<0, and (3) proves growth. Every
preceding step grows as well. The first descent below n_a, if it occurs,
therefore needs a later block for all a>=4. This proves arbitrarily
long delays inside the advertised valve row without producing a
divergent integer orbit.

## 6. Connection contract, replay, and next proof obligation

The map from a starting integer to (r,u,k) preserves exact first-reset
size through (1)--(3). Projection to (r,k) preserves the slope boundary
but loses the possible return equality; projection to a mod-nine state
loses even the slope. Restoring u and the ordered carry makes the
boundary exact. Counting (4) produces the mechanical density (5),
which preserves a source-frequency predicate and discards orbit chronology.

Run `python 04-computation/experiments/collatz_guards_20260921_valves.py`.
The stdlib [script](../../04-computation/experiments/collatz_guards_20260921_valves.py)
and [JSON](../../04-computation/experiments/collatz_guards_20260921_valves.json)
declare exhaustive source/cylinder universes and finite family controls.
All arithmetic is exact; explicit checks remain active under `python -O`.
The finite list of returns is not used as an all-cycle classification.

The productive remaining question concerns **sequences of reset blocks**:
prove that the actual blocks of every fixed n>1 eventually repay their
accumulated growth and carry. The spatial success rate in (5) does not
force that sequence to sample its success region. Neither a finite
mod-nine automaton nor the first-reset theorem resolves this chronology.
