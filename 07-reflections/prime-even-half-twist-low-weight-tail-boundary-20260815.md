# Prime-even half-twist low-weight tail boundary

**Status: exact elementary formulas for the three displayed tangent families,
plus FINITE-EXACT controls.  The uniform height-`<=7` zero/low-multiplicity
atlas is OPEN.  This is an unnumbered research companion, not an all-prime
classification and not an LRC(14) closure.**

Codex, 2026-08-15.

## Outcome

For `Q=2p`, with `p` an odd prime, put

```text
B_r={ell in Z/(2p)Z: ||r(2ell+1)/(4p)||<1/14}
```

and write `w(r,s)=|B_r intersect B_s|`.  The three height-eight tangent
families visible at the end of the bounded prime-even census admit exact
strict-endpoint formulas.  If

```text
N_AA(p)=#{odd x: 4p/21 < x < 2p/7, 3x == 2p (mod 5)},
N_35(p)=#{x:     2p/21 < x < p/7,  3x == p  (mod 10)},
N_53(p)=#{odd x: 4p/35 < x < p/7,  5x == p  (mod 6)},
```

then, for the representative tangent lifts described below,

```text
w_AA=2 N_AA(p),    w_BE(3/5)=4 N_35(p),
w_BE(5/3)=4 N_53(p).
```

The exact last boundaries at which these expressions can be at most eight
are

```text
A/A:       p=521 among odd p,
B/E, 3/5: p=601 among odd p,
B/E, 5/3: p=609 among odd p, but p=593 among primes.
```

Consequently every member of these three tangent families has weight at
least `12` at every prime `p>=607`.

This does **not** prove that these are all low-weight pairs.  That missing
exhaustion is precisely the remaining structural gap.

The independent direct-mask controls add three bounded facts:

1. the first prime beyond the previous `p<=599` atlas, `p=601`, has no
   seven-owner cover;
2. for every one of the `58` primes `607<=p<=997`, every normalized positive
   pair weight is at least `12`, every zero neighbor has raw rational height
   at most seven, and the zero degrees at roots `A=1,B=2,E=4` are
   `(19,31,23)`; and
3. one first-prime representative in each of the `96` invertible residue
   classes modulo `420` gives the same height-seven candidate clique number
   `6`.  Every maximum candidate clique is three complementary pairs with
   type profile `A^4BE`.

Items 1--3 are finite exact statements only.  In particular, the residue
controls are not an extrapolation to every larger prime in their classes.

## Inheritance and corrected near miss

The closest proved ordinary-interval mechanism is
[THM-3423](../01-canon/theorems/THM-3423-odd-interval-ratio-complement-and-dyadic-clique-law.md).
It turns disjoint odd intervals into a bounded-rational ratio graph and
classifies its complementary-pair equality case.  The dyadic source is
[THM-3435](../01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md),
which is **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED** for the exact
fibre decomposition and the target-free necessary gates used here: exactly
seven owners, at least two odd (`A`) owners, and at least one owner with
`v_2(r)>=2` (`E`).

The bounded predecessor is
[the prime-even cap-seven finite atlas](prime-even-half-twist-cap-seven-finite-atlas-20260815.md).
After the proved gates and the fixed-fibre invoice, each unresolved
seven-owner profile has cover excess

```text
Omega=sum_r |B_r|-2p <= 8.
```

For a literal cover, any one pair intersection is at most `Omega`: each point
in `B_r intersect B_s` already contributes at least one to
`sum_x(m_x-1)=Omega`.  Thus every pair in a remaining cover must have weight
at most eight.

MISTAKE-393 is the relevant audit repair.  A degree-`d` circle map is
bijective on each inverse-image component, not on their union, and the
two-sheet partner statement is valid only after canonical sign
normalization.  Here the associated corrected near miss is equally concrete:
a raw ratio modulo `p` does not determine an `A/A` intersection.  Its lift
modulo `4p`, equivalently an orientation bit, changes the intersection from a
linear overlap to a boundary tangency.  The derivation below retains this bit.

## Base intervals, sheets, and the orientation bit

Represent a literal sheet by the odd phase `t=2ell+1 (mod 4p)` and use
centered representatives after projection.  Odd-unit normalization gives
three owner types:

- `A`, an odd coefficient: one oriented sheet over the widened interval
  `|x|<2p/7`;
- `B`, a coefficient `2u` with `u` odd: two sheets over the odd interval
  `|x|<p/7`; and
- `E`, a coefficient `4v`: two sheets over `|y|<p/14`, including the
  reflection-fixed base fibre.

It follows that an `A/even` intersection is counted once over its base
intersection, an even/even intersection is counted twice, and only `A/A`
needs the lift/orientation bit.  Complementary `A` coefficients `r` and
`2p-r` select opposite sheets and are disjoint.

These projection statements organize the calculation; the companion's
controls do not trust them.  They reconstruct each literal mask directly
from

```text
14 min(r(2ell+1) mod 4p, 4p-r(2ell+1) mod 4p) < 4p.
```

## Exact `A/A` tangent formula

Normalize the first owner to `1`.  Consider the raw ratio `3/5` and choose
the odd lift `s` satisfying

```text
5s-3 == 2p (mod 4p).
```

The other signs and the reciprocal ratio `5/3` are obtained by coefficient
sign and owner interchange.  If `x,y` are the centered odd phases in the two
widened `A` intervals, then multiplication by the odd sheet coordinate gives

```text
3x-5y == 2p (mod 4p).
```

Because `|x|,|y|<2p/7`, the only possible representatives are

```text
3x-5y=+2p or -2p.
```

For the positive-`x` representative of the first equation,
`y=(3x-2p)/5`.  The strict condition `-2p/7<y<2p/7` is equivalent to

```text
4p/21 < x < 2p/7,
```

and integrality is `3x==2p (mod 5)`.  The negative solution is its central
reflection.  There is no deck factor for an `A/A` pair, so

```text
w_AA=2*#{odd x: 4p/21<x<2p/7, 3x==2p (mod 5)}.
```

The discarded lift `5s-3==0 (mod 4p)` would instead give a central, linearly
large intersection.  This is exactly the information destroyed by reducing
the coefficient ratio modulo `p`.

In integer endpoint form the one-sign count is

```text
L=floor(4p/21)+1, U=floor((2p-1)/7),
N_AA=#{x in [L,U]: x odd, 3x==2p (mod 5)}.
```

The parity and congruence select one residue modulo `10`.

## Exact `B/E` tangent formulas

Normalize the `B` owner to coefficient `2` and write an `E` owner as `4v`.
For raw ratio `3/5`, the base relation is `2v==3/5 (mod p)`.  If `x` is in
the centered odd `B` interval and `y` is in the centered `E` interval, then

```text
3x-10y == 0 (mod p).
```

The left side has magnitude below `8p/7`.  It is odd, whereas zero is even,
so the zero representative is impossible.  The only representatives are
`+p` and `-p`.  For positive `x` and `3x-10y=p`, the strict bounds
`0<x<p/7` and `|y|<p/14` reduce exactly to

```text
2p/21 < x < p/7,   3x==p (mod 10).
```

Central reflection supplies the opposite sign, and the even owner pullback
has two sheets.  Hence

```text
w_BE(3/5)=4*#{x: 2p/21<x<p/7, 3x==p (mod 10)}.
```

The congruence automatically makes `x` odd.  Its strict integer endpoints
are `floor(2p/21)+1` and `floor((p-1)/7)`.

For raw ratio `5/3`, the identical argument begins with

```text
5x-6y == 0 (mod p).
```

Parity again excludes zero and the size bound leaves only `+/-p`.  Substituting
`y=(5x-p)/6` gives

```text
4p/35 < x < p/7,   5x==p (mod 6),
w_BE(5/3)=4*#{odd x: 4p/35<x<p/7, 5x==p (mod 6)}.
```

The strict integer endpoints are `floor(4p/35)+1` and
`floor((p-1)/7)`.  Reversing the ordered edge reciprocates the raw ratio and
does not change its literal weight.

## Count recurrences and exact thresholds

All three one-sign counts are arithmetic-progression counts.  Shifting
`p` by `210` preserves their selected residue and changes the endpoints by

```text
A/A:       (L,U) -> (L+40,U+60),
B/E, 3/5: (L,U) -> (L+20,U+30),
B/E, 5/3: (L,U) -> (L+24,U+30).
```

The first shift exposes two new points in the selected class modulo `10`;
the other shifts expose one new point modulo `10` and modulo `6`,
respectively.  Therefore

```text
N_AA(p+210)=N_AA(p)+2,
N_35(p+210)=N_35(p)+1,
N_53(p+210)=N_53(p)+1.
```

An exhaustive check of the residue representatives gives the exact last
weight-at-most-eight boundaries `521`, `601`, and `609` stated above.  Since
`609` is composite and the last prime in that third family is `593`, every
one of the displayed tangent weights is at least `12` for prime `p>=607`.

Direct literal-mask values `(p,w_AA,w_35,w_53)` include

| `p` | `w_AA` | `w_35` | `w_53` |
|---:|---:|---:|---:|
| `401` | `8` | `8` | `8` |
| `521` | `8` | `8` | `12` |
| `547` | `10` | `8` | `12` |
| `571` | `12` | `12` | `8` |
| `593` | `12` | `12` | `8` |
| `601` | `12` | `8` | `12` |
| `607` | `12` | `12` | `12` |
| `1009` | `18` | `20` | `20` |

For each row the direct-mask path finds both sign/reciprocal presentations in
each applicable ordered type sector and agrees exactly with the interval
formula.

## Boundary prime and full-root scan

The independently frozen profile DFS decides `p=601` negatively with

```text
15 profiles, 65 recursive states, 50 candidate branches,
```

using no node or time cap.  This is the first exact negative beyond the
predecessor atlas's `p<=599` boundary.

For every prime from `607` through `997`, the new companion then compares
each normalized root `1,2,4` against every nonempty coefficient
`1<=r<2p`, using only direct literal masks.  This is `58` complete root rows,
not a sample of coefficients.  It finds

```text
zero degrees (A,B,E)=(19,31,23),
number of edges with 0<w<=8: 0,
minimum positive weight: 12,
zero neighbors outside raw height <=7: 0.
```

The histogram of the three root-wise positive minima is frozen in the output
companion.  This scan neither assumes the tangent formulas nor restricts
candidates to a rational atlas.

## Modulo-420 candidate graph

Let

```text
D7={+/-a/b mod p: gcd(a,b)=1, a,b>0, a+b<=7}.
```

For each root type, retain both coefficient lifts in `1<=r<2p` of every
ratio in `D7`, together with the root, and put an edge exactly when the two
direct literal masks are disjoint.  This retains coefficient type and the
`A` lift; it is not the bare Cayley graph on ratios.

The companion takes the first prime at least `607` in every invertible class
modulo `420`.  The `96` representatives range from `607` to `3319`.  On every
row it finds

```text
root degrees (19,31,23), maximum clique size 6,
four normalized maximum-clique presentations,
type profile A^4BE,
three complementary pairs in every maximum clique.
```

The four presentations are the two `A` root normalizations, one `B` root
normalization, and one `E` root normalization of the same equality orbit.
This is a finite residue-control atlas.  Its role is to show that the proposed
height-seven reduction has the expected six-clique obstruction in every
residue class; it does not prove that all actual zero or low-weight pairs have
entered the candidate set.

## Source, target, loss, and sidecars

| field | exact content |
|---|---|
| source | labelled literal masks `B_r` on all `2p` sheets |
| target | a typed coefficient-ratio graph with exact pair weight `w(r,s)` |
| map | normalize one owner by an odd unit, project coefficients modulo `p`, then record the reduced ratio |
| preserved | coefficient type under odd-unit normalization; literal pair weight when computed before quotienting |
| lost by raw ratio | the `A` lift/orientation, literal sheet locations, the fixed `E` fibre, complementary-pair labels, higher multiplicities, and the full union |
| required sidecars | ordered endpoint types `A/B/E`; the residue `(b s-a r) mod 4p` for an `A` endpoint; complementary partner `r <-> 2p-r`; fixed-fibre incidence |
| cover predicate | seven distinct owners cover every literal sheet |
| necessary graph predicate | every pair has `w<=Omega<=8`; after the desired atlas and the tangent thresholds, all seven vertices would have to form a `D7` zero clique |
| destroyed implication | a graph clique alone does not certify a literal cover, because it forgets union and higher overlap |

This is a weighted graph problem with ties, not an intrinsic tournament.

## Exact remaining gap

The missing statement is the following uniform typed near-complement lemma:

> For every sufficiently large odd prime `p` (with the desired threshold
> `p>=607`), every normalized typed coefficient pair with literal weight
> `w<=8` either has weight zero and raw ratio in `D7`, or belongs, with its
> endpoint types and `A` orientation lift attached, to the displayed
> `+/-3/5,+/-5/3` tangent shell.

The exact formulas prove that the tangent alternative has weight at least
`12` for every prime `p>=607`.  The desired lemma would therefore force every
pair in a putative cover to be a `D7` zero edge.  The modulo-420 candidate
atlas then exhibits the expected obstruction: its clique number is only six.

What remains OPEN is the uniform exhaustion, not another search for larger
controls.  In particular, the present work does not claim that `p=7,19` are
the only all-prime solutions.

## Reproduction and hashes

Run from the repository root:

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_prime_even_half_twist_low_weight_tail_probe_20260815.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_prime_even_half_twist_low_weight_tail_probe_20260815.py
```

Normal and optimized outputs are byte-identical.  The script audits that it
contains no Python `assert`, pins the predecessor companion by hash, and pins
its semantic payload.  LF-normalized SHA-256 values are

```text
script:   4bb0c853fa3bebabf1414d75aa2b43d4b00c22a2b49cbaaac5e06bde5986c4ba
output:   a9012981cef620a0eeac6883f6e6cb080272c3d99cb73703a0b45212f9cb9b80
semantic: 17411c3f578f369ca6ab424db0db6730b547db5deb598e87fa341da49dd6473b
base:     ed54b01f9bf155643b6407c6af8ee6f15c7a099f8ac3b60399dd49b496fb1d12
```

The frozen transcript is
[lrc_prime_even_half_twist_low_weight_tail_probe_20260815.out](../05-knowledge/results/lrc_prime_even_half_twist_low_weight_tail_probe_20260815.out),
and the executable companion is
[lrc_prime_even_half_twist_low_weight_tail_probe_20260815.py](../04-computation/lrc_prime_even_half_twist_low_weight_tail_probe_20260815.py).
