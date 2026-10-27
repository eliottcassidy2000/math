# Exact aliquot fibers and the odd-part carry under doubling

## Scope and notation

For a **positive integer** `n`, let

\[
\sigma(n)=\sum_{d\mid n}d,\qquad
s(n)=\sigma(n)-n=\sum_{\substack{d\mid n\\d<n}}d.
\]

Thus `s(1)=0`. A displayed chain ending at zero stops there; this note does
not define `s(0)`. Inverse fibers mean `s^{-1}({t})`, not an inverse function.
Every proof below is elementary and self-contained. Finite enumerations have
proved cutoffs, not empirical search limits.

## 1. Doubling retains an odd-part carry

Write `n=2^a m`, where `a>=0` and `m` is odd. Multiplicativity of the divisor
sum for coprime factors gives

\[
\sigma(n)=(2^{a+1}-1)\sigma(m),\qquad
s(n)=(2^{a+1}-1)\sigma(m)-2^a m.
\]

Substitution at `a+1` proves the exact law

\[
\boxed{s(2n)=2s(n)+\sigma(m).}
\]

The additive term is the divisor sum of the **odd part**, not a universal
constant and not a term that can be discarded. Equivalently,

\[
s(2^a m)=2^a\bigl(2\sigma(m)-m\bigr)-\sigma(m).
\]

For repeated doubling, the two-coordinate state

\[
(x,c)=\bigl(s(n),\sigma(\operatorname{oddpart}(n))\bigr)
\]

has the closed update `(x,c) -> (2x+c,c)`. This affine law describes
**doubling**, not iteration of `s`. Section 2 shows that the exact pair even
identifies `n`, but no simple affine update for aliquot iteration or globally
minimal state description is claimed.

| Fixed odd part | Carry `c` | Inputs | Aliquot outputs |
|---|---:|---|---|
| `5` | `6` | `5,10,20,40` | `1,8,22,50` |
| `25` | `31` | `25,50,100` | `6,43,117` |

For example, `43=2*6+31` and `117=2*43+31`. The change from odd part `5`
to odd part `25` changes the carry from `6` to `31`.

### The complete 13-fiber detects lost state

The sources `27` and `35` have the same aliquot sum:

\[
s(27)=1+3+9=13=s(35)=1+5+7.
\]

But their odd-part divisor sums are `40` and `48`, so

\[
s(54)=2\cdot13+40=66,\qquad
s(70)=2\cdot13+48=74.
\]

Consequently **no function of `s(n)` alone can produce `s(2n)` for all
positive integers `n`**. The odd-part carry repairs this failure exactly.
This is a genuine compositional obstruction, not a conclusion inferred from
the mere recurrence of the number 13 in unrelated problems.

## 2. Two exact observations recover the source itself

**Theorem (exact two-observation reconstruction).** The map

\[
n\longmapsto\bigl(s(n),s(2n)\bigr)
\]

is injective on positive integers. More precisely, for an observed pair
`x=s(n)`, `y=s(2n)`, set

\[
c=y-2x,\qquad d=y-x,\qquad a=v_2(d).
\]

Then

\[
\boxed{n=2^{a+1}c-d
      =2^{v_2(y-x)+1}(y-2x)-(y-x).}
\]

**Proof.** Write `n=2^a m` with `m` odd. The doubling law gives
`c=sigma(m)`. Adding this to the tower formula gives

\[
d=x+c=2^a\bigl(2\sigma(m)-m\bigr).
\]

The factor `2sigma(m)-m=m+2s(m)` is **positive and odd**. Consequently the
2-adic valuation of `d` is exactly the exponent `a` of the original source,
not merely a lower bound. Finally,

\[
2^{a+1}c-d=2^a m=n.
\]

This reconstructs a unique source from every actual pair. ∎

There is also an exact **image guard** for an arbitrary integer pair `(x,y)`:

1. Require `c=y-2x>0` and `d=y-x>0`.
2. Set `a=v2(d)` and `m=2c-d/2^a`; require `m>0`.
3. Require `sigma(m)=c`.

The recovered `m` is automatically odd. These conditions are necessary by
the proof and sufficient: for `n=2^a m`, the divisor-sum formulas give
`s(n)=d-c=x` and `s(2n)=2d-c=y`. Positivity alone is insufficient. For
example, `(x,y)=(1,4)` gives `c=2,d=3,m=1`, but `sigma(1)=1`, not `2`;
the candidate pair must be rejected.

Examples are `(0,1)->1`, `(6,43)->25`, `(43,117)->50`, and
`(117,265)->100`. In the last example the observation `265=s(200)` is a
**doubling observation**, not the next term `65=s(117)` of the forward
aliquot chain.

The doubling tower with odd part `25` makes the recovered exponent visible:

| Source `n` | Observation `x=s(n)` | Observation `y=s(2n)` | Carry `c` | Difference `d` | Recovered `v2(d)` |
|---:|---:|---:|---:|---:|---:|
| `25` | `6` | `43` | `31` | `37` | `0` |
| `50` | `43` | `117` | `31` | `74` | `1` |
| `100` | `117` | `265` | `31` | `148` | `2` |

Here `37=2sigma(25)-25` is the positive odd core of the difference, and
each subsequent doubling increments its 2-adic valuation exactly once.

### Exact integers cannot be replaced by residues

The complete 13-fiber already supplies the strongest directly relevant
counterexample. For every `a>=0`,

\[
\boxed{s(2^a\cdot35)-s(2^a\cdot27)=8(2^a-1).}
\]

Indeed both initial values are `13`, their carries are `48` and `40`, and
iteration of `(x,c)->(2x+c,c)` gives
`s(2^a n)=2^a s(n)+(2^a-1)c`. Their entire infinite doubling traces are
therefore identical **modulo 8**, despite different exact next values `66`
and `74`. More generally, within any fixed exact `s`-fiber, the doubled-trace
difference is `(2^a-1)` times the difference of the two odd-part carries.

For a separate small-prime check, for every `a>=0`, the prime cases `3` and
`7` give

\[
s(2^a\cdot7)-s(2^a\cdot3)=4(2^a-1).
\]

Thus their entire infinite doubling traces agree modulo `4`, although the
sources are different. Exact integer observations recover the source; even
all observations modulo this fixed modulus do not. This is not a claim
about unbounded-precision observations or about an aliquot-iteration clock.

### Corollary: any known prime dilation recovers the source

The role of `2` is not essential, but **primality is**. Fix a known prime
`p`, and write `n=p^a m` with `p` not dividing `m`. Then

\[
\sigma(n)=\frac{p^{a+1}-1}{p-1}\sigma(m),\qquad
\boxed{s(pn)=p\,s(n)+\sigma(m).}
\]

For exact observations `x=s(n)` and `y=s(pn)`, define

\[
c=y-px,\qquad d=y-x.
\]

The carry is `c=sigma(m)`, and

\[
d=p^a\bigl(pc-(p-1)m\bigr)
 =p^a\bigl(m+p\,s(m)\bigr).
\]

The last factor is positive and is congruent to `m`, a nonzero residue
modulo `p`. Hence `a=v_p(d)`, and the source is

\[
\boxed{n=\frac{p^{v_p(y-x)+1}(y-px)-(y-x)}{p-1}.}
\]

Thus `(s(n),s(pn))` is injective on positive integers for every fixed
**known prime** `p`. For a promised valid pair, reconstruction itself needs
only integer arithmetic and a `p`-adic valuation; no search over preimages
is required.

For an **arbitrary** integer pair, validity additionally requires the
divisor-sum guard. Require `c>0,d>0`, set `a=v_p(d)`, and form

\[
m=\frac{pc-d/p^a}{p-1}.
\]

Require `m>0`, `p` not dividing `m`, and `sigma(m)=c`. These conditions are
necessary and sufficient: the formulas re-encode the pair as
`x=(d-c)/(p-1)` and `y=(pd-c)/(p-1)`. The displayed division defining `m`
is automatically exact for integer `x,y`: `d-c=(p-1)x`, and division by
`p^a`, congruent to `1` modulo `p-1`, preserves the congruence `d=c`
modulo `p-1`. Its residue modulo `p` is nonzero, so the explicit `p`-free
condition is also automatic after integrality, but retaining it in a
validator makes the domain guard visible.

**Composite-dilation hostile check.** Taking multiplier `4` and source `3`
gives `x=s(3)=1`, `y=s(12)=16`. Substituting `p=4` into the prime-only
formula, with the largest power of `4` dividing `15` equal to `1`, would
give `(4*12-15)/3=11`, not `3`. This disproves extension of this formula
to arbitrary composite multipliers. The false candidate even has `sigma(11)=12=c`, so the divisor-sum guard cannot compensate for dropping primality. It does not assert that every different
reconstruction method fails for composite multipliers.

This corollary uses two labeled observations, one after a specified prime
dilation. It does not solve an inverse aliquot problem from one observation,
and it says nothing about termination of forward aliquot iteration.

## 3. A complete finite bound for every target greater than one

**Theorem.** Fix an integer `t>=2`. If `s(n)=t`, then either

1. `n=(t-1)^2` and `t-1` is prime; or
2. `n<=floor((t-1)^2/4)`.

**Proof.** Neither `n=1` nor prime `n` is possible, since those give `0`
and `1`. Let `p` be the least prime factor of composite `n`, and set
`q=n/p`. Then `q>=p`. If `q=p`, then `n=p^2` and `s(n)=1+p`, proving
the first alternative. If `q>p`, then `1,p,q` are three **distinct** proper
divisors. Therefore

\[
t\ge1+p+q\ge1+2\sqrt{pq}=1+2\sqrt n,
\]

which proves the second alternative. ∎

The prime-square branch must not be omitted: target `12` has the source
`121`, far beyond its ordinary bound `floor(121/4)=30`. The theorem also
gives the weaker universal bound `n<=(t-1)^2`. The program checks each
listed fiber both by the sharper theorem and, independently, by a divisor
sieve through that weaker universal bound.

For even `n` and `t>=2`, the two distinct proper divisors `1,n/2` also give
`n<=2(t-1)`. This is useful for hand checks but is not needed for completeness.

## 4. Complete inverse fibers

These are **all positive-integer sources**, not merely witness lists.

| Target `t` | Ordinary source bound | Prime-square exception | Complete fiber |
|---:|---:|---:|---|
| `5` | `4` | none (`4` not prime) | empty |
| `10` | `20` | none (`9` not prime) | `{14}` |
| `11` | `25` | none (`10` not prime) | `{21}` |
| `12` | `30` | `121=11^2` | `{121}` |
| `13` | `36` | none (`12` not prime) | `{27,35}` |
| `25` | `144` | none (`24` not prime) | `{95,119,143}` |
| `37` | `324` | none (`36` not prime) | `{155,203,299,323}` |
| `50` | `600` | none (`49` not prime) | `{40,94}` |
| `100` | `2450` | none (`99` not prime) | `{124,194}` |

For the small targets, completeness can be checked by hand. The composite
inputs through `36` have the following values:

```
n:     4  6  8  9 10 12 14 15 16 18 20 21 22 24 25 26 27 28 30 32 33 34 35 36
s(n):  3  6  7  4  8 16 10  9 15 21 22 11 14 36  6 16 13 28 42 31 15 20 13 55
```

Restricting this table to each proved ordinary bound gives the fibers for
`5,10,11,12,13`, with `121` added only for target `12`. All primes give `1`
and `n=1` gives `0`, so no omitted inputs qualify.

For distinct primes `p,q`, the exact semiprime identity is

\[
s(pq)=1+p+q.
\]

Hence the 37-fiber consists precisely of the products corresponding to

\[
36=5+31=7+29=13+23=17+19.
\]

The finite-bound check additionally excludes every other kind of source.
Likewise, the 25-fiber uses `24=5+19=7+17=11+13`. For the even targets,
`s(40)=50`, `s(94)=1+2+47=50`, `s(124)=1+2+4+31+62=100`, and
`s(194)=1+2+97=100`.

Target `5` is therefore untouchable by a complete elementary argument.
Nothing here classifies all untouchable numbers or proves a global assertion
about which odd targets have a preimage. Target `1` is excluded from the
finite-fiber theorem: every prime is a preimage of `1`.

## 5. Forward chains are a different operation

Repeated doubling and repeated application of `s` must not be conflated.
The exact forward chains are

```
  5 ->   1 -> 0
 10 ->   8 -> 7 -> 1 -> 0
 25 ->   6 -> 6 -> ...
 50 ->  43 -> 1 -> 0
100 -> 117 -> 65 -> 19 -> 1 -> 0
 74 ->  40 -> 50 -> 43 -> 1 -> 0.
```

In particular `25` reaches the fixed point `6`; it does not reach zero.
These explicit computations do not imply that arbitrary aliquot sequences
terminate or enter a cycle.

## 6. Exact links among 5, 13, 25, 37, 50, and 100

There are useful shared-factorization links, but their scopes should remain
visible:

1. **The carry 31 links 25 to a 37-preimage.** Since
   `sigma(25)=1+5+25=31` is prime and differs from `5`,
   \[
   s\bigl(5\sigma(25)\bigr)=s(155)=1+5+31=37.
   \]
   The same `31` is the carry in `25,50,100`. This is an exact special
   case of the semiprime identity, not an aliquot step from `25` to `155`.

   More generally, for an **odd prime** `p`, put
   \[
   q=\sigma(p^2)=s(p^3)=p^2+p+1,\qquad
   \kappa=2q-p^2=(p+1)^2+1.
   \]
   Here `q` is the carry and `kappa` is the odd slope on the doubling ray
   with odd part `p^2`. Since `p` does not divide `q`, multiplicativity gives
   the exact guarded-bridge formula
   \[
   \boxed{s(pq)-\kappa=(p+1)\bigl(s(q)-1\bigr).}
   \]
   Consequently `s(pq)=kappa` **if and only if `q` is prime**. The cases
   `p=3` and `p=5` give `(q,kappa,pq)=(13,17,39)` and `(31,37,155)`.
   The guard is essential: for `p=7`, `q=57=3*19`, `s(q)=23`, so
   `s(399)=65+8*(23-1)=241`, not the prime-only prediction `65`.

2. **The 100-chain really contains a 13-factor.** Factorization gives
   \[
   s(100)=\sigma(4)\sigma(25)-100=7\cdot31-100=117=9\cdot13.
   \]
   As `sigma(9)=13` and `sigma(13)=14`,
   \[
   s(117)=13\cdot14-117=65.
   \]

3. **The doubled 13-fiber joins the 37 lane.** The equality
   `s(35)=13` does not survive doubling as a scalar equality:
   \[
   35\xrightarrow{\times2}70\xrightarrow{s}74=2\cdot37
      \xrightarrow{s}40\xrightarrow{s}50\xrightarrow{s}43
      \xrightarrow{s}1\xrightarrow{s}0.
   \]
   Indeed `s(70)=sigma(2)sigma(5)sigma(7)-70=3*6*8-70=74` and
   `s(74)=1+2+37=40`.

No fourth-power representation theorem follows from these numerical links.
Transferring a conclusion to a fourth-power counting problem requires an
explicit common factorization, map, or preserved state; matching integers
alone are not such a transfer.

## 7. Reproduction and verification scope

Run from this directory:

```sh
python check_aliquot.py
python -O check_aliquot.py
```

The certificate compares three independent implementations (paired proper
divisors, prime-factor divisor sums, and a proper-multiples sieve) for every
`1<=n<=20000`; tests the source-bound alternatives on that interval; checks
the doubling law at `10000` inputs and the closed tower formula at `909`
inputs; reconstructs every source `1<=n<=10000` from its two observations;
checks the exact pair-image guard on all `7650` pairs with `0<=x<=50` and
`1<=y<=150` against an independent complete source sieve; certifies all nine
complete fibers with their proper-divisor witness lists; and verifies the
stated chains, lost-carry counterexample, and both modular-trace collisions.
The independent source-image sieve also checks two-observation injectivity
on every source through `45000`.
The guarded prime-square bridge and exact defect are checked at
`p=3,5,7,11,13,17,19`, including the composite-carry hostile example `p=7`.
For each dilation prime `p=2,3,5,7,11`, the checker additionally verifies
the carry law and reconstruction on `5000` sources, injectivity on `45000`
sources, and the arbitrary-pair guard on all `3100` pairs in
`0<=x<=30,1<=y<=100`. A composite multiplier `4` is explicitly rejected.

Checks use explicit exceptions, not `assert`, and therefore remain active
under `python -O`. The proofs establish general identities and completeness;
the finite computations independently check their stated arithmetic instances.
