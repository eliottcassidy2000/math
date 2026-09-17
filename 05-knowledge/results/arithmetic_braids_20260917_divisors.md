# Arithmetic braids: divisor balance and almost-prime sandwiches

**Status:** PROVED elementary classifications and local identities below;
FINITE-EXACT positive-center census; CITED Chen theorem. No claim here proves
twin primes, Collatz, Goldbach, or LRC(14). This is a research result note,
not a newly reserved canon theorem ID.

**Session:** arithmetic-braids-20260917, divisor/sandwich lane.

## Inheritance and live board

The closest proved mechanism is the multiplicative divisor fibre in
[THM-2422, operation fibres and twin-center ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
especially its equations (5)--(8a). Its canonical hostile is `4`: deleting
the repeated-factor diagonal makes a prime square look prime. The corrected
near miss is the exceptional twin center `4` in
[MISTAKE-268](../../01-canon/MISTAKES.md): only nonexceptional twin centers
belong to `6Z`. The least-used sidecar relevant here is the **whole prime
exponent profile**, rather than its sum or support. The related
[THM-2500, finite-hole relative atoms](../../01-canon/theorems/THM-2500-finite-hole-relative-atoms-prime-rays-and-cover-bounds.md)
likewise keeps valuation boxes when factor shadows lose witnesses.

The live concepts are:

| Object | Retained coordinate | Cheap hostile or decisive probe |
|---|---|---|
| Proper divisor lattice | exponent box and its Boolean subcube | `p`, `p^2`, `p^3`, `p^2q`, `p^2qr` |
| Fixed divisor defect | support size and multiplicative partitions | four distinct primes with one repeated |
| Ordered sandwich `(6k-1,6k+1)` | two separate Omega values | `k=4`, which breaks exact mixed symmetry |
| Local prime sieve | left/right/neutral residue label | complete residue systems mod `5,35,385` |
| Genuine prime predicate | unexposed cofactor and positive size | compare local reflection to positive prefix census |

The first two give exact classifications; the last three explain which part
of the apparent three-way pattern is an exact local law and which part needs
additional arithmetic information.

## DB1. The user's divisor equality has exactly the three claimed shapes

For `N>=2`, define

```text
F(N) = #{d: d|N, 1<d<N},
S(N) = #{d: d|N, 1<d<N, d squarefree},
U(N) = #{p: p|N, p prime, p<N}.
```

These definitions are necessary to reproduce the example `F(p)=S(p)=U(p)=0`.
In particular `U` is not the usual `Omega`, counting multiplicity, and at a
prime it is not even the usual distinct-prime count `omega(p)=1`.
Including `1` or `N` in `S` changes the equation. If `N=1` is admitted,
all three sets are empty and `1` is an additional trivial solution.

Every divisor counted by `U` is counted by `S`, so `S>=U`.

Write `N=prod_(i=1)^r p_i^{a_i}` with distinct primes and positive exponents.
Then

```text
F = prod_i(a_i+1)-2,
S = 2^r-1-[N squarefree],
U = r-[N prime].                                             (DB1)
```

Consequently, for `N>=2`,

```text
F(N)=S(N)+U(N)
iff N=p, or N=p^3, or N=p^2 q r,
```

where the primes in the last expression are pairwise distinct.

**Proof.** A squarefree composite has `F=S` and `U>0`, so fails. Primes
give `0=0+0`. If `N` is nonsquarefree, equality becomes

```text
prod_i(a_i+1)=2^r+r+1.                                      (DB2)
```

For `r=1`, this is `a_1+1=4`, giving `p^3`. For `r=2`, the right side is
`7`, which is not a product of two integers at least `2`. For `r=3`, the
right side is `12`; its only factorization into three integers at least
`2` is `3*2*2`. For `r>=4`, nonsquarefreeness forces

```text
prod_i(a_i+1)>=3*2^(r-1)>2^r+r+1,
```

because `2^(r-1)>r+1` at `r=4` and the difference increases afterward.
This proves both directions.

**Mechanism.** Divisors form an exponent box `prod_i {0,...,a_i}`;
squarefree divisors form its Boolean subcube `{0,1}^r`. For a nonsquarefree
number, `F-S` counts the box points outside that cube after removing `N`
itself. Equality says that this excess equals `r`, the number of prime
directions. This is a small-dimensional exponential balance, not an
unproved extrapolation from three numerical examples.

The examples are `p`: `0=0+0`; `p^3`: `2=1+1`; `p^2qr`: `10=7+3`.
The nearby failures `p^2`: `1!=1+1` and `p^2q`: `4!=3+2` locate the
boundary.

## DB2. Sharp defect and all-scale finite classification

Set `D=F-S-U`. At a nonsquarefree number with `r` distinct prime factors,

```text
D=prod_i(a_i+1)-2^r-r-1 >= 2^(r-1)-r-1.                    (DB3)
```

Equality in the inequality holds exactly at exponent profile
`(2,1,...,1)`, up to permutation. The proof is that every exponent factor
is at least `2`, and at least one is at least `3`; equality in that product
bound forces exactly the stated profile.

The nonsquarefree minima at support sizes `1,2,3,4,5,6` are
`-1,-1,0,3,10,25`. Squarefree composites have `D=-r`; primes have `D=0`.
Thus the only nonsquarefree numbers with `D<=0` are

```text
p^2, p^3, p^2q, p^2qr,
```

with distinct primes in the latter two patterns. The zeros are the DB1
patterns; `p^2` and `p^2q` have defect `-1`.

For every fixed integer `d`, there are only finitely many prime exponent
profiles with `D=d`, and they can be enumerated without a bound on `N`:

1. Include profile `(1)` if `d=0`; include the squarefree profile of length
   `-d` if `d<=-2`.
2. For each support `r` satisfying `2^(r-1)-r-1<=d`, factor
   `2^r+r+1+d` into `r` unordered integers `b_i>=2`.
3. Keep nonsquarefree profiles `a_i=b_i-1`.

The support bound is finite because its left side tends to infinity;
for each support the factorization list is finite. The companion implements
this complete inverse problem for `-6<=d<=12`. For example,

| Defect | All shapes, distinct prime letters |
|---|---|
| `-1` | `p^2`, `p^2q` |
| `0` | `p`, `p^3`, `p^2qr` |
| `1` | `p^4`, `p^3q` |
| `2` | `p^5`, `p^2q^2` |
| `3` | `p^6`, `p^4q`, `p^2qrs` |
| `4` | `p^7`, `p^3qr` |

This turns the observed isolated equality into an effective family of
classifications. No analytic number theory is needed.

## DB3. Almost-prime notation needs multiplicities and disjoint support

If `A,B,C` mean integers with exactly `1,2,3` prime factors counted with
multiplicity, then `p^2qr` is a **4-almost-prime**. In set-product notation
the three solution types can be written

```text
A,  {a^3:a in A},  {a^2 b:a in A, b=qr in B, q!=r, gcd(a,b)=1}.
```

The last two restrictions matter: unrestricted `A^2B` also contains
`p^4` and `p^2q^2`, both non-solutions. Literal `C^2B` consists of
8-almost-primes, since `Omega(c^2b)=2*3+2=8`, whether or not their supports
overlap. It cannot represent `p^2qr`.

This is the same information-loss mechanism as in the inherited divisor
fibre: an Omega class forgets which prime was repeated. Neither Omega alone
nor omega alone determines `D`.

## SW1. The sandwich object is an ordered matrix, with an exact local law

Let `A_j={n>=2:Omega(n)=j}` where Omega counts prime factors with
multiplicity, and define

```text
N_ij(K)=#{1<=k<=K:6k-1 in A_i, 6k+1 in A_j}.               (SW1)
```

The four prime/semiprime cells form the product `{1,2} x {1,2}`. Adding
3-almost-primes gives a `3 x 3` product. These are outcome tables. There is
no intrinsic pairwise winner relation, so they are not tournaments. Ranking
their frequencies would manufacture a transitive comparison and introduce
ties; it does not recover the discarded arithmetic.

There is, however, an exact source for a local **three-way split**. Let `Q`
be a finite set of primes at least `5`, `M=prod_(p in Q) p`, and let `h_-`
and `h_+` count the Q-primes dividing `6k-1` and `6k+1`, respectively.
Then the following polynomial identity is PROVED:

```text
sum_(k mod M) x^{h_-(k)} y^{h_+(k)}
   = prod_(p in Q)(p-2+x+y).                               (SW2)
```

**Proof.** Modulo each `p`, invertibility of `6` gives exactly one
left-hit residue, one distinct right-hit residue, and `p-2` neutral
residues. The Chinese remainder theorem makes these local choices
independent on the complete residue system modulo `M`; multiplying their
generating polynomials counts every residue once.

Thus the polynomial is a function of `x+y`. Its coefficient at `x^a y^b`
is `binom(a+b,a)` times the coefficient at `z^{a+b}` in
`prod_(p in Q)(p-2+z)`. Under the uniform measure on residues modulo `M`,

```text
E h_-=E h_+=sum_(p in Q) 1/p,
Cov(h_-,h_+)=-sum_(p in Q) 1/p^2.                          (SW3)
```

The covariance follows because one prime cannot divide both endpoints:
`gcd(6k-1,6k+1)=1`. Cross-prime contributions cancel by CRT independence;
at each prime the covariance is `0-(1/p)^2`.

This gives a rigorous microcosm-to-macrocosm transfer: the three labelled
local residue states generate the entire finite-prime joint distribution.
It also demonstrates why treating endpoints as independent is incorrect.

**Connection contract.** The source is the labelled local divisibility
state at each `p in Q`; the target is the exact joint count over `Z/MZ`;
the map is CRT; the preserved predicates are divisibility and side of the
center. The quotient forgets prime powers, factors outside `Q`, and the
positive magnitude of the integer representative. The restoration sidecar
is those valuations and the remaining cofactor. The cheapest decisive
tests are the exact tables at `Q={5}`, `{5,7}`, `{5,7,11}`; all pass in the
companion. This does not prove any prime or semiprime asymptotic.

## SW2. What symmetry survives for actual integers

Define Omega of a negative integer using its absolute value. The map
`k -> -k` exchanges the two endpoint types, because

```text
|6(-k)-1|=6k+1,   |6(-k)+1|=6k-1.
```

So the signed-center census on `-K<=k<=K, k!=0` is exactly symmetric.
Likewise `k -> -k` proves the symmetry in SW2. It does not preserve the
positive prefix `1<=k<=K`, and therefore cannot prove `N_12(K)=N_21(K)`.

The first failure is `K=4`: the centers `6,12,18` have prime pairs and
center `24` has `(23,25)`, a prime/semiprime pair. Thus
`N_12(4)=1`, `N_21(4)=0`. This is an explicit failure of exact equality,
not a refutation of possible asymptotic relations.

A useful missing coordinate is the number `b(n)` of prime factors
congruent to `5 mod 6`, counted with multiplicity. Since both endpoints
are coprime to `6`,

```text
n mod 6 = (-1)^{b(n)},
b(6k-1) is odd,   b(6k+1) is even.                        (SW4)
```

In particular a semiprime on the left has one prime of each reduced
residue class modulo `6`; a semiprime on the right has two in the same
class. A 3-almost-prime on the left has either one or three `5 mod 6`
factors, while on the right it has zero or two. This orientation sidecar
is completely invisible in the bare Omega table.

## SW3. Exact census and the Chen boundary

For the unfiltered universe of **1,000,000 centers** `W=6,12,...,6,000,000`,
the companion gives the following FINITE-EXACT table. Rows describe `W-1`;
columns describe `W+1`.

| Omega left / right | `1` | `2` | `3` | `>=4` |
|---|---:|---:|---:|---:|
| `1` | 37,915 | 78,689 | 59,706 | 30,192 |
| `2` | 78,277 | 157,420 | 112,416 | 52,182 |
| `3` | 59,992 | 112,305 | 72,363 | 28,755 |
| `>=4` | 30,161 | 52,125 | 28,736 | 8,766 |

The first `3 x 3` block covers `769,083` centers; `230,917` have at least
one endpoint with four or more prime factors. Probabilities on this finite
universe are exactly these entries divided by `1,000,000`; conditioning
on the displayed `3 x 3` block uses denominator `769,083` instead.

The four prime/semiprime cells have the suggested order at this scale,
except the mixed counts are unequal. Their difference is `412`. The
ordering is not universal across finite scales: through center `600`,
the four counts `(PP,PS,SP,SS)` are `(26,28,20,13)`, so semiprime/semiprime
is then the smallest of the four. These are observations at stated finite
scales, not asymptotic theorems.

**CITED literature audit.** Chen's original 1973
[*On the representation of a larger even integer as the sum of a prime
and the product of at most two primes*, Theorems I and II, pp.157--158](https://www.scribd.com/document/831214629/Chen-Prime-Paper)
states a sufficiently-large-even-integer prime-plus-at-most-two-primes
result and a shifted-prime result: for fixed even `h`, infinitely many
primes `p` have `p+h` with at most two prime factors, accompanied by lower
bounds. Theorem II with `h=2` concerns the union of prime and exact
semiprime outcomes. It does not give exact probabilities, a frequency
ranking, or equality of the two mixed cells. The primary research paper
by [Pintz, *An approximation to the twin prime conjecture and the parity
phenomenon*, introduction](https://arxiv.org/pdf/1004.1065)
explicitly explains that this union does not decide which exact
factor-count outcome occurs infinitely often.

There is an additional carrier detail: the bare shifted-prime statement
does not itself restrict centers to `6Z`, because `p=13`, `p+2=15=3*5`
has center `14`. Transport to the `6k` table requires an appropriately
rough or congruence-restricted formulation; this note does not import that
stronger formulation without a separate source audit. Likewise the twin
pair `(3,5)` has exceptional center `4`, excluded from SW1 by definition.

## Reproduction and audit scope

```bash
python 04-computation/experiments/arithmetic_braids_20260917_divisors.py
python -O 04-computation/experiments/arithmetic_braids_20260917_divisors.py --output /tmp/arithmetic_braids_divisors_optimized.json
```

On Windows replace the optional output path by a suitable absolute path.
The default output is
[the matching JSON](../../04-computation/experiments/arithmetic_braids_20260917_divisors.json),
which embeds the source SHA-256. The source is
[the dependency-free companion](../../04-computation/experiments/arithmetic_braids_20260917_divisors.py).
The ordinary and optimized Python runs both passed and produced byte-identical
LF-normalized JSON. Recorded SHA-256 values are

```text
source: b7a5161f87c6d1784e917678c455a587838a77651bd9356ad1c8a3c418dc5aab
output: 0491d8b49aafd3aab86d2ccd39a032343455940a3dd584c37ee02eb26ca7e60e
```

The exact audits comprise all `24,309` sorted exponent profiles of support
`1..9` with exponents `1..8`; direct divisor-set enumeration for every
`2<=N<=10,000`; a prime-power sieve through `6,000,001`; independent trial
factorization of its Omega values through `10,000`; an independent complete
ordered sandwich matrix through `k=1000`; all million unfiltered centers;
the modulo-six factor-parity law through `k=10,000`; and exact CRT tables
modulo `5,35,385` checked against the generating polynomial and covariance.
All load-bearing checks use explicit exceptions and remain active under
Python `-O`. Universal claims in sections DB1--DB3 and SW1--SW2 rest on the proofs above,
not on the finite profile or integer bounds.

**Stopping boundary / new questions.** The prime-exponent profile completely
solves the divisor-balance lane. The sandwich lane now has a precise local
transfer and a quantified loss ledger. A next substantive step would study
the signed discrepancy `N_ij(K)-N_ji(K)` while retaining the residue-factor
parity and the unsieved cofactor, or introduce a verified roughness window
before comparing with Chen's sieve. Repeating the same finite table at a
larger cutoff without a new invariant would not cross that boundary.
