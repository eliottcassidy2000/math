# Prime-power bias with its residue map retained

**Status.** PROVED: the finite residue-power transform and its inverse, the
exact mod-four square correction, its remainder bound, and the conditional
transfer of normalized logarithmic means. FINITE-EXACT: the complete first
lead through 26,861 and the declared arithmetic controls. CITED: classical
prime-race oscillation and the conditional logarithmic-density theorem.
No new analytic prime-distribution theorem or Collatz conclusion is claimed.

The first strict lead of primes `1 mod4` is at **26,861**, where the counts
are **1,473 versus 1,472**. Thus 1,473 is a class count, not the location of
the crossover. The useful reusable operation is a lossless conversion
between prime counts and weighted prime-power counts. Squaring alone loses
the prime's residue class; the full conversion retains the exponent and the
smaller argument at which each residue vector is read.

## 1. Inheritance, hostile, and proposed connection

The closest repo mechanism is Theorem 2.2 of
[signed sandwich discrepancy](collatz_mod6_20260917_sandwich_bias.md): an
unordered semiprime count equals half its ordered-pair count plus half the
prime-square diagonal. That note explicitly retains a separate pairing
residual; its prime/semiprime marginals do not determine the sign of the
original adjacent-pair count. Its pointwise claim first fails at `K=1453`.
That historical number is distinct from the present **1473** class count.

The board here is **prime base / power exponent / residue map / cutoff /
counting weight / analytic remainder**. The least-used sidecar is the map
`b -> b^k` on *all* residue classes, including nonunits. The cheap hostile is
`3^2=9`: class 3 contributes to the weighted count of class 1. Forgetting this
move breaks even the smallest nontrivial inverse calculation.

Source: the full residue count functions. Target: the full weighted
prime-power count functions. Map: (1) below. Preserved: every prime-count
value, once all needed lower cutoffs are retained. Lost by a single squared
snapshot: base residue and exponent. Restoration: the power-map and cutoff
indices in the finite Möbius inverse. This is classical prime-power
bookkeeping organized as an exact reusable interface, not a novelty claim.

## 2. An invertible transform despite noninvertible power maps

For an integer modulus `q>=2` and every residue `a`, define

\[
\pi_q(x,a)=\#\{p\le x:p\text{ prime},\ p\equiv a\pmod q\},
\qquad
J_q(x,a)=\sum_{\substack{p^k\le x\\p^k\equiv a\pmod q}}\frac1k.
\]

Here `p` is prime and `k>=1`. A prime power has its full weight at the
endpoint; no half-jump convention is used. Put both functions equal to zero
for `x<2`. On a residue vector `v`, let

\[
(P_kv)_a=\sum_{b^k\equiv a\ (q)}v_b.
\]

The maps can collapse classes, but `P_k P_l=P_(kl)`. Grouping prime powers by
their exponent gives the exact finite identities

\[
J_q(x)=\sum_{k\ge1}\frac1k P_k\pi_q(\lfloor x^{1/k}\rfloor),
\qquad
\pi_q(x)=\sum_{k\ge1}\frac{\mu(k)}k P_kJ_q(\lfloor x^{1/k}\rfloor).
\tag{1}
\]

Only `k<=floor(log_2 x)` can contribute. To prove the inverse, substitute the
first formula into the second. The coefficient at total exponent `m=kl` is
`(1/m) sum_(k|m) mu(k)`, which is 1 for `m=1` and zero otherwise. Integer
roots compose exactly: `floor(floor(x^(1/k))^(1/l))=floor(x^(1/(kl)))`.
Neither cancellation nor this root identity requires `b` to be a unit or a
power map to be injective.

The scope is the family of functions and their finite lower-cutoff queries.
A single vector `J_q(x)` is not asserted to determine all the unknown
lower-cutoff data. In particular, at `q=4,x=9`, using the scalar formula
without residue relabelling would give

\[
J_4(9,1)-\tfrac12J_4(3,1)=\tfrac32\ne\pi_4(9,1)=1.
\]

The correct square preimage also includes class 3 and subtracts
`J_4(3,3)/2=1/2`. Likewise `P_2 e_1=P_2 e_3`: the square map itself cannot
recover the base class. The inverse works across scales, not by inverting
that collapsed map.

## 3. Exact squares, and the term that must not be dropped

Write

\[
\Delta(x)=\pi_4(x,1)-\pi_4(x,3),\quad
D_J(x)=J_4(x,1)-J_4(x,3),\quad
\pi_o(y)=\#\{p\le y:p\text{ odd prime}\}.
\]

Odd powers preserve both odd classes; every even power lands in class 1.
Consequently

\[
D_J(x)=\Delta(x)+\frac12\pi_o(\sqrt x)+R(x),\tag{2}
\]

\[
R(x)=\sum_{\substack{k\ge3\\k\text{ odd}}}\frac{\Delta(x^{1/k})}{k}
 +\sum_{\substack{k\ge4\\k\text{ even}}}\frac{\pi_o(x^{1/k})}{k}.
\tag{3}
\]

All sums remain finite and exact. The script reports each summand as a
rational, together with the rigorous error bill

\[
|R(x)|\le B(x):=\sum_{k\ge3}\frac{\pi_o(x^{1/k})}{k}
\le x^{1/3}\log_2 x\qquad(x\ge2).\tag{4}
\]

Thus `D_J-square-B <= Delta <= D_J-square+B` is a valid interval when only
the aggregate bill is used. The omitted term is not identically zero:
`R(27)=-1/3`, whereas `R(1000)=2/15`.

The classical prime number theorem gives
`pi_o(sqrt x)/2 ~ sqrt(x)/log x`; (4) gives
`R(x)=o(sqrt(x)/log x)`. These statements identify the deterministic square
correction. **They do not give**
`pi_4(x,3)-pi_4(x,1) ~ sqrt(x)/log x`, because `D_J` is still present.

One can make the distinction exact. For `x>=2`, put

\[
E(x)=\frac{\log x}{\sqrt x}\,[\pi_4(x,3)-\pi_4(x,1)],\qquad
W(x)=\frac{\log x}{\sqrt x}D_J(x).
\]

Then (2)--(4) imply `E(x)+W(x)->1`. If the logarithmic mean
`(log X)^(-1) integral_2^X W(x) dx/x` exists and equals `m`, the same mean
of `E` exists and equals `1-m`. Proof: the logarithmic average of a function
tending to 1 tends to 1, by separating a fixed initial interval from its
tail. In particular a *zero mean* for `W` yields *mean one* for `E`, not
pointwise convergence of `E` to 1. The analytic zero-mean premise is not
proved by this finite checker.

## 4. The first reversal is a source-owned exact calculation

The sieve is independently checked by trial division for every integer
from 0 through 30,000. Scanning every prime up to 26,861 proves there is no
earlier strict lead; between consecutive primes the counts are constant.

| Cutoff `x` | `pi_4(x,1)` | `pi_4(x,3)` | Difference |
|---:|---:|---:|---:|
| 26,860 | 1472 | 1472 | 0 |
| 26,861 | 1473 | 1472 | 1 |
| 26,862 | 1473 | 1472 | 1 |
| 26,863 | 1473 | 1473 | 0 |
| 26,879 | 1473 | 1474 | -1 |

The prime 26,861 is the **2,946th prime overall**, including 2. At that
specific source the exact decomposition is

\[
\frac{50831}{2520}=1+\frac{37}{2}+\frac{1691}{2520}.
\]

The positive prime difference is not a negative square correction: the
weighted race has overcome that correction at this cutoff. At 100,000 the
same quantities are `Delta=-25`, square term `32`, and
`D_J=3029/360>0`. Weighted and unweighted race signs need not even agree.

`RaceReceipt` retains the actual cutoff and both count conventions.
`audit_receipt(packet,expected_bound)` rejects a valid packet borrowed from
another cutoff, even when its unweighted counts happen to be unchanged.
The production cap `x<=10^6` is explicit; the proofs of (1)--(4) have no such
cap. A failed finite `first_lead(limit)` query returns `None`, not a
universal sign theorem.

## 5. Connection to the earlier mod-six square diagonal

For either two-class unit group modulo 4 or 6, let `chi` have values `+1,-1`.
The same identity `chi(p^2)=1` appears in two different operations:

* In (2), one removes weighted prime powers to recover primes. The square
  term is subtracted from the weighted race.
* In the inherited unordered semiprime identity, with primes coprime to the
  modulus, `2 sum_(p<=q,pq<=x) chi(pq) = sum_(p,q,pq<=x) chi(p)chi(q)
  + pi_unit(sqrt x)`. The square term repairs the diagonal's multiplicity.

The checker independently enumerates both sides for `x<=200`. This common
diagonal mechanism is exact; it does not identify prime-power counts with
semiprime counts or fix the sign of the off-diagonal contribution.

Changing modulus also needs its excluded-prime sidecar. For `x>=2`,

\[
\pi_6(x,1)-\pi_6(x,5)
=\pi_3(x,1)-\pi_3(x,2)+1,
\]

because modulo 6 the negative class omits prime 2. There is no justified
map here from aggregate prime races to a prescribed Collatz source's native
valuation word, endpoint or ROOT certificate. The transferable method is
retaining the exact finite correction before interpreting an averaged bias.

## 6. Analytic status and primary sources

[Rubinstein--Sarnak, *Chebyshev's bias* (1994), authored publication record](https://collaborate.princeton.edu/en/publications/chebyshevs-bias/)
states the GRH and Grand Simplicity hypotheses. Under GRH plus the relevant
linear-independence hypothesis, the logarithmic density of
`pi_4(x,3)>pi_4(x,1)` is approximately `0.9959`. This weights cutoffs by
`dx/x`; it is not a claim about every cutoff.

[Granville--Martin, *Prime Number Races*, sections “The Mod 4 Race” and
“Where do these biases come from?”](https://arxiv.org/html/math/0408319v1)
explains the square correction and reports the classical oscillation:
both signs occur arbitrarily far out, on the normalized scale as well.
That rules out the proposed positive pointwise asymptotic. The
[authors' survey version](https://personal.math.ubc.ca/~gerg/teaching/613-Winter2011/PrimeNumberRaces.pdf)
explicitly reports the counts 1473 and 1472 at 26861; the computation here
reproduces them independently. No claim about Mertens sums or the
`pi(x)-Li(x)` Skewes race is imported into the mod-four identity.

## 7. Reproduction and scope

```text
python -B 04-computation/experiments/prime_power_bias_20261014d.py
python -B -O 04-computation/experiments/prime_power_bias_20261014d.py
```

Both runs reproduce **48,915 exact checks** and the saved
[JSON](prime_power_bias_20261014d.json), LF SHA256
`49ad4ecaacf187d4adaf0dbabc6fa7a49a38604f4cb98cbf4966db97b3a2392e`.
Declared universe: every `x=0..512`, every modulus `q=2..12`, all residue
classes, both transforms; power-semigroup controls for `k,l=1..7`; complete
trial division through 30,000; the finite crossing; and exact large snapshots
through `10^6`. Hostiles include omitted residue preimages, omitted higher
powers, wrong source, booleans/floats, a bounded query misread as universal,
and weighted/unweighted sign disagreement. No approximate roots enter any
check. Literature supplies the analytic claims; finite output supplies none
of their hypotheses or asymptotic conclusions.
