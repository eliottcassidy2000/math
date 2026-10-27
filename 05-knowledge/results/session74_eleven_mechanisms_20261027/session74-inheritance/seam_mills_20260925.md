# Mills prime shells: exact remainder decoding and the odd-power doubling seam

**Status: PROVED elementary identities and decoder statements;
FINITE-EXACT prime/carry controls; CITED Mills-type existence and
arithmetic-type theorems.** No universal Collatz route, proof of
Legendre's conjecture, or unconditional identification of the familiar
least Mills decimal is claimed. No novelty claim.

## Inheritance and board

Targeted searches of current navigation, canon, the corrections ledger,
reference files, and research results found no prior Mills-specific
result. The actual local predecessor is the missing remainder coordinate
in [arithmetic seams, operations](arithmetic_seams_20260921_operations.md),
and the distinction between a new prime divisor and an actually prime
value in [arithmetic seams, primes](arithmetic_seams_20260921_primes.md).
The same-coordinate decoder obligation was made exact for Collatz in
[root decoding and Mahler comparison](decoder_mahler_catalan_20260925.md).
Its inherited canon is
[THM-2228, Mahler carry tail and integral stabilization](../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md)
and [THM-4072, Mahler safe-terminal fibre product](../../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md).

Closest proved mechanism: a quotient becomes an exact decoder only with
the remainder and domain predicate retained. Canonical hostile: an
integer cube of a prime is composite, while a nearby prime can occupy
many different positions in its power shell. Corrected near miss: the
existence of one prime-generating constant does not establish the least
constant's digits or a prime in every square interval. Least-used
sidecar: the fractional part before exponentiation, and the integer gap
digit it generates.

Anchor: define the prime-shell decoder. Niche: show what n+2 does to its
dyadic valuation. Wildcard: the genuine arithmetic distinction between
cubic and quartic least prime-representing constants.

| Concept | Retained quantity | Decisive boundary |
|---|---|---|
| Power shell | parent p, exponent c, child q | p^c itself is composite |
| Floor decoder | fractional part and gap digit | truncating the binomial changes primality |
| Odd n+2 chain | valuation of the horizontal difference | odd and even exponents behave differently |
| Doubles/higher layers | v2 of the gap, plus its odd part | gaps 4,16,24 all occur at certified primes |
| Mills existence | an infinite nested prime chain | not every prime, not every square interval |
| Least constant | global minimality and infinite extension | standard cubic decimal remains conditional |

## 1. Primary source scope (CITED)

The user supplied [Mills' constant](https://en.wikipedia.org/wiki/Mills%27_constant).
Its technical claims were checked against primary papers, rather than
used as theorem dependencies by themselves.

Mills' existence theorem is unconditional: some A>1 has
`floor(A^(3^n))` prime for every positive integer n. The familiar
identification of the least A with `1.3063778838...` is a different
claim. [Caldwell--Cheng, 2005, Theorem 1 and Section 3](https://cs.uwaterloo.ca/journals/JIS/VOL8/Caldwell/caldwell78.pdf)
give its digit computation under RH, using the least-prime-above-the-cube
chain starting from 2. A finite list of certified terms is unconditional;
its all-future extension from that small start is the remaining issue
for identifying this specific least constant by that route.

[Saito, *Mills' constant is irrational*, Mathematika 71 (2025),
arXiv v2](https://arxiv.org/pdf/2404.19461v2), printed page 1, defines
`xi_c=min{A>1: floor(A^(c^n)) prime for all n>=1}`, for integer c>=2.
Existence is Corollary 3.4. **Theorem 1.1:** every integer c>=4 gives
transcendental xi_c. **Theorem 1.2:** xi_3 is transcendental or some
`xi_3^(3^m)`, m>=1, is a degree-three Pisot number; in particular xi_3
is irrational. These results do not assume RH. Question 1.6 on page 2
asks whether xi_2 is rational or irrational. The conclusions concern
the **least** constant, not every possible prime-representing A.

[Matomaki, 2010, *Prime-representing functions*](https://users.utu.fi/ksmato/wp-content/uploads/sites/1432/2025/04/Primerepfunc.pdf)
proves existence even for `floor(A^(2^n))`, with continuum many such A.
This selects successful infinite square-shell branches. It does not
prove that every interval `(m^2,(m+1)^2)` contains a prime.

Thus cubes are neither the only usable exponent nor equivalent to
quartics in every respect: there is a precise published distinction
between the current arithmetic-type results for xi_3 and xi_c, c>=4.

## 2. Exact prime-shell tree and nested decoder (PROVED)

Fix an **integer** c>=2. Given primes p,q, call q a child of p when

```text
p^c < q < (p+1)^c.                                      (1)
```

Then `p=floor(q^(1/c))`, so the parent, when it is prime, is unique.
The gap digit is `d=q-p^c`. Neither p^c nor `(p+1)^c-1` can be prime:
the first is a proper power, and the second factors as

```text
(p+1)^c-1 = p * (1+(p+1)+...+(p+1)^(c-1)).
```

Consequently a prime child also satisfies `q+1<(p+1)^c`.

Suppose an infinite prime chain `(p_n)_(n>=1)` satisfies (1). Set

```text
I_n = [p_n^(1/c^n), (p_n+1)^(1/c^n)].                    (2)
```

The lower endpoints strictly increase and the upper endpoints strictly
decrease. The derivative of `x^(1/c^n)` on x>=1 is at most `1/c^n`,
so `length(I_n)<=1/c^n ->0`. Hence these nested closed intervals have
a unique common point A. Because `I_(n+1)` has upper endpoint strictly
below that of I_n, the limit satisfies

```text
p_n <= A^(c^n) < p_n+1,
p_n = floor(A^(c^n))                                     (3)
```

for every n. Conversely, any A producing primes as in (3) supplies (1):
raising `p_n<=A^(c^n)<p_n+1` to c and taking the next floor gives the
interval; equality at the lower endpoint would make p_(n+1)=p_n^c
composite. This is an exact correspondence between infinite prime-shell
chains and their represented constants.

This construction has actual ancestry and a reconstruction map. It
encodes **one branch**, not all primes or all integers. For example,
prime 67 has cube-root floor 4, which is composite, so it has no prime
parent in the cubic tree. The prime children of 2 in the cubic tree
are exactly `11,13,17,19,23`. Their shared parent does not determine
the next branch.

**Connection card.** Source: an infinite prime chain with (1).
Target: a single real constant. Map: nested intervals (2). Preserved:
every labelled prime and its index via (3). Lost if only the current
parent is kept: the chosen child and all future branch choices.
Sidecar: gap digits or nested intervals. Cheapest test: the five
different cubic children of 2. This is a model of successful exact
decoding, but it has no supplied map to universal Collatz ancestry.

## 3. Why prime gaps make cubes work, and what squares change

Assume a uniform large-x prime-gap input of the form: a prime lies in
`(x,x+K*x^theta)` for all sufficiently large x, with theta<1. Applying
it at x=p^c finds a child inside the next power shell whenever

```text
c*theta < c-1,             equivalently c > 1/(1-theta).
```

This follows by comparing `K*p^(c*theta)` with
`(p+1)^c-p^c-1`, whose leading term is `c*p^(c-1)`.
Saito's section 2, printed page 3, uses the cited Baker--Harman--Pintz
input theta=21/40, yielding the threshold `40/19`. Every integer c>=3
is on this side. The same uniform argument does not give c=2.
Matomaki's existence result obtains selected infinite branches by a
stronger selection argument; it does not change this inequality into
a proof that every square shell works.

## 4. The fractional remainder is a nonlinear carry (PROVED)

Write `A^(c^n)=p_n+theta_n`, where `0<=theta_n<1`. Then

```text
d_n = floor(sum_(j=1)^c binom(c,j) p_n^(c-j) theta_n^j),
p_(n+1) = p_n^c + d_n,
theta_(n+1) = (p_n+theta_n)^c - p_n^c - d_n.             (4)
```

For c=3 the digit is the floor of
`3p^2 theta + 3p theta^2 + theta^3`. This is the exact counterpart of
retaining ordered carries in a Collatz block. Exponentiation without
the fractional sidecar does not determine the prime child.

Even a last term smaller than one can change the decision at a floor
boundary. Take p=2 and theta=9/10. Then

```text
2^3+3*2^2*theta+3*2*theta^2 = 1183/50 = 23.66,
(2+theta)^3 = 24389/1000 = 24.389.
```

Dropping theta^3 reports prime 23, while the actual floor is composite
24. This is a local decoder hostile, not a claim that this theta comes
from an infinite Mills prime chain. Approximate terms may be omitted
only with a certified error smaller than the distance to the relevant
integer boundary.

Ordinary halving has its own exact remainder law:

```text
(p+theta)/2 = floor(p/2) + ((p mod2)+theta)/2.             (5)
```

For odd prime p the new fractional part is in `[1/2,1)`, and the new
integer is `(p-1)/2`, not necessarily prime. Certified controls are
`11 ->5` and `1361 ->680`. Thus a Mills prime predicate is not closed
under the halving required by Collatz. The shared technique is retention
of the residue and fractional coordinate; it is not an established
equivalence of the two dynamics.

## 5. The {2,3,11} step and the higher-power gap seam (PROVED)

For prime p and child prime q=p^c+d, parity gives

```text
d odd iff p=2;           d even for every odd parent p.  (6)
```

In particular, the only prime pair satisfying `q=p^3+3` is
`p=2,q=11`. This supplies a genuine exact occurrence of `{2,3,11}`:
parent, correction, child. Its mechanism is the exceptional even prime,
not the earlier odd-square-shell criterion. It is specific to this
exponent and correction: the first quintic step is `37=2^5+5`.

The doubled gap row after the first even parent is real, but its deeper
valuation layers cannot be discarded. Complete trial division gives

```text
347   = 7^3  + 4,
79531 = 43^3 + 24,
300779= 67^3 + 16.
```

Each is the **first** prime above its displayed cube. In particular,
the multiple-of-four sea contains necessary prime choices.

There are sharper algebraic exclusions that retain the missing odd
part. A **pure** gap `d=2^(3t)` is impossible, since

```text
p^3+2^(3t)=(p+2^t)(p^2-p*2^t+2^(2t)),       t>=1.
```

Both factors exceed one. This does not ban valuation v2(d)=3: the
certified gap 24 above is allowed. At the quartic scale the specific
gap four is impossible by

```text
p^4+4=(p^2-2p+2)(p^2+2p+2),                p prime.
```

Thus exponent shape and the complete gap matter. A valuation-only
quotient cannot identify which candidates remain prime.

## 6. n+2 under odd powers stays on the doubles seam (PROVED)

For odd integer p>=1 and positive integer c,

```text
v2((p+2)^c-p^c) = 1                         if c is odd,
                 v2(c)+v2(p+1)+1            if c is even. (7)
```

For odd c, factor the difference by `(p+2)-p=2`. The remaining sum
has c odd summands and is odd. For even c write c=`2^s*d`, d odd.
The odd d factor does not change the valuation. The first squaring
contributes `(p+2)^2-p^2=4(p+1)`; every subsequent difference-of-squares
factor is the sum of two odd even powers, congruent to 2 modulo 8,
and adds one factor of 2. The total is `2+v2(p+1)+(s-1)`.

Consequently cubes and fifth powers take the odd n+2 chain to
differences of valuation exactly one. Squares give valuation at least
three; quartics give valuation at least four. This is an exact way the
same horizontal recursion enters different dyadic layers under a change
of representation. It concerns differences between power-shell bases;
it is distinct from the chosen prime gap d in (4).

## 7. Exact controls and stopping boundary

Run the [checker](../../04-computation/experiments/seam_mills_20260925.py):

```text
python 04-computation/experiments/seam_mills_20260925.py
python -O 04-computation/experiments/seam_mills_20260925.py
```

The [frozen output](seam_mills_20260925.out) agrees under both modes.
All primality decisions are complete trial division using a sieve
through 100000, valid for the declared universe n<=10^10. The sieve
is independently checked by divisor-by-divisor trial division through
100000; displayed cubic-prefix primes are independently trial divided
again. No probable-prime test or floating point is used.

The script checks 114 next-prime power shells with prime parent <=97,
integer exponent 2..6, and upper bound 10^10. Finite greedy prefixes are

| c | Certified prefix |
|---|---|
| 2 | 2,5,29,853,727613 |
| 3 | 2,11,1361,2521008887 |
| 4 | 2,17,83537 |
| 5 | 2,37,69343979 |
| 6 | 2,67 |

Any constant realizing the displayed cubic prefix lies in the exact
outer rational bracket

```text
1306377883863080 / 10^15 < A < 1306377883869479 / 10^15.
```

This finite-prefix bracket does not unconditionally identify the least
infinite cubic constant. Other controls cover 4608 rational binomial/
halving states, 5489 adjacent odd-power valuation cases, 150 cubic
factorizations, and 25 quartic factorizations.

The board after the probe: Mills provides a clean example of an exact
branch decoder whose fractional carry survives compression. The odd
chain has a provable valuation law, and the prime gap seam retains
essential higher-power layers. The remaining Collatz obligation is
still a map from tournament structure to a valid arithmetic certificate;
neither prime-shell existence nor exponentiation symmetry supplies it
without that map.
