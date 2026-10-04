# Three inverse divisions free joint ternary and19-adic home addresses

2026-10-04. **PROVED** for the construction, lifting law, and sharp uniform
odd-rank bound. **FINITE-EXACT** for the declared address and route controls.
Collatz convergence for an arbitrary supplied integer remains **OPEN**.
There is no historical priority claim for inverse words, valuation lifting,
or the Chinese remainder theorem.

[Script](../../04-computation/experiments/joint_home_address_20261004.py)
and [saved output](joint_home_address_20261004.out).

## 1. Inheritance and the changed coordinate

The [nineteen-ray bank](mod19_route_lifts_20261004.md) constructs a completed
route in every19-adic residue, with sharp uniform odd rank four. The
[joint one-row selector](mod19_recursive_observers_20261004.md), sections4--5,
instead fixes one inverse hub and row. Its ternary and19-adic addresses must
agree on their shared block digit; separately reachable addresses may be
incompatible. That theorem and its5/341 hostile remain valid.

The [shared-depth invariant](mod19_resonance_depth_20261004.md) identifies
the common clock factor, and the [torus map](inverse_ray_torus_clock_20261004.md)
retains the hub, orbit and direction quotient. The least-used coordinate here
is the number of inverse divisions between a changing exponent and the source.
Moving the exponent to the third inverse position changes which ternary
digits survive. It does not relax a guard on a fixed one-edge channel.

The board is **inverse depth / parameter period / division precision /
joint address / exact source / completed suffix**. The anchor is a checked
certificate constructor; the niche is two-prime address lifting; the wildcard
is using a denominator as an information channel. The canonical hostile is
still a different integer with the same finite residue data.

## 2. Uniform joint-address theorem

**PROVED.** For every a>=0, k>=1, every x modulo3^a, and every y modulo19^k,
there are infinitely many positive odd integers n such that

    n=x mod3^a,       n=y mod19^k,

with an explicit first-hit Collatz route to1 of at most four odd steps.
Four is the smallest uniform bound, even when a=0.

Write U(n)=oddpart(3n+1). Begin with an actual three-letter word
`(b1,b2,b3)` from a positive odd n0 to a fixed hub h, and put

    E=b1+b2+b3,
    A=h*2^E,
    B=9+3*2^b1+2^(b1+b2),
    n_t=(A*2^(18t)-B)/27,          t>=0.                 (1)

Assume every point on the three-edge prefix exceeds1 except possibly its
terminal hub, and assume h is prime to3 and19. The last exponent is absent
from B. Since `2^18=1 mod27`, all three inverse divisions remain integral.
Every inverse value increases with t, so positivity and the no-early-root
guard survive. The exact word at n_t is

    (b1,b2,b3+18t),

followed by the supplied first-hit route of h. Oddness is preserved at each
inverse step. No forward orbit search is needed to construct the certificate.

The following finite witness table chooses one ray for each source residue
modulo19. Each row is checked by three literal inverse divisions. The hub5
has its one-edge certificate5->1; the other hubs are already1.

| n0 mod19 | n0 | Prefix | Hub |
|---:|---:|---|---:|
|0|2261|(7,5,4)|1|
|1|9045|(9,5,4)|1|
|2|144725|(13,5,4)|1|
|3|17749|(12,3,4)|1|
|4|2417|(2,6,8)|1|
|5|36181|(11,5,4)|1|
|6|1507|(1,7,5)|5|
|7|1109|(8,3,4)|1|
|8|141|(3,5,4)|1|
|9|4835|(1,8,8)|1|
|10|4437|(10,3,4)|1|
|11|277|(6,3,4)|1|
|12|69|(4,3,4)|1|
|13|70997|(14,3,4)|1|
|14|565|(5,5,4)|1|
|15|283989|(16,3,4)|1|
|16|35|(1,5,4)|1|
|17|17|(2,3,4)|1|
|18|75|(1,2,8)|1|

These witnesses were selected from the complete box b1,b2,b3 in1..18;
no claim of globally smallest witnesses is made. The critical ray differs
from the earlier bank's parameter placement:

    n_t=(40960*2^(18t)-271)/27,
    word (1,7,5+18t,4),       n_0=1507.               (2)

Its variable is now the third source-side exponent. That position matters.

## 3. Two independent digit lifts

For distinct parameters t,s, subtracting (1) gives

    v3(n_t-n_s)=v3(t-s),
    v19(n_t-n_s)=1+v19(t-s).                           (3)

Indeed A is a unit at both primes,
`v3(2^18-1)=3`, and `v19(2^18-1)=1`. For an odd prime q and q|z-1,
binomial expansion first gives `v_q(z^q-1)=v_q(z-1)+1`; an exponent
prime to q leaves that valuation unchanged. Factoring t-s into a q-power
and a unit proves the needed lifting formula. Division by27 removes
exactly the three ternary powers and no factor19.

Consequently t modulo3^a bijects with all source residues modulo3^a,
while t modulo19^(k-1) bijects with the19-adic disk fixed by n0 modulo19.
The two parameter moduli are coprime. Each requested pair therefore has
exactly one parameter in

    0 <= t < 3^a*19^(k-1).                            (4)

Adding any nonnegative multiple of this period supplies the infinitely many
distinct positive certified representatives. This proves the theorem.

The digit compiler uses explicit nonzero derivatives. At ternary level a,
adding d*3^a to t changes the next source digit by `d*A mod3`. At19-adic
level k, adding d*19^(k-1) changes the next source digit by `d*A/9 mod19`.
The identities follow from

    (2^18-1)/27=9709=1 mod3,
    (2^18-1)/19=13797=3 mod19.

Both lifts retain division precision: the residue evaluator first works
modulo27 times the requested modulus, then divides the exact numerator
by27. Reducing modulo3^a before that division would destroy the address.

The compiler also builds the inherited inverse-ray AST. Independent AST
evaluators recover both requested residues without expanding the source.
At a=100,k=80 the tested parameter has495 bits, while the source has more
than10^140 binary digits. Its odd first-hit rank is still three.

## 4. Why this strengthens address coverage without contradicting the old gate

A useful general identity exposes the lost coordinate. For a valid family

    n_t=(A*2^(L t)-B)/3^d,

where L is positive even,3 does not divide A, and
`nu=v3(2^L-1)>=d`, one has

    v3(n_t-n_s)=nu-d+v3(t-s).                         (5)

Thus a shallower inverse prefix fixes nu-d initial ternary digits. For
L=18, nu=3: one inverse division fixes two digits; two fix one; three
free all ternary digits. This is an exact parameter-to-source map, not a
claim that arbitrary certified prefixes have interchangeable hubs.

The earlier one-row theorem fixes its hub and varies a six-exponent block.
Conditioning its source on a19-adic residue restricts that block modulo3.
Our constructor changes the retained inverse-depth data and uses a finite
bank of three-edge prefixes. It therefore has a different parameter map.
The old compatibility theorem is not retracted or bypassed on its own domain.

The uniform rank-four lower bound is inherited from the
[critical residue proof](mod19_route_lifts_20261004.md), section3: no
completed source in residue6 modulo19 has odd first-hit rank at most three.
An independent complete finite control here reduces each exponent modulo18
to1..18. Such changes preserve source residues modulo19 and all integrality
tests through three inverse divisions, because2^18=1 modulo both19 and27.
Allowing intermediate root self-returns only enlarges this test universe;
residue6 is still absent. Formula(2) attains rank four.

## 5. Remaining height and exact-source boundaries

The nineteen disjoint rays have counting function

    #(bank intersect[1,X])=(19/18)*log2(X)+O(1).

They are jointly dense in every declared3-adic/19-adic address while having
natural density zero as positive integers. A certificate at a requested
address belongs to the constructed n_t, not to every integer with that
address. The script explicitly retains a different odd integer with the
same two residues as a source-identity hostile.

The first-hit odd rank is bounded; ordinary time and height are not. With
the canonical parameter in(4), the displayed bank gives ordinary rank at
most `18*3^a*19^(k-1)+8`, since its largest base rank is26. This is an
address-construction bound, not a bound for a given integer's Collatz time.

Source -> target: a checked three-edge inverse prefix with rooted hub ->
joint source residues. Preserved: exact first-hit certificate and both
address coordinates. Lost under residues: integer magnitude and the lift
quotient of t. Required sidecar: the actual parameter, prefix, hub and AST.
The cheapest decisive tests are the critical residue6, all three source
rows in that ray, and an independently evaluated source alias.

## 6. Exact reproduction

    python -X utf8 -B 04-computation/experiments/joint_home_address_20261004.py
    python -O -X utf8 -B 04-computation/experiments/joint_home_address_20261004.py

The script constructs the full stated witness box, checks228 actual first-hit
routes, and verifies5700 pairs of parameter differences at both primes.
For every a=0,1,2,3 and k=1,2 it independently enumerates the complete
parameter language and checks all15200 joint addresses against the digit
compiler, including their next progression lift. The AST readers supply
an independent path to each requested residue. All words through odd rank
three with exponents1..18 are tested for the lower-rank obstruction.
Boolean, negative, out-of-range and zero-prime-precision inputs are rejected.
All tests use exact integers and remain active under optimization.
