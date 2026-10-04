# Primitive-prime layers measure the shared ternary depth of inverse-ray addresses

2026-10-04. **PROVED** elementary order, compatibility and all-depth existence
theorems below. **FINITE-EXACT** for the declared prime and address censuses.
No historical novelty claim. These results constrain address selection
inside a fixed guarded inverse family; they do not certify an arbitrary
supplied integer or prove Collatz convergence.

Artifacts: [script](../../04-computation/experiments/mod19_resonance_depth_20261004.py)
and [output](mod19_resonance_depth_20261004.out).

## 1. Inheritance and the connection that survives a failed prime match

The [inverse-ray codec](inverse_ray_ternary_addresses_20261004.md) proves
that one fixed source row has parameter b and ternary address period
`3^(a-1)` at source precision `3^a`. The
[mod19 selector, sections 4–5](mod19_recursive_observers_20261004.md) adds
the period `3*19^(k-v-1)` after peeling the hub valuation v, and proves that
the two addresses must agree on a shared block digit modulo 3.

The [cubic depth tower, section 2](cyclotomic_depth_towers_20261004.md)
already proves an all-height primitive-prime factorization of its missing
orbit index. Its exact small factors are 19, 87211 and then
`163*135433*272010961`; each factor carries a specified base-2 order.
Section 3 explicitly rejects equality with the golden norm tower at the
next layer: 87211 differs from 5779, and base 2 has order 5778 modulo 5779,
whereas the two golden roots have orders 27 and 54.

Closest mechanism: exact element order followed by noncoprime CRT.
Canonical hostile: equal shared ternary depth need not mean equal full
clocks. Corrected near miss: `9 | p-1` is necessary, not sufficient, for
the base-2 clock to supply a shared ternary digit. Least-used sidecar:
the actual block generator 64, with an anchored hub and source row.
The board is **primitive prime / element order / shared digit / hub
valuation / source guard / quotient loss**.

## 2. A prime's shared-depth invariant

For an odd prime `p!=3`, put

    d_p=ord_p(2),
    R(p)=max(0,v_3(d_p)-1).

**PROVED.** At every prime-adic precision `k>=1`,

    v_3(ord_(p^k)(64)) = R(p).                       (1)

First `ord_p(64)=d_p/gcd(d_p,6)`, so its 3-adic valuation is R(p).
Reduction from units modulo `p^k` to units modulo p has a kernel of order
`p^(k-1)`. The order at the higher level is consequently the base order
times a power of p. Since `p!=3`, that changes none of its 3-primary part.
This proof allows flat order lifts; it assumes no non-Wieferich condition.

For a fixed positive odd hub u with `3` not dividing u and a fixed row r,
the guarded inverse channel is

    n_b=(2^(kappa_r(u)+6b)*u-1)/3,       b>=0.

If the hub is certified, every allowed member carries that exact parent
certificate. Write `v=v_p(u)`. At source precision `p^k` its block period is

    N=1                              if k<=v,
    N=ord_(p^(k-v))(64)               if k>v.          (2)

Indeed equality of two addresses is equivalent, after subtracting and
canceling the unit factors and `p^v`, to `64^(b-c)=1 mod p^(k-v)`.
Thus each block residue modulo N gives exactly one reachable address.
The collapsed case in (2) is essential.

The ternary selector supplies a block residue `b3 mod M`, `M=3^(a-1)`;
the auxiliary selector supplies `bp mod N`. They can refer to one source
iff `b3=bp mod gcd(M,N)`. If `k>v`, (1) makes this exactly

    b3=bp mod 3^min(a-1,R(p)).                       (3)

When compatible, one progression modulo `lcm(M,N)` gives every block.
Among all pairs of separately reachable addresses, the compatible fraction
is precisely `3^(-min(a-1,R(p)))`. If `k<=v`, the index is 1. For `a=1`
the ternary address retains only the already selected row, so its index
is also 1. These are finite clock frequencies, not densities of positive
source integers.

This is the useful invariant: raising the **same** prime-adic precision
never adds shared ternary digits once the hub valuation has been passed.
Changing the prime can add them. The whole auxiliary address still retains
its other phase factors; (3) describes only its interface to the ternary
block address.

### Why 19 is first, and why the group size is not enough

The block clock has a factor 3 iff `9 | ord_p(2)`, which implies `9 | p-1`.
The first prime congruent to 1 modulo 9 is 19, and `ord_19(2)=18` confirms
that it works. Therefore 19 is the smallest prime distinct from 2 and 3
whose noncollapsed six-exponent block clock shares a ternary digit. No
smaller admissible prime develops this property at a higher power, by (1).

The converse congruence fails already at 127:
`9 | 126`, but `ord_127(2)=7` and R(127)=0. The units contain appropriate
torsion, but the specified generator does not visit it.

## 3. Every positive shared depth occurs in the inherited cubic tower

Define the integer

    C_r=(2^(2*3^(r-1))-2^(3^(r-1))+1)/3,    r>=1.

This is the earlier note's `Phi_(2*3^r)(2)/3`; we retain the divided
quantity in the current symbol. **PROVED all-depth supply:** for every
`r>=2`, C_r has a prime divisor p, and every such prime satisfies

    ord_p(2)=2*3^r,            R(p)=r-1.              (4)

Here is the inherited elementary proof, to expose the exact dependency.
Set `t=2^(3^(r-1))`. Since `t=-1 mod3`, writing `t=-1+3h` gives
`t^2-t+1=3(1-3h+3h^2)`. Its 3-adic valuation is exactly one. Therefore
`3` does not divide C_r; C_1=1, while C_r>1 for r>=2. Its prime divisors
are odd and differ from 3. For any such p, `t^2-t+1=0 mod p`, so
`t^3=-1`. The option `t=-1` would force p=3. Thus t has exact order 6.
If d is the order of 2, then `d/gcd(d,3^(r-1))=6`, which forces
`d=2*3^r`. This proves (4). Different levels have disjoint prime supports
because an element's order is unique.

No assertion that C_r is prime is needed, and no global primitive-divisor
theorem is imported. The exceptional C_1=1 records the absence of a fresh
order-six prime. For every desired shared depth D>=1, any prime divisor
of C_(D+1) supplies exactly that depth in (3).

| r | C_r | Shared ternary depth at every prime factor |
|---:|---:|---:|
|1|1|no auxiliary prime|
|2|19|1|
|3|87211|2|
|4|163 * 135433 * 272010961|3|

This is a source-aware use for the old missing-phase factorization: its
successive primitive-prime layers carry one, two, three, ... shared block
digits. It is a theorem about the clocks of guarded inverse families,
not a map from finite-field elements to positive Collatz sources.

## 4. A weaker invariant reconnects the unequal binary and golden layers

The earlier exact order certificates give

    ord_87211(2)=54,       ord_87211(64)=9,
    ord_5779(2)=5778,      ord_5779(64)=963=9*107.

Both primes therefore have R=2. In an anchored inverse ray, both impose
agreement on exactly two ternary block digits once `a>=3`. They agree on
the exponent quotient `b mod9`; they do not identify the full block clock.
The second clock retains another factor 107. The decisive old hostile
`2^54=2944 mod5779`, rather than 1, remains valid.

This is one verified common quotient after the stronger prime/order
identification failed. No assertion about R at all later golden norm
factors follows. The map's source is a prime observer with its specified
block generator and origin; the target is the cyclic ternary quotient
`Z/(3^R)`. It preserves shared block-phase compatibility, loses the other
order factors and all integer-route information, and requires the original
hub, row and certificate to recover a certified family.

There is also a flat-lift hostile to an unnecessary strengthening of (1).
Exact computation at p=3511 gives block orders

    585, 585, 2053935

at precisions p, p^2, p^3. The first order does not grow, but R=2 at
all three levels. Our proof needs only that possible order growth is a
power of p, so this case is included without an exception to the theorem.

## 5. Reproduction and exact universes

```text
python -X utf8 -B 04-computation/experiments/mod19_resonance_depth_20261004.py
python -O -X utf8 -B 04-computation/experiments/mod19_resonance_depth_20261004.py
```

The script uses exact trial division for primality, exact integer
factorization, and modular-return certificates excluding every prime
shortening of a claimed order. There are no probable-prime assumptions.

It checks all 166 primes from 5 through 1000, plus 3511, 5779 and 87211;
orders through p^3; ternary precisions a=1,...,6 and hub valuations v=0,...,3,
for 12,168 lifted gcd controls. The first four inherited C_r factorizations
and all their primitive orders are independently recertified. In the
prime census the R histogram is `{0:143,1:16,2:5,3:1,4:1}`; the least
primes in that complete range with exact positive depths 1,2,3,4 are
19,271,163,487. These need not increase with depth.

An independent address test uses the actual completed root ray
`n_b=(2^(4+6b)-1)/3`. For each p in
`{5,7,19,73,109,127,163,5779,87211}` and a=1,...,4 it enumerates the entire
joint block period and compares its actual address pairs with every pair
accepted by (3). All 36 full languages agree, including the no-resonance,
one-digit, two-digit and three-digit cases. They test actual source
addresses without using an orbit-search oracle or extrapolating coverage
of arbitrary integers.
