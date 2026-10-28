# All-order low-bit checksum for positive tournament channels

**Status: PROVED from the inherited Burnside channel, with FINITE-EXACT checks.**
This is an all-order congruence, not a new counting algorithm or a certificate
of every bit of a computed count. No literature-novelty claim is made.

For integers r>=1, s>=0, write

```
F_n(r,s) = sum_(odd lambda partitions n) 2^(r E(lambda)+s k)/z_lambda,
k = length(lambda),
E(lambda) = (n-k)/2 + sum_(i<j) gcd(lambda_i,lambda_j),
z_lambda = product_l l^(m_l) m_l!.
```

Let v2 be the exponent of 2 in a nonzero rational number and phi Euler's
totient. The inherited carrier interpretation ensures F_n(r,s) is a positive
integer. Then:

1. **Odd n>=1:** `v2(F_n(r,s)) = r(n-1)/2+s`.
2. **Even n>=4:** with `B=rn/2+2s`,
   `F_n(r,s)/2^B = phi(n)/2 (mod 2)`.
   In particular `v2(F_n(r,s))>=B`, with equality exactly when `phi(n)/2` is odd.
3. **Boundary cases:** `F_0(r,s)=1` and `F_2(r,s)=2^(r+2s-1)`.

The congruence in item 2 includes its assertion that the quotient is an
integer. When phi(n)/2 is even this theorem gives a strict lower bound, **not**
an exact higher valuation.

## Proof

All cycle lengths are odd, so

```
v2(z_lambda) = sum_l v2(m_l!) <= v2(k!).
```

The inequality follows because `k!/product_l m_l!` is an integer. Also

```
E(lambda) >= (n-k)/2 + binomial(k,2).
```

Work first in the localization `Z_(2)` of rationals with odd denominators,
where every odd denominator is 1 modulo 2. After establishing a valuation
bound for the rational Burnside sum, its already-proved integrality turns
that bound into ordinary integer divisibility.

### Odd n

The number k of odd parts is odd. The k=1 partition (n) contributes

```
2^A/n,   A=r(n-1)/2+s,
```

whose valuation is A and whose normalized residue is 1 modulo 2.

For odd k>=3, `v2(k!)<=k-2`: for example Legendre's identity
`v2(k!)=k-popcount(k)` applies, and an odd integer k>=3 has at least two 1s
in its binary expansion. Therefore every other term has valuation minus A
at least

```
r (k-1)^2/2 + s(k-1) - (k-2)
>= (k-1)^2/2 - (k-2)
= ((k-2)^2+1)/2 >= 1.
```

After dividing the sum by 2^A, only the single n-cycle contributes to the
low bit. The normalized sum is odd, proving item 1.

### Even n>=4

Now k is even and at least 2. For k>=4 use `v2(k!)<=k-1`. Every such term
has valuation minus `B=rn/2+2s` at least

```
r k(k-2)/2 + s(k-2) - (k-1)
>= k(k-2)/2-(k-1)
= ((k-2)^2-2)/2 >= 1.
```

**Hostile factorial bound:** the odd-k inequality `v2(k!)<=k-2` cannot be
used indiscriminately here: `v2(4!)=3>2`. The weaker valid bound above suffices.

It remains to examine two odd parts a,b with a+b=n.

- If a and b are distinct, `z=ab` is odd and
  `E=n/2-1+gcd(a,b)`. Coprime pairs have valuation exactly B, and their
  normalized term is `1/(ab)=1 (mod 2)`. If they are not coprime, their gcd
  is an odd integer at least 3, so their valuation exceeds B by at least 2r.
- If a=b, then a=n/2 is odd. Since n>=4, necessarily a>=3. Here
  `E=2a-1`, `z=2a^2`, and the valuation surplus over B is
  `r(a-1)-1>=1`.

Thus the only lowest-valuation partitions are the unordered distinct
coprime odd pairs `{a,n-a}`. For even n, all a coprime to n are odd and
`gcd(a,n-a)=gcd(a,n)`. There are phi(n) such a in 1,...,n-1. The involution
`a -> n-a` has no fixed point among them when n>=4, so the number of
unordered pairs is exactly phi(n)/2. Their residues add modulo 2, proving
item 2. The empty partition and the unique odd partition (1,1) of 2 give
item 3 directly. This also explains why n=2 must be separated.

## Consequence for the paired-tournament output

The aligned independent pair channel is `P_n=F_n(2,1)`. Consequently

```
odd n:       P_n = 2^n                    (mod 2^(n+1));
even n>=4:  P_n = (phi(n)/2 mod 2) 2^(n+2) (mod 2^(n+3));
P_0=1, P_2=8.
```

The check depends on n and a totient, not a full partition traversal. It can
reject an incorrect result's low bits cheaply, but cannot certify its higher
bits, distinguish every wrong answer, or justify a performance speedup.

## Bounded provenance check

The proof uses the same positive Burnside channel as Session 22 and does not
alter Session 89's frozen implementation. A targeted content search at
commit `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf` inspected the following known
counting sources for totient, valuation, 2-adic, low-bit, and factorial-two
statements:

- `05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/`;
- `05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/`;
- `05-knowledge/results/session87_four_vertices_enumeration_20261028/FAST_COUNTS.md`;
- `01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`;
- `01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md`;
- `04-computation/a000568_turbo.py` and `04-computation/a000568_gmp_enum_v4.c`.

Those sources contain the inherited factorial-two cancellation but the exact
all-order valuation/totient statement was not found in this bounded search.
That is **not** an exhaustive repository or literature search and establishes
no novelty priority.

## Verification

Run from this directory:

```sh
python3 check_low_bits.py
python3 -O check_low_bits.py
```

Both runs pass, without reading or requiring a committed raw JSON receipt.
The checker:

- computes six full channels `(1,0),(1,3),(2,0),(2,1),(3,2),(4,1)` through n=80
  using the unchanged exact fused counter;
- checks every asserted valuation/congruence and the two boundary cases;
- independently compares rational odd-partition sums through n=26;
- checks 6,402 individual summand valuation classifications;
- independently counts the coprime odd pairs and tests the invalid k=4
  factorial-bound shortcut.

Across the six channels it records 240 odd exact-valuation cases, 54 even
lower-bound equalities, 180 strict even cases, and 12 boundary cases. These
finite checks support the implementation; the proof above supplies the
all-order statement. No benchmark range was expanded.
