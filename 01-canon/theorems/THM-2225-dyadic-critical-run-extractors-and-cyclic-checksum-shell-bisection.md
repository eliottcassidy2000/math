---
id: THM-2225
title: "Dyadic critical-run extractors and cyclic-checksum shell bisection"
status: >
  PROVED + VERIFIED-EXACT. An unknown-bias Bernoulli coin admits two
  pathwise-total fair-bit extractors on every nonconstant input stream.
  If n is the initial constant-run length, dyadic block compression stops
  by 2n flips. A cyclic checksum on the second half of the first
  nonconstant dyadic prefix stops by max(2,2n-1). Fairness is exact for
  every 0<p<1: the first construction has a weight-preserving block-swap
  involution, while the second bisects every Hamming shell by a cyclic
  response whose image contains the antipode. The dyadic hypothesis is
  load-bearing for that checksum mechanism.
source: codex-gauss-2026-07-24-dyadic-coin-extractors
related:
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/dyadic_coin_critical_run_checksum_thm2225.py
output: 05-knowledge/results/dyadic_coin_critical_run_checksum_thm2225.out
script_sha256: 7b094c491b2e6fccdf740b57a78979c44a00c6a883894ff4629cd3e4074249c8
output_sha256: c1a66ce97c49f7b6566d2ee9e83a7b2c52cdb47a78d8d5d4e2762ac29a54aa15
hash_basis: working-tree bytes (LF)
---

# THM-2225 -- dyadic critical-run fair-coin extractors

Let the input bits be independent, with

```text
P(0)=p,       P(1)=q=1-p,       0<p<1.              (1)
```

For a nonconstant stream, let `n` be the length of its initial constant
run. Define

```text
M=2^ceil(log_2(n+1)),       m=M/2.                   (2)
```

Equivalently, `M` is the unique power of two such that

```text
m<=n<=M-1.                                           (3)
```

The first `m` input bits are constant and the first `M` are not. Thus the
sets

```text
S_M={w in {0,1}^M:
     w_1=...=w_m and w is nonconstant}               (4)
```

partition the nonconstant streams according to their first applicable
dyadic prefix. Constant infinite streams have probability zero under (1).

## 1. The `2n` block-compression extractor

Read through bit `M`. Starting with the word `w` of length `M`, perform
the following deterministic compression.

1. Pair adjacent symbols.
2. If some pair is unequal, choose the leftmost unequal pair and declare
   `01` heads and `10` tails.
3. If every pair is equal, replace `00` by `0` and `11` by `1`, then
   repeat on the compressed word.

A nonconstant power-of-two word cannot compress all the way to one symbol
without exposing an unequal pair, so the rule is pathwise total on (4).

Suppose the decisive pair occurs after `r` compressions. It represents two
adjacent homogeneous dyadic blocks of length `2^r`. Swap those two blocks.
This operation is a fixed-point-free involution on `S_M`:

- all lower-level pairs remain equal;
- all same-level pairs to the left remain unchanged;
- the decisive pair remains the leftmost unequal pair, with its order
  reversed;
- the numbers of zeros and ones are preserved.

It also preserves membership in `S_M`. Below the top compression level,
aligned block pairs cannot cross the midpoint `m`; at the top level the
operation merely swaps two constant halves. Hence the involution pairs
every heads word with a tails word of exactly the same Bernoulli weight

```text
p^(number of zeros) q^(number of ones).              (5)
```

It follows that the output is fair for every `p`. Finally, (3) gives

```text
M=2m<=2n.                                            (6)
```

This proves part (a).

## 2. The cyclic-checksum extractor

The sharper construction uses a different carrier. If `M=2`, declare
`01` heads and `10` tails. Now assume `M>=4`, so `m>=2`.

Write the second half of the `M`-bit word as

```text
y_i=w_(m+i),       1<=i<=m,
```

and define the cyclic checksum

```text
s(y)=sum_(i=1)^m i y_i mod m,
       represented in {0,...,m-1}.                  (7)
```

Declare

```text
heads iff s(y)<m/2,
tails iff s(y)>=m/2.                                 (8)
```

The last coefficient in (7) is zero modulo `m`. Consequently, if the
initial transition has occurred by flip `M-1`, the decision is already
known after `M-1` flips; the unobserved bit `y_m` cannot change it. Only
when `n=M-1` is flip `M` needed to reveal the transition and select this
dyadic shell.

## 3. Exact Hamming-shell bisection

It is enough to prove that (8) splits every total Hamming-weight class in
`S_M` equally. All words in one such class have the same probability (5).

First fix an interior tail weight

```text
1<=j=sum_i y_i<m.                                    (9)
```

Let `rho` cyclically move the bit in position `i` to position `i+1`,
with `m` moving to `1`. Then

```text
s(rho y)=s(y)+j mod m.                              (10)
```

Put `g=gcd(j,m)`. Along any cyclic-rotation orbit, the checksum visits
every residue in one coset modulo `g` equally often. Indeed, the checksum
advances by `j`, whose additive order modulo `m` is `m/g`; if the word
has a longer stabilizer period, each residue is merely repeated the same
number of times.

Because `m` is a power of two and `j<m`,

```text
g divides m/2.                                      (11)
```

Translation by `m/2` therefore preserves each checksum coset and pairs
its residues below `m/2` with those at least `m/2`. Every rotation orbit,
and hence every tail-weight class (9), has equally many heads and tails.

For total word weights other than `m`, the constant first half determines
whether the tail is appended to `0^m` or `1^m`, so the preceding argument
applies directly. At total word weight `m`, there are only two words:

```text
0^m 1^m,       1^m 0^m.                             (12)
```

Their tail checksums are

```text
m(m+1)/2 = m/2 mod m,       0,                      (13)
```

respectively, so (8) separates them. This proves shell-by-shell equality
of the two output probabilities.

Stopping at `M-1` does not disturb the proof. Expanding an early terminal
prefix by its two possible last bits gives two members of `S_M` with the
same checksum and output. Thus the variable-length prefix code is exactly
the full-shell coloring above after the identity `p+q=1`.

## 4. The sharp pathwise bound

If `n<=M-2`, the rule stops at `M-1`, and (3) gives

```text
M-1=2m-1<=2n-1.                                     (14)
```

If `n=M-1` and `M>=4`, it stops at `M`, while

```text
M<=2M-3=2n-1.                                       (15)
```

For `n=1`, it stops at two flips. Therefore every nonconstant input
satisfies

```text
tau<=max(2,2n-1).                                   (16)
```

This proves part (b).

## 5. Response images, tournaments, and the LRC warning

On the weight-`j` tail shell, (10) is the response homomorphism

```text
Z/mZ -> Z/mZ,       delta |-> j delta.              (17)
```

An output-reversing rotation exists exactly when the antipode `m/2` lies
in its image:

```text
j delta=m/2 mod m
iff gcd(j,m) divides m/2.                           (18)
```

For every `1<=j<m`, condition (18) holds exactly when even `m` is a power
of two. If `m` has an odd factor, write `m=2^a u` with odd `u>1` and take
`j=2^a`; then `gcd(j,m)` does not divide `m/2`. The smallest checksum
hostile example is `m=6,j=2`: the lower/upper counts are `7` and `8`.

This distinguishes two notions that look similar in a Gram shadow.
THM-2195's cyclic Hamming derivative and subgroup variance detect that a
nonconstant word moves under rotation. Fair bisection needs the stronger
Smith/image statement that the required antipode is actually attained.

The distinction is also load-bearing in THM-2222. Its identity `LR=-R`
is a transfer-operator eigenline, but `L` is thirteen-to-one rather than
an invertible orthogonal group action; after normalization the eigenvalue
is `-1/13`, not a sign character. Therefore no orbit-balance conclusion
follows. The cover domination and even/odd checkpoint inclusions are the
necessary sidecar. The coin construction succeeds because (17) retains
the exact response image, not merely a positive variance or eigenvalue.

## 6. Exact referee

Run

```bash
python3 04-computation/dyadic_coin_critical_run_checksum_thm2225.py
python3 -O 04-computation/dyadic_coin_critical_run_checksum_thm2225.py
```

The companion exhausts `S_M` through `M=16`, checking both decisions,
the block-swap involution, every Hamming-shell count, last-bit merging,
and every pathwise bound. It separately audits (11), (13), and (16) for
all dyadic `m<=4096`, and freezes the `m=6,j=2` hostile example. The
ordinary and optimized outputs are byte-identical. QED.
