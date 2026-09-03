---
id: THM-4371
title: "LRC(14) minimal wall-LCM adaptive owner clock"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 AND THM-2061/2066/4070/4349 +
  VERIFIED-EXACT. For an eleven-element positive body H of height h, let
  D_H be the lcm of its doubled speeds and all positive odd a<12h. Universal
  odd two-tail safety is equivalent to emptiness of the single owner relation
  at clock 7D_H. Equivalently D_H is the odd-tail lcm times
  2^(1+max_{c in H} v_2(c)); at fixed height the body dependence is only this maximal
  2-adic scale. The height-only clocks 7D_h and
  7 lcm(1,...,12h-1) are also complete. This sharpens only THM-4349's
  explicit clock and is not LRC(14).
source: codex-2026-09-03-LRC-minimal-wall-lcm-clock
depends_on:
  - LRCUpTo13
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-4070-lrc14-d2-mod14-two-bank-affine-ray-firewall
  - THM-4349-lrc14-adaptive-owner-clock-completeness-and-twelve-height-tail-cutoff
related:
  - THM-4075-lrc14-divisor-complete-dyadic-owner-word-closure-through-30
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization
  - THM-4365-lrc14-cofinite-828-quotient-fibre-and-centered-residue-exit-law
  - THM-4367-lrc14-active-first-exit-scale-collision-classification
script: 04-computation/lrc14_adaptive_clock_lcm_sharpening_probe_adaptive_clock_20260903.py
output: 05-knowledge/results/lrc14_adaptive_clock_lcm_sharpening_probe_adaptive_clock_20260903.out
script_sha256: 6eb8b09a135aa98b241562f2480d8dad9ef282a77abf536674e27ec07d3d47a9
output_sha256: b591050f7ebd533de8e12f69a0721b08ea62cea6323432a252b6eb4d6f674efe
hash_basis: raw LF bytes
---

# THM-4371 -- one minimal wall LCM is a complete adaptive owner clock

**PROVED RELATIVE TO CITED `LRCUpTo13` AND THM-2061/2066/4070/4349 +
VERIFIED-EXACT.**  Let `H` be an eleven-element set of distinct positive
integers and put

```text
h=max H,
O_h={a:1<=a<12h and a is odd},
D_H=lcm({2c:c in H} union O_h),
Q_H=7D_H.                                                (1)
```

Then the following are equivalent.

1. `R_(Q_H)(H)=empty`, for the labelled complementary-owner relation of
   THM-2066.
2. There is some positive clock `N` with `R_N(H)=empty`.
3. `2H union {a,b}` is `1/14`-lonely for every two distinct positive odd
   integers `a,b`.
4. The same assertion holds just for distinct positive odd
   `a,b<12h`.
5. No distinct positive odd pair has THM-4070's continuum containment
   `G(H) subset Sigma_2({a,b})`.

Thus the explicit clock in THM-4349,

```text
28(42h+1)^2 lcm(1,...,12h),                             (2)
```

may be replaced by `Q_H`.

There are two body-independent versions.  Define

```text
D_h=lcm({2,4,...,2h} union O_h),
Q_h=7D_h,
Q_tilde_h=7 lcm(1,...,12h-1).                          (3)
```

The five conditions are also equivalent to `R_(Q_h)(H)=empty`, and to
`R_(Q_tilde_h)(H)=empty`.

These are completeness clocks, not assertions that their relations are
always empty.  The theorem does not close the dyadic seam or LRC(14).

## 1. Exact normal form and inherited safe block

Put

```text
L_h^odd=lcm(O_h),                nu_H=max_{c in H} v_2(c). (4)
```

The odd part of every `c in H` is itself an element of `O_h`, and therefore
divides `L_h^odd`.  The only new prime-power contribution of `2c` is its
`2`-part.  Hence

```text
D_H=2^(1+nu_H)L_h^odd,
D_h=2^(1+floor(log_2 h))L_h^odd.                        (5)
```

In particular `D_H|D_h|lcm(1,...,12h-1)`.  Formula `(5)` also shows that, at
fixed height, the body dependence in `(1)` is exactly its largest `2`-adic
core scale.

The equivalence of conditions 2--5 is THM-4349, using THM-2061/4070 for the
finite metric box and continuum criterion.  We prove the new implication
`4 => 1`; THM-2066 gives `1 => 3`.

Cited `LRCUpTo13` supplies a point `y_0` satisfying

```text
||cy_0||>=1/12                         for every c in H.
```

Since `y -> ||cy||` is `c`-Lipschitz and
`1/12-1/14=1/84`, the closed arc

```text
I={y:dist(y,y_0)<=1/(84h)},             |I|=1/(42h),    (6)
```

lies in

```text
G(H)={y:||cy||>=1/14 for every c in H}.                 (7)
```

At any clock `Q`, the `Q`-grid points in `I` contain a cyclic block of `K`
consecutive residues with

```text
K>=Q/(42h)-1.                                           (8)
```

Every residue in this block belongs to THM-2066's labelled safe packet
`A_Q(H)`.  Closed core safety and strict tail danger are retained at the two
endpoints of `I`.

## 2. Doubling cancels one factor two in every rational wall

Assume condition 4.  Fix distinct positive odd `a,b<12h`.  The safe-time set
of the physical row

```text
S(H;a,b)=2H union {a,b}                                (9)
```

is a nonempty closed finite union of intervals and points.  It is not the
whole circle, since time zero is unsafe.  Some component endpoint is
therefore a safety wall for a physical speed `m in S(H;a,b)`:

```text
x=(14k+/-1)/(14m).                                    (10)
```

In the owner packet the body phase is `y=2x`, so

```text
y=(14k+/-1)/(7m).                                     (11)
```

Its reduced denominator divides `7m`, not merely `14m`.  There are exactly
two possible speed species:

```text
m=2c with c in H,                 or                m in O_h.
```

It follows from `(1)` that every bounded pair has a witness

```text
y_(a,b)=r_(a,b)/Q_H in G(H) minus Sigma_2({a,b}).       (12)
```

Indeed

```text
Q_H=lcm({14c:c in H} union {7a:a in O_h})=7D_H.        (13)
```

The selected endpoint `x` is one of the two lifts over `y` and is safe for
both tails.  Thus the witness in `(12)` is literal, including its lift label.
If a component is a singleton, its endpoint is still valid: safety uses
weak inequalities, while danger and `Sigma_2` use strict inequalities.

## 3. This wall LCM already exceeds the long-block threshold

Let `Q=Q_H`, and let `z mod 2Q` be any eligible odd tail class.  On the block
in `(8)`, subtracting the strict `Q/7` eligibility inequality at the first
grid point from the next `K-1` inequalities gives

```text
|jz|_Q<2Q/7,                         1<=j<=K-1.         (14)
```

For completeness, THM-4349's no-wrap lemma says that if `q` is the centered
representative of `z modulo Q`, then `(14)` implies

```text
|q|<2Q/[7(K-1)]
   <=2Q/[7(Q/(42h)-2)].                                (15)
```

The final expression is below `12h+1` whenever

```text
Q>1008h^2+84h.                                         (16)
```

The clock `(1)` satisfies this strictly without auxiliary padding.  The two
largest elements of `O_h` are coprime, and `D_H` is even:

```text
gcd(12h-1,12h-3)=1,
D_H is divisible by 2(12h-1)(12h-3).
```

Consequently

```text
Q_H
 >=14(12h-1)(12h-3)
 =(1008h^2+84h)+42(24h^2-18h+1)
 >1008h^2+84h.                                         (17)
```

The final inequality holds for every integer `h>=1`; here necessarily
`h>=11`.  Also `Q_H` is even, so centering an odd residue modulo `Q_H`
preserves odd parity.  Every eligible class therefore has a nonzero odd
centered representative satisfying

```text
|q|<12h+1,                         hence |q|<=12h-1.    (18)
```

Moreover `(8)` and `(16)` give `K>24h+1`, so the core-safe block contains two
adjacent residues.  This will recover the one owner bit apparently lost by
centering modulo `Q`.

## 4. Adjacent owner rigidity forces the top bits to agree

Suppose for a contradiction that

```text
(u,v) in R_Q(H),                         Q=Q_H.          (19)
```

Both classes are individually eligible.  Let `a,b` be their centered
representatives modulo `Q`.  By `(18)`, these are nonzero odd integers with

```text
|a|,|b|<12h.                                            (20)
```

The original tails live modulo `2Q`, so there are unique bits
`epsilon,delta in {0,1}` with

```text
u=a+epsilon Q mod 2Q,
v=b+delta Q mod 2Q.                                    (21)
```

We first prove the local rigidity needed to compare these bits.  Extend a
packet label cyclically to an integer `s`.  Replacing a canonical label by
`s+Q` adds an odd integer to each odd tail's nearest integer, so it flips
both owner bits and preserves complementarity.  Choose consecutive lifted
labels `s,s+1` from the block `(8)`.  For `q=a` or `q=b`, put

```text
n_q(t)=nint(qt/Q),
e_q(t)=qt/Q-n_q(t).                                    (22)
```

Eligibility gives

```text
|e_q(s)|<1/7,                         |e_q(s+1)|<1/7.
```

Since `q` is centered, `|q|/Q<=1/2`, and hence

```text
n_q(s+1)-n_q(s)=q/Q-e_q(s+1)+e_q(s),
|n_q(s+1)-n_q(s)|<1/2+2/7=11/14<1.                    (23)
```

The left side is an integer, so it is zero.  Thus each centered owner is
constant across these two adjacent labels.

On the other hand, `(21)` gives the exact owner identities

```text
omega_u(t)=n_a(t)+epsilon t mod 2,
omega_v(t)=n_b(t)+delta t mod 2.                        (24)
```

The two owner words are complementary at both `s` and `s+1`.  Subtracting
the two complementarity equations modulo two and using `(23)` yields

```text
epsilon+delta=0 mod 2,                  so epsilon=delta. (25)
```

Thus the apparently missing top bit is common to the two tails and cancels
at every packet label.  Eligibility depends only on reduction modulo `Q`, so

```text
(a,b) mod 2Q belongs to R_Q(H).                         (26)
```

If `|a|!=|b|`, apply the synchronized witness `(12)` to the distinct positive
tails `|a|,|b|`.  Changing a tail's sign preserves its literal two-lift
danger mask, so `(a,b)` fails to cover that witness, contradicting `(26)`.
If `|a|=|b|`, the two signed tails have identical masks and cannot cover both
lifts over any nonempty packet point.  The block `(8)` is nonempty, giving
the same contradiction.  Repeated residue classes have therefore been
handled without assuming that actual tails or clock residues are distinct.

We conclude that `R_(Q_H)(H)=empty`, proving `4 => 1`.

## 5. Height clocks and divisor transport

The divisibilities following `(5)` give

```text
Q_H | Q_h | Q_tilde_h.                                 (27)
```

THM-2066 divisor transport maps a relation at a multiple clock into the
relation at every divisor clock.  Hence emptiness at `Q_H` propagates upward
to both clocks in `(3)`.  Conversely, emptiness at either larger clock
certifies condition 3 by THM-2066.  This proves all claimed equivalences.

The especially simple universal clock is therefore

```text
boxed: 7 lcm(1,...,12h-1),                             (28)
```

while `(1)` and `(5)` retain the strongest common-wall bound proved here.

## 6. Hostiles, audit, and exact scope

The top-bit issue in Section 4 is genuine at one point.  For the height-`11`
clock `Q=7D_11`, take `c=8` and `r=Q/(14c)`.  Then `r` is odd.  At phase
`r/Q=1/112`, the classes `1` and `1+Q` have the same centered residue but
opposite literal masks:

```text
at r:       (danger,safe),              (safe,danger). (29)
```

At the adjacent phase their masks agree:

```text
at r+1:     (danger,safe),              (danger,safe). (30)
```

Thus one-point centering really loses the bit, while the adjacent safe block
recovers it exactly.  This is a literal tail-mask control; the displayed
`r` is not asserted to belong to `A_Q(H)` for every body.  It is the same
controlled-forgetting warning seen in THM-4363/4365/4367, but those theorems
are analogies rather than dependencies.

The max-`72` body

```text
{1,5,11,13,17,19,23,37,41,70,72}
```

is still the canonical practical hostile: all clocks `15,...,43` fail, while
clock `47` is empty.  Its new body clock has `378` decimal digits.  It is
`1,171,280,000` times smaller than `(2)`, but remains astronomically larger
than the exhibited clock `47`; this theorem preserves, rather than erases,
the need for adaptive small-clock search.

The exact companion verifies:

- the wall divisibilities, normal forms, parity, and strict threshold for
  every `1<=h<=256`;
- all `31,824` eleven-bodies of height `11<=h<=18`;
- `40,600` exact doubled-wall fractions;
- `2,566,368` finite cases of the no-wrap lemma;
- `880,594` finite cases of adjacent-owner rigidity; and
- the literal top-bit hostile `(29)--(30)`.

Normal, optimized, and fixed-hash-seed replays byte-match the frozen output.
The hashes in the frontmatter use raw LF bytes.

This theorem sharpens only the explicit completeness clock in THM-4349.  It
does not prove `R_(Q_H)(H)` empty for every body, improve the inherited
`12h` tail cutoff, classify unsafe bodies, or lower the LRC(14) frontier.
The phrase "minimal wall LCM" refers to the exact common multiple `(13)`
furnished by the rational-endpoint proof; it does not assert global numerical
optimality among all possible proof clocks.  **QED.**
