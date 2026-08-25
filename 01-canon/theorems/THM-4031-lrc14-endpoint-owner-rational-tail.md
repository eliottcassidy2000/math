---
id: THM-4031
title: "LRC(14) endpoint-owner rational tail on odd-dilation rays"
status: >
  PROVED ARITHMETIC + EXACT FINITE-PHASE TAIL + VERIFIED-EXACT +
  INDEPENDENTLY REFEREED. On every
  THM-4025 odd-dilation ray, the owner-relaxed margin tends to
  (4U-3t)/(126U). Hence every old-strip ray has an eventual all-survivor
  arithmetic tail, while every ray on the other side closes at every odd
  scale. A 126t-periodic residue numerator gives exact phase thresholds and
  a closed universal bound. This is a stopping theorem for one necessary
  gate, not an LRC counterexample or an all-height LRC conclusion.
source: root + incoming_signal_scan + endpoint_tail_referee / OEIS finite-head rational-tail continuation, 2026-08-24
audit: >
  PASS. Fraction evaluation of the full owner minimum dominates an independent
  integer endpoint numerator in 280 small cells. Every one of the 63t odd
  phases is checked for six hostile/control rays, including periodicity,
  first success and last failure. The canonical ray's complete exact head is
  audited through its proved tail. An independent engine checked 2,640 exact
  gate cells and 12,537 phase residues. Normal and optimized streams byte-match.
depends_on:
  - THM-4025-lrc14-owner-residue-odd-dilation-semigroup
  - THM-4003-lrc14-scale-two-component-erosion-boundary-strip
related:
  - THM-4029-lrc14-ap-cover-twelve-owner-rational-tail
script: 04-computation/lrc14_owner_endpoint_rational_tail_thm4031.py
output: 05-knowledge/results/lrc14_owner_endpoint_rational_tail_thm4031.out
script_sha256: 974e384184ac36e74437c334b0050005cce27613a7ecc7b275739b33667d85b7
output_sha256: c103136dc4c217a81eb48f166de195e0b2aa0ea8c797c48e09c7f389eca8c504
hash_basis: raw LF bytes
---

# THM-4031 -- the endpoint-owner rational tail

**PROVED ARITHMETIC + EXACT FINITE-PHASE TAIL + VERIFIED-EXACT + INDEPENDENTLY
REFEREED.** This theorem
identifies the ordinary-order tail that THM-4025 deliberately left open. It
shows that THM-4003's owner-relaxed residue gate eventually loses all power
on every physical old-strip ray. An arithmetic survivor is only an unresolved
cell, not an LRC counterexample. LRC(14) remains open.

## 1. Inherited gate and theorem

Retain THM-4025's arithmetic objects for odd `t` and `1<U<=t`:

```text
epsilon_i(t,U)=min_(1<=u<=U)e_i(t,u),                 i=1,2,
B(t,U)=max(0,2/63-epsilon_1(t,U)-epsilon_2(t,U)),
D(t,U)=t(2U-1)/(84U(U-1)).                            (1)
```

For every positive odd `k`, put

```text
M_k=(B-D)(kt,kU).                                     (2)
```

Then

```text
lim_(k->infinity,k odd) M_k
 =2/63-t/(42U)
 =(4U-3t)/(126U).                                    (3)
```

The convergence is from below.

Consequently:

- if `4U>3t`, every sufficiently large odd scale is an owner-relaxed
  arithmetic survivor;
- if `4U<=3t`, every positive odd scale is an arithmetic closure.

In particular, every physical old-strip ray of THM-4025 lies on the first
side and has an ordinary-order all-survivor tail.

The closest proved mechanism is THM-4025's divisibility-order monotonicity.
Its canonical hostile is the ordinary-order revival on the ray `(5,4)`.
The least-used sidecar is the target endpoint owner one unit below the scaled
maximum.

## 2. Endpoint owner and quantitative squeeze

At scale `(kt,kU)`, choose

```text
u_k=kU-1.                                             (4)
```

Since `gcd(k,kU-1)=1`,

```text
g_k=gcd(kt,kU-1)=gcd(t,kU-1)<=t.                     (5)
```

The directed numerators are odd and their moduli even, so the positive
residues never vanish. Therefore

```text
0<e_i(kt,u_k)<g_k/u_k<=t/(kU-1),       i=1,2.         (6)
```

The target minima include `u_k`, so `epsilon_i(kt,kU)->0` and
`B(kt,kU)->2/63`. Direct algebra gives

```text
D(kt,kU)=t/(42U)+t/[84U(kU-1)].                      (7)
```

More quantitatively,

```text
(4U-3t)/(126U)-[2t+t/(84U)]/(kU-1)
 < M_k
 <=(4U-3t)/(126U)-t/[84U(kU-1)].                     (8)
```

This proves `(3)` and the survivor side. On the other side,
`t/(42U)>=2/63`, so `(7)` gives `D>2/63>=B` at every scale.

The sign wall is sharp. Equality `4U=3t` cannot occur for odd `t`, because
its two sides have opposite parity. As a purely algebraic gate observation
outside that parity domain, `(7)` would put equality on the closure side;
this is not an extension of the theorem's domain.

## 3. Exact periodic endpoint numerator

The endpoint witness is stronger than the coarse squeeze. Define

```text
r_1(k)=[k(3t-4U)+4]^+_(42g_k),
r_2(k)=[k(9t+16U)-16]^+_(126g_k),
h(k)=3r_1(k)+r_2(k).                                  (9)
```

Then

```text
e_1(kt,u_k)+e_2(kt,u_k)=h(k)/[126(kU-1)],             (10)
```

The minima obey `epsilon_1+epsilon_2<=e_1(kt,u_k)+e_2(kt,u_k)`. Therefore
the truncation in `B=max(0,2/63-epsilon_1-epsilon_2)` is harmless:

```text
B>=2/63-e_1(kt,u_k)-e_2(kt,u_k).
```

Consequently the full margin satisfies the exact sufficient bound

```text
M_k >=
 [2(4U-3t)(kU-1)-2Uh(k)-3t]/[252U(kU-1)].             (11)
```

The tuple `(g_k,r_1(k),r_2(k),h(k))` has valid period

```text
P=126t                                                   (12)
```

on odd `k`; least period is not claimed. Indeed `P` preserves parity and
is divisible by `t`, so `g_(k+P)=g_k`. The change in the first numerator is
`P(3t-4U)`, divisible by `42g_k`; the change in the second is
`P(9t+16U)`, divisible by `126g_k`. Thus both positive residues repeat.

## 4. Exact finite-head/phase-tail compiler

Assume `a=4U-3t>0`. For every odd representative
`s in {1,3,...,P-1}`, set `h_s=h(s)` and

```text
n_s=max(0,
  ceil((2Uh_s+3t-2a(sU-1))/(2aUP))),
K_s=s+n_s P.                                         (13)
```

Along the phase `k=s+nP`, the numerator in `(11)` is affine with positive
slope `2aUP`. Hence its endpoint certificate is nonnegative exactly for
`n>=n_s`. If `n_s>0`, its last failed point is
`L_s=s+(n_s-1)P`. Put

```text
L=max({L_s:n_s>0} union {-1}),
K_phase=1 if L<1, and L+2 otherwise.                  (14)
```

Then every positive odd `k>=K_phase` is certified as an arithmetic survivor.
The infinite tail has compiled into a finite table of `63t` odd phases and
one affine inequality per phase.

There is also a closed all-phase threshold. Since positive residues are
strictly below their even moduli,

```text
h(k)<=3(42g_k-1)+(126g_k-1)
    <=252t-4=:H.                                      (15)
```

Let `OddCeil(A/B)` be the least positive odd integer at least `A/B`. Then

```text
K_univ=OddCeil((2a+2UH+3t)/(2aU))                     (16)
```

certifies every positive odd `k>=K_univ`. This bound is explicit, not
claimed sharp.

## 5. Treating the n-th odd multiplier as n

Under the ordinal chart

```text
k=2n-1,                                               (17)
```

the valid period `126t` in `k` becomes a `63t`-phase clock in the natural
index `n`. The theorem-bearing state is

```text
(n mod 63t, g_k, r_1(k), r_2(k), kU-1).               (18)
```

The first four coordinates are finite periodic data; the last is an
unbounded affine height denominator. This is the OEIS-style finite-head plus
rational-tail representation. A unary closure/survivor bit word loses both
the reason for the tail and its proof.

The offset in `u_k=kU-1` is load-bearing. Choosing the scaled endpoint
`u=kU` instead gives the exact invariant gaps of THM-4025, not `O(1/k)`
decay.

## 6. Exact controls and stopping boundary

The canonical ray `(t,U)=(5,4)` has:

```text
exact endpoint-certificate tail K_phase=1071,
universal bound K_univ=1259,
complete arithmetic closure set={1,3,5,7,11,13}.      (19)
```

The closure set is **FINITE-EXACT**: the audit evaluates every odd scale below
`1071`, after which `(13)--(14)` proves the tail. Thus the ordinary-order
revival in THM-4025 is a finite-head phenomenon, while the endpoint bound is
deliberately not sharp.

Further exact controls include endpoint tails `777` on `(11,9)`, `19` on
`(11,11)`, `2775` on `(13,10)`, and `3627` on `(17,13)`. The hostile
ray `(3,2)` lies on the other side and closes at every odd scale.

This is a genuine stopping theorem: the THM-4003 owner-relaxed gate cannot
close the arithmetic tail of any old-strip ray. It does not say that the
underlying LRC row exists, much less that it is a counterexample. A survivor
means only `D<=B`, so this necessary obstruction failed to close the cell.

THM-4003 gives the LRC failure implication only for physical target heights
`11<=kU<91^6`. The arithmetic limit and tail remain valid beyond that bound,
but receive no new LRC interpretation here. In particular, if
`K_phase U>=91^6`, the asymptotic theorem supplies no in-scope LRC target on
that ray.

## 7. Replay

```text
python3 -B 04-computation/lrc14_owner_endpoint_rational_tail_thm4031.py
python3 -B -O 04-computation/lrc14_owner_endpoint_rational_tail_thm4031.py
```

Both runs reproduce the frozen raw-LF output. **QED.**
