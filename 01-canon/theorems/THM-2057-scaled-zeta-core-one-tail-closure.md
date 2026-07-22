---
id: THM-2057
title: "The scaled missing-clock sieve closes two AP one-tail planes"
status: >
  PROVED. General sieve: if a core C contains no multiple of N, 2<=N<=14,
  then a scaled one-tail counterexample aC union {w} must have Na|w. Hence it
  pays the lcm of all missing clocks. For C={1,...,11,13}, the missing clocks
  are exactly 12 and 14. Their lcm forces 84a|w, where the explicit scaled
  affine binding witness below is strict. Thus every row
  {a,2a,...,11a,13a,w} satisfies LRC(14). For C={1,...,12}, the missing
  clocks 13 and 14 force 182a|w; the explicit deep-well phase below closes
  every row {a,2a,...,12a,w}. This closes two rank-two planes, not LRC(14).
source: codex-2026-07-21-LRC-scaled-zeta-core
script: 04-computation/lrc_kelvin_farey_scaled_core_codex_20260721.py
result: 05-knowledge/results/lrc_kelvin_farey_scaled_core_codex_20260721.out
script_sha256: 9f4d2780c6d9d24dfcdf8e99ae1a1fe327ad98e97e21044effd95769c0e4c309
result_sha256: 35c25efa343bf489f79c9593379378699483c254915f9d9c50a080d119c9b520
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2047
related:
  - THM-724
  - THM-1014
  - THM-2053
  - THM-2055
  - THM-2056
  - THM-2059
  - HYP-2896
  - HYP-8846
---

# THM-2057 -- scaled missing-clock sieve and two one-tail closures

Let

```text
C={1,2,...,11,13},
S(a,w)=a C union {w},          a,w in Z_(>0).              (1)
```

Repeated speeds, if any, impose no extra constraint, so the statement covers
both sets and labelled rows.

**Theorem.** For every positive integers `a,w`,

```text
M(S(a,w))>=1/14.                                         (2)
```

Equivalently, the whole integral rank-two plane

```text
a(1,2,...,13)+b e_12,       a>0, 12a+b=w>0,              (3)
```

satisfies LRC(14).

## 1. A modular unit-orbit lemma

We use the following elementary lemma.

**Lemma.** Let `2<=N<=14`, let `a,w` be positive integers, and put `q=Na`.
If `q` does not divide `w`, there is an integer `k` such that

```text
gcd(k,N)=1,                 ||k w/q||>=1/14.              (4)
```

*Proof.* Put

```text
g=gcd(w,q),       h=q/g>=2,       c=w/g,
```

so `gcd(c,h)=1`. Choose a unit `r mod h` away from the two ends of the residue
interval as follows:

```text
h<=14:             r=1;
h>14 odd:          r=(h-1)/2;
h>14, h=0 mod 4:  r=h/2-1;
h>14, h=2 mod 4:  r=h/2-2.                                (5)
```

In every case `gcd(r,h)=1` and

```text
min(r,h-r)>=h/14.                                        (6)
```

Choose `k_0 mod h` with `c k_0=r mod h`. It is a unit modulo `h`. We may
replace `k_0` by `k_0+h t` so that it is also coprime to `N`: if a prime
`p|N` divides `h`, every representative is already nonzero modulo `p`; if
`p` does not divide `h`, exactly one class of `t mod p` is forbidden. Choose
the allowed classes simultaneously for the primes dividing `N`.

For the resulting `k`,

```text
||k w/q||=||k c/h||=min(r,h-r)/h>=1/14,
```

which proves the lemma. QED.

The point of using a *unit* `r` is that no Jacobsthal, class-group, or
equidistribution input is hidden here. Formula (5) is an explicit central
unit for every modulus.

**Missing-clock corollary.** Let `C_0` be any finite set of positive integers
with no element divisible by `N`, where `2<=N<=14`. Then

```text
Na does not divide w  implies  M(a C_0 union {w})>=1/14.  (7)
```

Indeed, use the lemma's `t=k/(Na)`. Every core phase is
`||c k/N||>=1/N>=1/14`, because `k` is a unit modulo `N` and `c` is nonzero
modulo `N`; the tail satisfies (4).

Consequently, if

```text
Q(C_0)=lcm{N:2<=N<=14 and C_0 contains no multiple of N}, (8)
```

then any counterexample of the form `a C_0 union {w}` must pay the exact
divisibility tax

```text
a Q(C_0) divides w.                                      (9)
```

This is the transferable statement. It converts several surviving safe
clocks into one high-codimension divisibility ray before any geometry or
Euler argument is used.

## 2. The `12a` clock

Suppose first that `12a` does not divide `w`. Apply the lemma with `N=12` and
write `t=k/(12a)`. For every `i in C`,

```text
||a i t||=||i k/12||>=1/12,                              (10)
```

because `k` is a unit modulo `12` and none of the residues
`1,2,...,11,13=1 mod 12` is zero. The tail has distance at least `1/14` by
(4). Thus `t` proves (2).

This strictly extends HYP-2896's first branch. The old `a=1` proof only had
to observe `12 does not divide w`; the unit-orbit lemma supplies the missing
numerator when the core has arbitrary scale `a`.

## 3. The `14a` clock

Now suppose `12a|w` but `14a` does not divide `w`. Apply the lemma with
`N=14` and put `t=k/(14a)`. Since `k` is a unit modulo `14` and no element of
`C` is zero modulo `14`,

```text
||a i t||=||i k/14||>=1/14       for every i in C,        (11)
```

while the tail again satisfies (4). Hence (2) follows.

Notice that the route decision is divisibility by the *scaled* clocks
`12a,14a`, not by `12,14`. Projectivizing the parameter plane too early would
erase exactly this datum.

## 4. The double-kill binding family

It remains that both `12a` and `14a` divide `w`. Their least common multiple
is `84a`, so

```text
w=84a m,             m>=1.                               (12)
```

For the unscaled row `C union {84m}`, put

```text
t_m=(35m+2)/(84m+5),
D=84m+5.
```

The phase computation needed here is self-contained. In the order
`1,2,...,11,13,84m`, the numerator distances `D||v t_m||` are

```text
35m+2, 14m+1, 21m+1, 28m+2, 7m, 42m+2,
7m+1, 28m+1, 21m+2, 14m, 35m+3, 35m+1, 7m.
```

Every entry is at least `7m`, with equality at speeds `5` and `84m`.
Consequently the exact margin of this witness is

```text
min_(v in C union {84m}) ||v t_m||=7m/(84m+5)>1/14,      (13)
```

where the strict inequality is `14m>5`. This is a lower-bound certificate
for `M`, not a claim that the witness is globally maximizing.

For `S(a,w)`, use `t=t_m/a`. Then every phase is identical to its phase in
(13): `ai t=i t_m` and `w t=84m t_m`. Therefore the last branch is not merely
safe but strict. This completes the proof of (2). QED.

## 5. The adjacent AP-tail plane

The same sieve closes a second full rank-two plane. Put

```text
C'={1,2,...,12},
S'(a,w)=a C' union {w}.                                  (14)
```

The missing clocks of `C'` in `{2,...,14}` are exactly `13` and `14`. Thus
the missing-clock corollary proves `M(S'(a,w))>=1/14` unless both `13a` and
`14a` divide `w`. In the remaining branch,

```text
w=182a m,                 m>=1.                           (15)
```

Use the explicit phase

```text
t=14m/[a(182m+1)].                                      (16)
```

Write `Q=182m+1=13(14m)+1`. For `1<=j<=6`, the least residue of
`j(14m)` modulo `Q` is at least `14m`. For `7<=j<=12`, its complementary
residue is

```text
Q-j(14m)=(13-j)14m+1>=14m+1.                            (17)
```

Hence every core speed `aj` has distance at least `14m/Q`. Since
`182m=-1 mod Q`, the tail has distance exactly `14m/Q`. Finally,

```text
14m/(182m+1)>1/14  iff  14m>1.                           (18)
```

Therefore every row `S'(a,w)` satisfies LRC(14), with a strict witness on the
double-kill ray. This recovers the deep-well family identified in THM-724 and
THM-1014, but the displayed argument needs only the explicit phase and proves
the threshold statement for every scale and every tail. QED.

## 6. Exact route certificates

The proof is algorithmic and has only three leaves:

```text
12a does not divide w  -> central unit on the 12a clock;
12a divides w, 14a does not -> central unit on the 14a clock;
84a divides w          -> scaled affine binding phase (13).               (19)
```

The adjacent AP-tail plane has the parallel certificate

```text
13a does not divide w  -> central unit on the 13a clock;
13a divides w, 14a does not -> central unit on the 14a clock;
182a divides w         -> deep-well phase (16).                           (20)
```

No finite-height search is involved. The route also explains why the very
large THM-2053 determinant residual on this plane is misleading: the phase
certificate is controlled by two scaled clock orbits and one divisibility
ray, not by Euclidean parameter size.

## 7. Assumption challenge and transfer target

The useful vertices here are neither runners nor normal-fan sectors. They are

```text
safe numerator orbits,
killed clocks,
divisibility sublattices,
binding rays.
```

This quotient preserves the actual LRC witness margin and the route modulus;
it forgets endpoint-owner topology away from the one-tail core. The transfer
target for a general THM-2052 star is therefore: find a lower-rank core whose
safe residues contain a complete unit orbit on one or two clocks, then prove
that simultaneous clock killing forces a finite collection of affine binding
families. HYP-8871 records that sidecar program alongside the Kelvin/Farey
address from THM-2056.
