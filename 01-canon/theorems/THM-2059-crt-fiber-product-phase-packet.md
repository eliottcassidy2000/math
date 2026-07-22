---
id: THM-2059
title: "CRT fiber-product formula for one-tail safe-phase packets"
status: >
  PROVED. For every scaled core aC with one tail w and every clock N, the
  safe conditions at t=k/(Na) form a CRT fiber product of a core packet modulo
  N and a tail packet modulo Na/gcd(w,Na). A residue-histogram dot product
  counts compatible classes modulo their lcm; multiplying by the lift factor
  Na/lcm counts safe grid indices. The dot product has an exact zero-mode plus
  finite-Fourier-fluctuation split and a Cauchy certificate. This gives
  an exact phase-packet carrier at arbitrary moduli and packages the
  missing-clock step used in THM-2057; it is a certificate generator, not
  LRC(14).
source: codex-2026-07-21-LRC-CRT-packet
script: 04-computation/lrc_crt_phase_packet_codex_20260721.py
result: 05-knowledge/results/lrc_crt_phase_packet_codex_20260721.out
script_sha256: dd42a0c3369bfede13ddda1871ab6e1db3280a879f2e735093e4713c0010ae7c
result_sha256: 922a12e188a35b607d8e760fd03b001e0054159a26b6e7cb2c20a52252b2c26e
hash_basis: normalized repository blobs (LF)
depends_on: []
related:
  - THM-1002
  - THM-2047
  - THM-2057
  - THM-2053
  - THM-2058
  - HYP-8846
  - HYP-8871
  - MISTAKE-233
---

# THM-2059 -- CRT fiber products are the exact one-tail clock carrier

Let `C` be a finite set of positive integers, let `a,w,N` be positive
integers, and let

```text
S=aC union {w},        Q=Na,
g=gcd(w,Q),            h=Q/g,
u=w/g.                                                     (1)
```

Thus `gcd(u,h)=1`. For an integer `x`, write

```text
|x|_m=min(x mod m, m-(x mod m)).                           (2)
```

Define the exact threshold packets

```text
A_N(C)={r mod N: 14|cr|_N>=N for every c in C},
B_h(u)={s mod h: 14|us|_h>=h}.                             (3)
```

Put `d=gcd(N,h)` and, for `j mod d`, define the reduction histograms

```text
alpha_j=#{r in A_N(C): r=j mod d},
beta_j =#{s in B_h(u): s=j mod d}.                         (4)
```

Then the number of packet types modulo `lcm(N,h)` is

```text
P_N(C;a,w)=sum_(j mod d) alpha_j beta_j.                   (5)
```

The number of safe grid indices `k mod Q` at phases `t=k/Q` is exactly

```text
(Q/lcm(N,h)) P_N(C;a,w).                                  (6)
```

In particular, a positive dot product in (5) proves `M(S)>=1/14`. It is an
LRC(14) instance certificate when `S` is a thirteen-speed relative row.

## Proof

At `t=k/Q`, a core speed has phase

```text
act=ack/(Na)=ck/N,
```

so all core speeds are at distance at least `1/14` precisely when
`k mod N` belongs to `A_N(C)`. The tail has phase

```text
wt=wk/Q=guk/(gh)=uk/h,
```

so it is safe precisely when `k mod h` belongs to `B_h(u)`.

The two residue requirements `k=r mod N` and `k=s mod h` are compatible
exactly when `r=s mod d`. By the generalized Chinese remainder theorem,
every compatible pair has one solution modulo `lcm(N,h)`, and every solution
arises from one such pair. Grouping compatible pairs by their common residue
`j mod d` gives (5). Since both `N` and `h` divide `Q`, their least common
multiple divides `Q`; each solution type has exactly `Q/lcm(N,h)` lifts
modulo `Q`. This proves (6). Every counted index supplies the displayed safe
phase, proving the certificate claim. QED.

## Relation to the missing-clock sieve

The missing-clock step inside THM-2057 is the especially rigid small-clock
case. If `2<=N<=14` and no
member of `C` is divisible by `N`, then every unit `r mod N` lies in
`A_N(C)`: each `cr` is a nonzero residue and hence has centered distance at
least `1`, while `1/N>=1/14`. THM-2057's central-unit lift proves that the
histograms in (4) overlap whenever `Na` does not divide `w`. Thus the
missing-clock conclusion `Na|w` for a hypothetical counterexample is a
nonemptiness application of the exact dot product (5). THM-2057's two full
plane closures additionally require their separate affine binding phases.

The gain is that (5) remains exact for `N>14`, when neither all nonzero
residues nor all units are automatically safe. Instead of discarding those
clocks, one computes the actual core packet and asks whether its reduction
histogram overlaps the tail packet. This is the natural join with THM-2058's
proved primitive phase-order packets and labelled owner intervals; THM-2059
adds compatibility, while neither theorem forces every overlap to be positive.

A zero overlap rejects only that clock grid, not the row.

## Exact bulk/fluctuation split

Write `A=sum_j alpha_j`, `B=sum_j beta_j` and

```text
alpha^0=alpha-(A/d)1,       beta^0=beta-(B/d)1.
```

Then (5) has the exact orthogonal decomposition

```text
P_N=AB/d + <alpha^0,beta^0>.                             (8)
```

Consequently

```text
P_N >= AB/d-||alpha^0||_2 ||beta^0||_2.                 (9)
```

Thus the clock is certified whenever the right side is positive. The test is
integer-exact after squaring:

```text
(AB)^2 > (d sum alpha_j^2-A^2)(d sum beta_j^2-B^2).      (10)
```

With the Fourier convention
`alpha_hat(r)=sum_j alpha_j exp(-2 pi i rj/d)`, Parseval gives

```text
P_N=(1/d) sum_(r mod d) alpha_hat(r) conjugate(beta_hat(r)),
||alpha^0||_2^2=(1/d) sum_(r!=0)|alpha_hat(r)|^2,         (11)
```

and similarly for `beta`. Hence (8) has a nonnegative zero mode, positive when
both packets are nonempty, plus a signed/complex nontrivial-mode sum. This is
the rigorous remnant of the
incoming “bulk plus cusp” intuition: it is finite Fourier analysis on the
actual CRT carrier, not a modular-form identification. If `P_N=0` with
`A,B>0`, the two nonnegative histograms have disjoint support, and the
nontrivial modes cancel the entire zero mode exactly.

## Assumption challenge and tournament analysis

The vertices are not runners, modular cusps, or primes. They are residues in
two different finite packets, and the decisive relation is the symmetric
compatibility edge

```text
r ~ s  iff  r=s mod gcd(N,h).                             (7)
```

Quotienting either packet to its total cardinality destroys the obstruction:
two nonempty packets can occupy disjoint reduction classes and have dot
product zero. The lossless quotient is the pair of histograms `(alpha,beta)`.

A tournament is not a clean carrier here. Relation (7) is bipartite and
symmetric, with unavoidable ties inside each reduction class. Artificially
orienting it would discard exactly the CRT compatibility being measured.
For comparing candidate clocks one may rank the scalars `P_N`, breaking ties
by strict margin or denominator cost, but that ranking is a search scheduler,
not the proof object. The exact proof object is the compatibility graph and
its histogram dot product.

## Computational audit

The companion uses integer arithmetic only. It compares (6) with direct grid
enumeration on `53760` rows from four unrelated core families, checks the
small missing-clock specialization, audits (10), and exhibits strict `N>14`
packet certificates. The frozen output ends in `PASS`.
