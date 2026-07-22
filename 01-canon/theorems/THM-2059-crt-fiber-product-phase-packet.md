---
id: THM-2059
title: "CRT fiber-product formula for one-tail safe-phase packets"
status: >
  PROVED. For every scaled core aC with one tail w and every clock N, the
  safe phases t=k/(Na) are exactly a CRT fiber product of a core packet modulo
  N and a tail packet modulo Na/gcd(w,Na). Their number is an explicit dot
  product of residue-class histograms. This extends the phase-packet carrier
  to arbitrary moduli and recovers THM-2057 as a special nonemptiness lemma;
  it is a certificate generator, not LRC(14).
source: codex-2026-07-21-LRC-CRT-packet
script: 04-computation/lrc_crt_phase_packet_codex_20260721.py
result: 05-knowledge/results/lrc_crt_phase_packet_codex_20260721.out
script_sha256: 814deeb54a10e8632871015f6b4b84a7a982356c2525a66f8973e90fb0f53af0
result_sha256: f0d5303e07ea3c38cbe5e019821bc8d96dabc650f8e335fa3a72683e42aca822
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-1002
  - THM-2047
  - THM-2057
related:
  - THM-2053
  - THM-2058
  - HYP-8846
  - HYP-8871
  - MISTAKE-226
---

# THM-2059 -- CRT fiber products are the exact one-tail clock carrier

Let `C` be a finite set of positive integers and let

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

In particular, a positive dot product in (5) is an exact LRC(14)
certificate for the row `S`.

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

THM-2057 is the especially rigid small-clock case. If `2<=N<=14` and no
member of `C` is divisible by `N`, then every unit `r mod N` lies in
`A_N(C)`: each `cr` is a nonzero residue and hence has centered distance at
least `1`, while `1/N>=1/14`. THM-2057's central-unit lift proves that the
histograms in (4) overlap whenever `Na` does not divide `w`. Its conclusion
`Na|w` for a hypothetical counterexample is therefore a nonemptiness theorem
for the exact dot product (5).

The gain is that (5) remains exact for `N>14`, when neither all nonzero
residues nor all units are automatically safe. Instead of discarding those
clocks, one computes the actual core packet and asks whether its reduction
histogram overlaps the tail packet. This is the natural join with the
primitive phase-order counts claimed in THM-2058.

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
small missing-clock specialization, and exhibits strict `N>14` packet
certificates. The frozen output ends in `PASS`.
