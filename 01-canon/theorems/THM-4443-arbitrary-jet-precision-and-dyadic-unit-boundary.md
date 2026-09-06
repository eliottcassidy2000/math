---
id: THM-4443
title: "Arbitrary-jet precision and the dyadic unit boundary"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 fifth research wave
depends_on:
  - THM-4439-all-node-twojet-metric-precision-by-terminal-clusters
related:
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
  - THM-4435-four-node-metric-blindness-and-universal-hermite-precision
proof: 05-knowledge/results/hermite-higher-jet-unit-boundary-overnight-hexagon-sep05.md
script: 04-computation/hermite_higher_jet_precision_overnight_hexagon_sep05.py
output: 05-knowledge/results/hermite_higher_jet_precision_overnight_hexagon_sep05.out
script_sha256: 6263c2b9b5fc65036e7b5bd3bf2d751d529defb63635f07aa8beba9ee1a52a97
output_sha256: 0fa167e65476acbe873812aae7e7a52dd48a9451cd49d9a7ec0d3d5475cf11b6
independent_script: 04-computation/hermite_higher_jet_referee_overnight_hexagon_sep05.py
independent_output: 05-knowledge/results/hermite_higher_jet_referee_overnight_hexagon_sep05.out
independent_script_sha256: f14c3ac44f34cbbe275408895de781eb6363d7586fd03e4fffae6c86a095f697
independent_output_sha256: c0c03411b4c78f6101d660930eb2a8e5bc3c7467abc34864d913a675c3fed1e7
hash_basis: raw LF bytes
---

# THM-4443 -- Arbitrary-jet precision and the dyadic unit boundary

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[complete proof, failure mechanisms and independent referee](../../05-knowledge/results/hermite-higher-jet-unit-boundary-overnight-hexagon-sep05.md)
are part of this theorem. No external priority is claimed.

For distinct integer nodes x_i with arbitrary positive multiplicities m_i,
observe all Hasse jets of orders0,...,m_i-1 on degree<sum_i m_i polynomials.
Set

```text
Q_i(X)=product_(j!=i)(X-x_j)^(m_j),
a_(i,l)=[T^l]Q_i(x_i+T)^(-1), 0<=l<m_i.
```

The largest integer Smith factor is the least common multiple of the
reduced denominators of all a_(i,l). Its p-exponent is the sharp extra
data precision L_p=max_(i,l)(-v_p(a_(i,l))). The cardinal Hermite inverse
column for datum(i,r) has top coefficient exactly a_(i,m_i-r-1), so every
proposed denominator is attained, not only an upper bound.

For three nodes with uniform multiplicity3 at p=2, write the pair depths
as(e,e,e+d), d>=1. If d>=2, L_2=8e+5d-1. If d=1, choose closest pair
x,y and outsider z, and put t=2(z-x)/(y-x) modulo16. Then

```text
L_2=max(7e+4,8e+Gamma(t)),
Gamma=3 for t=3 mod4; 2 for t=5 mod8;
      1 for t=1 mod16; 0 for t=9 mod16.
```

Exchanging closest endpoints replaces t by2-t and preserves Gamma; unit-
affine changes preserve it as well. Simultaneous cancellation of
t^2+3t+4 and t^2-7t+14 has minimum valuation1,2,3,4 respectively.
All four branches are attainable. The maximum hides some distinctions
at shallow depth; the residue is not claimed necessary at every fixed e.

In particular the isometric triples2^e(0,1,2) and2^e(0,1,3) have losses
max(7e+4,8e+1) and max(7e+4,8e+3), different for every e>=2. The first
displayed full9x9 Smith lists have largest exponents18 and19. Uniform
k=1 and k=2 are metric-only at every node count, while two-node uniform
observers are metric-only at every k. Thus uniform k=3,n=3 is the first
possible multiplicity/node-count boundary, not a smallest-diameter claim.

Dropping uniformity fails earlier: weighted isometric sets2^e(0,2,1) and
2^e(1,3,0), each with multiplicities(2,2,1), have losses max(3e+2,4e+1)
and max(3e+2,4e), respectively9 and8 at e=2. This does not contradict
[THM-4439](THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md)'s
complete uniform two-jet hypothesis.

Primary controls: complete3,160-row diameter80 head, all four unit classes
at17 scales,1,000 signed/unit-affine pairs,49 literal full inverse/Smith
matrices,112 two-node controls,5,947 explicit gates. The independent
no-primary-import referee compares all largest integer factors on126
literal matrices, with1,262 explicit gates. Both normal/optimized pairs
match and the root replayed the referee. General higher-node full
partitions and other-prime higher-jet metric classifications remain OPEN.
