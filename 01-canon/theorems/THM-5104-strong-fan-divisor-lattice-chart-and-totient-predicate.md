---
id: THM-5104
title: Strong fan divisor-lattice chart and totient predicate
status: PROVED elementary synthesis of inherited fans and classical repunit arithmetic; independently audited. No historical-priority claim.
source: 05-knowledge/results/session241_prime_relations_20261110/repunit_chart/PROOF.md
session: Session241, 2026-11-10
depends_on:
  - THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity
  - THM-4801-collatz-meets-the-constructible-polygons-trunk-five-initial-segment-principle-collatz-closed-family
---

# THM-5104: a divisor-lattice chart with actual strong tournaments

Let Q_1=C3. For n>=2 take the transitive chain1,...,2n+1 and add vertex0,
which loses exactly to n,n+1,2n+1. The inherited Fibonacci fan theorem
THM-4748 gives h_n=H(Q_n)=4^n−1. Q_n is strong; its order is3 for n=1
and2n+2 for n>=2.

For all positive a,b,n:

- gcd(h_a,h_b)=h_gcd(a,b), and h_a divides h_b iff a divides b.
- phi(n)=#{1<=a<=n:gcd(h_a,h_n)=3}.
- h_(a+b)=h_a*h_b+h_a+h_b.
- For every odd prime p, p divides h_n iff ord_p(4) divides n. If
  n=ord_p(4)*t, then v_p(h_n)=v_p(4^ord_p(4)−1)+v_p(t).

Thus the image with numerical divisibility is exactly the positive-integer
divisor lattice, with bottom3 and ordinary-gcd meet. Its internal join is
h_lcm(a,b), and its incidence Möbius function and divisor counts transport
from the index lattice. The map does not require factoring n to construct Q_n.

The Euclidean exponent step
(4^b−1)−4^(b−a)(4^a−1)=4^(b−a)−1 proves the gcd identity. Strict monotonicity
gives the divisibility equivalence; the totient predicate follows at once.
The other identities are expansion and classical order/LTE. The inherited
THM-4801 initial-run criterion also gives phi(h_n) a power of2 exactly for
n in{1,2,4,8,16}, unconditionally. This last function is not phi(n).

## Losses are part of the theorem

The normalized R_n=h_n/3 also preserves divisibility and gcd, but is only a
coordinate: R_3=21 is not a Hamilton-path count. Ordinary lcm is not preserved:
lcm(h_2,h_3)=315, while h_6=4095. Ordinary multiplication is not preserved:
h_2*h_3=945. Numerical primality is not preserved: prime index3 maps to63.
An explicit induced Q_2 inside Q_3 shows raw containment is not divisibility.
The chart has only floor(log_4(X+1)) values up to real X>=0. The shifted addition
identity is not a local gluing operation for arbitrary tournaments.

[Complete proof, inheritance and hostile witnesses](../../05-knowledge/results/session241_prime_relations_20261110/repunit_chart/PROOF.md)
and [independent audit](../../05-knowledge/results/session241_prime_relations_20261110/embedding_audit/AUDIT.md)
separate all-index reasoning from finite exact controls. This is an arithmetic
representation inside the spectrum, not Chan–Pak surjectivity. Chan–Pak OPEN.
