---
id: THM-4501
title: "Recursive positive Collatz motif families with exact natural frequency and dyadic scaling"
status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED
source: crossroads-family-20260926, kind-pasteur
depends_on:
  - 01-canon/theorems/THM-4497-coprime-graph-rank-compression.md
related:
  - 05-knowledge/results/crossroads_family_20260926_orbits.md
  - 05-knowledge/results/crossroads_family_20260926_fractal.md
scripts:
  - 04-computation/experiments/crossroads_family_20260926_orbits.py
outputs:
  - 05-knowledge/results/crossroads_family_20260926_orbits.out
audit: "Root and geometry independently checked affine realization, CRT support, valuation disjointness, infinite-union density, exact phase period, and dyadic closure losses. Ordinary and optimized Python outputs match."
---

# THM-4501 — Recursive motif families and their occurrence rate

**PROVED.** Full proofs and exact controls are in
[the orbit-family note, sections 2–7](../../05-knowledge/results/crossroads_family_20260926_orbits.md).
T denotes the shortcut Collatz map n/2 or (3n+1)/2. U denotes its odd
accelerated map. These are genuine positive integer orbit families;
they do not assert an infinite positive no-descent orbit.

## Explicit family

Put a=4347, P=7*23*29=4669, M=P*2^12=19124224, and

    r_k=4348*3^(-k) mod M,
    F_k={2^k u-1: u>=4348, u=r_k mod M}, k>=0.        (1)

Every n in F_k has its first k+11 shortcut iterates strictly above n.
After the first k odd steps it enters a+M Z_(>=0). Its next ten odd
nodes retain a prime-support triangle on positions 0,3,6: the three
edge primes 7,23,29 each divide exactly the corresponding pair and no
other node in those ten positions. The base ten-node word is

    4347,6521,4891,7337,5503,8255,12383,18575,27863,41795.

The families F_k are pairwise disjoint. For every fixed K>=0,

    #((union_(k>=K) F_k) intersect [1,X])
        =2^(1-K) X/19124224+O(log X).                 (2)

The implied constant can depend on the fixed base modulus. The least
period of the phase sequence r_k is tau=59136, and

    F_(k+tau)={2^tau(n+1)-1: n in F_k}.               (3)

Thus (1) gives an exact arithmetic recursion, an unbounded range of
finite growth lengths, and an exact occurrence frequency.

## Proof and retained coordinates

For g=a+2^12 t the i-th odd node, through shortcut time S_i<=11, is
m_i(g)=m_i(a)+3^i 2^(12-S_i)t. Taking t a multiple of P preserves all
residues modulo the triangle primes. Every nonempty base prefix has
multiplier at least 9/8, so it stays above its positive source. Exact
support is certified on the listed finite base; the affine formula
propagates it to every member of the progression.

For n=2^k u-1, the first j<=k shortcut iterates are
3^j 2^(k-j)u-1. They are odd where a next step is required and strictly
increase. The residue condition makes their landing point 3^k u-1 a
member of a+M Z_(>=0). Also v_2(u)=2, so v_2(n+1)=k+2, proving
disjointness. If u_k is the least admissible u, then
4348<=u_k<=4347+M and the exact k-th count is

    max(0,1+floor((X-(2^k u_k-1))/(2^k M))).

Only O(log X) levels contribute; summing bounded rounding errors
and the geometric comparison tail proves (2). Summing densities without
this tail control would be insufficient.

The residue clock has period ord_(M/gcd(4348,M))(3), namely
ord_(2^10*7*23*29)(3)=lcm(256,6,11,28)=59136.
The common cutoff on u then proves the exact set equality (3).

## Closure is not natural density

The dyadic closure of F_k is

    C_k=2^k r_k-1+2^(k+12) Z_2.

The closure of the tail union is {-1} union union_(k>=K) C_k, with Haar
measure 2^(1-K-12), **4669 times its natural density**. The odd CRT
factor vanishes in dyadic closure. Each such closure has Hausdorff
dimension 1, while their intersection as K tends to infinity is {-1},
of dimension 0. The integer sets themselves have empty intersection.

Thus arbitrarily long finite realization, exact recurrence, positive
fixed-horizon density, and large closure dimension do not supply one
positive integer realizing every horizon.

## Relation to 27

The literal positive ancestors of 27 are only {27*2^r:r>=0}; an odd
predecessor of a multiple of 3 is impossible. The more useful common-tail
comb consists of all positive odd direct predecessors of 41:

    n_j=(82*4^j-1)/3=27,109,437,1749,... .

It obeys n_(j+1)=4n_j+1 and has zero natural density. Its longer total
stopping time does not imply longer initial growth: every j>=1 descends
within two shortcut steps. A separate version of (1) based on 27 retains
52 positive-slope steps and the selected triangle 91,175,325; its primes
are exclusive only within that selected triple, not the whole prefix.

For every positive target y coprime to 3 its basin is dense in Z_2:
solve 2^(L+r)y=C_w mod 3^e for any prescribed parity word, using that
2 generates the units modulo 3^e. The resulting positive source
(2^(L+r)y-C_w)/3^e realizes the word and then reaches y. This is a
closure statement with no natural-density conclusion. In particular,
known convergent positive integers realize every finite parity prefix.
