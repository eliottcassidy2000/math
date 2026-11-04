# Certified leading digits of tournament class counts

**PROVED finite remainder bound and exact prefix certification / FINITE-EXACT
implementation controls. This computes requested leading bits of A000568,
not its entire integer at large orders and not the Hamiltonian spectrum.**

The existing full-integer response DP is retained as the recommended engine
when the complete integer is needed. The new interface answers a different
query: determine the exact binary length and first B binary digits of a(n),
with an explicit rational certificate and work governed by the requested
precision. It avoids constructing a(n) or its quadratic-length identity power.
For n=10000 it certifies a binary length of49,876,542 and256 leading bits
from the identity permutation alone.

## 1 Inheritance and the live concepts

The closest proved mechanism is [THM-514 unit-tail kernel](../../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).
[THM-516 divisor features](../../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md)
and the [Session87 degree-bucket engine](../../session87_four_vertices_enumeration_20261028/FAST_COUNTS.md)
already supply future-gcd response compression, exact integral arithmetic,
and full prefixes through200. Their successor
[Session89 fused channels](../../session89_consecutive_tournament_pairs_20261028/counting/README.md)
already shares multiple counting queries. No such mechanism is claimed new.

The least-used relevant sidecar was the explicit boundary in
[Session98 all-poles expansion](../../session98_tournament_growth_20261029/limits/ALL_POLES_CORE_EXPANSION.md):
asymptotic truncation does not license exact rounding without an absolute
error bound after restoring2^binom(n,2)/n!. The present bound supplies that
missing finite certificate. Its cycle-core expansion is inherited, as is
[THM-4845 Gaussian and meromorphic synthesis](../../../../01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md).
The public sequence and classical Davis formula were checked at
[OEIS A000568](https://oeis.org/A000568) on2026-11-03. No external novelty or
priority claim is made for this elementary bound.

The concept board has four entries:

| Concept | Operation and retained information | Hostile or decisive test |
|---|---|---|
| Exact full integer | Inherited response DP retains all future gcd interactions | Compare every requested prefix against its exact output |
| Leading-bit query | Keep low moved-mass cores and a rigorous positive remainder | Same normalized error may be enormous in the raw integer |
| Modular tail query | Keep an independently proved residue and the interval width | Two integers separated by the modulus must not both fit |
| Hamiltonian spectrum | Retain H on each actual tournament, not only total orbit mass | Many asymmetric classes do not force a prescribed prime H |

The corrected near miss is rounding the identity contribution to get the full
count. At n100 its error exceeds1 by an enormous margin, even though the
identity alone nearly determines the leading digits. The checker rejects
that inference explicitly.

## 2 A finite positive remainder bound

Let Q=binom(n,2), a_n=A000568(n), and f_n=n!a_n/2^Q. Write an odd cycle type
as a non-singleton core mu and a tail of fixed points. For the core let

    m=sum(parts),  t=number of parts,
    E=sum_i (ell_i-1)/2 + sum_(i<j) gcd(ell_i,ell_j),
    z=product_ell ell^multiplicity * multiplicity!.

All parts of mu are odd and at least3. The empty core is permitted.
The inherited exact formula is

    f_n = sum_(m<=n, cores mu)
              (n)_m/z * 2^(-D_n(mu)),
    D_n(mu)=binom(m,2)-E+(m-t)(n-m).                 (1)

Here (n)_m=n!/(n-m)! and the empty-core term is1. Let S_(n,M) be the sum
restricted to m<=M, with0<=M<=n.

**Proposition C1.** Set c=floor((n-1)/3), r=n/2^c. If r<1 and M<n, then

    0 <= f_n-S_(n,M) <= r^s/(1-r),                 (2)

where s=max(3,M+1), improved to s=5 when this would give s=4.
For M=n the remainder is exactly zero. In particular r<1 for every n>=13.

**Proof.** A fixed-point-free odd permutation has every unordered-pair orbit
of size at least3. A fixed pair would require either two fixed vertices or
a2-cycle, and all its other orbit sizes are odd. Therefore

    E <= binom(m,2)/3,       t <= m/3.

Consequently

    D_n(mu) >= m(m-1)/3 + 2m(n-m)/3
              = m(2n-m-1)/3
             >= m(n-1)/3 >= mc.                  (3)

For each mass m, m!/z counts permutations of that cycle type. These types
are a subset of S_m, so sum_(mu of mass m)1/z<=1. Since(n)_m<=n^m, the
entire mass-m contribution in(1) is at most r^m. Summing the geometric series
over m>=s proves(2). No such core has mass1,2 or4, which proves the stated
improvement. For r<1 at n>=13, check n13,14,15; increasing n by3 doubles the
denominator while n+3<2n. QED.

This is intentionally conservative. It uses no assumed asymptotic regime,
no unverified decimal estimate, and no cancellation. For arbitrary n the
finite endpoint M=n remains an exact fallback.

## 3 Exact prefix certification

For an integer shift e>=0, form the rational interval

    L=2^(Q-e) S_(n,M)/n!,
    U=2^(Q-e) [S_(n,M)+R_(n,M)]/n!,               (4)

where R is(2), or0 when M=n. Then L<=a_n/2^e<=U. If floor(L)=floor(U)=q,

    q = floor(a_n/2^e)                            (5)

is exact. If q>0, its bit length plus e is the exact bit length of a_n.
Discarding excess low bits of q gives any smaller requested leading prefix.

`heads.py` takes B requested bits and initially uses

    e=max(0,Q-bitlength(n!)+1-B).

It grows M only until(5) holds. This choice makes the explicitly constructed
power2^(Q-e) have at most bitlength(n!)+B bits when e>0, instead of Q bits.
Small n use the exact finite endpoint. A supplied core-mass cap raises an
error if no certificate has been obtained; it never silently returns an
approximation. No floating-point logarithm occurs.

At large n, M=0 can suffice even for a long requested prefix, provided the
identity quotient is farther than the certified remainder from the next integer. The program
checks that distance; it does not infer it from small error alone. The
rational endpoints, binary head, cutoff and arithmetic sizes are returned.

**Connection contract.** The source is the exact odd-cycle Burnside sum;
the map forgets cores of mass>M and replaces their weight by a proved positive
interval. It preserves any binary prefix separated from the next integer by
that interval. It destroys the remaining low digits and all class identities.
The indispensable sidecar is the error interval after scaling, not merely
its normalized size. The cheapest test is equality of the two integer floors.

## 4 Exact reconstruction with an independent residue

**Proposition C2.** Suppose a rational interval[L,U] is proved to contain an
integer a and a residue a=r mod m is independently proved. If U-L<m, exactly
one integer with that residue can lie in the interval. It is

    a = r + m ceil((L-r)/m).                      (6)

The candidate must lie in[L,U]. This follows because distinct integers in
the same residue class differ by at least m. The implementation rejects an
interval of width>=m and rejects an inconsistent candidate.

This connects the head certificate to the inherited exact modular response
engine, but no faster full-integer hybrid is claimed. The bounded controls
supply residues from an independently computed exact integer. Obtaining
those residues at larger orders remains genuine work. A residue guessed
from the spectrum conjecture is not allowed.

## 5 Exact computation and measured work

The following are exact head certificates. The128-bit rows with n<=200 are
also compared against independently produced full integer counts.

| n | Exact binary length of a_n | Requested bits | Maximum retained core mass | Core types retained |
|---:|---:|---:|---:|---:|
|32|379|128|25|142|
|40|621|128|16|32|
|64|1721|128|8|6|
|100|4426|128|3|2|
|128|7412|128|3|2|
|160|11775|128|0|1|
|200|18655|128|0|1|
|1000|490971|256|0|1|
|10000|49876542|256|0|1|

At n10000 the certified256-bit prefix is

    e7d782997a96b31a0cb2523c988f010034622aef43f5aa98b08af4148e3c84b7

in hexadecimal. The exact n! has118,459 bits and the identity power used in
the scaled calculation has118,715 bits. This is not the full49,876,542-bit
class count and is not an extension of the OEIS full-term table.

`benchmark.json` records five alternating CPU measurements, with100 head
queries per timing batch at moderate orders. The baseline is the current
Session87 exact full-prefix engine, followed by extracting the same head.
These are different output workloads: the baseline also produces every
complete integer through n. No speedup ratio for complete A000568 computation
is asserted. At n1000 and10000, only the certified head is timed; no full-count
baseline is claimed feasible or run there.

Median CPU seconds from that same-host run:

| n | Full prefix CPU seconds | Certified128-bit query CPU seconds |
|---:|---:|---:|
|32|0.000792|0.00376769|
|64|0.005514|0.00013303|
|100|0.029107|0.00003509|
|160|0.330404|0.00002563|
|200|1.470206|0.00003249|

At n32 the head interface is slower: requesting128 of379 total bits makes
the conservative bound retain many core types. The full-prefix engine remains
the stronger choice there. At n10000 the256-bit query took median
0.032301 CPU seconds; this is a measured finite query, not a
full-count computation or a runtime guarantee.

## 6 A complete integer checkpoint at order256

The inherited degree-bucket engine was also rerun through256, and the entire
value a(256) is saved in hexadecimal in[full256.json](full256.json), with a
SHA-256 of its minimal big-endian binary encoding. It has30,957 binary digits
and two-adic valuation131. Its independently certified256-bit head agrees.
A separate inherited bounded-cycle calculation gives

    a(256)/2^128 = 4264 mod65536.

That check visits542,134 odd partitions having2,4 or6 cycles, of which97,444
survive the16-bit valuation cutoff. It uses the
[October25 fixed-precision formula](../../port_residual_integral_20261025/precision/README.md),
not the future-response DP. The separate explicit four-bit even-order formula
also gives8 modulo16. Normal and optimized complete-checkpoint receipts agree.

This is a reproducible full-integer checkpoint, not an enumeration record:
the[earlier computation history](../../strong_attachment_tail_compiler_20261028/HISTORY.md)
already documents complete exact scalar computation through300. No speedup
for full-value computation is claimed. Reproduce with

```sh
python3 -B 05-knowledge/results/orbit_spectrum_session_20261103/counting/full_checkpoint.py > /tmp/full256.json
```

The script has a120-second wall-time cap and verifies the source hash of the
inherited full engine. Its full integer SHA-256 is

    2b415fefb721acef1bd9e740d4d53b57b4d6aa42db0921b33a0a43d0f99a8ef3

## 7 Why this does not prove spectrum completeness

For a uniformly labelled tournament, f_n equals the expected automorphism
group order, by Burnside's double count. The truncation controls aggregate
symmetry mass. Since nontrivial tournament automorphism groups have odd order
at least3, the probability of nontrivial symmetry is at most(f_n-1)/2.
This elementary consequence concerns a random labelled tournament.

The inherited prime-H asymmetry theorem says that every witness with prime
H>3 is asymmetric. Its converse is false, and an abundance of asymmetric
classes does not imply a witness at a prescribed H. In particular(1)--(6)
supply neither the prime/7p constructor nor the missing pointwise support
needed to prove that7 and21 are the only positive odd forbidden values.
For that task a count refinement must retain the Hamiltonian statistic.
The asymmetry input is credited to [Session142 prime-lane proof](../../session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md), including
its sharper7p symmetry classification; it was recovered rather than rederived
as a new result here.

## 8 Reproduction and controls

From the repository root:

```sh
python3 -B 05-knowledge/results/orbit_spectrum_session_20261103/counting/check.py > /tmp/heads-normal.json
python3 -B -O 05-knowledge/results/orbit_spectrum_session_20261103/counting/check.py > /tmp/heads-optimized.json
cmp /tmp/heads-normal.json /tmp/heads-optimized.json
python3 -B 05-knowledge/results/orbit_spectrum_session_20261103/counting/heads.py 10000 --bits 256 --max-core-mass 0
```

The deterministic receipt covers804 exact head queries at n0..200 with
B=1,16,64,128; independent rational Davis sums through40;973 remainder
checks;25 literal unordered-pair cycle-type controls through8;68 exact
head/residue reconstructions; and the n1000/10000 rational certificates.
Positive controls include the exact full endpoint at small orders. Hostile
controls reject insufficient mass caps, identity-only full rounding, and
nonunique modular intervals. Runtime checks remain enabled under`-O`.
The separate benchmark receipt contains nondeterministic CPU measurements.
