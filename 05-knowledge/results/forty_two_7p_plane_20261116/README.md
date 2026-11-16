# The 41/42/43 seam and the 7p plane

The 7p obstruction is not a lack of numerical formulas for multiplying
by seven. It is the difference between an invertible formula and its
physical range. Below are exact positive constructions, scope-limited
negative results, and the remaining target-dependent step.

## 1. Induced lifts across the center 42

The archived [ready restart packet](../ready_restarts_rounding_20261116/README.md)
supplies one synchronized 41/43 pair at center 42. The [birth-fibre
atlas](../birth_fibers_binary_growth_20261116/README.md) supplies another
matched-root pair. Here are two **other** order-nine outputs. In each
case the prime parent is literally induced, so these are two-vertex
attachments without any change to the old arcs:

| Outgoing masks of strong output | H | common E | g |
|---|---:|---:|---:|
|256,273,35,279,357,395,303,351,4|287|105|77|
|256,269,281,97,299,327,407,319,8|301|105|91|

The root is vertex 8. O=E+g at both roots. The induced order-seven
parent counts 41 and 43, respectively; the replay prints eligible
subset masks. These are actual outputs at 7(42−1) and 7(42+1), centered
at 294=7·42. Appending a common near-sink word of all zeros gives
the entire exact pair of rays

    7(30·2^k+11), 7(30·2^k+13), k>=0.

The quotients are automatically coprime to 30. At k=1 they are the
twin primes 71 and 73; at k=6 they are 1931 and 1933. This does not
prove infinitude of twin primes. Since g differs, this is not an
arbitrary fixed-root binary continuation. See [another induced pair](RAIL.md)
with common E=217, which starts at the same 41/43 but follows
quotient rails 62·2^k−21 and 62·2^k−19. Those induced parents themselves
have common E=27, not 31. If parent and child continue in parallel the
literal sevenfold scaling holds at the anchor only; its drift is
56(2^k−1). This distinguishes an exact starting lift from an
operator intertwining for every height.

### Exact inverse with a changed right source

Append an arbitrary uniform near-sink block B to the first pair.
With its actual two-count profile (u,v)=(H(B),P_B(2)),
the normalized right targets are

    p=30v+11u, or p=30v+13u.

For prescribed p, solve u≡11p or 7p modulo 30 respectively
(here 7p on the right is multiplication in the congruence for u,
not a second graph target), then test v=(p−11u)/30 or
(p−13u)/30. Since v≥2u, the unknown u is at most p/71 or p/73.
This is a target-decreasing, physical-profile obligation with
modulus 30 instead of 108·2^s in the old 756 route. For p>5
a prime, u must be coprime to 30; it does not have to be prime.
For example the seven-vertex source (2,20,17,7,105,15,47)
has profile (49,872). With one transitive pad (u,v)=(49,1744),
the two quotient primes are 52859 and 52957.
The E=217 pair similarly reads (62v−21u,62v−19u);
the literal example with (23,272) yields prime factors
16381 and 16427. These specify exact witnessed profiles;
they do not assert that the progressions always intersect the
physical profile set.

## 2. Two exact barriers in the old reader

**Same-block prime-gap obstruction.** On the old
756·2^s v+Iu row with I in {−33,−31,−29}, divisibility by
7 forces 7|u, since all three intercepts are units modulo 7.
If the normalized quotient is prime >3 then 3∤u, since
the leading term is divisible by 3. A positive such count u
must be at least 35 because 7 is forbidden. Consequently
two simultaneous prime quotients from the same block and
height cannot be twins: the smallest possible difference is
10. For adjacent intercepts the difference is 2u/7; for
outer intercepts it is 4u/7. The bound is sharp. With
(u,v)=(35,676) from literal rows (16,17,3,71,12,31,55),
s=2 and adjacent intercepts −31 and −29 give the prime
quotients 291877 and 291887. The 41/43 constructions above
change the core and so are not subject to this barrier.

**No uniform sevenfold amplification after normalization.**
Suppose a count p≥5 prime was itself produced by a nonempty
near-sink step and is at the same root, in its normalized
profile (p,g,E,E+g). For an active normalized root the
[domination bound](../birth_fibers_binary_growth_20261116/README.md)
is p≤4E−1. Also E<p after a near-sink step: from the previous
profile H'−E'=H−(O−E)+g≥H, since O−E is the original full
root-end count, which injects into the root deletion paths.
If E=0, amplifying needs H(B)=7, impossible. Otherwise a
nonempty uniform near-sink block (u,v) could give 7p only if

    (7−u)p = 2(v−u)E.

Because v≥2u, u is 1, 3 or 5.
- u=1 gives a transitive block with v=2^m. The four possible
  divisors 2^m−1 of 3p give E=3p,p,3,1. The first two
  violate E<p, the last violates p≤4E−1. If E=3, the only
  prime p≥5 with p=2^m−1 and p≤11 is 7, not a source.
- u=3 means a triangle and transitive padding, by the
  small strong count bound and SCC factorization.
  Then v=12·2^a and 3 divides v−3, contrary to
  (v−3)E=2p.
- u=5 implies the unique order-four strong H5 core and
  padding, so v=32·2^a. Since (v−5)E=p and v−5>1,
  E=1, again impossible.

Any sequence of such uniform blocks is one ordinal
concatenated block. This rules out a whole all-order
triangular architecture, not the nonuniform induced lifts.

## 3. From 21 to 63 and back to primes

63=3·21=9·7 is not a target 7p with p prime. Its nontrivial
SCC cannot be built as a product of smaller counts: every
nontrivial factorization uses 7 or 21. Yet 63 has strong
witnesses, first born on eight vertices. The two birth cores
are a converse pair, replayed explicitly here.

For one of these cores (2,16,3,23,165,79,31,111), use root
1. The actual profile is (63,11,26,37). After k singleton
near-sinks the count is 63+52(2^k−1). Since 2 has order three
modulo seven, the 7-lane returns exactly at k=3j. The quotient

    P_j=9+52(8^j−1)/7

is a 7-adic isometry on j: v_7(P_i−P_j)=v_7(i−j), by
lifting the exponent on 8^d−1. It gives a permutation of
the residues modulo every 7^a and a unique modular inverse
at every level. At the first return it gives the actual prime
quotient 61, i.e. the physical count 427. Nevertheless it
never gives quotient 2 modulo 3: for even j it is 0, for
odd j it is 1. So all 7-adic phases coexist with a missing
ordinary residue class and exponentially spaced absolute
outputs. The second known normalized 63 root (E=38) has
the same local features, with first return 7·85.

**Selective channel, exact up to 301.** Below 630, any
strong extension of either 63 core can be built by single
mixed attachments, each keeping the current graph strong.
If no remaining single vertex has mixed contact with a
current strong S, S is a module. Its contraction inside a
strictly larger strong graph lies on a directed triangle, so
the larger graph contains C3[S,1,1]. This costs at least
3H(S)+Q2(S)≥(|S|+2)H(S)≥630: cutting each Hamilton path
once witnesses Q2(S)≥(|S|−1)H(S).
That proves the mixed-only grammar in this extended window.
A tournament with count 7p, p≥5 prime, has a single strong
nontrivial SCC: the other possible factorization needs H7.
Neutral singleton pads do not change the count, and the old
strong 63 core lies in that SCC. Thus the grammar suffices
for all such extensions at every order.

Replay mixed cuts to 301, using the weighted even-path cut
law of [the reader-price packet](../exact_reader_price_catalan_20261116/README.md).
One core has labelled layer sizes 1,58,266,0, and converse
gives the other core. The exact prime quotients attainable
from these cores up to 43 are 19,23,29,31,37,41,43.
In particular 41 requires two mixed attachments here, 43
just one, and 11,13,17 cannot be recovered from these cores
at *any* order. This is not an exclusion of those totals
from the global tournament spectrum.

## 4. The planar versus 3D Jacobian comparison, precisely scoped

The complex planar conjecture is still open.
[Gao's November revision](https://arxiv.org/abs/2608.00222v2)
gives explicit constant-Jacobian noninjective maps in three
and higher dimensions and describes the tangent-sweep
construction and nonproperness. The repo has the exact
degree-seven three-to-one example in THM-1300, the non-slice
test in THM-4562, and the twisted planar quotient with
a ramification factor in THM-4727. The latter is not a map
of the planar Keller class; the former does not descend
through a simple slice. Earlier [quintic source flow work](../planar_quintic_s5_petersen_resolvent_20261019.md)
likewise needed the labelled fibre, not merely a scalar invariant.

The usable analogy is a range issue. The map
(u,v)↦(30v+11u,30v+13u) has nonzero determinant (absolute
value 60); it can be inverted algebraically. But an inverse
lattice point need not be a physical two-colour profile.
Neither algebraic invertibility nor local residue closure
implies physical surjectivity. The nonuniform two-vertex
lifts are a concrete route outside the uniform triangular
restriction. No implication about JC2 is being claimed.

Extra number trail: [the modulus-43 six-step clock](../six_reset_profile_progressions_20261116/README.md),
where 43 is the only odd prime permitting a certain period-six
collapse; [the (2,3,7) Hurwitz deficit 1/42](../ogg_triangular_chowla_20260927.md);
and [the earlier 41/59 Frobenius edge defect](../session32_square_monster_20261023/SYNTHESIS.md).
These are separate marked-fibre phenomena, not a derivation
of one number from another. Run check_rail.py and check_63.py
with standard Python.
