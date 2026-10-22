# Independent audit: converse packing congruences

Status: ACCEPTED. No corrections required to the frozen producer report.

Reviewed source: `session29-converse-congruences/REPORT.md`, SHA256 `66be30171d977ff4114f74db099eb7745f368ed5d16a0e08884122aa115696e6`; checker SHA256 `5d86197856c452b6a884579073025a580d0940ed1d0ce8dc077c1ee8014534ab`. This audit covers that report, not any later half-path extension.

## Proof review

1. The primary source was independently opened and Corollary20 inspected: Irving–Omar, *Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra*, https://arxiv.org/html/2412.10572v3 . It gives the Hamiltonian count of the complement as a sum over compatible permutations with only odd cycles, weighted by two per nontrivial cycle. The complement's loops cannot be used by a Hamiltonian path; its nonloop arcs are the converse tournament. Reversing paths gives the same H. Removing permutation fixed points yields exactly the odd-cycle packing formula used here. This is a cited classical input, not a new theorem of this packet.

2. The anti-involution argument is correct. Tournament automorphism groups have odd order because an involution would interchange endpoints of an arc. If a is anti and a² has odd order m, a^m is involutory and anti. The explicit conjugator for two anti-involutions is valid: with g=tau sigma and odd order m, h=g^((m+1)/2) obeys h sigma=tau h. There is at most one fixed vertex, so odd invariant cycle supports must contain that vertex.

3. A fixed packing can contain at most one individually fixed odd cycle. Every remaining member occurs with its partner, and such a pair is eligible precisely when its supports are disjoint. This proves Phi(x)=I(Q,x²)+x sum_C I(Q_C,x²). The producer correctly states J=Phi+2M with nonnegative integral M, rather than claiming the quotient graph Q determines all packing counts.

4. If A is the number of nontrivial cycle orbits, D counts internally disjoint ones, and I counts overlapping ones, then alpha_1=F+2A and alpha_2≡D mod2. Substitution in H≡1+2alpha_1+4alpha_2 mod8 gives H≡1+2F+4I. No fixed-point assumption is silently used for odd order; it is precisely what prohibits two individually fixed cycles in a two-cycle packing. Conjugacy makes F and I independent of the chosen anti-involution.

5. For even n, F=0 implies H≡1 mod4. Converse pairs all ordinary classes with H≡3 mod4, proving the stated multiplicity parity, including for any converse-invariant restricted family. The merged-node tiling count 2H/|Aut(T)| uses the free action on ordered Hamiltonian paths, and is correctly restricted to a non-self-converse pair. This does not exclude H=21 or resolve the general spectrum.

6. The generic mirror seam formula is correct, including singleton X. The four patterns using the seam are disjoint and exhaustive: XY, YXY, XYX, YXYX. Their factors give (E+Q)². Without the seam, the baseline H² loses exactly S² paths. Rédei parity then forces E+Q≡S mod2. Equal-parity squares agree modulo4, giving the even mirrored target's H≡1 mod4.

## Independent controls

`independent.py` imports no producer or repository modules. It generates cycles by directed DFS, computes packing coefficients by a vertex-pivot polynomial recurrence, and counts paths by memoized endpoint recursion. Fixed packings are separately checked by feasible combinations and compared with the quotient formula.

-167 self-converse orientations: all fixed-reversal orientations for n≤5, plus40 deterministic samples at each of n6 and n7.
-All anti-involutions and automorphisms are enumerated for the n≤5 controls; explicit conjugators and invariance of F/I are verified.
-285 portal cases across every labeled tournament of orders1 through4: both port parity and the full mirror-seam formula are compared with direct path counting on the constructed tournament, up to order8.
-Normal and optimized-mode outputs are identical. Checks use explicit exceptions, not removable assertions.

The finite tests are controls. The proof supplies the all-order statements. Novelty beyond the sources and targeted repository lineage is not asserted. No repository or Git files were modified.
