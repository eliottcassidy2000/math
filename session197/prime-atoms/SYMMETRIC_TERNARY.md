# A stronger amplitude barrier for a self-converse fixed clone root

Session 197. This is a separate restriction on a symmetry-preserving subclass of the mechanism in `MINIMUM_TERNARY.md`.

**Theorem.** Suppose a rooted tournament `(T,r)` has an anti-involution fixing `r`, and transitive cloning of `r` gives `H_m=A*3^m-1` for every positive `m`, or merely for all sufficiently large `m`. Then

    A = 2 modulo4,  and A>=78.

No attainment at 78 is claimed. The general asymmetric mechanism has the smaller sharp amplitude 60. This theorem does not cover an arbitrary root of an unmarked self-converse tournament unless an anti-involution fixes that root.

## 1. The symmetry residue

An anti-involution is an arc-reversing permutation whose square is the identity. It fixes at most one vertex, so the source here has odd order and `r` is its unique fixed vertex. Replacing `r` by `m` transitive clones extends the anti-involution by reversing the clone order. Every even `m` therefore gives an even-order self-converse tournament.

The companion, independently audited `session197-orbit-wildcard/EVEN_CONVERSE.md` proves that every even-order self-converse tournament has Hamilton count one modulo four. Briefly, the odd-cycle collection formula gives `H=1+2c modulo4`, where `c` is the number of directed odd cycles. A fixed-point-free anti-involution pairs all such cycles: an invariant odd cycle would have an invariant odd-sized vertex set. Thus `c` is even.

Putting `m=2` in the proposed clone identity gives `9A-1=1 modulo4`, hence `A=2 modulo4`. The eventual-identity case is identical, using any sufficiently large even `m`, or the exact run expansion from the preceding packet.

## 2. An all-order strong-core cutoff

The run and SCC arguments in `MINIMUM_TERNARY.md` apply without modification. An eventual pure ternary identity is an exact one, its run matching number is two, and every strong component outside the root component has count one because its product divides `gcd(3A-1,9A-1)=1`. Delete those singleton components.

The remaining root component is strong and invariant under the anti-involution. Its restriction still fixes `r`, so this component has **odd** order. For a hypothetical amplitude `A<=74`, its count is `3A-1<=221`.

The exact order-eleven branch of Busch's Theorem 1 is `H>=9*5^2=225`. All orders at least twelve are also excluded by `5^11>221^3` and the general bound `5^(n-1)<=H^3`. Consequently the odd strong root component must have order at most nine.

The source is [Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), DOI `10.37236/1141`. The sharper order-eleven branch, rather than the weaker uniform exponential estimate at that order, is needed for this cutoff.

## 3. Complete enumeration with the root retained

Fix a source anti-involution `tau`. The induced involution on Hamilton words, obtained by reversal followed by `tau`, has an odd nonempty fixed set. Choose one fixed word and relabel its successive vertices `0,1,...,n-1`. Then `tau(i)=n-1-i`, and its fixed root is exactly the middle vertex `r=(n-1)/2`.

This is a root-preserving normalization, not merely an unrooted isomorphism argument. All consecutive arcs are forward; each remaining pair-orbit under reflection has one independent orientation bit. Thus every possible marked source appears in the complete universe of `2^floor((n-1)^2/4)` normalized graphs at its order. No source can be missed because its fixed root is different from a chosen canonical representative's middle vertex.

The exact matching-rank inventories are:

| Order | Presentations | Rank 0 | Rank 1 | Rank 2 | Rank 3 | Rank 4 |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| 3 | 2 | 1 | 1 | 0 | 0 | 0 |
| 5 | 16 | 1 | 2 | 13 | 0 | 0 |
| 7 | 512 | 2 | 6 | 154 | 350 | 0 |
| 9 | 65536 | 8 | 32 | 2544 | 10016 | 52936 |

Only rank two can give a pure ternary clock, by the exact matching/run theorem proved in `MINIMUM_TERNARY.md`. Every rank-two source is tested by exact Hamilton counts after one, two and, when needed, three clones. None satisfies both

    H_2=3H_1+2,  H_3=9H_1+8.

Hence no possible source of odd order at most nine has the required clock. The strong-core cutoff excludes every amplitude at most 74 at arbitrary order. The residue `A=2 modulo4` then forces `A>=78`, as claimed.

## 4. Verification and scope

Run `verify_symmetric_ternary.py` normally and under `python -O`. It uses the exact Hamilton, augmenting-path matching and physical cloning routines in `verify_minimum_ternary.py`, and independently generates all 66067 normalized sources displayed above. All physical anti-involutions, complete rank histograms, even-clone congruences and surviving clone equations are checked by explicit exceptions. A separate discovery probe used a forward Hamilton recurrence and subset-state matching and obtained the same rank histograms and empty pure-clock set.

The finite bank is complete only through order nine. Its conversion to an all-order amplitude bound uses the proved SCC reduction and exact strong-order cutoff. Neither absence at larger amplitude nor the existence of a symmetric clock at the lower bound is asserted.
