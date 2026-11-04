# Ambient module reversals: the complete order-ten congruence sieve

**FINITE-EXACT complete order10 universe**, combined with the **PROVED**
all-order ambient descent in [B3-B4](../algebra/PROOF.md). Every9733056
tournament isomorphism class of order10 was checked. Exactly12670 satisfy
S=E modulo2; **none** also satisfy C=C transpose modulo4. This is a
necessary-condition sieve for module-converse collisions, not a complete
order11 reconstruction census and not a sufficiency test for a collision.

## Consequence in arbitrary ambient quotients

If reversing a nontrivial homogeneous module M leaves the full named
two-port response unchanged, it leaves H,S,E,C unchanged. B3-B4 force

    S_u(M)=E_u(M) mod2,
    C_uv(M)=C_vu(M) mod4.

Odd modules of size at least3 violate the second condition by the inherited
odd-current theorem, at every size. The complete inherited endpoint census
has no parity-balanced class at orders2,4,6 (class counts1,4,56). At order8
its6880 classes contain precisely two parity-balanced cases, codes39887362
and72392840. Both fail the cover-skew condition, as independently checked in
the algebra packet; one witness has C_02-C_20=150=2 modulo4. The present
complete order10 sieve has no combined survivor.

Therefore **reversing any nontrivial module of size at most11 is visible
to two ports in every ambient tournament**, regardless of ambient order.
An invisible single-module converse would require an even module of size
at least12, not self-converse, satisfying both modular conditions. This
statement does not exclude arbitrary changed-arc patterns, nor simultaneous
interacting reversals beyond the proved disjoint odd-module result.

## Complete universe and exact implementation

The source is the inherited, independently authenticated complete
[order10 bank](../../consecutive_flags_order10_20261028/THEORY.md):191536
order9 parent representatives and9733056 child addresses in increasing
delta-varint order. Each address specifies its parent and all nine birth
incidences. Both compressed and decompressed files are SHA-256 checked
against that bank's manifest before computation; source hashes are frozen
before and after. All addresses, counts and full-input EOF are checked.
There is no strongness, connectivity, regularity or symmetry filter.

The first pass uses an exact bitset endpoint DP modulo2. For every support
X and v in X, the bit at v is the parity of the old end bits on X-v at
predecessors of v; a singleton contributes1. Thus, by induction on support,
it equals E_v(T[X]) modulo2. Running on the converse gives S. There is no
probabilistic or hash filter. Every one of the12670 passing classes is then
recomputed with the inherited integer endpoint DP and full ordered C matrix.
All endpoint bits are compared again; every ordered skew is tested modulo4.
Unsigned subtraction is valid here because2^32 is divisible by4, and all
integer endpoint/cover counts themselves fit in32 bits at order10.

For a class failing endpoint parity, B3 excludes it. For a passing class
failing any ordered skew, B4 excludes it. Isomorphism-class enumeration
suffices because both conditions are invariant under every relabelling;
no sidecar canonicalization or untested quotient representative is used.

The full receipt is [full10.json](full10.json). It reports an empty list
of combined survivors. An empty list is not silently promoted to a theorem
about unenumerated orders or about other mechanisms of response equality.

## Independent controls and reproduction

The independent integer audit compiles the **actual** fast parity function
with a separate exact-DP comparison harness. It checks all1099 labelled
tournaments of orders1..5 embedded into order10 by an ordered transitive
tail, and1024 specified pseudorandom order10 tournaments:2123 comparisons
of all ten start and end bits. The recurrence proof excludes false negatives
outside this finite control family. The first10000 complete-bank classes
also have normal/-O matching receipts. The full census is one compiled-O3
run; the bounded controls are not advertised as its independent repetition.

```text
python3 run.py
python3 run.py --limit 10000
python3 -O run.py --limit 10000
python3 audit.py
python3 -O audit.py
```

The full source and interpretation are reviewed in
[INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md). Lower even-order counts and
both order8 hostiles are pinned to the prior endpoint/parity receipts and
this packet's algebra checks. These results strengthen an all-order ambient
theorem along one specific operation; they do not answer universal two-port
reconstruction or imply unit-distance colour forcing.
