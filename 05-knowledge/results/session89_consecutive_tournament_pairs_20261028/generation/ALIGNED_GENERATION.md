# Exact generation of consecutive, aligned tournament pairs

**PROVED:** the two quotient decompositions below, their full stabilizers, and
completeness of the finite algorithms. **FINITE-EXACT:** all representative keys,
named rooted-child classes and complete stabilizer groups agree among three
implementations through core order `n=5`; an independent literal `2^(n²)`-input
orbit census verifies all those fields through `n=4`. No order-six pair bank is
materialized, and no asymptotic speedup is claimed.

## Object and output contract

An object is `(T,S,m)`: two ordinary tournaments on the **same** `n` core vertices,
plus an `n`-bit root-outneighborhood `m`. Its second consecutive tournament is
`U=S+root`. The root is distinguished and fixed; one permutation relabels `T`,
`S`, and `m` together. There is no converse merging or swapping of components.
The labelled universe is exactly `2^(2 binom(n,2)+n)=2^(n²)`.

The full returned dictionary is

```
(T_representative_code, S_aligned_code, mask)
    -> (named_S88_rooted_child_key, complete_pair_stabilizer)
```

`T_representative_code` follows the byte-pinned Session87 degree-order convention.
Within that fixed representative, `(S_aligned_code,mask)` is lexicographically
least under `Aut(T)`. The stabilizer is the **entire concrete permutation group**
on the core, not just its size; every permutation fixes the distinguished root
when extended to `U`. A rooted-child key points into the retained Session88
rooted dictionary and its complete certificates. Matrices are packed with the
Session87 edge convention; a mask bit means `root -> core vertex`.

This is a labelled binary-matrix encoding, **not** ordinary simultaneous
row/column permutation equivalence of binary matrices. See the encoding lane
for the distinct action. Multiplying unrelated unlabelled class counts discards
the shared alignment and does not implement this contract.

## Two exact decompositions — proof

Use the forward relabelling convention `g.A` here. (The code stores inverse,
destination-to-source orderings, so its direct action pulls back coordinates.)

Fix representatives `T,S` with full groups `G=Aut(T)`, `H=Aut(S)`. Every core
pair whose components have these classes can be normalized to `(T,g.S)` for
some `g in S_n`. Two alignments `g,g'` describe isomorphic normalized core pairs
iff `g'=a g b` for some `a in G`, `b in H`: an isomorphism must preserve `T`,
and the remaining difference is an automorphism of `S`. Thus alignments are
exactly `G \ S_n / H`. At an alignment `g`, the surviving symmetry is

```
K_g = Aut(T,g.S) = G intersection g H g^-1.
```

For a fixed aligned core pair, two masks extend to equivalent pairs iff they
lie in one `K_g`-orbit. The full stabilizer of the extended pair is
`(K_g)_m={k in K_g : k.m=m}`. Different core-double-coset blocks cannot merge,
because deleting the distinguished root preserves the shared core pair.
Therefore one representative per double coset followed by one representative
per `K_g` mask orbit is a full, duplicate-free transversal with exact groups.

Alternatively, first fix an ordinary class `T` and a **pointed** child `(U,r)`.
The common core is `U-r`. Let `H=Aut(U,r)` restricted to this core. All shared
alignments are exactly `G \ S_n / H`, with full pair stabilizer
`G intersection g H g^-1`. Here there is no remaining mask stage because its
attachment has already been incorporated in the pointed child. This is exact,
but it often repeats alignment work across many attachments of the same core.

The proof includes zero vertices: `S_0` and all groups are trivial and there is
one empty-core/root pair.

## Three different meanings of “consecutive pairs”

### Independent tournaments on an aligned shared core — benchmarked here

`T_n` and `U_(n+1)` are independent tournament structures; the old vertices are
identified, but `U` restricted to them need **not** equal `T`. This is the full
`2^(n²)` universe and requires the alignment-plus-residual-mask decomposition.

### Nested one-vertex extensions — a proper subproblem

If the requirement is `U-root=T`, there is no independent second core tournament
or relative alignment. Fix an ordinary `T` and enumerate root masks modulo
`Aut(T)`; this is exactly the named pointed-child bank already computed by
Session88. Its class counts for `n=0,...,5` are **1, 2, 4, 12, 48, 296**. The
checker independently extracts the nested subset `T_code=S_code` from each
full aligned-pair bank and verifies equality with that pointed-state count.
The labelled nested universe has only `2^(binom(n,2)+n)` elements.

### An overlapping stream — reuse the middle object, retain the relative port

Suppose an incoming pair ends at a middle tournament `M`, and an outgoing pair
starts at the same ordinary tournament `M`. Normalize both copies to one chosen
canonical `M` and retain the actual normalization transports. Cache:

- its representative, the complete concrete `Aut(M)`, and Session88 named
  pointed row/root transporters;
- the incoming pair's full stabilizer **embedded in `Aut(M)`**, denoted `A`;
- the outgoing pair's full stabilizer embedded in `Aut(M)`, denoted `B`;
- the actual relative map whenever shared labels are already specified.

If the relative map is not prescribed, compatible abstract gluings are
`A \ Aut(M) / B`; they are not all of `S_|M|` and not a product of independent
pair-class counts. Indeed every admissible identification of the two canonical
middle copies is an automorphism of `M`; changing either pair representative
changes that identification on the left by `A` or on the right by `B`.
A gluing represented by `g` has chain stabilizer `A intersection g B g^-1`.

The incoming group fixes its distinguished deletion root and preserves the
previous tournament; its permutations extend from the incoming shared core to
all of `M` before the canonical transport is applied. After gluing, propagate
the resulting **residual chain stabilizer**, not the larger unconditioned
`Aut(M)`. Compose the actual normalization/inclusion transports when moving to
the next middle object. This is an exact cache/gluing contract, not a newly
benchmarked arbitrary-length stream implementation.

A designated ordinary class and node identify a **pointed isomorphism type**,
not a unique labelled orientation or transporter. For the explicitly tested
seven-vertex tournament `i -> j` iff `(j-i) mod 7` belongs to `{1,2,4}`, the full
automorphism group has order 21 and the stabilizer of vertex zero has order 3.
There are `6!/3=240` labelled realizations with that root in a fixed position;
three distinct automorphisms still fix the root. The checker verifies these
groups directly over all `7!` permutations. A representative plus a coherent
transport coset is the missing sidecar; a class name and anonymous node alone
do not replace it.

## What is reused, and why it saves work

1. Run pinned Session88 block filling through order `n+1` once per setup, retaining
   its ordinary banks, named pointed states and complete root transporters.
2. Normalize each ordinary core once, transporting its **certified full group**
   into the chosen representative. No pair-size canonicalizer is invoked.
3. Enumerate the distinct labelled images `S_n/H` of each ordinary `S` once.
4. For each first class `T`, mark its `G`-orbits on those images. The fixed
   elements of `G` simultaneously supply the complete residual group `K`.
5. Compute mask representatives **and their stabilizers** once per identical
   concrete subgroup `K`. Reuse that recipe for every core alignment having it.
6. Emit the full bank and attach the named Session88 pointed-child key.

At `n=5` there are 8,784 aligned core pairs but only **five** distinct concrete
residual subgroups in this representative convention. Their entire reusable
mask recipes need 216 mask-action evaluations. The final 279,968 output records
still have to be emitted; that output cost is not hidden or removed.

This implementation uses exhaustive permutations to obtain right-coset images,
not a scalable permutation-group package. One could replace that step with
stabilizer-chain double-coset machinery without changing the proof, but that
speedup is **not implemented or measured here**.

## Matched baselines and evidence

All three methods return exactly the same full dictionary. The baseline already
fixes the first ordinary class, reuses its automorphism group, precomputes mask
actions, reuses edge images per second tournament, and rejects at the first
smaller joint image. It then extends **every** labelled second tournament by
every mask. The successful change is hoisting alignment classification out of
the mask loop and reusing complete subgroup recipes, not comparing against an
unnecessarily literal `2^(n²) * n!` canonicalization baseline.

Every timed trial independently rebuilds and pays for the same Session88 setup
through `n+1`, including the named rooted dictionary and labelled lookup. Full
output equality is checked after timing. Benchmarks alternate execution order;
setup and total timings are recorded separately. Timings are local Python
observations, not complexity bounds. Very small orders have overhead and no
consistent speedup.

| n | Pair classes | Full stabilizer-order histogram | Core double cosets | Residual recipes |
|---:|---:|---|---:|---:|
| 0 | 1 | `1:1` | 1 | 1 |
| 1 | 2 | `1:2` | 1 | 1 |
| 2 | 8 | `1:8` | 2 | 1 |
| 3 | 88 | `1:84, 3:4` | 12 | 2 |
| 4 | 2,752 | `1:2720, 3:32` | 176 | 3 |
| 5 | 279,968 | `1:279448, 3:512, 5:8` | 8,784 | 5 |

At `n=5`, baseline work is 393,216 joint candidates, 658,680 joint comparisons
and 24,576 core-edge actions. Core-first uses 12,432 core-edge actions, 216
mask actions, and 1,440 right-coset permutation visits. These are different
operation types, deliberately **not** combined into a misleading speed factor.
The exact enumeration independently satisfies
`sum_outputs n!/|Aut(pair)| = 2^(n²)`.

Seven alternating, fresh-setup-inclusive trials in `results_n5.json` gave:

| n | Reference median | Core-first median | Rooted-first median | Reference/core-first |
|---:|---:|---:|---:|---:|
| 4 | 0.007833 s | 0.005727 s | 0.017767 s | 1.37× |
| 5 | 0.546617 s | 0.263014 s | 1.644683 s | 2.08× |

These are illustrative same-machine measurements under the recorded full-output
contract, not a universal performance guarantee. The baseline is faster at
orders 0–2 in this run; the order-three difference is small.

The alternative rooted-child-first construction has 279,968 double cosets
already at the alignment stage. It is correct but slower in the measured
implementation. A designated child class and root do not, by themselves,
remove the need to enumerate shared alignments.

## Hostile controls and boundaries

- **Class multiplication:** at `n=2`, there is one ordinary core class and four
  pointed-child classes, but **eight** aligned pairs, not four.
- **Wrong mask group:** take a cyclic first triple (code 2) and transitive second
  triple (code 0). `|G|=3`, but the residual intersection has order one. Masks
  1 and 2 are `G`-equivalent and are nevertheless distinct aligned pairs.
  Quotienting masks by `G` instead of `K` loses valid outputs.
- **Anonymous subgroup order:** two order-three groups at `n=4` have different
  concrete vertex actions and different mask orbits. Recipes are cached by the
  complete subgroup, never its order alone. A conjugacy-based cache would need
  the missing conjugating transport.
- The literal independent audit explicitly tests all labelled inputs only
  through `n=4`. At `n=5`, all three full banks agree and satisfy exact mass;
  no full 33,554,432-input independent audit is claimed.
- The source uses explicit runtime checks, not `assert`; verification is also
  run under `python -O`. A count match alone does not pass the verifier.

## Independent self-converse lift cross-check

`check_lift.py` independently implements the encoding lane's lift on
`{0,p_i,q_i}`: internal `p_i -> q_i`, plus layer `T`, minus layer converse `T`,
cross-layer symmetric discrepancy graph `R=T XOR S`, and root incidences from
`m`. Its explicit anti-involution fixes zero and exchanges each `p_i,q_i`.

For **every** pair class through `n=3`, the checker verifies that anti-involution
entry by entry, canonicalizes the full lift, and checks that its image is a
bijection onto the self-converse subset of the **entire** Session88 ordinary
bank. Lift orders 1, 3, 5, 7 have self-converse counts **1, 2, 8, 88**; the
order-seven test filters all 456 ordinary classes rather than a handpicked
subfamily. Full-order canonicalization is used only for this audit.

For each lift, the embedded complete pair automorphism group is checked equal
to the centralizer of the designated anti-involution in the full lift group.
These groups need not be the same as the full lift group: at lift order seven,
one class has pair group order 3 and full lift group order 21. The checker also
counts every involutive anti-automorphism and verifies the conjugacy-orbit
identity `#anti-involutions = |Aut(lift)| / |Aut(pair)|` in all tested cases.
No order-nine full tournament bank is claimed. The general bijection/conjugacy
proof belongs to the encoding lane; these checks independently audit it.

## Reproduce

From a repository containing the unchanged Session87/88 source:

```sh
python3 -B 05-knowledge/results/<session89-packet>/generation/verify.py
python3 -B -O 05-knowledge/results/<session89-packet>/generation/verify.py
```

From staging (replace the source path if necessary):

```sh
python3 -B generation/verify.py --baseline-root ../session88-symmetry/verified-repo
python3 -B generation/check_pairs.py \
  --baseline-root ../session88-symmetry/verified-repo --max-n 5 --benchmark 7
```

`verify.py` reconstructs results from source and checks a pinned semantic hash;
it does not require the timing receipt. `results_n5.json` records the measured
seven-repeat run. Its timing fields are excluded from semantic hashes.
The loader verifies the exact two inherited implementation files and never
writes into their checkout. Provenance pin: Session88
`03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`.
