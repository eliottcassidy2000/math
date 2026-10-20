# Blue extrema after the hollow-frame repair

Status: **PROVED** elementary mean, bounds, and explicit family;
**CITED** general regular-subset reduction; **FINITE-EXACT** maxima through
order eleven. No general exact maximum or asymptotic equivalent is claimed.

## 1. Inheritance, objects, and the corrected question

Use the hollow reflection frame of
[`ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md`](../odd_blue_spectrum_hollow_cubes_20261020/ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md),
Sections 2–4. Its pair index set has size `k`, its tournament has order
`n=2k+1`, and its external coordinates are:

- a tournament `A` on the `x_i` that beat the fixed center `f`;
- a symmetric graph with `b_ij=1` precisely when `x_i -> y_j`;
- an independent matching vector `a`, with `a_i=1` precisely when
  `x_i -> y_i`.

The external frame has `k(k-1)` independent bits: `binom(k,2)` tournament
bits and `binom(k,2)` cross-graph bits. Write

    h(D) = number of rho-symmetric directed Hamiltonian words,
    H = Aut(A,b),    H_a = stabilizer of the matching vector,
    B(T_a) = h(D)/|H_a|.

The path set depends on the external frame, not on the matching bits.
It is positive and odd. The finite group `H` embeds faithfully into the
odd-order tournament automorphism group `Aut(A)`.

The closest proved mechanism is precisely that centralizer quotient.
The hostile example is the inherited order-nine `B=39` witness, which
refutes the earlier `3^floor((n-2)/2)` cap even at odd order.
The corrected near miss is replacing `B` by `h` without retaining `H_a`.
The useful sidecar is the action on binary matching completions.

The bounded concept board for this continuation is:

1. External-frame cube: uniform labelled moments of `h`.
2. Matching cube: attainability of `B=h`, not an assumed generic rigidity.
3. Path blocks: elementary universal upper bounds.
4. Cyclic three-block tournaments: explicit superexponential examples.
5. Exact small frames: independent complete censuses, not sampled maxima.

The honest extremal object is

    M_k = max {B(T): T self-converse of order 2k+1}.

## 2. Every frame attains its full symmetric-path count

**CITED input.** David Gluck, *Trivial Set-Stabilizers in Finite
Permutation Groups*, Canadian Journal of Mathematics 35 (1983), 59–67,
Corollary 1, p. 66, says that every finite permutation group of odd order
has a subset with trivial setwise stabilizer. The property is defined on
p. 60. This includes intransitive actions; faithfulness matters.

Primary sources:

- DOI: <https://doi.org/10.4153/CJM-1983-005-2>
- Publisher PDF:
  <https://www.cambridge.org/core/services/aop-cambridge-core/content/view/9667B85D9DDD0FE911C8382287B86167/S0008414X00022550a.pdf/div-class-title-trivial-set-stabilizers-in-finite-permutation-groups-div.pdf>

The original theorem and definition were checked, not inferred from a
secondary abstract. Its proof uses solvability of odd-order groups; this
is not advertised as an elementary substitute for Feit–Thompson.

Apply it to the faithful permutation group `H` on the pair indices.
The characteristic vector of the regular subset is an admissible matching
completion with `H_a=1`. Since all other stabilizers have order at least one,

    max_a B(T_a) = h(D),             M_k = max_D h(D).       (2.1)

Thus the old maximum-blue question is exactly an external-frame path-count
question. The group quotient has not been silently dropped: a cited theorem
supplies the missing regular completion.

The finite maxima below do not require this citation: `B<=h` gives the upper
bound, and explicit trivial-stabilizer completions attain each reported
maximum.

## 3. Exact labelled mean and factorial lower bound

**PROVED.** If the `k(k-1)` external bits are uniform and independent, then

    E_D h(D) = k!,
    sum_D h(D) = 2^(k(k-1)) k!.                              (3.1)

A symmetric Hamiltonian word is determined by its first half

    z_(p_1), ..., z_(p_k), f,

where `p` is a permutation of the pair indices, one vertex is chosen from
each pair, and the last chosen vertex is an `x`. There are
`2^(k-1) k!` potential prefixes.

Each of the first `k-1` successive edges constrains exactly one external
bit: an `A` bit if both endpoints have the same polarity, and a `b` bit
if they have opposite polarity. The unordered index pairs are all
distinct, since they are edges of a simple index path. Consequently no
external edge orbit is reused, even if some edge types agree. These
`k-1` constraints are independent, with success probability `2^-(k-1)`.
The last edge to `f` is already forced. Summing their indicators proves
(3.1).

Together with (2.1) and oddness of `h`, this gives

    M_1=1,          M_k >= k!+1 for every k>=2.              (3.2)

In particular, no fixed-base exponential bound `C^k` can bound all `M_k`.
The historical `3^(k-1)` cap already fails at `k=4`; for every `k>=5`,
the factorial mean alone forces its failure because `k!>3^(k-1)`.
The mean is over labelled external frames, not over tournament classes
and not over matching-completion orbits.

## 4. Elementary upper bounds from disjoint path blocks

**PROVED.** For every `k>=1`,

    M_k <= 2^((k-1)/2) k!                 if k is odd,
    M_k <= 3 * 2^(k/2-2) k!               if k is even.     (4.1)

These inequalities also hold for every `h(D)`, without the cited input.

For two fixed pair indices, there are eight choices of their order and
polarities. Exactly four give a directed two-vertex word: each of the
four polarity choices selects two tournament vertices and exactly one
of the two orders is directed.

If the last vertex must be an `x`, there are four potential choices.
The two all-`x` orders contribute exactly one path. The two mixed choices
both end at an `x`; they contribute either zero or two, according to the
cross bit. Hence there are at most three directed choices.

For `k=2r+1`, divide the index word into `r` adjacent two-index blocks and
a final singleton. There are `k!/2^r` ordered choices of the unordered
blocks. Ignoring the edges between blocks leaves at most `4^r` choices
inside them; the last singleton is forced to be an `x`. This yields
`2^r k!`.

For `k=2r`, take `r` adjacent blocks. The first `r-1` have at most four
choices each and the last has at most three. This yields
`(k!/2^r)4^(r-1)3`, as stated. The `k=2` bound is sharp. No sharpness
claim is made at larger `k`.

## 5. An explicit elementary superexponential family

The growth conclusion is not only a nonconstructive consequence of
Gluck plus averaging.

Take `k=3m` and let

    A = C3[TT_m,TT_m,TT_m],

the cyclic tournament on three blocks, each inducing a transitive
`m`-vertex tournament. Put `b_ij=1` for all distinct indices.

Now no `y_i` can be followed by an `x_j` of a different pair. A valid prefix
ends at an `x`, so it contains no `y` at all. Therefore

    h(D) = HP(A).

There are at least `3(m!)^3` Hamiltonian paths in `A`: visit the three
blocks in cyclic order `m` times, choose one of three starting blocks,
and choose any permutation of the vertices in each block.

Here `Aut(A)=C3`. To verify that assertion, the vertices at each internal
rank have a common outdegree `m+r` and form a directed triple across the
three blocks. An automorphism must act on each rank triple by a cyclic
rotation. Between any two ranks the adjacency matrix is a cyclic
permutation matrix or its sum with the identity, so the rotations at the
two ranks must agree. Thus only the three simultaneous block rotations
remain. The argument also covers `m=1` directly.

Consequently the frame group is `C3`, and a matching vector supported at
one index has trivial stabilizer. Its completion satisfies

    B(T) = h(D) >= 3(m!)^3,       |T|=6m+1.                 (5.1)

This is a fully specified elementary family. For `m=5`, its certified
lower bound `5,184,000` exceeds `3^14=4,782,969`; the inequality continues
for every `m>=5`. Its growth exceeds every fixed-base exponential.
No exact count for this entire family is claimed.

## 6. Complete external-frame census through k=5

Two independent implementations agree at every labelled frame:

1. `enumerate_frames.cpp` sums each potential prefix's indicator over
   the subcube of frames that support it.
2. `replay_subset_dp.cpp` rebuilds every frame's adjacency and computes
   the count by a forward pair-subset dynamic program.

Both enumerate precisely all `2^(k(k-1))` external frames, not the
`2^(k^2)` matching-completed frames and not all tournaments.

| k | n | external frames | max h=max B | E h | Var h |
|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 1 | 1 | 1 | 0 |
| 2 | 5 | 4 | 3 | 2 | 1 |
| 3 | 7 | 64 | 9 | 6 | 6 |
| 4 | 9 | 4096 | 39 | 24 | 75 |
| 5 | 11 | 1048576 | 195 | 120 | 1500 |

The variances in this table are **FINITE-EXACT** moments, not an asserted
all-`k` formula.

The first implementation takes `125,829,120` subcube updates at `k=5`.
The independent replay compares all `1,048,576` outputs pointwise. On
the recorded run they took about `0.12` and `1.03` seconds, respectively.
No order-thirteen exhaustive computation was attempted.

`verify_extrema.py` also checks the first four counts against every
class in the inherited census: `2+8+88+2752=2850` class checks.
All first and second moments, histograms, maximum counts, and oddness
checks agree. These are complete finite maxima, not sample records.

### Explicit order-eleven witness

With vertices `x_i=2i`, `y_i=2i+1`, and `f=10`, the witness has external
mask `7208` and matching mask `2`. External bits are ordered by lexicographic
pairs `i<j`, first the `A` bits and then the `b` bits.

Its adjacency row masks are

    [1448,425,1099,880,1031,980,1043,860,1108,343,682].

Bit `j` in row `i` means `i -> j`. The anti-involution swaps consecutive
vertices and fixes `10`.

The frame group is generated by the pair-index cycle `(1 2 3)` and has
order three. The matching vector has its only `1` at index `1`, so its
stabilizer is trivial. Directly enumerating all `1920` potential signed
prefixes finds exactly `195` full symmetric Hamiltonian words. Thus

    (h, |H|, |H_a|, B, c3) = (195,3,1,195,52).

`n11_witness.json` retains the complete adjacency, all 195 words, the frame
automorphisms, and the matching-stabilizer certificate.

### Unique maximizing hollow type and eight maximum-blue classes

Exactly forty labelled external frames have `h=195`. They are precisely
the full `S5` orbit of the witness, of size `5!/3=40`. This is checked
by transforming the witness under every pair-index permutation and
comparing the resulting set to the exhaustive maximum list.

Within this unique maximizing hollow type, the group cycles indices
`1,2,3` and fixes `0,4`. Of its 32 matching vectors, eight are fixed by
the group and 24 have trivial stabilizer. After quotienting, there are
eight classes with `B=65` and eight classes with `B=195`.

Therefore the global maximum `195` occurs in exactly **eight**
self-converse tournament classes of order eleven, four on each triangle
parity shore. They form one auxiliary `Q3` packet under the inherited
odd-orbit toggle theorem. This statement concerns the auxiliary class
cube, not adjacency in the actual antipodal blue graph.

## 7. Reproduction and remaining frontier

Run from this packet's directory:

    g++ -O3 -std=c++17 enumerate_frames.cpp -o /tmp/enumerate_frames
    g++ -O3 -std=c++17 replay_subset_dp.cpp -o /tmp/replay_subset_dp
    for k in 1 2 3 4 5; do
      /tmp/enumerate_frames "$k" "h_k${k}.bin" > "exact_k${k}.json"
      /tmp/replay_subset_dp "$k" "h_k${k}.bin" > "replay_k${k}.json"
    done
    python3 verify_extrema.py
    python3 build_witness.py

All core checks run in an isolated copy of this packet. The inherited
class check is optional: `verify_extrema.py` automatically uses the sibling
`odd_blue_spectrum_hollow_cubes_20261020` when present, or the sibling
staging directory `session14-odd-blue-pairing` when present. Otherwise it
reports zero inherited-class checks and still completes all core checks.
An explicit inherited location can be supplied with
`--inherited /path/to/odd_blue_spectrum_hollow_cubes_20261020`. That directory
must contain `pairing_n3.json`, `pairing_n5.json`, `pairing_n7.json`, and
`pairing_n9.json`.

Binary arrays use native little-endian unsigned 16-bit integers on the
recorded platform; JSON summaries and source are portable. No nonstandard
packages are required. The command intentionally does not use `-DNDEBUG`;
the independent comparisons are assertions.

The main open questions are now:

- Determine the exact maximum of `h(D)` for general `k`.
- Determine the asymptotic growth of `M_k/k!`.
- Characterize maximizing external-frame types and their matching-orbit
  distributions.

The exact labelled mean does not determine the maximum, and the finite
variance sequence does not supply an asymptotic law. No literature-priority
or novelty claim is made. These results do not alter the still-open
all-order antipodal incidence problem or any LRC claim.
