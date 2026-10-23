# Delayed opposite blocks become temporary pulses

October 22, 2026. **PROVED:** the source-disjoint delayed-pulse replacement
at every finite native-depth cut. **FINITE-EXACT:** selected native prefixes
and improved all-future bounds for the same three frozen annuli. Uniform
153N for every N>=4096 and arbitrary-source frontier entry remain **OPEN**.

## 1. Inheritance and the new coordinate

The closest mechanism is synchronous two-channel cancellation in
[the linked tail envelope](../future_tail_blocks_20261021b/README.md).
The canonical hostile is two opposite vectors at different current nodes:
their sum being zero does not authorize erasure before both arrive.
The corrected near miss is confusing a primitive-edge delay with the
integer clock determining its sign. The least-used coordinate here is an
**authenticated delay between two current nodes**, rather than equality of
those nodes.

Board: immutable source fibres / current native node / two charge channels /
native delay / actual integer clock / first-arrival terminal. The anchor is
153N coverage, the niche is finite delay cancellation, and the wildcard is
the repeated affine-word connection to powers of three.

The setup, Q, delta, kappa, stopped frontier V and complete frozen prefixes
are exactly those of the linked envelope and
[angular packet](../tail_geometry_20261021c/README.md). Each pending current
node owns a disjoint set of original sources with their mod15/mod35 registers
and complete integer-clock prefix. Its signed four-vector is b=(w,h), with
old and new channels separate. Nothing below changes those sources.

## 2. The pulse theorem

Fix a native snapshot depth D. Suppose distinct pending nodes u,v have a
supplied literal native prefix

    P^ell(u)=v,  ell>0,

whose entire node list avoids V, and whose actual summed integer clock is a.
Assume the separate-channel identity

    b_v + (-1)^a b_u = 0.                              (1)

If v first reaches V at terminal t after tau further native edges and clock
c, then u first reaches that same terminal after tau+ell edges and clock
c+a. If v never enters V, neither does u. Therefore their combined four-vector
at any future cut D+r is exactly

    0                         r<tau,
    (-1)^c b_v                tau<=r<tau+ell,
    0                         r>=tau+ell.              (2)

For tau=infinity the vector is identically zero. No value of tau, c or t is
assumed, and the unknown suffix is never queried. Equation(2) follows by
adding the two first-arrival indicators and using(1). It is an identity in
all four coordinates before taking the indefinite energy delta.

Thus replace the pair, for purposes of an all-future upper bound, by one
**optional pulse** of vector b_v. At any fixed cut this pulse is either absent
or contributes one signed vector to one literal terminal. The generic linked
assignment bounds therefore apply. The replacement loses the restriction
that presence lasts precisely ell cuts; losing it enlarges the choice space
and is safe for an upper bound.

For several such pairs, require disjoint pending owners. Their native paths
may intersect, but their original source fibres may not be spent twice.
Summing(2) preserves the exact total old and new field at every cut, including
all cross terms with other pulses, ordinary pending blocks and arrived ports.
One must evaluate a complete bound on this replacement carrier; its saving
cannot simply be added to an independently optimized saving on the old carrier.

This is **not** a synchronous native state and is not numerical deletion of
both sources. In particular, the monotonicity theorem for ordinary prefix
extension does not automatically apply to a later update of this relaxation.
Keep the original prefix, pair paths and pulse types; a sound earlier gate
can always be retained in the running minimum.

## 3. Exact native improvement

Discovery explored at most16 native edges from each supplied pending node,
stopping at V. It greedily selected disjoint opposite pairs. This is not
an optimal matching claim. The consumer reads only the selected literal
paths, independently checks each native edge and its stripped 4u+1 tower,
checks avoidance of V, and verifies the actual summed clock and(1).

| N | D | Ordinary blocks before | Pulses | Blocks after relaxation | Selected prefix edges | Old angular bound | Pulse angular bound | 153N |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4096 | 105 | 195 | 39 | 156 | 158 | 551844 | 521233 | 626688 |
| 8192 | 110 | 348 | 63 | 285 | 237 | 1118881 | 1039528 | 1253376 |
| 32768 | 130 | 760 | 144 | 616 | 573 | 4142214 | 3856740 | 5013504 |

Each row retains the previous exact checks at all cuts through D. The pulse
bound covers every cut after D, including nonarrival. The original source,
prefix-edge and state commitments are unchanged. New selected paths add
authenticated child relations, not a claim that their endpoints are grounded
in V. Source-address decoding still recovers every mod15/mod35 guard and clock.

One literal new relation is1865 ->1399: native delay1, integer clock2.
At N8192 the two original pending vectors are opposite, so their future
sum is a pulse lasting one native cut. Using delay1 as its sign clock would
give the wrong sign and is explicitly rejected by the checker.

N7 remains a hostile: no selected pair, bound1153>1071, and an earlier cut
already fails. The intended all-N target starts at4096, not at1.

[Weighted destination coupling](../weighted_tail_20261022/README.md) can be
applied to the typed pulse relaxation itself. Its independent follow-on
certificate, when present, pins this packet's receipt and recomputes the
whole objective rather than adding improvements.

## 4. Why 3^k+1 suggests this operation

The relevant prior source is
[the 2,10,42 neighbour macros, section3](../collatz_2_10_42_neighbour_macros_20261005.md).
For p_k=(2*4^k+1)/3 and k>=2 it proves the exact odd-step word

    p_k --(1, 2 repeated k-1 times, 2+epsilon_k)--> y_k,
    y_k=(3^k+1)/2^epsilon_k=U(3^(k-1)),
    epsilon_k=1 for even k, 2 for odd k.

The source and smaller child have the same endpoint, but their odd-step
depths differ by k and their summed halving clocks differ by2k+1. This is
a genuine delayed, opposite-sign common future at the appropriate Nyquist
specialization. It is not immediate cancellation: their arrival cuts differ.
Nor does the endpoint identity alone prove that the child reaches1.

For comparison the formal repeated shortcut (3x+1)/2 sends2^k-1 to3^k-1
in k legal odd steps, whereas (3x-1)/2 sends2^k+1 to3^k+1. These identities
keep the affine intercept and clock; they do not identify the positive and
negative Collatz graphs. The native P certificate above is checked separately
and does not borrow an odd-step delay from this analogy.

## 5. Remaining uniform obligation and reproduction

A usable sufficient target is: for every N>=4096, find a finite D, check its
earlier cuts, authenticate a disjoint pulse bank, and make a sound all-future
gate on that relaxation at most153N. The present theorem supplies a new
finite operation; it does not prove that such a bank exists at every scale.
Repeatedly finding more pairs in finite cases is not that uniform theorem.

Run from the repository root:

    python -B 05-knowledge/results/lagged_tail_pulses_20261022/check.py
    python -B -O 05-knowledge/results/lagged_tail_pulses_20261022/check.py

`paths.json` contains the selected paths. `relaxations.json` explicitly marks
the output as optional-pulse relaxation, not a synchronous state. `receipt.json`
pins the complete old carrier, new paths and exact bounds. Controls include
arbitrary unknown finite tails and nonarrival, every cut in the declared tiny
model, independent edge arithmetic, duplicate owners, and corrupted clocks.

Receipt LF SHA256:
`c0be85b63ccb3c00d7263b7180c63a88f3c4fcc6d5c4225728578a189c90a4ef`.
