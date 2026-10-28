# Session88: targeted review of concurrent symmetry packets

2026-10-28. **Source-read + independently replayed small witnesses.**
This is an integration review, not a second order-nine research run. Session88's baseline pin remains
`7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51`.

## Exact sources and evidence levels

1. [Strong attachment/tail packet at 6a9503c](https://github.com/eliottcassidy2000/math/tree/6a9503cc7c1eb6adeae7b4e3b1db11f5d2f75105/05-knowledge/results/strong_attachment_tail_compiler_20261028):
   README, THEORY N1/N2/N4/E3, AUDIT, native.py/native-normal.json,
   root_switch_skip.py/root-switch-skip.json, benchmark-roots.json.
2. [Kernel125 packet at e7ac99a](https://github.com/eliottcassidy2000/math/tree/e7ac99aad4b4bf9d43348bc800c030928c3d98a0/05-knowledge/results/kernel125_closed_carriers_20261028):
   README, THEORY K4/K5, AUDIT, root_minimum.py,
   root_degree.json/root_neighbor.json.

The listed files were read at their exact pinned commits. Both source audits explicitly call their tests
**same-agent** alternate-path verification, not independent review. Their
reported full finite results remain source-reported here unless listed in the
independent replay section below. No literature-priority claim is made.

## Constructive results: what composes and what does not

| Producer / output obligation | Reported finite scope | Consequence |
|---|---|---|
| Strong-tail E3: bank, Aut-order and block-size histograms | All classes through 9; n=9 rows 1,716,608 → 191,536, classifications 15,449,472 → 1,723,824 | Same ordered 191,536-member bank; complete numeric rows are computed only once/class. |
| Session88 full blockfill | Through 8; root classifications 434,048 → 55,040 | Complete named Q, every indexed root transporter set, full Aut permutations, parent data and raw bank. |
| Kernel125 K5 degree / neighbor-degree selectors | Through 9; competitor queries 278,461 → 24,005 | One representative/class and Aut-order histograms, with exact bank bijection; not full Q or all root responses. |

**Implementation detail:** `root_switch_skip.generate` does not retain its
individual row supports in the returned result. It marks `assigned[]`, keeps
the children and aggregate histograms, and discards the row. The theorem allows
Q recovery *if* supports and one row/block are stored; the saved current output
requires reclassifying representatives or adding this storage to recover Q.
It also does not retain Session88's full vertex-indexed transport sidecar.

Strong-tail's matched order-eight 5.31× is 7.495934/1.411203 seconds, three
alternating samples, parent setup included and extra Hamilton/prime audits
omitted on both sides. Session88's 3.69× is 12.254597/3.323452 seconds for its
larger full-output obligation; its separate independent audits are also outside
timing. Neither ratio transfers to the other's output contract, and neither is
an asymptotic factor. Kernel125's 19.877/9.842 generation seconds are descriptive
runs before full-bank bijection, not a matched cross-session speed claim.

The selectors use a **minimum** invariant/key order; Session87 used maximum
signatures and a different pointed-key order. Counts and literal representative
identities cannot be copied between them. Their exact bank bijections, not
matching aggregate counts alone, certify the same ordinary classes.

Safe composition is: select one representative/class, compute its **complete**
root response (or infer every omitted root by certified automorphisms), then
fill its named equivalence block and transport the sidecar. A higher invariant
can dismiss a competitor for a minimum-root decision but cannot identify that
root's missing pointed type or transporters. All parent-mask orbits are still
generated. No multiplicative combination of timing/query ratios is justified.

## Mathematical scope and independent checks

### P7 pair-versus-joint marker hostile — independently replayed

On Z/7, orient by steps {1,2,4}. The ordered triples (0,1,3) and (0,1,5)
have identical pair-orbit data but distinct triple orbits, with 0 versus 1
common outneighbor. Marker j beats exactly its specified old port; every other
old vertex beats it, and the markers form the specified transitive order.
For order (0,1,2), the resulting uncoloured ten-vertex tournaments have
H=4281 versus 4315. All six marker orders separate (the other value pair is
4287/4325). Different H is an uncoloured invariant, so marker colours need not
be retained in the observation. This is not claimed to be modular substitution.

`check_small.py` independently exhausts all 5,040 S7 permutations, checks all
six ordered pair observations, and computes every marker H by both subset DP
and nonmemoized path DFS. No source module is imported.

### Free-frame tree transport — scoped proof; P7 finite replay

For one fixed free transitive G-orbit of frames, s sites modulo simultaneous
gauge have exactly |G|^(s−1) alignments: normalize the first chart to identity,
leaving s−1 arbitrary group elements. Tree gains recover all other relative
positions by ordered path multiplication. This is an all-s elementary proof
only under that fixed-free-orbit/query-language hypothesis; it does not count
arbitrary port choices, changes of atom, or complement/reset closure.

Strong-tail N1 uses **left gauge**, gains g_i^−1 g_j. Kernel125 K4 uses
**right gauge**, gains g_j g_i^−1. These are consistent conventions but require
opposite path product orders; Paley21 is nonabelian. The new checker obtains
Aut(P7)=21 literally, verifies a minimum free frame (0,1) of size 2, and checks
both conventions on all 21^3=9,261 four-site relative configurations and
148,176 ordered endpoint path checks (including identities). The other five
atoms' reported minimum frame sizes are source-read only.

### Native complement/reset closure — source-reported finite evidence

Strong-tail N4 uses every path presentation X_n, of size
2^binomial(n−1,2), with consecutive path arcs fixed. For P7, H=189 and
|Aut|=21 give nine path orbits; these meet seven immediate partner classes.
The source exhausts native complement plus fresh unnamed reset at n=7 and
reports all **456 ordinary classes and 32,768 presentations** reachable.
The number 32,768 is path presentations at order 7, not all labelled order-7
tournaments. Native complement is not global converse. Strongness is not an
imposed closure filter. The full closure was not independently replayed here.

### Bounded observational separation — source recurrence plus small replay

Observations are ordinary isomorphism class and blue fixed-point flag. The
menu is native complement, fresh unnamed Hamilton-path reset, and ordered
singleton joins on either side, with finite horizon M. Downward signatures
record the observation and both successor classes, paired with the signature
of the complement. Induction gives the coarsest carrier for this typed finite
menu; class-constant reset rows introduce no extra split. This is not the
full named-language theorem and not a universal joint-attachment minimum.

Source horizon 7 has 1,1,2,8,64,1019,24754 cells at orders 1..7; precisely five
two-element cells remain at order 6. The source supplies order-eight contexts
separating all five and therefore obtains a discrete horizon-8 carrier for
n≤6. The new checker independently reconstructs the five stated terminal mask
pairs and exhausts all 40,320 S8 permutations for each, finding no isomorphism.
It does **not** independently establish that these are the only residual pairs:
the exhaustive horizon-7 partition remains source-reported. All-order
faithfulness for this exact singleton/unnamed menu remains OPEN.

## Reproduction and stopping boundary

```sh
python3 -B check_small.py > /tmp/concurrent-normal.json
python3 -B -O check_small.py > /tmp/concurrent-optimized.json
cmp /tmp/concurrent-normal.json /tmp/concurrent-optimized.json
```

Both runs PASS with byte-identical receipts:
`8189ae9ab16cc7498ef32863ee734404b523b3285a1b71ef1c62b99b136c366a`.
All guards use explicit exceptions, not removable assertions. The two saved
receipts record the finite universe and omitted replays. No order-nine run,
full n7 closure rerun, new all-order classification, or new performance claim
was needed for this integration review.
