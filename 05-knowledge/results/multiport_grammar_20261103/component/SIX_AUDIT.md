# Independent audit: strong-six supplier, named rewrite and complete scalar tree

2026-11-03. **INDEPENDENT ANALYTIC / SOURCE / REPLAY AUDIT PASS.**
This audit covers [SIX_PROBE.md](SIX_PROBE.md), its complete six-atom and
two-target receipt, and the supplementary [SIX_TREE.md](SIX_TREE.md) with all
28 successor decisions. All producer artifacts were preserved. No correction
remains. Native compilation and its finite-clock patches are separate from
this certificate.

## 1. Exact atom coverage, independently reconstructed

I read the complete producer source. Its augmentation argument generates all
isomorphism classes: deleting any vertex reduces to a previously generated
class, and its transported incidence vector is among all tested augmentations.
The canonicalization examines every permutation, before any strongness or
port filter. The orbit-stabilizer calculation is used with actual stabilizers.

The separate verifier
[`six_audit.py`](../../../../04-computation/multiport_grammar_20261103/six_audit.py)
imports **no producer implementation**. It independently enumerates all
labelled tournaments through order six and tests strongness by breadth-first
reachability from every vertex. For every proposed class representative it
generates its full labelled orbit, using the opposite relabelling convention
from the producer. It checks the canonical minimum, stabilizer size, pairwise
disjointness of the orbits, and exact equality of their union with the complete
labelled universe. The strong-orbit union is separately compared with the
independently enumerated strong graphs.

At order six this gives **32,768 labelled tournaments, 56 classes, 22,320
labelled strong tournaments and 35 strong classes**. In particular the 35
classes are not accepted merely from a published class count or aggregate
orbit mass. Their disjoint actual labelled orbits are verified to be the
entire independent strong universe.

For each strong representative and every port subset, the verifier counts
ordered two-slot covers by one full permutation and a cut. A cut removes only
the edge across the slot boundary; either empty slot is allowed. Endpoint
vectors and a joint endpoint matrix then give every port signature. This
differs both from the producer's complementary-subset permutation calculation
and from its inherited subset DP. All **2,240 raw marked order-six cases**
produce 1,950 distinct signatures. The combined bank has exactly 2,120
signatures, and every stored actual representative is independently recounted.
The complete sorted signature set has compact-JSON SHA256
`e2283995f5ec7b288f0883605035d4530adfc2e15f0d3929986fea69b3cf0b75`.

Relabelling transfers every port subset of a labelled graph to a subset of
its class representative. Therefore testing all 64 subsets of each strong
class is exhaustive for the whole labelled marked universe. Full-signature
replacement is sound for ordered assembly and its scalar observers, but does
not identify named vertices or native guard data.

The finite order-six statistics independently agree: minimum h=15, maximum
h=45, and minimum positive f=15. Neither requested target is supplied by one
strong-six atom. The two-component `(6,1)` supplier is consequently a precise
hostile to extrapolating this single-atom failure to arbitrary words.

## 2. Independent all-length, all-fixed-shift inverse

The audited bounded-intercept proof applies without modification to this bank.
For a normalized word the first and final atoms are marked. With at least
three atoms, the actual cross term `g_first*f_tail` gives `j>=2h`. With one
or two atoms, `h<=45^2`. Thus `h<=max(2025,|c|)`, and retaining that same joint
term bounds f by the displayed `L(c)`. Initial unmarked prefixes are quotiented
only by `(h_prefix,odd(p_prefix))`, with one actual representative retained.
The lost dyadic factor is absorbed by a fixed clock shift; other observer
coordinates are not identified.

The source's grouping by `(H,P,F)` does not quotient the remaining coordinates:
it expands every actual representative after each divisibility calculation.
Its recursion terminates because the positive tail marginal decreases by at
least a factor of two. Normalization removes words whose first atom is
unmarked only after the complete recursion. The right-inverse controls retain
zero-marginal prefixes, as required.

The independent verifier uses counters of **full signatures**, retaining their
word multiplicities rather than the producer's word-index lists. It builds its
own prefix monoid and Hamilton-count monoid from the independently recovered
bank. It deliberately uses the strictly larger envelope

```
B_audit(c)=2*45^2+|c|
```

and the corresponding larger f bound. For the two original targets it makes
33 normalized fixed-marginal queries rather than the producer's 30, reproduces
every producer query's candidate multiplicity inside its tighter bound, and
finds exactly the same accepted full-signature multiplicities:

| Address | Producer queries | Producer normalized visits | Independent result |
|---|---:|---:|---|
|24|20|1,610|none|
|44|10|805|one representative word/signature|

The prefix pairs are independently exactly `(1,1),(5,1)` for 24 and `(1,1)`
for 44. In particular the prefix of Hamilton count five is not omitted from
the negative 24 decision. These are complete all-word decisions, not bounded
word-length scans. Uniqueness concerns the chosen signature-representative
bank and prefix quotient; it is not uniqueness of labelled physical suffixes.

## 3. Actual positive and old-name-preserving adjacency path

Literal permutation/cut counting confirms the strong-six atom

```
rows=(2,4,57,3,11,27), ports={4},
sigma=(17,224,28,44,0,2,4).
```

Appending a marked singleton gives

```
rows=(66,68,121,67,75,91,0), ports={4,6},
sigma=(17,448,252,105,46,17,4).
```

Both cross terms in the joint law are retained: the new joint is `2+44=46`.
The resulting family satisfies

```
-29+252*2^(t-2) = 3*(-11+21*2^t)+4,  t>=3.
```

The physical refinement also passes. In the original 4-sibling, move the last
two prefix names into the suffix description, without changing a single arc.
The new suffix is literally TT7 with ports `{2,4,6}` and signature
`(1,128,84,21,12,1,0)`. The displayed ten internal reversals are exactly the
inversions of the proposed first-six tournament. The eleventh reversal removes
port 2 by reversing its old edge to the same physical root. All old names,
the root and total order `t+6` are retained.

Vertex 6 remains a marked terminal singleton throughout, including for any
ordering of these eleven independent edge reversals. The S182 terminal-SCC
criterion therefore proves strongness of every complete intermediate graph
for **all t>=3**. This is an analytic all-clock assertion; the literal checks
at t=3,4,5 are controls. A second independent endpoint-removal DP recounts all
36 complete intermediate graphs at those clocks and verifies their named
root defects. All twelve suffix signatures, all actual adjacency states and
all 36 root/count rows agree with the producer receipt. The root defect changes
from 12 to 46, with its intermediate changes explicitly retained.

The nonempty-prefix restriction matters. The same bulk target at the original
clock t=2 has an empty reblocked prefix and actual count **219**, whereas
`3*73+4=223`. This is included as an independent hostile control. The bulk
formula is not asserted at t=1 or t=2. Any separate finite-clock patches and
native helper programs require their own proofs; they cannot change this
fixed-suffix, shift-minus-two domain.

## 4. Complete supplementary fourteen-state scalar tree

I read the whole supplementary source and `SIX_TREE.md`. The BFS applies the
same exhaustive successor decision to both children of every accepted state.
Its depth cap is a stopping safeguard; an incomplete nonempty frontier would
be reported explicitly. Here the frontier is empty at depth six, before the
requested cap eight.

The independent `--tree` mode rechecks **all 28 decisions** with the signature
counter and larger envelope above. It makes **324 fixed-marginal queries**,
compared with 288 in the producer, reproduces every tighter-query candidate
multiplicity (43,939 total visits), and agrees with every accepted full-signature
multiplicity. All stored accepted witnesses are independently recounted by
endpoint DP. It then checks the initial state, both arithmetic successors of
every declared state, the exact accepted-state list, and the empty frontier.

The complete scalar address set is therefore

```
empty,2,4,22,42,44,422,442,444,4222,4442,4444,42222,44444.
```

Induction on address length pays the all-depth conclusion: there is no
unexamined accepted state from which a longer address could start. This is
not an inference from six successful sampled levels. Equality of `(c,odd(f))`
is used only for scalar existence under arbitrary replacement and fixed clock
shifts. It supplies no equality of actual sources or their changing defects.

Consequently the common-entry language with SCCs through order six is still
finite. The negative 24 decision requires an SCC of order at least seven for
a repair **within this uniform attachment architecture**, or a departure from
that architecture/uniform requirement. It neither proves seven sufficient nor
excludes other entry states, separately chosen integer witnesses, clock-varying
suffixes, or different constructions. No reusable two-branch physical class or
global spectrum theorem follows. Only the 4-to-44 edge has the named adjacency
path audited above; native programs for any edges are a separate interface.

## 5. Replays and frozen pins

Independent normal and optimized full `atom_six_probe.py --shift` runs pass
**2,999 always-active gates** and reproduce both JSON and stdout byte-for-byte.
An independent optimized supplementary tree run reproduces its frozen receipt.
The independent verifier's normal and optimized `--tree` runs both pass
**3,443 gates** and produce identical 16,624-byte receipts.

Commands, from the worktree root:

```
python3 04-computation/multiport_grammar_20261103/atom_six_probe.py --shift --output /tmp/atom_six_normal.json
python3 -O 04-computation/multiport_grammar_20261103/atom_six_probe.py --shift --output /tmp/atom_six_optimized.json
python3 -O 04-computation/multiport_grammar_20261103/atom_six_probe_tree.py --output /tmp/six_tree_audit.json
python3 04-computation/multiport_grammar_20261103/six_audit.py --tree
python3 -O 04-computation/multiport_grammar_20261103/six_audit.py --tree --output /tmp/six_full_audit_optimized.json
```

| Artifact | SHA256 |
|---|---|
| `SIX_PROBE.md`, pre-promotion | `39a21eb2b716d3ae79e1b858a1594ff4957c43f57bb7efa14749e0381aaf081a` |
| `atom_six_probe.py` | `ed295d2f7ce1a7815a7c377ede9544f3118efe3736458d1fc921d460c21900ba` |
| `atom_six_probe.json` | `1a57397d99deacf460d6e7c01f80da705b0304240d0fec363108cedb5bc4d951` |
| `SIX_TREE.md`, pre-promotion | `525c4aeef5795740a520a8a4e4345da86f1b5550c5d9e906127697eef21f5cbf` |
| `atom_six_probe_tree.py` | `cfb8e65193a9db580ab0613ec57f07439e3a236e2b2dd1b14b54a29621881d45` |
| `six_tree_probe.json` | `53e6eaf99e78acff793c24525338e25776c65f90a24ed147ca067caee87d0dc0` |
| `six_audit.py` | `9e4973431b1dd7c3df19dc57e8137bade858c4d1a2434e5bd0b03cae06beab7b` |
| `six_audit.json` | `1771c82c68794d3e74e1268c09c835773af66146b80b6263384dc188e682f34c` |

**PASS.** All mathematical and implementation claims covered above are accepted
with their stated scopes. The producer sources, previous audit bundles and
Git state were not modified. Parent owns status promotion and integration.
