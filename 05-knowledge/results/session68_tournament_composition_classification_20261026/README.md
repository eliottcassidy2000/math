# Session68 — classify the decoration, retain the action

**Status:** independently audited proofs for named composition families,
with separately scoped exact controls. The canonical structural inputs and
normal-basis theorem are **CITED**. There is no claim of a closed atlas of all
prime tournaments, all metagraphs, or a simple-group classification analogue.

The 12-vertex example is not an isolated anomaly. It begins a complete family
of cyclic-triple realization fibres with an exact automorphism spectrum.
The arithmetic bridge is literal for a cyclic quotient action, but fails
when that action is replaced by the larger hierarchical wreath group.

## The classification obtained

| Carrier or operation | Exact result | Proof and boundary |
|---|---|---|
| Canonical strong-module tree | Ordered linear children and child-coloured oriented prime quotients classify ordinary tournaments; coherent child models split the automorphism extension | [Classification §§1–3](session68-classification/CLASSIFICATION.md); inherited/cited mechanism, not an atlas of prime quotients |
| Iterated C3 tree, with A/A-converse terminal blocks | Every ordinary realization of one fixed C3 hypergraph at `n=4·3^d` has a palette representative; all orders `3^(3^d+e)`, `0<=e<=(3^d-1)/2`, occur | [Full fibre and polynomial](session68-classification/CLASSIFICATION.md); palette representatives do not exhaust labelled realizations |
| Modular-prime carousel R_m, m odd | Ordinary types are binary necklaces; least period d gives `Aut=C3^m ⋊ C_(m/d)`; all labelled lifts number `2^(2m+1)` | [Arithmetic](session68-arithmetic/REPORT.md), [ordinary-fibre completion](session68-root/FIBRE_COMPLETION.md); concurrent carousel classification credited below |
| Chosen normal-basis coordinates | Rotation is Frobenius, least word period is field degree, and the surviving cyclic stabilizer is the relative Galois group | [Field dictionary and hostiles](session68-arithmetic/REPORT.md); no identification of tournament composition with field multiplication, and not a classifier for the larger wreath action |
| Two aligned necklaces of periods d,e | There are `gcd(d,e)` relative-phase gluings, each with joint period `lcm(d,e)` | [CRT/phase proof](session68-arithmetic/REPORT.md); separately quotienting the two necklaces loses the phase |
| Literal Hamilton-path-complement line | Every uniform member of either regular family exits its fixed cyclic-triple count; exact flux is `score(first)-score(last)-1` | [Line boundary](session68-root/METAGRAPH_BOUNDARY.md); counts and ordinary classification alone do not reconstruct target incidence |

Here `A=1→C3`. The semidirect products retain the actual factor-permutation
action; an abstract group order is not a substitute for that action.

### The full ternary spectrum

Let `P_d(z)` count ordinary realization types by the quotient-stabilizer
exponent e. Then

    P_0(z)=2,
    P_(d+1)(z)=(P_d(z)^3-P_d(z^3))/3+z P_d(z^3).

Every coefficient from degree 0 through `(3^d-1)/2` is positive and even;
both extreme coefficients are 2. The two extreme groups are elementary
abelian and iterated cyclic wreath groups, respectively. Every path-cover
coefficient is constant across the fibre, while ordinary marked-Hamilton-path
fibre sizes vary inversely with the automorphism order, by total ratio
`3^((n-4)/8)`. At n=12 the common H is 462915 and the two ordinary fibre
sizes are 5715 and 17145. These are not automatically converse-merged degrees.

### Same size and triangle count, different classification carriers

Both following families have 36 vertices and **1929 cyclic triples**,
but their cyclic-triple hypergraphs are not isomorphic: the canonical maximal
strong blocks have sizes 12 and 4, respectively.

| Family | Ordinary types | Full labelled lifts | Automorphism orders and numbers of types |
|---|---:|---:|---|
| Depth-two C3 tree | 24 | `2^22` | `3^9,…,3^13`: `2,10,8,2,2` |
| Prime carousel R9 | 60 | `2^19` | `3^9,3^10,3^11`: `56,2,2` |

The tree's 512 standard colourings meet all 24 ordinary orbits, not all
`2^22` labelled lifts. The carousel's cyclic action admits the normal-basis
field model. Two explicit degree-nine field elements nevertheless have
hierarchical stabilizers 9 and 1: [independent field controls](session68-root/FIELD_CONTROLS.md).
This is the precise failure of transferring a cyclic classifier to a wreath action.

## Random recovery and provenance

Seed **2026102668** selected six paths without replacement from each frozen,
sorted stratum: 3976 theorem files, 3337 reflections, 5580 results, and
18224 source files. This is stratified sampling, not uniform sampling of
the entire repository. The [native dispositions](session68-arithmetic/sample_dispositions.md)
distinguish useful leads, historical errors and bounded non-use.

The frozen full source blobs and the independent `verify_random_recovery.py`
reproduce all 24 draws and the native first-140-line excerpt extraction.
No sampled historical program is executed. Source integrity is not a claim
that every historical statement was true or independently audited.

The original research pin is `fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8`.
Native authored files are preserved byte-for-byte; their input pins are not
publication-parent pins. The surrounding [Session67](../session67_metagraph_enumeration_information_20261026/README.md)
enumeration/correction work is preserved independently.

### Concurrent results explicitly credited

- `d761e4cf0ff92264b043586a6160b274dae0b29c`:
  [decorated classification, sharp 4p onset, Sylow towers and marked hinges](../tournament_composition_classification_20261026/README.md).
- `c9a8823436ed8c49a3effafd40353880a9909a76`:
  [complete carousel fibres, necklace/Möbius strata, converse pairing, modular-flip spectra and pointed growth](../cyclictriple_composition_20261026/README.md).
- `993de01394dea261c50d4bd9c992f8bfd1325ecc`:
  [concurrent normal-basis/Frobenius dictionary](../module_shattering_lines_carousel_finite_field_20261026.md).
  That commit was labelled **pre-audit**. This packet's dictionary is independently
  proof-audited and hostile-tested; it does not inherit all of that commit's
  module-shattering, F1 or line-census claims by association.

The present continuation contributes the all-depth ternary spectrum and
full ordinary-fibre normalization, independent field/hierarchy hostiles,
the aligned relative-phase bridge, and the independently verified literal-line
flux boundary. No global historical-priority claim is made. The modular-flip
operation and the literal path-complement operation remain different carriers.

## Verification and remaining boundary

[Independent audit](session68-audit/AUDIT.md) and its exact reviewed-input
manifest separate mathematical acceptance from source integrity and finite
replay. The independent checker imports no producer code: it uses complete
128-realization orbit action, W2 orbit BFS, direct order-12 Hamilton-word DFS,
and all 385 covers of a six-vertex module-reversal control. No formal-assistant
verification or all-order finite census is claimed.

From this packet directory:

    python3 -I -B verify_packet.py --replay
    python3 -I -B -O verify_packet.py --replay
    python3 -I -B verify_mutations.py

The envelope checks the exact file universe, sizes and SHA256 hashes, copies
the packet to a temporary location, and compares isolated normal/-O outputs
against frozen receipts without rewriting them. It needs only Python's
standard library, no network or repository checkout. The native optional
Git-based random replay can separately be run with an explicit repository
path containing the pinned source commit.

**OPEN:** a practical atlas of all modular-prime tournaments, and a closed
exact run/allocation/port grammar for arbitrary literal metagraph composition.
Neither these family classifications nor their number-theoretic dictionary
settles those obligations or the unrelated LRC(14) frontier.
