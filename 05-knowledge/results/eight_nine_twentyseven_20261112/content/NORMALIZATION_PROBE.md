# A bounded same-count normalization probe

**FINITE-EXACT.** No context-content improvement occurs in the explicitly
bounded physical operation family below. This is neither a global
counterexample to normalization nor an exhaustive fixed-H tournament census.

## Sources and operation universe

Start with the named order-seven minimum witnesses in `receipt.json`
for H55,171,175,189, whose context contents are5,9,7,63 respectively.
For each source T and every a,b>=0 with a+b<=3, form the literal ordinal
padding `TT_a→T→TT_b`, of order7 through10. On every such padded base:

1. Reverse each single arc.
2. Reverse each unordered pair of distinct arcs.
3. For every Hamilton path P of the padded base, complement every arc
   except the spine arcs of P.
4. If the base has order at most9, add one vertex in each of its2^n
   possible incidence orientations.

All candidates are counted directly. In particular, native spine
complementation is **not assumed to preserve H**. No neutral output is
recursively searched, and the list is not quotiented by isomorphism.
Exact duplicate neutral labelled outputs are deduplicated only when
computing their full context content.

## Outcome and stopping gate

| initial H | candidate attempts | neutral attempts | distinct neutral labelled outputs | context content of every neutral output |
|---:|---:|---:|---:|---:|
| 55 | 9907 | 86 | 82 | 5 |
| 171 | 11067 | 32 | 28 | 9 |
| 175 | 11107 | 32 | 28 | 7 |
| 189 | 11247 | 32 | 28 | 63 |

There are43328 operation attempts in total. The5900 native-spine attempts
give **no** neutral output. The one-arc, two-arc, and one-vertex families
give42,60,80 neutral attempts, respectively; none changes the content.
The full per-source breakdown and content histograms are recorded in
`normalization_receipt.json`.

The stopped implication is precise: the successful native H513 repair
cannot be generalized into a universal neutral operation, and neither
padding followed by one or two arc reversals nor one further vertex
supplies a repair for these four sources in this range. Larger changes,
different starting isomorphism classes, and sequences passing through other
counts remain untested. The conjectured unrestricted minimum
`min_(H(T)=N) g(T)=gcd(N,3)` remains **OPEN**.

## Controls and reproduction

The driver first reproduces the actual positive control C513→N513:
complement outside the spine `(1,0,4,5,6,2,7,3)` and check the reduction
27→3. An independent Python cover computation recounts this resulting
graph and verifies its literal transformation. Any improvement would
receive the same independent recount; there were none.

The finite driver is embedded in `normalization_probe.py`. It compiles
inside a temporary directory under this owned packet, removes its binary
and temporary source on completion, and has a180-second execution limit.
Run from the repository root:

    python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/content/normalization_probe.py

The script outputs both the explicit search receipt and the positive
control's complete ordered-cover profile. No compiled artifact is retained.
