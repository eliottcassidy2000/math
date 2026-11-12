# Higher-order exact certificates by complementary family elimination

**Session283. PROVED algorithms and threshold grammar; FINITE-EXACT computation.**
A000568 counts are ordinary unlabelled tournaments. These results do not solve
Chan–Pak and do not partition tournament isomorphism classes themselves:
Burnside sectors are permutation cycle types whose fixed-tournament sets can
overlap. Their weighted aggregate counts classes exactly.

## What changed

Session279's strongest degree-four certificates needed a moment packet that its
producer obtained by listing **all odd partitions**. The new producer lists only
long residual cores and sums all short/backbone completions behind a finite gcd
response interface. An integer factorial basis makes that recursion practical.
It outputs exactly the same richer `(M,delta,R)` packets, retaining the data
that the earlier efficient unstratified jet omitted.

[Factorization and exact compiler proof](PROOF.md) · [Sharp periodic boundary and adaptive schedule](THRESHOLD_PROOF.md)

The first failure of the degree-four certificate is now classified for **every
initial odd backbone top b>=27**. It is

    N_b = 4(b+2)+t(b+2),

where t is an explicitly classified periodic function on odd integers, with
period6930 and values12,16,18,20,22,24,26,28. Only53 bounded four-cycle patterns
are needed. This yields a constant-width search for the smallest sufficient
backbone in this certificate family at every n, but **not** a constant-time
calculation of the moments or the huge exact count.

## Exact higher-order results

The consumer never sees a stored exact answer. When its interval is not a
singleton it combines that interval with the inherited independent modular
engine. Only then are its recovered integers compared with the inherited full
counter and cache. All agree.

| n | backbone top | residual cores listed | exact residue bits needed | packet seconds | inherited full-counter seconds |
|---:|---:|---:|---:|---:|---:|
|64|11|266|0|0.023|0.0027|
|100|19|608|58|0.122|0.0118|
|150|27|2210|92|1.242|0.126|
|150|35|774|0|1.168|0.145|
|200|49|1240|0|7.763|0.812|
|250|59|2539|0|42.546|5.656|

The n250 output has29,489 bits. These are not newly discovered numerical
values: the point is that the tighter structural certificate recovers them
without an independently supplied residue. **The full optimized C++ counter
is faster than this Python certificate pipeline in these matched measurements.**
The new producer is faster than the old all-partition **packet producer** at
n64/n70; that is a different workload and must not be advertised as a full-
counter speedup. Timings are machine-specific single observations, not a
statistical benchmark or asymptotic guarantee.

The exact adaptive schedule selects b33 atn150 (first rational ambiguity156),
b45 atn200 (210), and b57 atn250 (260). Fresh minimal-backbone runs at150/200
also return singleton intervals. The b57 selection at250 follows from the
proved schedule; the saved actual n250 packet run used the slightly larger59.

## Reproduction

From the repository root:

```bash
P=05-knowledge/results/session283_grammar_20261112/enumeration
python3 "$P/residual_compiler.py" > /tmp/packet-controls.json
python3 "$P/thresholds.py" > /tmp/thresholds.json
python3 "$P/complete.py" --cases 64:11,100:19,150:27,150:35,200:49,250:59 > /tmp/completion.json
python3 "$P/complete.py" --cases 150:33,200:45 > /tmp/minimal-completion.json
```

Requires Python3 and, for `complete.py`, g++ plus GMP. Run with `python3 -O`
as well; mathematical data agree after excluding timing fields. The compiler
controls compare integer-short, raw-rational, projective-rational and direct
full-partition packets at n0..24 for four backbones, then compare the full
packets against the old producer at64/70. The threshold receipt stores all
3465 odd residue classes compactly as indices into its53-pattern table.

Saved `packet_n150_B35.json` is a standalone packet accepted by the inherited
Session279 `ground_moments.certify` consumer. `completion.json` records the
recovered integers, residues, intervals, statistics and dependency hashes;
`minimal_completion.json` records the smaller scheduled backbones.

### Exploratory failure retained in scope

The first n200 runs computed their recursion but failed when Python's default
4300-decimal-digit string limit rejected packet serialization. They were not
reported as successful certificates. The final producer explicitly enables
arbitrary-size integer rendering and the full final cases were rerun. No
mathematical count or bound changed.

## Remaining limits

At fixed backbone the residual-core list still grows superpolynomially.
Increasing the backbone transfers complexity to its response recursion; the
finite periodic schedule controls the *certificate obstruction*, not that
state-space cost. Exact class counting and isomorphism classification remain
different tasks. The concrete Chan–Pak analogy is the need for a sufficient
continuation interface and a singleton intersection of independent constraints;
there is no theorem here that turns a class-count residue into a path-count
realizer.
