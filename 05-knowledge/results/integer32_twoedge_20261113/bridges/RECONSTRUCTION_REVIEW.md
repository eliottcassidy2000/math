# Independent review of the exact wedge image and repair

**Review verdict: PASS in the stated scopes.** This is an audit, not a status
promotion. The producer's R1–R3 remain RESERVED until the root's decision.
Reviewed 2026-11-13 against the following frozen SHA256 pins:

| Producer artifact | SHA256 |
|---|---|
| `reconstruction/PROOF.md` | `f2ba72de4e691ff25f46a1849109a4db30fa0a6917c31a3363b0190b5998647c` |
| `reconstruction/INHERITED_AUDIT.md` | `089421c85674284d820f12b2b22f116aaee29a05775d0f139e5ff58754ed448a` |
| `reconstruction/check.py` | `2f6991f470496b2b802fe8a20760fd459562148678b9015be2fab1272cc8ab05` |
| Both producer receipts | `bdd9fbc57328a637c1153cde5bea2c7c6f29d650d5744b6eb1d758055dbfcab1` |

## Proof and type review

**R1 — PASS.** With outgoing incidence bits a_vu, each wedge is their pairwise
XOR at v. The star equations construct those local incidences modulo one bit
per vertex. Their edge mismatches `1+a_uv+a_vu` have zero sum on every triangle
exactly because the three beta values sum to one. Complete-graph integration
therefore supplies the vertex gauge correcting all edge mismatches at once.
The line-graph kernel is exactly the constants, so realization has exactly
two opposite orientations. This proves sufficiency as well as necessity;
the parity matrix may be redundant but has the claimed codimension.
The formal jet claim keeps its fixed constant/linear/disjoint-pair/repeated
coefficients and certifies only the quadratic truncation.

**R2 — PASS.** A changed arc subset F changes precisely the line-graph cut,
with weight `sum_v d_F(v)(n−1−d_F(v))`. The four displayed path families between
adjacent line-graph vertices have distinct edges, including their middle edges
and terminal edges. There are exactly `2+2(n−3)=2n−4` paths. Every nonempty
proper cut separates an adjacent pair, so each path crosses and yields the
lower bound. A single changed arc attains it. Global converse gives the
excluded zero word; arbitrary relabelling is not part of this metric.

**R3 — PASS.** One erroneous wedge affects at most one vote in each of those
edge-disjoint path families. At most n−3 errors leave a strict majority out
of2n−4 votes. Each path uses at most three observed bits, giving the stated
O(n^4) abstract local-XOR algorithm. The implementation additionally checks
realizability and the final radius, so an out-of-promise input cannot receive
a false radius certificate. The midpoint of two distance2n−4 codewords gives
the sharp ambiguity at radius n−2. Translation by a valid cut shifts every
vote by its endpoint difference, justifying reference-word error enumeration.

The two-slot bit is a local mixed-incidence test because a directed triple
can occupy one slot and the complementary Hamilton path the other. A zero
one-slot bit has a different meaning. The transitive order-four hostile and
the proof's slot qualifications prevent exporting the repair guarantee to
one-slot support. No unrestricted incidence-matrix realization is inferred.

**Inherited one-slot review — PASS relative to its stated dependency.** The
common-edge argument for distinct legal relocations is valid, including its
special order-three and order-four gap cases. Propagation recovers directed
positive support. Completion of dark arcs remains the inherited BJR-based
directed-support theorem, as the producer explicitly says. This audit does
not independently reprove that external theorem or turn the small census
into its all-order replacement. The module-run reversal and endpoint
uniformity arguments preserve full-tournament per-slot first marginals;
they do not preserve all induced endpoint marginals or all joint moments.

## Independent exact checks

[check_reconstruction_audit.py](check_reconstruction_audit.py) imports the
frozen producer only as the implementation under test. Its independent
oracles build directed incidence sets from each arc orientation, compute
cut weights from changed-arc degrees, test star/triangle equations directly,
and search the complete code image for every nearby word. No producer encoder
or decoder is used to construct those expected answers.

Both [audit receipts](reconstruction.normal.json) and
[optimized receipts](reconstruction.optimized.json) agree, with **68,665
always-active gates**:

- All converse states at orders3..6 give image sizes4,32,512,16384 and minimum
  distances2,4,6,8 by the independent directed-incidence/cut calculations.
- Every formal field at orders3 and4 is checked:8 and4096 fields. The parity
  image, exact decoder and exhaustive nearest-word oracle agree. Exactly4
  and416 fields lie within the guaranteed radius.
- Every one of466 reference errors of weight at most2 at order5 is compared
  against the complete512-word code. At order6, the first256 patterns ordered
  by weight then lexicographically are compared against all16384 words; this
  is explicitly a sample, not all errors or all formal fields.

Separately rerunning the full producer in normal and optimized modes recovered
its exact saved57,186-gate byte stream and pinned receipt hash. The bridge
checker also supplies an independent literal permutation/cut enumeration of
two-slot supports at all labelled tournaments through5. The producer's
one-slot/dark-arc census shares its pinned ordinary representative bank,
but rebuilds paths and directed supports rather than importing old algorithms.

Reproduce the independent audit from the repository root:

```
python3 -B 05-knowledge/results/integer32_twoedge_20261113/bridges/check_reconstruction_audit.py
python3 -B -O 05-knowledge/results/integer32_twoedge_20261113/bridges/check_reconstruction_audit.py
```

No correction to R1–R3 or the frozen checker was required. The inherited
external dependency and the finite enumeration bounds remain explicit.
