# Minimal corrections for root integration — not applied by auditor

## 1. Historical THM-209 truth-surface banner

Target: `01-canon/theorems/THM-209-independence-polynomial.md`.

Pinned commit: `63f52936f8dd8592d4f2fee869748baad2e1997e`.
Git blob: `983d8655b66627fbd39661c69a3a324f2fe261b6`.
SHA256: `b6bf2c296da08e6f16bee6a1e04a48bf0a38ee3b48c57054298172e3fc1bfd2f`.

The same bytes were independently checked in local HEAD
`97f0fde6f26d7ad732442b78593b598eb5484df4` and the worktree. This is a
local HEAD check, not a new remote fetch or claim about subsequently arriving
commits. Provenance: `THM209-correction-provenance.json`.

Original status line:

```
**Status:** VERIFIED (exhaustive n=3..7, sampled n=8: 5K, n=9: 5K, n=10: 200), PROOF NEEDED
```

Original definition:

```
For a tournament T on n vertices, let OddCyc(T) denote the set of all directed odd cycles in T (as vertex sets). Define the **odd-cycle disjointness graph** G(T):
```

Suggested banner, leaving every existing byte below it as historical provenance:

> **CURRENT CORRECTION (2026-10-26).** The correctly defined odd-cycle
> collection formula is **PROVED for all tournament orders**; see
> [THM-002 / OCF](THM-002-ocf.md), external proof 1. Conflict-graph vertices
> are literal directed odd cycles modulo cyclic rotation, not their vertex
> supports: different cycles on one support remain distinct vertices. The
> historical “PROOF NEEDED” status and “as vertex sets” wording below are
> superseded. A five-vertex witness has four directed triangles and three
> directed 5-cycles: its literal conflict graph is `K7` and `H=15`, whereas
> collapsing equal supports gives `K5` and `I(K5,2)=11`. The body below is
> preserved as historical provenance.

Certificate, with vertex set `0,...,4` and rows giving adjacency:

```
0 0 0 0 1
1 0 0 1 0
1 1 0 0 0
1 0 1 0 0
0 1 1 1 0
```

The triangles, modulo cyclic rotation, are `(0,4,1)`, `(0,4,2)`,
`(0,4,3)`, `(1,3,2)`. The 5-cycles are `(0,4,1,3,2)`,
`(0,4,2,1,3)`, `(0,4,3,2,1)`. On five vertices all these odd cycles
intersect, so the literal conflict graph is `K7`. All three 5-cycles have
the same support, so support collapse gives `K5`. The independent checker
computes `H=15` by Hamiltonian-path recursion, separately from OCF.
The witness is minimal in order: up to four vertices an odd cycle has
length three, and each 3-support carries at most one directed cycle.

## 2. Producer report decoder contract clarification

Target: `session66-graph-pairs/REPORT.md`, section 8, after the reproduction
commands or before the verification list. Theorem J/I proofs need no change.

Suggested text:

> The supplied `decode_joint` and `decode_graph` routines are **valid-image
> decoders, not image-membership validators**: their contract assumes input
> was produced by the corresponding codec, up to relabelling. Malformed
> inputs require separate validation and need not be rejected by these routines.

No producer API change is necessary for its stated valid-image theorem.
The independent audit's `strict_decode` additionally recognizes each full
image and rejects all tested malformed graphs.

## 3. Minor inherited convention, not a new theorem defect

`THM-002-ocf.md` proof 1 writes `Dbar = D^op`. Irving–Omar's digraph
complement includes diagonal loops, so strictly `D^op` is `Dbar` with
loops removed. Hamiltonian paths use distinct vertices and no loop edge;
therefore their counts agree, and the inherited OCF is unaffected. If that
source sentence is edited, say “`Dbar` with its loops removed equals
`D^op`, so `ham(Dbar)=ham(D^op)=ham(D)`.” The Session66 report does not
assert the problematic literal complement equality.
