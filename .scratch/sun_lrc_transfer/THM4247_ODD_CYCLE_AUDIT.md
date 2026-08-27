# Audit of THM-4247 and the updated LRC residual odd-cycle artifacts

Verdict: mathematical/count/hash PASS. No maintained file was edited.

## THM-4247

All six normalized source/output SHA-256 values equal the theorem
frontmatter. The three normal replays reproduce the frozen output. The exact
projection rows sum correctly:

```text
degree 34: 36,288 - 1,536 - 576 = 34,176
degree 42: 16,992 -   672 - 192 = 16,128.
```

The deleted rows are disjoint. The projection identity forces the companion
visible degrees `(124,12)`, `(4,132)`, `(156,12)`, and `(0,168)`, so the
degree-four fibre and pure-hidden exclusions are applied to the correct rows.
The 24-vector degree-twelve hidden shell is exactly two free 12-element
orbits under target units and `T`; `T^2=-omega` makes the displayed two-coset
orbit enumeration complete.

The attachment logic has the needed direction only: common full-map image
implies hidden projection zero on all twelve attachments. The two symbolic
orbit calculations retain `t(t^2-1)` in the reduced denominator; the inherited
degree-three denominator bound removes the raw orbit-zero degree-eight
factor. The reciprocal common roots are therefore only `t=+/-1`, which give
`Z/U=0` and are outside the gate. No converse or wall closure is claimed.

Concrete non-consequence-bearing source corrections:

1. `jc23_w0_hidden_degree12_attachment_audit_thm4247.py` lines 20--21 still
   say the script exhausts degree-42 vectors. `main()` actually enumerates the
   24 degree-12 vectors and two representatives.
2. The same file's finite-field comment says `d=tD(t^2), deg D=6`; in this
   degree-12 lane the asserted denominator degree is three, hence `deg D=1`.
3. `denominator_degree_budget()`, `enumerate_degree_42()`, and
   `rank_one_resultant_certificate()` are unused degree-42 remnants; the
   live degree-12 output is hard-coded rather than obtained from the helper.
   The error message `degree-42 vector became zero` is likewise stale.

These should be cleaned to prevent a future auditor from mistaking unused
degree-42 code for a dependency. They do not change any current assertion or
theorem conclusion.

## LRC residual and odd-cycle scripts

The cross-postprocessor replays exactly and its normalized hashes match
THM-4242:

```text
181,162 - 36 = 181,126
removed ray = {(50,r):590<=r<=625}
36*C(30,9)=515,057,400 body obligations
fixed-50 remainder = 592 - 36 = 556, max other endpoint 589.
```

Both odd-cycle programs replay their frozen outputs. The NetworkX and
standard-library paths agree on:

```text
739 nonisolated vertices, one connected component
181,126 edges
30,912,074 triangles
737 triangle-covered vertices
triangle-free edges = bridges = {(616,756),(616,760)}
leaves = {756,760}; sole articulation = 616
remaining 737-vertex block connected with no articulation
fixed-50 degree 556; neighbour-induced edges 150,861
neighbour-induced triangles 26,807,003
edge/triangle digest
d682ec7b70b7eecc641559e26d0afa1e23b0ede3ff90c025a637e6c02e208623.
```

The bridge conclusion is valid: an edge in a triangle cannot be a bridge,
and the only triangle-free edges terminate at the two leaves. The
737-vertex remainder is biconnected because it is connected, has more than
two vertices, and has no articulation point. The scripts correctly label
this as a proof-obligation graph, not a safety graph.

One presentation clarification would help: `max_core_size=512` in the
NetworkX output means the number of vertices in the maximum `k=500` core,
not that the degeneracy is 512. The preceding `degeneracy=500` is correct.
