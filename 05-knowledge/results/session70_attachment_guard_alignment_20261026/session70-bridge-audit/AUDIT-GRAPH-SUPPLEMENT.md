# Independent graph supplement: exact pairing and overwrite semantics

2026-10-26. **ACCEPTED** for Theorems A/B and the explicit Petersen/K5/K3,3
and gracefulness statements described below. This supplements, and does not
modify, `AUDIT-INITIAL.md` or its frozen evidence.

## Reviewed final producer artifacts

Directory: `session70-graph-labeling`.

- `REPORT.md` SHA-256:
  `3859b1c1e588da3945764ec840870981af9e7b5bd2f95114dedccb053f23208b`.
- `verify_graph_alignment.py` SHA-256:
  `824a72a986719d53983d2fb2bd12cc0a1f6fa43f2ba945db27d014a2c508d61a`.
- `receipt.json` SHA-256:
  `346900a91d5f5cc0392c8b4f41a8ba0811be7979c68c4bc22f175398b9bfa4ff`.

The final producer checker was independently invoked normally and with `-O`;
both outputs match that receipt byte-for-byte. The `--hostile` invocation
returns exit code 1 in both modes, with the expected error that bare K5
cycles do not all lift through Petersen ports. Saved stdout/stderr and
`graph-hostile.validation.json` record the actual outcomes.

## 1. Theorem A: exact named endpoint pairing

The final statement has the required hypotheses: r nontrivial vertex-disjoint
paths; 2r distinct named endpoints; mandatory old edges; fresh internally
disjoint connectors touching only exposed endpoints; intermediate linear
forest validity; final spanning-path acceptance. Its exact quotient is the
endpoint perfect matching, with `(2r-1)!!` states.

The independent all-order proof and 11,136 exact distinct-pair contexts are
already frozen in the initial audit. The final prose correctly excludes
singleton endpoint coincidences and unrestricted-to-fixed-Q_n transfer. It
also distinguishes terminal cycle closure and hidden interiors/counting from
the particular Boolean forest language. No additional restriction is needed.

## 2. Theorem B: signs really are unobservable under the stated requests

The theorem concerns the inherited factorial graceful family with
`m>=k^2+2k`, isolated leaves `b_i=i(k+1)`, and fixed colors `d>k`.
The fixed forest consists of one tree H and k isolated b_i. This claim can
also be checked directly: for each i<=k, the replacement fixed edge
`i--(b_i+i)` attaches i to the star core, while b_i loses its original star
edge. The index `b_i+i=i(k+2)` cannot equal any `b_j` for 1<=i,j<=k, since
`i(k+2)=j(k+1)` would force k+1 to divide i. The upper bound on these vertices
is `k(k+2)<=m`. Thus no isolated b_i is accidentally reattached by a fixed
edge, and the rest is a tree.

For every free color d<=k, the two edges incident with b_i are
`(b_i-d)--b_i` and `b_i--(b_i+d)`. Their other endpoints belong to H, because
distinct b_i differ by at least k+1. Both are inside the allowed interval.
All other color-d edges join two vertices of H and create a cycle if selected.
Consequently a completion is exactly one attachment of each color to a
different isolated leaf: `(pi,epsilon) in S_k x {-,+}^k`.

The accepted operation language is precise:

- overwrite one or two distinct color coordinates by specified actual edges;
- no-ops are legal;
- observe only resulting validity and future validity;
- success advances to the overwritten selection, failure leaves it unchanged.

An overwrite is valid iff its resulting leaf assignments are a permutation.
The signs do not enter the criterion. On success, the new assignment pi is
determined by the old pi and the request, independently of old signs; on
failure it is unchanged in both states. Induction proves equal pi states
indistinguishable under every request word. Conversely, overwriting color d
by a specified edge at leaf i is legal iff `i=pi(d)`, so one-step guards
recover pi. Hence there are exactly k! behavioral classes, each containing
2^k literal sign states. This is an all-order proof, not extrapolation from
finite testing.

No-op semantics are essential: a request equal to the current edge is a
no-op in one sign state and a genuine change in the other. Observing the
change indicator, forbidding no-ops, or observing literal output edges
distinguishes those states. Such protocols are correctly excluded. This
theorem does not prove graceful labelings exist for every prescribed tree
shape, nor equate reconfiguration of varying shapes with that conjecture.

### Independent actual-edge computation

`audit_graph_supplement.py` imports no producer code. It includes all actual
candidate edges of each free color, including the always-invalid H--H
contacts. It checks:

- k=1,2,3 at both m=k^2+2k and m=k^2+2k+2: all actual initial choices,
  all one-/two-color requests on every valid state;
- k=4, m=24: all 384 formula states and all one-color actual-edge requests;
- actual forest validity against the assignment criterion; accepted successor
  assignments; exact one-color profile classes; a change-indicator hostile.

The independent total is **112,122 explicit checks**, including **105,328
actual-edge overwrite requests**. Normal and `-O` receipts are byte-identical.
The finite universes demonstrate controls and boundary behavior; the proof
above covers arbitrary k,m and all future words.

## 3. Petersen contraction, suppression, and connectedness

The final report correctly distinguishes one-run quotient cycles from the
complete two-port response. Its local equations
`d_O+x=d_I+x=2`, `x in {0,1}` are necessary and sufficient for unique
2-factor lift. They allow one run with a spoke or two separate singleton
runs without it. Keeping only a simple K5 Hamilton cycle would omit the
latter case and would not prove full Petersen non-Hamiltonicity.

The report repairs that logical gap explicitly: classify perfect matchings
by the number of spokes, yielding exactly six, and verify each complementary
2-factor is two 5-cycles. Its explicit one-spoke matching and two cycles
agree with independent enumeration. The native graph has 120 unoriented
Hamilton paths in the replay; Hamilton path and Hamilton cycle targets are
not conflated.

The K3,3 suppression marks are also now explicit. The initial independent
audit's exact six-to-two refinement was incorporated and separately checked
by the producer: exactly two of the six quotient Hamilton cycles use all
three marked subdivided edges, and their expansions equal precisely the two
spanning cycles in Petersen-minus-0. This is an exact cycle-family
correspondence, not merely equal cardinalities. The marks are essential for
spanningness; planarity/minor predicates do not supply them.

## 4. K5 and K3,3 are also distinct graceful-labeling controls

The final added elementary comparison is correct. K3,3 with one part labeled
0,1,2 and the other 3,6,9 has distinct edge differences 1,...,9.
For K5, a graceful labeling would require exactly five odd edge differences
among 1,...,10. If a of its five labels are odd, the odd-difference count is
`a(5-a)`, which can only be 0,4,6. Thus K5 is not graceful. The independent
checker also exhausts all 84 normalized possible label sets containing 0
and 10. This has no implication that planarity decides gracefulness, nor
does it concern the Graceful Tree Conjecture's prescribed tree shapes.

## 5. Audit scope and replay

The full final producer program was replayed, including its Q_n and inherited
graceful examples. Independent second implementations here target Theorems
A/B, the native tournament rank family, Petersen port/subdivision behavior,
and K5/K3,3 gracefulness. External citations, the full old corpus, Gerbicz's
all-n glue construction, and the full Kuratowski characterization were not
independently re-audited here. The producer marks those distinctions rather
than promoting them from these finite controls.

```
python3 -B audit_graph_supplement.py > graph-audit.normal.json
python3 -O -B audit_graph_supplement.py > graph-audit.optimized.json
cmp graph-audit.normal.json graph-audit.optimized.json
```

`SOURCES-GRAPH.json` freezes the reviewed producer identities.
`SHA256SUMS-GRAPH` freezes this supplement's evidence. The independently
frozen initial manifest remains unchanged and is verified separately.
