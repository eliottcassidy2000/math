# Independent audit of the Hamilton-support image packet

**PASS, 2026-10-26, for I1–I5 with the stated inherited reconstruction
dependencies and finite scopes.** No mathematical correction is required.
This review does not establish all-order black-core connectivity or a
polynomial-time support-image recognition algorithm.

Reviewed frozen producer:

| File | SHA256 |
|---|---|
| `../image/README.md` | `291067c75a06f87a71540c94fda13a16cfb735e882707229933e5c0377e6b8f5` |
| `../image/check.py` | `3df577b1b1d859eaa692a175cb435b379ac1145a062d94e06a23d7288425ce06` |
| `../image/normal.json`, `optimized.json` | `b514815e682eb2765da02637873dd34f7b52beb21e5125eab1c3aee404fb3870` |

## Mathematical review

**I1.** The declared domain is a loopless oriented graph on a fixed vertex
set. In a genuine support, the ordered strong-component decomposition is
recoverable; prescribed endpoints in successive tournament components show
that every consecutive cross-component arc is traced, while no path can
use a skipping cross-component arc. The missing-pair rule is precisely the
inherited G1 decoder on the original support, without feeding newly inferred
arcs back into its counts. G1's cited structural dependency is retained;
the new note need not reprove it. THM-4807 separately certifies uniqueness
and endpoint reconstruction.

The converse direction of the image test is elementary and exact: D is a
subgraph of its completed tournament T, so HP(D) is a subset of HP(T).
Equal finite counts imply equal languages. The additional U(D)=D gate then
gives U(T)=D. Neither implication assumes that a locally decoded graph was
already an image. Nonnegative coefficients make equality at0 and1
equivalent to absence of every positive-degree coefficient. The empty graph
unit and singleton conventions are consistent.

**I2.** Cutting a Hamilton word at exactly its missing arcs gives one
ordered, vertex-disjoint spanning family of nonempty D-paths. Its order,
all endpoints and actual missing connectors uniquely recover the word.
Conversely any such compatible cover concatenates legally. Thus the
coefficient counts exactly the covers with the stated number of connectors.
The empty cover is handled separately. The four-cycle example has one new
word requiring both missing chords; individual chord tests are insufficient.
The example is explicitly outside the accepted decoder domain, avoiding
an unsupported claim about that more restricted class.

**I3.** The five-vertex example passes the local decoder and support
saturation but has defect polynomial7+2z. Its two added arcs really become
visible. The first-failure claim is finite and is supported by the complete
producer census through five, independently replayed here; the audit also
checks all smaller partial graphs by a separate implementation. Deleting
the middle vertex of the transitive triple gives an edgeless pair, while
every two-vertex tournament has a nonempty support. This proves failure of
induced heredity and rules out a forbidden-induced-subgraph description of
the full image class. The corrected deletion U(C(D)−v) is well-defined on
the image and commutes with successive deletions and relabelling by the
inverse identities. The insertion examples show both visibility directions.

**I4.** For two labelled tournaments, the complement operation agrees on
exactly the preserved path arcs and disagrees on every other pair. Hence
agreement being one spanning directed path is an iff, and that path is
unique. The blue condition is the inherited literal positional-reflection
anti-automorphism condition, not arbitrary self-converseness. Complement
preserves it. For each labelled target isomorphic to S, exactly |Aut(S)|
permutations transport S to it; uniqueness of the agreement path proves
the displayed division. Inverting the transport gives symmetry after
division by |Aut(T)|. Paths, presentations, literal lines, and the later
converse-merged graph remain distinct carriers; no alignment is inferred
from separate ordinary isomorphism types.

**I5.** Directly complementing the padded source leaves inner relation h,
forward arcs a→s and t→b, all other outer cross arcs backward, and b→a.
This is the stated framed hub and depends only on h and the ordered ports.
The preserved Hamilton path closes via b→a, proving strong connectivity.
Positional reflection swaps the outer singleton markers and restricts to
the inner positional reflection, so the colour test is exact. Two black
paths with the same ports therefore give a two-line walk; neither distinct
sources nor a route between different ports is asserted. The terminal-basin
contraction theorem was checked in its inherited source and is correctly
used only as a remaining connectivity reduction.

## Reproduction and independent controls

Both producer modes replay **76,618 gates** and match the frozen receipt
byte for byte. The producer exhausts all59049 oriented graphs at order5,
obtaining1624 locally decoded candidates,1504 support-saturated candidates,
480 false local images and1024 genuine images. These remain FINITE-EXACT.

`check_independent.py` imports no producer code. It uses edge sets, literal
vertex permutations, graph searches and separately enumerated ordered path
covers. Its **1682 explicit gates** pass normally and with `-O`, yielding
identical receipts. It covers:

- Every oriented graph through order4, compared with the actual support
  images of all tournaments, including empty/singleton cases.
- The five-vertex false image and its exact7+2z polynomial, independently
  by literal Hamilton words and by all compatible ordered path covers.
- The four-cycle joint-only polynomial4+z² and its unique new word,
  individual inert chords and the nonhereditary deletion control.
- All20 ordinary tournament-type pairs at orders3 and4, with every
  vertex permutation, verifying the agreement test and transporter factors.
- All16 Hamilton words of those representatives for the padded hub,
  strong-connectivity and colour-preservation checks.

Independent checker SHA256:
`1e346d4edc098b636f0579b5089152a06dd56aae87500f44e7b5adec99761b31`.
Independent normal/optimized receipt SHA256:
`9d0d1d4759e07c4fa98c84e0a4fd8f65320fe6664610fa511230045279ec5185`.

Reproduce from this directory:

```sh
python3 -I -B check_independent.py > replay.json
python3 -I -B -O check_independent.py > replay.optimized.json
cmp replay.json normal.json
cmp replay.json replay.optimized.json
python3 -I -B ../image/check.py > producer.replay.json
cmp producer.replay.json ../image/normal.json
```

A subsequent status/link-only promotion of the reviewed proof is covered
by this audit, with its hash lineage recorded by the producer. Any change
to a mathematical statement, proof, checker or finite universe needs renewed
review. No producer file was edited during this audit.

## Authorized status promotion

After PASS, the root requested promotion of the header, the I1 candidate
status label, and the audit link only. The resulting proof SHA256 is
`4c89b9c42fc0072e2aecacac2daf28fb0e30edd5bf8026519ce391b597331db8`.
Reversing precisely those status/link substitutions recovers the reviewed
candidate hash `291067c75a06f87a71540c94fda13a16cfb735e882707229933e5c0377e6b8f5`.
This was checked programmatically; no mathematical body, checker or receipt
changed. The producer manifest records this promotion separately.
