# Independent audit of saturated bounded-profile membership

## Verdict

**PASS for the mathematical corollary and the guarded finite controls.**
The proof correctly gives exact membership of every fixed positive rectangle
`1 <= B <= X, 1 <= s <= K`, conditional only on the cited, independently
audited finite-profile-interface theorem. It is an all-order membership
decision, not merely an order-bounded census. No proof of F(8,137), numerical
tail coverage, semilinearity, or polynomial binary-input running time follows.

Two scope clarifications should accompany promotion:

1. The stated search complexity covers bounded rectangle discovery and its
   shared predecessor DAG. Expanding all instruction words or graph adjacency
   data has output-dependent cost. In particular, the demonstration script's
   extra literal subset-counting validation is not covered by that bound.
2. The implemented state-limit stop returns `INCOMPLETE`. An unexpected
   exception, killed process, or uncontrolled memory exhaustion need not
   produce JSON and supplies no absence certificate. Only a normally exhausted
   search reporting `COMPLETE` certifies all missing pairs absent.

Neither clarification requires a change to the mathematical algorithm.
The author agreed to make both explicit in the reviewed release delta.

## Frozen object and dependency

The original candidate is `frozen_source/PROOF.md`, SHA-256
`82a43f3d6464dc05c3912842311ee1080aab28ae0e95dfd1714c2a2ad471eeba`.
The original implementation hashes are:

- `check_membership.py`:
  `497ed0e235dbf79b68f9c8c840e7c5005f23de314838c5ac1840ae3a8eb4c5eb`.
- `forest_runtime.py`:
  `571dc0a668a10f67ad9889bfe5d2079f4dde60089e82fc90cc9ffc7728703979`.

`INPUT_SHA256.json` records the original candidate and receipt hashes.
The dependency reviewed here is the published finite-profile interface at
commit `bf716a495ec2d49dbeb46e22798e37d43e4b0c26`, path
`session203/finite-profile-interface/PROOF.md`, with local proof SHA-256
`cb2de007a4de0a2446764e3b162024fa9f3cd5ea749c54342acab322ed1c6534`.
Its register-representative theorem, positive forest transfer, and weak-tail
preservation are cited proved dependencies, not inferred from this audit's
finite tests.

## Mathematical audit

### Saturation and merging

For nonnegative integer coefficients, saturation at `c` commutes with any
linear sum: a contributing coordinate at least `c` makes both saturated sums
equal to `c`; otherwise the contributions are identical. Coefficient zero
erases either representation exactly. Thus every legal transition and each
nonnegative scalar observation respects the quotient. No monotonicity of an
individual coordinate along a construction is needed.

The state key includes active labels. Equal keys have identical legal next
instructions and equal capped successor observations under every continuation.
This justifies merging even if their uncapped forest counts, tournament
orders, or hidden identities differ. It does not turn arbitrary formal
vectors into physical states: only vectors reachable from the empty graph
are admitted.

The bound `c = max(X,K)+1` puts every requested value strictly below the
saturation sentinel. A capped positive readout in the target rectangle is
therefore exact, while a readout equal to `c` is rejected as overflow.

### Positive B readout

The coefficient for each one-component forest not starting at the root is
one, and the coefficient for each one- or two-component forest with an end
at the root is one. Their overlap legitimately has coefficient two. These
are the literal `h-s` and `L` contributions, both counted without subtraction.
The implementation's `any` is correct because a physical root can be the end
of at most one vertex-disjoint path component.

The signed formula `h+L-s` must not be used after independently capping its
terms. The audit includes the arithmetic hostile `7-6=1` versus
`min(5,7)-min(5,6)=0`. The implementation uses the positive forest expression
for its capped search and uses the signed expression only for uncapped exact
literal witness validation.

### Active controls and roots

An introduction uses a currently free register and selects only active old
vertices as out-neighbors. Every forgotten old vertex beats every later new
vertex, matching the interface convention. Forget operations use only active
labels; labels may subsequently be reused for new vertices.

The search need not retain a root marker before its final observation. Choose
an active root, mark it, then forget all other ports if desired. Every forest
contributing to B or s has at most two components, so those additional forgets
cannot discard it through the two-bulk-start restriction. Root endpoint labels
are unchanged. Consequently observing with other active ports present is
exact, and the final root was never silently forgotten.

### Completeness and witnesses

The interface theorem provides a W(K)-register representative for every
actual positive bounded-start profile, possibly after its proved weak-tail
replacement. The search exhausts all legal words modulo a continuation-safe
congruence, so it reaches a state with that representative's capped readout.
Conversely each first-discovery predecessor trace is a legal construction.
Appending a root mark gives an actual rooted tournament, regardless of which
uncapped representative first reached the capped state.

The explicit state guard is conservative: it may stop before processing
other queued states or newly found successors, losing some positive answers,
but reports `INCOMPLETE`. It never upgrades an unfinished run to absence.
Readouts already recorded have independent physical witness traces.

### State and complexity bounds

There are `(W+1)^2` endpoint-pair types. At most W components start at active
vertices and at most two start in bulk, hence at most `W+2` components.
Counting all multisets of up to that size gives the claimed admissible-basis
upper bound `D <= binomial((W+1)^2+W+2,W+2)`. Invalid patterns only enlarge
this upper bound. The complete capped state has at most
`N = 2^W (c+1)^D` possibilities.

Per state, introduction actions number at most `W 2^W`, forget actions at
most W, and root observations at most W. Each source pattern has at most
`(W+3)^2` introduction choices and one forget image. Arithmetic operands
have `O_K(log(c+1))` bits, including bounded intermediate sums. Standard
indexed search stores coordinates, queue references, and predecessor indices
within `O_K(N D log(c+1))` bits. A deterministic dictionary adds only
fixed-K polynomial factors in `log(c+1)` to the claimed time bound.

Thus the core decision search is polynomial in numerical X at fixed K,
with enormous fixed-K degree and constants. It is not polynomial in the
binary length of X. A first-discovery witness word has fewer than N
instructions, but materializing its tournament or independently counting its
paths requires additional work, as addressed by the promotion clarification.

## Independent finite verification

`independent_audit.py` does not import the producer's algorithms. Its search
implements introduction through aggregate endpoint-type multiplicities, not
the producer's component-slot loops. Its literal forest check enumerates
directed edge subsets and rejects degree violations and cycles. Its scalar
witness check recursively enumerates all directed simple paths, rather than
using the producer's subset dynamic program.

The exact independent searches give:

| K | X | State limit | Result | States | Transitions | Profile count |
|---|---|---|---|---|---|---|
| 1 | 1 | 1,000,000 | COMPLETE | 145 | 340 | 1 |
| 1 | 5 | 1,000,000 | COMPLETE | 447 | 1,046 | 2 |
| 1 | 13 | 1,000,000 | COMPLETE | 1,084 | 2,534 | 5 |
| 1 | 137 | 1,000,000 | COMPLETE | 26,760 | 62,470 | 41 |
| 2 | 13 | 1,000 | INCOMPLETE | 1,000 | 2,015 | 6 verified |
| 1 | 1 | 1 | INCOMPLETE | 1 | 1 | 0 |

The complete 137 run agrees on every returned pair with both saved producer
receipts, not merely their totals. All 1,099 labelled tournaments of orders
1 through 5 were checked against literal forests, with 5,405 complete forest
vectors at full or successively forgotten active controls and 5,405 rooted
readouts. Every small literal bounded profile occurs in the completed search.

All 88 witness instances across normal/optimized 137 and both guard receipts
were independently replayed and recounted. This includes 330 legal register
reuse events and witnesses up to order 8. The optimized 137 receipt duplicates
the 41 normal witnesses; 88 is an instance count, not a distinct-source count.

Additional controls check 6,030 scalar saturation identities, runtime normal
and `-O` guards, exact 144/145-state stopping boundaries, all three invalid
nonpositive CLI arguments, and fresh normal/optimized full 137 receipts.
All assertions in the audit are explicit exceptions and stay enabled with
Python optimization.

The finite checks corroborate the general proof and implementation. They
are not the premise for the all-order completeness theorem. K=2,X=13 remains
an intentionally incomplete run, and no larger-budget absence is claimed.

## Reproduction and scope of workspace actions

From this audit directory:

```
python3 -B independent_audit.py > RESULT.json
python3 -B -O independent_audit.py > RESULT_OPTIMIZED.json
```

The two result objects should agree after removing `elapsed_seconds`.
`frozen_source/` makes the audit self-contained. No producer input or shared
math checkout file was edited. The shared checkout was already dirty and
historically behind the supplied published dependency, so it was neither
pulled nor rebased. The read-only startup was run; the targeted message
processor reported no `.machine-id` and was not registered by this audit.
No commit or publication was performed.
