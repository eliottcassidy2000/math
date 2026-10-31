# Independent audit: every bounded eight-core with five arbitrary tails

**Status: PASS. All1,287 complete trees independently reconstructed and
bound to the producer receipt in normal and optimized Python. Both mode
receipts are byte-identical, as are their control receipts.**

## Scope and prerequisites

The claim is for all1,287 eight-element subsets C of {1,...,13}, with
five arbitrary positive integer speeds, at radius1/14. Attachments may
repeat or coincide with core elements. The incoming nine-core/four-tail
certificate is a proved dependency, checked independently in ../audit_lrc.
No unrestricted LRC(14), optimal loneliness radius, or literature priority
is asserted.

## Independent conductor argument

Let J be a closed interval of length W. A speed v has x=vW periods on J.
Its bad-set measure, divided by W, is at most

    M(x) = [floor(x)/7 + min(x-floor(x),1/7)] / x.

This follows by taking floor(x) full periods and one residual interval.
On every piece M is increasing up to the residual peak n+1/7, then
decreasing. For k=2,3,4 and x>=k/7, the initial piece is at most1/k;
every later peak has value(n+1)/(7n+1)<=1/k because
(7-k)n>=k-1 for n>=1. For k=5 and x>=10/7, the first allowed n=1
piece decreases from1/5; all later peaks satisfy(n+1)/(7n+1)<=1/5
because n>=2. Hence the coefficients2/7,3/7,4/7,10/7 suffice.

The equality case is not discarded. If k>=2 bad sets covered J while
each had measure<=W/k, their summed measure would be<=W. A finite
relative-open cover of a connected interval, with no one member covering
it, must have two members meeting in a nonempty relative-open set.
Such an intersection has positive measure. Thus a cover would have
summed measure strictly greater than W, a contradiction. Closed safe
endpoints are essential. Positive surviving measure is not asserted.

For k=1 the coefficient1/7 follows separately from the open bad component
length1/(7v); a closed interval of that length retains safe endpoints.
The sharper general conductor is independently audited elsewhere in the
session. This argument verifies exactly the constants used here.

## Stronger inheritance rule, including repeated tails

If an assigned tail v belongs to {1,...,13} and is absent from C, promote
v to a ninth bounded core speed. If v is already in C, choose any missing
bounded speed d and enlarge C to C union{d}. In the latter case the
modified row has one additional constraint and contains the dropped
repeated v already. In both cases, the other four original tails remain
arbitrary positive speeds, including any fixed prefix. The nine-core
four-tail theorem supplies a witness for a row containing every original
constraint. This is a sufficient strengthening, not an equivalence.

No speed above13 can trigger this rule. Every noninherited finite branch
is visited with its exact sorted lower bound. All internal widths must
be positive, and every finite final leaf receives an original-clock
witness. The inherited bounded branch is not an omission of future
multiples or an argument from marginal safety.

## Independent computational path

The auditor imports only its own `geometry.py`, reusing the independently
written four-tail audit's full-comb intersection. Producer geometry uses
clipped tooth indices; the auditor generates the entire comb on[0,1]
and intersects ordered closed interval families. Finite leaf witnesses
are first full-intersection endpoints, rather than producer phase arrivals.

Every one of1,287 roots must appear in lexicographic combination order.
The auditor independently determines every cutoff, width, child branch,
inherited marker, endpoint witness, count and canonical full-tree hash.
The compact receipt does not determine the audit search universe.

Run:

    python3 -B check.py ../lrc/five_tail_pruned_receipt.json > result.json
    python3 -B -O check.py ../lrc/five_tail_pruned_receipt.json > result.optimized.json
    cmp result.json result.optimized.json
    python3 -B controls.py ../lrc/five_tail_pruned_receipt.json > controls.json
    python3 -B -O controls.py ../lrc/five_tail_pruned_receipt.json > controls.optimized.json
    cmp controls.json controls.optimized.json

Use the **pruned** five-tail receipt named above. Earlier exploratory
unpruned receipts are not inputs to this audit.

## Completed exact audit

-1,287 original eight-core families.
-37,610 four-tail,91,872 three-tail,315,277 two-tail and373,687 one-tail nodes.
-24,899 finite leaves, with original-speed witnesses.
-16,731 valid nine-core inherited branches, including padded repeated tails.
-All1,287 full-tree commitments checked, covering41,646,941 canonical bytes.
-Maximum cutoffs for1,2,3,4,5 remaining tails:84,80,81,79,100.
-3,012 exact closed-conductor controls and4,000 bad-measure envelope checks.
-Seven hostile/endpoint controls and eight rejected certificate mutations.

The receipt SHA256 is
`2dfbc2cb3e7e2969ffa8dedb5283637e4cd1db7bbe939d158dfed570770130cd`.

To overlap independent work with producer execution, the completed audit
used the equivalent split path:

    python3 -B rebuild.py /tmp/independent-families.json
    python3 -B bind_rebuilt.py ../lrc/five_tail_pruned_receipt.json /tmp/independent-families.json

Both commands were also run with `-O`, with byte-identical intermediate
reconstructions and final result receipts. `rebuild.py` calls the same full
independent `reconstruct` routine used by `check.py`; `bind_rebuilt.py` compares
every reconstructed record and aggregate. Binding alone is insufficient:
the reconstruction must be run first. The single-command `check.py` route
above reruns both stages without any saved intermediate dependency.
No external prover, Lean check, or human review is claimed.

## Lossless receipt transport

The full compact JSON is distributed as a gzip/base64 text archive to keep
repository API transfers small. From this directory, first run
`python3 -B ../materialize_receipts.py` to recover the exact named JSON files.
The decoder checks byte counts and SHA256. The proof still requires full
regeneration and independent replay; compression and hashes are not proofs.
The packet-level `verify.py` decodes into memory automatically and compares
the regenerated receipts byte for byte.
