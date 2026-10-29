# Exact contraction of forced-zero cycle stages

**PROVED recurrence preservation / FINITE-EXACT checks / measured
implementation experiment.** No new response quotient or asymptotic bound.

## Source recurrence and schedule

The inherited [Session87 prefix evaluator](../session87_four_vertices_enumeration_20261028/FAST_COUNTS.md)
uses the current degree d, future gcd responses h_j, and an integer weight
scaled by odd(d!). It processes odd lengths p in decreasing order. A
multiplicity k changes d to d+kp, adds k gcd(p,j) to each remaining response,
and applies the exact integer cycle-selection and orientation factor.
The final unit stage is already folded into a polynomial evaluation.

For target prefix through N, put r=N-d. If p>r, every valid continuation
has k=0 at that stage. The weight is unchanged; only an inaccessible
response coordinate is discarded. After inserting k copies at length p,
send the new state directly to

\[
q=\max(1,\text{largest odd integer at most}\min(p-2,N-d-kp)).
\]

An agenda stores buckets by q and degree. It merges equal remaining
response keys before their next active stage. The unit stage at q=1
handles the boundary r=0 as well, with only zero multiplicity available.

**Theorem.** This schedule gives exactly the same complete prefix, for
every N>=0 and every positive literal vertex-color count, as the inherited
descending recurrence.

**Proof.** Every skipped stage has forced zero multiplicity, hence identity
weight transfer. The successive coordinate deletions compose to direct
restriction to j<=q and j<=r. All intervening paths are unique; none adds
a multiplicity factor. Histories reaching the same degree and surviving
response at q have identical allowed continuations and factors, even if
they arrived from different earlier stages. Adding their weights commutes
with every later linear transfer. Induct on descending q, and use the
unchanged exact unit-tail fold. QED.

The schedule contracts proven identity paths. It does not use a guessed
future-zero test or discard any composite cycle length that can fit.
Response packing and odd-factorial normalization are inherited, not new.

## Matched experiment

Both methods return **all coefficients0 through N**. Three runs per method,
alternating order; no concurrent mathematical workload was deliberately
launched. Other system activity is uncontrolled. Python's Windows CPU timer
has visible quantization, so these small differences are not strong timing
evidence.

| N | Reference transfers | Scheduled transfers | Reference median CPU s | Scheduled median CPU s |
|---:|---:|---:|---:|---:|
|100|90219|87000|0.03125|0.03125|
|150|601461|582239|0.23438|0.21875|
|200|2994814|2905899|1.32812|1.28125|

The exact reduction at200 is88,915 transfers, about2.97 percent. Exploratory
single runs at300 were essentially tied at about32 wall seconds; they
are not included as a repeat-controlled benchmark or used to claim a gain.
`benchmarks.json` retains both clocks, all samples, source hash and output
hashes. No memory reduction is measured; peak *active* scheduled states
would not measure all pending agenda states.

This is a modest exact work reduction, not a replacement for the current
production method on performance grounds. The more useful new recursion
in this packet is [precision truncation](PRECISION.md), which changes the
requested output and must not be benchmarked as if it produced full counts.

## Structural boundary

The schedule is safe because the numeric capacity proves a unique zero
continuation. A native complement or re-marking cannot be skipped merely
because a scalar statistic is unchanged. The81/27 transporter and the
earlier six-vertex interval-family exit remain contrary controls. The
source group action and ordered ports must be retained for structural
enumeration. This scheduling theorem applies to the Burnside arithmetic
graph, not automatically to the black/blue metagraph.
