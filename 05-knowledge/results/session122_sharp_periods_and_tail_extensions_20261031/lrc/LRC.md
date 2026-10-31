# Four arbitrary tails over every bounded nine-core

**Status:** PROVED by an analytic conductor and a FINITE-EXACT complete certificate, independently replayed in normal and optimized Python. Not a claim of unrestricted LRC(14), optimal radius, or literature priority.

## Inheritance and scope

The closest proved mechanism is Session 121 T3: every ten-element set
C contained in {1,...,13} admits three arbitrary positive integer attachments
at radius 1/14. The canonical hostile remains the radius-1/4 core {1,2}:
attachments 3 and 12 work individually but not jointly. The indispensable
sidecar is the complete common closed interval set, including isolated points.
A safe interval's width controls a sufficient all-large-tail regime, not the
entire common-time predicate.

The new target is every nine-element set in the same bounded pool, with
**four** arbitrary positive integer attachments. There are 715 such cores.
Repeated attachments and attachments already present in the core are allowed.
This is a genuine next layer rather than a count of dilated copies.

## P1. Partial-period conductor

Let a core's common safe set contain a closed interval of width W>0. Let
k>=2, 0<r<=1/(4k-2), and positive real speeds u=u_1<=...<=u_k. Set

    c_k = 2kr / (1-2r(k-1)).

If u W >= c_k, all tails are simultaneously safe on a positive-measure
subset of that interval. For integer speeds, the sufficient inclusive cutoff
is ceil(c_k/W).

**Proof.** The radius assumption gives 0<c_k<=1. Select an interval J of
length c_k/u from the core interval. The bad set for u on J has measure at
most 2r/u because J is no longer than one u-period. For any other speed v,
partition J into full v-periods and a residual interval. Its bad measure is
strictly less than 2r|J|+2r/v: for a positive residual the full-period length
is strictly less than |J|, and for a zero residual the extra term is unused.
There is at least one such other speed. Thus the union of all bad sets has
measure strictly less than

    2r/u + (k-1)(2r c_k/u + 2r/u) = c_k/u = |J|.

The surviving set has positive measure. The proof retains actual time and
works with coincident speeds. QED.

At r=1/14 this gives the exact constants used by this certificate:

    k=1: 1/7; k=2: 1/3; k=3: 3/5; k=4: 1.

The k=1 statement is separate: every closed interval of length at least
2r/u meets the closed safe set of speed u, since each open bad component has
length exactly 2r/u. Equality can produce only endpoints, so the one-tail
case claims existence, not positive measure.

These constants are sufficient, not optimal. Any sharper conductor proved
elsewhere in this session does not alter the constants or finite universe of
this sealed certificate.

## P2. Complete branch theorem

Fix a nine-core C and sort four attached speeds nondecreasingly. Maintain
its exact common cover J after the already assigned prefix, the number k of
unassigned tails, and the lower bound lo equal to the last assigned speed
(initially k=4, lo=1). At every node verify W=max width(J)>0, and put
B=ceil(c_k/W).

Every next speed v>=max(lo,B) and every subsequent speed are handled by P1
(or its single-tail version). For every integer lo<=v<B:

1. If v is in {1,...,13} but absent from C, the enlarged set C union {v}
   has ten **distinct** bounded speeds. The other three attached speeds,
   both previously assigned and still unassigned, are arbitrary positive
   integers. Session 121 T3 discharges this entire branch. The receipt marks
   it explicitly as `inherit_ten_core`.
2. Otherwise, when k>1 compute the full exact intersection with speed v's
   closed safe comb and recurse with k-1 and lower bound v.
3. When k=1 find the first safe arrival in the complete cover, compare with
   a separate literal full-comb intersection, and check that the returned
   rational time satisfies every original speed.

Empty and point-only internal covers are hard failures, not silently pruned.
The inherited branch is a logically proved enlargement, not a statistical
heuristic or a missing-speed filter claiming future irrelevance. Every finite
integer branch is present in the tree. Therefore a successful complete replay
proves the all-positive-integer four-tail theorem.

## P3. Result and evidence

The producer is `four_tail_certificate.py`. It uses `Fraction` throughout,
explicit exceptions rather than removable assertions, and publishes one
compact record per core. Each record retains the original core, all exact
safe endpoints, initial width/cutoff, every level's node counts and largest
cutoff, inherited-branch count, full-tree byte count and SHA256 commitment.

A hash is not a proof. The deterministic generator reconstructs every branch
and checks every leaf. An independent auditor must reconstruct the complete
trees using independently implemented full-comb geometry and bind all 715
hashes. The producer's optional `--full-directory` writes all committed trees
for direct inspection; they are not required large repository blobs.

    python3 -B four_tail_certificate.py --output regenerated.json
    python3 -B -O four_tail_certificate.py --output regenerated.opt.json
    cmp regenerated.json four_tail_receipt.json
    cmp regenerated.opt.json four_tail_receipt.json

The copied source helper `source_lrc_banks.py` is unchanged from Session121's
attributed helper, ultimately `clique_palette_joint_ports_20261030/lrc_banks.py`.
Its full-comb routines supply initial cores and a second finite-leaf decision;
the independent auditor must import neither helper nor producer.

The conclusion strictly contains Session121's ten-core/three-tail families:
a row containing exactly nine bounded speeds and four speeds above 13 is
now admitted by the completed certificate and independent audit. This does
not classify arbitrary nine-speed core shapes or settle unrestricted LRC(14).

## Completed verification

All 715 cores passed in normal and optimized production; the compact receipts
are byte-identical. The separate auditor reconstructed 34,077,244 bytes of
canonical full trees independently and bound all 715 compact commitments.
Both independent modes agree. Exact totals:

- 715 roots and 23,281 three-tail width nodes;
- 218,081 two-tail and 407,207 one-tail width nodes;
- 44,133 finite leaf witnesses and 52,347 inherited branches;
- maximum cutoffs for one through four remaining tails: 72,94,104,92.

The independent controls pass 4,428 conductor-boundary cases, seven
hostile/endpoint checks, and reject eight certificate mutations. See
`../audit_lrc/AUDIT.md` and its exact receipts.

Compact receipt SHA256:
`dca3d9be25f87c8e9cf60417f823efc526defe65b8ef003f276f12636da1d920`.

## Lossless receipt transport

The full compact JSON is distributed as a gzip/base64 text archive to keep
repository API transfers small. From this directory, first run
`python3 -B ../materialize_receipts.py` to recover the exact named JSON files.
The decoder checks byte counts and SHA256. The proof still requires full
regeneration and independent replay; compression and hashes are not proofs.
The packet-level `verify.py` decodes into memory automatically and compares
the regenerated receipts byte for byte.
