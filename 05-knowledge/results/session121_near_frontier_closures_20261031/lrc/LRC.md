# Three arbitrary attachments: a genuine next finite-core layer

**Status:** PROVED by an analytic conductor plus a FINITE-EXACT complete
certificate, independently replayed in normal and optimized Python.
No claim of unrestricted LRC(14), or literature priority, is made.

## Inheritance and nearest frontier

The incoming `clique_palette_joint_ports_20261030/LRC.md` L4 proves all positive
pairs over all 78 eleven-element subsets of `{1,...,13}`. Its parent
`precision_lifts_ramsey_deficits_20261030/LRC.md` L2 supplies the two-tail
conductor and the exact one-tail first-arrival mechanism. The closest actual
extension is not dilation of those cores: it is **three independently varying
positive speeds over ten-element cores**. The canonical hostile remains
`{1,2}` at radius `1/4`: speeds 3 and 12 attach individually but not jointly.
The necessary sidecar is the complete common closed interval set, never three
marginal success bits.

This packet neither changes earlier maximum-radius statements nor replaces
any primitive bank by a dilation count. All earlier endpoints, primitive core
lists, and their source certificates remain intact.

## T1. A many-tail conductor

Write `G_C(r)={t in [0,1]: ||ct||>=r for every c in C}`. Suppose it contains a
closed interval of width W>0. Let k>=2, `0<r<=1/(4k-2)`, and let the k positive
integer speeds be sorted `u_1<=...<=u_k`. If `u_1>=ceil(1/W)`, all k speeds
attach simultaneously, in the original clock. In fact there is a positive
measure set of successful times inside that core interval.

**Proof.** Choose a subinterval J of length `1/u_1`. The first speed is bad
on a set of measure exactly `2r/u_1`. For any speed v and interval of length
ell, partition from the left into `floor(v ell)` full periods and one residual
piece. Each full period has bad measure `2r/v`; the residual has at most
`2r/v`. Consequently its bad measure is strictly less than
`2r ell+2r/v`: if the residual length is positive, the complete periods have
length strictly below ell; if it is zero, there is no residual bad measure.
Apply this to each `u_j`, j>=2, on J. The union of all bad sets has measure
strictly below

    2r/u_1 + sum_(j=2)^k (2r/u_1 + 2r/u_j)
      <= (4k-2)r/u_1 <= 1/u_1.

The complement therefore has positive measure. Every point of that complement
satisfies all actual core and attached inequalities. Repeated speeds do not
invalidate a union bound. Equality in the stated radius threshold is allowed
because the intermediate estimate is strict. This is a sufficient conductor,
not a claimed optimal radius threshold. QED.

At r=1/14 it covers k=2,3,4. Only the k=3 finite-core families are certified
here; the k=4 conductor alone is not a four-tail family certificate.

## T2. Complete finite certificate for three tails

Fix r=1/14 and a core C with largest safe interval width W>0. Sort any three
attachments as u<=v<=w. Put B=ceil(1/W).

1. If u>=B, T1 finishes all three.
2. For every integer 1<=u<B retain the **complete** common cover
   `J_u=G_(C union {u})(r)`, including closed endpoints and isolated points.
   If its largest width W_u>0, put V_u=ceil(1/W_u). If v>=V_u, the inherited
   two-tail conductor finishes v,w. Alternatively T1 with k=2 gives the same
   cutoff here, since 1/14<=1/6.
3. For each u<=v<V_u retain the complete cover
   `K_uv=G_(C union {u,v})(r)`. If its largest width W_uv>0, put
   `T_uv=ceil(2r/W_uv)`. Every w>=T_uv succeeds because the longest open bad
   gap for speed w has length `2r/w`; a closed interval at least this long
   meets its closed safe set, including equality.
4. Exactly the finite triples `u<=v<=w<T_uv` remain to be checked literally.

All three positivity conditions are checked, not presumed. A point-only or
empty branch is an explicit stopping failure for this certificate and must
not be silently omitted. Thus successful finite checks prove the all-positive-
integer statement, not merely a sampled prefix.

## T3. The bounded structural family

The supplied complete enumeration tests every

    C subset of {1,...,13}, |C|=10,

and all finite branches in T2. When its result `all_triples_accepted` is true,
it proves, for **all positive integers u,v,w**, a common time t satisfying

    ||s t|| >= 1/14  for every s in C union {u,v,w}.

There are exactly 286 labelled cores; each is primitive (ten selected integers
from this range cannot all share a divisor greater than one). Overlaps and
repeated added speeds are included. This is not an assertion for arbitrary
ten-speed cores. It strictly enlarges the stated prior bounded-core coverage:
for example C={1,...,10} and attachments 101,103,107 are in this class but
contain only ten members of {1,...,13}, so the eleven-core two-tail theorem
cannot directly apply.

The safe intervals in the receipt are complete at the exact claimed radius,
not substitutes for maximum-radius optimization. No maximum-radius theorem
is asserted or weakened.

## Exact computation and independent controls

`triple_certificate.py` constructs a restricted joint cover by enumerating only
tooth indices meeting an actual core interval. It retains rational endpoints
with `Fraction` and uses explicit exceptions, not removable assertions.
Every first-attachment cover is independently compared with the source full
comb compiler. Every finite third-attachment decision is independently compared
with a literal full-comb intersection, and every returned first-arrival time
is checked against **every original speed**, not merely the last speed.
The separate auditor independently replayed all second-attachment widths and
thresholds without importing producer or copied source geometry.

Reproduce with:

    python3 triple_certificate.py --output triple_receipt.json
    python3 -O triple_certificate.py --output triple_receipt.optimized.json
    cmp triple_receipt.json triple_receipt.optimized.json

The copied `source_lrc_banks.py` is unmodified source from
`clique_palette_joint_ports_20261030/lrc_banks.py`, read at source snapshot
0578f16726e4e535189fb6074cc2ddb024c95d82. Earlier source attribution and checks
remain valid; no relabelling as new mathematics is intended.

**Remaining frontier.** A four-tail family over nine-element cores has a
conductor from T1 but still needs its full finite tree. The more important
unrestricted obstruction remains a structural uniform bound on the complete
phase geometry of arbitrary core shapes. Neither this finite core census nor
its ordinary dilation closure supplies that bound.

### Exact older-family comparison

`THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion.md` already
certifies the same 286 *base sets*, but only physical rows `3C union T` with
all three tails prime to 3. `THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays.md`
propagates those chart families to `3hC union T`, still with ternary-unit tails.
The present physical row is **C union {u,v,w}**, with no divisibility
restriction. Its quantifiers therefore are not another corollary that counts
the same dilation rays. The older statements and their sharp endpoints remain
unchanged. Their primitive-body reduction is still relevant to the unrestricted
clock-three entry frontier, which this bounded undilated theorem does not solve.

### Completed producer results and compact evidence

All 286 cores passed. The full deterministic finite tree has 13,693 first
joint covers, 574,336 second-joint width branches, and 67,585 final finite
queries. Maximum cutoffs are 137 (three-tail), 314 (two-tail), and 82
(one-tail). The minimum initial width is `9/1232`, attained at
`(1,2,3,5,6,7,8,9,11,13)`, agreeing with THM-4442's independently inherited
geometry without using that geometry as an input.

The separate controls pass 480 literal complete-comb comparisons, 1,827
all-large tests exactly at the conductor's radius boundary for k=2,3,4,
and four hostile/closed-point tests. Normal and optimized controls agree.

The full 12MB receipt is deliberately regenerable rather than required as a
large checked-in blob. `triple_compact_receipt.json` retains **every original
core and all its original safe endpoints**, each initial cutoff, counts, and
SHA256 commitments to each complete finite branch tree and the entire ledger.
A hash alone proves nothing: run the producer and independent auditor to
reconstruct the committed evidence. To regenerate the compact receipt:

    python3 compact_receipt.py triple_receipt.json triple_compact_receipt.json

The full local ledger remains retained as research evidence. The publication
should include the deterministic generator, compact receipt, normal/-O
agreement, and independent audit rather than silently omitting the finite
branches from its proof obligation.

An equivalent thirteen-distinct-speed formulation is: **every row having at
least ten distinct speeds in `{1,...,13}` is safe at radius `1/14`.** Consequently this
single three-tail theorem contains the earlier 78 eleven-core all-pair
families as special cases, while admitting rows with exactly ten bounded
speeds that those families did not directly cover.

The completed normal and `python3 -O` full ledgers are byte-identical, as are
their regenerated compact receipts and the separate controls. `verification.json`
records script/source/receipt hashes; independently audited status should only
be promoted after the separate auditor's completed replay, not from these
producer comparisons alone.

For a compact-only replay, the producer still constructs and checks the full
finite branch tree in memory, then writes only the compact receipt:

    python3 -B triple_certificate.py --compact --output regenerated.compact.json
    python3 -B -O triple_certificate.py --compact --output regenerated.optimized.compact.json
    cmp regenerated.compact.json triple_compact_receipt.json
    cmp regenerated.optimized.compact.json triple_compact_receipt.json

Expected compact SHA256:
`e264a1ab7b31a3c7e1ca082c46ebddedd212594df568a638b4d8960a3c620bb0`.
No full ledger is printed to stdout in either mode; stdout contains progress
and summary only. Omitting `--compact` preserves full-ledger reproduction.

The companion `../audit_lrc/AUDIT.md` reports independent PASS in normal and
optimized Python: all 574,336 second-width full-comb reconstructions and
67,585 finite leaves, 1,827 conductor-boundary and six hostile/endpoint
controls, five rejected certificate mutations, and all 286 compact branch-tree
hash bindings. Its normal/optimized summaries are byte-identical.
