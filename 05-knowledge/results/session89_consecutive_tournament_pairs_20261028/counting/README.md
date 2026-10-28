# Session 89: shared exact channels for consecutive tournament pairs

**Status: PROVED carrier/transfer specialization and FINITE-EXACT implementation
checks.** The formula and response quotient are inherited from Sessions 21/22;
the degree-bucket layout and terminal fold are inherited from Session 87. The
new implementation shares the response traversal and transfer metadata while
retaining independent positive channel weights. This is exact **counting**, not
class generation, transporter construction, or a new enumeration record.

## First distinguish the paired object

Fix an old vertex set V of size n and a named new root *. An **independent,
aligned pair** consists of a tournament T on V and a rooted tournament U on
V union {*}. Quotient by one simultaneous permutation of V, fixing *. The
restriction U|V is a second tournament layer, and the incidences with * are
one binary vertex attribute. Thus the object has two tournament layers and
one binary vertex bit; its count is F_n(2,1).

This is not the product of two separately unlabelled class counts, because
independent quotients forget their alignment. Already n=2 gives

```
A_2 = 1, rooted_A_3 = 4, A_2 * rooted_A_3 = 4,
aligned independent pair count F_2(2,1) = 8.
```

If the task instead requires **U|V=T**, there is only one independent
tournament layer. The count is F_n(1,1), i.e. the rooted order-(n+1) count.
There is no extra factor A_n. A concrete parent T also has its own extension
orbits under Aut(T); a global Burnside count does not produce those orbits.

## Inherited multichannel formula

For integers r>=1 and s>=0, take r ordered tournament layers on V and s
ordered binary attributes on every vertex. Every isomorphism acts on all
layers and attributes together. For an odd partition lambda of n, define

```
t = length(lambda)
E(lambda) = sum_i (lambda_i-1)/2 + sum_(i<j) gcd(lambda_i,lambda_j)
z_lambda = product_l l^(m_l) m_l!
F_n(r,s) = sum_(odd lambda partitions n) 2^(r E(lambda)+s t)/z_lambda.
```

A permutation with an even vertex cycle fixes no tournament: a suitable
power reverses the endpoints of an antipodal pair. On an odd cycle type,
there are E independent orientation bits per layer and t invariant vertex
bits per attribute. Burnside gives the formula. F_n(1,0)=A000568(n), while
F_n(1,1) counts rooted next-order tournaments. **r=0 is rejected:** without
any tournament layer even cycle types need not vanish. For example, binary
vertex-color orbits at n=2 are 3, whereas the odd-only sum is 2.

The exact F_n(r,s) carrier and integer transfer are already documented in

- `05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/README.md`;
- that packet's `channels.py` (the historical `a,b` names mean `r,s` here).

The same prior packet proves the self-converse reduction
`C_(2n+1)=F_n(2,1)` and `C_(2n)=F_n(2,0)`. Therefore the aligned-pair totals
coincide with an inherited self-converse channel. We do not claim a new
sequence, a new literature formula, or a newly constructed object bijection.

## Share states, not first moments

At descending odd cycle length l, the inherited state is the used degree d
and future responses `R_j=sum_(processed L) m_L gcd(L,j)`. Adding k l-cycles
changes the energy by

```
DeltaE = k(l-1)/2 + l*k(k-1)/2 + k R_l.
```

Let `O(d)=odd(d!)` and `v(k)=v_2(k!)`. For channel (r,s), its independent
normalized weight is `W_(r,s)=O(d) sum_histories 2^(rE+st)/z`. The transfer is

```
d' = d+lk
B = O(d') / [O(d) l^k O(k)]
W'_(r,s) += B W_(r,s) 2^(r DeltaE + s k - v(k)).
```

B is an integer, being the odd part of the labelled cycle-selection integer
`d'!/[d! l^k k!]`. The shift is nonnegative since r>=1, s>=0 and
`v_2(k!) <= binomial(k,2)`. Response update, response expiration, destination
key and B are all independent of r,s. All weights are strictly positive on
reachable states, so every allowed channel has exactly the same state
support. Signed-core channels are **not** supported by this argument/API.

`fused_prefix` uses one degree-bucket dictionary per stage with a vector of
independent channel weights at each response key. Its two-channel path is
specialized to avoid a Python channel loop on every transition. It computes
common degree metadata, response decoding, destination keys and hash lookups
once rather than separately. The l=1 terminal fold evaluates each channel's
own polynomial, sharing the degree multiplier. For k=0, response expiration
still runs; different old keys can merge even though arithmetic is identity.

No state quotient shrinks. For two channels the number of graph traversals
is halved, but big-integer arithmetic for both weights remains. The fused
weights use more peak memory than sequential scalar calls. This gives an
implementation speedup for a declared **batch**; it does not show that asking
for A_n in addition to F_n(2,1) is faster than asking only for F_n(2,1).

### Exact hostile merge / weighted Frobenius defect

The cores gamma=(3,7) and eta=(5,5) have the same degree 10 and terminal
response t=2 (take bound >=11 so the 1-coordinate survives). Their energies
and centralizers differ:

| Core | E | z |
|---|---:|---:|
| (3,7) | 5 | 21 |
| (5,5) | 9 | 50 |

For the contribution of these two histories, put
`M_j = 2^(5j)/21 + 2^(9j)/50`. Then

```
M_0 = 71/1050
M_1 = 6176/525
M_2 = 2778112/525
M_0 M_2 - M_1^2 = (2^5-2^9)^2/(21*50) = 1536/7 > 0.
```

Thus even the normalized shortcut `M_2=M_1^2/M_0` fails at a genuine merge.
The additional root bit contributes the common factor `2^t=4`, which does
not repair the loss. Squaring integer weights also squares their centralizer
weights incorrectly. The second moment must be maintained independently
along the shared graph, not reconstructed by squaring a completed ordinary
counter or its aggregated intermediate weights.

## Do not confuse the matrix carrier

An arbitrary n-by-n binary matrix under simultaneous row/column conjugation
is a different object, allowing loops and independently chosen opposite arcs.
For a permutation of cycle type lambda, its fixed matrix has
`sum_(i,j) gcd(lambda_i,lambda_j)` independent bits. For **odd** lambda this
number equals `2 E(lambda)+length(lambda)`, but matrix Burnside uses **all**
partitions. At n=2 the identity contributes 8 to the average and the
transposition contributes 2: matrix count 10 versus tournament-pair count 8.
The shared odd-type exponent is not an equality of the full orbit problems.

## Verification and reproducibility

All code is standard-library Python. Historical sources are read directly or
through one exact Git blob at commit
`03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`, with SHA-256 checks before execution;
no historical sources are modified or copied into this packet.

```sh
python3 verify.py --receipt verification.json
python3 -O verify.py --receipt verification_optimized.json
python3 benchmark.py --receipt benchmarks.json
```

In an external staging directory, add `--baseline-root /path/to/repository`.
The loader first tries the frozen source's path, and otherwise reads that
exact commit/path using `git show`. It does not fetch or modify the repository.

Both normal and optimized runs use explicit checks, not removable `assert`
statements. The controls include:

- Independent rational partition sums for (1,0),(1,1),(2,0),(2,1),(3,2),
  every n through 30.
- Fused two-channel agreement with S87 ordinary and generalized S87 paired
  scalar calls through selected bounds up to 160; exact matching per-stage
  state and transition telemetry. This is a moderate prefix, not a new record.
- S87 scalar color agreement for s=0,1,2,3 through 80, and inherited Session 22
  channel agreement through 30. The stronger S87-based scalar comparator is
  used for the central performance claim.
- Terminal-fold ablation, empty/one-channel/arbitrary-pair/duplicate-channel
  checks, and rejection of r=0, negative/nonliteral parameters.
- Literal orbit censuses for (1,0),(1,1),(2,1) at every n=0..4 and (3,2) at
  n=0..3; the paired n=4 universe has 65,536 labelled objects.
- Independent literal binary-matrix conjugacy through n=4, all-partition
  agreement, and the odd-type exponent identity through n=15.
- Product-of-classcounts and weighted-Frobenius-defect hostile controls.

Paired totals n=0..10 are

```
1, 2, 8, 88, 2752, 279968, 95458560, 111698291584,
457509297625088, 6662951988432581120, 349330527429800077778944.
```

See `verification.json` for complete literal-universe sizes, orbit-size
histograms, checked channels, and the hostile rational data.

## Same-output finite benchmark

The benchmark rotates/reverses execution order, performs seven repetitions
per method, verifies every returned value, and records both process CPU and
wall time. Both compared methods return **the entire prefixes 0..N for
(1,0) and (2,1)**. The strongest comparator calls the byte-identical Session 87
`degree_bucket_dp` for (1,0) plus its minimally r,s-generalized scalar version
for (2,1); a second comparator repeats both generalized scalar calls. The
older Session 22 channel engine is context, not the main speedup denominator.
Untimed `tracemalloc` probes measure Python allocations, not process RSS.

The following tables are generated from the completed `benchmarks.json` receipt.

Median process CPU seconds, seven repetitions per method:

| N | Repeated strongest S87-based | Repeated generalized S87 | Repeated historical S22 | Fused two-channel | Gain vs strongest |
|---:|---:|---:|---:|---:|---:|
| 80 | 0.02775 | 0.02783 | 0.04078 | 0.01967 | 1.41x |
| 120 | 0.15402 | 0.15772 | 0.24194 | 0.12053 | 1.28x |
| 160 | 0.74864 | 0.75702 | 1.17781 | 0.62711 | 1.19x |

| N | Repeated peak Python bytes | Fused peak Python bytes | Fused / repeated |
|---:|---:|---:|---:|
| 80 | 395,068 | 647,152 | 1.64x |
| 120 | 2,549,796 | 4,492,100 | 1.76x |
| 160 | 11,962,224 | 21,472,272 | 1.80x |

At N=160, every scalar channel has 844,358 mathematical
transitions, 27,506 peak states and 31,449 degree-transfer
metadata calculations. The fused method traverses that support once and computes
those metadata once; two repeated calls traverse it twice and compute
62,898 metadata. It still updates two independent big integers
per graph transition. These exact work counts do not themselves determine the
wall/CPU ratio. The diminishing measured gain is reported, not extrapolated.

A generic four-channel control at N=80, channels [[1, 0], [1, 1], [2, 0], [2, 1]], took
0.07101 CPU seconds for repeated generalized S87 calls and
0.05439 seconds for the generic fused implementation
(1.31x ratio). The dedicated two-channel
implementation is the main result; batching does not promise an arbitrary
number of channels will improve by the same factor.

## Practical consequence and stopping boundary

When both the ordinary and aligned-pair prefixes are required, use the fused
two-channel interface if the extra memory is acceptable. When only one count
family is needed, use the scalar interface. Reuse here amortizes shared
symmetry metadata; it does not reconstruct individual classes, their stabilizers,
or their alignment. Canonical generation remains a different algorithmic task.
No larger prefix or new record was pursued.
