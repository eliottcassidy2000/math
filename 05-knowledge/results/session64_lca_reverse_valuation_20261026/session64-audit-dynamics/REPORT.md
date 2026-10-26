# Independent audit: punctured valuation coding, integer rays, and guarded debt

Session64, 2026-10-26. **ACCEPTED WITH ONE PROSE CLARIFICATION.** The arithmetic
proofs in the reverse-valuation packet are sound on their stated domains. The
independent checker passes in normal and optimized Python with identical bytes.
The root `CONNECTIONS.md` is accepted in the audited domains: sections 1–3 and
the dynamics-related limitations in sections 4 and 6. Tournament claims and
the finite-DAG counting theorem were not re-audited here.

## Useful result and correction

The parity-driven parameter replay independently recovers all 10,000 odd
returns, 19,939 parameter macros, final counter 9,884, final x bit length
11,669, and both producer SHA-256 receipts. The **odd-section minimum is 6**;
the minimum over **all parameter macros starting at the first odd section is
5**, attained at the very first odd macro `(a,b)=(6,23)->(5,36)`. If the initial
blank-derived parameter state `(0,4)` is included, its counter is of course 0.
These are different observation cuts, not contradictory computations.

At the parent's explicit authorization, only the last paragraph of producer
`REPORT.md` was changed to name the odd-section convention of JSON `min_a` and
record the independently checked macro minimum. Producer code and JSON were
not changed. Its three-file `SHA256SUMS` scope was preserved and refreshed.

- Original producer report SHA-256:
  `8be89a25ae2d91334f90f0f36dd244b89d5f610d7ac5db25cc76a889875f4f9e`.
- Corrected producer report SHA-256:
  `4c5aceeb0253ec258feda8bfa016d1322b6483fce07eb4fcb21bca9e4e4362ed`.
- Unchanged producer code SHA-256:
  `97743626055018e7d8f877753bafb2d2b1a333eb13e3dee23bd5aed616954465`.
- Unchanged producer JSON SHA-256:
  `e9fc0d6f0956bf7c10cef00a6a86652f0ec15ea76b482d59b529d520807e344d`.
- Audited root connections SHA-256:
  `987df5a36abfabd3e279145e47850a196ac8b50b94453620e1c9440217219dd1`.

## Proof critique and closed boundaries

1. **Exact coding domain.** The branches have respective multipliers
   `2^k/3` and `(2/3)^k`, valuation exactly k, and map odd units bijectively to
   the exact-first-letter cylinders. The complete finite preimage trees of
   `-1/3` and `+1/3`, including empty words, must be deleted. A numerator zero
   has infinite valuation, never a letter. In any cylinder the deleted set is
   still the complete countable exceptional set lying there, not merely its
   singular endpoint. The producer phrase “finite exceptional tail” is read
   as a tail reaching a singularity after finitely many steps, not a claim
   that the cylinder loses only finitely many points.

2. **Countable one-sided shift.** For every infinite positive word, A_n>=n
   forces both affine tails to vanish 2-adically. Each infinite tail is odd,
   so the first-branch numerator has its prescribed finite valuation; no
   coded point hits a singularity. Conversely nonsingular odd orbits provide
   the convergent identities and uniqueness. Length-n cylinders have modulus
   `2^(A_n+1)`; that extra oddness bit is essential. These cylinders prove
   continuity both ways with the product topology of the discrete countable
   alphabet. This is not an accelerated endomorphism of every odd 2-adic, or
   the different compact binary parity-shift statement.

3. **Isometry and LCA.** At the first unequal letters p,q, the two nonconstant
   terms have distinct valuations p,q, hence their difference has valuation
   min(p,q); the common prefix adds A. Thus the same-code conjugacy preserves
   `v2(x-y)=A+min(p,q)`, including equality/infinite depth. This preserves the
   labelled low-bit trie, not integer order or merely the number of shared
   valuation letters. Equality of shifted infinite tails is instead actual
   forward coalescence. The `3,13->5` and `1,-5/3->1,-1` controls separate
   these predicates. A rooted dynamical LCA needs its own basin and cycle
   convention. No universal positive-integer root basin is inferred.

4. **Integer inverse and rays.** Since `3x-1` is prime to 3, an integer image
   forces `v3(y)=k-1`; solving gives the unique candidate and its modulo-3
   integrality criterion. For x>1 the odd integer map strictly increases;
   for x<0 it strictly decreases; 1 is fixed. A backwards integer path stays
   on its sign side and strictly approaches a bounded integer endpoint, so
   terminates. Nontrivial components are unbranched one-sided rays. Common
   futures on one ray require comparability by iteration; Collatz integer
   sibling mergers are not transported as integer mergers.

5. **Eventual periodicity versus integer locus.** A periodic H code has its
   rational affine fixed point in `(1/3,1]`, by the real geometric sum;
   finite h-prefixes preserve that interval. Thus an eventually periodic H
   code cannot encode any integer except 1. This implication is valid because
   periodic tails make both completion sums equal the same rational fixed
   point. It does NOT identify the real and 2-adic sums for arbitrary words,
   or prove rational points must be eventually periodic. Constant words
   show ordinary point 1 goes to twisted 3/5, and ordinary -1 to twisted 1.

6. **Barrier guards.** `R(w)=1-min D_j` over both endpoint prefixes is the
   exact positive integer threshold. Its reverse recurrence and chronological
   monoid follow directly from the split prefix minimum. Executing n blocks
   asks only for j<n, while readiness at the returned odd section also asks
   j=n; a terminal counter zero executes the last block and then halts. Even
   parameter states may have counter zero without halting. Internal even
   macros only add resource, so odd-section guards suffice. Every finite
   valuation word has a positive odd integer realization >=5 by its residue
   cylinder; this strengthens the explicit all-1 bad-prefix family without
   making any assertion about one prescribed seed. Repeated formal words
   with nonnegative total drift need only one period's R; negative total
   drift eventually defeats every finite resource. Integer realization of
   an *infinite* periodic word is a separate restriction.

7. **Affine ledgers.** The H normalized height decreases but stays above
   x_0-1 for x_0>1: it is not a well-founded termination rank and does not
   control the counter. The C ledger retains a positive affine carry, giving
   the stated exact descent inequality; the coefficient sign alone fails
   at fixed point 1. No positive-integer Collatz convergence follows.

8. **Haar harmonic scope.** Normalized odd-cylinder mass is `2^-A`, hence
   iid letters with mass `2^-k`. With q^2+q=1 the kernel is exactly
   `1/[q(2-q^2)]=1`. Stop the bounded martingale at 0 or >=M; M consecutive
   ones guarantee exit from the finite live strip with uniform positive
   probability. Overshoot preserves the q^M bound. The increasing lower-exit
   events exhaust finite ruin paths, proving exact eventual ruin q^a. No
   asymptotic drift assumption or informal optional-stopping limit is needed.
   The fixed ordinary integer 27 has Haar mass zero: q^6 is not a probability
   assigned to its actual halting statement. Infinite geometric tails are
   retained, not silently truncated/renormalized, in both exact horizon checks.

## Independent finite controls

`audit.py` does not import, call, or copy producer functions. It uses explicit
exceptions, integers, and Fractions, with these alternative routes:

| Universe / route | Checked count |
| --- | ---: |
| All odd residues modulo 4096, direct forward iteration in both maps | 4,096 sources |
| Words of length <=4, alphabet 1..4, valuation sum <=11, against full residue-class census | 548 map/word cylinders |
| Closed-series singular preimages, replayed to numerator zero | 548 |
| Periods length 1..3, alphabet 1..4, both maps | 168 closures |
| All unordered pairs of those 84 periodic words, including equivalent periods | 3,486 paired cross-map metric checks |
| Odd signed targets -2047..2047, branches 1..12 | 24,576 inverse candidates / 2,048 targets |
| Odd signed sources -4999..4999 | 5,000 injective images |
| Backward paths of those targets, to root or fixed 1 | 2,028 integer steps |
| Words length 0..5, alphabet 1..3; counters 0..10 | 364 words / 4,004 guard comparisons |
| Same words, direct parameter replay including execution-only and returned guards | 4,004 native-parameter replays |
| All chronological cuts of those words | 2,005 |
| Odd positive sources 1..127, 12 steps, both maps | 1,536 affine ledger steps |
| Backwards horizon ruin recursion, horizons 1,6,12,24,48 | All reported masses reproduced |
| Parity-driven parameter iteration from (0,4), detecting rather than preselecting return lengths | 10,000 odd returns |

The inverse branch bound is exhaustive for the target universe: v3(|y|)<=6
there, so the only possible k is <=7; enumeration through 12 adds hostile
impossible branches. The horizon recursion is backwards dynamic programming
for `P_a(ruin by t)`, with the exact bound `P_a=0` for a>t; unlike the producer
it does not propagate forward live probability masses or safe sinks.

Replay receipts reproduced:

- Valuation CSV SHA-256:
  `99b0e9e1fcb53c64d24da5363dffc6d6a653be6d37ccc61d25a65fd88eacfb2e`.
- Final integer hexadecimal-encoding SHA-256:
  `4f24d6f0690f5d19c6da5c96f15a7af324d10f024bf8a75d16fb2b1498c06886`.

These are **FINITE-EXACT** implementation controls; the proofs, not finite
enumeration, carry the unbounded statements. The blank-orbit result remains
**FINITE_PREFIX_ONLY**.

## Sources, native-machine boundary, and reproduction

The closest proved prior mechanisms were inspected via targeted `git show`
at checkout HEAD `8941edc2fbde192cd9b1242dcb70ca87d09414f3` because sparse files
were absent. Their byte snapshots are `prior-valuation-lca.md` and
`prior-debt-range.md`, from
`05-knowledge/results/session56_rmq_lca_connections_20261025/`.
`prior-source-audit.md` records the inherited printed-clock discrepancy.
No repository file was edited.

The parameter targets and blank start match
[the primary paper, section 2.2](https://arxiv.org/html/2610.22603v2#S2.SS2),
inspected directly in this audit. **CITED, not rebuilt:** their identification
with the native Turing machine and blank tape. The displayed halt-family
clock `3b^2+9b-1` is negative at b=0 and is excluded. No corrected all-b native
clock or Lean rebuild is claimed. A parameter macro is not one native
Turing-machine transition. Antihydra blank-tape halting remains **OPEN**.

From the scratch parent directory:

```sh
python3 -B session64-audit-dynamics/audit.py > session64-audit-dynamics/normal.json
python3 -O -B session64-audit-dynamics/audit.py > session64-audit-dynamics/optimized.json
cmp session64-audit-dynamics/normal.json session64-audit-dynamics/optimized.json
(cd session64-audit-dynamics && sha256sum -c SHA256SUMS)
```

When the packet is relocated, pass `--packet PATH_TO_PRODUCER_PACKET`.
Both runs completed with exit code 0 and identical output. No denied or
uncertain actions, pending permission requests, or unresolved audit blockers.
