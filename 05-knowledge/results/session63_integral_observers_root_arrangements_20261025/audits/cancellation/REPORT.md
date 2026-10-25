# Independent cancellation-depth and tournament audit

2026-10-25. **PASS** for the producer's scoped Theorems 1–4, the concrete
50-vertex realization, and the independent order-join family. No mathematical
correction is required. Finite tests support the proofs; they do not replace
the all-parameter arguments. No novelty, black-move, HN, or LRC claim is made.

## Provenance correction and recovered identity

The initial audit inherited the producer's statement that the historical
50-vertex source had not yet been located. That statement is now superseded.
Incoming checkpoint `a971f9b4fbdf409c31ba61ad9e8b3b00648536cf` contains:

- `05-knowledge/results/unprotected_ports_cancellation_20261025/valuation/README.md`:
  the complete regular elementary-abelian spectrum and its integral proof;
- `05-knowledge/results/unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md`,
  section 4: the concrete 50-vertex witness;
- `05-knowledge/results/unprotected_ports_cancellation_20261025/valuation/INDEPENDENT_REVIEW.md`:
  the prior independent audit.

I read these sources. The source connection set is the negative of the
producer's set. Thus `u -> -u` identifies the cores and `(u,b) -> (-u,b)`
identifies their ordered-pair inflations. This is the recovered witness up
to explicit isomorphism, not merely another example with matching numbers.
The regular spectrum and witness are **RECOVERED CONCURRENT RESULTS**,
independently rederived and verified here, not new discoveries.

`provenance_audit.py` independently parses the source's literal connection
set from its J7 formula, rather than importing the producer's copied set.
It verifies all 625 core and 2500 inflated adjacency entries, all six source
file SHA256 pins, and the source convention `Q(u)=-P(-1-u)` including its
modulo-five coefficients `[4,3,3,4,1]` in degrees 21 through 25.
Normal and optimized outputs agree byte-for-byte.

The producer's all-odd-action parity/prime-support interface, general
integral averaging precision criterion, and synchronized diagonal spectrum
are separate scoped extensions not stated in the inspected incoming packet.
This is a bounded source comparison, not a claim of historical novelty.
The cyclic-action comparison is likewise already in the incoming source.

## Mathematical audit

1. For an odd-order action, all cycle lengths are odd. Every cycle factor
   has a positive linear coefficient at `t=-1`, so the rational root order
   is the minimum cycle count, without leading cancellation.
2. Reduction cannot lower this order. It raises it exactly at prime
   divisors of the first nonzero integral shifted coefficient. Every
   subset orbit has odd size, giving `P(t)=(1+t)^n mod 2` coefficientwise.
   Complement reciprocity forces the remaining degree parity in odd
   characteristic. Empty and identity actions have the stated boundaries.
3. Writing the numerator coefficient as `|A|` times the integer orbit
   coefficient proves the valuation-e criterion. Keeping numerator data
   modulo `ell^(e+1)` is sufficient. Division by a nonunit in a field is
   nowhere used.
4. For `E=(C_p)^k`, every nonidentity element has `p^(k-1)` cycles of length
   p. In the integral shifted expansion the normalized j-term has
   p-valuation `j-1-v_p(j)`. For odd p only j=1 survives modulo p, giving
   exact depth `p^k-p+1`. At k=1 this equals the rational depth 1: no jump.
5. In s synchronized copies the lower bound becomes
   `j-1+v_p(s)-v_p(j)` and the j=1 coefficient is s. Therefore characteristic
   p depth is `sp^k-p+1` when p does not divide s, and `sp^k` otherwise.
   Independent product actions instead multiply orbit polynomials and add
   root orders. Equal point-orbit counts do not identify these actions.
6. The core is regular and hence pair-free: a homogeneous directed pair
   would have outdegrees differing by one. The inflation has only its
   intended pairs. Every larger transitive module would contain overlapping
   adjacent homogeneous pairs, so these are the canonical transitive blocks.
7. Rooted refinement is discrete and translations are transitive, proving
   the entire core automorphism group, not merely a subgroup. In the
   inflation the degree-25 layer induces the core; unique homogeneous mates
   determine the degree-24 layer. No hidden internal automorphisms remain.
8. The regular core is strongly connected; its inflation is strongly
   connected. In an order join, strong components are fixed setwise, giving
   the independent product automorphism action. An outneighbor of a vertex
   inside the earlier component separates it from a vertex in a later
   component, ruling out cross-component homogeneous pairs.
9. The invariant-XOR rank is the number of gap orbits, not a root
   multiplicity. Polynomial coefficient reduction does not recount actual
   invariant-space dimensions. The companion integral specialization
   framework owns the invariant-lattice base-change statement.

## Independent exact evidence

`audit.py` imports no producer code. It uses direct convolution of shifted
cycle factors and exact integer division; independently it constructs the
literal action's cycle index and takes binomial moments. Its refinement uses
incoming-neighbor block counts, independently of the producer's ranked
outneighbor colors.

- 47 arithmetic cases: `(p,k)` in `(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),
  (7,2),(11,1),(13,1),(17,1),(19,1)`, `s` in `{1,2,3,p,p^2}`, deduplicated,
  with `sp^k <= 800`.
- Every case checks characteristic zero and all 14 primes through 43:
  **705 exact depth checks**, plus the complete characteristic-p polynomial.
- Nine odd-action controls: cyclic orders `1,3,5,7,9,15,21,25`, and the
  nonabelian order-21 affine action on seven points. These test the dyadic
  law, prime support, and odd-characteristic reciprocity.
- The even transposition has polynomial `1+t+t^2` and characteristic-two
  depth zero, demonstrating the oddness boundary.
- Three independent C3 copies have characteristic-three depth 3, while
  three synchronized copies have depth 9.
- The full concrete 25- and 50-vertex adjacency objects are checked, including
  every homogeneous pair and all 25 translations. Root partition counts
  are `2,3,15,21,25`. No full `2^25` subset census is claimed.

The concrete result is

```
full automorphism group: (C5)^2
number of gaps: 25
invariant-XOR rank: 1
subset orbits: 1,342,208
root depths (Q,F5,F2,F3): (5,21,25,25)
rational leading coefficient: 3000
```

The cyclic order-25 hostile has depths `(1,5)` over `(Q,F5)`, not `(5,21)`.
The nonabelian order-21 action has leading coefficient 2 and no jump at
its bad-denominator primes 3 or 7. Thus neither group order nor denominator
primes alone determine cancellation depth.

Fresh independent normal/optimized replays agree byte-for-byte. Fresh
producer normal/optimized replays also agree with each other and with the
producer's saved 80-case output. All tests use explicit exceptions, not
Python assertions removed by optimization.

## Reproduction

Python standard library only:

```sh
python3 -B audit.py > replay_normal.json
python3 -B -O audit.py > replay_optimized.json
cmp replay_normal.json replay_optimized.json
python3 -B provenance_audit.py --repo <repo-root> > provenance_normal.json
python3 -B -O provenance_audit.py --repo <repo-root> > provenance_optimized.json
cmp provenance_normal.json provenance_optimized.json
```

For the producer rerun, use the producer's unchanged `check.py`; the audit
manifest pins that input and the audited result text. Producer source-recovery
checker results are additional evidence, not substitutes for our source-parsed
isomorphism check. `SHA256SUMS` pins this frozen audit packet.

## Execution incident, reconciliation, and stopping state

The first shell attempts stalled. A pure-V8 BigInt computation nevertheless
completed independently, with 47 arithmetic cases, the tournament checks,
nine odd-action controls, and the even hostile. Its numerical verdict was
valid; the accompanying initial file write failed after 259.5 seconds.
No earlier stalled shell execution was claimed to have completed.

After environment recovery, the failed destination was confirmed absent and
no matching pending processes remained. The delayed `audit.py` existed;
two old outputs were empty. Fresh normal/optimized runs were then performed
and verified. A late-discovered complete producer output matched the fresh
rerun but was not retroactively treated as a completed contemporaneous audit.
The interrupted output files were removed after replacement by fresh verified
receipts. The single bounded persistence retry succeeded. `INDEPENDENT_V8.json`
is explicitly a summary transcribed from the completed tool output, not an
invented recovered in-memory store; `audit.py` and the fresh replay files are
the portable reproducible evidence.

No mathematical blocker remains. All artifacts are scratch files; no
repository edits or pushes were made by this audit. Publication is owned by
the root/publisher task.
