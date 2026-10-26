# Independent hostile audit: signed modular ports and clone defect

**Date:** 2026-10-26. **Verdict: PASS, with the stated scope boundaries.**
No mathematical correction to the sealed producer report is required.
This is an independent audit, not an arithmetic Frobenius identification or
a proof of Collatz convergence. No producer or repository files were edited;
no publication was performed by this auditor.

## 1. Exact audited input and provenance

Audited packet: `session66-modular/{REPORT.md,check.py,normal.json,
optimized.json,MANIFEST.json,SOURCE_MANIFEST.json}` and its seven archived
source blobs. An unchanged snapshot is in `input/`.

- Producer report SHA-256:
  `30ab70d0dbdaaf0ac722008f5d99cce371fb872034173247e7c047cb884b8753`.
- Producer checker SHA-256:
  `d0dd6d29b12ac645cd540d011db7dd012939852d613e7321ef6022b884fa029e`.
- Source commit: `63f52936f8dd8592d4f2fee869748baad2e1997e`.
- All 12 entries in the producer manifest match their declared byte lengths
  and SHA-256 hashes. Each of the seven archived sources additionally matches
  the exact Git blob at its declared path at the source commit, including its
  recorded Git blob ID. Thus this is not merely a self-consistent local hash
  check.
- `initial-git-verified.json` preserves the first Git-verified run. The final
  checker now defaults to a standalone offline replay: SHA-256, byte lengths
  and independently recomputed Git blob hashes authenticate the archived
  sources. Optional `--repo PATH` additionally checks the declared paths at
  the pinned commit; `git-verified-final.json` records that final check.
- Both producer executions were replayed using `-B`, once normally and once
  with `-O`. Both exactly reproduce both sealed output files, byte for byte:
  261, 902 authored checks.

The shared repository was used read-only. Its HEAD was already newer than the
source pin; no current-branch content was silently substituted for a pinned
source. The read-only startup packet is `startup.txt`. `AGENTS.md`, navigation,
the relevant tournament frontier, protocol and guardrails were read. The
targeted message command cannot run because this checkout has no `.machine-id`;
no registration or read-state change was attempted. That administrative limit
does not obstruct the mathematical audit. The parent owns repository freshness
and publication.

## 2. All-parameter arithmetic proof audit

### Signed inverse fibres: accepted

For every odd integer `n` and `b=±1`, `3n+b` is a nonzero even integer.
Valuation must mean `v2(abs(3n+b))`; negative sources are not exceptions.
Every target `q=(3n+b)/2^a` is odd and not divisible by three.

Conversely, for every odd `q` prime to three and every integer `a>=1`,
exactly one `b=±1` represents `2^a q mod3`. The numerator `2^a q-b` is an
odd multiple of three, so `n=(2^a q-b)/3` is an odd integer. Since `q` is
odd, `3n+b=2^a q` has exact valuation `a`, not merely valuation at least `a`.
The source and target have the same sign, including `q=±1`.

Multiplication of that equation by two gives
`(n,b)->(2n+b,-b)`, raising the valuation by one while retaining the exact
integer target. Iteration for **integer `j>=0`** gives
`(2^j n+b(2^j-(-1)^j)/3,(-1)^j b)`; the carry is integral because
`2=-1 mod3`. This is an inverse-fibre operation, not forward orbit conjugacy.
Negation instead gives `U_(-b)(-n)=-U_b(n)` at the same valuation and same
forward time. The report keeps these different operations separate.

### The signed carry and extra binary bit: accepted as an iff

Fix any odd `m>=3`, residues `r,s`, sheet `b`, and `a>=1`, with
`2^a s=3r+b modm`. Put `c=(3r+b-2^a s)/m`, which can be negative.
For **every integer height `t`**, not only positive heights,

`3(r+mt)+b = 2^a s + m(3t+c)`.

Because `m` is odd, division by `2^a` is integral exactly when
`2^a | 3t+c`. After that division, its parity is odd exactly when
`(3t+c)/2^a = 1-s mod2`. These are jointly equivalent to

`3t+c = (1-s)2^a mod2^(a+1)`.

The resulting odd target forces the source odd, so there is no missing
source-parity hypothesis in the converse. The corresponding direct source
condition is `3n+b=2^a mod2^(a+1)`, equivalently the report's inverse-of-three
guard. Using modulus `2^a` instead would only certify divisibility.

The hostile is uniform: if `3n+b=2^a y` with `y` odd, then
`3(n+m2^a)+b=2^a(y+3m)` and `y+3m` is even. This raises the valuation while
preserving both coarse source congruences. The same argument applies for
negative `n`; zero numerators cannot occur at an integer source.

CRT with any odd `m` joins every residue `r` to every prescribed exact
valuation guard, with arbitrarily large positive representatives. In
particular the two branches at residue `r=0` have target residues `b/2`
and `b/4`, distinct modulo every odd `m>1`. A deterministic fully odd map
on the bare odd-modulus residue table therefore does not exist. The report's
replacement of an unqualified odd-table alternation by the precise
inverse-exponent alternation is mathematically necessary.

## 3. The 14/15 gluing and cycle classification are proved, not extrapolated

An odd residue modulo 210 is uniquely determined by its residue modulo 105.
Consequently `(n,b)->(2n+b,-b)` is bijective: two is invertible modulo 105,
and the output has a unique odd lift.

The map `Phi(n,b)=3n+b mod315` is well-defined: changing `n` by 210 changes
`3n+b` by 630. Every image is prime to three. From its residue `z mod315`,
recover the unique `b=±1` with `b=z mod3`, then `n=(z-b)/3 mod105`, then
the unique odd lift modulo 210. This proves a bijection onto all 210 residues
prime to three, with explicit inverse. The identity
`Phi(L(n,b))=2 Phi(n,b)` proves conjugacy of these finite permutations.

CRT modulo `315=9*5*7` now gives an exact classification. The mod-nine
coordinate is a unit and has doubling period six. A nonzero mod-five
coordinate has period four; the zero coordinate has period one. A mod-seven
coordinate has period one or three, both dividing six. Thus every state has
period six if `5|z`, and twelve otherwise. The state counts are
`6*1*7=42` and `6*4*7=168`, giving **seven 6-cycles and fourteen 12-cycles**.
This is an all-state CRT proof with a confirming full census, not a guessed
period formula.

The analogous conjugacies for the factors are doubling on residues prime to
three modulo 21 and 45. They give `2^1 6^2` on the signed odd mod14 factor and
`6^1 12^2` on the signed mod15 factor. The correct product is a fibre product
over the **same sheet**. Of the 420 naive factor pairs, 210 have mismatched
sheets and cannot represent one signed source. Mismatch persists after the
lift. For example, `(1,+1) mod14` and `(1,-1) mod15` cannot be glued merely
because their integer residues CRT to1.

The distinction between finite and integer cycles is load-bearing. For any
nonzero odd integer source, an integer lift of length `j>0` raises its exact
valuation by `j`; it cannot return to the original signed source. All 21
finite-cycle representatives also have their nonzero integer wrap carries
explicitly checked. No integer cycle or terminal claim follows.

## 4. Signed pair-orbit proof and an actual unrestricted minimizer

For a permutation `sigma`, put a bit on each unordered pair, using a fixed
ordering of its endpoints. Transport acts by a permutation of these bits,
with a flip when the endpoint ordering reverses. For one pair-orbit of length
`ell`, XOR of its failure bits equals the XOR of these transport flips.
Conversely any failure pattern of that parity determines exactly two edge
assignments, by choosing its first bit and propagating. Different orbits use
disjoint edge variables. This proves the full signed orbit enumerator in the
report, as well as attainability of its minimum.

Negative twist can occur only when the two endpoints lie in one vertex
cycle and return swapped. If its length is `d`, swapping requires separation
`d/2`; hence `d` is even. Its antipodal pairs form exactly one negative orbit.
No cross-cycle pair can swap endpoints. Therefore

`min_T D_sigma(T) = number of even vertex cycles`.

This is valid for every finite permutation, including the empty set, fixed
points, odd cycles and mixed types. An explicit general constructor chooses
one arbitrary first edge bit per orbit and propagates with no failures until
the closing step. The last comparison fails precisely for negative twist.

Applying this constructor to the actual 210 signed ports gives a concrete
tournament, with all 21, 945 unordered pairs checked directly, of defect **21**.
It is encoded in `constructor-210.json`; no unenumerated optimum is being
reported only as a search expectation.

## 5. Actual antipodal pair modules: sharp restricted minimum 77

Let all vertex cycles have even lengths `2r`. Within each cycle the antipodal
pairs are disjoint and partition its vertices. Sigma takes pairs to pairs,
inducing a quotient permutation `pi` with cycles of lengths `r`.

Requiring each pair to be homogeneous means every exterior vertex sees its
two members identically. Since two distinct blocks are both homogeneous,
their four cross edges all have the same direction. Thus such tournaments
are exactly arbitrary quotient tournaments together with an independently
chosen internal direction in every pair. Existence and independence are not
extra unproved assumptions.

Under sigma every quotient failure contributes four original failures.
The internal edges form one orbit per original vertex cycle, returning with
their endpoints swapped after `r` steps. Its failure parity is odd, its
minimum is one, and the signed-orbit propagation construction attains it.
These internal variables are independent of all quotient variables. Hence

`D_sigma(T)=4 D_pi(Q)+I_sigma(T)`, and

`min_(pair modules) D_sigma(T)=c+4e(pi)`

for every such permutation, where `c` is its number of vertex cycles. The
quotient has an even cycle exactly when the corresponding original cycle
length is divisible by four. This proves the all-parameter sharp formula,
not just a lower bound or census pattern.

For the 210-state map, `c=21` and quotient cycle type is `3^7 6^14`. The exact
restricted minimum is therefore **21+4*14=77**. The audit independently
constructs an optimal quotient, combines it with optimal internal edges,
and directly checks:

- 105 disjoint antipodal blocks;
- all 21, 840 exterior clone equalities;
- quotient defect 14;
- internal defect 21;
- external defect 56;
- full original defect 77.

The certificate has a separate read-only verifier,
`verify_constructor.py`, which does not import the audit generator or any
producer code. It rebuilds adjacency matrices from the saved bits, checks
the actual affine permutation, antipodality, all clone equalities, quotient
lifting, and direct failure counts. Certificate SHA-256:
`20a795a6f4ec840f864dddbdcc5a06ea2c0def18a4fba8961e7a0ad30d57db67`.

The four-cycle hostile is decisive: unrestricted minimum 1 but module minimum 5.
The six-cycle control has module minimum 1 because its quotient is a three-cycle.
This identifies divisibility by four, not evenness alone, as the extra tax.
It does not assert that these blocks are maximal canonical clones or that
the native Collatz/Sierpinski objects possess them.

## 6. Ternary and arithmetic-Frobenius boundaries

For every `k>=0`, signed ports modulo `3^k` are conjugate to doubling on
the units modulo `3^(k+1)` by the same inverse formula. Binomial induction
gives `v3(4^(3^j)-1)=j+1`, hence the order of two modulo `3^(k+1)` is
`2*3^k`. The action is one even cycle.

The first nonzero ternary digit of a nonzero difference changes sign under
doubling. Reducing modulo `3^(k+1)` does not alter that first nonzero digit.
Thus the native restricted ternary tournament has every arc reversed and
defect `binom(2*3^k, 2)`, whereas optimization over all tournaments has minimum 1.
The auditor checked every pair for `k=0,..., 5`, including wrap boundaries.

The report correctly imports only the elementary permutation mechanism from
Session 65. An arithmetic Frobenius/genus identification would additionally
need the supersingular carrier, its actual Frobenius action and its geometric
theorem. None is supplied or claimed for this modular doubling permutation.
Likewise no preferred arithmetic orientation, positive orbit certificate,
descent theorem, integer-cycle classification, or LRC(14) conclusion is proved.

## 7. Independent computational universe and reproduction

The independent `check.py` imports no producer code. Normal and `-O` outputs
are byte-identical, with **595, 164 explicit checks** (exceptions, not asserts).
Their declared universes include:

- all odd sources from -5001 through 5001 on both sheets; negation and ordered
  lifts `j in {0, 1, 2, 3, 7, 16}`;
- all odd moduli 3 through 35, all residues and sheets, exponents 1 through 6,
  and every height `-2^(a+1) <= t < 2^(a+1)`: 325, 584 guard-iff tests,
  including 2, 703 negative carry cases and 7, 752 missing-top-bit hostiles;
- all 210 signed ports, all 420 naive factor pairs, and all finite cycles;
- all 132, 496 native ternary pair comparisons for `k=0,..., 5`;
- all 30 permutation conjugacy types on zero through six vertices, and every
  tournament in each type: 367, 966 tournament/permutation cases, comparing
  the **complete defect histograms** to the signed-orbit polynomial;
- every all-even cycle type on zero through eight vertices, and every
  quotient/internal assignment: 12 types and 5, 331 module orientations;
- actual 210-vertex free and restricted minimizing constructions, followed
  by a separate saved-certificate verification.

Conjugating a vertex permutation transports tournaments bijectively and
preserves defect, so one representative per cycle type is complete for the
invariant. The finite tests do not replace the all-parameter proofs above.
No search through all `2^21945` tournaments, no external code, and no
Lean/Coq formal verification is claimed.

Run from this audit directory; no repository, Git executable, or network is
required for the default mathematical and certificate replays:

```sh
python3 -B check.py --write-constructor > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
python3 -B verify_constructor.py > constructor-normal.json
python3 -B -O verify_constructor.py > constructor-optimized.json
cmp constructor-normal.json constructor-optimized.json
python3 -B input/check.py > producer-replay-normal.json
python3 -B -O input/check.py > producer-replay-optimized.json
cmp producer-replay-normal.json producer-replay-optimized.json
cmp producer-replay-normal.json input/normal.json
# Optional pinned-path provenance recheck, when a repository is available:
# python3 -B check.py --repo /path/to/math > git-verified-final.json
```

An additional frozen-location replay outside the workspace, with no sibling
`session66-repo`, reproduced the normal and `-O` receipts and the constructor
verification byte for byte; see `PORTABILITY.json`.

Audit concepts: exact binary guard; signed modular gluing; pair-orbit twist;
weighted quotient; internal/external independence. The relevant inherited
research pattern is **controlled forgetting requires a sidecar**: the common
sheet is required at CRT gluing, and weights plus internal holonomy are required
at clone merging. No new meta-pattern or novelty claim is needed.
