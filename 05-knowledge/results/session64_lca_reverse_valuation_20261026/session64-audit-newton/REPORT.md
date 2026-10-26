# Independent audit: native Newton clone-excess compression

2026-10-26. **Verdict: PASS within the stated scope.** No mathematical or
implementation correction is required in the reviewed Newton packet. A prose-only
clarification requested by the synthesis reviewer was applied, as recorded below. The
observation-image theorem and the prime-power LCA bridge are **PROVED** by the
arguments reviewed below. Native-engine validation here is **FINITE-EXACT**,
with limits recorded explicitly. No repository files were edited by this audit.

## Reviewed inputs and inheritance

Main production source: `../session64-counting/newton_modular.py`, SHA-256
`0342e4682a78de8505fde96b71f7d34c9038dbde3ace163793a3375c7c93ab6a`.
Reviewed proof: `../session64-counting/NEWTON_QUOTIENT.md`.
LCA proof/checker: `../session64-root/newton-lca/{REPORT.md,check.py}`.
Exact hashes of these and inherited source/proof files are in
`INPUTS_SHA256SUMS`. This audit pins the producer source before executing it.

Closest mechanism: the inherited integral gauge and constant-tail modular
quotient in `05-knowledge/results/port_residual_integral_20261025/counting`.
Canonical hostile: independent integral moment coordinates need not lie in
the polynomial evaluation image. Corrected near miss: polynomial degree
slots and arbitrary sparse exponents have different image lattices.
Least-used sidecar: the actual source exponents `m-b(D)` are nonnegative.
The audit board is evaluation kernel / integral transition / signed support /
terminal node / ordinary versus core / absolute precision. Anchor is native
arithmetic; niche is the exact lattice and LCA proof; wildcard is the scalar
substitution hostile below. No literature or priority claim is made.

## Proof review

1. **Full image, not only an upper bound.** The monic polynomials
   `q_j(X)=product_(i=1)^j (X-2^i)` give a unitriangular integral basis change
   from monomials. At row `r`, `q_j(2^r)=0` for `r<=j`; for `r>j` its
   valuation is `tau_j=j(j+1)/2`. Sequentially inspecting rows `1,2,...`
   proves that the *entire* kernel consists of coefficients divisible by
   `2^max(K-tau_j,0)` for `j<S`, together with unrestricted coefficients for
   `j>=S`. Earlier terms already vanish at every later row, so there is no
   cancellation gap in this induction. Thus the image really is the stated
   direct sum of cyclic groups of orders `2^(K-tau_j)`. Reducing that group
   modulo its doubles shows its minimum number of generators is the number
   of nontrivial factors; multiplying their orders gives the bit count.
   For degree `<h`, retain only `j<min(h,S)`.

2. **Sparse exponents are a different sublattice.** The degree-`<h` formula
   is not asserted for `h` arbitrary monomials. Independently exhaustive
   images at `K=S=4` have sizes `16` for support `{0}`, `8` for `{1}`, `1`
   for `{4}`, `128` for `{0,1,4}`, and `64` for `{0,2,5}`. Consequently
   `min(h_support,S,L(K))` is a valid adaptive storage bound, not a universal
   sparse-image invariant-factor formula. The producer makes this distinction.

3. **Native operation preserves the exact kernel.** The actual batch is
   `P'(X)=T*2^J0*X^delta*P(2^a X)`, with `a=ell*k`, integer signed `T`,
   and `J0,delta>=0` by the inherited integral gauge. At each destination
   node it observes a source node `a+r` in `1..S` and multiplies its value by
   an integer. This proves kernel stability for every legal batch, including
   excess, skips, coordinate deletion and incoming linear merges. It does
   not rely on positivity or on reconstruction of a discarded coefficient.
   The source and destination `b` are fixed by their exact keys.

4. **Terminal boundary and forbidden inverses.** A completing batch from
   `S>0` has `k>=1`, so the terminal scalar reads `P(2^a)` with `a>0`.
   It never asks the source quotient for `P(1)`. An already terminal state
   remains its separate scalar. No division by `2` or `X` occurs after
   truncation: modulo `4`, constants `4` and `2`, or monomials `X^2` and
   `X`, respectively witness why either proposed inverse would be invalid.
   The factor computation and final normalization use only odd inverses.

5. **Semantic support and merges.** At fixed committed degree `u` and
   cycle count `m`, the exact sign is `(-1)^(u-m)`. Thus same-`m` incoming
   integer weights do not cancel. Each legal nonzero batch shifts support
   injectively, and the permitted dense length decreases with the horizon.
   Keeping modular-zero sparse keys and zero nonterminal dense vectors is
   therefore essential to the stated adaptive-support argument. The code
   does so; independent intermediate-state checks verify actual instances.
   This changes weight coordinates, not literal clone ownership or ports.
   The clone-cut proof remains the semantic justification of pair-free
   core counts; no additional literal-clone or port quotient is introduced.

6. **Prime-power LCA bridge.** For prime `p` and `r>i`, the least-significant
   digit trie gives depth `v_p(p^r-p^i)=i`. Summing these depths yields the
   same triangular losses; the induction in item 1 works over `p^K`.
   Named nodes, their order and depths are retained; unit residues and
   native counting semantics are not supplied by a bare tree. The reviewed
   root proof does not claim a recurrence at other primes or a formula for
   arbitrary node configurations. Its checker was replayed separately,
   not substituted for validation of the native clone recurrence.

## Independent native validation

`check_independent.py` imports the actual `newton_target` implementation but
**does not call producer `controls` or reuse producer `audit_stage`**. It
replaces the stage callback with its own oracle. Its independent paths are:

- a fresh original `f(d)=sum gcd(c,d)` recurrence, using rational weights
  and products of geometric factors rather than the producer beta routine;
- an unmerged odd-partition Burnside oracle, with one geometric factor
  `(1+z^d)^(-1)` per cycle and no defect/moment/Newton carrier;
- exact integer divided differences at `2,4,...` for mixed Newton
  coordinates, not production `multiply_x` or `sparse_to_newton`;
- direct observations at every available positive node, not production
  `evaluate`.

Normal and optimized Python runs give byte-identical PASS receipts:

| Check | Exact finite coverage |
| --- | --- |
| Full native core totals | `N=0..48`; `K=1,2,3,4,6,9,12,16,24`; 441 comparisons |
| Unmerged partition oracle | 24,661 partitions, producing all ordinary/core totals through 48 |
| Actual intermediate nonterminal states | 105,201 state visits across those runs |
| Positive-node observations | 1,315,755 |
| Semantic zero sparse keys / zero dense vectors | 79,222 / 6,943 |
| Pre-final terminal state visits | 1,234, checked as separate scalars |
| Reachable legal batches | 4,044 across `N=1..24`, every distinct reachable `(S,D)` and legal `(k,e)` at each stage |
| Signed-excess / terminal batches among them | 2,293 / 2,132 |
| Native transfer / kernel node checks | 13,295 each |
| Unrestricted ordinary totals, separately | 49 exact totals through 48 |

The legal-batch test uses a seeded choice of precision `1..12`, arbitrary
signed polynomial coefficients on exponents `{0,1,3,8}`, exact factorial
quotients, and one seeded invisible Newton generator per batch. It checks
all destination nodes for each such batch. This is not a claim to enumerate
all arbitrary coefficient inputs or all kernel generators at all precisions.
The mathematical argument above supplies the universal closure statement.

At ordinary order 50, a fresh partition calculation and fresh `f` recurrence
agree. The absolute residue modulo `4` is `0`, but after the known division
by `2^25`, the normalized residue modulo `4` is `2`. Absolute 27-bit
arithmetic retains those two normalized bits; absolute 2-bit arithmetic
cannot. This is a hostile, not a normalized-Newton precision theorem.

## Literal semantic control and an additional hostile

`check_literal.py` imports no counting or clone implementation. It enumerates
all labelled tournament masks through order 6, forms their complete vertex-
relabeling orbits, and tests homogeneous pairs directly from the arcs. Normal
and optimized outputs agree:

- ordinary classes: `1,1,1,2,4,12,56`;
- pair-free core classes: `1,1,0,1,2,5,26`.

This extends literal coverage one order beyond the producer's own literal
control. It also exposes why ordinary counts must not be substituted for the
cycle inventory: the naive scalar series `A(z/(1+z))` first disagrees with
these core counts at order 4 (`0`, not `2`); order 6 gives `20`, not `26`.
Automorphism-weighted vertex inventories, and hence the exact signed
clone-excess construction, cannot be recovered by that naive scalar step.
The producer does not make this false claim.

## Exact prose-only amendment

At the parent publisher's request, section 2 of `NEWTON_QUOTIENT.md` was
clarified from “lower triangular unit matrix” to “lower triangular matrix
with unit diagonal entries.” A unit here means an odd invertible residue,
not necessarily the number 1: for example `q_2(8)/2^3=3`. This removes an
ambiguity without changing any mathematical statement, code, test input or
receipt. The amended proof is hash-pinned by `INPUTS_SHA256SUMS` and now
explicitly included in the producer's `SHA256SUMS`. The independent proof
review above already used “unit” in the invertible sense.

## Replays, reproduction and limits

The producer's normal and optimized receipts both replay byte-for-byte,
including its 405 complete counts and 75,304 intermediate checks through
order 44. Both LCA receipts also replay byte-for-byte: 588 diagonal/LCA
checks, 58,968 kernel-transition evaluations and 72 exhaustive image cases.
These replays are distinguished from the independently written native checks.

From the workspace root:

```sh
python3 -B session64-audit-newton/check_independent.py \
  --repo session64-repo --producer session64-counting/newton_modular.py \
  > session64-audit-newton/independent-normal.json
python3 -O -B session64-audit-newton/check_independent.py \
  --repo session64-repo --producer session64-counting/newton_modular.py \
  > session64-audit-newton/independent-optimized.json
cmp session64-audit-newton/independent-{normal,optimized}.json
python3 -B session64-audit-newton/check_literal.py > session64-audit-newton/literal-normal.json
python3 -O -B session64-audit-newton/check_literal.py > session64-audit-newton/literal-optimized.json
cmp session64-audit-newton/literal-{normal,optimized}.json
```

For another checkout/location, substitute its actual paths; no producer
source modification is needed. SHA manifests cover source, receipts and
this report. No runtime, memory, complexity, reachable-fibre optimality,
normalized-Newton precision, new unrestricted ordinary-state reduction,
port transport, or black-core connectivity result is certified here.
There is no unresolved audit blocker.
