# Prefix-free proof scheduling with an atomic integer coverage ledger

2026-10-05. **PROVED:** the input code and normalized atomic weights;
conditional fairness; exact certificate transport and finite-coverage tests.
**FINITE-EXACT:** the declared scheduler experiment. **OPEN:** decay of the
full unresolved mass to zero, equivalently universal Collatz termination
for this complete literal search. No optimal universal-machine or
algorithmic-complexity claim is made about the implemented checker.

Artifacts: [program](../../04-computation/experiments/atomic_prefix_certificate_20261005.py)
and [saved output](atomic_prefix_certificate_20261005.out).

## 1. Inheritance and the coordinate that must survive

The closest proved mechanism is already explicit in
[the fixed-integer audit, section 4A](crossroads_poset_20260926_integer.md):
a full-support atomic mass retains each integer, whereas Haar averaging
or Abel normalization can erase it. Its exact decay-to-zero reduction is
not a proved decay estimate. The canonical hostile there is a surviving
singleton, and the plateau \(Z_2=Z_3\) forbids assuming strict decrease at
every clock tick.

The practical inheritance is the retained literal evidence and grounded
closure in [adaptive observation union](adaptive_observation_union_20261004.md),
conditional fairness in [fair frontier extension](fair_frontier_extension_20261004.md),
and explicit point versus schema obligations in
[the finite-seed compiler](finite_seed_receipt_compiler_20261005.md).

[Reimann's survey, arXiv:2408.05121v1](https://arxiv.org/abs/2408.05121v1)
is the source of the present prefix-complexity/dimension question. The
results below are elementary proved constructions, not applications of a
paper theorem to Collatz. The exact dyadic measure and local-dimension
calculation are developed separately in
[algorithmic Collatz measure](algorithmic_collatz_measure_20261005.md).

The board is **input atom / proof code / exact source / reusable suffix /
fair service / unresolved obligation**. Its key type boundary is that a
short proof's code weight and its input's coverage weight are different.

## 2. A normalized, computable weight on every odd input

Index positive odd sources by \(n=2m-1\), \(m\ge1\). If
\(k=\lfloor\log_2m\rfloor\), let

\[
\gamma(m)=0^k\operatorname{binary}(m),\qquad
\ell(m)=2k+1,\qquad
\mu(\{2m-1\})=2^{-\ell(m)}. \tag{1}
\]

The initial zeros determine the length of the following binary payload,
so these codewords are prefix-free. Block \(2^k\le m<2^{k+1}\) has
\(2^k\) entries of weight \(2^{-(2k+1)}\), hence block mass
\(2^{-(k+1)}\). Summing gives exactly one. The tempting weights
\(4^{-(k+1)}\) sum to one-half and need a factor of two.

Every positive odd integer has a computable, strictly positive atom. For
example \(\mu(\{27\})=1/128\). Removing that singleton preserves natural
density one and removes a Haar-null set, but still removes exactly
\(1/128\) of this input measure. No density-one or dimension statement
can substitute for accounting for that mass.

For the first \(M\ge1\) indexed inputs, put \(k=\lfloor\log_2M\rfloor\).
The exact unlisted tail is

\[
\mu\{2m-1:m>M\}
=2^{-k}-(M-2^k+1)2^{-(2k+1)}. \tag{2}
\]

In particular, the first \(2^K-1\) inputs leave mass \(2^{-K}\).
These formulas concern the **input** space, independently of any asserted
Collatz behavior.

## 3. Prefix-free receipts are a separate space

An actual first-hit certificate for \(n=2m-1\) is a finite tuple of exact
positive valuations \(a_1,\ldots,a_r\). Encode it as

\[
\gamma(m)\,\gamma(r+1)\,\gamma(a_1)\cdots\gamma(a_r). \tag{3}
\]

The parser first reads the source and the number of exponent fields,
then exactly those fields. It rejects trailing bits, malformed fields,
wrong valuations, a nonroot endpoint, and any suffix after the first
visit to one. Before the arithmetic validity filter, the whole syntax
is prefix-free with Kraft sum

\[
\left(\sum_m2^{-\ell(m)}\right)
\sum_{r\ge0}2^{-\ell(r+1)}
\left(\sum_{a\ge1}2^{-\ell(a)}\right)^r=1.
\]

This does not mean proof-code mass measures input coverage. Even the ROOT
receipt for input one is `11`, of code weight \(1/4\), whereas input one
has atom \(1/2\). Other receipt languages may give several descriptions
of the same source. The ledger therefore verifies a receipt, extracts its
exact source, and counts that source atom once. Resubmitting the same
proof never increases coverage. This parser is a concrete decidable
certificate checker, not an optimal universal prefix-free machine, and
its lengths are not asserted to be shortest description lengths.

## 4. A fair, explicitly budgeted scheduler

For each indexed input \(m\), keep a resumable literal search task for its
first visit to one. Stage \(T\ge1\) activates every task with
\(\ell(m)\le T\), and reserves

\[
2^{T-\ell(m)}
\]

quanta for each unfinished task. A quantum observes at most one odd
Collatz edge, preserving its exact valuation; if its current value
already has a verified route, it splices that suffix instead. Every
reserved stage is finite and

\[
\sum_{\ell(m)\le T}2^{T-\ell(m)}<2^T. \tag{4}
\]

This is a bound on logical odd-edge quanta, not constant bit-time
arithmetic: numbers and retained words can grow. The saved output also
counts actual queries, distinct observed edges, and separate certificate
verification work. Skipped reservations are not charged as executed
queries.

**Conditional completeness.** If a source has a finite odd first-hit rank
\(r\), its literal task terminates after at most \(\max(1,r)\) quanta.
Through stage \(T\ge\ell(m)\) it has been offered

\[
2^{T-\ell(m)+1}-1
\]

quanta. Thus, writing \(D=\max(1,r)\), it is certified no later than
\(T=\ell(m)+\lceil\log_2(D+1)\rceil-1\). Suffix reuse can only finish
it earlier. Every finite batch of actually terminating sources therefore
finishes at a finite stage. No assumption about the other tasks is needed.

The same Kraft scheduling argument applies to an effectively advertised
prefix-free set of finite-action proof searches. It does not require
knowing in advance which search will succeed. In particular, it does
not silently assume access to the halting domain of a universal machine.
The implementation uses the explicit gamma-coded source tasks above;
it neither constructs a universal-machine halting probability nor gives
a termination bound independent of the still-unknown rank \(r\).

Verified suffixes are proved descendants of a completed route. Unfinished
observations are retained separately. Applying the inherited reverse
reachability closure from ROOT to their exact edge graph can recover
additional certificates, but an ungrounded cycle never certifies itself.

## 5. What a mass bound actually proves

Let \(C_t\) be the finite set of distinct verified sources and

\[
c_t=\sum_{n\in C_t}\mu(\{n\}),\qquad R_t=1-c_t.
\]

This is an exact rational lower bound on covered input mass and an exact
mass of inputs not yet in the ledger. It is not an assertion that those
remaining inputs fail to converge.

For the first \(N\) indexed inputs, the smallest atom is

\[
\delta_N=2^{-(2\lfloor\log_2N\rfloor+1)}.
\]

**Finite coverage certificate.** If \(R_t\le\delta_N\), all odd inputs
\(1,3,\ldots,2N-1\) are certified. Otherwise any missing one contributes
at least \(\delta_N\), and the finite ledger leaves infinitely many other
positive atoms unlisted. Its residual would therefore be strictly greater
than \(\delta_N\). Equality is sufficient for this finite-ledger theorem.
If an arbitrary **possibly infinite** certified set were allowed, only
the strict inequality is sufficient in general: the complement of one
minimum-weight input has residual exactly \(\delta_N\) and still misses it.

The scalar bound can be inefficient. If the known sources are precisely
the first \(2^K-1\) inputs, using this scalar criterion to infer the first
\(2^J-1\) needs \(2^{-K}\le2^{-(2J-1)}\), or \(K\ge2J-1\).
That is a quadratic prefix overhead. Direct knowledge of the verified
prefix already certifies its larger range; the scalar test does not
improve that fact. A localized ledger instead computes the missing mass
inside the requested interval, which is zero iff all its sources are
certified.

**Audit correction.** The first draft used the conservative strict test
and incorrectly used an infinite complement-of-one set to claim strictness
was necessary for a finite ledger. The repaired API uses \(\le\); a
ROOT-only ledger is an explicit equality control. The dyadic overhead
boundary is correspondingly \(K\ge2J-1\), not \(K\ge2J\).

For this sound conditionally complete search, monotone convergence on a
countable probability space gives

\[
\lim_t R_t=\sum_{n\text{ has no finite ROOT route}}\mu(\{n\}). \tag{5}
\]

One can prove (5) without a measure-theory theorem: first truncate the
input sum at \(M\), where only finitely many indicators occur, and then
use the computable tail (2). Positivity of every atom gives

\[
\lim_tR_t=0\quad\Longleftrightarrow\quad
\text{every positive odd input reaches one}.
\]

This is the inherited full-support reduction in a prefix-code coordinate.
No bound tending to zero has been proved here. In particular, neither
Kraft normalization, fair service, nor observed finite mass gains proves
that all tasks terminate.

Finite verified seeds and universal family schemas must still be kept
distinct. A receipt for one seed adds its one atom. A proved infinite
constructor may supply many further individual receipts or a separately
justified family mass; the ledger does not count every member merely
because a family name occurs in a finite assumption manifest.

## 6. Exact finite experiment and the cost comparison

Run from the repository root:

```text
python -B 04-computation/experiments/atomic_prefix_certificate_20261005.py
python -B -O 04-computation/experiments/atomic_prefix_certificate_20261005.py
```

The scheduler runs stages one through fourteen, initially knowing only
ROOT. Stage fourteen activates the first 127 indexed inputs, namely the
odds through 253. Two runs use the same task code and reservations:
independent literal tasks, and tasks sharing all checked suffixes.

| stage-14 outcome | independent tasks | suffix memory |
|---|---:|---:|
| actual odd-edge queries | 853 | 234 |
| distinct observed edges | 175 | 175 |
| distinct certified sources, including outside the active interval | 54 | 149 |
| covered input mass | \(501/512\) | \(2079263/2097152\) |
| certificate replay edges, before final format round trips | 478 | 5720 |

Both runs make 32,385 reservations, but skip unused budgets after a task
finishes. The suffix run uses 294 quanta and the independent run 854;
some quanta close an already known endpoint without an edge query.
The exact observed-edge dictionaries agree. Hence the improvement is
proof reuse and avoiding repeat queries, not discovery of new distinct
orbit edges or a demonstrated total runtime speedup. Literal certificate
checking costs more in the simple suffix implementation.

Applying grounded reverse closure to those same 175 observed edges adds
the certificates for 219 and 329, for 151 total certified sources and
mass \(2079583/2097152\). This postprocessing replays 2,640 certificate
edges; it makes no new odd-edge query. It proves 116 of the 127 active
sources; the least still missing is 127. In particular all odds through
125 are directly certified, while the single global residual only
guarantees the first seven indexed inputs, through thirteen. This is a
concrete distinction between the exact ledger and its coarse scalar.

Other controls include 1,023 gamma round trips, 1,024 exact finite-tail
identities, 203 complete receipt-code round trips, the source/proof-weight
hostile, singleton/equality boundaries, and eleven malformed/type/root
padding controls. These are finite exact checks of the implementation.
The infinite normalization, fairness, and implication (5) have the
separate proofs above.
