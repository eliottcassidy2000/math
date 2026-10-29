# Operation-closed synthesis: a recurrence becomes a falsifiable grammar

Session `operation_closed_synthesis_20261029`, starting revision
`4f53ff2859d8c79ec417d441b74777d9e6f3b8bc`.

**Status.** PROVED elementary interface statements and regular-language
enumeration; CITED classical unique-Hamiltonian-cycle construction;
FINITE-EXACT tournament determinant criterion through order 12. The extension
of that criterion to every order is OPEN. LRC(14), general Collatz convergence,
and the general Hamiltonian-cycle/torsion conjecture remain OPEN. All audits
here are same-agent checks using independently implemented algorithms, not an
independent referee report. No literature-first claim is made.

## What the synthesis actually produces

The most productive common principle is **equivalence under the next legal
operation**, rather than equality of present invariants. It explains why
useful recurrences appear after adding one missing state bit, and why a small
formula need not give a small closed state space. But this principle has two
additional requirements: keep the native integer lattice when torsion matters,
and use a positive functional that actually sees the target sources when
existence or convergence matters. These are distinct obligations, not a single
universal entropy or rank.

The application is a direct, small generator for the unique-Hamiltonian-cycle
tournaments in [HYP-9445](../../hypotheses/HYP-9445-run-observer-torsion-and-hamiltonian-cycles-star-law-fermat-routes-and-gaussian-gf-questions.md).
The classical construction replaces a whole tournament census by two or three
legal insertions per parent. Testing its two-run observers gave:

| Order | All unique-cycle classes | Unimodular two-run observers |
|---:|---:|---:|
|4|1|1|
|5|3|3|
|6|8|7|
|7|21|16|
|8|55|36|
|9|144|81|
|10|377|182|
|11|987|409|
|12|2584|919|

Orders10--12 are the new extensions of the inherited order 9 census. In all
4,180 classes, the determinant is 0 or+/-1. This does **not** classify the
Smith torsion of singular observers. In particular the first singular case
at order 6 is torsion-free; singular and torsion-bearing are not synonyms.

More informative than the three new numbers is the following exact finite
observation. Give the five classical insertion ports the names a,b,c,d,e
defined in [THEORY.md](THEORY.md). Then, throughout the census,

    det M2(T_w) != 0  iff  w contains neither be nor ee.       (E)

This is a precise all-order conjecture, with a next proof obligation: transport
the determinant or its native kernel through each insertion. Counting more
orders alone will not prove it.

## A candidate explanation, with its exact boundary

Pattern avoidance itself is an all-order regular language. Its three productive
states have transition-count matrix

    C = [[1,0,1], [1,1,1], [1,1,0]],

in the order A,B0,B1, starting in B0. If v_l counts admissible avoiding words of
length l, then

    sum v_l x^l = (1+x)/(1-2x-x^2+x^3),
    v_0,v_1,v_2 = 1,3,7,
    v_l = 2v_(l-1)+v_(l-2)-v_(l-3).

**PROVED for this explicitly defined tournament subfamily** via its injective
growth words; **CONDITIONAL on(E)** for the unimodular class count at arbitrary
order. It predicts2,065 at order 13, among6,765 unique-cycle classes. That next
count is not computed here.

There are three different correct minima for this model:

| Observation task | Minimum |
|---|---:|
|Scalar number of avoiding words by length, stationary linear realization|3 dimensions|
|Accept/reject the avoiding language, ignoring why a word fails|4 deterministic states|
|Distinguish illegal words, legal singular candidates, and accepted words|6 deterministic states|

The minima are proved for the language model, not claimed as all-order minima
for tournament observers. One extra live bit distinguishes B0 from B1. Its
necessity is already a literal tournament fact at order 5: the sources encoded
by `1000100111` and `1000100101` have typeB and det M2=1. After the same e-port,
the former still has determinant1 and the latter determinant0. Keeping the
current type and determinant is therefore not a lawful recurrence state.

This is the same *failure mechanism* as the inherited denominator101 pair
whose common appended letter produces35 versus175. It is not a bijection
between the tournament matrices and the arithmetic cycles.

## What transfers to the open problems

**A000568.** Finite-state structure is useful on a genuine subclass. It cannot
produce all tournament counts by one fixed finite stationary linear transfer
matrix: those counts are at least `2^(n(n-1)/2)/n!`, which grows faster than
exponentially. This elementary obstruction explains why the full cycle-index
and Gaussian-field recurrences need growing indices or nonstationary weights.
It does not rule out faster algorithms or more general concise recurrences.

**LRC(14).** The old tight-AP obstruction distinguishes pointwise information
from information only up to measure zero. We give an exact endpoint atomic
functional, including owner multiplicities, which detects closed-threshold
safety even at an isolated point. For speeds1,...,13 at threshold1/14 there
are 182 owner-tagged endpoint candidates, six distinct safe times and12 safe
owner tags. The integer182 is also the order 10 tournament count above, but
these two derivations provide no mathematical identification of the objects.
The first-good-period27 hostile also passes the endpoint test:646 distinct
safe endpoint locations,736 safe owner tags. No uniform positivity proof for
arbitrary13-speed rows follows from this exact reformulation.

**Collatz.** The older Elias-weighted source measure assigns positive mass to
every positive odd integer. Under this measure the nonroot `5x+1` cycle
`13 -> 33 -> 83 -> 13` leaves survivor mass at least `69/2048` at every horizon.
This gives a concrete rejection test for a proposed convergence argument:
reparametrizing the word or proving generic cylinder contraction cannot erase
that atomic mass. For `3x+1`, proving survivor mass tends to zero would prove
convergence; no such decay theorem is obtained here.

**Integral tournament composition.** The T5 four spurious Z/6 factors and the
prime-power clock carry show why characteristic-zero decomposition, reduction
modp, and exact integer composition answer different questions. The new
determinant calculations use the native integer matrix as a control. They do
not infer Smith forms by concatenating rational sign sectors.

The full proofs, operation definitions and failure boundaries are in
[THEORY.md](THEORY.md). [RECOVERY.md](RECOVERY.md) gives the source/target/map/loss
ledger and the older work that generated the synthesis.

## Verification and next decisive tests

Reproduction from the repository root:

```text
python -X utf8 -B 05-knowledge/results/operation_closed_synthesis_20261029/tournament_probe.py --max-n 12
python -X utf8 -B 05-knowledge/results/operation_closed_synthesis_20261029/grammar_probe.py
python -X utf8 -B 05-knowledge/results/operation_closed_synthesis_20261029/audit.py
```

Python packages: NumPy, Numba, SymPy, python-flint. All integer path counts
are bounded by `12! * 11 < 2^63`; FLINT determinants use arbitrary precision.
There are no removable assert statements in the new proof gates.

The audit compares120 complete native matrices from independent prepend and
append recursions, tests literal permutation cutting on all classes through
order 6, checks native determinants against the sign-sector product, and checks
converse closure on all 4,180 classes. Normal and optimized replay receipts
are recorded in the manifest. Source snapshots have hashes and exact Git
revisions; incoming work under audit is not silently used as a proved input.

Next tests, in order of value:

1. Prove or refute(E) using an insertion-compatible kernel/Schur-complement
   invariant. Minimal bad seeds are `ee` and `cbe`; explain why every later
   allowed insertion preserves their degeneracy.
2. Prove that avoiding words keep both rational sectors unimodular, with
   integral row/column operations. Determinant agreement alone is not that
   proof. The three live states provide a candidate induction interface.
3. For the general HYP-9445 torsion question, locate the first obstruction
   beyond this unique-cycle family. The present census says nothing new about
   all strong tournaments with two or more Hamiltonian cycles.
4. In LRC, use the endpoint functional to formulate a source-located positive
   certificate for one of the actual open entry strata. Do not replace that
   obligation by full algebra rank, a topology computation, or positive mean.

