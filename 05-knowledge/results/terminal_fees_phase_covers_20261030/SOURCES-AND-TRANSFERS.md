# Algorithmic improvements: exact interfaces before numerical records

**Status:** source-checked research synthesis, with the precise read scope below.
No external submission, posting, Lean build or independent audit of the long
external proofs was performed. External instructions are not session actions.

## Circuit lower bounds: a proved barrier to retuning an unchanged relaxation

The requested [Circuit Breaker page](https://www.yukon.org/circuitbreaker) was
read, together with its public verifier rules and
[baseline proof](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/master/submission/proof.md)
(Git blob1c739b148ffc706654374b337f512df8a490a8e2). The page's stronger evaluated
coefficient is provisional pending human review; it is not imported here as
a proved lower bound. The baseline itself claims no new theorem: it combines
published circuit and affine-disperser results. Those underlying long proofs
were not independently rechecked in this session.

One short exact obstruction *can* be checked locally. The baseline's unchanged
linear program includes

    u <= 3+phi,
    u <= 5-Q,
    u <= (5-2phi+Q)/2,
    c = u-phi/2.

Multiply the three bounds by7/10,1/10,1/5. The Q coefficients cancel, the
constant is31/10, and the phi coefficient is1/2. Therefore c<=31/10.
This excludes an improvement by retuning that LP, not by changing the proof
interface or adding a genuinely valid new case constraint.

That is the exact analogy with the earlier positive-envelope Collatz hostile:
searching harder cannot repair an impossible system of potential inequalities.
The useful operation is to retain another multiplier digit, not to interpret
an unsuccessful optimization as evidence against actual Collatz content.

Incoming Session117 supplies a more concrete local model:
[point-defect circuits](../session117_ramsey_circuit_interfaces_20261030/circuits/POINT_DEFECT_CIRCUITS.md)
have exact B2 costs controlled by support and branching contexts. The proof
uses two different downstream contexts to expose information a one-bit
intermediate cannot retain. Its explicit upper bound also stops this small
construction from being advertised as a large circuit lower-bound advance.
For our applications, the corresponding test is whether two states identified
by a quotient can have different terminal fees or attachment verdicts.

## Big O Speedrun: preserve the computational model and the input quantifiers

The supplied `/like` page returned404. The homepage, about page, problem pages,
[verifier README](https://github.com/bigospeedrun/verifier/blob/main/README.md),
[target guide](https://github.com/bigospeedrun/verifier/blob/main/docs/targets.md),
and complete [3SUM/APSP statement](https://github.com/bigospeedrun/verifier/blob/main/challenges/3sum-apsp/EndStatement.lean)
were inspected instead. No record on that site is used as a mathematical
dependency here, and no claim of practical speed follows from a proof label.

The fixed statement gives one program and a bit-width bound for **all** inputs
of each allowed size, and counts actual word-RAM steps. For APSP it retains
reachability as well as distances; for min-plus products it requires a witness
attaining each reported minimum, as well as the lower inequalities. The file
explicitly says that it defines targets and proves nothing.

This gives two useful tests of our own proposed improvements:

* The phase-cover compiler really eliminates enumeration proportional to the
  *added speed*, but preprocessing the core and exact-integer bit costs remain.
* Faster min-plus multiplication of the25-state excursion kernel does not
  compute its unknown native edge prices. The bottleneck is an all-word price
  certificate, not the tiny final convolution. This prevents applying a major
  algorithmic improvement to a step that was never the proof obstruction.

An extension worth testing is adaptive certificate refinement: first solve a
cheap quotient; split states only when a retained terminal predicate differs;
then validate the resulting weighted transitions independently. Our multiplier
experiment implements one fully specified refinement of this kind. A general
polynomial bound on how many splits suffice is **OPEN**, not implicit in the
procedure. The previous exponential exact-observer lower bounds still apply.

## Quasi-Riemann: retained scale ratios and genuinely replicated targets

Primary source examined:
[OpenAI, The Quasi-Riemann Hypothesis, October5,2026](https://github.com/openai/math/tree/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026),
specifically `build/paper2.tex`, blob2000faddbbebac5de0ecfe0b962534ea61a852d5.
Read scope: introduction, five-step proof outline, canonical induction and its
completion, and the relevant theta/deconvolution formulas. This is a **CITED
mechanism**, not an independent proof audit of the entire zero-free theorem.
The source claims an11/12 zero-free half-plane in its stated setting; that is
not the Riemann hypothesis and is not a theorem proved in this repository by
this session. The separate public QuasiRiemannTracker README was also read;
its verification labels do not replace checking hypotheses and statements.

There are three particularly useful operations in the source.

1. It embeds the desired fixed character sum into a family indexed by u.
   At suitable prime sixth powers u=p^6 the character becomes an explicit
   coprimality indicator, so many family members reproduce the target with a
   controlled error. A mean-square bound can then constrain the fixed target.
   The copies are established by an arithmetic identity: they are not a claim
   that an arbitrary fixed point behaves like an average point.
2. It completes a cubic sum by adding cube indices, then uses Möbius inversion
   to recover the original object. It does not discard terms that could cancel
   the desired summand. This is directly relevant to the repo's signed-source
   and orientation-trace failures: the inverse is part of the representation.
3. Its recursive estimates retain a scale ratio while both scales contract.
   A strict separation parameter bounds the number of induction rounds, and
   derivative norms and small exponent losses are budgeted across those rounds.
   Changing a modulus without retaining this quantitative separation would not
   transport the estimate.

The source and target of each proposed transfer must be separated:

| source operation | our target | preserved predicate | missing obligation |
|---|---|---|---|
| many exact arithmetic copies of one sum | fixed-source ROOT energy | actual marked source has positive mass | construct equally faithful copies with controlled cost |
| cube completion plus Möbius inversion | quotient/controller completion | recover original signed contribution | prove a bounded inverse in the relevant weighted norm |
| joint scale contraction with a retained ratio | changing Collatz precision | original integer clock and paid terminal fees | a uniform separation or summable debt across unbounded depth |
| exact residue law | LRC interval phase covers | the actual added speed and clock | arbitrary-core nonempty coverage, rather than selected residues |

In particular, dividing one scale by5 repeatedly takes an unbounded number
of rounds as the input grows. A constant per-round error would accumulate.
The paper's bounded-depth argument therefore cannot simply be renamed a
Collatz proof. Our first-return potential instead seeks a **telescoping**
payment: terminal fees pay actual new-layer returns with no separate overhead
for each controller change. A uniform proof of those fees is still required.

## Cheap hostile probes for the next general theory

* **Completion/inversion:** append two opposite signed copies. Any rule claiming
  that the completed total bounds the chosen summand without an inverse fails.
* **Scale recursion:** E_j=E_(j-1)+1 has bounded one-step error but unbounded
  depth loss. A theorem must give contraction, bounded depth, or telescoping.
* **Target replication:** identical quotient labels with different exact target
  phases are not replicas. Test the retained multiplier endpoint and LRC point
  components before averaging.
* **Algorithmic acceleration:** count the costs of constructing the kernel,
  reading exact integers and preserving the actual target, not only its final
  small matrix calculation.

These are connections between proof mechanisms, not isomorphisms between
Collatz, circuit complexity, L-functions and lonely runner. Their value is to
specify which extra coordinate or inverse would make a proposed transport valid.
