# Proof and computation audit

2026-10-25. Internal adversarial proof review plus separate exact algorithms.
No independent human/agent audit; no formal theorem-prover build.

## Proof gates

- **Types:** base coefficients, character indices, root coordinates and
  actual rational points are separate objects. The conductor is an ideal
  of both B and A, principal on A; no claim that it is principal on B.
- **Hypotheses:** B1 uses a UFD with distinct prime radicands times units.
  Domain/normality is invoked explicitly for normalization. Hourglass
  normality is inherited from the identified root-cover theorem. The
  geometric statements are characteristic zero. No characteristic-two
  sign representation is smuggled into them.
- **Necessity and sufficiency:** every multiplication row is tested in the
  conductor definition. Each required prime exponent has an attaining
  hostile row. Conversely those exponents make all products land in B.
- **Symmetry and charts:** all seven punctures of all 120 labelled
  complementary Fano pairs are checked. Character labels and twist ratios
  remain fixed; triple-overlap transitions and conductor transition are
  tested separately. A switch of a chosen basis does not relabel a port.
- **Quotient loss:** the tensor algebra keeps eight components per factor.
  Reduced geometric sheet counts are distinguished from rank-64 scheme
  fibres. Unit squareclasses remain after restoring integral monomials.
- **Generalization:** the composite counterexample refutes the support-only
  extension of the prime formula; the corrected C1 theorem keeps the full
  coordinate images. The n=1 coprime example recovers an independently
  recognizable conductor and its exact gap count.
- **Consequence:** no rational-point existence, finite torification of the
  whole surface, universal Collatz convergence or Busy Beaver termination
  is inferred from an algebraic normalization or a finite test.

## Exact computation

All six scripts use explicit exception gates, which remain active with
`python -O`. Each normal/optimized receipt pair is byte-identical.

| Script | Explicit universe / result | Gates |
| --- | --- | ---: |
| check.py | 120 Fano pairs x 7 charts; 800 ordered F_2^4 decompositions; all 22 actual boundary strata; twist and chart checks | 21,817,748 |
| audit_monoid.py | Independent additive-monoid BFS, all 120 pairs, every exponent in [0,2]^6; index valuation 16 at each branch | 87,480 membership checks |
| check_prime.py | Prime and composite decompositions listed in CYCLIC C4; sharp composite hostile and coprime gap counts | 39,643 |
| check_native.py | 576 literal sources, all ten simplified rule branches; blank entry after 162 native steps | 1,719 |
| check_histogram.py | 2,672 inverse histories, sources 1,5,13, depths 1..6, valuations 1..8, moduli 3,9,27,81 | 31,963 |
| check_stack.py | 972 complete low-dimensional complement/inertia cases; 23 norm-seven strata including the open; literal source decoder | 2,313,956 |

The monoid audit imports no producer code and generates reachable exponent
vectors by addition before comparing with the predicted residue ideals.
It confirms the displayed algebra in a finite box, not the all-order
proof. Native-machine verification uses literal tape sets and the raw
transition table, not the arithmetic macro to simulate itself.

Positive controls include trivializing the residual quadratic twist,
coordinate decompositions with no shared ramification, empty correction
terms, and the exact blank-tape entry. Hostiles include a lost squareclass,
three independent residual classes, omitted conductor factors, non-surjective
composite coordinate maps, duplicate source owners, and a printed negative
native clock. Their failed implications and repaired scopes are recorded in
MISTAKES.md and the proof files.

Run from repository root:

```sh
python3 -B 05-knowledge/results/normseven_boundary_20261025/check.py
python3 -B 05-knowledge/results/normseven_boundary_20261025/audit_monoid.py
python3 -B 05-knowledge/results/normseven_boundary_20261025/check_prime.py
python3 -B 05-knowledge/results/normseven_boundary_20261025/check_native.py
python3 -B 05-knowledge/results/normseven_boundary_20261025/check_histogram.py
python3 -B 05-knowledge/results/normseven_boundary_20261025/check_stack.py
python3 -B 05-knowledge/results/normseven_boundary_20261025/verify.py
python3 agents/check_docs.py
```

Repeat the first six commands with `-O` to reproduce their optimized
receipts. The packet verifier checks source hashes at pinned Git revisions,
receipt/script consistency, both replay modes, and the manifest. Attached
PDF hashes are recorded; absence of a user's local Download files on
another machine does not prevent checking the repo packet.

`agents/check_docs.py` passed after the route update. The current frontier
entry is a replacement, keeping the bounded-startup budget intact. No
chronological session log was used as a current mathematical authority.

The final sync received f093b00695 and 2b69da9a76. Their mathematical
result files (cancellation, specialization, twist mass, root stack and
Fourier effectivity) were read, and their root-stack result led to the
additional exact 2-product theorem in STACK_GLUE. The frontier conflict
was resolved by retaining both S63 and the new chart route. No unrelated
source text was discarded.

## Remaining limitations

The explicit genus/normality input is inherited, with a pinned proof file;
this session did not formally reprove it. C1's algebraic proof was reviewed
internally, not by an independent agent. The attached papers were read and
selected mechanisms tested, but their full Lean developments were not
downloaded or built. Mazur's large density theorem remains source-cited.
The displayed Antihydra clock discrepancy is an exposition-level witness,
not an audit result about its formal code. Universal fixed-source ROOT
positivity and six-colour forcing remain open.
