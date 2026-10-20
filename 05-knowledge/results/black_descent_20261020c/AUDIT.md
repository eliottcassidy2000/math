# Audit and reproduction ledger

Status: VERIFIED checks supporting the separately stated mathematical
proofs and FINITE-EXACT results. No numerical survival is treated as an
all-order proof.

## Proof review

The cycle agent derived the prime-orbit and even terminal-block proofs.
The root read both in full and checked the boundary signs, sumset
threshold, endpoint inequalities and the exclusion of the last two
positions when choosing the even terminal block. A third agent
independently checked both proofs, including the pigeonhole step
comparing at least delta candidate predecessors with at most delta-1
bad successors. The inherited all-blue even-cycle construction and
repaired non-affine black closed-path theorem remain explicit dependencies.

The root derived the SCC minima, terminal count, basin contraction and
four-step macro. An independent agent checked the symbolic arcs and
every score witness for blackness. That audit corrected the prose arc
ledger: ac and xb change, while xc cancels. The actual target and all
computed certificates already used the correct arcs.

The snark agent supplied the rooted carriers and cut-preservation lemma.
The root checked the support-containment direction, independent branch
path cut counting, small exterior cases and the role of the five
independent port owners. The fixed-gadget boundary predicate is
exhaustively certified, rather than assumed from a graph name.
The root separately replayed the double-star verifier, all 224 closure
counts and the 32,647-subset NN gate.

The root read the 105/17 maps, determinant cancellation, explicit extra
blindspace, normalization conventions and spread/code argument. Exact
spectra, modular rank and explicit null vectors independently certify
the rational ranks. The final spread walk is a bounded hostile probe
inside one symplectic orbit, not a classification of all spreads.

## Exact universes and controls

| Checker | Universe | Positive / hostile controls |
|---|---|---|
| cycles.py | All nonempty cross-set systems for (p,ell)=(3,3),(3,4),(5,3), every sign choice; 2,823 strong self-converse n8/n9 classes | 279,488 additive certificates, literal constructed cycles, actual triangle recounts; fixed-cycle and equal-endpoint failures |
| scc_normalization.py | All labelled B of orders0..5 and every HP of B; complete class/line tables n3..7 | 7,888 guarded four-black macros; one-singleton last move blue; exact SCC minimum and count comparison |
| snark checker | Named graphs; all one-edge triples up to automorphism with full orbit mass; every cut deletion size0..4; every NN subset of size>=3 | Petersen supports, looped dumbbell, multiedge closures, 100 bad fixed double-star interfaces, positive nonphysical remainder |
| 105_17_probe.py | Entire 105/105/91 carriers; all relevant symmetric-group cycle types; explicit 17-spreads | Exact coefficient kernel, transposition fixed counts15/25, rank84 versus70, physical closure normalization |
| spread_walk_probe.py | 10,000 deterministic transvection-walk visits from one spread | Residual C10 and C4+C6, all full first witnesses and operation words; no universal evenness inference |

The inherited n8/n9 tournament universes were authenticated in the
previous packet. The present class scans do not independently re-enumerate
all labelled tournaments. Non-self-converse cases in the theorem use
the all-black cycle-average argument.

## Commands

From the repository root:

```text
python3 05-knowledge/results/black_descent_20261020c/cycles.py
python3 05-knowledge/results/black_descent_20261020c/scc_normalization.py
python3 04-computation/black_descent_20261020c_snarks.py
python3 05-knowledge/results/black_descent_20261020c/105_17_probe.py
python3 05-knowledge/results/black_descent_20261020c/spread_walk_probe.py
```

Each was also run with `python3 -O`; checks use explicit failures.
The agent's two complete snark transcripts are byte-identical, with
SHA256 `14ea4d7088ee7718537c1dc1411ad4da43ae1fca02cad2c7fc9e26634a983121`.
The two algebra receipts are byte-identical between normal and optimized
execution. Source hashes are recorded within those receipts.

Repository checks: `git diff --check`, the local Markdown-target audit,
and `python3 agents/check_docs.py`. The documentation check concerns the
maintained truth surface; it is not a mathematical theorem prover.

## Remaining limits

All-order black descent is proved; complete black connectivity is not.
The SCC macro has a real two-singleton guard. The full degree polynomial
has inherited exact finite formulas, but its arbitrary-order extremizers
and leading asymptotic factor are not evaluated. Strict seven-port
support containment works at a chosen double-star interface and fails
at some other interfaces. No proof makes a Petersen reduction
unavoidable in every cyclically five-connected snark.

Concurrent commits cb2874b2e7 and c75e354bdd were read for their
boundary-conditioned black hubs, trivalent attachments and partial Kempe
linkages. They strengthen the bridge reduction and retain the missing
simultaneous-disjointness coordinate. Their old open descent statement
is superseded by the new proof, while their connectivity caveat remains.
