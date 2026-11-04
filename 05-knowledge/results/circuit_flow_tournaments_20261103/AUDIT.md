# Scope and reproduction

The all-size claims CF1-CF7 have explicit combinatorial/algebraic proofs
in PROOFS.md. CF8 is a complete finite-state induction proof: the full
reachable mod8 and mod27 evaluation vectors, transitions and reaching
words are in clocks.json. The minimized19-state CRT output observer is
proved minimal because its19 current outputs are distinct and reachable.
This does not assert that the full arithmetic carrier has19 states.

Exact controls:

| Universe | Check |
|---|---|
| All64 labelled four-vertex seeds, each of4 roots, five modules | 1280 compiler cases; slots1,2,3 each |
| All75 nonempty labelled inputs of orders1-4, both actual gates | 150 full-profile cases; slots1,2,3; strong outputs |
| Both five-vertex gates with singleton and two-vertex modules | 8 independent permutation-and-cut counts |
| Five inputs: singleton, transitive2/3/4, cyclic3 | Literal two-gate commutator versus `-72 pc3` |
| Six literal gate words through order13 | Separate recursive Hamilton counter and birth replay |
| Actual H611 pair | Separate Boolean-matrix native replay; all9 root profiles; two-rooted-birth counts6467/6379 |
| Actual12-vertex hostile pair | Full automorphism orders81/27;67 common labelled triples; all12 nonzero profile evaluations equal; cut291 reads2450493/2753631 |
| Mod8 / mod27 unary grammar | Complete banks of5/162 vectors; independent compiler transitions; exact minimized clock |

`flow_operator.py` enumerates exterior words and mandatory gap occupancy;
it does not fit the context operator from numerical samples. The checker
uses interpolation only to evaluate tiny *input* polynomials, whose
degree is already bounded by the input order; slot counts are computed
independently from actual graph masks. The final audit uses a separately
written recursive path counter and Boolean-matrix native complement.

All required conditions use explicit exceptions, not removable assertions.
`reproduce.py` runs check.py, clocks.py, native_reset.py, audit.py both
normally and optimized and compares receipt bytes. Hashes are recorded in
reproduction.json. No external datasets or network access are needed.

`clock_explore.py` is a bounded discovery tool, not a completeness oracle:
the initial searches modulo11/13 hit their20000-state cap; searches modulo81
and648 hit a2500-state cap. No classification is inferred from those
truncations. The proof uses only the fully closed8/27 banks.

Validity review: CF4 requires nonempty input; the empty-input polynomial
has a nonzero zeroth value. Arithmetic costs exclude initial-profile and
bit costs. Re-marking is not free. All named source rows and maps remain
available. The neutral reset changes root responses and need not preserve
F49/Q4 phases. Its explicit birth/native clocks and common root frames
are retained; no arithmetic clock is relabelled as physical event time.

No new global A000568 bound, global spectrum coverage bound, or ready return
at35207 is claimed. The inherited failed55->111 and previous source-return
receipts are pinned by hash in audit.json, not reinterpreted as unrestricted
impossibility statements.
