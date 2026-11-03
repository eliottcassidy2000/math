# Verification and proof review

**Scope:** the declared constructions, operator identities and conditional
compiler. Full spectrum coverage, all native metagraph components and a
universal target-selector are not certified by this packet.

## Exact computational controls

| Job | Universe and independent path |
|---|---|
| `junction.py` |13 literal wrapper sources,126 rooted singular seals,206 complete native/common-automorphism cases through order4; independent permutation counts at orders8/9; named17-query reset, slot-divisor contrast, H31 cover tuning and ten controller entry classes |
| `gluing/check.py` |All74 labelled sources of orders2–4; complete marked covers before/after G versus a separate Hamilton DP;126 abstract multiplicity controls,312 operator entries, exact weighted gcd and prime escape; modular CA and mod49 jets checked against literal polynomials |
| `axiom/counter_check.py` |All16 counter pairs with0<=a,b<=3; complete marked covers at a+b<=3; literal instruction states, native endpoint edits, source strongness, priority hostile, address-map and inherited core-budget controls |
| `machines/check.py` |Input0–31 for register/CA/Turing/marked-source agreement; literal big-integer FRACTRAN and independent Hamilton DP for0–15; fixed program, named native exits, graph-word recovery and blue/black criterion |
| `machines/hankel.py` |All31 strong word sources behind the specified7-by-7 matrix, orders5–13; independent subset DP and odd-cycle collection enumeration; exact rational determinant |
| `incoming_checks.py` |The incoming15-vertex input41 priority gadget; actual91-chord complement identity and its two reflected unequal chords certify a black line |

Every verification gate uses explicit runtime checks that remain enabled
under Python optimization. `verify_packet.py` reruns each frozen job normally
and under `-O`, checks byte equality against its receipt, and records hashes
in `verification.json`. The machine Hankel job intentionally uses its local
`check.py`; the other jobs import no repository implementation. The source
pin file records the inherited mathematical inputs. None of these controls
is an all-order scalar-spectrum census.

## Proof review

- Root independently reviewed the gluing kernel bijection, uniform-context
  divisor, unitriangular weighted-gcd proof, prime escape, CA and jet formulas.
- The gluing agent reviewed the counter carrier's full-polynomial collision,
  strongness, masks, zero-test distinction and inherited native compiler cost.
- The gluing and axiom agents independently derived the singular seal
  identity and checked the five-vertex native junction. The gluing review
  found a genuine draft error: the mod3/gcd guard did not imply eight-state
  controller entry. The k2 ratio10 hostile and the bounded ratio0 bridge
  repair it; J4 and the exact checker retain that guard explicitly.
- Root checked the conditional coverage theorem against the actual inherited
  finite core bound and single-arc compiler, including strongness of prime/
  7p targets, physical-name preservation and temporary-interface scope.
- The axiom agent's completed [independent review](axiom/REVIEW.md) covers
  the fixed machine phases, source programs, colours, bounds and Hankel
  scope, as well as the gluing gcd and mod49 jets. Its one machine indexing
  clarification was applied; code and receipts were unchanged.
- Incoming work through `ab860360b8` was read as mathematical input and
  reconciled in [INCOMING](INCOMING.md). A demonstrated colour error in the
  priority scanner was repaired in the source note, router and output
  labels. Its arithmetic counts are unchanged; the new one-example
  correction job is separate from the five research jobs.

No Lean formalization or proof-assistant certification is claimed.
The finite machine simulators have explicit safety caps; mathematical
all-input termination follows from the stated decreasing-register proof,
not from a finite run cap. The coverage axiom is never used to justify
the unconditional formulas or the computed witnesses.

## Filing checks

The repository documentation checker and direct local Markdown-link
validation passed; all34 direct inherited source pins were checked against
published revisions, with the one repaired prose source explicitly marked.
Whitespace checks passed. `manifest.json` records final packet files other than itself;
`sources.json` records direct inherited source files, not an entire
transitive dependency graph. A successful push is required before the
session is reported complete.
