# LRC(14): a divisor-complete owner-word escape from the fixed clock bank 15..43

Status: **FINITE-EXACT COUNTEREXAMPLE TO A UNIFORM FIXED-BANK EXTENSION**.
This is not an unsafe LRC row and not an arbitrary-entry theorem.

## Signal

The divisor-complete eleven-body

```text
H=(1,5,11,13,17,19,23,37,41,70,72)
```

escapes every fixed owner-word clock `15<=N<=43` in the precise sense relevant
to the dyadic-owner-word closure theorem: for each clock, either its core-safe
packet `A_N(H)` is empty or its literal complementary-tail relation `R_N(H)` is
nonempty.  Thus no clock in this bank meets the sufficient certificate

```text
A_N(H) nonempty and R_N(H)=empty.
```

This is a sharp boundary for extending the through-height-40 census merely by
reusing its same clock bank.  It says nothing adverse about the row itself:
the adaptive clock `N=47` repairs this body, with

```text
A_47(H)=(12,20,21,26,27,35),   R_47(H)=empty.
```

## Exact scope and mechanism

The body has eleven distinct positive speeds and is divisor-complete through
14.  At every clock, the audit enumerates the labelled packet and all unordered
odd residue pairs modulo `2N`, including repeated classes.  For every nonempty
packet at clocks 15..43 it prints a least pair whose two strict literal danger
masks cover both lifts over every packet label.  Empty packets fail the
certificate's nonemptiness antecedent.  All arithmetic is integral, including
the closed body-safe test and strict tail-danger test.

The local search was heuristic only; the displayed witness and all claims about
it are independently rechecked by exhaustive finite enumeration.  Searches at
maximum bounds 41 through 71 found no witness, but that negative search is not
exhaustive and is not a minimality claim.  The exact positive claim is only the
single displayed maximum-72 witness and its clocks 15..43; the first-repair scan
is exact for this witness through clock 500.

## Reproduction

```bash
clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_owner_escape_local_search_owner_escape_arbitrary_entry_20260902.cpp \
  -o /tmp/lrc14_owner_escape
/tmp/lrc14_owner_escape 72 5000000 23

python3 \
  04-computation/lrc14_owner_escape_exact_audit_owner_escape_arbitrary_entry_20260902.py
```

The second command deterministically reproduces
`05-knowledge/results/lrc14_owner_escape_exact_audit_owner_escape_arbitrary_entry_20260902.out`.

## Concept-board update

- **Fixed pool:** fails uniformly once the body reaches this maximum-72 layer.
- **Adaptive clock:** immediately survives; `N=47` closes the witness.
- **AP/divisor suppliers:** divisor completeness alone does not force a
  certificate inside a preassigned bounded bank.
- **Deletion/refactorization:** the small entries force many empty low-clock
  packets, while the sparse high entries arrange coverable packets at the
  remaining clocks.  That is the escape mechanism, not an LRC obstruction.
- **Multiple anchors:** remain a plausible entry route because changing the
  reference changes the body and hence its packet geometry.

No theorem ID is assigned: this is a finite exact hostile for one proof bank,
not a positive closure theorem.
