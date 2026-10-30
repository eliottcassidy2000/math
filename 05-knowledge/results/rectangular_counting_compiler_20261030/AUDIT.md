# Verification and limitations

2026-10-30. Same-agent review. The formal check and native computations
are different verification paths; no independent human review is claimed.

## Formal result

`check_lean.py` verifies the external and Mathlib revisions, requires a
clean external checkout, invokes Lean on this packet's root module and
checks every requested axiom report. The successful output and exact
source hashes are in `lean_output.txt` and `lean_receipt.json`.

```text
python -X utf8 -B check_lean.py --external <pinned-external-checkout>
```

External revision: `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
Lean:4.35.0. Mathlib:`56a005855ea17e0b0958f73ec801ae140b0a1679`.
Command inside that checkout: `lake env lean <packet>/ReducedEngine.lean`.
Exit0. All eight printed theorems use exactly the standard axioms
`propext`, `Classical.choice`, `Quot.sound`; no `sorryAx` is present.

The imported external engine already had a successful full build in
the prior packet. The present formal root imports EngineGroup and
contains the new adapter, all-depth proofs, boundary-family identity,
and exact leaf recurrence. It is intentionally a packet-local overlay,
not installed in the external project's Solution root or the math
repository's general Lean root. The old Solution theorem has not been
modified to use this producer. The new claims are reachable from the
root command above, not from an assumed unrun monorepo build.

An initial tactic attempt left an equality unsolved and consequently
printed `sorryAx` for dependent declarations. That failed run was not
accepted. Normalizing and closing the equality repaired it; the stored
receipt is from the subsequent successful complete module check.

Hypotheses are satisfiable: K0>=2 permits K0=2 or the actual48;
arbitrary constant clique tables give examples at every declared shape;
all moduli are allowed. No proof comes from an impossible constraint.

## Native controls

From this directory:

```text
python -X utf8 -B verify.py --output certificate.json
python -O -X utf8 -B verify.py --output certificate.optimized.json
```

The outputs match byte-for-byte. Explicit checks raise exceptions and
remain active under `-O`. Total local checks:1,106,254. The inherited
native81/27 helper's internal checks are not included in that total.

- Full ambient enumeration: all65536 binary assignments on a2x2 grid,
  with one actual rectangle and32768 members in each fibre. Three
  signed response coordinates agree after the exact multiplicity factor.
- All dense choices at G=0,1,2 for the listed q, and G=3,q=2,3;
  deterministic samples at G=4,q=2,5. Certificate records each universe,
  sample size and exhaustiveness flag. Dense colex rank and zero-pad
  extension have exact inverse checks. This Python coordinate order is
  its own explicit implementation; we do not assert that Lean's
  phiEquiv gives the same numerical code order.
-72 full/reduced effective-modulus cases, with every count through
  twice the bank lcm, plus inherited dyadic cutoff checks b=0,...,24.
  Their all-size justification is the already-proved CRT theorem.
- Actual uniform depth1 and2 thresholds: exact big integers for n,
  sOf, K/G/width, q, rectangle and unused exponents. Only bounded
  modular powers are evaluated; no enormous branch set or input N
  is materialized. Integer receipts give bit length, low64 bits and
  SHA256 of the unsigned big-endian integer bytes, plus decimal text
  for small enough values. The small odd bank is diagnostic, not a
  sufficient zero bank for these vast n.
- Hostiles: an observer of an erased coordinate; equal totals with
  different boundary vectors; nonuniform fibre multiplicities; original
  mod3 readings losing Boolean counts;81/27 cyclic-triple collision.

The g2,lambda4 witness is labelled generic and is never substituted
for the actual schedule. Its semantic identity comes from the general
rectangle counting formula; it is not a valid EngineGroup K>=2 grid.

## Proof/model audit

- **Types:** the reduced object is a count/table in ZMod m; the formal
  result does not turn it into a running-time estimate. A finite count
  bound is supplied separately when modular zero testing is used.
- **Directions:** rho equality is sufficient for child equality. We do
  not assert the converse or that all rectangle choices are necessary
  on every T. The observer threshold Q>B is an iff on the stated full
  interval; its necessity witness is the positive count Q.
- **Quotients:** only artificial constant-size padding fibres are
  removed here. Genuine object multiplicities remain weighted.
- **Ports/signs:** named boundary values are preserved pointwise;
  negative internal coefficients and the row/column grouping remain.
- **Schedules:** width and exponent depend on size and depth. Adapting
  the schedule to a boundary assignment requires keeping that dependence.
- **Circuit model:** unrestricted shared B2 DAGs, actual primary input
  length, and all-length language quantifiers remain separate from
  table size, one-pass memory and template feasibility.
- **Scope:** external current badges are cited reports. Our successful
  Lean audit is only for the explicitly pinned engine and new overlay.

## Reproduction seal

```text
python -X utf8 -B seal.py
```

`seal.py --write` creates the manifest; `--staged` also checks exact
staged bytes. Packet .gitattributes forces LF, avoiding platform line-end
changes that would invalidate the receipt. No generated Python cache is
part of the packet. Source pins include the inherited proof documents
and runtime hostile helper. The attached PDF and public HTML captures
are identified by hash in SOURCES; private source files are not copied
into the publication packet.
