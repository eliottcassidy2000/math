# Validation and concurrent integration

Status: **PROVED** for the scoped mechanisms in the linked notes;
**FINITE-EXACT** for the universes below. No floating-point tolerance is
used. The controls support, and do not replace, the quantified arguments.
Unrestricted black-core connectivity and a small aggregate A000568
recurrence remain **OPEN**.

## Reproduction

From the repository root, with Python, NetworkX and SymPy available:

```text
python 05-knowledge/results/deficiency_three_20261020/verify.py --replay
python agents/check_docs.py
```

The verifier checks the normalized SHA256 manifest, rejects removable
`assert` statements in the five experiment sources, and runs every source
both normally and with `-O`. Each stdout must equal its committed `.out`
transcript after CRLF-to-LF normalization. The manifest covers this
packet's notes, sources and experiment outputs, plus their direct inherited
code and census inputs. It excludes itself and `verify.out`. Regenerating
the manifest is an explicit maintenance operation (`--record`), not part
of ordinary verification. The verifier is an integrity/replay check, not
a formal proof assistant.

## Exact universes and controls

| Source | Universe and decisive controls |
|---|---|
| [ownership_check.py](ownership_check.py) | All 1099 labelled tournaments of orders one through five for the clone predicate; complete deficiency-zero through three menus at feedback budgets one, two and three, compared bijectively with independent canonicalization. Full compiler comparison at budgets one and two. Shared endpoints, light-to-heavy mergers, changed endpoint/gap roles, and missing external port information are hostile controls. |
| [black_check.py](black_check.py) | 5505 R4-endcap certificates over all fixed-path middles through order six; 960 guarded larger-ear certificates; inherited nonstrong blue sources at orders eight and nine. Every step tests direction, blue/black colour, and final labelled matrix. F5 is a failed boundary-arc guard, not a connectivity counterexample. |
| [ear_probe.py](ear_probe.py) | All 17280 first/last black-path pairs for the ten-vertex F5 target, yielding 18 detours; 6606 full-tiling endcap/middle controls (endcap orders four through nine, all fixed-path middles through order six); 822 general guarded endcaps at orders four through six with middles of order zero through two. The directed triangle and regular pentagon test the order and rotational-symmetry boundaries. |
| [boundary_check.py](boundary_check.py) | All 2187 ternary boundary words with the parity filter stated explicitly, all 105 junction relabellings, 19 connected simple-interior five/seven-vertex poles and a connected 17-vertex control. Independent colouring and near-perfect-matching engines agree. Nonphysical positive tensors, connected golden-moment collisions, and crossing/ternary detectors test the scope of each reduction. |
| [centers_check.py](centers_check.py) | Every ordinary census class of orders three through nine, direct strong-connectivity filtering, and separate Hamiltonian/automorphism replay of every deficit-two survivor; 89 padding controls and a nonrigid-half hostile. This is not the all-order center proof. |

Full command outputs are the adjacent `.out` files. No sampled population
is described as exhaustive outside these universes.

## Proof audits

The root independently checked the ownership gap bound, the exact clone
predicate and ownership quotient, the five R4 detour ledgers, and the
uniform full-tiling endpoint charges. The boundary certificate received
a separate [independent audit](boundary_audit.md), reconstructing its
support cover without importing the producer and checking the physical
Kempe-bijection and matching-factor arguments.

The general endcap detour was also independently audited: the closing
arc and chord supply directedness, the two final endpoint charges force
blackness, and the first move is blue exactly at the rotation-automorphism
boundary. The regular pentagon demonstrates that this last hypothesis is
essential. Any finite checks described only in that audit are supporting
controls; the retained executable universes are listed above.

An invalid draft inference that transferred a leaf's degree constraint to
its neighbor was removed before publication. Its mechanism and repaired
scope are recorded in [centers.md](centers.md) and the mistakes ledger.

## Incoming mathematical signal

The isolated worktree began this phase at `7d21b9fdcf`, incorporating
the all-order black descent and cyclic-five Petersen reductions. During
the session the following additional commits arrived on `origin/main`:

* `2e8ff23dc3`: odd blue spectrum pairing and hollow toggle cubes;
* `71f9f0bcba`: paired-black leaves and the center reduction;
* `b3d59a8495`: the complete deficit-two center classification and exact
  leaf series;
* `760a29ae33`: sharper prime-batch width bounds, the exact synchronized
  triangle sector, and the remaining four-endpoint recognition grammar.

The [incoming center packet](../deficit_two_centers_20261020/README.md)
proves `D(z)=z^4` by an all-order reduction to a finite necessary symmetry
space and complete certificates. Its root README's audited status
supersedes historical pending-audit text in frozen manuscripts. The root
read the cutoff argument, hand local-configuration supplement, complete
census proof, and independent full-census audit. This session does not
claim to rerun the incoming 2121984-certificate census. Our small-order
checker and divisor decomposition are separately identified in
[centers.md](centers.md).

The last packet's README and connection contracts were read before
integration. Its four-endpoint grammar supplies a new target for local
ownership, subject to retaining triangle atom types and relative cyclic
orientation. Transitive-chain weights cannot silently replace those
triangle modules. The packet has an explicit same-agent audit and does
not claim an independent-agent audit; this session does not alter that
status or rerun its experiments.

The targeted message-intake command was attempted during startup but
could not identify an agent because the local `.machine-id` registration
was absent. No identity or registration was fabricated. Current
mathematics was recovered from navigation, proved sources, exact outputs,
and incoming commits rather than from a chronological log.
