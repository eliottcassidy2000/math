# Audit, source pins and reproduction

**VERIFIED computational packet; scoped proofs in the linked notes;
Chan–Pak and universal odd-coefficient coverage OPEN.** This is same-agent
mathematical review using independent algorithms, not independent-agent
peer review or formal proof-assistant verification.

## What was checked

The complete ordinary order-ten catalogue has 9,733,056 tournaments.
The Python endpoint engine and an independently implemented C++ engine
agree on its entire Hamilton-count histogram, strong connectivity,
unit-source counts and selected-mask counts. Marked orbit counts use
actual automorphisms, with converse types kept separate. A second C++
pass checks H, defect and U individually for all 478,684 retained strong
marked interfaces on 465,023 source graphs; the U checksum is 19,177,986,462.
This is stronger than checking only aggregate totals.

The complete stable census uses a proved finite run bound. The mixed
certificates use entire tensors, including higher coefficients. Independent
split-cover kernels check both obstruction examples and all 1,379 single
birth-clock slope kernels. Literal Hamilton counts check the small
inflation grids. The complete order-eleven descendant census independently
recovers the unique two-root state predicted by the cylinder intersection.
It does not enumerate order-eleven births.

All 786 selected lower-head interfaces in the completed binary compiler
are independently recounted. Its audit includes 1,572 exact affine checks
and 1,571 whole constructed tournament counts; one larger control is
algebraic. All-length conclusions follow from the constructor identities,
not from these finite checks. The remaining finite head entries use the
explicitly pinned inherited complete order-eleven spectrum. Two newly
needed finite heads have literal adjacency witnesses.

There are 37 ordinary/optimized Python replay pairs. Thirty-six compressed
receipts are byte-identical. The order-nine two-block pair has an identical
complete mathematical payload and only a historical scope spelling change:
“order-nine” versus “order-9”, and “order eleven” versus “order 11”.
The manifest checker permits exactly those replacements. Its own normal
and optimized receipts are also byte-identical. Executable validity gates
use explicit exceptions and survive optimized Python.

The source and dependency manifest is [packet_audit.json.gz](packet_audit.json.gz).
It records 38 source files, 22 loaded repository dependencies, six inherited
data files, 77 output receipts and the inspected figure. Text hashes
normalize CRLF to LF. It also checks local documentation links, Python
syntax and agreement between the major receipt totals. It is a consistency
check on finished computations, not a replacement for rerunning them.

## Pinned inputs and provenance

The public input is McKay's [ordinary tournament catalogue](https://users.cecs.anu.edu.au/~bdm/data/tourn10.txt.gz),
linked from his [directed graph data page](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html).
The compressed file is 37,522,053 bytes, SHA256:

```
6e45e6d53d48938f9e1dd016a6da954a41f707a0ac287223d260b389b9cb3f32
```

The inherited complete spectrum is
[least_order_11.txt](../path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt),
with LF-normalized SHA256:

```
42f2fe3f6ff6d198baab67d526b6e727bd7d62ad51d4b1ebdc0ec34900d5f928
```

The source manifests also pin the inherited order-nine catalogue, earlier
ordinary tournament bank, contraction DAG, hostile cyclic-port receipt
and prime wheel. The logarithmic-core cutoff, prime/two-prime-lane
reduction, minimum-order component recursion and path-tail identities are
credited in [FORBIDDEN.md](FORBIDDEN.md). Coefficient 1039 already had an
order-eight seed; this packet adds a renewable stable family. No improved
general bounded-coverage record is claimed.

## Reproduction

Run from the repository root with Python 3, NumPy, Numba and SymPy.
The independent C++ wrappers require `g++` with C++17 support. The optional
figure uses Matplotlib; it is not a numerical input. The wrappers currently
use the Windows `TEMP` directory. Place the public compressed catalogue at
`$env:TEMP/math-tourn10-mckay.txt.gz` after verifying the hash above.
`census_ten.py` also accepts an explicit catalogue path.

In the following commands, set the packet prefix once in PowerShell:

```powershell
$packet = '05-knowledge/results/joint_stability_forbidden_20261112'
python3 -X utf8 "$packet/audit_packet.py"
```

That first command checks the committed packet quickly. Full reproduction
is substantially more expensive. Run each dependency chain in order:

```powershell
python3 -X utf8 "$packet/census_nine.py" unit
python3 -X utf8 "$packet/census_nine.py" stable
python3 -X utf8 "$packet/joint_nine.py"
python3 -X utf8 "$packet/birth_union.py" 9
python3 -X utf8 "$packet/census_ten.py"
python3 -X utf8 "$packet/audit_catalogue.py"
python3 -X utf8 "$packet/audit_unit_slopes.py"
python3 -X utf8 "$packet/stable_census.py" 10
python3 -X utf8 "$packet/joint_ten.py"
python3 -X utf8 "$packet/birth_union.py" 10
python3 -X utf8 "$packet/descendants_eleven.py"
```

The family and coefficient chains are:

```powershell
python3 -X utf8 "$packet/obstruction.py"
python3 -X utf8 "$packet/contractible_obstruction.py"
python3 -X utf8 "$packet/crossing_family.py"
python3 -X utf8 "$packet/coefficient_band.py"
python3 -X utf8 "$packet/affine_births.py"
python3 -X utf8 "$packet/audit_affine_births.py"
python3 -X utf8 "$packet/bulk_head_interfaces.py"
python3 -X utf8 "$packet/head_target_ten.py"
python3 -X utf8 "$packet/finish_bulk_heads.py"
python3 -X utf8 "$packet/audit_bulk_heads.py" --completed
python3 -X utf8 "$packet/higher_odd_families.py"
python3 -X utf8 "$packet/coefficient_intersections.py"
python3 -X utf8 "$packet/complete_rays.py"
python3 -X utf8 "$packet/natural_recursion.py"
python3 -X utf8 "$packet/library_boundary.py"
```

For a literal graph from either illustrative complete family:

```powershell
python3 -X utf8 "$packet/complete_rays.py" --family crossing --parameter 3 --exponent 2
python3 -X utf8 "$packet/complete_rays.py" --family 1039 --parameter 2 --exponent 1
```

The scoped negative probes are `two_block_probe.py 8` and `9`,
`rank_raise.py 1` and `2`, `rank_raise_exterior.py`,
`reflection_crossings.py 13`, `decorate_obstruction.py` with and without
`--all-marks`, `double_decoration.py`, `extreme_crossings.py`, and
`extreme_ten.py`. Their receipts declare the exact universes. A failed
template search is not a nonexistence theorem outside that template.
The baseline interface audit is `audit_bulk_heads.py` without arguments;
the illustrative discovery search is `head_repair.py`.

Add `-O` immediately after `python3` for optimized replays; the common
writer uses separate `.optimized.json.gz` outputs. The two complete C++
audits provide independent algorithmic verification in place of duplicating
the expensive ordinary-catalogue Python pass. Finally rerun
`audit_packet.py` in both modes, `git diff --check`, and the repository's
`agents/check_docs.py`.

## Remaining boundary

The negative target is pointwise convergence of the forbidden sets to
{7,21}; no finite order has only those two omissions. The packet gives a
finite decision cutoff for each candidate through an inherited theorem,
not a uniform proof that every allowed candidate eventually disappears.
Stable-port acyclicity, strong confluence modulo converse, all-order birth
generation, and uniform lower-head closure for the 478 even higher-mode
axes remain OPEN. The finite family libraries have proved sparse output
images and cannot by themselves establish the universal conjecture.

The ordinary confluence failure and joint-stability failure are now
counterexamples, not open questions. The repaired local criterion requires
whole adjacent-line identities; checking the double-contraction corner
alone is insufficient. See [CROSSINGS.md](CROSSINGS.md).
