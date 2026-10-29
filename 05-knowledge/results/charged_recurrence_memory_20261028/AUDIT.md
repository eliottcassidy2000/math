# Scope, reproduction and proof audit

**Same-agent review with alternate computational routes.** No independent
researcher or Lean audit is claimed. Packet-local proved statements remain
in this result packet; no global canon theorem IDs were reserved.

## Exact universes and controls

| receipt | declared universe | checks |
|---|---|---:|
| charged.json | full charged cycle inventories for C and Hol at p=3,7,11, depths 1,2; literal C_2 permutations at p=3,7; literal base affine actions | 139 |
| jets.json | rational moments through arity p+1, depths 0..6, p=7,11,19; p=7 arity 49 at depths 3,4 | 455 |
| memory.json | affine support p=7,11,19 through depth 4; ternary supports through depth 6; exact Hankel certificates at p=7,11 depth 2 and ternary C depth 3 | 186 |
| renewal.json | every fixed-path presentation of ordinary tournaments at n=3,4,5,6; literal 81/27 pair | 130 |
| frobenius.json | all 1053 specified F_27 semilinear maps, every adjacency; support through depth 3 | 2143 |
| rooks.json | all 74 fixed-path presentations at n=3,4,5; all forbidden-arc subsets; directed triangle hostile | 149 |

There are **3,202 explicit guards** per full replay. Counts describe checks,
not independent theorems or statistical samples. Exhaustive universes are
only those in the table. All arithmetic is integral or rational; modular
rank checks are accompanied by mathematical upper bounds, so a full modular
minor certifies the rational rank. No floating-point tolerance is used.

## Alternative paths

- Charged wreath cycle substitution versus enumeration of all 83,349 actual
  C_2 permutations for p=7, matching each charge/cycle cell, not only totals.
- Fourier jet projection versus the complete charged cycle inventory.
  p=19 additionally checks the composite multiplier group h=9, so the
  formulas are not tested only on prime-order multiplier groups.
- Boolean branching support versus independent coin reachability; positive
  histograms versus modular Hankel minors and explicit annihilating windows.
- Actual complement classification versus the odd-partition Burnside count;
  class presentation mass versus Hamilton dynamic programming divided by Aut.
  Direct hidden-state evolution over Q checks the memory kernel at n=4.
- Cycle-pack removal versus literal injective acyclic arc subsets; final
  path counts versus a separate Hamilton dynamic program.

## Proof checks and boundaries

F1 keeps paired-copy multiplicity separate from cycle length. F2 uses
character orthogonality only after identifying uniform child products.
F3 proves positivity of the surviving blocks explicitly; it does not repeat
THM-1770's invalid atomwise noncancellation argument. Exponents E_j in F3
are integers: 1+p+...+p^(j-1) has parity j.

M1's minimal order depends on the r=0 convention 0^0=1. Omitting arity zero
removes a zero-support root. M2 excludes p=3, whose root fixed support lacks
1 and is handled separately. S2's interval proof supplies the all-depth
length bound; finite coin checks alone would not establish that bound.

N1 preserves the original presentation and path throughout the complement
word. M^t instead describes resampling. N3's set carrier and linear span
are distinct mathematical objects; only complement, not re-marking, is in
their operation menu. All census classes are ordinary, not converse-merged.

F_27 is represented by an irreducible cubic over F_3. The displayed
semilinear group is checked as an automorphism subgroup; no undocumented
full-Aut result is used. Its nonabelian phase rule is explicitly outside the
prime-field multiplier-kernel theorem. Binary normal-basis tests are not
silently substituted for this odd-field action.

The 81/27 pair survives; the cycle-index transporter hostile survives;
the directed triangle rejects the unguarded rook formula. No compression
here supersedes the inherited all-context ordered-run reconstruction or
claims closure for arbitrary general strong atoms.

## Provenance and runnable entry points

The startup and archaeology base was d3741fe6ed02. RECOVERY records which
older mechanisms were selected, including their correction status and
exact file paths. External classical sources are linked at their claims.
The seal pins packet files, imported repo Python dependencies and principal
proof inputs using LF-normalized SHA256, so Windows line endings do not
invalidate the seal on another host. It excludes only seal.json itself.

`replay.py` recomputes all six receipts without rewriting them. Its JSON
normalization removes Python tuple/list display differences, not data.
`python -O` leaves every explicit RuntimeError guard active. `seal.py`
checks hashes; `seal.py --write` deliberately reseals after an intended edit.
