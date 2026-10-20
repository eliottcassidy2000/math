# Independent review of the six end-atom theorem and inventory connection

**AUDIT PASS.** Reviewed `five_atoms.md`, `five_atoms_verify.py`, and
`inventory.md`. No correctness repair is required. This review covers
the written all-order mechanisms and independently reconstructs their
small controls; it does not promote an interior-context theorem.

Reviewed SHA256 values:

```text
five_atoms.md
b38baa0616f87ba5dc21f5b0fe3e067015df67639499407904c524c36dbbd4d9
five_atoms_verify.py
4b91e444c07800859a68be2617a11f3b7af64ce3eb1fd759173432ad4846b4a1
inventory.md
f856ed678a230a3015d7fe86176ea13ee3e2db8bee20b2f1a015ce6cc6c56fca
```

## Six-atom coverage and the finite base

An independent enumeration of all1024 labelled tournaments of order5,
using reachability from every vertex and a brute minimum over all120
permutations for canonicalization, gives544 strong labelled tournaments,
50 strong standard-Hamiltonian-path presentations, and six ordinary
isomorphism classes. The six atom rows in `five_atom_probe.json` cover
these classes exactly. No canonicalization or strongness helper from
the main scripts was imported in this check.

For a source context with its Hamiltonian path relabelled `(0,...,m-1)`,
the free pairs are precisely `1<=i<j<m` with `j>=i+2`. This yields
`1,1,1,2,8,64` complete presentations for m1..6, totaling77. An independent
path expander, adjacency toggle, and positional anti-reversal predicate
replayed all six templates on all77 contexts. Each template restored A
literally and produced a local tournament with degree multiset
`{4,3,1,1,1}`, hence exactly `1->1->C3`.

The sole blackness failure was exactly the recorded one:

```text
template0, A=(30,12,24,16,2), P=(0,1,2,3,4),
blue flags=(false,true,false,false).
```

The replacement path `(0,2,3,4,1)` is directed in that same labelled A,
makes all four steps black, and restores the same literal A. Thus the
theorem's path-choice quantifier is correct; it does not falsely claim
that every A path works.

## Uniform proof: no free context arc is used

All four paths traverse A contiguously, in alternating path directions.
Every internal path edge is toggled four times; every other internal
context edge is untouched. Cross arcs depend only on the first and last
path vertices, and their labelled boundary ledger cancels. Their finite
local targets and port directedness can therefore be checked with two
distinct context ports; m1 is separately in the complete finite base.
No source hypothesis is needed for this literal identity. The source
hypothesis is used for blackness.

For the two middle paths, the written positional image table checks
directly. For templates0/1, the two occurrences of u lie at the end of
the reversed context or beginning of the forward context. Reversal pairs
`(u,a_(m-5))` with `(b4,a_(m-6))`. For template4 it pairs
`(u,a_(m-1))` with `(b4,a_(m-2))`. Templates2/3 at the third move and
template5 at both middle moves give the analogous `(m-3,m-4)` pair.
At m>=7 each source arc is a known true nonpath arc, and each displayed
backward cross arc is false and outside the previously toggled boundary.

The two exceptional local image rows also hold directly:

* Template2, second move: `(u,b3)` maps in reversed order to `(b2,b0)`.
  The first is an untouched forward cross arc; the second is absent in
  the initial atom and is not toggled by the first path.
* Template3, second move: `(u,b1)` maps to `(b4,b3)`, with the same true /
  false ledger.

Thus every proposed anti-reversal would equate a known arc with a known
nonarc. This is a uniform proof, not extrapolation from large transitive
contexts. As an additional hostile audit, an independent symbolic matrix
left **all other context arcs unset** and replayed every path for m7..100.
All1128 middle-step index/orientation checks were determined by the known
arcs alone and passed; no comparison consulted an unset entry. This is
stronger than merely evaluating the same formulas in TT_m, while the
all-order conclusion remains based on the index and boundary argument.

The first and fourth moves are black because their corresponding endpoint
tournaments have singleton first SCC and respectively strong-five / C3
last SCC. They are nonselfconverse. Positional anti-reversal persists when
the same path edge is traversed backwards, so nonselfconverseness at the
final endpoint legitimately certifies the fourth move.

## Code audit

`replay` expands the two proxy labels into the chosen context path,
relabels only the context positions, records blue flags before moving,
and compares the final adjacency matrix to the literal restored context
plus the stored local target. Its `move` helper checks path legality.
The `need` checks are explicit exceptions and remain active under `-O`.
The small-base loop checks the entire77-presentation universe per atom,
and gates the only repair on the exact atom, context, and blue-flag vector.

The TT_m checks in `five_atoms_verify.py` are correctly described as
checks of index identities and certificate values. By themselves they
would not prove arbitrary-context blackness, but the written proof
identifies every consulted arc and supplies the missing universal step.
The independent unset-arc audit above confirms that dependency claim.

The persisted fully independent C5 implementation remains available in
`five_atom_audit.py`: it imports none of the main proof helpers and checks
all1100 labelled source contexts and7888 Hamiltonian-path presentations
through m6. Its all-size proof is in `five_atom_audit.md`.

## Inventory theorem and empty-context boundary

The `inventory.md` four-path ledger is correct for every nonempty source
context, including m1 where first and last context vertices coincide.
Its surviving arc-toggle set is exactly `{xz,ab,bc,ca}`; all context and
cross arcs cancel. The blackness score certificates are respectively
`m+6,3,5,m+6` against the required blue sum `n-1=m+5`. The third uses
the internal opposite-position pair b,y, so it still excludes blue at
the endpoint-score exceptional value m3.

The file now correctly excludes the empty context. At m0, merely erasing
the context from the four words makes the third path invalid (its c->y
arc is absent after the first two moves). Moreover, the proposed target
`C3->C3` is pure blue, so no all-black path to that target is possible.
The source `TT3->C3` is active; the text does not mistakenly call it pure.

The inventory-joining corollary correctly retains a nonempty source
context: for adjacent inventories t,t+1 it uses
`r=n-3(t+1)>=1` and `A=1^r C3^(t-1)`. Its endpoints have singleton first
SCC and triangle last SCC, so both are active. Together with the inherited
fixed-inventory theorem, this joins every admissible inventory in the
three-atom language. The conclusion is explicitly about one component
of the **ambient** black metagraph, not an induced graph of terminal
SCC words, and it claims neither a global feedback invariant nor global
black-core connectivity.

No edits outside this new review file were made for this audit.
