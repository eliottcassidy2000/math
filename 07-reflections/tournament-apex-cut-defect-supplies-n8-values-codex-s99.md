# Tournament Apex Cut Defect Supplies n=8 Values -- Codex S99

The user's clue was right, but it needed one refinement.  Raw cut size is not
the invariant.  The useful invariant is the insertion-weight distribution that
the apex cut induces on old Hamiltonian paths.

In the fixed-path model, add vertex `7` after an old rooted `n=7` tournament.
The edge `6 -> 7` is fixed, and the six free apex-row bits decide which of
`0..5` the apex beats.  A Hamiltonian path of the old tournament can accept the
apex before its first vertex, after its last vertex, or between consecutive
vertices `a,b` when `a -> 7 -> b`.  Summing those legal insertion slots over
old paths gives the new `H` exactly.

The exact audit covers all `2^21` rooted `n=8` rows by this apex-extension
description.  Balanced raw `w=3` cuts are broad suppliers: among strong
extensions they produce `463` distinct `H` values up to max `513`, and in
`[1,100]` they miss only `1,5,7,20,21`.  At this rooted level `w=3` does not
uniquely explain `49` and `75`; those values also appear under the
single-unbalanced `w=1` row.  This refines, rather than refutes, HYP-2879:
after quotienting to non-isomorphic strong-ear targets, balanced `w=3` misses
exactly `49,75` and boundary `w=1` fills them.  The rooted ledger explains what
raw cut cardinality forgets.

The sharper signal is the defect form:

```text
49 = 2*25 - 1,  with insertion distribution {1:1, 2:24};
75 = 2*39 - 3,  with insertion distribution {1:3, 2:36}.
```

Both examples use strong `n=7` bases.  So a single apex win almost doubles the
old atom and subtracts a small odd defect.  This fits the S530 tiling model:
the apex tile is the source-sink boundary closer, not an ordinary interior
diagonal.  It should act like a boundary defect term.

The KPS S31f pull adds the necessary guardrail: the literal longest tile
`(n-1,0)` is not required for `H=49` or `H=75` at `n=7`.  The "apex" object
here is therefore the boundary-row insertion defect, not a single named tile
in every fixed-base tiling.  The same pull confirms the 7-adic atom side:
`49=7^2` and `75` are irreducible `m=7` atoms, while `7` and `21` stay
forbidden.

Balanced examples for the same values are less diagnostic:

```text
49 via w=3: {2:4, 3:11, 4:2};
75 via w=3: {2:2, 3:13, 4:8}.
```

Those are centered around insertion weight `3`, as expected for balanced cuts.
They supply many values, but they do not expose the near-doubling law.

This suggests a constructive strong-atom route after HYP-2877 and HYP-2879:

```text
old labelled H-atom
  + apex cut with near-constant insertion profile
  -> c*H - d.
```

For the live `w=1` clue, `c=2` and `d` is a small odd source-sink defect.  That
looks like a real operation, not just a census coincidence.

The assumption challenged here is that balanced/unbalanced cut cardinality is
the proof object.  It is only an address.  The proof object is the path-weight
profile, together with the strong-ear quotient, because that pair preserves the
`H` operation while retaining the boundary exception.

Next I would classify all `c*H-d` apex profiles in the full n=7 to n=8 ledger
and compare the possible odd defects with the E_7 odd-hole data.  If the small
defects line up with the pentagon/heptagon split, this becomes a concrete
bridge between HYP-2877's strong atoms and HYP-2878's apex-prime odd holes.
