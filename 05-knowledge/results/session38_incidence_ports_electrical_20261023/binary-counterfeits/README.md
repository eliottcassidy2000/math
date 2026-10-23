# Every binary quotient has a literal mirrored counterfeit

23 October 2026. **FINITE-EXACT, complete order-eight universe.** All 255
partitions of the fixed nine E33 symbols into two nonempty color blocks
admit a genuine alternative endpoint-port matching satisfying the declared
constraints. The proof is not inferred from class collisions or failed
forest tests. There are 61 distinct literal four-port switches covering
these 255 menus; no minimum witness-cover claim is made.

Together with the independently certified three-color forest menu
(c_head+c_edge) mod 3, this proves that **three colors are necessary and
sufficient among quotients of this fixed nine-symbol alphabet for the
converse-authenticated order-eight wiring problem**. The bound is not
universal over arbitrary invariants or symbol-assignment rules.

## Alphabet and exhaustive scope

For an outgoing edge v→w, the base symbol is

    z = 3(c_w mod 3) + (c_vw mod 3),  in {0,…,8},

where c_w counts cyclic triangles through w and c_vw counts cyclic
triangles through the unordered edge. A vertex retains the multiset of its
outgoing symbols, and the tournament retains the multiset of these rows.
The row boundaries and direction are retained.

A binary partition is represented by a nine-digit string starting with 0.
The remaining eight digits range over all nonzero binary words, so both
color blocks are nonempty. Complementing the color names gives the same
unlabelled partition; requiring the first digit to be 0 removes exactly
that duplication. Thus there are exactly 2^8−1=255 cases, and the checker
verifies their exact set, not merely their count.

## What each witness preserves

For each menu the receipt gives two actual closed source masks m1,m2.
Exchange their open partners m1^F,m2^F and do the reflected exchange,
where F=2^21−1. All other literal ports keep their native partners.
The four source masks are distinct and black. Their four target masks are
permuted bijectively.

Every certificate independently verifies:

- the source and target are literal closed/open ports in the complete
  mask-to-class ledger;
- every target retains the actual recolored signature required by its
  particular source, both before and after reflection;
- endpoint-score compatibility and the complete target-port inventory;
- commutation of the alternate matching with path reflection;
- nonnegative ordinary transport entries and full class-level converse
  reflection;
- unchanged blue wiring and unchanged loop multiplicities;
- even literal black multiplicities, including reflection-fixed cells;
- the pointwise cyclic-triangle transport law and the entire merged
  c3-current vector;
- nonzero changes to both directed and undirected merged wiring.

Each witness uses four distinct merged nodes and changes four merged
black multiplicities by ±2. Hence these are actual counterexamples to
uniqueness in the stated relaxed port model, even with parity, electrical
current and converse information retained.

They are not different complement-realizable tournament wirings. Indeed,
each altered source is deliberately no longer paired with its actual
bitwise complement. The missing information is literal partner ownership,
not the listed margins or scalar current.

## Why this establishes the lower bound

A successful two-color menu would be one of the 255 enumerated partitions,
but each has a verified counterfeit. One color is a coarsening of each
such menu, so the same witness also excludes it. This establishes the
necessary lower bound of three colors without assuming that every failed
forest criterion means a different feasible transport.

The separate complete quotient census supplies a positive three-color
example with partition `012120201`, equivalently
(c_w+c_vw) mod 3. Its ordinary/converse signature-class graph is a forest,
so its two source-margin systems uniquely force every ordinary target-class
count by leaf elimination. The positive proof and the literal negative
witnesses concern the same converse-authenticated model.

The three-color success is finite at order eight. It does not contradict
the exact-E collisions or the all-size finite-alphabet obstruction.

## Artifacts and reproduction

- `candidates.json`: the 255 partition strings and deterministic source/
  target class choices used to locate the witnesses. These choices are
  not trusted by the verifier.
- `check.py`: recomputes all 6,880 classes' E33 symbols directly from their
  canonical adjacency words, finds actual donor masks, and checks every
  condition above without importing another research producer.
- `receipt.json`: all 255 literal matchings, exact ordinary cell changes,
  directed and undirected merged changes, 61-switch deduplication count,
  check counts and input hashes.
- `check.normal.out` and `check.optimized.out`: frozen replay logs.

Run with the complete wiring packet as input:

    python -B check.py --wiring PATH
    python -O -B check.py --wiring PATH

PATH must contain `classes.tsv.gz`, `ports.tsv.gz`, `edges.tsv.gz`, and
`class_of_mask.u16le.gz`. The default workspace path is the sibling
`session38-wiring`. No raw binary, external solver, or exploratory probe
script is needed. Default runs compare the frozen receipt; `--write`
intentionally regenerates it. Every check survives optimized Python.
