# Independent audit of the two-height encoding

2026-10-23. Reviewer: port_rigidity subagent. **PASS** after the wording
clarification below and literal position-label guards.

I read README.md and check.py independently. The telescoping identity
d_ij=1+2(h_ij-h_i+1,j) yields the stated suffix capacity condition and
the unique inverse. The boundary complement value 4-d, versus 2-d in
the interior, is necessary because adjacent path arcs stay fixed. The
reflection formula includes converse, so it preserves the directed path.
The commuting involution action has precisely the stated two- and
four-element orbits: a fixed antidiagonal free bit excludes any
complement-reflection fixed point for n>=3. The count formulas concern
line multiplicities, not endpoint support.

The order-four hostile has two individually realizable triangle values
but an impossible shared arc. Its minimality should be stated as the
first **cross-coordinate** capacity obstruction with individually legal
triangles; a bare invalid array already exists at n=3. The author adopted
that clarification. The proof does not require a new global metagraph
determination or connectivity claim.

I reran the final consumer with

    python -B 05-knowledge/results/metagraph_two_height_20261023/check.py

It passed 1,512,370 checks, including the saved receipt comparison, with
LF-normalized receipt SHA256
`36fdc58da3459ddf5d140df6bf11fbdace5154b1422a355c4c610ce6d00748bd`.
The complete fixed-path universes n=3..7 and the malformed value/address
controls match the declared scope. No floating-point tolerance or
assert-only acceptance is involved. The separate root-curl theorem
gives an independent reconstruction perspective; both coordinate maps
retain the full fixed-path state and therefore respect the incoming
all-labelled-flips memory obstruction.
