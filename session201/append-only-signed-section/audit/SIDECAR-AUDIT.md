# Independent audit of the separately scoped sidecars

**Verdict: PASS under the nontrivial source domain S of the main proof.**
This approval is separate from the main theorem audit. It covers the exact
first-growth gap, the scaled signed-33 append-only hostiles, the three
one-letter collision loci, and the all-order parity exclusion of O/Z
collisions from the formal boundary orbit. It does not cover longer-word
collision classification or any all-target occupancy claim.

Audited text SHA256:
`1fd38bf4651fbf35bf7e7cce033680e22d8e955ed4fb83c62dbbba2d0237254e`.
Audited verifier SHA256:
`20db30b03dbfc3fe4bfbb6a988882cc4505338a32a5b4a70ff82c3b2f0825752`.

## Domain qualification

The phrase “physical source” in the first-growth statement is read in the
main proof's explicitly declared domain: the nontrivial source orbit S of
`t=(1,0,2,1)`. It must not be extended to the optional singleton
`e=(1,1,1,1)` used in the older frozen dependency. At that singleton,
`F(e)=D(e)=1`, so they are not two distinct values. The main proof excludes
that singleton, and every physical source used by the sidecar verifier
belongs to S. Making this inherited domain explicit in the release is
recommended; no mathematical repair is needed within S.

## All-height argument

### First growth and the represented hostile

Before the first appended Z or J, only O letters can occur and D stays
fixed. With no preceding O, the first increase is exactly F(s). With r>=1
preceding O letters, E has been reset to zero and U is `2^r U_s`, giving
`D_s+2^r U_s>F(s)`. Subsequent D values cannot decrease. The empty word and
Z attain `D_s` and `F(s)`, respectively. Since actual states in S have
`U>=2d>=2E` and U positive, `U>E`; hence `F(s)>D_s`.

This proves the exact first two D-tail values at all heights. Their
identification with the first two section targets uses the main theorem
only in the positive-f chamber; it is not claimed for general f-zero
sources, whose deletion and short-completion targets may lie below D.

For `t.ZOOZ=(15,1,23,1)`, direct arithmetic gives f=8 and signed target 33,
strictly between D=15 and F=37. This proves no append-only repair even
though its target exceeds D. The separately displayed representative
`t.OOOOZ=(33,1,63,1)` is correct.

Scaling the initial balanced b by `(OJ)^r` yields the claimed hostile
`3^r(15,1,23,1)`. For the positive representation, `t.OOOOJ` is balanced
with D=U=33, and OJ multiplies a balanced state's D and U by 3. Its E,d
need not scale by that same factor on the first iteration, but the
sidecar does not assert that they do. Thus the given balanced witness
has D=U=`33*3^r`, and W gives the required section witness.

### Collision loci and parity

Solving the four component equalities gives exactly:

- O=J: `D=E=U=d`.
- Z=J: `E=0`, `d=D`, `U=2D`.
- O=Z: `d=0`, `E=U`, with D free.

The independent checker also proves completeness by exact rational matrix
ranks. The O-J and Z-J difference matrices have rank 3 and the respective
claimed one-dimensional kernels. O-Z has rank 2 and kernel basis
`(1,0,0,0)`, `(0,1,1,0)`. These are exact linear identities, not an
extrapolation from a finite lattice census.

The parity conditions `D=d=1 mod 2`, `U=E mod 2` comprise exactly two
residue states. Both are preserved by all three operations, and the
initial z lies in one of them. This proves the invariant for every word
length by induction. Its odd deletion count excludes the O=Z kernel.
The argument is entirely formal and assigns no physical graph to z or
its descendants.

## Independent checks and author replay

`independent_sidecar_check.py` imports only the independent audit
implementation, not the author's code. Normal and optimized runs are
byte-identical PASS results. They include:

- Exact rational rank/nullspace tests of all three collision loci.
- All six transitions from the two parity residue states.
- All 555 actual sources with U<=128, first-growth O prefixes of lengths
  0..10, and direct D-tail inventories through height `8 U_s`.
- Direct first-two-section-target checks for all 307 positive-f sources
  in that source bank, using the reviewed main theorem's complete bound.
- Scaled hostile formulas for r=0..50, with complete no-repair target
  searches for the first three hostiles and independently constructed
  positive section representatives.

Both author sidecar runs were replayed and match both frozen receipts
byte-for-byte, SHA256
`29e0418c6230ea64c8b6ad5de8d555b4c9ba376a8ddb2b53c6cf9ccf32505a96`.
They cover the complete U<=3000 physical bank (59,723 profiles), 62,764
formal boundary states reached through length 11, and 2,401 small-lattice
collision controls. Those finite controls support, but do not replace,
the all-order algebraic and parity proofs above.

Original input bytes are preserved in `original-inputs/`. Their
fingerprints are recorded in `SIDECAR-INPUT-HASHES.json` and
`ORIGINAL-INPUT-HASHES.json`. No original input was edited during this audit.
