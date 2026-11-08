# Independent audit: Session228 recursive physical polynomial grammar

**Verdict: PASS for the stated local theorems and finite-exact compilation. No mathematical correction required.** This does not establish Chan–Pak completeness, unrestricted minimum depth/order, or literature priority.

## Scope and independence

Audited the frozen `PROOF.md`, `README.md`, `verify.py`, `CHECK.json`, and `CHECK.optimized.json` in the original `session228-recursive` packet, now bundled byte-for-byte under `inputs/`. Their exact SHA-256 hashes are in both independent receipts and the portable manifest. No producer or repository file was edited. The repository was initially identified by the requester at `1a178a44f13f9602a3be9e84714be91cbfc7d7a6`; it was concurrently advanced by other work during this audit (observed head `9f9bac37c83c839cb5fd3fe2220991232a87b2b2`). The mathematical packet itself is identified by hashes, not that moving checkout.

The standalone standard-library checker `audit.py` does not import or execute the producer. It independently constructs seed profiles using canonical set partitions and endpoint Hamilton recursion, verifies all seed Hamilton counts by permutations, represents complete polynomials by Newton forward differences, compiles both ordered operand directions, and computes the scalar closure by forward product relaxation. It reads the producer receipt only as expected values and literal witnesses. Separate replay files execute the producer unmodified.

## Results reproduced

- All 1,099 labelled seeds through order five produce 15 order/profile pairs.
- All 121 pairs of labelled nonempty inputs through order three were tested for both ordered join and singleton-cycle construction: 242 full literal-profile controls.
- Profile counts for orders 1–13 are `1,1,2,3,8,14,29,63,139,304,489,709,1032`.
- Exactly 577 attainable scalars through 4095. Minimum-depth histogram: 79 at depth 0, 390 at depth 1, 107 at depth 2, 1 at depth 3.
- The first scalars at depths 1, 2, 3 are respectively 17, 241, 4001.
- Sorted scalar/minimum-depth digest: `7c895b5f3a40b78ab702a2d63b0cd43daccc6fe39f4d353006f98c331b44fa50`.
- Direct literal Hamilton counts for the TT4 collar chain are 1,17,241,4001. Both stored depth witnesses are verified independently.
- Both order-nine collision polynomials were verified completely by ten integer evaluations, rather than checking only the three displayed moments. The two successive collars give literal counts 150831 and 155151.
- Literal C3 exterior singleton padding at every order 3–14 retains H=3, including graphs beyond the global cutoff order 13.
- A separate general-seed-bank test includes the physical order-14 seed C3=>TT11, polynomial t^14+2t^12. At target 3 it correctly remains a depth-zero scalar generator despite exceeding N(3)=3.
- Compiling to order seven with and without target pruning gives exactly the same retained profiles and minimum depths at H<=21.
- Independent normal and optimized receipts are byte-identical. Producer normal and optimized replays are each byte-identical to their frozen expected receipts.

## Mathematical proof audit

### Polynomial law

The single-colour layout A_pre,B_pre,z,A_post,B_post exhausts legal possibilities, since no B-to-A transition can avoid z. The two forbidden events and their intersection have exactly the asserted bin counts. Their inclusion-exclusion formula simplifies to (1). The argument works even for j=1: terms involving zero available bins correctly use Q(0)=0. Nonempty factors are essential and explicitly required. The join identity and scalar specialization are correct.

### Strong-node bound and pruning

Every induced tournament has at least one Hamilton path, giving r_T>=2^|T|-1. The algebraic second inequality has gap 2(2^(a-1)-1)(2^(b-1)-1). A singleton and a transitive opposite factor attain it at every cyclic output order. Thus N(X) is sharp as a uniform bound on cyclic-node order, not on arbitrary grammar outputs.

H is positive and monotone along both grammar operations: joins multiply positive H values, and the cyclic scalar law is at least either child count. Hence pruning an intermediate profile at H>X cannot discard a successful target. This remains true when a large exterior transitive pad costs no additional H: the pad is treated in scalar closure, and only cyclic nodes receive the order bound.

Flattening outer joins leaves seed or cyclic factors. Large seeds remain explicit scalar generators; a large seed cannot be hidden inside a successful small-order cyclic node, because graph order only increases. Every factor's H is at most X. These facts establish completeness for arbitrary finite seed banks containing K1, including disconnected seeds, not only the particular seed-five calculation.

### Minimum depth and deduplication

At fixed order every child has smaller positive order, so there is no cyclic dependency in the profile dynamic program. Equal complete polynomials have identical future polynomial responses under both operations. Replacing either witness by a lower-depth witness changes neither scalar feasibility nor order and cannot worsen max-based composition depth. Thus profile deduplication is valid for minimum depth, not just membership. Seeds retain depth zero even if the same graph/profile also admits nested constructions. In scalar closure a generator 1 cannot improve a positive depth below the minimum already obtainable from the non-unit factors. K1 provides the depth-zero unit.

The producer omits reversed size pairs. This is legitimate because both polynomial operations are symmetric in their operands, even though the physical ordered join graphs can differ. The independent implementation includes both directions and obtains the same profile counts and depth digest.

### Linear-summary dimension

For every k>=0, upward-shift paths give the top coefficient k!, and terms reaching evaluation at zero vanish. Consequently the response-to-moment matrix is triangular and invertible over Q. On the physical transitive tournaments TT_n, n=1,...,k+1, the moment matrix j^n is nonsingular by the usual Vandermonde argument. Therefore any assignment into Q^d admitting rational linear readouts for all f_0,...,f_k has d>=k+1. No linearity assumption on the assignment itself is needed for this rank obstruction. The moment vector realizes equality. This proves the quantified theorem for all k; the finite rank checks through k=8 are controls, not its justification.

This statement concerns simultaneous linear readouts over arbitrary physical nonempty inputs. It does not rule out nonlinear one-coordinate encodings, specialized finite seed domains, or alternate constructions. The packet correctly states these boundaries.

## Attribution and novelty audit

The cited THM-3121 proves the full quotient path-cover substitution kernel. THM-3134 proves the endpoint-jet/full-profile calculus. THM-3166 supplies the falling-factorial colour coordinate and ordered-join multiplication. The packet accurately treats the displayed binary formula as a specialization, not a new general substitution theorem.

The related THM-3202 **C3 repeated-join moving-jet formula and C-finite obstruction** treats moving factor jets and a non-C-finite obstruction in the repeated-join parameter. It explicitly disclaims an information-theoretic minimality claim. The present linear-response rank theorem is a different, tightly scoped assertion for collar depth, while the need for growing full-profile information is inherited context. There is a separate unrelated THM-3202 filename in the repository, so title-qualified citations should be used if this context is added. No omission here invalidates the stated result or requires attributing this exact rank theorem to that source.

The target-bounded complete grammar compiler and exact minimum-depth calculation are the main scoped deliverables. No exhaustive novelty search was performed and no priority claim is certified.

## Reproduction

Run `python3 -I -B audit.py > AUDIT.json` and `python3 -I -B -O audit.py > AUDIT.optimized.json`; compare receipts byte-for-byte. The checker uses explicit exception gates, not optimization-removable assertions. `REPLAY.json` and `REPLAY.optimized.json` are separate producer replays. Recreate them with `python3 -I -B inputs/verify.py > REPLAY.json` and `python3 -I -B -O inputs/verify.py > REPLAY.optimized.json`; compare each to its corresponding `inputs/CHECK*.json`. The independent checker locates bundled inputs relative to its own file, regardless of the current working directory. `SHA256SUMS` uses paths relative to this audit directory and can be checked with `sha256sum -c SHA256SUMS`. `AUDITED_INPUTS.json` records repository-relative paths and hashes of outside inherited theorem sources, which are context rather than runtime dependencies.
