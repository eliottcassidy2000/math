# Independent audit: append-only binary children

**PASS, without a mathematical correction to the frozen packet.** The all-order source-slice classification, arbitrary-rebasing quantifiers, append-length reduction and exact intersection are valid in the stated physical domain. No eventual-coverage theorem follows.

## Frozen inputs

`INPUTS.json` records the source packet. The audited proof hash is
`b21204b9de7160c9b94aad791ae78cf1f4d587b4d2704d8233d239737e85f54a`;
the audited author verifier hash is
`a9aabb98e8bc3cfb02df9e5f07efeca78c24863cef14fa4f4afec2ebf27fcfb9`.
The source proof and verifier remained unchanged throughout the audit. The statement may be promoted from pending audit without changing its mathematical content.

## All-order reasoning

1. **Domain and inherited dependencies.** The proof uses actual nontrivial rooted start-one sources. The retained dependencies explicitly state completeness of the physical grammar, including its weak-tail decomposition, and the graph-level balanced inverse. Even without importing that structural completeness, every conclusion is valid for the explicitly defined seed/O/Z/J orbit. No cone point is treated as reachable merely because its coordinates are nonnegative.

2. **Growth is uniform and independent of rebasing.** The three displayed readout rows are correct exact linear identities. Each is coefficientwise at least three halves of the F row on the nonnegative physical cone. Thus any extension of length at least two has readout at least 9n/4, strictly greater than 2n+1 when n>=5. This excludes every longer word, not merely words up to a tested length. Choosing a different actual representative of the same n does not affect the inequality. The O/J offset is a negative positive-odd integer; it cannot supply the plus child.

3. **Negative Z slice.** Removing final Z letters is legitimate because Z preserves h. A remaining O cannot have h=-1. A remaining J must have had D=1, whose nontrivial sources are precisely the transitive O chains. Their J images and the seed give exactly the stated numerical profiles t.Z^m. Restoration of the Z tail is exhaustive.

4. **Positive Z slice.** After removing the Z tail, the final operation must be O. Its predecessor has f=0 and a+p=1. The two cases (a,p)=(1,0) and (0,1) are exhausted correctly. The first balanced inverse forces e=1 by U_W>=D_W. In the second case, when e>0, the actual final J and balanced inverse force e=2^j. Crucially j=1 comes from the singleton inverse W, so that boundary is included. The three displayed families give exactly every actual profile with h=1, with the correct formulas.

5. **O/J minus slice.** The E=0 case gives D=1 directly from D-2E=1. For E>0 the last-letter alternatives are exhaustive. The proposed even-D predecessor of Z is rightly excluded. The J case reduces to f=0,d=1; the supplementary E=0 argument correctly proves that this forces the seed, instead of incorrectly inferring it from E=0 alone. The balanced E=1 branch yields precisely J(t.Z^j).

6. **Exact intersection.** The powers-of-two valuation calculations prove that the first minus family meets the three plus families only at 5,9,17 respectively. The second minus family can meet only the plus family 10*2^k-1 by reduction modulo three. The remaining valuation cases give only 9. Therefore the two independently rebased branches are jointly available exactly at 5,9,17.

7. **Hostile 25 and scope.** Count25 is represented, but every actual F=25 source has only one-letter outputs45 or63. Growth excludes all longer words from reaching49 or51. The theorem therefore forbids this prescribed-parent operation, not the child values themselves. Both49 and51 have literal C1 witnesses. Different surgeries, terminal replacement, changes of parent readout, and C2/C3 interfaces remain outside the theorem.

## Boundary cases

The theorem deliberately assumes odd n>=5. The only F=3 source is the seed, whose three one-letter outputs all have F=5; any second step is at least 15/2 and hence cannot reach7. Thus3 does not create a missing both-child exception. The singleton is outside the nontrivial parent-source domain and is used only where the inherited balanced inverse permits it. Count1 is not silently admitted as a source parent. Empty append words are excluded and could not reach either child for n>=5 anyway.

## Independent exact controls

The standalone `check.py` imports no author code or data. It verifies the linear transfer/readout identities on a coordinate basis, then independently generates all126,910 physical numerical profiles with U<=5000 by increasing cover buckets. It recovers exactly:

- 13 profiles with h=-1;
- 76 profiles with h=1;
- 24 profiles with D-2E=1;
- the complete finite parent-count sets through5001;
- the arithmetic intersection {5,9,17} through100,000,000;
- 1,926 two-letter hostile checks;
- both F=25 source profiles and their child set {45,63}.

The cutoff is complete for every target<=5001 because F>=U+1. Literal memoized path/cycle counts independently check the small parent controls, hostile25 and witnesses49/51, including the frozen source word JZJJZ. The author separately counts all six actual binary child controls. Its normal and optimized replays agree with both frozen author receipts exactly.

The independent normal and optimized outputs are identical. These finite controls support the proof but do not substitute for its all-order arguments.

Reproduce from this audit directory:

    python3 check.py > check-normal.json
    python3 -O check.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json NORMAL.json

The optional D-interface extension discussed during coordination is not included in this audited frozen proof and requires its own frozen addendum.
