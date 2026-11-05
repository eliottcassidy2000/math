# Independent review: sharp dyadic self-converse C2 image

Session197. **PASS, no correction.** This review binds the separate corollary in `session197-prime-atoms/dyadic-profile-universality/`:

- SELF_CONVERSE_C2.md: `fc8a488c39f7786e9cc74e1e2d05d551a8e5ab3493fbe31b1fb878b3fb55dfd2`.
- verify_self_converse_c2.py: `d3f9e94960b886d0e35721999c46a57f5fd534bcfff2044c657917e46ac48e82`.
- Both frozen receipts: `df5cb3f850ab424cb79ed8180b73d14f8867b18f962f0878e1be447f7fd92260`.

## Dependencies and all-order argument

The joint-profile proof is the exact `49bfc941828bf680dcaca1d5e4ec448f3b689d0f72c71f164bcb8b6bf0088c64` snapshot, independently reviewed in `session197-lowcycle/dyadic-profile-audit/REPORT.md` (hash `1ee6acd8cffa9d4bbe42edde98f13bc1e4db7233a77e8e43ee81a333f367d7e4`) and published at commit33a63e57e68b818d2282721ef40b34aca420c8ab. Its exact endpoint/deletion coordinates and quadratic order bound have the needed scope.

The self-converse structural dependency is the exact `3e6d6e7344413504f299329fec05304600deef9432911f7aa9fecbd8fd4cb908` snapshot, separately independently reviewed and published at c62aad43732ca55d1c10a3372d7ed9da7473218a. This audit does not infer that theorem from the present residue tests.

The odd-square bit lift is correct at j=3 and at every subsequent bit: changing x by2^(j-1) changes the square by2^j modulo2^(j+1), while preserving every lower bit. In the5-mod8 case U=2 shifts the desired square residue by4; in the1-mod8 case U=4 shifts it by24. Both leave an odd-square target. D=x-U stays odd, so the joint-profile source theorem applies.

One or two universal sinks give actual weak start-one sources with E=0,d=1, prescribed even U, and unchanged odd D. Weak sources are expressly allowed by the diagonal theorem. The use of the same actual source graph on both sides supplies a physical anti-involution, not only equal scalar data. The order calculation is exactly2(4k²+3k+2+2)+1=8k²+6k+9 in the larger padding case. Small moduli are correctly handled by lifting to k=3.

Neutral loops give unbounded source D, not merely longer presentations with a fixed count. The positive diagonal expression is at least D², so it gives infinitely many distinct output counts in each permitted nonregular dyadic residue. The exceptional residue comes from the single R5 graph and is not claimed to have infinitely many counts. Necessity and the all-dyadic obstruction criterion keep this exception separate. In particular, satisfying every congruence test is not claimed to imply an exact integer realization.

## Source replay and separate controls

Normal and optimized Python source runs, performed in a copy, reproduced the frozen receipt exactly. No frozen source was changed.

The independent standalone `dyadic-self-converse/check.py` imports no source. It maintains all four odd square roots while extending the modulus, instead of using the source's single-root lift. Its modular word constructor reverses both the correction-loop order and each two-loop pair. It computes the resulting actual rooted one-apex Hamilton profile by assigning each transitive deletion vertex to the path before or after the root, not by the source's halving recurrence. Thus the source residue checks compare physical assignment counts, including endpoint and deletion data.

Normal and `-O` independent runs agree on:

- all510 allowed nonregular residues at moduli2³ through2¹⁰;
- twelve additional exact targets through modulus2⁶⁴;
- 510 independent physical source-profile checks;
- 30,377 literal anti-involution edge checks;
- sixteen actual final graph H/C counts by unmemoized directed Hamilton-word recursion, including both weak-padding choices.

Checker hash: `8f910216130671ff5ee737816950cad0e79785184430cc6f2c4d65a268549e84`.

Both independent receipt hashes: `8596524638f38500d2f6e17a586d1395430833b21e11e156f0e5e0ee275ff0b0`.

Run `python3 dyadic-self-converse/check.py` and its `-O` form from this review directory. The checker uses only the standard library and writes its JSON to stdout. The source replay copy is in `dyadic-self-converse/source/`.

The conclusion is the exact dyadic image of the specified structural class. There is no exact-target sufficiency, cofinite coverage, prime-infinitude or unrestricted-spectrum conclusion.
