# CP03 independent publication and replay review

**PASS for the frozen CP03 candidate. Remote publication and CI remain separate checks.**

Reviewed candidate `session190-publish/CP03-candidate/package.json`, SHA-256 `1d7dc7a9da557dbd45ca5bbfc987782f1bd208ffd7c579d3cf069bf72c4264ec`.

## Exact file selection and lossless storage

The candidate contains exactly 64 intended originals: 21 producer spectrum files, 28 spectrum audit files, 6 interface producer files, 7 interface audit files, the updated S190 navigation note, and the CP03 reproduction note. The selections were derived independently from the four frozen manifests, not trusted merely because the publisher's map listed them.

All 74 stored objects satisfy the 48-KiB bound. Each object's byte count, SHA-256 and Git blob SHA-1 was independently checked. Every direct file and every concatenated gzip stream was decoded independently, without importing the publisher's packing or restoration code. Every reconstructed original equals its exact frozen source byte-for-byte. Original bytes total 2,668,296; stored bytes total 740,691.

The package is bound to the final spectrum proof `d1e5debdf607c3eb83f8054c08f94ec1137ef0b4b17aac775973f878429d5479`, its independent manifest `b74f98385607e79bc4ee326bb4162cafcacc655a56a1b696efcb24c98534f7d5`, interface proof `f33466652c3d8a88b437185930da9034869e95327ba480c157d1360e4f3c7566`, and interface audit manifest `c55bbe61f68bf7484e1cf301ce590cc15b9e4638f955eee0731861f6679239ac`. All constituent pins appear in the receipts.

## Portable dependency layout and exact replays

An isolated copy contained only the 62 approved mathematical files in the intended layout:

- `spectrum/spectrum/` and `spectrum/spectrum_audit/`;
- `interface/session190-interface/` and `interface/session190-planar-parity-audit/`.

The spectrum audit's fallback resolves the producer sibling correctly. Its complete inherited order-eight cover is included in the audit directory, so no unpublished sibling is needed. Interface scripts have no external repository imports or runtime file dependencies.

All eight scripts were replayed normally and with optimized Python: spectrum producer verifier, producer census, independent local verifier, independent wrapper verifier, independent census, final source binder, interface producer verifier, and interface independent audit. Each normal/optimized pair agrees byte-for-byte and matches the corresponding frozen output. Both complete spectrum census implementations were rerun in both modes. After all 16 runs, every one of the 62 frozen original files still matches its original hash. Generated executables and caches are not publication files.

## Composite restoration

Independently assembled stored objects from CP01-final, CP02-final-v2 and CP03-candidate were restored with the current v2 restorer. All 133 current original files match independently decoded bytes. Only administrative paths are superseded: the S190 navigation note, RECEIPTS.md, RESTORE_ARCHIVES.py, and the paper-patch reproduction instructions. No frozen mathematical file is overwritten by a later checkpoint.

Normal and optimized verification agree. A second ordinary restoration creates no additional files, and no stored object changes. The checked package pins and superseded records are retained in `PACKAGE_REVIEW_RECEIPT.json`.

All 16 relative links in the updated navigation note resolve to files in the restored composite result. The reproduction instructions correctly require ordered composite manifests and temporary copies for scripts that write outputs.

## Narrative scope

The navigation preserves the substantive mathematical limits: the four-arc result is an exact common-exterior construction, not arbitrary successor or spectrum closure; the common scalar F is not the full rooted response. The two parent exclusions retain the fixed-parent/common-clock hypotheses. The six-port gadget asserts parity equality in arbitrary permitted exteriors, but the constant integer-difference identity only for planar exteriors. Its explicit 2-versus-0 hostile prevents an existence-equivalence reading. The all-order Walsh comparison remains tied to its credited inherited basis theorem.

This review establishes source selection, local object integrity, lossless reconstruction, dependency portability and exact replay. It does not yet assert that these objects exist at an immutable remote commit or that any CI has passed. Those require the separate post-publication evidence.
