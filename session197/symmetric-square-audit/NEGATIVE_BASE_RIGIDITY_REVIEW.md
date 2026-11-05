# Independent review: negative-unit base rigidity

Status: **PASS**, no correction. This review concerns the exact marked, root-fixed self-converse transitive-cloning mechanism, not unrestricted Hamilton-path values.

## Frozen source

- `../session197-prime-atoms/NEGATIVE_BASE_RIGIDITY.md`: SHA256 `bea7513641935a9ad0b5af15551dac626ec981667bddc62dbae054bed3339b9d`.
- `../session197-prime-atoms/verify_negative_base_rigidity.py`: SHA256 `9c90530171f672d9fed1e0546edad25e9ec6f4d946110d0f4422ae72a5b18b00`.
- `../session197-prime-atoms/receipt-negative-base-rigidity.json`: SHA256 `1d63e74124e4d9c2b59c2d95bec91a73c0780f409c307065de24f94c5c2de354`.

## Proof audit

Every internal fragment between consecutive root symbols is nonempty, starts at an out-neighbor of the root, and ends at an in-neighbor. Any permutation of those fragments remains a valid root-run word. Their disjoint nonempty labelled vertex sets make the permutation action free, proving `(k-1)! | a_k`. No symmetry assumption is needed here.

For the largest nonzero run index K, the exact exponential expansion has leading base K and positive coefficient a_K. Therefore an eventual real-base identity `H_m=A b^m-1`, with A positive and b greater than one, forces b=K and A=a_K. The remaining distinct integer exponentials are linearly independent, so the identity holds at every positive integer m. In particular K is at least two. No m=0/deletion assertion is made.

A fixed-root anti-involution has no other fixed vertex, hence the original source is odd-order. Reversing each even clone chain extends it to an even-order self-converse graph. The inherited `H=1 mod 4` theorem applies. An even b gives `H=-1 mod 4`, while odd b requires `A=2 mod 4`. For odd b at least five, `(b-1)!` is divisible by four, contradicting that residue. Thus b=3 exactly.

The imported maximal-run/matching equality then gives crossing matching rank two. The already independently audited square theorem yields `A=2q^2` or `4q^2`; the residue selects `2q^2` with q odd. Its independently audited finite exclusion and strong-core reduction give `A>=242`. Attainment at 242 is correctly not asserted. Those inherited claims are not disguised as consequences of the small new control script.

## Independent computation

`check_negative_base_rigidity.py` uses a root-revisiting subset dynamic program, rather than the source's exterior-permutation and root-gap enumeration. It independently counts all run coefficients for all 5,405 labelled rooted tournaments through order five. All 26,705 coefficient divisibility tests pass, including 10,435 nonzero coefficients. This is a control, not the universal proof.

- Independent checker SHA256: `b2b64169ddd88ebe61831b69da80f9ca327f68c141df97b04a9161609b2fa46b`.
- Independent receipt `negative-base-independent.json` SHA256: `d8888397cca208b8e4f8dc9f4e70642f4dc651f87c3d072bd962d7b9c0340e7a`.

Both source and independent checker were replayed normally and under Python `-O`; their respective outputs match byte for byte, and the source reproduces its frozen receipt. Frozen source files were not changed.
