# Session 143: recovered cube/octonion structure

**Result:** an explicit inclusion-lattice and character-action correspondence,
with the exact sign obstruction; plus an octonionic C4 lift of the 180° symmetry
of square-magic equations. No conjecture is solved by these bridges.

* `PROOF.md`: proofs, source recovery, maps, sidecars, hostile controls and scope.
* `verify.py`: portable stdlib checker, independent recursive and polynomial paths.
* `certificate.json`: exact results; `certificate.optimized.json` is identical.
* `negative.log`: expected rejected one-sign mutation.
* `SOURCE_PINS.json` and `sources/`: frozen inherited evidence at repository commit
  `8286d710498a3b4e9b9f1ac683c6eb77716f3f4b`, with blob and content hashes.
* `search_hits.json`: bounded local recovery receipt, not a repository-wide claim.
* `PORTABLE_MANIFEST.json` and `SHA256SUMS`: exact handoff inventory.

## Main conclusions

1. Octonion multiplication is signed XOR; its multiplicative associator is
   `(-1)^det_F2(a,b,c)`. The additive associator is `0` or `±2` times a basis unit.
2. The binary-supported subalgebras form a 16-vertex subspace lattice. Choosing
   three independent labels selects an 8-vertex Boolean cube of R/C/H/O
   subalgebras. There are 28 such selected cubes, not a single canonical one.
3. Multiquadratic Galois sign changes act by the same eight explicit diagonal
   automorphisms on octonions. The representation/lattice correspondence is real;
   the fields and octonions remain different algebras.
4. An explicit left multiplication realizes a signed order-four lift of a
   square's half-turn after squaring its perimeter roots. It preserves the actual
   eight magic-line equations, but does not solve them or change the squared-value
   multiset. Other XOR translations fail; a norm/opposite-pair hostile is included.
5. The nine-cycle is `+4 mod9`. It cannot be faithfully transported to eight
   basis labels. The intended “lambda cube” remains unresolved: type-system cube,
   Lucas cube and scalar λ³ uses are kept separate.

## Reproduce

From any working directory, using this directory's absolute path as `$PACKET`:

```sh
python3 "$PACKET/verify.py" > normal.json
python3 -O "$PACKET/verify.py" > optimized.json
cmp normal.json optimized.json
cmp normal.json "$PACKET/certificate.json"
# Must fail, with cochain-versus-recursion diagnostic:
python3 "$PACKET/verify.py" --mutate-sign
```

All intentional changes are confined to this packet, not the source checkout.
