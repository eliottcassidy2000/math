# Second proof review: affine lift holonomy and character blocks

23 October 2026. **PASS within the stated complete-bijection, fixed-fibre
interface.** No mathematical correction requested. This is a proof and
interface review independent of the producer; the supplemental finite
Fourier test imports no producer code. It does not independently validate
an exhaustive native move inventory, which the producer correctly leaves
open.

## Frozen producer files reviewed

Directory `session44-lift-holonomy`, exact SHA256:

- README.md: `d664710c472570b59fc8b35db9d2e1316d7a651b18c9989af684d438074b6c03`
- check.py: `1dd8ada55c3859943f5d48a75f882ad9b697a536a03606ad1a8ccabf3701cce0`
- receipt.json: `9585908464baa1413e9603c738257f25a427211bbb9489ebb8e21433b1d000bb`
- source_pins.json: `0480da94234c9c9f54b473d9e69587b57b40f988016c3d2d74d6ac1d68dea26d`
- manifest.json: `ca8118269f5b2fd45806cfbe7e9d27064ba43133547c746801034969ed2aa160`

The immutable source commit is
`a201515762059411a745360467afff224f7db4c1`. Both producer normal and optimized
replays passed and were byte-identical (41,229 recorded checks). The
producer's final receipt-equality check is outside that recorded count.

## Proof audit

1. **Both directions of the component theorem are valid.** Tree transport
   `S_v` converts every lifted edge to its based affine map. Therefore
   paths cannot leave an orbit. Conversely each based generator is realized
   by a root-to-edge-to-root excursion, and reverse edges realize inverses.
   A word in generators followed by a tree path realizes every claimed
   component state. Loops, repeated edges, internal maps, and r=0 do not
   require an omitted argument.

2. **Translation rank is used only when valid.** When all jointly based
   generators are translations, the generated subgroup is exactly their
   F3 span. Cosets give equal-size components; annihilating characters
   separate them. The affine reflection counterexample correctly prevents
   importing equal ternary orbit sizes to general affine holonomy. Section
   existence is distinguished from orbit transitivity.

3. **The Markov hypotheses are sufficient.** Fibre-independent edge
   probabilities and detailed balance with positive base weights imply
   lifted stationarity `pi_v/3^r` and self-adjointness. Zero-probability
   edges must be discarded before computing support holonomy, as the text
   now explicitly says. Poincare and absolute-value gaps are not confused.
   The native C9 example need not be positive semidefinite; positivity is
   asserted only for the separate epsilon family, where it is true.

4. **Character signs and transposes are correct.** Acting on a target
   frequency k by the Markov pullback yields
   `chi_k(c_e) chi_(M_e^T k)` at the source. Reverse affine maps provide the
   adjoint coefficients. Consequently full dual-frequency orbits, not
   individual characters, are the invariant sectors for general affine
   maps. Tree gauge is necessary in the stated holonomy interpretation and
   is retained. In the translation case every annihilating character gives
   the base block; the multiplicity of eigenvalue1 also follows directly
   from the Dirichlet-form equality, independent of a spectral guess.

5. **Converse and native types are honest.** Reversing a path in the
   converse tournament preserves the partner square and sends clock j to
   n-1-j. The order20 construction retains fixed literal source/sink
   owners; its internal owners are transported, not falsely said to stay
   fixed. All guard claims there follow from actual source/target
   isomorphisms for any supplied class potential. No nonexistent numeric
   order20 electrical atlas is asserted.

## Critical interface boundaries

- A partial relation on a fibre is not an affine-cover edge, even if it is
  a partial identity. On two F3 fibres, only the edge0--0 gives five
  components; completing it to the identity matching gives three. A
  nonlinear guard can create exactly this boundary. The theorem excludes
  such completion rather than silently performing it.
- Components apply to the supplied vertex carrier and supplied complete
  transitions. Adding other legal moves may merge them; unseen native
  paths may enlarge the carrier. The text correctly requires completeness
  before upgrading this to a full-native component classification.
- The new two-state selected observer packet in the coupled-moves packet
  is **not itself an F3^r fibre**: two is not a power of three. Its legal
  two-slot move can be represented by a two-vertex base with singleton
  fibres, but it does not license a bijection on arbitrary larger ternary
  charts. This matters when combining the packets.
- Per-slot target preservation, one frozen electrical character, and
  aggregate current preservation remain distinct conditions. The native
  order20 automorphisms preserve the first and hence the other class-based
  guards. Other packet moves need their own exact guard certificates.

## Independent finite interface test

`independent_affine_blocks.py` builds a reversible two-base-state chain
with F3² fibres, a shear and its inverse at one vertex, and an affine
reflection with nonzero translation at the other. It evaluates the full
state operator and independently Fourier-transforms every coefficient in
`Q[w]/(w²+w+1)`. It compares these with the affine transpose/phase formula,
checks that different dual-frequency orbits do not mix, and checks every
Hermitian coefficient pairing. The frequency orbit sizes are1,1,1,6.
It separately verifies the partial/full-matching component boundary above.

Normal and optimized independent outputs agree exactly: **818 checks**.
The producer replays and these independent outputs are preserved beside
this note. No floating arithmetic or optimization-sensitive assertions are
used. Parent's separate native/spectral audit is complementary; this note
does not replace it or claim its checks as our own.
