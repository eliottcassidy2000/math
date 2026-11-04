# Independent audit: preceding-context terminal block laws

**PASS, no corrections required.** Reviewed `mining/PROOF.md`, `verify.py`, and `FRONTIER.md`, and compared the local pinned copy of `transport/TWO_DIGIT_ENTRY.md` with the inherited signatures and construction used here.

## Mathematical verification

The inherited B2 signature has endpoint coordinate e=1, as correctly retained by the author. Its difference from3B0 is (0,24,0,18,−2,−2,0); B4−3B0 is (−2,8,0,0,−6,−2,−2). Substitution into the actual ordered-sum interface gives exactly

- Δf2=0 and Δ(h−j)2=2h;
- Δf4=−2f and Δ(h−j)4=4h+2j+2e.

The stated full-attachment formulas therefore follow without assuming transitivity or strongness of the preceding context. The actual marked terminal SCC proves all three full attachments strong. The empty-context tuple is used only in the declared coordinates and is explicitly not claimed to be a seven-coordinate identity.

For outer TT clocks, uniform+4 at two distinct clocks forces f=0. For a nonempty context any port supplies a two-slot cover consisting of that singleton and a Hamilton path of the complement, so f=0 exactly when the context is unmarked. Then j=e=0, leaving h=1. The tournament equivalence h=1 iff transitive is justified by the path-insertion argument: every path extends while retaining its old relative order, so the three distinct Hamilton orders of an induced triangle remain distinct after extension. The paired-recipe iff statement is correct and explicitly does not exclude other +4 recipes.

The shortened two-final-port guard has the stated actual signatures and Δj=−2, with Δh=Δf=Δs=0. Thus its broader+2 extension follows by the same proof. This is an additional actual family, not a returning grammar. The direct two-birth schedule preserves every old named edge and has a marked terminal singleton after its first step and a marked terminal triangle after its second.

## Independent computation

`verify.py` in this audit directory does not import the author's checker. It obtains each seven-coordinate signature by enumerating literal directed path permutations of every induced subset, independently of the author's endpoint recurrence and of the ordered-sum formulas. Full graphs are counted by forward path-extension DP. It checks:

- All75 labelled preceding context/port states through order3, including the empty context.
- All675 full attachments with outer TT1, TT2 and C3.
- All225 literal name-preserving two-birth schedules, including exact equality to the claimed final adjacency, preservation of every old edge, and intermediate/final strongness.
- All75 paired uniform-locus tests at two distinct TT clocks.
- Both three-vertex and shortened two-vertex entry signatures.

Normal and optimized audit runs byte-match. The author's own normal and optimized checker runs also byte-match, with75 contexts and675 full graphs. See `receipt.json` and saved stdout controls.

## Attribution and scope

The local pinned common-entry proof contains exactly the inherited three signatures and literal birth schedule used by the author. The claimed new contribution is their arbitrary preceding marked-context law, its shortened two-port source, and the exact paired-recipe obstruction, rather than priority for the entry itself. `FRONTIER.md` distinguishes targeted mining from an exhaustive novelty survey and explicitly defers actual missing44/native-lineage certification to that separate lane. Neither the new proof nor this audit supplies an unrestricted spectrum theorem, a native guarded program, or recursive two-digit closure.
