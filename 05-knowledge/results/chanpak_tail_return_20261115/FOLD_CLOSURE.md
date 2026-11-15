# The folded windows cannot close by repainting an ordinal tail

The [fold overlays](TAIL_FOLDS.md) produce genuinely wider absolute windows. The following exact continuation obstruction says what must change if one tries to turn the two overlaid arms into a single matched window. There is a useful analogy with [tail-biting trellises](https://www.conferences.hu/mtns2010/proceedings/Papers/224_Gluesing-Luersen.pdf): a closed trajectory must meet a condition on the actual state, not only agree on an output symbol. Here the state condition is an integer address identity.

## The right state remembers the whole final ordinal word

Let a tournament spine consist of an arbitrary nonempty painted core R followed ordinally by a transitive painted word w of width L. Let h=H(R); it equals the Hamilton count of the entire spine. For any initial marked start weight S_0, the final weighted start is exactly
```
S_final=2^L S_0+h val_L(w),
```
where val_L uses the usual binary order with the first chronological mark most significant. Since h is odd, (h,S_final,L) read the **entire** tail backwards without knowing R or its apex mask:
```
val_L(w) = (h^(-1) S_final) mod 2^L,   in [0,2^L).
```
This follows directly from the singleton cut transfer and invertibility modulo 2^L. It is not an estimate or a bound on the number of bits required by a stochastic decoder.

Thus two different L-bit ordinal tails on different old cores of the same Hamilton count cannot have the same right state, regardless of how old contacts are changed. If the lengths differ, the two words must at least have the same final min(L,L') chronological marks, by applying the same inverse at that shorter length. In particular the distinct addresses in an exact fold cannot be made invisible by tweaking only the older masks while retaining the same ordinal boundary and length. One must change the physical boundary, recontact across it, or use a different kind of continuation.

## Appending more ordinary bits after the fold cannot merge it either

For the exact even-to-odd fold, the two suffix addresses are a=2u and a+1, with u>0. In chronological marks they differ only in their first bit (the least significant address bit is the complement of that mark); their width L is the same. Hence starting from one common state with S_0≥h, their final start readouts differ by exactly h2^(L−1). Each is at least 2^L S_0, so the larger is strictly less than twice the smaller.

Appending different numbers of transitive marks cannot merge these two S's: to match while the appended lengths differ, the larger of the two start values must be at least twice the smaller. Indeed after removing excess steps on the longer side, any common-length equality forces the intervening S's equal: the residual width-r bit corrections differ by less than h2^r, and h,2^r are relatively prime. Each step at least doubles S. At equal appended lengths, common-length equality already requires the same initial S. These arguments cover *all* transitive append lengths; they do not assume a common continuation word.

There is therefore a genuine separation between widening a coverage interval and delivering its arms to one state. The wider folded windows have exact counts but different terminal states; the [one-edge recontact](SINGLE_EDGE_AND_TARGET.md) and [two-sided recontact](TWO_SIDED_RECONTACT.md) show exact state-preserving switches for different structures. The next construction must make a change of the latter kind available at the seam of the former kind. The obstruction does not apply to nonordinal arc changes that invalidate the binary memory formula.

Replay: `python replay_fold_closure.py` checks the terminal-word observer on uniformly sampled literal painted spines, and verifies the distinct exact terminal starts in the three widened fold examples.
