# A reflected three-prime window at a single absolute height

This complements [the dihedral construction](DIHEDRAL_HEIGHT.md). It is a numerical reflection of the reader outputs, not a claimed vertex action on the three frames.

On outgoing masks `(254,92,248,240,96,194,128,18,255)`, switch the cuts of masks 52, 77, 155. The literal `(H,B)` profiles are `(2033,34288),(2039,34288),(2045,34288)`. This attains the sharp step 6 for same-order equal-height switches, as in [FRAME_REPAIR](FRAME_REPAIR.md), but at both successive steps in one lattice.

The existing slope-96, intercept-minus-one paired readers in [replay.py](replay.py), after k ordinal singletons on the source, produce `p=96·34288·2^k−H` and `7p`. At k=10, the three exact p values are 3,370,645,507; 3,370,645,513; and 3,370,645,519. Their upper counts are 23,594,518,549; 23,594,518,591; and 23,594,518,633. The ordinary values are prime by exact trial division. A request in this menu selects exactly one of the signed offsets by its residue modulo 96 because the three distinct H values have distinct residues. It then selects k by testing whether the resulting positive quotient by 96·34288 is a power of two.

Independent replay (run this alongside the earlier files):

```python
from frame_replay import cut
from replay import prime, profile
A=(254,92,248,240,96,194,128,18,255)
for mask,h in zip((52,77,155),(2033,2039,2045)):
    assert profile(cut(A,mask)) == (h,34288)
    assert prime(96*34288*1024-h)
```

This is an absolute decoder for one three-offset menu, not a placement theorem for arbitrary targets.
