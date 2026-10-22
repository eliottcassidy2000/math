# AUDIT: are the two non-Paley DRT(31) classes of the twice-doubled Paley construction converses of each other?
import sys
sys.path.insert(0, '.')
from aud_newclaims_lib import *
SP31 = paley_S(31)
fam31 = drts_from(doubled(doubled(skew_hadamard_from_drt(paley_S(7)))))
cf = canon(fam31)
reps = {}
for c, D in zip(cf, fam31): reps.setdefault(c, D)
Ds = list(reps.values())
cc = canon([Ds[0], Ds[1], -Ds[0], -Ds[1]])
print(f"n = 31 doubled classes: {len(Ds)}; class(D0^T) == class(D1): {cc[2] == cc[1]}; class(D1^T) == class(D0): {cc[3] == cc[0]}")
