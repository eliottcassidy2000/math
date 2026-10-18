import sys
from fractions import Fraction
from small_n4 import SmallN4, E, H, Fm, G1, G2, B1, B2
from lie_voa_ss import *
c = Fraction(sys.argv[1]); Dmax = int(sys.argv[2]); P = int(sys.argv[3]) if len(sys.argv) > 3 else 2147483647
alg = SmallN4(c, P)
sv = find_singular(alg, 2, (0,))
print(f"small N=4 c={c}: singular vectors at weight 2 charge 0: {len(sv)}")
ss = SingularSupport(alg, Dmax, sv)
h = Fraction(1, 2)
ss.run(gen_modes=[(h, G1), (h, G2), (h, B1), (h, B2), (1, E), (1, H), (1, Fm)], extra_ops=[(0, Fm), (0, E)])
print("generators (D, Li, charge, #):", ss.gens_found)
