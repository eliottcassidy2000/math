from fractions import Fraction
from lie_voa_ss import *
from max_submodule import simple_quotient_dims, per_weight
from small_n4 import SmallN4
from n2_sca import N2SCA, L as NL, J as NJ, GP, GM
import kw_char
P = 2147483647
# Virasoro (3,4) vs Rocha-Caridi
from vir_ss import vacuum_character
alg = Virasoro(Fraction(1, 2), P)
d = per_weight(simple_quotient_dims(alg, 14, [(-1, 0), (-2, 0)]))
print("Vir(3,4) max-submodule dims:", list(d.values())); print("Rocha-Caridi           :", vacuum_character(3, 4, 14))
# sl2 k=-1/2 vs KW
alg = AffineSL2(Fraction(-1, 2), P)
d = per_weight(simple_quotient_dims(alg, 9, [(0, 0), (-1, 2)], height=lambda ch: ch[0]))
print("sl2 k=-1/2 :", list(d.values())); print("KW         :", kw_char.kw(3, 2, 9)[0])
# small N=4 c=-9 vs Schur index |coeffs|
alg = SmallN4(-9, P)
h = Fraction(1, 2)
from small_n4 import E, Fm, G1, G2, B1, B2, L as SL, H as SH
d = per_weight(simple_quotient_dims(alg, 5, [(0, E), (-1, Fm), (-h, G1), (-h, G2), (-h, B1), (-h, B2), (-1, SL), (-2, SL), (-1, SH)], height=lambda ch: ch[0]))
print("small N=4 c=-9:", list(d.values()), "(Schur |coeffs|: 1,0,3,4,9,12,22,36,60,88,135)")
# N=2 minimal models: compare with V/<singular vector>
for kk in (1, 2, 3):
    c = Fraction(3 * kk, kk + 2)
    alg = N2SCA(c, P)
    nplus = [(-h, GP), (-h, GM), (-1, NL), (-2, NL), (-1, NJ)]
    d = per_weight(simple_quotient_dims(alg, 6, nplus))
    print(f"N=2 c={c}: max-submodule dims", [d[k] for k in sorted(d)])
