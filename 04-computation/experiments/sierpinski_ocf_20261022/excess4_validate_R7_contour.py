# Contour check for excess4_validate.out: R_7 showed relative errors ~2e-9 at K = 32, r = 0.5 (vanishing z^2, z^4
# coefficients ~1e-10).  If that is Cauchy-sum aliasing, K = 64 points on r = 0.35 must remove it.
# Output: excess4_validate_R7_contour.out.  mac-mini-2026-10-22-loop4.
import sys
sys.path.insert(0, '.'); sys.path.insert(0, 'audit2')
import loop_engine5 as LE5
from b_mscaling import mscaling
from common import rot_A, skew
S = skew(rot_A(7)); Sh = S / 7
X = mscaling(S, 14, verbose=False)
co = LE5.series_general(Sh, 4, 14, K=64, r=0.35)
for L in range(2, 15, 2):
    ex = float(X[L].get(4, 0))
    print(f"L={L:2d}: X_4 engine {co[L]:+.14e} exact {ex:+.14e} diff {co[L] - ex:+.1e}", flush=True)
