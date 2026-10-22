# Audit2 (D): size of the worst residual y - y_pred in the populations (uniformity of HYP-9365), using H values verified
# in e_populations_*.out (regular n = 11 by Held-Karp here; circulant n = 25 from circ25_H.txt, all recomputed there).
import sys, re, subprocess
import numpy as np
sys.path.insert(0, '.')
from e_predictor import two_loop_dense, y_exact, rho_circ_exact
from common import code_to_A, skew, circulant_A, regular_classes
HC = sys.argv[1]
codes = regular_classes(11)
Hs = [int(x) for x in subprocess.run([HC, 'hk'], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split()]
D = []
for c, Hv in zip(codes, Hs):
    t = two_loop_dense(skew(code_to_A(c))); y, _ = y_exact(Hv, 11, t['rho']); D.append(y - t['y_pred'])
D = np.array(D)
print(f"regular n=11: y - y_pred in [{D.min():+.4f}, {D.max():+.4f}], n*(y - y_pred) in [{11 * D.min():+.3f}, {11 * D.max():+.3f}]")
rows = []
for line in open('../../lattice_bridges_20261021/circ25_H.txt'):
    m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if m: rows.append((tuple(int(x) for x in m.group(1).split(',')), int(m.group(2))))
D = []
for conn, Hv in rows:
    t = two_loop_dense(skew(circulant_A(25, conn))); y, _ = y_exact(Hv, 25, rho_circ_exact(25, conn)); D.append(y - t['y_pred'])
D = np.array(D)
print(f"circulants n=25: y - y_pred in [{D.min():+.4f}, {D.max():+.4f}], n*(y - y_pred) in [{25 * D.min():+.3f}, {25 * D.max():+.3f}]")
