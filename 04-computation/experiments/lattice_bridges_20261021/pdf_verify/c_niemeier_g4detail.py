import sys; sys.path.insert(0, '.')
exec(open('c_niemeier.py').read().split("print(\"\\nGenus-4 test")[0].replace('print(', 'pass; (lambda *a, **k: None)('))
import numpy as np
for G in C4:
    d30 = count_L(LAT['E8^3'], G) - count_L(LAT['D16E8'], G)
    det = round(np.linalg.det(G.astype(float)))
    print("det", det, "offdiag", [int(G[i,j]) for i in range(4) for j in range(i+1,4)], " F_30 coeff =", d30,
          " counts E8^3 / D16E8:", count_L(LAT['E8^3'], G), count_L(LAT['D16E8'], G))
