#!/usr/bin/env python3
# audit K2: explicit order-13 sandwiches C3 => Z => C3^op for 711 = 9*79 and 855 = 9*95 (Z = the order-7 strong
# self-converse census classes with H = 79, 95; codes in THM-4922's census convention). Reuses audit K1's helpers.
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "..", "audit_K1", "aud_275.py")).read().split("B5 = decode(")[0]
ns = {}; exec(src, ns)
decode, H, converse, join, self_converse = ns["decode"], ns["H"], ns["converse"], ns["join"], ns["self_converse"]
C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
for code, e in ((51714, 79), (90561, 95)):
    Z = decode(code, 7)
    T = join(C3, Z, converse(C3))
    print(f"Z (code {code}): H = {H(Z)}, self-converse {self_converse(Z)}; C3 => Z => C3^op: order {len(T)}, H = {H(T)} "
          f"(= 9 * {e} = {9 * e}), self-converse {self_converse(T)}")
