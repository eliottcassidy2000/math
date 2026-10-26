# Audit extra: the explicit reduced form of 2^21 TT7(P_q), q = 3 (mod 4), in the documents' labellings, with the
# octahedron (T_oct = (q-5) T_W4) and K7 - C7 (T = eps (q-6) T_(K6-P6), THM-4780 (4)) reduced; eps fixed numerically.
# Output: A7, B7, C7 and the 37 genuine seven-vertex survivors with their multiplicities (canonical labellings).
import sys
sys.path.insert(0, '.')
import audit_census_symbolic as acs
from audit_graph_sums import S_norm
res = acs.census(7)
lab = acs.express(res, 7)
A7 = lab['1']; B = lab['T_W4']; C = lab['T_K6-P6']; O = lab['T_oct']
B7 = acs.padd(B, acs.pmul(O, (-5, 1)))
reps, info = acs.class_info(7)
c7 = None
genuine = []
for key, v in res.items():
    if key == 'one' or key[1] != 7: continue
    mask = key[2]; deg = acs.degrees(7, mask)
    if all(d == 4 for d in deg):
        c7 = (mask, v)
    else:
        genuine.append((mask, v))
mask, v = c7
E = [acs.EDGES[7][n] for n in acs.bits(mask)]
Ekp = [e for e in acs.EDGES[6] if e not in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]]
r = S_norm(11, 7, E) / ((11 - 6) * S_norm(11, 6, Ekp))
r2 = S_norm(7, 7, E) / ((7 - 6) * S_norm(7, 6, Ekp))
assert r == r2 and abs(r) == 1
eps = int(r)
C7 = acs.padd(C, acs.pmul(acs.pscale(v, eps), (-6, 1)))
out = [f'2^21 TT7(P_q) = A7(q) + B7(q) T_W4 + C7(q) T_(K6-P6) + sum over the 37 genuine survivors m_H T_H   (q = 3 mod 4)',
       f'  A7 = {acs.pstr(A7)}',
       f'  B7 = {acs.pstr(B7)}   (from non-dense graphs {acs.pstr(B)}, octahedron {acs.pstr(O)} times (q-5))',
       f'  C7 = {acs.pstr(C7)}   (from non-dense graphs {acs.pstr(C)}, K7 - C7 with m = {v[0]}, eps = {eps}, times (q-6))',
       f'  A7 divisible by q(q-1): {acs.peval(A7, 0) == 0 and acs.peval(A7, 1) == 0}',
       f'  37 genuine survivors (canonical labelling: vertices 0..6, order orientation; complement edges listed):']
for mask, v in genuine:
    comp = [e for e in acs.EDGES[7] if not mask >> acs.EIDX[7][e] & 1]
    out.append(f'    m = {v[0]:+4d}  complement {comp}')
print('\n'.join(out))
open('audit_tt7_explicit.out', 'w').write('\n'.join(out) + '\n')

# the genuine seven-vertex part G(q) = 2^21 TT7 - A7 - B7 T_W4 - C7 T_(K6-P6) at small q (TT7 from the source-sink
# decomposition), to show it does not vanish identically
from audit_graph_sums import W4
from audit_tt7_structure_helpers import tt7
lines2 = []
for p in (7, 11, 19, 23, 31, 43):
    TW = p * (p - 1) * S_norm(p, 5, W4); TK = p * (p - 1) * S_norm(p, 6, Ekp)
    known = acs.peval(A7, p) + acs.peval(B7, p) * TW + acs.peval(C7, p) * TK
    G = (1 << 21) * tt7(p) - known
    lines2.append(f'  q = {p}: TT7 = {tt7(p)}; genuine part G(q) = 2^21 TT7 - A7 - B7 T_W4 - C7 T_(K6-P6) = {G}')
print('\n'.join(lines2))
open('audit_tt7_explicit.out', 'a').write('\n'.join(lines2) + '\n')
