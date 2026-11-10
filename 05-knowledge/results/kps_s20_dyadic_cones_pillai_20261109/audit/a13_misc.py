# a13: miscellaneous checks.
#  (a) the unique strong classes with H = 5 (order 4) and H = 13 (order 5): Q and the vertex responses B_w;
#  (b) fixed-port tails do NOT have bounded root responses: port response R along tails from small cores
#      (odd part of R at the port after each letter), contradicting "finite banks have bounded root responses";
#  (c) least-squares exponent of N(r0) against r0 over r0 = 63..2047 (HYP-9656 (B): "roughly r0^1.8").
import sys, math
from a_lib import *
S19 = sys.argv[1]
for H0, n in ((5, 4), (13, 5)):
    codes = [int(l.split()[0]) for l in open('%s/cls_%d.txt' % (S19, n)) if l.strip()]
    tours = [from_code(c, n) for c in codes]
    for c, A, r in zip(codes, tours, run(tours)):
        if r['H'] == H0 and is_strong(A):
            print('(a) strong H=%d class code %d (order %d): Q = %d (odd part %d); B_w = %s; odd parts %s'
                  % (H0, c, n, r['Q'], oddpart(r['Q']), r['R'], sorted(set(oddpart(x) for x in r['R']))))

def add_sink(A):
    n = len(A); B = [row + [1] for row in A]; B.append([0] * (n + 1)); return B
def add_bridge(A, u):
    n = len(A); B = [row + [0 if i == u else 1] for i, row in enumerate(A)]
    last = [0] * (n + 1); last[u] = 1; B.append(last); return B
for name, T0, u in (('singleton core', TT(1), 0), ('3-cycle core', [[0, 1, 0], [0, 0, 1], [1, 0, 0]], 0)):
    A = T0; seq = []
    for L in 'SBSBSSBSBSSSB':
        A = add_sink(A) if L == 'S' else add_bridge(A, u)
        r = run([A])[0]; seq.append((L, r['H'], r['R'][u], oddpart(r['R'][u])))
    print('(b) %s, port u fixed, letters SBSBSSBSBSSSB: (letter, h, R_u, oddpart R_u) =' % name, seq)

r0 = [63, 127, 255, 511, 1023, 2047]; N = [127, 431, 2303, 5153, 20921, 55343]
x = [math.log(a) for a in r0]; y = [math.log(b) for b in N]
mx = sum(x) / len(x); my = sum(y) / len(y)
slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)
print('(c) least-squares exponent of N(r0) on r0 = 63..2047: %.3f ; endpoint exponent log(55343/127)/log(2047/63) = %.3f'
      % (slope, math.log(55343 / 127) / math.log(2047 / 63)))
