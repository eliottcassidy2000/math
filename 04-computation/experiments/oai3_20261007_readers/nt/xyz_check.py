# Lane nt: check the elementary characterizations tied to idoneal numbers, n <= N.
import numpy as np, sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
euler=[1,2,3,4,5,6,7,8,9,10,12,13,15,16,18,21,22,24,25,28,30,33,37,40,42,45,48,57,58,60,70,72,78,85,88,93,102,105,112,120,130,133,165,168,177,190,210,232,240,253,273,280,312,330,345,357,385,408,462,520,760,840,1320,1365,1848]
def unrepresented(strict):
    rep = np.zeros(N + 1, dtype=bool)
    x = 1
    while True:
        y0 = x + 1 if strict else x
        if x * y0 + (x + y0) * (y0 + (1 if strict else 0)) > N: break
        y = y0
        while True:
            z0 = y + 1 if strict else y
            n0 = x * y + (x + y) * z0
            if n0 > N: break
            rep[n0::(x + y)] = True
            y += 1
        x += 1
    return [n for n in range(1, N + 1) if not rep[n]]
u_strict = unrepresented(True)      # ab+ac+bc, 0<a<b<c   (Rains)
u_weak = unrepresented(False)       # xy+yz+zx, x,y,z>=1  (Borwein-Choi)
print("N =", N)
print("not ab+ac+bc (0<a<b<c):", len(u_strict), "== Euler's 65:", u_strict == euler)
print("not xy+yz+zx (x,y,z>=1):", u_weak)
print("  == {1,4} U {idoneal = 2 mod 4}:", u_weak == sorted([1, 4] + [n for n in euler if n % 4 == 2]))
