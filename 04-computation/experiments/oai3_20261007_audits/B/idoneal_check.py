from math import gcd, isqrt
from sympy import factorint
Ns = [int(x) for x in open('exp2_sieve_2p1e11.txt').read().split()]
print("count", len(Ns))
idoneal = sorted(N//4 for N in Ns if N % 4 == 0)
print("idoneal (N/4 for 4|N):", len(idoneal), idoneal)
known65 = [1,2,3,4,5,6,7,8,9,10,12,13,15,16,18,21,22,24,25,28,30,33,37,40,42,45,48,57,58,60,70,72,78,85,88,93,102,105,112,120,130,133,165,168,177,190,210,232,240,253,273,280,312,330,345,357,385,408,462,520,760,840,1320,1365,1848]
print("equals Euler's 65:", idoneal == known65)
def sqfree(n): return all(e == 1 for e in factorint(n).values())
m19 = [n for n in idoneal if n % 4 == 1 and sqfree(n)]
print("squarefree idoneal = 1 mod 4:", len(m19), m19)
def classno(D):  # D<0, count primitive reduced forms
    N = -D; h = 0
    a = 1
    while 3*a*a <= N:
        for b in range(-a+1, a+1):
            if (b*b + N) % (4*a): continue
            c = (b*b + N)//(4*a)
            if c < a: continue
            if gcd(gcd(a, abs(b)), c) != 1: continue
            if c == a and b < 0: continue
            h += 1
        a += 1
    return h
for m in m19:
    t = len(factorint(m)) if m > 1 else 0
    h = classno(-4*m)
    print(m, "h(-4m) =", h, " 2^omega(m) =", 2**t, " genus-field degree 2^(omega+1) =", 2**(t+1), "match" if h == 2**t else "MISMATCH")
# Borwein-Choi: which idoneal n are not xy+yz+zx (x,y,z>=1)?
LIM = 3_000_000
rep = bytearray(LIM+1)
# n = xy+yz+zx with 1<=x<=y<=z: n + x^2 = (x+y)(x+z) ... direct loop
x = 1
while 3*x*x <= LIM:
    y = x
    while x*y + y*y + x*y <= LIM:   # z>=y
        z = y
        base = x*y
        while base + z*(x+y) <= LIM:
            rep[base + z*(x+y)] = 1
            z += 1
        y += 1
    x += 1
non = [n for n in range(1, LIM+1) if not rep[n]]
print("not xy+yz+zx (x,y,z>=1) up to", LIM, ":", non)
print("all idoneal:", all(n in known65 for n in non), " count", len(non))
