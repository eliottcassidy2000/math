# Independent check of #158 Lemma (seven-point certificate) and the Moser gadget facts.
import itertools
from mpmath import mp, mpf, sqrt, mpc, sin, acos, pi, fabs
import sympy as sp

mp.dps = 60
s3, s11 = sqrt(3), sqrt(11)
A = mpc(s3/2, mpf(1)/2); B = mpc(s3/2, -mpf(1)/2); T = mpc(s3, 0)
u = mpc(mpf(5)/6, s11/6)
G = {'0': mpc(0,0), 'A': A, 'B': B, 'T': T, 'uA': u*A, 'uB': u*B, 'uT': u*T}
rot = mpc(mpf(35)/37, mpf(12)/37); shift = mpc(mpf(-290)/250, mpf(149)/250)
Z = {g: rot*(z+shift) for g, z in G.items()}

# exact edge check with sympy (unit distances preserved by rotation+translation)
S3, S11, I = sp.sqrt(3), sp.sqrt(11), sp.I
Ae = (S3+I)/2; Be = (S3-I)/2; Te = S3; ue = (5+I*S11)/6
Ge = {'0':0,'A':Ae,'B':Be,'T':Te,'uA':ue*Ae,'uB':ue*Be,'uT':ue*Te}
names = list(Ge)
edges = []
for a, b in itertools.combinations(names, 2):
    d = sp.nsimplify(sp.expand((Ge[a]-Ge[b])*sp.conjugate(Ge[a]-Ge[b])))
    d = sp.simplify(d)
    if d == 1: edges.append((a, b))
print("unit edges (exact):", len(edges), edges)

# 3-colourability by brute force
def colorable(k):
    for col in itertools.product(range(k), repeat=len(names)):
        c = dict(zip(names, col))
        if all(c[a] != c[b] for a, b in edges): return True
    return False
print("3-colourable?", colorable(3), " 4-colourable?", colorable(4))

# region test: 0<l<4, l!=1, P<P' ; also check the sector-sign product directly
ok = True
for g, z in Z.items():
    xi, up = z.real, z.imag
    l = xi**2 + up**2
    P = up**2*(3*xi**2-up**2)**2
    Pp = l**3*(1-l/4)*(l-1)**2
    r = sqrt(l); phi = mp.atan2(up, xi); dl = acos(r/2)
    prod = sin(3*(phi-dl))*sin(3*(phi+dl))
    print(f"{g:3s} xi={float(xi):+.6f} ups={float(up):+.6f} l={float(l):.6f} P={float(P):.6g} P'={float(Pp):.6g} margin={float(Pp-P):.4g} sinprod={float(prod):+.4g}")
    ok &= (0 < l < 4) and (P < Pp) and (prod < 0)
print("all seven placed vertices strictly in the three-label region:", ok)
# smallest margin
print("min (P'-P):", min(float((lambda z: (z.real**2+z.imag**2)**3*(1-(z.real**2+z.imag**2)/4)*((z.real**2+z.imag**2)-1)**2 - z.imag**2*(3*z.real**2-z.imag**2)**2)(z)) for z in Z.values()))
# gamma facts used in Lemma 6.x (singleton-pair alternation)
gam = 2*mp.asin(1/(2*s3))
print("cos(gamma) =", mp.nstr(mp.cos(gam), 30), " |1-e^{i gamma}| =", mp.nstr(abs(1-mp.expjpi(gam/pi)), 30), " 1/sqrt3 =", mp.nstr(1/s3, 30))
