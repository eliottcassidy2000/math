"""S(4) = sum_n C(2n,n)C(4n,2n)/((4n+1)64^n):  an INDEPENDENT route to the
lemniscatic reduction of THM-3012 addendum 5, via Jacobi sn/dn FIRST MOMENTS over
a quarter period, plus the mechanism of the k=4 wall, the CM curve behind it, and
a bounded PSLQ exclusion with live positive controls.

Chain (each link refereed below):
  C1  the controlling 2F1 of THM-3012 (5) at k=4 is 2F1(1,1;3/2;z)
      = asin(sqrt z)/sqrt(z(1-z)) -- ELEMENTARY.  Schwarz angles (1/2,1/2,0):
      k=4 is the ONLY k>=4 with a cusp (nu=|1-4/k|=0).  So the obstruction at
      k=4 is not in the inner 2F1 but in the OUTER theta-integration.
  C2  pi S(4) = 2 J,  J = int_0^1 [asin v + asinh v]/sqrt(1-v^4) dv
      (equivalently the owner's (2 sqrt2/pi) int_0^1 K(k)/(2-k^2) dk).
  C3  J = JA + JB with
        JA = int_0^{pi/2} th dth/sqrt(1+sin^2 th)      = (1/sqrt2) int_0^K t dn(t,1/2) dt
        JB = int_0^{K(-1)} log(sn(u,-1)+dn(u,-1)) du   = (1/2)     int_0^K t sn(t,1/2) dt
      Fourier series of dn and sn at the CM nome q = e^{-pi} give
        S(4) = varpi/4 + (2 varpi/pi^2)(2V - T),
        T = sum_{n odd} sech(pi n)/n^2,
        V = sum_{n>=0} (-1)^n/((2n+1)^2 sinh((2n+1)pi/2)).
  C4  bridge to addendum 5: Lambda = 2G + D = pi^2/8 + 2V - T, and 2V = T - U/2 + pi^2/12.
  C5  the curve is E32: y^2 = x^3 - x (CM by Z[i]); control L(E32,1) = varpi/4.
  C6  bounded PSLQ exclusion for W = 2V - T = Lambda - pi^2/8.
"""
import mpmath as mp, itertools, time
mp.mp.dps = 170
PI = mp.pi; r2 = mp.sqrt(2); G = mp.catalan
OK = []
def ck(name, val, tol=mp.mpf(10)**-150):
    ok = abs(val) < tol; OK.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {mp.nstr(val,5)}", flush=True)

print("C1  the controlling 2F1 at k=4 is ELEMENTARY (the wall mechanism)")
for z in [mp.mpf('0.3'), mp.mpf('0.77'), -mp.mpf('0.4')]:
    ck(f"2F1(1,1;3/2;{z}) = asin(sqrt z)/sqrt(z(1-z))",
       mp.hyp2f1(1,1,mp.mpf(3)/2,z) - mp.asin(mp.sqrt(z))/mp.sqrt(z*(1-z)))
for k in [1,2,3,4,5,6]:
    lam = mp.mpf(2)/k; nu = abs(1-mp.mpf(4)/k)
    print(f"    k={k}: Schwarz (lam,mu,nu)=({lam},{lam},{nu})  sum={2*lam+nu}"
          f"  {'SPHERICAL' if 2*lam+nu>1 else 'EUCLIDEAN'}{'  CUSP' if nu==0 else ''}")

print("C2  the owner's elliptic representation")
Iell = mp.quad(lambda k: mp.ellipk(k**2)/(2-k**2), [0,1])
J  = mp.quad(lambda p: (PI/2 - p + mp.asinh(mp.cos(p)))/mp.sqrt(1+mp.cos(p)**2), [0, PI/2])
S4 = 2*J/PI
print("  S(4) =", mp.nstr(S4,100))
ck("S(4) = (2sqrt2/pi) int_0^1 K(k)/(2-k^2) dk", S4 - 2*r2/PI*Iell)
mp.mp.dps = 120
ck("S(4) = 3F2(1/4,3/4,1/4;1,5/4;1) [dps120]",
   S4 - mp.hyper([mp.mpf(1)/4,mp.mpf(3)/4,mp.mpf(1)/4],[1,mp.mpf(5)/4],1), mp.mpf(10)**-110)
ck("S(4) = int_0^1 2F1(1/4,3/4;1;u^4) du [dps120]",
   S4 - mp.quad(lambda u: mp.hyp2f1(mp.mpf(1)/4,mp.mpf(3)/4,1,u**4),[0,1]), mp.mpf(10)**-110)
mp.mp.dps = 170
# v-form has a (1-v)^{-1/2} endpoint singularity: tanh-sinh attains ~1e-100 here,
# so it is checked at 1e-90.  The phi-form (v = cos phi) is analytic -> full precision.
ck("S(4) = (2/pi) int_0^1 (asin+asinh)/sqrt(1-v^4) dv",
   S4 - 2*mp.quad(lambda v:(mp.asin(v)+mp.asinh(v))/mp.sqrt(1-v**4),[0,1])/PI, mp.mpf(10)**-80)

print("C3  Jacobi first-moment form and the Lambert reduction")
vp = 2*mp.ellipk(-1); Kl = mp.ellipk(mp.mpf(1)/2); Km = mp.ellipk(-1)
JA = mp.quad(lambda t: t/mp.sqrt(1+mp.sin(t)**2), [0, PI/2]); JB = J - JA
ck("varpi = sqrt2 K(1/sqrt2) = Gamma(1/4)^2/(2 sqrt(2 pi))",
   vp - mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*PI)))
ck("JA = (1/sqrt2) int_0^K t dn(t,1/2) dt",
   JA - mp.quad(lambda t: t*mp.ellipfun('dn',t,mp.mpf(1)/2),[0,Kl])/r2)
ck("JB = (1/2) int_0^K t sn(t,1/2) dt",
   JB - mp.quad(lambda t: t*mp.ellipfun('sn',t,mp.mpf(1)/2),[0,Kl])/2)
ck("JB = int_0^{K(-1)} log(sn(u,-1)+dn(u,-1)) du   [regulator form]",
   JB - mp.quad(lambda u: mp.log(mp.ellipfun('sn',u,-1)+mp.ellipfun('dn',u,-1)),[0,Km]).real)
T = mp.nsum(lambda j: mp.sech(PI*(2*j+1))/(2*j+1)**2, [0, mp.inf])
V = mp.nsum(lambda j: (-1)**j/((2*j+1)**2*mp.sinh(PI*(2*j+1)/2)), [0, mp.inf])
W = 2*V - T
print("  T =", mp.nstr(T,60)); print("  V =", mp.nstr(V,60)); print("  W = 2V-T =", mp.nstr(W,60))
ck("JA = pi varpi/8 - (varpi/pi) T", JA - (PI*vp/8 - vp*T/PI))
ck("JB = (2 varpi/pi) V",            JB - 2*vp*V/PI)
ck("S(4) = varpi/4 + (2 varpi/pi^2)(2V - T)", S4 - (vp/4 + 2*vp*W/PI**2))

print("C4  bridge to addendum 5 (independent confirmation of that reduction)")
chi = lambda n: 0 if n%2==0 else (1 if n%4==1 else -1)
D = 4*mp.nsum(lambda j: chi(2*j+1)/((2*j+1)**2*(mp.exp(PI*(2*j+1))-1)), [0, mp.inf])
U = mp.nsum(lambda m: mp.sech(PI*m)/m**2, [1, mp.inf])
Lam = 2*G + D
print("  D =", mp.nstr(D,60)); print("  U =", mp.nstr(U,60)); print("  Lambda =", mp.nstr(Lam,60))
ck("addendum-5  S(4) = (2sqrt2/pi^2) K(1/sqrt2) Lambda", S4 - 2*r2/PI**2*Kl*Lam)
ck("BRIDGE      Lambda = pi^2/8 + 2V - T",               Lam - (PI**2/8 + 2*V - T))
ck("BRIDGE      2V = T - U/2 + pi^2/12",                 2*V - (T - U/2 + PI**2/12))
ck("addendum-5  Lambda = 5 pi^2/24 - U/2",               Lam - (5*PI**2/24 - U/2))

print("C5  the CM curve behind it: E32 = [y^2 = x^3 - x], CM by Z[i]")
from sympy import primerange
def ap_cm(p):
    if p == 2 or p % 4 == 3: return 0
    for b in range(0, int(mp.sqrt(p))+2):
        a2 = p-b*b; a = int(round(mp.sqrt(a2)))
        if a*a == a2 and a % 2 == 1:
            for aa in (a,-a):
                for bb in (b,-b):
                    if (aa+bb) % 4 == 1: return 2*aa
    raise ValueError(p)
NM = 800; av = [0]*(NM+1); av[1] = 1
for p in primerange(2, NM+1):
    ap = ap_cm(p); av[p] = ap; pk = p*p; pr, pr2 = ap, 1
    while pk <= NM:
        cur = ap*pr - (0 if p == 2 else p)*pr2; av[pk] = cur; pr2, pr = pr, cur; pk *= p
for n in range(2, NM+1):
    if av[n] == 0:
        for p in primerange(2, n+1):
            if n % p == 0:
                q = p
                while n % (q*p) == 0: q *= p
                m = n//q
                if m > 1: av[n] = av[q]*av[m]
                break
def Lval(s, N=32, eps=1):
    c = 2*PI/mp.sqrt(N); tot = mp.mpf(0)
    for n in range(1, NM+1):
        if av[n] == 0: continue
        x = c*n
        tot += av[n]*(x**(-s)*mp.gammainc(s,x) + eps*x**(-(2-s))*mp.gammainc(2-s,x))
    return tot*(2*PI)**s*mp.mpf(N)**(-mp.mpf(s)/2)/mp.gamma(s)
mp.mp.dps = 120
L1 = Lval(1); L2 = Lval(2); Lp0 = 32*L2/(4*PI**2)
ck("CONTROL  L(E32,1) = varpi/4  (BSD, rank 0)", L1 - vp/4, mp.mpf(10)**-95)
print("  L(E32,2)  =", mp.nstr(L2,50)); print("  L'(E32,0) =", mp.nstr(Lp0,50))
mp.mp.dps = 150

print("C6  bounded PSLQ exclusion for W = Lambda - pi^2/8, with live controls")
L = mp.log(1+r2); lg2 = mp.log(2); z3 = mp.zeta(3)
B = {'1':mp.mpf(1),'pi':PI,'pi^2':PI**2,'pi^3':PI**3,'1/pi':1/PI,'G':G,'G/pi':G/PI,'G*pi':G*PI,
     'L':L,'L^2':L**2,'L*pi':L*PI,'L/pi':L/PI,'log2':lg2,'log2^2':lg2**2,'log2*pi':lg2*PI,
     'r2':r2,'r2*pi':r2*PI,'r2/pi':r2/PI,'r2*G':r2*G,
     'varpi':vp,'varpi^2':vp**2,'varpi/pi':vp/PI,'pi/varpi':PI/vp,
     'varpi^2/pi^2':vp**2/PI**2,'pi^2/varpi^2':PI**2/vp**2,'pi^3/varpi^2':PI**3/vp**2,
     'varpi^2/pi':vp**2/PI,'pi/varpi^2':PI/vp**2,'G*varpi':G*vp,
     'z3':z3,'z3/pi^2':z3/PI**2,'z3/pi^3':z3/PI**3,
     'L2E32':L2,'L2E32*pi^2/varpi^2':L2*PI**2/vp**2,"L'(E32,0)":Lp0}
nm = list(B); vv = [B[k] for k in nm]
CORE = ['1','pi','pi^2','pi^3','G','varpi','varpi^2','G*varpi','G/pi','G*pi','pi^3/varpi^2',
        'varpi^2/pi','z3','L2E32',"L'(E32,0)",'L2E32*pi^2/varpi^2','L^2']
TOL = mp.mpf(10)**-130
def sweep(tgt, idx, maxsize, maxcoeff):
    hits = []
    for sz in range(1, maxsize+1):
        for c in itertools.combinations(idx, sz):
            r = mp.pslq([tgt]+[vv[i] for i in c], maxcoeff=maxcoeff, maxsteps=20000, tol=TOL)
            if r and r[0] != 0:
                rec = -sum(mp.mpf(r[i+1])*vv[j] for i,j in enumerate(c))/r[0]
                if abs(rec-tgt) < mp.mpf(10)**-140: hits.append((r,[nm[i] for i in c]))
    return hits
ALL = list(range(len(nm))); CIDX = [nm.index(x) for x in CORE]
print(f"  tier A: all {len(nm)} atoms, subsets <= 2, |c| <= 10^8, 150 dps, tol 1e-130")
print(f"  tier B: {len(CORE)}-atom core, subsets <= 3, |c| <= 10^6")
print(f"  core = {CORE}")
t0 = time.time()
for cn, cv in [('CONTROL Li2(1/2) = pi^2/12 - log2^2/2', PI**2/12 - lg2**2/2),
               ('CONTROL 7pi^3/180 = sum coth(pi n)/n^3', 7*PI**3/180),
               ('CONTROL varpi/4 = L(E32,1)', vp/4)]:
    hA = sweep(cv, ALL, 2, 10**8); hB = sweep(cv, CIDX, 3, 10**6)
    print(f"    {cn:40s}: A {len(hA)} / B {len(hB)}  {(hA+hB)[:1]}", flush=True)
for tn, tv in [('W = 2V - T', W), ('T', T), ('V', V), ('D', D), ('U', U),
               ('Lambda', Lam), ('S(4)', S4), ('pi^2 S(4)/varpi', PI**2*S4/vp)]:
    hA = sweep(tv, ALL, 2, 10**8); hB = sweep(tv, CIDX, 3, 10**6)
    print(f"    {tn:40s}: A {len(hA)} / B {len(hB)}  {(hA+hB)[:1]}", flush=True)
print(f"  sweep {time.time()-t0:.0f}s")
print("\n" + ("ALL S4-JACOBI-ROUTE REFEREE CHECKS PASSED" if all(o for _,o in OK)
               else "FAILURES: " + str([n for n,o in OK if not o])))
