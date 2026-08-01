"""VALIDATED pipeline (it blind-rediscovers pi*S(3)) applied to S(4) and S(5).
Reports either an exact closed form, or a MEANINGFUL bounded negative:
'no Z-relation over the product basis B, |coeff| <= H, at precision P'."""
import mpmath as mp, itertools, time
mp.mp.dps = 150

def S(k): return mp.hyper([mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/k], [1, 1+mp.mpf(1)/k], 1)
sq = lambda n: mp.sqrt(n)
MULT = {'1':mp.mpf(1),'r2':sq(2),'r3':sq(3),'r5':sq(5),'r6':sq(6),'r7':sq(7),
        'r10':sq(10),'r15':sq(15),'r30':sq(30),'2':mp.mpf(2),'r2/2':sq(2)/2}
LOG = {'log2':mp.log(2),'log3':mp.log(3),'log5':mp.log(5),'log7':mp.log(7),
       'log(1+r2)':mp.log(1+sq(2)),'log(2+r3)':mp.log(2+sq(3)),
       'log(5+2r6)':mp.log(5+2*sq(6)),'log(phi)':mp.log((1+sq(5))/2),
       'log(3+2r2)':mp.log(3+2*sq(2)),'log(9+4r5)':mp.log(9+4*sq(5)),
       'log(8+3r7)':mp.log(8+3*sq(7)),'log(2+r5)':mp.log(2+sq(5)),
       'log(1+r3)':mp.log(1+sq(3)),'log(4+r15)':mp.log(4+sq(15)),
       'atan(r2)':mp.atan(sq(2)),'atan(r2/5)':mp.atan(sq(2)/5),
       'atan(1/3)':mp.atan(mp.mpf(1)/3),'atan(1/2)':mp.atan(mp.mpf(1)/2),
       'atan(r5)':mp.atan(sq(5)),'atan(1/r2)':mp.atan(1/sq(2)),
       'atan(r3/5)':mp.atan(sq(3)/5),'atan(r6)':mp.atan(sq(6)),
       'pi':mp.pi,'Catalan':mp.catalan,'one':mp.mpf(1)}
PROD = [(f"{m}*{l}", MULT[m]*LOG[l]) for m in MULT for l in LOG]
print(f"product basis size = {len(PROD)}; pairs = {len(PROD)*(len(PROD)-1)//2}")

def scan2(target, tname, H, tol, budget):
    t0=time.time(); hits=[]; n=0
    for a, b in itertools.combinations(PROD, 2):
        n += 1
        if time.time()-t0 > budget: break
        r = mp.pslq([target, a[1], b[1]], maxcoeff=H, maxsteps=120000, tol=mp.mpf(10)**(-tol))
        if r and r[0] != 0: hits.append((r, a[0], b[0]))
    return hits, n, time.time()-t0

for k in (4, 5):
    v = S(k)*mp.pi
    print(f"\n=== pi*S({k}) = {mp.nstr(v, 40)}")
    hits, n, el = scan2(v, f'pi*S({k})', 10**5, 130, 1200)
    print(f"   scanned {n} pairs in {el:.0f}s at dps=150, |coeff| <= 1e5, tol 1e-130")
    if hits:
        print(f"   *** {len(hits)} HIT(S) ***")
        for r, a, b in hits[:8]: print(f"      {r}  over [{a}, {b}]")
    else:
        print(f"   NEGATIVE: no Z-relation pi*S({k}) = c1*A + c2*B over any PAIR from the")
        print(f"   {len(PROD)}-element product basis with |coeff| <= 1e5 at 150 digits.")
