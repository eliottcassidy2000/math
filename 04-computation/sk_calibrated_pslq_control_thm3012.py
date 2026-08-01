"""CALIBRATED search.  The 4-element basis {1,sqrt3} x {log(5+2sqrt6), atan(sqrt2/5)}
FINDS pi*S(3); the 40-element basis DOES NOT find the same relation at dps=240.
So wide-basis PSLQ nulls are meaningless.  Scan SMALL bases systematically instead.

Shape suggested by k=1,2,3:   pi*S(k) = alpha*L1 + q*L2,  alpha algebraic, q rational.
"""
import mpmath as mp, itertools, time
mp.mp.dps = 150

def S(k): return mp.hyper([mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/k], [1, 1+mp.mpf(1)/k], 1)

sq = lambda n: mp.sqrt(n)
MULT = {'1':mp.mpf(1),'r2':sq(2),'r3':sq(3),'r5':sq(5),'r6':sq(6),'r7':sq(7),
        'r10':sq(10),'r15':sq(15),'r30':sq(30),'r2/2':sq(2)/2,'r3/3':sq(3)/3}
LOG = {'log2':mp.log(2),'log3':mp.log(3),'log5':mp.log(5),'log7':mp.log(7),
       'log(1+r2)':mp.log(1+sq(2)),'log(2+r3)':mp.log(2+sq(3)),
       'log(5+2r6)':mp.log(5+2*sq(6)),'log(phi)':mp.log((1+sq(5))/2),
       'log(3+2r2)':mp.log(3+2*sq(2)),'log(9+4r5)':mp.log(9+4*sq(5)),
       'log(8+3r7)':mp.log(8+3*sq(7)),'log(2+r5)':mp.log(2+sq(5)),
       'atan(r2)':mp.atan(sq(2)),'atan(r2/5)':mp.atan(sq(2)/5),
       'atan(1/3)':mp.atan(mp.mpf(1)/3),'atan(1/2)':mp.atan(mp.mpf(1)/2),
       'atan(r5)':mp.atan(sq(5)),'atan(r3)':mp.atan(sq(3)),
       'pi':mp.pi,'Catalan':mp.catalan,'one':mp.mpf(1)}

def scan(target, tname, H=10**5, tol=130, size=2, budget=900):
    """all bases of `size` products alpha*L, plus optionally pi and 1"""
    prods = [(f"{m}*{l}", MULT[m]*LOG[l]) for m in MULT for l in LOG]
    t0 = time.time(); hits = []
    for combo in itertools.combinations(prods, size):
        if time.time() - t0 > budget: print("   (budget)"); break
        vals = [target] + [c[1] for c in combo]
        r = mp.pslq(vals, maxcoeff=H, maxsteps=200000, tol=mp.mpf(10)**(-tol))
        if r and r[0] != 0:
            hits.append((r, [c[0] for c in combo]))
    return hits

print("CONTROL: rediscover pi*S(3) by the SAME scan (size 2).")
h = scan(S(3)*mp.pi, 'pi*S(3)', size=2, budget=600)
print(f"   hits: {len(h)}")
for r, nm in h[:4]:
    print(f"      {r}  over {nm}")
