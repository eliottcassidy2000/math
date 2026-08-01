"""WEIGHT HYPOTHESIS.  pi*S(1) = 8sqrt2/3 (weight 0), pi*S(2) = 4log(1+sqrt2) and
pi*S(3) = sqrt3 log(5+2sqrt6) - 2atan(sqrt2/5) (weight 1).  Since S(4) = (2sqrt2/pi)
int_0^1 K/(2-k^2) dk and int_0^1 K dk = 2G, S(4) plausibly sits at WEIGHT 2
(Catalan, Li_2, pi^2, log^2).  A pure log/arctan basis would then miss it by
construction -- consistent with every prior negative.  TEST IT."""
import mpmath as mp, itertools, time
mp.mp.dps = 130
def S(k): return mp.hyper([mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/k], [1, 1+mp.mpf(1)/k], 1)
sq = lambda n: mp.sqrt(n)
G = mp.catalan; PI = mp.pi
l2 = mp.log(1+sq(2)); l3 = mp.log(2+sq(3)); l6 = mp.log(5+2*sq(6))

W2 = {  # weight-2 constants
 'G':G, 'pi^2':PI**2, 'log2^2':mp.log(2)**2, 'l2^2':l2**2, 'l3^2':l3**2, 'l6^2':l6**2,
 'pi*l2':PI*l2, 'pi*l3':PI*l3, 'pi*l6':PI*l6, 'pi*log2':PI*mp.log(2),
 'Li2(1/2)':mp.polylog(2, mp.mpf(1)/2), 'Li2(r2-1)':mp.polylog(2, sq(2)-1),
 'Li2(-r2+1)':mp.polylog(2, 1-sq(2)), 'Li2(1/4)':mp.polylog(2, mp.mpf(1)/4),
 'Li2(3-2r2)':mp.polylog(2, 3-2*sq(2)), 'log2*l2':mp.log(2)*l2,
}
MULT = {'1':mp.mpf(1),'r2':sq(2),'r3':sq(3),'r6':sq(6)}
PROD = [(f"{m}*{w}", MULT[m]*W2[w]) for m in MULT for w in W2]
print(f"weight-2 product basis: {len(PROD)} elements")

for name, tgt in [('S(4)', S(4)), ('pi*S(4)', PI*S(4)), ('pi^2*S(4)', PI**2*S(4)),
                  ('S(5)', S(5)), ('pi*S(5)', PI*S(5))]:
    t0=time.time(); hits=[]
    for a, b in itertools.combinations(PROD, 2):
        if time.time()-t0 > 260: break
        r = mp.pslq([tgt, a[1], b[1]], maxcoeff=10**5, maxsteps=100000, tol=mp.mpf(10)**-115)
        if r and r[0] != 0: hits.append((r, a[0], b[0]))
    # also singletons
    for a in PROD:
        r = mp.pslq([tgt, a[1]], maxcoeff=10**6, maxsteps=100000, tol=mp.mpf(10)**-115)
        if r and r[0] != 0: hits.append((r, a[0], '-'))
    print(f"  {name:10s}: {len(hits)} hit(s)" + (f"  {hits[:3]}" if hits else "  (none)"), flush=True)

print("\nCONTROL: the same weight-2 machinery must NOT spuriously fire on pi*S(2)")
t0=time.time(); h=[]
for a, b in itertools.combinations(PROD, 2):
    if time.time()-t0 > 120: break
    r = mp.pslq([PI*S(2), a[1], b[1]], maxcoeff=10**5, maxsteps=100000, tol=mp.mpf(10)**-115)
    if r and r[0] != 0: h.append((r,a[0],b[0]))
print(f"  pi*S(2) against weight-2 basis: {len(h)} hit(s)  (expect 0 -- it is weight 1)")
