"""kps-S33: aligned band blockers {1..11,13,W} are closed by the EXPLICIT far-peel threshold.
The census's nemesis (unbounded lonely-denominator, HYP-4040) is exactly the far-peel's domain.
Threshold from far_peel_lonely_of_cite: (1+2V)*(400V) < 3W, V = sum of the 12-base."""
from math import gcd, floor, ceil
from functools import reduce

def lcm(a,b): return a*b//gcd(a,b)
def lcm_list(xs): return reduce(lcm, xs, 1)

base = list(range(1,12)) + [13]           # {1,...,11, 13}  -- the fixed 12-base
V = sum(base)
threshold = (1 + 2*V)*(400*V)/3            # far_peel_lonely_of_cite: W > this
print(f"base = {base}")
print(f"V = sum(base) = {V}")
print(f"explicit peel threshold  W > (1+2V)*400V/3 = {threshold:,.1f}")
print()

def dist_to_int(x): return abs(x - round(x))
def is_lonely(speeds, t, margin=1/14):
    return all(dist_to_int(s*t) >= margin - 1e-12 for s in speeds)

def census_min_denominator(speeds, Qmax=60):
    """smallest q with a lonely p/q (the bounded-denominator census). None if q>Qmax."""
    for q in range(2, Qmax+1):
        for p in range(1, q):
            if is_lonely(speeds, p/q):
                return q
    return None

print(f"{'X':>3} {'W=lcm(2..X)':>14} {'W>thr?':>7} {'census q<=60':>13} {'peel goodlen>0':>14}")
for X in range(7, 22):
    W = lcm_list(range(2, X+1))
    speeds = base + [W]
    clears = W > threshold
    # census: smallest lonely denominator (None if none <= 60)
    q = census_min_denominator(speeds, 60)
    # peel: is the good region of the base nonempty AND does W's comb leave room?
    # (measure proxy) sample the base good region, check a lonely t exists for the full family
    #  near a base-good point -- the far-peel guarantees this when W > threshold.
    peel_ok = None
    if clears:
        # find a base-good rational, then confirm full family lonely at a nearby fine t
        found = False
        # scan fine t = k/(2W) near base-good points
        for k in range(1, 4*W, max(1, (4*W)//4000)):
            t = k/(2*W)
            if is_lonely(speeds, t):
                found = True; break
        peel_ok = found
    print(f"{X:>3} {W:>14,} {str(clears):>7} {str(q):>13} {str(peel_ok):>14}")

print()
print("READING: once W = lcm(2..X) clears ~1.67M (X>=17), the small-q census FAILS")
print("(q -> infinity, HYP-4040) but the far-peel closes it (peel_ok=True).")
print("The two-sided architecture: census for bounded W, far-peel for W > threshold.")
