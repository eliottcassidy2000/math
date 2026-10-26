# Exact test: is TT5(P_p) a polynomial (or a quasi-polynomial by residue class) in p?  Lagrange interpolation in
# exact rationals through some primes, prediction checked exactly on all the others.
import re
from fractions import Fraction
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT5=(\d+)', l) for l in open('paley_tt5_data.out')) if m]
def interp(pts):
    def f(x):
        tot = Fraction(0)
        for i, (xi, yi) in enumerate(pts):
            w = Fraction(yi)
            for j, (xj, _) in enumerate(pts):
                if j != i: w *= Fraction(x - xj, xi - xj)
            tot += w
        return tot
    return f
for deg in (5, 6):
    pts = rows[:deg + 1]
    f = interp(pts)
    misses = [(p, t, f(p)) for p, t in rows[deg + 1:] if f(p) != t]
    print(f"degree {deg} through the first {deg + 1} primes: {len(rows) - deg - 1 - len(misses)} of {len(rows) - deg - 1} other primes exact")
    if misses: print("   first misses:", [(p, t, float(fp)) for p, t, fp in misses[:3]])
for mod in (8, 3, 24):
    for r in sorted({p % mod for p, t in rows}):
        sub = [(p, t) for p, t in rows if p % mod == r]
        if len(sub) < 8: continue
        f = interp(sub[:6])
        ok = all(f(p) == t for p, t in sub[6:])
        print(f"  p = {r} mod {mod}: {len(sub)} primes; degree-5 interpolant through the first 6 exact on the rest: {ok}")
