from math import gcd

def totient(n):
    return sum(1 for k in range(1,n+1) if gcd(k,n)==1)

def in_band(n, q, p, i):
    # {i*p/q} in [1/n, (n-1)/n]  <=>  q <= n*((i*p)%q) <= (n-1)*q  (matches Lean bandOK)
    r = (i*p) % q
    return q <= n*r <= (n-1)*q

def live_set(family, n, q):
    return [p for p in range(1, q) if all(in_band(n, q, p, v) for v in family)]

print("=== Tight AP {1,...,n-1}, threshold 1/n:  liveCount(q) =? phi(n)*[n|q] ===")
for n in range(3, 25):
    fam = list(range(1, n))
    ok = True; firstfail=None
    for q in range(2, 6*n+2):
        c = len(live_set(fam, n, q))
        expected = totient(n) if q % n == 0 else 0
        if c != expected and firstfail is None:
            ok=False; firstfail=(q,c,expected)
    l1 = live_set(fam, n, n)
    units = [o for o in range(1,n) if gcd(o,n)==1]
    print(f"n={n:2d} phi={totient(n):2d} law holds<=q={6*n+1}: {ok}"
          + (f"  FAIL q={firstfail[0]} got={firstfail[1]} exp={firstfail[2]}" if not ok else "")
          + f"   live@q=n: {l1} == units {units}? {l1==units}")
