from mpmath import mp, mpf, hyp3f2, pi, sqrt, log, atan, nstr, pslq
import itertools
mp.dps = 150
def S(k): return hyp3f2(mpf(1)/4, mpf(3)/4, mpf(1)/k, 1, 1+mpf(1)/k, 1)
s2, s3, s6 = sqrt(2), sqrt(3), sqrt(6)
C = {"pi": pi, "log(5+2r6)": log(5+2*s6), "r3 log(5+2r6)": s3*log(5+2*s6),
     "atan(r2/5)": atan(s2/5), "log(1+r2)": log(1+s2), "r3 log(1+r2)": s3*log(1+s2),
     "log(2+r3)": log(2+s3), "r3 log(2+r3)": s3*log(2+s3), "atan(r2)": atan(s2),
     "atan(r2/3)": atan(s2/3), "atan(r2/7)": atan(s2/7), "atan(r6/5)": atan(s6/5),
     "r2": s2, "r3": s3, "1": mpf(1), "log2": log(2), "r3 log2": s3*log(2),
     "log3": log(3), "r3 log3": s3*log(3), "atan(1/r2)": atan(1/s2)}
names = list(C)
for k in (6, 12):
    T = pi*S(k); hits = 0
    print(f"--- pi*S({k}) = {nstr(T,30)}")
    for size in (2, 3):
        for combo in itertools.combinations(names, size):
            r = pslq([T]+[C[n] for n in combo], tol=mpf(10)**-120,
                     maxcoeff=800, maxsteps=50000)
            if r and r[0] != 0:
                terms = " + ".join(f"{c}*{n}" for c, n in zip(r[1:], combo) if c)
                print(f"    {r[0]}*T + {terms} = 0"); hits += 1
                if hits > 4: break
        if hits > 4: break
    print(f"    relations: {hits}")
# sanity: the SAME machinery must re-find the known k=3 relation
T3 = pi*S(3); hits3 = 0
print("--- control: does the sweep re-find k=3?")
for size in (2, 3):
    for combo in itertools.combinations(names, size):
        r = pslq([T3]+[C[n] for n in combo], tol=mpf(10)**-120,
                 maxcoeff=800, maxsteps=50000)
        if r and r[0] != 0:
            terms = " + ".join(f"{c}*{n}" for c, n in zip(r[1:], combo) if c)
            print(f"    {r[0]}*T + {terms} = 0"); hits3 += 1
            if hits3 > 3: break
    if hits3 > 3: break
print(f"    k=3 relations found: {hits3} (control must be > 0)")
