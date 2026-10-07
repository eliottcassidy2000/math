import sympy as sp, itertools
a = [2, 40, 544, 6912, 87552, 1094144, 13534208, 165978112, 2022215680, 24520433664, 296325611520, 3572784594944, 43009509490688, 517207712333824]
ks = list(range(1, 15))
# try bases with optional k-polynomial factors
cands = []
for b in [12, 10, 9, 8, 6, 5, 4, 3, 2, 1, -2, -4]:
    for d in range(0, 3):
        cands.append((b, d))
def try_set(S):
    syms = sp.symbols('c0:%d' % len(S))
    eqs = []
    n = len(S)
    for idx in range(n):
        k = ks[idx]
        eqs.append(sum(syms[j] * (S[j][0] ** k) * (k ** S[j][1]) for j in range(n)) - a[idx])
    sol = sp.solve(eqs, syms, dict=True)
    if not sol: return None
    sol = sol[0]
    if len(sol) < n: return None
    for idx in range(n, len(ks)):
        k = ks[idx]
        if sum(sol[syms[j]] * (S[j][0] ** k) * (k ** S[j][1]) for j in range(n)) != a[idx]:
            return None
    return sol
found = False
for r in range(1, 6):
    for S in itertools.combinations(cands, r):
        if not any(s == (12, 0) for s in S): continue
        sol = try_set(list(S))
        if sol:
            print("FOUND", S, sol); found = True
    if found: break
if not found: print("no closed form with up to 5 terms of these shapes")
