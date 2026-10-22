# The Camion residue as a function of the base residue: rho_C(a) = (sigma_p(a) - 1)/p mod p, a in F_p^x (lift-independent).
# Look for symmetries: a -> 1/a, a -> -a, a -> a^k; and whether rho_C(a) is a low-degree polynomial in a (Lagrange degree).
# mac-mini-2026-10-21-fanopsp.
import importlib.util, io, contextlib
spec = importlib.util.spec_from_file_location("cl", "camion_lift.py")
buf = io.StringIO()
with contextlib.redirect_stdout(buf): cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
def lagrange_degree(vals, p):
    # degree of the interpolating polynomial over F_p through the points (a, vals[a]), a = 1..p-1 (value at 0 unknown: try all)
    best = None
    for v0 in range(p):
        pts = [(0, v0)] + [(a, vals[a]) for a in range(1, p)]
        # Newton divided differences mod p
        xs = [x for x, _ in pts]; ys = [y % p for _, y in pts]
        coef = ys[:]
        for j in range(1, p):
            for i in range(p - 1, j - 1, -1):
                coef[i] = (coef[i] - coef[i - 1]) * pow((xs[i] - xs[i - j]) % p, -1, p) % p
        deg = max([i for i in range(p) if coef[i] % p] or [0])
        best = deg if best is None else min(best, deg)
    return best
for p in (11, 13, 17, 19):
    cl.K = 2
    rho = {a: (cl.sigma_lim(a, p, 2)[0] - 1) // p % p for a in range(1, p)}
    inv = all(rho[a] == rho[pow(a, -1, p)] for a in range(1, p))
    neg = all(rho[a] == rho[p - a] for a in range(1, p))
    antineg = all((rho[a] + rho[p - a]) % p == 0 for a in range(1, p))
    print(f"p={p}: rho_C(a) for a=1..{p-1}: {[rho[a] for a in range(1, p)]}")
    print(f"   rho(1/a)=rho(a): {inv}  rho(-a)=rho(a): {neg}  rho(-a)=-rho(a): {antineg}  min Lagrange degree (value at 0 free): {lagrange_degree(rho, p)} (generic p-2 = {p-2})")
