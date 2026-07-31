"""DECISIVE TEST 1 of HYP-9070: is K = 1 (leading form a pure power of ONE
linear form) forced for polynomial automorphisms of the affine plane?

Build automorphisms as composites of affine maps and triangular maps
(x,y) -> (x, y + p(x)) / (x + p(y), y), verify Jac = const, then factor the
leading forms P_n, Q_m and count DISTINCT linear factors K."""
import random
import sympy as sp

x, y = sp.symbols('x y')


def jac(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))


def hom_top(F):
    F = sp.Poly(sp.expand(F), x, y)
    d = F.total_degree()
    return sp.expand(sum(c * x ** a * y ** b for (a, b), c in F.terms() if a + b == d)), d


def K_of(form):
    """number of DISTINCT linear factors of a binary form"""
    if sp.simplify(form) == 0:
        return None, None
    fl = sp.factor_list(form)
    facs = [f for f, e in fl[1] if sp.total_degree(f) > 0]
    return len(facs), [(sp.factor(f), e) for f, e in fl[1] if sp.total_degree(f) > 0]


def rand_affine(rng):
    while True:
        a, b, c, d = [rng.randint(-3, 3) for _ in range(4)]
        if a * d - b * c != 0:
            return (a * x + b * y + rng.randint(-2, 2), c * x + d * y + rng.randint(-2, 2))


def rand_tri(rng, deg):
    coeffs = [rng.randint(-3, 3) for _ in range(deg + 1)]
    coeffs[-1] = rng.choice([1, -1, 2])
    if rng.random() < 0.5:
        p = sum(coeffs[i] * x ** i for i in range(deg + 1))
        return (x, sp.expand(y + p))
    p = sum(coeffs[i] * y ** i for i in range(deg + 1))
    return (sp.expand(x + p), y)


def compose(F, G):
    """F after G"""
    return (sp.expand(F[0].subs({x: G[0], y: G[1]}, simultaneous=True)),
            sp.expand(F[1].subs({x: G[0], y: G[1]}, simultaneous=True)))


rng = random.Random(9070)
rows = []
for trial in range(10):
    F = (x, y)
    for _ in range(rng.randint(1, 2)):
        F = compose(rand_affine(rng), F)
        F = compose(rand_tri(rng, rng.randint(2, 3)), F)
    F = compose(rand_affine(rng), F)
    J = sp.simplify(jac(F[0], F[1]))
    if not J.is_number or J == 0:
        continue
    Pn, n = hom_top(F[0]); Qm, m = hom_top(F[1])
    if n < 2 or m < 2:
        continue
    Kp, fp = K_of(Pn); Kq, fq = K_of(Qm)
    g = sp.gcd(n, m)
    rows.append((n, m, g, n // g, m // g, Kp, Kq))
    if len(rows) <= 8:
        print(f"  deg=({n},{m}) g={g} (a,b)=({n//g},{m//g}) Jac={J}  "
              f"K(P_n)={Kp} K(Q_m)={Kq}   P_n={sp.factor(Pn)}")

print(f"\n{len(rows)} automorphisms sampled")
allK = [(kp, kq) for *_, kp, kq in rows]
print("  distinct (K(P_n),K(Q_m)) observed:", sorted(set(allK)))
print("  K=1 for BOTH in all samples:", all(kp == 1 and kq == 1 for kp, kq in allK))
ab = sorted(set((a, b) for _, _, _, a, b, _, _ in rows))
print("  (a,b) observed:", ab)
print("  every sample has a=1 or b=1 (Jung-van der Kulk):",
      all(a == 1 or b == 1 for a, b in ab))
