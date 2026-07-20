"""opus-2026-07-20-S422 (HYP-8190): working Targets 1/2/6 concurrently.

T6: unit-distance census -- generate G_4, G_5 arc-flip metagraphs from scratch,
    K4-freeness, then numeric unit-distance embedding (Gauss-Newton, validated on
    the Moser spindle control).  "Even a clean negative upgrades the cards."
T1: (C1) the torus-compressed degree-<=3 presentation of F in dim 6 via gadget
    vars s=xy, w=x^2 z, t=s^2 -- SYMBOLIC det test: does the naive substitution
    gadget preserve Keller-ness?  (Either answer is a result: constant -2 gives a
    dim-6 cubic-degree Keller presentation; non-constant DEMONSTRATES the naive
    gadget fails and the BCW/dBvdE constructions must be pinned from the text --
    correcting S421's implication that the transport is mechanical.)
    (C2) toy validation of the symmetric/HN stage: P = y1*x2^3 quartic; Hessian
    nilpotency; Zhao vanishing Delta^m(P^m) = 0; Delta^m(P^(m+1)) behavior.
"""
import itertools, numpy as np, sympy as sp
from math import comb

# ================= T6: metagraphs and unit distance =================
def metagraph(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    perms = list(itertools.permutations(range(n)))
    def canon(mask):
        best = None
        for p in perms:
            v = 0
            for b, (i, j) in enumerate(pairs):
                pi, pj = p[i], p[j]
                if pi < pj:
                    bit = (mask >> pairs.index((pi, pj))) & 1
                else:
                    bit = 1 - ((mask >> pairs.index((pj, pi))) & 1)
                v |= bit << b
            if best is None or v < best:
                best = v
        return best
    canon_cache = {}
    def C(mask):
        if mask not in canon_cache:
            canon_cache[mask] = canon(mask)
        return canon_cache[mask]
    classes = {}
    for mask in range(1 << m):
        c = C(mask)
        classes.setdefault(c, []).append(mask)
    idx = {c: i for i, c in enumerate(sorted(classes))}
    E = set()
    for mask in range(1 << m):
        a = idx[C(mask)]
        for b in range(m):
            nb = idx[C(mask ^ (1 << b))]
            if nb != a:
                E.add((min(a, nb), max(a, nb)))
    return len(idx), sorted(E)

def has_k4(V, E):
    adj = [[False] * V for _ in range(V)]
    for a, b in E:
        adj[a][b] = adj[b][a] = True
    for q in itertools.combinations(range(V), 4):
        if all(adj[a][b] for a, b in itertools.combinations(q, 2)):
            return True
    return False

def embed_unit(V, E, tries=300, iters=400, seed=1):
    rng = np.random.default_rng(seed)
    best = 1e18
    for _ in range(tries):
        p = rng.normal(size=(V, 2)) * 1.5
        for _ in range(iters):
            # Gauss-Newton on edge residuals
            r = []
            J = np.zeros((len(E), 2 * V))
            for k, (a, b) in enumerate(E):
                d = p[a] - p[b]
                q = d @ d
                r.append(q - 1.0)
                J[k, 2 * a:2 * a + 2] = 2 * d
                J[k, 2 * b:2 * b + 2] = -2 * d
            r = np.array(r)
            if np.abs(r).max() < 1e-12:
                break
            try:
                step = np.linalg.lstsq(J, -r, rcond=None)[0]
            except np.linalg.LinAlgError:
                break
            p = p + 0.8 * step.reshape(V, 2)
        res = max(abs((p[a] - p[b]) @ (p[a] - p[b]) - 1) for a, b in E)
        mind = min(np.hypot(*(p[a] - p[b])) for a in range(V) for b in range(a + 1, V))
        if res < best:
            best, bestp, bestmind = res, p.copy(), mind
        if best < 1e-10 and bestmind > 1e-4:
            break
    return best, bestmind

print("=" * 80)
print("T6: UNIT-DISTANCE CENSUS OF THE ARC-FLIP METAGRAPHS")
print("=" * 80)
# Moser spindle control (7 vertices, 11 edges, known unit-distance, chi=4)
SP_E = [(0,1),(0,2),(1,2),(1,3),(2,3),(3,4),(4,5),(4,6),(5,6),(0,5),(0,6)]
res, mind = embed_unit(7, SP_E, tries=150)
print(f"  control Moser spindle: max|d^2-1| = {res:.2e}, min pair dist = {mind:.3f} "
      f"-> embedder {'VALID' if res < 1e-9 and mind > 1e-3 else 'FAILED'}")
for n in (4, 5):
    V, E = metagraph(n)
    k4 = has_k4(V, E)
    print(f"  G_{n}: V = {V}, E = {len(E)}, K4-free: {not k4}", end="")
    if k4:
        print("  -> NOT unit-distance (K4 obstruction)")
    else:
        res, mind = embed_unit(V, E)
        verdict = "YES realizable" if res < 1e-9 and mind > 1e-3 else f"no embedding found (residual {res:.1e})"
        print(f"  -> numeric embed: {verdict} (min pair dist {mind:.3f})")

# ================= T1 (C1): the dim-6 gadget test =================
print()
print("=" * 80)
print("T1/C1: TORUS-COMPRESSED DEGREE-<=3 PRESENTATION IN DIM 6 -- KELLER OR NOT?")
print("=" * 80)
x, y, z, s, w, t = sp.symbols('x y z s w t')
F1h = z + 3*z*s + 3*z*t + z*s*t + 4*y**2 + 7*y**2*s + 3*y**2*t
F2h = y + 3*x*z + 6*x*z*s + 3*x*z*t + 12*s*y + 9*y*t
F3h = 2*x - 3*x*s - x*w
G = sp.Matrix([F1h, F2h, F3h, s - x*y, w - x**2*z, t - s**2])
degs = [sp.total_degree(sp.Poly(g, x, y, z, s, w, t)) for g in G]
print(f"  component degrees: {degs} (all <= 3: {max(degs) <= 3})")
# section check: G restricted to s=xy, w=x^2 z, t=(xy)^2 reproduces (F, 0, 0, 0)
sub = {s: x*y, w: x**2*z, t: (x*y)**2}
u_ = 1 + x*y
F_orig = [u_**3*z + y**2*u_*(4+3*x*y),
          y + 3*x*u_**2*z + 3*x*y**2*(4+3*x*y),
          2*x - 3*x**2*y - x**3*z]
sec = [sp.expand(G[i].subs(sub) - F_orig[i]) for i in range(3)] + \
      [sp.expand(G[i].subs(sub)) for i in range(3, 6)]
print(f"  section reproduces F and kills gadget rows: {all(e == 0 for e in sec)}")
detG = sp.expand(sp.Matrix(G).jacobian([x, y, z, s, w, t]).det())
isconst = detG.free_symbols == set()
print(f"  det J(G) constant? {isconst}")
if isconst:
    print(f"  det = {detG}  -> DIM-6 DEGREE-3 KELLER PRESENTATION ACHIEVED")
else:
    p = sp.Poly(detG, x, y, z, s, w, t)
    print(f"  det J(G) = {sp.factor(detG.subs({w: 0, t: 0}))} ... (truncated view at w=t=0)")
    print(f"  #monomials in det: {len(p.terms())}; det on the section s=xy,w=x^2z,t=s^2:"
          f" {sp.expand(detG.subs(sub))}")
    print("  -> NAIVE SUBSTITUTION GADGET DOES NOT PRESERVE KELLER-NESS (demonstrated).")
    print("     The BCW/dBvdE reductions are genuinely non-naive; S421's 'pin the dimension'")
    print("     understated this: pinning the CONSTRUCTION from van den Essen Ch. 6 is the task.")

# ================= T1 (C2): toy HN / Zhao-Delta machinery =================
print()
print("=" * 80)
print("T1/C2: TOY SYMMETRIC-QUARTIC STAGE -- HESSIAN NILPOTENCY + ZHAO DELTA MACHINERY")
print("=" * 80)
v = sp.symbols('v1 v2 v3 v4')          # (x1, x2, y1, y2)
P = v[2] * v[1]**3                      # <y, H(x)> with H = (x2^3, 0): quartic
H = sp.hessian(P, v)
Hk = H
nilp = None
for k in range(1, 5):
    if Hk == sp.zeros(4, 4):
        nilp = k
        break
    Hk = sp.expand(Hk * H)
print(f"  P = y1*x2^3: Hessian nilpotent of index {nilp} (HN: {nilp is not None})")
Gm = sp.Matrix([v[i] + sp.diff(P, v[i]) for i in range(4)])
dg = sp.expand(Gm.jacobian(v).det())
print(f"  det J(id + grad P) = {dg} (Keller: {dg.free_symbols == set()})")
def lap(f):
    return sp.expand(sum(sp.diff(f, vi, 2) for vi in v))
def lappow(f, m):
    for _ in range(m):
        f = lap(f)
    return f
vanish_mm = [sp.expand(lappow(P**m, m)) == 0 for m in range(1, 5)]
van_m1 = [sp.expand(lappow(P**(m+1), m)) for m in range(1, 7)]
print(f"  Zhao vanishing Delta^m(P^m) = 0 for m = 1..4: {vanish_mm}")
print(f"  Delta^m(P^(m+1)) zero-pattern m = 1..6: {[e == 0 for e in van_m1]}")
print("  (invertible toy: eventual vanishing expected -- this is the machinery the")
print("   transported F-witness must BREAK at infinitely many m)")
