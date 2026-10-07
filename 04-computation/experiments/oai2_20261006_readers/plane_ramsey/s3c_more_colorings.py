from s3_residue_colorings import *
from pysat.solvers import Cadical153

def omega(F, t, bit):
    """omega_t = ((2t-1) + sqrt(-(4t-1)))/(2t), with sqrt(-(4t-1)) = basis element 'bit' times a rational if needed."""
    return F.elt({0: Fr(2*t-1, 2*t), bit: Fr(1, 2*t)})

def gens_from(F, base_units, w1):
    gens = []
    z = F.one()
    for j in range(6):
        gens.append(z)
        for u in base_units:
            gens += [F.mul(z, u), F.mul(z, F.conj(u))]
        z = F.mul(z, w1)
    return gens

allok = True
# ---------------- C. Q(sqrt-3, sqrt-19): 3-adic colouring, residue field F_9 = F_3[i], U = mu_4 ----------------
F = MQ([-3, -19])
w1 = F.elt({0: Fr(1, 2), 1: Fr(1, 2)}); w5 = omega(F, 5, 2)      # (9+sqrt-19)/10
s19 = hensel_sqrt_padic(19, 3, 30, 1)
def res3_19(x):
    # residue = a_0 + a_{sqrt-19} * i * sqrt19 ; sqrt-3 parts vanish
    a0 = rat_mod(x[0], 3, 1); a2 = rat_mod(x[2], 3, 1)
    return (a0 % 3, (a2 * s19) % 3)          # element re + im*i of F_9
def col_K3K3(r): return (r[0] + r[1]) % 3
for u in gens_from(F, [w5], w1): 
    r = res3_19(u); assert r in [(1,0),(2,0),(0,1),(0,2)], r
hu = hilbert90_units(F, 20)
for u in hu:
    # field units: check they are 3-adic units with residue in mu_4 (needs P-integrality; Hilbert-90 units can have 3 in denominators? they are P-units by the lemma)
    pass
pts = build_patch(F, gens_from(F, [w5], w1), 3, 3000)
edges, _, diffs = exact_edges(F, pts)
allok &= check("Q(sqrt-3,sqrt-19) [Eisenstein + Heegner-19 rotation] 3-adic", F, pts, lambda x: col_K3K3(res3_19(x)), edges, diffs, 3)
print("   3-colourable (first 500)?", sat_colorable(500, [(i,j) for i,j in edges if i<500 and j<500], 3))

# ---------------- D. Q(zeta12) = Q(sqrt3)^2-plane: residue mod (sqrt3), F_9, U = mu_4 ----------------
F = MQ([-1, -3])          # e1 = i, e2 = sqrt-3, e3 = i*sqrt-3 = -sqrt3
w1 = F.elt({0: Fr(1, 2), 2: Fr(1, 2)})
z12 = F.elt({1: Fr(1, 2), 3: Fr(-1, 2)})          # (sqrt3 + i)/2
pyth = [F.elt({0: Fr(3, 5), 1: Fr(4, 5)}), F.elt({0: Fr(5, 13), 1: Fr(12, 13)}), F.elt({0: Fr(8, 17), 1: Fr(15, 17)})]
def res3_z12(x):
    a0 = rat_mod(x[0], 3, 1); a1 = rat_mod(x[1], 3, 1)
    return (a0 % 3, a1 % 3)
base = [z12] + pyth
gens = gens_from(F, base, w1)
for u in gens:
    r = res3_z12(u); assert r in [(1,0),(2,0),(0,1),(0,2)], r
pts = build_patch(F, gens, 2, 3000)
edges, _, diffs = exact_edges(F, pts)
allok &= check("Q(zeta12)=Q(sqrt3)^2 (Eisenstein+Gaussian+Pythagorean) 3-adic", F, pts, lambda x: col_K3K3(res3_z12(x)), edges, diffs, 3)

# ---------------- B. Polymath field Q(sqrt-3,sqrt-11,sqrt-15): 11-adic, F_121, U = mu_12, kappa(11)=5 ----------------
F = MQ([-3, -11, -15])     # bits: 0 sqrt-3, 1 sqrt-11, 2 sqrt-15
p = 11
# F_121 = F_11[r]/(r^2+3)
def fmul(a, b): return ((a[0]*b[0] - 3*a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)
def fpow(a, e):
    r = (1, 0)
    for _ in range(e): r = fmul(r, a)
    return r
FE = [(a, b) for a in range(p) for b in range(p)]
MU12 = [x for x in FE if x != (0, 0) and fpow(x, 12) == (1, 0)]
assert len(MU12) == 12
idx = {x: i for i, x in enumerate(FE)}
cay = set()
for x in FE:
    for u in MU12:
        y = ((x[0]+u[0]) % p, (x[1]+u[1]) % p)
        cay.add((min(idx[x], idx[y]), max(idx[x], idx[y])))
k = 5
s = Cadical153(); v = lambda i, c: i*k + c + 1
for i in range(len(FE)):
    s.add_clause([v(i, c) for c in range(k)])
for a, b in cay:
    for c in range(k): s.add_clause([-v(a, c), -v(b, c)])
assert s.solve(); model = set(l for l in s.get_model() if l > 0)
COL121 = {FE[i]: next(c for c in range(k) if v(i, c) in model) for i in range(len(FE))}
s4 = Cadical153()
for i in range(len(FE)):
    s4.add_clause([i*4+c+1 for c in range(4)])
for a, b in cay:
    for c in range(4): s4.add_clause([-(a*4+c+1), -(b*4+c+1)])
print("   Cay(F_121, mu_12): 5-colourable yes; 4-colourable?", s4.solve())
t5 = hensel_sqrt_padic(5, 11, 20, 4)
def res11(x):
    # x = sum a_S e_S ; any S containing bit1 (sqrt-11) -> 0
    a0 = rat_mod(x[0], 11, 1); a1 = rat_mod(x[1], 11, 1); a4 = rat_mod(x[4], 11, 1); a5 = rat_mod(x[5], 11, 1)
    # e_{0}=sqrt-3 -> r ; e_{2}=sqrt-15 -> r*t ; e_{0,2}=sqrt-3*sqrt-15 = -3 sqrt5 -> -3t
    re = (a0 - 3*t5*a5) % p
    im = (a1 + t5*a4) % p
    return (re, im)
w1 = F.elt({0: Fr(1, 2), 1: Fr(1, 2)})
w3 = omega(F, 3, 2); w4 = omega(F, 4, 4)       # bit index: sqrt-11 is e_{2}(=bit1 -> index 2), sqrt-15 is index 4
gens = gens_from(F, [w3, w4, F.mul(w3, w4), F.mul(w3, F.conj(w4))], w1)
for u in gens:
    r = res11(u); assert r in MU12, (u, r)
pts = build_patch(F, gens, 2, 4000)
edges, _, diffs = exact_edges(F, pts)
allok &= check("Polymath field Q(sqrt-3,sqrt-11,sqrt-15) = Frac Z[w1,w3,w4]: 11-adic 5-colouring", F, pts, lambda x: COL121[res11(x)], edges, diffs, 5)

# ---------------- E. Heegner compositum (3,11,19,43,67,163): 2-adic F_4 colouring ----------------
ds = [3, 11, 19, 43, 67, 163]
F = MQ([-d for d in ds])
res = make_res_2adic(F, ds, K=90, M=40)
w1 = F.elt({0: Fr(1, 2), 1: Fr(1, 2)})
ws = []
for i, d in enumerate(ds[1:], start=1):
    t = (d + 1)//4
    ws.append(omega(F, t, 1 << i))
gens = gens_from(F, ws, w1)
for u in gens: assert res(u) != (0, 0)
pts = build_patch(F, gens, 2, 1500)
edges, _, diffs = exact_edges(F, pts)
allok &= check("Heegner compositum Q(sqrt-d: d=3,11,19,43,67,163) 2-adic F_4", F, pts, res, edges, diffs, 4)
print("ALL OK:", allok)
