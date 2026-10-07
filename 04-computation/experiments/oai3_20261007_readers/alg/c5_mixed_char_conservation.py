"""C5: mixed-characteristic conservation of number for THM-1300's map at p = 2.

For P in Z_2^3, w = F(P), the fibre scheme Z = Spec Z_2[x,y,z]/(F - w) is a complete
intersection of dimension 1 in the 4-dimensional (unramified) regular local ring
R = Z_2[x,y,z]_(2, x-xb, y-yb, z-zb) at the closed point Pb = P mod 2.  If the mod-2 fibre
F_bar^{-1}(w_bar) is finite at Pb, then O_Z is CM, 2-torsion free, and Serre's
chi^R(O_Z, R/2R) = length(O_Z/2) = mu_bar(Pb) (local length of the mod-2 fibre),
which equals the number of Qbar_2-points of F^{-1}(w) in the residue disc of Pb.

Here: (a) mu_bar(Pb) and the fibre dimension at each Pb in F_2^3 (Groebner over GF(2));
(b) the distribution of N_disc(P) = #(Q_2-rational integral preimages of F(P) congruent to
P mod 2), exact from the mod-8 census (Hensel: v_2(det JF) = 1, so mod 8 decides).
"""
import itertools
import numpy as np
import sympy as sp

x, y, z = sp.symbols('x y z')
u = 1 + x*y
F = [sp.expand(u**3*z + y**2*u*(4 + 3*x*y)),
     sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)),
     sp.expand(2*x - 3*x**2*y - x**3*z)]

def Fint(P, m):
    X, Y, Z = P
    U = 1 + X*Y
    return ((U**3*Z + Y**2*U*(4 + 3*X*Y)) % m, (Y + 3*X*U**2*Z + 3*X*Y**2*(4 + 3*X*Y)) % m,
            (2*X - 3*X**2*Y - X**3*Z) % m)

def std_count(G, gens, B=14):
    lms = [sp.Poly(g, *gens, modulus=2).monoms(order='grevlex')[0] for g in G.exprs]
    cnt = 0
    for e in itertools.product(range(B), repeat=3):
        if not any(all(e[i] >= l[i] for i in range(3)) for l in lms):
            cnt += 1
    return cnt  # == B^3 - ... ; if it equals the cap behaviour, the ideal is not 0-dimensional

print("(a) mod-2 fibres at the points of F_2^3")
local = {}
for Pb in itertools.product(range(2), repeat=3):
    wb = Fint(Pb, 2)
    # translate Pb to the origin: x -> x + xb etc.
    sub = {x: x + Pb[0], y: y + Pb[1], z: z + Pb[2]}
    I = [sp.expand(f.subs(sub) - c) for f, c in zip(F, wb)]
    G = sp.groebner(I, x, y, z, order='grevlex', modulus=2)
    zero_dim = G.is_zero_dimensional
    out = "positive-dimensional"
    mu = None
    if zero_dim:
        # local length at the origin: dim F_2[x,y,z]/(I + m^K), K large (stabilises)
        vals = []
        for K in (6, 10, 14):
            mK = [x**a*y**b*z**c for a in range(K + 1) for b in range(K + 1) for c in range(K + 1) if a + b + c == K]
            GK = sp.groebner(list(G.exprs) + mK, x, y, z, order='grevlex', modulus=2)
            vals.append(std_count(GK, (x, y, z), B=K + 2))
        mu = vals[-1]
        out = "0-dim; global fibre length over F_2-bar = %d; local length at Pb (K=6,10,14) = %s" % (
            std_count(G, (x, y, z), B=20), vals)
    else:
        # local dimension: is the origin on a positive-dimensional component?
        out = "positive-dimensional fibre"
    local[Pb] = mu
    print("   Pb =", Pb, " w_bar =", wb, ":", out)

print("\n(b) exact disc counts from the mod-8 census (N_disc = # integral Q_2-preimages = P mod 2)")
m = 8
pts = list(itertools.product(range(m), repeat=3))
img = {}
for P in pts:
    img.setdefault(Fint(P, m), []).append(P)
for Pb in itertools.product(range(2), repeat=3):
    hist = {}
    for P in pts:
        if tuple(c % 2 for c in P) != Pb:
            continue
        fib = img[Fint(P, m)]
        # each integral preimage occupies 2 classes mod 8 (|det| = 1/2), distinct preimages are distinct mod 4
        same_disc = sum(1 for Q in fib if tuple(c % 2 for c in Q) == Pb) // 2
        total = len(fib) // 2
        key = (same_disc, total)
        hist[key] = hist.get(key, 0) + 1
    print("   Pb =", Pb, " mu_bar =", local[Pb], "  (N_disc, N_total): share of the disc ->",
          {k: "%d/64" % v for k, v in sorted(hist.items())})
