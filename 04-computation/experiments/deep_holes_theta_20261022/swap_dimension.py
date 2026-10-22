"""Is the swap defect forced by representation theory? (session opus-2026-10-22-S53)

Setting.
- D = (Z/h)^2 is the hyperbolic discriminant form, with q(u, v) = uv/h and b((u,v),(u',v')) = (uv' + u'v)/h.
- rho_D is the Weil representation of signature 24 (Borcherds' convention):
    rho(T) e_g = e(q(g)) e_g,
    rho(S) e_g = (1/h) sum_d e(-b(g, d)) e_d.
- The torus t_a(u, v) = (a u, v/a), a in (Z/h)^*, and the swap sigma(u, v) = (v, u) both commute with rho_D.
- V^- is the space of vectors that are torus-invariant and sigma-antisymmetric.
- The antisymmetric invariant vectors are 1_{S_d} - 1_{S_{h/d}}, where S_d = dZ x (h/d)Z runs over the maximal
  isotropic subgroups; d and h/d are paired by sigma.

Swap theorem (would follow). If dim S_12(rho on V^-) equals the number of these pairs, then every cusp form in
S_12(V^-) is Delta times an invariant vector. In particular E = Theta_M - sigma^* Theta_M is such a form. With the
torus symmetry of Theta_M this proves HYP-9363.

Dimension formula (Borcherds, Reflection groups of Lorentzian lattices, sec. 7; Skoruppa). For k > 2, with
alpha(A) the sum of beta in [0,1) over the eigenvalues e(beta) of A:

    dim M_k(rho) = d + dk/12 - alpha(e^(pi i k/2) rho(S)) - alpha((e^(pi i k/3) rho(ST))^(-1)) - alpha(rho(T)),
    dim S_k = dim M_k - dim ker(rho(T) - 1).

The implementation is validated on the trivial representation (dim M_12 = 2) and on the invariant subspace.
Prints a table; ALL CHECKS PASSED means the representation-theoretic proof goes through for every h tested."""
import sys, math, cmath
import numpy as np

HS = [int(x) for x in sys.argv[1:]] or [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 18, 22, 25, 30, 46]
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def alpha(eigs):
    s = 0.0
    for lam in eigs:
        b = (cmath.phase(lam) / (2 * math.pi)) % 1.0
        if abs(b - 1.0) < 1e-9 or abs(b) < 1e-9:
            b = 0.0
        s += b
    return s


def dims(Smat, Tmat, k):
    d = Smat.shape[0]
    ST = Smat @ Tmat
    eS = np.linalg.eigvals(cmath.exp(1j * math.pi * k / 2) * Smat)
    eST = np.linalg.eigvals(np.linalg.inv(cmath.exp(1j * math.pi * k / 3) * ST))
    eT = np.linalg.eigvals(Tmat)
    dimM = d + d * k / 12 - alpha(eS) - alpha(eST) - alpha(eT)
    kerT = sum(1 for lam in eT if abs(lam - 1) < 1e-9)
    return dimM, dimM - kerT, kerT


# validation: trivial representation
dm, ds, _ = dims(np.eye(1, dtype=complex), np.eye(1, dtype=complex), 12)
check(abs(dm - 2) < 1e-9 and abs(ds - 1) < 1e-9, f'trivial rep: dim M_12 = {dm:.3f}, dim S_12 = {ds:.3f}')

for h in HS:
    els = [(u, v) for u in range(h) for v in range(h)]
    idx = {g: i for i, g in enumerate(els)}
    N = len(els)
    q = np.array([((u * v) % h) / h for (u, v) in els])
    Tfull = np.diag(np.exp(2j * np.pi * q))
    Sfull = np.zeros((N, N), dtype=complex)
    for i, (u, v) in enumerate(els):
        for j, (u2, v2) in enumerate(els):
            Sfull[j, i] = cmath.exp(-2j * math.pi * ((u * v2 + u2 * v) % h) / h) / h
    units = [a for a in range(1, h) if math.gcd(a, h) == 1]
    # torus orbits
    orb_id = [-1] * N
    orbits = []
    for i, (u, v) in enumerate(els):
        if orb_id[i] >= 0:
            continue
        O = set()
        for a in units:
            ainv = pow(a, -1, h)
            O.add(((a * u) % h, (ainv * v) % h))
        oid = len(orbits)
        orbits.append(sorted(O))
        for g in O:
            orb_id[idx[g]] = oid
    # basis of V^- : 1_O - 1_{sigma O} over pairs O != sigma O
    sig_orb = []
    for O in orbits:
        u, v = O[0]
        sig_orb.append(orb_id[idx[(v, u)]])
    basis_vecs = []
    seen = set()
    for oid, O in enumerate(orbits):
        so = sig_orb[oid]
        if so == oid or oid in seen:
            continue
        seen.add(oid); seen.add(so)
        vec = np.zeros(N)
        for g in O:
            vec[idx[g]] += 1
        for g in orbits[so]:
            vec[idx[g]] -= 1
        basis_vecs.append(vec / np.linalg.norm(vec))
    Bm = np.array(basis_vecs).T          # N x m, orthonormal columns (orbits disjoint)
    m = Bm.shape[1]
    Sres = Bm.T @ Sfull @ Bm
    Tres = Bm.T @ Tfull @ Bm
    inv_err = np.linalg.norm(Sfull @ Bm - Bm @ Sres) + np.linalg.norm(Tfull @ Bm - Bm @ Tres)
    dm, ds, kerT = dims(Sres, Tres, 12)
    divisors = [d for d in range(1, h + 1) if h % d == 0]
    pairs = sum(1 for d in divisors if d < h // d)
    # invariant vectors inside V^-: 1_{S_d} - 1_{S_{h/d}}
    invs = []
    for d in divisors:
        if d < h // d:
            vec = np.zeros(N)
            for (u, v) in els:
                if u % d == 0 and v % (h // d) == 0:
                    vec[idx[(u, v)]] += 1
                if u % (h // d) == 0 and v % d == 0:
                    vec[idx[(u, v)]] -= 1
            invs.append(vec)
    inv_ok = all(np.linalg.norm(Sfull @ w - w) < 1e-8 and np.linalg.norm(Tfull @ w - w) < 1e-8 for w in invs)
    print(f'     h = {h:2d}: |D| = {N}, torus orbits = {len(orbits)}, dim V^- = {m}, invariance error {inv_err:.1e}; '
          f'dim M_12(V^-) = {dm:.3f}, dim S_12(V^-) = {ds:.3f}, #antisymmetric invariants = {pairs} (verified invariant: {inv_ok})', flush=True)
    check(abs(dm - round(dm)) < 1e-6 and inv_ok and inv_err < 1e-8, f'h = {h}: integral dimension, invariants verified')
    check(abs(ds - pairs) < 1e-6, f'h = {h}: dim S_12(V^-) = #antisymmetric invariants = {pairs} (so S_12(V^-) = Delta x invariants)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
