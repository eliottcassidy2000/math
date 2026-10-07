# NUMERICAL test of the CONDITIONAL corollary of #090:
# for a lattice Lam of covolume N in R^2 and N points on R^2/Lam, the periodic Gaussian energy
#   E(X) = (1/N) sum_{i,j} sum'_{lam} exp(-pi*al*|x_i-x_j+lam|^2)  >=  E(A) = theta_A(al) - 1,
# with equality for X = A/Lam when Lam is a sublattice of the covolume-1 triangular lattice A.
import numpy as np, sys
from scipy.optimize import minimize
rng = np.random.default_rng(20261006)
c = np.sqrt(2/np.sqrt(3)); wq = np.exp(1j*np.pi/3)
def EA(al, K=12):
    s = 0.0
    for j in range(-K, K+1):
        for l in range(-K, K+1):
            if j == 0 and l == 0: continue
            z = c*(j + l*wq); s += np.exp(-np.pi*al*abs(z)**2)
    return s
def make_torus(B):  # B: 2x2 columns = basis of Lam
    return np.array(B, dtype=float)
def energy_grad(flat, N, B, al, Kbox):
    X = flat.reshape(N, 2)
    Binv = np.linalg.inv(B)
    D = X[:, None, :] - X[None, :, :]                       # N,N,2
    co = D @ Binv.T; co -= np.round(co); D = co @ B.T       # reduce mod Lam
    ks = np.arange(-Kbox, Kbox+1)
    LL = np.array([[m, n] for m in ks for n in ks]) @ B.T   # M,2
    T = D[:, :, None, :] + LL[None, None, :, :]              # N,N,M,2
    r2 = (T**2).sum(-1)
    W = np.exp(-np.pi*al*r2)
    i0 = np.where((LL**2).sum(1) < 1e-18)[0][0]
    idx = np.arange(N); W[idx, idx, i0] = 0.0
    E = W.sum()/N
    G = (-2*np.pi*al*W[..., None]*T).sum(axis=2)           # dE/dD_ij (per ordered pair)
    grad = (G.sum(1) - G.sum(0))/N
    return E, grad.ravel()
def run(name, N, B, al, starts=60, Kbox=None, lattice_pts=None):
    B = make_torus(B)
    if Kbox is None:
        lmin = min(np.linalg.norm(B[:, 0]), np.linalg.norm(B[:, 1]))
        Kbox = int(np.ceil(4.5/ max(lmin, 0.3))) + 2
    ref = EA(al)
    best = np.inf; below = 0; vals = []
    for s in range(starts):
        x0 = (rng.random((N, 2)) @ B.T).ravel()
        res = minimize(energy_grad, x0, args=(N, B, al, Kbox), jac=True, method='L-BFGS-B',
                       options={'maxiter': 3000, 'gtol': 1e-11, 'ftol': 1e-15})
        vals.append(res.fun); best = min(best, res.fun)
        if res.fun < ref - 1e-9: below += 1
    out = f"{name:34s} N={N:2d} al={al:4.2f}  E(A)={ref:.12f}  best={best:.12f}  best-E(A)={best-ref:+.2e}  below={below}/{starts}"
    if lattice_pts is not None:
        El, _ = energy_grad(lattice_pts.ravel(), N, B, al, Kbox)
        out += f"  E(A/Lam)-E(A)={El-ref:+.1e}"
    print(out); sys.stdout.flush()
def cvec(z): return np.array([z.real, z.imag])
def sub_tri(a, b):  # Lam = alpha*A, alpha = a + b*wq
    al_ = a + b*wq
    B = np.column_stack([cvec(al_*c), cvec(al_*c*wq)])
    N = a*a + a*b + b*b
    # coset reps of A/alpha A
    pts = []; seen = set()
    Binv = np.linalg.inv(B)
    for j in range(-N, N+1):
        for l in range(-N, N+1):
            p = cvec(c*(j + l*wq)); co = Binv @ p; co -= np.floor(co + 1e-9)
            key = tuple(np.round(co, 6))
            if key not in seen: seen.add(key); pts.append(B @ co)
    assert len(pts) == N, (len(pts), N)
    return N, B, np.array(pts)
for al in (1.0, 2.0):
    N, B, P = sub_tri(2, 1);  run("hex torus alpha=2+w (N=7)", N, B, al, lattice_pts=P)
    N, B, P = sub_tri(3, 1);  run("hex torus alpha=3+w (N=13)", N, B, al, starts=40, lattice_pts=P)
    N, B, P = sub_tri(4, 1);  run("hex torus alpha=4+w (N=21)", N, B, al, starts=30, lattice_pts=P)
    # thin index-21 sublattice Z u + Z 21 v
    u = cvec(c+0j); v = cvec(c*wq); B = np.column_stack([u, 21*v]); P = np.array([k*v for k in range(21)])
    run("thin torus Zu+Z21v (N=21)", 21, B, al, starts=30, lattice_pts=P)
    # incompatible tori: square of area 7 and 21
    run("square torus area 7", 7, np.diag([np.sqrt(7)]*2), al)
    run("square torus area 21", 21, np.diag([np.sqrt(21)]*2), al, starts=30)
