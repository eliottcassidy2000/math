# Is the Camion lift kappa: F_p^x -> Z_p^x (kappa(a) = omega(a) sigma_p(a)^(-2), sigma_p(a) = lim (a/p)^k s_a(p^k)) multiplicative like
# the Teichmueller lift?  Computes kappa(a) mod p^K for every residue a, and tests kappa(a)kappa(b) = kappa(ab) and kappa(a)^(p-1) = 1.
# mac-mini-2026-10-21-fanopsp.
import importlib.util, io, contextlib
spec = importlib.util.spec_from_file_location("cl", "camion_lift.py")
buf = io.StringIO()
with contextlib.redirect_stdout(buf): cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
for p, K in ((7, 3), (11, 3), (13, 3), (17, 3)):
    q = p**K; cl.K = K
    kap = {}
    for a in range(1, p):
        sig, eps = cl.sigma_lim(a, p, K)
        kap[a] = cl.teich(a, p, K) * pow(sig * sig % q, -1, q) % q
    mult = sum(1 for a in range(1, p) for b in range(1, p) if kap[a] * kap[b] % q == kap[a * b % p] % q)
    mult2 = sum(1 for a in range(1, p) for b in range(1, p) if (kap[a] * kap[b] - kap[a * b % p]) % (p * p) == 0)
    teich_eq = [a for a in range(1, p) if (kap[a] - cl.teich(a, p, K)) % (p * p) == 0]
    print(f"p={p}: kappa(1) = {kap[1]}; multiplicative mod p^{K}: {mult}/{(p-1)**2} pairs, mod p^2: {mult2}/{(p-1)**2};"
          f" kappa(a) = omega(a) mod p^2 for a in {teich_eq}; kappa(-1) = {kap[p-1]} (omega(-1) = {q-1})")
# Teichmueller rigidity for other residues: sigma_p(a) = sigma_p(a + t p)?
print("rigidity of sigma_p for residues a = 3, 5, p-1 and their lifts a + t p (t = 1, 2, 7):")
for p, K in ((7, 3), (11, 3), (13, 3)):
    cl.K = K
    for a in (3, 5, p - 1):
        base = cl.sigma_lim(a, p, K)[0]
        print(f"  p={p} a={a}: equal for all lifts: {all(cl.sigma_lim(a + t * p, p, K)[0] == base for t in (1, 2, 7))}")
