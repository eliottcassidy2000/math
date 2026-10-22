# The Camion lift kappa(a) = omega(a) * sigma_p^(-2) is the unique p-adic base a' = a (mod p) for which EVERY power p^k is a
# Camion pseudoprime to base a' (s_a'(p^k) = a'^((p^k-1)/2) mod p^k), just as the Teichmueller lift omega(a) is the unique base
# making every p^k a Fermat pseudoprime.  Check: take the integer a' = kappa(2) mod p^K and test v_p(D_a'(p^k)) for k = 1, 2, 3,
# where D_a'(n) = s_a'(n) - a'^((n-1)/2); a generic lift of 2 has h = 1.  mac-mini-2026-10-21-fanopsp.
import importlib.util, io, contextlib, sys
spec = importlib.util.spec_from_file_location("cl", "camion_lift.py")
buf = io.StringIO()
with contextlib.redirect_stdout(buf): cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
for p, K in ((7, 3), (11, 4), (13, 4)):          # interpolation needs p > 2K - 1
    q = p**K
    cl.K = K
    sig, eps = cl.sigma_lim(2, p, K)
    w = cl.teich(2, p, K); kappa = w * pow(sig * sig % q, -1, q) % q
    for a, label in ((kappa, "kappa(2)"), (2, "2"), ((2 + p) % q, "2+p")):
        vals = []
        for k in (1, 2, 3):
            n = p**k
            s = cl.smod(a, n, q)
            D = (s[n] - pow(a, (n - 1) // 2, q)) % q
            vals.append(cl.vp(D, p, K))
        print(f"p={p} base {label:9s} (= {a} mod {p}^{K}): v_p(D_a(p^k)) for k=1,2,3: {vals}  (capped at {K})")
