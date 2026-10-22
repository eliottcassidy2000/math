# Wormald's Conjecture 1 for regular base tournaments: among regular tournaments on d vertices the blow-up gain
# rho = prod_k (d^2 + s_k^2)/(d^2 - s_k^2)  (s_k: skew spectrum) is maximised by the rotational R_d.
# Exhaustive over all regular d-tournaments (gentourng -d m -D m d), d = 5..13.  Also reports the minimiser and
# whether the minimum is the doubly regular value ((d+1)/(d-1))^((d-1)/2).  mac-mini-2026-10-22-sierp.
import subprocess, sys
import numpy as np
def rho_from_code(code, d):
    A = np.zeros((d, d)); p = 0
    for i in range(d):
        for j in range(i + 1, d):
            if code[p] == '1': A[i, j] = 1
            else: A[j, i] = 1
            p += 1
    S = A - A.T
    ev = np.linalg.eigvals(S)
    s2 = np.sort(np.abs(ev.imag)) ** 2           # each sigma^2 appears twice (conjugate pair) plus one 0
    s2 = s2[1:][::2]                              # one per pair
    return float(np.prod((d * d + s2) / (d * d - s2))), A
for d in [int(x) for x in sys.argv[1:]] or [5, 7, 9, 11]:
    m = (d - 1) // 2
    codes = [c for c in subprocess.run(['gentourng', f'-d{m}', f'-D{m}', str(d)], capture_output=True, text=True).stdout.split() if not c.startswith('>')]
    rows = [(rho_from_code(c, d)[0], c) for c in codes]
    rows.sort()
    # rotational R_d value
    A = np.zeros((d, d))
    for i in range(d):
        for r in range(1, m + 1): A[i, (i + r) % d] = 1
    S = A - A.T; s2 = np.sort(np.abs(np.linalg.eigvals(S).imag)) ** 2; s2 = s2[1:][::2]
    rR = float(np.prod((d * d + s2) / (d * d - s2)))
    dr = ((d + 1) / (d - 1)) ** ((d - 1) / 2)
    top = rows[-1][0]; second = rows[-2][0] if len(rows) > 1 else None
    nmax = sum(1 for r, _ in rows if abs(r - top) < 1e-9)
    sec_str = ('%.10f' % second) if second is not None else '-'
    print(f"d={d}: {len(rows)} regular classes; max rho = {top:.10f} (attained by {nmax} class(es)); rho(R_d) = {rR:.10f}; R_d is the maximiser: {abs(top - rR) < 1e-9}; "
          f"second = {sec_str}; min rho = {rows[0][0]:.10f}; doubly-regular bound = {dr:.10f}", flush=True)
