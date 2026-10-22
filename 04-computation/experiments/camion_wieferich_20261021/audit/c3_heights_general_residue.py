# c3_heights_general_residue.py -- second-pass audit of THM-4750 5'(c) for residues a != 2 (mod p): with
# kappa(a) = omega(a) sigma_p(a)^(-2) computed from DIRECT values (c1_rigidity_direct.c), the integer bases
# a_h = kappa(a) + 3 p^h must have Camion tower height exactly h (initial segment + freezing for a general base),
# and kappa(a) mod p^K must pass every level k <= K-1.  Also prints v_p(a - kappa(a)) = v_p(Lambda_p(a)) for the residue itself.
# Usage: python3 c3_heights_general_residue.py /path/to/c1_binary      Second-pass auditor, 2026-10-21.
import sys, subprocess, re
C1 = sys.argv[1]
def run(p, kmax, *bases):
    out = subprocess.run([C1, str(p), str(kmax), "explicit"] + [str(b) for b in bases], capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        m = re.match(r'p=\d+ eps=([+-]1) base a=(\d+)', line)
        if m: res[int(m.group(2))] = {int(k): (int(sv), int(dv), int(v)) for k, sv, dv, v in re.findall(r'k=(\d+): eps\^k s\(p\^k\)=(\d+), eps\^k D\(p\^k\)=(\d+) \[v=(\d+)\]', line)}
    return res
def teich(a, p, K):
    q = p ** K; x = a % q
    for _ in range(K + 3): x = pow(x, p, q)
    return x
def vp(x, p, cap):
    if x == 0: return cap
    v = 0
    while x % p == 0 and v < cap: x //= p; v += 1
    return v
allok = True
for p, kmax, residues in ((5, 5, (3, 4)), (7, 4, (3, 5, 6)), (11, 3, (3, 5, 10)), (13, 3, (3, 5, 12))):
    K = kmax + 1; q = p ** K
    for a0 in residues:
        w = teich(a0, p, K); sig = run(p, kmax, w)[w][kmax][0]
        kap = w * pow(sig * sig, -1, q) % q
        lam0 = run(p, kmax, a0)[a0][kmax][1]
        bases = [kap] + [(kap + 3 * p ** h) % q for h in range(1, kmax + 1)]
        res = run(p, kmax, *bases); line = []
        for idx, b in enumerate(bases):
            vals = [res[b][k][2] for k in range(1, kmax + 1)]
            pred = [k + 1 for k in range(1, kmax + 1)] if idx == 0 else [(k + 1 if k < idx else idx) for k in range(1, kmax + 1)]
            if vals != pred: allok = False
            line.append(f"{'kappa' if idx == 0 else 'h=' + str(idx)}:{''.join(map(str, vals))}{'' if vals == pred else '(!)'}")
        print(f"p={p:2d} a={a0:2d}: kappa(a) mod p^{K} = {kap:7d}; v(a - kappa) = {vp((a0 - kap) % q, p, K)}, v(Lambda_p(a)) = {vp(lam0, p, K)};  "
              f"valuation patterns v(D_b(p^k)), k=1..{kmax}: " + " ".join(line))
print("all designed-height patterns as predicted:", allok)
