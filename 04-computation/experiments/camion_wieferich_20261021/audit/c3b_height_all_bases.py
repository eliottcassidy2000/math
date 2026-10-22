# c3b_height_all_bases.py -- second-pass audit of THM-4750 5'(c), broad version: for EVERY integer base 1 <= a < p^K, p | a excluded,
# compare the Camion tower height read from direct values (v_p of eps^k D_a(p^k) mod p^(k+1), k = 1..kmax, via the initial-segment
# and freezing lemma) with v_p(a - kappa(a)), kappa computed per residue from the Teichmueller base.  Heights >= kmax are capped.
# Usage: python3 c3b_height_all_bases.py /path/to/c1_binary        Second-pass auditor, 2026-10-21.
import sys, subprocess, re
from collections import Counter
C1 = sys.argv[1]
def run(p, kmax, bases):
    res = {}
    for i in range(0, len(bases), 60):
        out = subprocess.run([C1, str(p), str(kmax), "explicit"] + [str(b) for b in bases[i:i + 60]], capture_output=True, text=True, check=True).stdout
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
for p, kmax, K in ((5, 4, 4), (7, 3, 3)):          # bases a < p^K; direct values mod p^(kmax+1)
    Q = p ** (kmax + 1)
    kap = {}
    for a0 in range(1, p):
        w = teich(a0, p, kmax + 1); sig = run(p, kmax, [w])[w][kmax][0]; kap[a0] = w * pow(sig * sig, -1, Q) % Q
    bases = [a for a in range(1, p ** K) if a % p]
    res = run(p, kmax, bases)
    bad = 0; hist = Counter()
    for a in bases:
        vals = [res[a][k][2] for k in range(1, kmax + 1)]          # v(D_a(p^k)) capped at k+1
        h = next((k for k in range(1, kmax + 1) if vals[k - 1] < k + 1), None)   # first k with v < k+1: then h = v there
        h_dir = vals[h - 1] if h is not None else kmax              # (capped at kmax)
        h_kap = min(vp((a - kap[a % p]) % Q, p, kmax + 1), kmax)
        if h_dir != h_kap: bad += 1
        hist[h_kap] += 1
    print(f"p={p}: {len(bases)} integer bases a < {p}^{K}: tower height from direct D_a(p^k) (k <= {kmax}) vs v_p(a - kappa(a)) "
          f"(both capped at {kmax}): {bad} mismatches; height histogram {dict(sorted(hist.items()))}")
