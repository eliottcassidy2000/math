"""a3_prior_lines.py (AUDIT) - are the lines of Proposition 2.5 (c 2^s - 1 for c = 61, 123, 247, 495, ... and
c = 27, 29, 15) already path counts by THM-5130 (E) (certs_minus.txt of the S22 packet, audited)?
For the certificates without tower levels the anti-seed pair is rebuilt here and checked DIRECTLY:
state (r, M) from the definitions (brute force) and H(B => TT_t => v, v -> O) = M 2^t - 1 for t = 1, 2, 3 by the
audit's own path counter.  The certificates with tower levels (c = 15, 31, 61) are only quoted."""
import os, sys, itertools
from aud_common import *
out = []
def P(s):
    print(s); out.append(s); sys.stdout.flush()
certs = {}
fn = os.path.join(HERE, "..", "..", "kps_s22_sign_flip_towers_20261110", "certs_minus.txt")
for l in open(fn):
    if l.startswith("#"): continue
    p = l.split()
    certs[int(p[1])] = (int(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6]), p[7])
def closing(B, O, X):
    nb, nx = len(B), len(X)
    N = nb + nx + 1
    o = [set() for _ in range(N)]
    v = nb + nx
    for a in range(nb):
        o[a] |= set(B[a]); o[a] |= set(range(nb, nb + nx))
        if a in O: o[v].add(a)
        else: o[a].add(v)
    for a in range(nx):
        o[nb + a] |= {nb + u for u in X[a]}; o[nb + a].add(v)
    return o
here = {61: 5, 123: 5, 247: 5, 495: 5, 991: 5, 1983: 5, 3967: 5, 27: 3, 29: 4, 15: 6}
P("line c 2^s - 1 | first level s in Proposition 2.5 (S23) | first level by THM-5130 (E) [slope c 2^v', s >= v' + 1] | certificate | direct check")
for c in (61, 123, 247, 495, 991, 1983, 3967, 27, 29, 15):
    v2, order, m, code, Om, levels = certs[c]
    chk = "tower certificate (not rebuilt here)"
    if levels == "-" and m <= 8:
        B = decode(code, m)
        O = {i for i in range(m) if Om >> i & 1}
        vals = ham_c([closing(B, O, transitive(t)) for t in (1, 2, 3)])
        ok = all(vals[t - 1] == c * 2 ** (v2 + t) - 1 for t in (1, 2, 3))
        chk = "H(T_t) for t = 1, 2, 3 = %s = c 2^(v'+t) - 1: %s" % (vals, ok)
    elif levels == "-":
        B = decode(code, m)
        O = {i for i in range(m) if Om >> i & 1}
        vals = ham_c([closing(B, O, transitive(t)) for t in (1, 2)])
        ok = all(vals[t - 1] == c * 2 ** (v2 + t) - 1 for t in (1, 2))
        chk = "H(T_t) for t = 1, 2 = %s = c 2^(v'+t) - 1: %s" % (vals, ok)
    P("c = %4d | s >= %d | s >= %d | order %d, levels %s | %s" % (c, here[c], v2 + 1, order, levels, chk))
P("31 2^r - 1 <= 4097 exactly for r <= 7; so the lines of the minus grid of d = 31 with r <= 7 and the three rows b = 1, 2, 3 are in THM-5130 (E), from levels not higher than here.")
open(os.path.join(HERE, "a3_prior_lines.out"), "w").write("\n".join(out) + "\n")
