# Switching families of doubly regular tournaments.  A DRT D of order n bordered by a source is the skew-Hadamard
# matrix K = I + S' of order n + 1; renormalising K at row r (switch so that r is a source) and deleting r gives the
# DRT D_r.  Prints H(D_r) and the nauty class of D_r for every r.  Usage: drt_switch.py codes_file hk_big labelg.
import sys, subprocess
import numpy as np

def code_to_S(c):
    m = len(c); n = int((1 + (1 + 8 * m) ** 0.5) / 2)
    S = np.zeros((n, n), dtype=int); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = 1 if c[k] == '1' else -1; S[i, j] = s; S[j, i] = -s; k += 1
    return S

def S_to_code(S):
    n = S.shape[0]
    return ''.join('1' if S[i, j] == 1 else '0' for i in range(n) for j in range(i + 1, n))

def d6(S):
    n = S.shape[0]; adj = (S == 1).astype(int)
    bits = [adj[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    out = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        out += chr(63 + v)
    return out

def family(S):
    n = S.shape[0]
    B = np.zeros((n + 1, n + 1), dtype=int); B[1:, 1:] = S; B[0, 1:] = 1; B[1:, 0] = -1
    m = n + 1
    assert ((np.eye(m, dtype=int) + B) @ (np.eye(m, dtype=int) + B).T == m * np.eye(m, dtype=int)).all()
    fam = []
    for r in range(m):
        d = np.array([B[r, j] if j != r else 1 for j in range(m)])
        Br = d[:, None] * B * d[None, :]
        keep = [j for j in range(m) if j != r]
        Dr = Br[np.ix_(keep, keep)]
        assert (Dr @ Dr.T == n * np.eye(n, dtype=int) - np.ones((n, n), dtype=int)).all()
        fam.append(Dr)
    return fam

if __name__ == '__main__':
    codes = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    HK, LG = sys.argv[2], sys.argv[3]
    for ci, c in enumerate(codes[:int(sys.argv[4]) if len(sys.argv) > 4 else len(codes)]):
        fam = family(code_to_S(c))
        inp = '\n'.join(S_to_code(D) for D in fam) + '\n'
        Hs = subprocess.run([HK], input=inp, capture_output=True, text=True).stdout.split()
        canon = subprocess.run([LG, '-q'], input='\n'.join(d6(D) for D in fam) + '\n', capture_output=True, text=True).stdout.split()
        classes = {}
        for H, cf in zip(Hs, canon): classes.setdefault(cf, set()).add(H)
        print(f"seed DRT {ci}: family of {len(fam)} renormalisations -> {len(classes)} nauty classes; H values by class: "
              + '; '.join(f"[{len([1 for x in canon if x == cf])}x] H = {','.join(sorted(v))}" for cf, v in classes.items()), flush=True)
