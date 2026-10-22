# |Aut(Sigma_N)| via nauty (dreadnaut, digraph mode) for the Sierpinski interval tournaments, against v_3(H(Sigma_N)).
# Free action of Aut on Hamiltonian paths gives |Aut| | H.  mac-mini-2026-10-22-sierp.
import subprocess, re, math
def arc(x, y):
    z = y - x
    while z % 3 == 0: z //= 3
    return z % 3 == 1
H = {}; V3 = {}
for line in open('sierpinski_interval_90.out'):
    m = re.match(r'N=\s*(\d+) H=(\d+) v3=(\d+)', line)
    if m: H[int(m.group(1))] = int(m.group(2)); V3[int(m.group(1))] = int(m.group(3))
def autsize(N):
    lines = [f"n={N} $=0 d g"]
    for x in range(N):
        outs = [str(y) for y in range(N) if y != x and arc(x, y)]
        lines.append(' '.join(outs) + (';' if x < N - 1 else '.'))
    script = '\n'.join(lines) + '\nx\nq\n'
    out = subprocess.run(['dreadnaut'], input=script, capture_output=True, text=True).stdout
    m = re.search(r'grpsize=([\d.e+]+)', out)
    if not m: return None, out[-300:]
    g = m.group(1)
    if 'e' in g or '.' in g:
        v = round(math.log(float(g)) / math.log(3))
        assert abs(math.log(float(g)) / math.log(3) - v) < 1e-6, g
        return 3 ** v, None          # nauty prints large orders in floating point; all orders here are powers of 3
    return int(g), None
res = []
for N in range(1, 91):
    a, err = autsize(N)
    if a is None: print("parse fail", N, err); break
    va = 0; t = a
    while t % 3 == 0: t //= 3; va += 1
    res.append((N, a, va, V3[N], H[N] % a == 0))
for N, a, va, vh, div in res:
    print(f"N={N:3d} |Aut|=3^{va}*{a // 3**va}  v3(H)={vh}  excess={vh - va}  Aut|H: {div}")
