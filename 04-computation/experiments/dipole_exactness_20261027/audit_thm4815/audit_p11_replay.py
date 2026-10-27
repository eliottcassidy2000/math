"""audit_p11_replay.py -- independent replay of Session74's 15-move P_11 certificate (supplement/certificate15.json):
every path must be a Hamiltonian path of the current tournament; L must drop by one per move to 0; the end must be
transitive; reports F (forward arcs) along the way.  Only the 'path' fields are used (rows are compared, not trusted).
Run: python3 -I audit_p11_replay.py"""
import json
CERT = '/Users/e/Documents/GitHub/math-wt-pyth-cuboids/05-knowledge/results/session74_eleven_mechanisms_20261027/session74-paley/supplement/certificate15.json'
d = json.load(open(CERT))
n = 11; QR = {(x * x) % n for x in range(1, n)}
A = [[1 if i != j and (j - i) % n in QR else 0 for j in range(n)] for i in range(n)]
def L(A):
    t = sorted(sum(r) for r in A); return sum(abs(t[i] - i) for i in range(n)) // 2
def F(A): return sum(A[i][j] for i in range(n) for j in range(i + 1, n))
ok = True; Ls = [L(A)]; Fs = [F(A)]
for st in d['steps']:
    Q = st['path']
    if sorted(Q) != list(range(n)) or not all(A[Q[i]][Q[i + 1]] for i in range(n - 1)): ok = False; print('illegal path at depth', st['depth']); break
    B = [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)]
    for i in range(n - 1): B[Q[i]][Q[i + 1]] = 1; B[Q[i + 1]][Q[i]] = 0
    A = B
    rows = [sum(A[u][v] << v for v in range(n)) for u in range(n)]
    if 'rows' in st and rows != st['rows']: ok = False; print('row mismatch at depth', st['depth'])
    Ls.append(L(A)); Fs.append(F(A))
trans = sorted(sum(r) for r in A) == list(range(n))
print('moves', len(d['steps']), 'legal and rows match:', ok, '; L sequence', Ls, '; transitive at end:', trans, '; F sequence', Fs)
print('F(P_11) =', Fs[0], ', F(TT_id) =', n * (n - 1) // 2, ', parity of labelled distance P_11 -> TT_id:', (Fs[0] + n * (n - 1) // 2) % 2)
