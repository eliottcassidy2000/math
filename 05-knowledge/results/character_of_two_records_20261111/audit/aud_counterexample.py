"""AUDIT: an explicit counterexample to HYP-9667 (bS) and (bK).

U = ((source) + D23_0) switched at the vertex set X, where D23_0 is the first doubly regular tournament of order 23
of 04-computation/experiments/sierpinski_ocf_20261022/drt23_codes_all.txt (class D23_0 of drt23_classes.json,
H = 15759953146384333) and X is given below (vertex 0 is the source, vertices 1..23 are those of D23_0).

Checks, all with the audit's own programs:
  - D23_0 is doubly regular; U is a tournament; S S^T = 23 I (so S^2 = -23 I);
  - H(U) mod 2^32 by subset dynamic programming (aud_ham.c), on U and on a relabelled converse of U;
  - H(U) mod 2^64 by inclusion-exclusion over vertex sets (aud_ham_ie.c), a second method without shared code;
  - the numbers of directed cycles of every length mod 2^32 (aud_cyc.c global), and 1 + 2 c_odd = H (mod 4).
usage: python3 aud_counterexample.py /path/aud_ham /path/aud_ham_ie /path/aud_cyc REPO_ROOT [write|nowrite [XSTR TAG]]
(default X as below; the second recorded run uses XSTR = 10000000000000000000000, the single vertex 1, TAG = N24_single)
"""
import sys, os, random, subprocess
from aud_common import from_upper_code, is_doubly_regular, check_tournament

HAM, HIE, CYC, ROOT = sys.argv[1:5]
lane = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")
codes23 = [l.strip() for l in open(os.path.join(lane, "drt23_codes_all.txt")) if l.strip() and not l.startswith("#")]
T = from_upper_code(codes23[0])
n = len(T)
N = n + 1
XSTR = sys.argv[6] if len(sys.argv) > 6 else "10110111011011111111001"   # X: vertex v (1..23) is in X iff XSTR[v-1] == '1'
TAG = sys.argv[7] if len(sys.argv) > 7 else "N24"
X = [v for v in range(1, N) if XSTR[v - 1] == "1"]

U0 = [[0] * N for _ in range(N)]
for j in range(1, N):
    U0[0][j] = 1
for a in range(n):
    for b in range(n):
        U0[a + 1][b + 1] = T[a][b]
inX = [v in X for v in range(N)]
U = [[(1 - U0[a][b]) if (a != b and inX[a] != inX[b]) else U0[a][b] for b in range(N)] for a in range(N)]
check_tournament(U)
S = [[(1 if U[a][b] else -1) if a != b else 0 for b in range(N)] for a in range(N)]
sh = all(sum(S[a][x] * S[b][x] for x in range(N)) == ((N - 1) if a == b else 0) for a in range(N) for b in range(N))


def text(V):
    return str(len(V)) + "\n" + "\n".join("".join(str(x) for x in r) for r in V) + "\n"


def run(cmd, V):
    return subprocess.run(cmd, input=text(V), capture_output=True, text=True, check=True).stdout


print("D23_0 doubly regular: %s; code: %s" % (is_doubly_regular(T), codes23[0]))
print("switching set X (vertices 1..23 of D23_0, as a 0/1 string): %s  (|X| = %d)" % (XSTR, len(X)))
print("U is a tournament on N = %d vertices; S S^T = (N-1) I: %s; scores: %s" % (N, sh, sorted(sum(r) for r in U)))
h1 = int(run([HAM, "0"], U))
random.seed(24)
perm = list(range(N))
random.shuffle(perm)
W = [[U[perm[b]][perm[a]] for b in range(N)] for a in range(N)]
h2 = int(run([HAM, "0"], W))
print("H(U) mod 2^32, subset DP: %d; on a relabelled converse: %d; H mod 4 = %d, H mod 8 = %d; (bS) predicts 1 + N/2 = %d (mod 4)" % (
    h1, h2, h1 % 4, h1 % 8, (1 + N // 2) % 4))
h3 = int(run([HIE], U))
print("H(U) mod 2^64, inclusion-exclusion: %d; mod 2^32: %d (agrees with the subset DP: %s); H mod 4 = %d" % (
    h3, h3 % 2 ** 32, h3 % 2 ** 32 == h1, h3 % 4))
c = {}
for l in run([CYC, "global"], U).split("\n"):
    if l.strip():
        _, a, b = l.split()
        c[int(a)] = int(b)
vec = "".join(str(c[l] % 2) for l in range(3, N, 2))
codd = sum(c[l] for l in range(3, N, 2))
print("cycle counts mod 2^32 by length 3..24: %s" % [c[l] for l in range(3, N + 1)])
print("parities of (c_3, c_5, ..., c_23): %s; (bK) predicts 00000000000; 1 + 2 c_odd = %d (mod 4)" % (vec, (1 + 2 * codd) % 4))
if len(sys.argv) > 5 and sys.argv[5] == "write":
    with open("counterexample_%s_adjacency.txt" % TAG, "w") as f:
        f.write(text(U))
    with open("counterexample_%s_outmasks.txt" % TAG, "w") as f:
        f.write("%d\n" % N + "".join("%d\n" % sum(1 << j for j in range(N) if U[i][j]) for i in range(N)))
    print("written: counterexample_%s_adjacency.txt (rows of 0/1), counterexample_%s_outmasks.txt (N, then N decimal out-neighbour masks)" % (TAG, TAG))
