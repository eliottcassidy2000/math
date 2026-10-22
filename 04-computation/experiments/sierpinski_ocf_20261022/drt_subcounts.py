# Subtournament census of P_19 and the non-Paley doubly regular D_19: 4-vertex types (predicted equal for all DRTs
# of one order) and 5-vertex types (where differences may start).  Type = sorted score sequence + 3-cycle count,
# refined by H of the subtournament.  mac-mini-2026-10-22-sierp.
import sys, itertools
from collections import Counter
sys.path.insert(0, '.')
from drt_switch import code_to_S
def subtype(S, X):
    sc = tuple(sorted(sum(1 for w in X if w != v and S[v][w] == 1) for v in X))
    H = 0
    for p in itertools.permutations(X):
        if all(S[p[i]][p[i + 1]] == 1 for i in range(len(p) - 1)): H += 1
    return (sc, H)
q = 19; sq = {(x * x) % q for x in range(1, q)}
P = [[0 if i == j else (1 if (j - i) % q in sq else -1) for j in range(q)] for i in range(q)]
D = code_to_S(open(sys.argv[1]).readline().strip()).tolist()
for k in (4, 5):
    cP = Counter(subtype(P, X) for X in itertools.combinations(range(q), k))
    cD = Counter(subtype(D, X) for X in itertools.combinations(range(q), k))
    print(f"k={k}: P_19 {dict(sorted(cP.items()))}")
    print(f"     D_19 {dict(sorted(cD.items()))}   equal: {cP == cD}", flush=True)
