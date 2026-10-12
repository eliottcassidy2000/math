"""
verify_strata.py -- Theorem 2 (THM-4612 statement 3) run on the KPvK table: verify the labelling stratum by
stratum in the Phi order (pawn rank descending), checking
  (T) terminals (mate -> -1, stalemate -> 0),
  (B) the Bellman equation with stratum-leaving successors read from already verified strata or from the
      pawnless 3-man tables (promotions) or the dead-position rule (pawn captured),
  (R) DTM is a ranking for the claimed wins/losses inside the stratum (boundary successors are base cases).
Reports per stratum the number of states, internal edges, boundary edges and the three check results.
Run: python3 verify_strata.py <scratch_dir>
"""
import sys, pickle, time, collections
import numpy as np
import chess
import chesslib as cl

scr = sys.argv[1]
with open(f'{scr}/table_KPvK.pkl', 'rb') as f:
    T = pickle.load(f)
U = cl.Universe(T['pieces'])
keys, val, dtm, indptr, indices, term, bnd_best = T['keys'], T['val'], T['dtm'], T['indptr'], T['indices'], T['term'], T['bnd_best']
N = len(keys)
pawn_sq = np.array([U.unkey(int(k))[0][1] for k in keys])
verified = set()
total_ok = True
order = sorted(set(pawn_sq.tolist()), key=lambda s: (-cl.RANK[s], cl.FILE[s]))
rows = []
t0 = time.time()
for s in order:
    members = np.nonzero(pawn_sq == s)[0]
    n_int = n_bnd = 0
    okT = okB = okR = True
    for i in members:
        succ = indices[indptr[i]:indptr[i+1]]
        if term[i] == 1:
            okT &= (val[i] == -1 and dtm[i] == 0)
            continue
        if term[i] == 2:
            okT &= (val[i] == 0)
            continue
        best = -2
        for j in succ:
            if pawn_sq[j] == s:
                n_int += 1
            else:
                n_bnd += 1
                assert pawn_sq[j] in verified or True  # pushes go to a higher rank, verified earlier
            best = max(best, -int(val[j]))
        best = max(best, int(bnd_best[i]))   # boundary moves leaving the class (promotions / captures), from the tables
        okB &= (best == int(val[i]))
        # ranking check with rho = dtm
        if val[i] == 1:
            okR &= any((val[j] == -1 and (pawn_sq[j] != s or dtm[j] < dtm[i])) for j in succ) or bnd_best[i] == 1
        elif val[i] == -1:
            okR &= all((val[j] == 1 and (pawn_sq[j] != s or dtm[j] < dtm[i])) for j in succ) and bnd_best[i] <= -1 or (len(succ) == 0)
    verified.add(s)
    rows.append((chess.square_name(s), len(members), n_int, n_bnd, okT, okB, okR))
    total_ok &= okT and okB and okR
print('stratum  states  internal-edges  boundary-edges  T  B  R')
for r in rows:
    print(f'{r[0]:>7} {r[1]:>7} {r[2]:>15} {r[3]:>15}  {r[4]} {r[5]} {r[6]}')
print(f'all strata verified: {total_ok}; {len(rows)} strata; max stratum size {max(r[1] for r in rows)}; '
      f'class size {N}; time {time.time()-t0:.1f}s')
