"""
solve_small.py -- exact WDL/DTM tables for KQvK, KRvK, KBvK, KNvK, KPvK by retrograde analysis
(chesslib.solve), with the classical checks, and the per-stratum decomposition of KPvK:
the 48 pawn squares are the strata; inside a stratum only the two kings move.

Writes tables to the scratch directory given as argv[1] (pickle) and prints a census.
Run: python3 solve_small.py <scratch_dir>
"""
import sys, pickle, time, collections
import numpy as np
import chess
import chesslib as cl

out_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
tables = {}
summary = {}
for pieces in (['K', 'Q', 'k'], ['K', 'R', 'k'], ['K', 'B', 'k'], ['K', 'N', 'k'], ['K', 'P', 'k']):
    U, keys, index, G, val, dtm = cl.solve_class(pieces, tables)
    wtm = (keys & 1).astype(bool)
    rec = dict(N=len(keys),
               wtm_win=int(((val == 1) & wtm).sum()), wtm_draw=int(((val == 0) & wtm).sum()), wtm_loss=int(((val == -1) & wtm).sum()),
               btm_win=int(((val == 1) & ~wtm).sum()), btm_draw=int(((val == 0) & ~wtm).sum()), btm_loss=int(((val == -1) & ~wtm).sum()),
               max_dtm_plies=int(dtm[val != 0].max()) if (val != 0).any() else -1,
               edges=int(len(G['indices'])))
    summary[U.name] = rec
    print(' ', rec, flush=True)
    with open(f'{out_dir}/table_{U.name}.pkl', 'wb') as f:
        pickle.dump(dict(pieces=pieces, keys=keys, val=val, dtm=dtm, indptr=G['indptr'], indices=G['indices'],
                         term=G['term'], bnd_best=G['bnd_best'], nmoves=G['nmoves']), f, protocol=4)

# ---- KPvK by strata: the pawn square is the stratum coordinate ----
U, index, val = tables['KPvK']
with open(f'{out_dir}/table_KPvK.pkl', 'rb') as f:
    T = pickle.load(f)
keys, dtm = T['keys'], T['dtm']
pawn_sq = np.array([U.unkey(int(k))[0][1] for k in keys])
wtm = (keys & 1).astype(bool)
print('\nKPvK by stratum (pawn square): positions, white-to-move wins, black-to-move draws')
strata = collections.OrderedDict()
for s in sorted(set(pawn_sq.tolist()), key=lambda s: (cl.RANK[s], cl.FILE[s])):
    m = pawn_sq == s
    strata[chess.square_name(s)] = (int(m.sum()), int(((val == 1) & wtm & m).sum()), int(wtm[m].sum()),
                                    int(((val == 0) & ~wtm & m).sum()), int((~wtm)[m].sum()))
for rank in range(1, 7):
    row = []
    for file in range(8):
        s = chess.square(file, rank)
        n, w, nw, d, nb = strata[chess.square_name(s)]
        row.append(f'{chess.square_name(s)}:{w}/{nw}|{d}/{nb}')
    print('  ', '  '.join(row))
# stratum structure: inside a stratum all moves are king moves (reversible); every pawn move leaves it
# count internal edges that stay in the stratum vs leave
indptr, indices = T['indptr'], T['indices']
src = np.repeat(np.arange(len(keys)), np.diff(indptr))
same = (pawn_sq[src] == pawn_sq[indices])
print(f'KPvK internal edges: {len(indices)}; king moves (stay in stratum): {int(same.sum())}; pawn pushes (leave stratum): {int((~same).sum())}')
# pawn pushes always go to a higher rank: check the potential
assert (cl.RANK[0] == 0)
push_ok = all(cl.RANK[int(pawn_sq[indices[e]])] > cl.RANK[int(pawn_sq[src[e]])] for e in np.nonzero(~same)[0])
print('every stratum-leaving internal move advances the pawn:', push_ok)
# longest wins
for name in ('KQvK', 'KRvK', 'KPvK'):
    with open(f'{out_dir}/table_{name}.pkl', 'rb') as f:
        T = pickle.load(f)
    U2 = cl.Universe(T['pieces'])
    i = int(np.argmax(np.where(T['val'] != 0, T['dtm'], -1)))
    print(f'{name}: longest DTM {int(T["dtm"][i])} plies at {U2.fen(int(T["keys"][i]))} (KPvK plies stop at the promotion boundary)')
with open(f'{out_dir}/summary_small.pkl', 'wb') as f:
    pickle.dump(summary, f)
print('\nSUMMARY', summary)
