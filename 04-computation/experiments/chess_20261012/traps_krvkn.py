"""
traps_krvkn.py -- a central fortress in KRvKN (rook versus knight): the largest trap inside the region
{Black king on d4, e4, d5 or e5; knight anywhere; White king and rook anywhere; both sides to move}.
Boundary: RxN or KxN -> KRvK table (Black to move); NxR or KxR -> dead position.  Everything certified.
Control: the same with the Black king in a corner square (a1) -- the knight's side should keep much less.
Run: python3 traps_krvkn.py <scratch_dir>
"""
import sys, pickle, collections
import numpy as np
import chess
import chesslib as cl

scr = sys.argv[1]
tables = {}
for name in ('KQvK', 'KRvK', 'KBvK', 'KNvK', 'KPvK'):
    with open(f'{scr}/table_{name}.pkl', 'rb') as f:
        T = pickle.load(f)
    U = cl.Universe(T['pieces'])
    tables[name] = (U, {int(k): i for i, k in enumerate(T['keys'])}, T['val'])

U = cl.Universe(['K', 'R', 'k', 'n'])
CENTRE16 = {chess.square(f, r) for f in range(2, 6) for r in range(2, 6)}
def near_knight(sq):
    return sq[2] in CENTRE16 and chess.square_distance(sq[2], sq[3]) <= 2
def box4(sq):
    return sq[2] in {chess.D4, chess.E4, chess.D5, chess.E5}
def corner(sq):
    return sq[2] == chess.A1
REGIONS = [('centre d4 e4 d5 e5', box4), ('corner a1 (control)', corner), ('centre c3-f6 with knight within 2 of the king', near_knight)]
only = sys.argv[2] if len(sys.argv) > 2 else None
for lab, region in REGIONS:
    if only and only not in lab:
        continue
    kk, ii = cl.enumerate_positions(U, square_filter=region)
    GG = cl.build_graph(U, kk, ii, cl.oracle_from_tables(tables), in_region=True)
    aa, ss = cl.extract_trap(U, kk, GG, chess.BLACK, label=lab)
    bad = cl.verify_trap(U, kk, GG, chess.BLACK, aa)
    tw = (kk & 1).astype(bool)
    print(f'  [{lab}] trap {int(aa.sum())} of {len(kk)} ({100*aa.mean():.1f}%); Black to move {int((aa & ~tw).sum())}, White to move {int((aa & tw).sum())}; closure violations {len(bad)}')
    # quiet attacker-to-move members (knight not attacked, rook not attacked)
    quiet = 0; att = np.nonzero(aa & tw)[0]
    for i in att:
        sq, t = U.unkey(int(kk[i])); b = U.board(sq, t)
        if GG['term'][i] == 0 and not b.is_attacked_by(chess.BLACK, sq[1]):
            quiet += 1
    print(f'   quiet attacker-to-move members (rook not attacked): {quiet}')
    de = np.nonzero(aa & ~tw)[0]
    c2 = collections.Counter()
    indptr, indices, term, bnd_best = GG['indptr'], GG['indices'], GG['term'], GG['bnd_best']
    for i in de:
        succ = indices[indptr[i]:indptr[i+1]]
        internal_alive = int(aa[succ].sum()) if len(succ) else 0
        if term[i] == 2:
            c2['stalemate'] += 1
        elif internal_alive == 0:
            c2['boundary move only'] += 1
        else:
            c2['internal safe move'] += 1
    print('   defender-to-move members:', dict(c2))
    w = ss[aa & ~tw & (GG['term'] == 0)]
    if len(w):
        print(f'   width: min {int(w.min())}, median {int(np.median(w))}, max {int(w.max())}')
    # knight distance to the Black king among quiet attacker-to-move members
    dist = collections.Counter()
    for i in att:
        sq, t = U.unkey(int(kk[i]))
        dist[chess.square_distance(sq[2], sq[3])] += 1
    print('   attacker-to-move members by king-knight Chebyshev distance:', dict(sorted(dist.items())))
