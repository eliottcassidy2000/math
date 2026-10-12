"""
traps_analysis.py -- mechanism analysis of the two KRvKB traps (safe corner a8 vs wrong corner a1).
Rebuilds the regions (same predicates as traps.py), extracts the maximal traps again, and classifies:
  * attacker-to-move members: is the rook en prise (attacked by the bishop or the king, undefended)? is White
    stalemated? how many are 'quiet' (no boundary move at all)?
  * defender-to-move members relying only on a certified boundary move (width 1 and bnd_best >= 0 with no
    alive internal successor) versus members with an internal safe move;
  * the surviving region as a function of the Black king square and of the bishop's distance to the king.
Also checks the closure of the extracted trap with verify_trap (must report no violation).
Run: python3 traps_analysis.py <scratch_dir>
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

U3 = cl.Universe(['K', 'R', 'k', 'b'])
def near(sq, corner, d=2):
    return max(abs(cl.FILE[sq] - cl.FILE[corner]), abs(cl.RANK[sq] - cl.RANK[corner])) <= d

for corner, lab in ((chess.A8, 'safe corner a8'), (chess.A1, 'wrong corner a1')):
    def region(sq, corner=corner):
        return cl.IS_DARK[sq[3]] and near(sq[2], corner)
    kk, ii = cl.enumerate_positions(U3, square_filter=region, verbose=False)
    GG = cl.build_graph(U3, kk, ii, cl.oracle_from_tables(tables), in_region=True, verbose=False)
    aa, ss = cl.extract_trap(U3, kk, GG, chess.BLACK, label=lab, verbose=True)
    bad = cl.verify_trap(U3, kk, GG, chess.BLACK, aa)
    print(f'  closure re-verified: violations {len(bad)}')
    tw = (kk & 1).astype(bool)
    indptr, indices, term, bnd_best = GG['indptr'], GG['indices'], GG['term'], GG['bnd_best']
    # attacker (White) to move members
    att = np.nonzero(aa & tw)[0]
    c = collections.Counter()
    for i in att:
        sq, t = U3.unkey(int(kk[i])); b = U3.board(sq, t)
        wk, wr, bk, bb = sq
        rook_attacked = b.is_attacked_by(chess.BLACK, wr)
        rook_defended = b.is_attacked_by(chess.WHITE, wr)
        if term[i] == 2:
            c['White stalemated'] += 1
        elif term[i] == 1:
            c['White mated'] += 1
        elif rook_attacked and not rook_defended:
            c['rook en prise, undefended'] += 1
        elif rook_attacked:
            c['rook attacked but defended'] += 1
        else:
            c['quiet (rook safe)'] += 1
    print(f'  attacker-to-move members {len(att)}:', dict(c))
    # defender-to-move members: boundary-only vs internal
    de = np.nonzero(aa & ~tw)[0]
    c2 = collections.Counter()
    for i in de:
        succ = indices[indptr[i]:indptr[i+1]]
        internal_alive = int(aa[succ].sum()) if len(succ) else 0
        if term[i] == 2:
            c2['stalemate'] += 1
        elif internal_alive == 0:
            c2['boundary move only (captures the rook / forces a dead position)'] += 1
        elif bnd_best[i] >= 0:
            c2['internal safe move and a boundary draw'] += 1
        else:
            c2['internal safe moves only'] += 1
    print(f'  defender-to-move members {len(de)}:', dict(c2))
    # quiet attacker members by Black king square
    bk = np.array([U3.unkey(int(k))[0][2] for k in kk])
    quiet = np.zeros(len(kk), dtype=bool)
    for i in att:
        sq, t = U3.unkey(int(kk[i])); b = U3.board(sq, t)
        if term[i] == 0 and not b.is_attacked_by(chess.BLACK, sq[1]):
            quiet[i] = True
    cq = collections.Counter(chess.square_name(int(s)) for s in bk[quiet])
    print('  quiet attacker-to-move members by BK square:', dict(sorted(cq.items())))
    # region census by BK square: fraction surviving
    for s in sorted(set(bk.tolist())):
        m = bk == s
        print(f'   BK {chess.square_name(int(s))}: region {int(m.sum())}, trap {int((aa & m).sum())} ({100*(aa&m).sum()/m.sum():.1f}%), '
              f'attacker-to-move in trap {int((aa & m & tw).sum())} of {int((m & tw).sum())}')
    # a few quiet examples
    ex = [U3.fen(int(kk[i])) for i in np.nonzero(quiet)[0][:5]]
    print('  quiet examples:', ex)
