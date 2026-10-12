"""
traps_krvkb_exceptions.py -- what the safe-corner KRvKB trap excludes.
Loads the pickles written by traps.py (keys, alive, safe) for the a8 region and tabulates simple features of
the excluded states versus the trap members, then measures how much of the complement a few human predicates
cover ("the bishop hangs", "Black is in check and the bishop hangs", "the king is in check on the edge", ...).
Run: python3 traps_krvkb_exceptions.py <scratch_dir>
"""
import sys, pickle, collections
import numpy as np
import chess
import chesslib as cl

scr = sys.argv[1]
U = cl.Universe(['K', 'R', 'k', 'b'])
with open(f'{scr}/krvkb_{chess.A8}.pkl', 'rb') as f:
    S = pickle.load(f)
keys, alive = S['keys'], S['alive']
N = len(keys)
tw = (keys & 1).astype(bool)
feat = collections.defaultdict(lambda: [0, 0])   # feature -> [in trap, excluded]
preds = collections.OrderedDict()
preds['BTM: bishop attacked and undefended'] = [0, 0]
preds['BTM: in check'] = [0, 0]
preds['BTM: in check and bishop attacked'] = [0, 0]
preds['WTM: rook attacks bishop and bishop undefended'] = [0, 0]
preds['WTM: White can give check next move'] = [0, 0]
preds['any: bishop attacked by rook or king and undefended'] = [0, 0]
b = chess.Board(None)
for i in range(N):
    sq, t = U.unkey(int(keys[i]))
    wk, wr, bk, bb = sq
    b.clear_board()
    for p, s in zip(U.pobj, sq):
        b.set_piece_at(s, p)
    b.turn = t
    col = 0 if alive[i] else 1
    bishop_att = b.is_attacked_by(chess.WHITE, bb)
    bishop_def = b.is_attacked_by(chess.BLACK, bb)
    in_check = b.is_check()
    side = 'BTM' if t == chess.BLACK else 'WTM'
    feat[(side, 'bishop attacked', bishop_att)][col] += 1
    feat[(side, 'bishop undefended', not bishop_def)][col] += 1
    feat[(side, 'in check', in_check)][col] += 1
    feat[(side, 'BK square', chess.square_name(bk))][col] += 1
    feat[(side, 'WK-BK distance', chess.square_distance(wk, bk))][col] += 1
    feat[(side, 'bishop-BK distance', chess.square_distance(bb, bk))][col] += 1
    if t == chess.BLACK:
        if bishop_att and not bishop_def: preds['BTM: bishop attacked and undefended'][col] += 1
        if in_check: preds['BTM: in check'][col] += 1
        if in_check and bishop_att: preds['BTM: in check and bishop attacked'][col] += 1
    else:
        rook_att = b.is_attacked_by(chess.WHITE, bb)
        if rook_att and not bishop_def: preds['WTM: rook attacks bishop and bishop undefended'][col] += 1
        can_check = False
        for m in b.legal_moves:
            if b.gives_check(m):
                can_check = True; break
        if can_check: preds['WTM: White can give check next move'][col] += 1
    if bishop_att and not bishop_def: preds['any: bishop attacked by rook or king and undefended'][col] += 1
n_in, n_out = int(alive.sum()), int((~alive).sum())
print(f'region {N}: trap {n_in}, excluded {n_out}')
print('\nfeatures (in-trap count, excluded count):')
for k in sorted(feat, key=lambda k: (k[0], k[1], str(k[2]))):
    v = feat[k]
    print(f'  {k}: in {v[0]}, out {v[1]}, excluded share {100*v[1]/(v[0]+v[1]):.1f}%')
print('\nhuman predicates (members of the predicate inside the trap / in the complement; coverage of the complement):')
for k, v in preds.items():
    print(f'  {k}: in-trap {v[0]}, excluded {v[1]}, covers {100*v[1]/n_out:.1f}% of the complement, false positives (in-trap) {v[0]}')
