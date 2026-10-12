"""
audit_A_kbpvk_resolve.py -- audit A: INDEPENDENT exact re-solve of the KBPvK sub-class (dark bishop, a-pawn,
1,182,440 states) under the note's flagged promotion assumption (an uncapturable promotion is lost for Black),
with a retrograde solver that also resolves states whose every legal move is a losing boundary move (the case
chesslib.solve() never resolves).  Compares with kbpvk_sub.pkl and reports the corrected counts.
Usage: python3 audit_A_kbpvk_resolve.py <session_tables_dir> <my_tables.pkl>
"""
import sys, time, pickle, collections
import numpy as np
import chess
from audit_A_solver import solve, fen_of, decode_session_key, W, B, RANK, FILE, kings_adjacent

sdir, mypath = sys.argv[1], sys.argv[2]
with open(mypath, 'rb') as f:
    MY = pickle.load(f)
KP = MY['KPvK']; KP['index'] = {st: i for i, st in enumerate(KP['states'])}
DARK = [((FILE[s] + RANK[s]) % 2 == 0) for s in range(64)]
pieces = ['K', 'B', 'P', 'k']
pobj = [chess.Piece.from_symbol(c) for c in pieces]

# enumeration of the sub-class
t0 = time.time()
states = []
b = chess.Board(None)
for wk in range(64):
    for wb in [s for s in range(64) if DARK[s]]:
        for wp in [chess.square(0, r) for r in range(1, 7)]:
            for bk in range(64):
                sq = (wk, wb, wp, bk)
                if len(set(sq)) < 4 or kings_adjacent(wk, bk):
                    continue
                b.clear_board()
                for p, s in zip(pobj, sq):
                    b.set_piece_at(s, p)
                for turn in (W, B):
                    b.turn = turn
                    if b.is_attacked_by(turn, bk if turn == W else wk):
                        continue
                    states.append((sq, turn))
index = {st: i for i, st in enumerate(states)}
N = len(states)
print(f'enumerated {N} states ({time.time()-t0:.0f}s)', flush=True)


def kpvk_white_value(board):
    pm = board.piece_map(); sq = []
    for c in ['K', 'P', 'k']:
        for s, p in pm.items():
            if p.symbol() == c:
                sq.append(s); break
    return int(KP['val'][KP['index'][(tuple(sq), board.turn)]])


def escape_rule(board):
    moves = list(board.legal_moves)
    if not moves:
        return -1 if board.is_check() else 0
    for m in moves:
        if board.is_capture(m):
            board.push(m); dead = board.is_insufficient_material(); board.pop()
            if dead:
                return 0
    return None


# graph
t0 = time.time()
indptr = np.zeros(N + 1, dtype=np.int64); chunks = []
term = np.zeros(N, dtype=np.int8); bbest = np.full(N, -2, dtype=np.int8); bunc = np.zeros(N, dtype=np.int32)
n_assumed = 0
for i, (sq, turn) in enumerate(states):
    b.clear_board()
    for p, s in zip(pobj, sq):
        b.set_piece_at(s, p)
    b.turn = turn
    moves = list(b.legal_moves)
    if not moves:
        term[i] = 1 if b.is_check() else 2
        indptr[i + 1] = indptr[i]
        continue
    succ = []; best = -2
    for m in moves:
        if b.is_capture(m) or m.promotion:
            b.push(m)
            if b.is_insufficient_material():
                v = 0
            elif m.promotion:
                v = escape_rule(b)
                if v is None:
                    v = -1; n_assumed += 1        # the flagged ASSUMPTION
            else:
                v = kpvk_white_value(b)          # Kxb -> KPvK, White to move
            b.pop()
            best = max(best, -v)
        else:
            nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
            succ.append(index[(tuple(nsq), not turn)])
    bbest[i] = best
    chunks.append(np.array(succ, dtype=np.int32))
    indptr[i + 1] = indptr[i] + len(succ)
indices = np.concatenate(chunks)
print(f'graph: {len(indices)} internal edges, promotions resolved by the assumption {n_assumed} ({time.time()-t0:.0f}s)', flush=True)
G = dict(indptr=indptr, indices=indices, term=term, bbest=bbest, bunc=bunc)
val, dtm, nbol = solve(N, G)
tw = np.array([st[1] == W for st in states])
nl = np.where(tw, val != 1, val != -1)
print(f'MINE: win {int((val==1).sum())}, draw {int((val==0).sum())}, loss {int((val==-1).sum())}; max DTM {int(dtm[val!=0].max())}; '
      f'boundary-only losses resolved {nbol}; not-lost for Black {int(nl.sum())} ({100*nl.mean():.2f}%)')

# compare with the session pickle
with open(f'{sdir}/kbpvk_sub.pkl', 'rb') as f:
    S = pickle.load(f)
skeys, sval = S['keys'], S['val']
mism = []; seen = np.zeros(N, dtype=bool); missing = 0
for k, v in zip(skeys.tolist(), sval.tolist()):
    st = decode_session_key(k, 4)
    j = index.get(st)
    if j is None:
        missing += 1; continue
    seen[j] = True
    if int(val[j]) != v:
        mism.append((fen_of(pieces, *st), v, int(val[j])))
print(f'vs kbpvk_sub.pkl: session {len(skeys)} states, mine {N}; missing {missing}, unseen {int((~seen).sum())}; mismatches {len(mism)}')
cnt = collections.Counter((('BTM' if ' b ' in f else 'WTM'), v, mv) for f, v, mv in mism)
print('   mismatch types (side, session label, true label):', dict(cnt))
bk_of = collections.Counter(chess.square_name(chess.parse_square(f.split()[0] and chess.Board(f).king(B) and chess.square_name(chess.Board(f).king(B)))) for f, _, _ in mism) if False else None
bksq = collections.Counter(chess.square_name(chess.Board(f).king(B)) for f, _, _ in mism)
print('   mismatches by Black king square:', dict(sorted(bksq.items())))
for f, v, mv in mism:
    print('   ', f, 'session', v, 'mine', mv)
if len(sys.argv) > 3:
    with open(sys.argv[3], 'wb') as f:
        pickle.dump(dict(states=states, val=val, dtm=dtm, mismatches=mism), f, protocol=4)
# corrected by-square list and corner count
bk_sq = np.array([st[0][3] for st in states])
c = collections.Counter(chess.square_name(int(s)) for s in bk_sq[nl])
print('   not-lost by BK square (top 12, corrected):', c.most_common(12))
CORNER = {chess.A8, chess.B8, chess.A7, chess.B7}
in_corner = np.array([s in CORNER for s in bk_sq])
print(f'   corner states {int(in_corner.sum())}, not-lost among them {int((nl & in_corner).sum())}')
print(f'   corner share of the not-lost set: {100*in_corner[nl].mean():.2f}%')
# examples of unsafe Kxb captures from the corner region (BK in corner, Black to move)
ex = []
for i, (sq, turn) in enumerate(states):
    if turn != B or sq[3] not in CORNER or len(ex) >= 3:
        continue
    b.clear_board()
    for p, s in zip(pobj, sq):
        b.set_piece_at(s, p)
    b.turn = B
    for m in b.legal_moves:
        if b.is_capture(m) and b.piece_type_at(m.to_square) == chess.BISHOP:
            b.push(m); v = kpvk_white_value(b); b.pop()
            if v == 1:
                ex.append((b.fen(), m.uci())); break
print('   examples of corner states where Kxb lands in a lost KPvK position:', ex)
