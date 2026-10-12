"""
audit_A_krvkb_sample.py -- audit A: re-read the two saved KRvKB traps (krvkb_56.pkl = safe corner a8,
krvkb_0.pkl = wrong corner a1), recount the headline numbers, and check closure of a random sample of trap
members with OWN move generation (no chesslib): membership = the saved 'alive' flags, boundary = dead
positions and MY independent KRvK table (audit_A_solver.py).  Also recomputes the width of each sampled
defender-to-move member and compares with the saved 'safe' value.
Usage: python3 audit_A_krvkb_sample.py <session_tables_dir> <my_tables.pkl> [n_sample]
"""
import sys, pickle, random, collections, time
import numpy as np
import chess

W, B = chess.WHITE, chess.BLACK
RANK = [chess.square_rank(s) for s in range(64)]
FILE = [chess.square_file(s) for s in range(64)]
sdir, mypath = sys.argv[1], sys.argv[2]
n_sample = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
with open(mypath, 'rb') as f:
    MY = pickle.load(f)
KR = MY['KRvK']; KR['index'] = {st: i for i, st in enumerate(KR['states'])}
pieces = ['K', 'R', 'k', 'b']


def decode(key, n=4):
    turn = W if (key & 1) else B
    key >>= 1
    sq = []
    for _ in range(n):
        sq.append(key % 64); key //= 64
    return tuple(sq), turn


def encode(sq, turn):
    k = 0
    for s in reversed(sq):
        k = k * 64 + s
    return k * 2 + (1 if turn == W else 0)


def near(sq, corner, d=2):
    return max(abs(FILE[sq] - FILE[corner]), abs(RANK[sq] - RANK[corner])) <= d


def krvk_black_value(board):
    """board: KRvK with Black to move after a capture of the bishop; returns Black's (mover) value from my table."""
    if board.is_insufficient_material():
        return 0
    pm = board.piece_map()
    sq = []
    for c in ['K', 'R', 'k']:
        for s, p in pm.items():
            if p.symbol() == c:
                sq.append(s); break
    return int(KR['val'][KR['index'][(tuple(sq), board.turn)]])


for corner, fname, label in ((chess.A8, 'krvkb_56.pkl', 'safe corner a8'), (chess.A1, 'krvkb_0.pkl', 'wrong corner a1')):
    with open(f'{sdir}/{fname}', 'rb') as f:
        P = pickle.load(f)
    keys, alive, safe = P['keys'], P['alive'], P['safe']
    tw = (keys & 1).astype(bool)
    print(f'[{label}] region {len(keys)}; alive {int(alive.sum())}; black-to-move alive {int((alive & ~tw).sum())}, '
          f'white-to-move alive {int((alive & tw).sum())}')
    wh = collections.Counter(int(x) for x in safe[alive & ~tw])
    print(f'   width histogram (saved): {dict(sorted(wh.items()))}; mode {max(wh, key=wh.get)}')
    # region sanity: every key has BK within distance 2 of the corner and a dark bishop
    bad_region = 0
    for k in keys[:200000].tolist():
        sq, t = decode(k)
        if not near(sq[2], corner) or (FILE[sq[3]] + RANK[sq[3]]) % 2 != 0:
            bad_region += 1
    print(f'   region predicate violated on first 200k keys: {bad_region}')
    adict = {int(k): bool(a) for k, a in zip(keys.tolist(), alive.tolist())}
    sdict = {int(k): int(s) for k, s in zip(keys.tolist(), safe.tolist())}
    random.seed(7)
    members = [k for k, a in adict.items() if a]
    sample = random.sample(members, n_sample)
    viol = []; width_mismatch = 0; n_def = 0; kinds = collections.Counter()
    t0 = time.time()
    for k in sample:
        sq, turn = decode(k)
        b = chess.Board(None)
        for c, s in zip(pieces, sq):
            b.set_piece_at(s, chess.Piece.from_symbol(c))
        b.turn = turn
        moves = list(b.legal_moves)
        if turn == B:
            n_def += 1
            if not moves:
                if b.is_check():
                    viol.append(('defender mated', b.fen()))
                else:
                    kinds['stalemate'] += 1
                    if sdict[k] != 0:
                        width_mismatch += 1
                continue
            internal = 0; bnd = 0
            for m in moves:
                if b.is_capture(m):
                    b.push(m); dead = b.is_insufficient_material(); b.pop()
                    if dead:
                        bnd = 1
                    else:
                        viol.append(('unexpected non-dead capture by defender', b.fen(), m.uci()))
                else:
                    nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
                    k2 = encode(nsq, W)
                    if adict.get(k2, False):
                        internal += 1
            w = internal + bnd
            kinds['internal only' if (internal and not bnd) else ('boundary only' if (bnd and not internal) else 'both')] += 1
            if w == 0:
                viol.append(('defender has no safe option', b.fen()))
            if w != sdict[k]:
                width_mismatch += 1
        else:
            if not moves:
                kinds['attacker no moves'] += 1
                continue
            for m in moves:
                if b.is_capture(m):
                    b.push(m); vb = krvk_black_value(b); b.pop()
                    if vb < 0:
                        viol.append(('attacker capture wins', b.fen(), m.uci())); break
                else:
                    nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
                    k2 = encode(nsq, B)
                    if not adict.get(k2, False):
                        viol.append(('attacker leaves the trap', b.fen(), m.uci())); break
    print(f'   sample {n_sample} members: closure violations {len(viol)}; defender members {n_def}, width mismatches vs saved {width_mismatch}; '
          f'kinds {dict(kinds)}; {time.time()-t0:.0f}s')
    for v in viol[:5]:
        print('   ', v)
