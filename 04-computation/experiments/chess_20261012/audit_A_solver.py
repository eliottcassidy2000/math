"""
audit_A_solver.py -- INDEPENDENT retrograde WDL/DTM solver for KQvK, KRvK, KBvK, KNvK, KPvK (audit A).
Does NOT import chesslib.  Own enumeration, own legality filter, own boundary oracle, own
bucket-by-DTM retrograde.  Then compares WDL (and DTM) with the session's pickled tables.

Usage: python3 audit_A_solver.py <session_tables_dir> <my_tables_out.pkl>
"""
import sys, time, pickle, itertools, collections
import numpy as np
import chess

W, B = chess.WHITE, chess.BLACK
RANK = [chess.square_rank(s) for s in range(64)]
FILE = [chess.square_file(s) for s in range(64)]


def kings_adjacent(a, b):
    return max(abs(FILE[a] - FILE[b]), abs(RANK[a] - RANK[b])) <= 1


def enumerate_class(pieces):
    """Tablebase-legal placements: distinct squares, pawns on ranks 2-7, kings not adjacent,
    side NOT to move not in check; both sides to move.  Returns list of (squares, turn)."""
    wk = pieces.index('K'); bk = pieces.index('k')
    ranges = [[s for s in range(64) if 1 <= RANK[s] <= 6] if c in 'Pp' else list(range(64)) for c in pieces]
    pobj = [chess.Piece.from_symbol(c) for c in pieces]
    states = []
    b = chess.Board(None)
    for sq in itertools.product(*ranges):
        if len(set(sq)) < len(sq):
            continue
        if kings_adjacent(sq[wk], sq[bk]):
            continue
        b.clear_board()
        for p, s in zip(pobj, sq):
            b.set_piece_at(s, p)
        for turn in (W, B):
            b.turn = turn
            opp_king = sq[bk] if turn == W else sq[wk]
            if b.is_attacked_by(turn, opp_king):
                continue
            states.append((sq, turn))
    return states


def class_name_of(board):
    w = ''.join(sorted((p.symbol() for p in board.piece_map().values() if p.color == W), key='KQRBNP'.index))
    bl = ''.join(sorted((p.symbol().upper() for p in board.piece_map().values() if p.color == B), key='KQRBNP'.index))
    return w + 'v' + bl


def boundary_value(board, tables):
    """Value for the side to move AFTER the boundary move: 0 if dead, else from my tables, else None."""
    if board.is_insufficient_material():
        return 0
    name = class_name_of(board)
    ent = tables.get(name)
    if ent is None:
        return None
    pieces, index, val = ent['pieces'], ent['index'], ent['val']
    pm = board.piece_map()
    sq = []
    used = set()
    for c in pieces:
        for s, p in pm.items():
            if p.symbol() == c and s not in used:
                used.add(s); sq.append(s); break
    j = index.get((tuple(sq), board.turn))
    return None if j is None else int(val[j])


def build(pieces, states, index, tables):
    N = len(states)
    pobj = [chess.Piece.from_symbol(c) for c in pieces]
    indptr = np.zeros(N + 1, dtype=np.int64)
    chunks = []
    term = np.zeros(N, dtype=np.int8)
    bbest = np.full(N, -2, dtype=np.int8)
    bunc = np.zeros(N, dtype=np.int32)
    nmoves = np.zeros(N, dtype=np.int32)
    b = chess.Board(None)
    for i, (sq, turn) in enumerate(states):
        b.clear_board()
        for p, s in zip(pobj, sq):
            b.set_piece_at(s, p)
        b.turn = turn
        moves = list(b.legal_moves)
        nmoves[i] = len(moves)
        if not moves:
            term[i] = 1 if b.is_check() else 2
            indptr[i + 1] = indptr[i]
            continue
        succ = []
        best = -2
        for m in moves:
            if b.is_capture(m) or m.promotion:
                b.push(m)
                v = boundary_value(b, tables)
                b.pop()
                if v is None:
                    bunc[i] += 1
                else:
                    best = max(best, -v)
            else:
                nsq = list(sq)
                nsq[sq.index(m.from_square)] = m.to_square
                succ.append(index[(tuple(nsq), not turn)])
        bbest[i] = best
        chunks.append(np.array(succ, dtype=np.int32))
        indptr[i + 1] = indptr[i] + len(succ)
    indices = np.concatenate(chunks) if chunks else np.zeros(0, dtype=np.int32)
    return dict(indptr=indptr, indices=indices, term=term, bbest=bbest, bunc=bunc, nmoves=nmoves)


def solve(N, G):
    indptr, indices, term, bbest, bunc = G['indptr'], G['indices'], G['term'], G['bbest'], G['bunc']
    assert not bunc.any(), 'uncertified boundary moves'
    val = np.zeros(N, dtype=np.int8)
    dtm = np.full(N, -1, dtype=np.int32)
    res = np.zeros(N, dtype=bool)
    src = np.repeat(np.arange(N, dtype=np.int32), np.diff(indptr))
    order = np.argsort(indices, kind='stable')
    pred = src[order]
    pcount = np.bincount(indices, minlength=N)
    pptr = np.zeros(N + 1, dtype=np.int64)
    pptr[1:] = np.cumsum(pcount)
    nun = np.diff(indptr).astype(np.int64)
    nonlosing_bnd = bbest >= 0
    buckets = collections.defaultdict(list)
    for i in np.nonzero(term == 1)[0]:
        val[i] = -1; dtm[i] = 0; res[i] = True; buckets[0].append(int(i))
    for i in np.nonzero(term == 2)[0]:
        val[i] = 0; dtm[i] = 0; res[i] = True
    for i in np.nonzero((term == 0) & (bbest == 1))[0]:
        val[i] = 1; dtm[i] = 1; res[i] = True; buckets[1].append(int(i))
    n_boundary_only_losses = 0
    for i in np.nonzero((term == 0) & (nun == 0) & (bbest == -1))[0]:
        if not res[i]:
            val[i] = -1; dtm[i] = 1; res[i] = True; buckets[1].append(int(i)); n_boundary_only_losses += 1
    d = 0
    while buckets:
        if d not in buckets:
            d += 1
            if d > 10000:
                break
            continue
        layer = buckets.pop(d)
        for j in layer:
            if val[j] == -1:
                for p in pred[pptr[j]:pptr[j + 1]]:
                    if not res[p]:
                        val[p] = 1; dtm[p] = d + 1; res[p] = True; buckets[d + 1].append(int(p))
            else:
                for p in pred[pptr[j]:pptr[j + 1]]:
                    if res[p]:
                        continue
                    nun[p] -= 1
                    if nun[p] == 0 and not nonlosing_bnd[p]:
                        val[p] = -1; dtm[p] = d + 1; res[p] = True; buckets[d + 1].append(int(p))
        d += 1
    return val, dtm, n_boundary_only_losses


def fen_of(pieces, sq, turn):
    b = chess.Board(None)
    for c, s in zip(pieces, sq):
        b.set_piece_at(s, chess.Piece.from_symbol(c))
    b.turn = turn
    return b.fen()


def decode_session_key(key, n):
    turn = W if (key & 1) else B
    key >>= 1
    sq = []
    for _ in range(n):
        sq.append(key % 64); key //= 64
    return tuple(sq), turn


if __name__ == '__main__':
    sdir = sys.argv[1]
    out = sys.argv[2]
    tables = {}
    for pieces in (['K', 'Q', 'k'], ['K', 'R', 'k'], ['K', 'B', 'k'], ['K', 'N', 'k'], ['K', 'P', 'k']):
        name = ''.join(c for c in pieces if c.isupper()) + 'v' + ''.join(c.upper() for c in pieces if c.islower())
        t0 = time.time()
        states = enumerate_class(pieces)
        index = {st: i for i, st in enumerate(states)}
        N = len(states)
        t1 = time.time()
        G = build(pieces, states, index, tables)
        t2 = time.time()
        val, dtm, nbol = solve(N, G)
        t3 = time.time()
        turns = np.array([1 if st[1] == W else 0 for st in states], dtype=bool)
        rec = dict(N=N, edges=int(len(G['indices'])),
                   wtm_win=int(((val == 1) & turns).sum()), wtm_draw=int(((val == 0) & turns).sum()), wtm_loss=int(((val == -1) & turns).sum()),
                   btm_win=int(((val == 1) & ~turns).sum()), btm_draw=int(((val == 0) & ~turns).sum()), btm_loss=int(((val == -1) & ~turns).sum()),
                   max_dtm=int(dtm[val != 0].max()) if (val != 0).any() else -1, boundary_only_losses=nbol)
        print(f'[{name}] MINE: {rec}  (enum {t1-t0:.0f}s, build {t2-t1:.0f}s, solve {t3-t2:.0f}s)', flush=True)
        if rec['max_dtm'] >= 0:
            i = int(np.argmax(np.where(val != 0, dtm, -1)))
            print(f'   longest DTM {rec["max_dtm"]} plies at {fen_of(pieces, *states[i])}')
            # longest white-to-move win
            m = np.where((val == 1) & turns, dtm, -1)
            i = int(np.argmax(m))
            print(f'   longest white-to-move win {int(m[i])} plies at {fen_of(pieces, *states[i])}')
        tables[name] = dict(pieces=pieces, index=index, val=val, dtm=dtm, states=states, term=G['term'], bbest=G['bbest'])
        # compare with the session table
        try:
            with open(f'{sdir}/table_{name}.pkl', 'rb') as f:
                T = pickle.load(f)
        except FileNotFoundError:
            print('   (no session table)'); continue
        assert list(T['pieces']) == pieces, (T['pieces'], pieces)
        skeys, sval, sdtm = T['keys'], T['val'], T['dtm']
        missing = 0; wdl_bad = 0; dtm_bad = 0
        seen = np.zeros(N, dtype=bool)
        for k, v, dd in zip(skeys.tolist(), sval.tolist(), sdtm.tolist()):
            st = decode_session_key(k, len(pieces))
            j = index.get(st)
            if j is None:
                missing += 1; continue
            seen[j] = True
            if int(val[j]) != v:
                wdl_bad += 1
                if wdl_bad <= 5:
                    print('   WDL MISMATCH', fen_of(pieces, *st), 'mine', int(val[j]), 'session', v)
            elif v != 0 and int(dtm[j]) != dd:
                dtm_bad += 1
                if dtm_bad <= 5:
                    print('   DTM MISMATCH', fen_of(pieces, *st), 'mine', int(dtm[j]), 'session', dd)
        print(f'   vs session table: session has {len(skeys)} states, mine {N}; session states not in mine {missing}; '
              f'mine not in session {int((~seen).sum())}; WDL mismatches {wdl_bad}; DTM mismatches (nonzero values) {dtm_bad}', flush=True)
    with open(out, 'wb') as f:
        pickle.dump({k: dict(pieces=v['pieces'], states=v['states'], val=v['val'], dtm=v['dtm'], term=v['term'], bbest=v['bbest']) for k, v in tables.items()}, f, protocol=4)
    print('saved', out)
