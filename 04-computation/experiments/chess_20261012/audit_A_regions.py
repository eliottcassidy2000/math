"""
audit_A_regions.py -- INDEPENDENT closure checks for audit A.  Does NOT import chesslib.
 (a) KPvK rook-pawn corner predicate {WP on the a-file, BK in {a8,b8,a7,b7}}: enumerate the region with own code,
     check the closure clauses of the note (defender: stalemate / move into F / certified boundary draw;
     attacker: no move / every move into F or a certified boundary move with defender value >= 0), widths;
     h-file mirror.  Promotions certified by MY independent 3-man tables (audit_A_solver.py), and, as a second
     run, by the one-ply escape rule only.
 (b) KBPvK wrong-bishop corner predicate {WB dark, WP a-file, BK in corner}: closure with the stated rules
     (Kxp -> dead; Kxb -> my KPvK table; promotion -> one-ply capture into a dead position / stalemate; else unsafe).
 (c) probe of a latent chesslib.solve() gap on kbpvk_sub.pkl: Black-to-move states labelled 0 whose every legal
     move is a losing boundary capture (solve() never resolves such states).  Then a full (T)+(B) Bellman
     re-check of the whole sub-class pickle with own move generation and the flagged promotion assumption.
 (d) 2,000 random not-lost states of table_KPvK.pkl satisfy the closure clauses w.r.t. the not-lost set.
Usage: python3 audit_A_regions.py <session_tables_dir> <my_tables.pkl> [skip_full]
"""
import sys, pickle, random, time, collections
import numpy as np
import chess

W, B = chess.WHITE, chess.BLACK
RANK = [chess.square_rank(s) for s in range(64)]
FILE = [chess.square_file(s) for s in range(64)]
DARK = [((FILE[s] + RANK[s]) % 2 == 0) for s in range(64)]
assert DARK[chess.A1] and not DARK[chess.A8]

sdir, mypath = sys.argv[1], sys.argv[2]
skip_full = len(sys.argv) > 3 and sys.argv[3] == 'skip_full'
with open(mypath, 'rb') as f:
    MY = pickle.load(f)
for name, ent in MY.items():
    ent['index'] = {st: i for i, st in enumerate(ent['states'])}


def kings_adjacent(a, b):
    return max(abs(FILE[a] - FILE[b]), abs(RANK[a] - RANK[b])) <= 1


def my_value(board):
    """value for the side to move at board: 0 if dead, my table value, else None."""
    if board.is_insufficient_material():
        return 0
    pm = board.piece_map()
    w = ''.join(sorted((p.symbol() for p in pm.values() if p.color == W), key='KQRBNP'.index))
    bl = ''.join(sorted((p.symbol().upper() for p in pm.values() if p.color == B), key='KQRBNP'.index))
    ent = MY.get(w + 'v' + bl)
    if ent is None:
        return None
    sq, used = [], set()
    for c in ent['pieces']:
        for s, p in pm.items():
            if p.symbol() == c and s not in used:
                used.add(s); sq.append(s); break
    j = ent['index'].get((tuple(sq), board.turn))
    return None if j is None else int(ent['val'][j])


def make_board(pieces, sq, turn):
    b = chess.Board(None)
    for c, s in zip(pieces, sq):
        b.set_piece_at(s, chess.Piece.from_symbol(c))
    b.turn = turn
    return b


def legal_state(pieces, sq, turn):
    wk = sq[pieces.index('K')]; bk = sq[pieces.index('k')]
    if len(set(sq)) < len(sq) or kings_adjacent(wk, bk):
        return None
    b = make_board(pieces, sq, turn)
    if b.is_attacked_by(turn, bk if turn == W else wk):
        return None
    return b


def escape_rule(board):
    """Black (bare king) to move after a promotion: 0 if stalemated or a capture reaches a dead position,
    -1 if mated, None otherwise (uncertified)."""
    moves = list(board.legal_moves)
    if not moves:
        return -1 if board.is_check() else 0
    for m in moves:
        if board.is_capture(m):
            board.push(m); dead = board.is_insufficient_material(); board.pop()
            if dead:
                return 0
    return None


# ------------------------------------------------------------------ (a) KPvK corner
def kpvk_corner(file_idx, corner, promo_oracle, label):
    pieces = ['K', 'P', 'k']
    n_states = n_btm = 0; widths = collections.Counter(); viol = []; n_promo = 0
    for wk in range(64):
        for wp in [chess.square(file_idx, r) for r in range(1, 7)]:
            for bk in corner:
                for turn in (W, B):
                    b = legal_state(pieces, (wk, wp, bk), turn)
                    if b is None:
                        continue
                    n_states += 1
                    moves = list(b.legal_moves)
                    if turn == B:
                        n_btm += 1
                        if not moves:
                            if b.is_check():
                                viol.append(('defender mated', b.fen()))
                            else:
                                widths[0] += 1
                            continue
                        internal = 0; bnd = 0
                        for m in moves:
                            if b.is_capture(m):
                                b.push(m); dead = b.is_insufficient_material(); b.pop()
                                if dead:
                                    bnd = 1
                            elif m.to_square in corner:
                                internal += 1
                        wdt = internal + bnd
                        widths[wdt] += 1
                        if wdt == 0:
                            viol.append(('defender has no safe option', b.fen()))
                    else:
                        if not moves:
                            continue   # attacker stalemated: fine
                        for m in moves:
                            if m.promotion:
                                n_promo += 1
                                b.push(m); v = promo_oracle(b); b.pop()   # value for Black to move
                                if v is None or v < 0:
                                    viol.append(('attacker promotion not certified safe', b.fen(), m.uci())); break
                            elif b.is_capture(m):
                                viol.append(('unexpected attacker capture', b.fen(), m.uci())); break
                            # other moves keep BK in the corner and the pawn on its file: successor in the region
    print(f'[{label}] states {n_states} (black-to-move {n_btm}); widths {dict(sorted(widths.items()))}; '
          f'promotions examined {n_promo}; closure violations {len(viol)}')
    for v in viol[:5]:
        print('   ', v)
    return n_states, viol


CORNER_A = {chess.A8, chess.B8, chess.A7, chess.B7}
CORNER_H = {chess.H8, chess.G8, chess.H7, chess.G7}
kpvk_corner(0, CORNER_A, my_value, 'KPvK a-file corner, promotions certified by my tables')
kpvk_corner(0, CORNER_A, escape_rule, 'KPvK a-file corner, promotions certified by the one-ply escape rule only')
kpvk_corner(7, CORNER_H, my_value, 'KPvK h-file mirror, my tables')

# ------------------------------------------------------------------ (b) KBPvK corner
def kbpvk_corner():
    pieces = ['K', 'B', 'P', 'k']
    n_states = n_btm = 0; widths = collections.Counter(); viol = []
    promo_cert = collections.Counter()
    kxb = collections.Counter()
    t0 = time.time()
    for wk in range(64):
        for wb in [s for s in range(64) if DARK[s]]:
            for wp in [chess.square(0, r) for r in range(1, 7)]:
                for bk in CORNER_A:
                    for turn in (W, B):
                        b = legal_state(pieces, (wk, wb, wp, bk), turn)
                        if b is None:
                            continue
                        n_states += 1
                        moves = list(b.legal_moves)
                        if turn == B:
                            n_btm += 1
                            if not moves:
                                if b.is_check():
                                    viol.append(('defender mated', b.fen()))
                                else:
                                    widths[0] += 1
                                continue
                            internal = 0; bnd = 0
                            for m in moves:
                                if b.is_capture(m):
                                    b.push(m)
                                    if b.is_insufficient_material():
                                        bnd = 1                      # Kxp -> KBvK dead
                                    else:
                                        v = my_value(b)              # Kxb -> KPvK, White to move; my table
                                        kxb['safe' if (v is not None and v <= 0) else 'unsafe'] += 1
                                        if v is not None and v <= 0:
                                            bnd = 1
                                    b.pop()
                                elif m.to_square in CORNER_A:
                                    internal += 1
                            wdt = internal + bnd
                            widths[wdt] += 1
                            if wdt == 0:
                                viol.append(('defender has no safe option', b.fen()))
                        else:
                            if not moves:
                                continue
                            for m in moves:
                                if m.promotion:
                                    b.push(m); v = escape_rule(b); b.pop()
                                    promo_cert['certified draw' if v == 0 else ('mate' if v == -1 else 'uncertified')] += 1
                                    if v is None or v < 0:
                                        viol.append(('attacker promotion not certified', b.fen(), m.uci())); break
                                elif b.is_capture(m):
                                    viol.append(('unexpected attacker capture', b.fen(), m.uci())); break
    print(f'[KBPvK dark bishop, a-pawn, BK corner] states {n_states} (black-to-move {n_btm}); widths {dict(sorted(widths.items()))}; '
          f'promotions {dict(promo_cert)}; Kxb boundary {dict(kxb)}; closure violations {len(viol)}; {time.time()-t0:.0f}s')
    for v in viol[:5]:
        print('   ', v)


kbpvk_corner()

# ------------------------------------------------------------------ (c) kbpvk_sub.pkl probe and full Bellman re-check
def decode(key, n):
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


with open(f'{sdir}/kbpvk_sub.pkl', 'rb') as f:
    S = pickle.load(f)
skeys, sval = S['keys'], S['val']
pieces4 = ['K', 'B', 'P', 'k']
sidx = {int(k): i for i, k in enumerate(skeys.tolist())}
print(f'[kbpvk_sub.pkl] {len(skeys)} states; not-lost for Black (BTM val != -1 or WTM val != 1): '
      f'{int(sum(1 for k, v in zip(skeys.tolist(), sval.tolist()) if ((k & 1) and v != 1) or ((not (k & 1)) and v != -1)))}')
t0 = time.time()
probe_candidates = 0; mislabelled = []
for k, v in zip(skeys.tolist(), sval.tolist()):
    if (k & 1) or v != 0:
        continue
    sq, turn = decode(k, 4)
    b = make_board(pieces4, sq, turn)
    moves = list(b.legal_moves)
    if not moves or any(not b.is_capture(m) for m in moves):
        continue
    probe_candidates += 1
    safe = False
    for m in moves:
        b.push(m)
        vv = 0 if b.is_insufficient_material() else my_value(b)
        b.pop()
        if vv is not None and vv <= 0:
            safe = True; break
    if not safe:
        mislabelled.append(b.fen())
print(f'[solve() gap probe] black-to-move draws whose every legal move is a capture: {probe_candidates}; '
      f'of these with every capture losing (would be mislabelled by solve()): {len(mislabelled)} {mislabelled[:3]}; {time.time()-t0:.0f}s')

if not skip_full:
    t0 = time.time()
    badT = badB = 0; examples = []
    for k, v in zip(skeys.tolist(), sval.tolist()):
        sq, turn = decode(k, 4)
        b = make_board(pieces4, sq, turn)
        moves = list(b.legal_moves)
        if not moves:
            exp = -1 if b.is_check() else 0
            if v != exp:
                badT += 1; examples.append(('T', b.fen(), v, exp))
            continue
        best = -2
        for m in moves:
            if b.is_capture(m) or m.promotion:
                b.push(m)
                if b.is_insufficient_material():
                    vv = 0
                elif m.promotion:
                    vv = escape_rule(b)
                    if vv is None:
                        vv = -1          # the flagged ASSUMPTION of the note
                else:
                    vv = my_value(b)     # Kxb -> KPvK (my table)
                b.pop()
                if vv is None:
                    raise RuntimeError('uncertified ' + b.fen())
                best = max(best, -vv)
            else:
                nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
                j = sidx.get(encode(nsq, not turn))
                if j is None:
                    raise RuntimeError('successor missing ' + b.fen())
                best = max(best, -int(sval[j]))
        if best != v:
            badB += 1
            if len(examples) < 5:
                examples.append(('B', b.fen(), v, best))
    print(f'[kbpvk_sub.pkl full (T)+(B) re-check with own move generation + assumption] terminal errors {badT}, '
          f'Bellman violations {badB}; examples {examples[:5]}; {time.time()-t0:.0f}s')

# ------------------------------------------------------------------ (d) 2,000-sample closure check on table_KPvK.pkl
with open(f'{sdir}/table_KPvK.pkl', 'rb') as f:
    T = pickle.load(f)
assert list(T['pieces']) == ['K', 'P', 'k']
tkeys, tval = T['keys'], T['val']
tdict = {int(k): int(v) for k, v in zip(tkeys.tolist(), tval.tolist())}
pieces3 = ['K', 'P', 'k']


def not_lost(key):
    v = tdict[key]
    return (v != 1) if (key & 1) else (v != -1)


nl_keys = [k for k in tdict if not_lost(k)]
print(f'[table_KPvK.pkl] not-lost states {len(nl_keys)} of {len(tdict)}')
random.seed(2026)
sample = random.sample(nl_keys, 2000)
viol = []; agree = 0; myidx = MY['KPvK']['index']; myval = MY['KPvK']['val']
for k in sample:
    sq, turn = decode(k, 3)
    if int(myval[myidx[(sq, turn)]]) == tdict[k]:
        agree += 1
    b = make_board(pieces3, sq, turn)
    moves = list(b.legal_moves)
    if turn == B:
        if not moves:
            if b.is_check():
                viol.append(('mated', b.fen()))
            continue
        ok = False
        for m in moves:
            if b.is_capture(m):
                b.push(m); dead = b.is_insufficient_material(); b.pop()
                if dead:
                    ok = True; break
            else:
                nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
                if not_lost(encode(nsq, W)):
                    ok = True; break
        if not ok:
            viol.append(('defender no safe option', b.fen()))
    else:
        for m in moves:
            if m.promotion:
                b.push(m); vv = my_value(b); b.pop()
                if vv is None or vv < 0:
                    viol.append(('promotion unsafe', b.fen(), m.uci())); break
            elif b.is_capture(m):
                viol.append(('unexpected capture', b.fen())); break
            else:
                nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
                if not not_lost(encode(nsq, B)):
                    viol.append(('attacker escapes to a lost-for-defender state', b.fen(), m.uci())); break
print(f'[2,000-sample closure of the not-lost set of table_KPvK.pkl] violations {len(viol)}; '
      f'pickle value agrees with my independent table on {agree}/2000 sampled states')
for v in viol[:5]:
    print('   ', v)
