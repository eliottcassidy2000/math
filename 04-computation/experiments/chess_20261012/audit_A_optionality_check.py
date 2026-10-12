"""
audit_A_optionality_check.py -- audit A: independent re-simulation of the optionality witnesses (section 6 / THM-4613
statement 3) with OWN code on MY independent tables (audit_A_solver.py), reproducing the session's tie-break
protocol: White picks the candidate maximising the secondary score, ties broken by the smallest session key
(key = 2*sum(sq_i 64^i) + turn, pieces in table order); Black replies DTM-maximising with ties broken by
python-chess move-generation order; a cycle = a repeated white-to-move state.
Usage: python3 audit_A_optionality_check.py <my_tables.pkl> [KRvK|KQvK] [O1|O2|nonincreasing|optimal|all]
"""
import sys, pickle, time, collections
import numpy as np
import chess

W, B = chess.WHITE, chess.BLACK
mypath = sys.argv[1]
which = sys.argv[2] if len(sys.argv) > 2 else 'KRvK'
variant = sys.argv[3] if len(sys.argv) > 3 else 'all'
with open(mypath, 'rb') as f:
    MY = pickle.load(f)
ent = MY[which]
pieces, states, val, dtm = ent['pieces'], ent['states'], ent['val'], ent['dtm']
index = {st: i for i, st in enumerate(states)}
N = len(states)


def key_of(sq, turn):
    k = 0
    for s in reversed(sq):
        k = k * 64 + s
    return k * 2 + (1 if turn == W else 0)


def board_of(i):
    sq, turn = states[i]
    b = chess.Board(None)
    for c, s in zip(pieces, sq):
        b.set_piece_at(s, chess.Piece.from_symbol(c))
    b.turn = turn
    return b


# successor lists (internal, non-capture) in python-chess generation order
t0 = time.time()
succ = [None] * N
nmoves = np.zeros(N, dtype=np.int32)
for i in range(N):
    sq, turn = states[i]
    b = board_of(i)
    out = []
    moves = list(b.legal_moves)
    nmoves[i] = len(moves)
    for m in moves:
        if b.is_capture(m) or m.promotion:
            continue
        nsq = list(sq); nsq[sq.index(m.from_square)] = m.to_square
        out.append(index[(tuple(nsq), not turn)])
    succ[i] = out
print(f'[{which}] successor lists built in {time.time()-t0:.0f}s', flush=True)
keys = np.array([key_of(*states[i]) for i in range(N)], dtype=np.int64)
tw = np.array([st[1] == W for st in states])
nwin = np.zeros(N, dtype=np.int32)
for i in np.nonzero(tw & (val == 1))[0]:
    nwin[i] = sum(1 for j in succ[i] if val[j] == -1)


def black_reply(j):
    s = succ[j]
    if not s:
        return None
    best = s[0]
    for x in s[1:]:
        if dtm[x] > dtm[best]:
            best = x
    return best


def score(j, kind):
    s = succ[j]
    if not s:
        return 10 ** 6
    return int(min(nmoves[x] for x in s)) if kind == 'O1' else int(min(nwin[x] for x in s))


def run(kind, restrict, max_plies=200):
    t0 = time.time()
    starts = np.nonzero(tw & (val == 1))[0]
    mated = cycled = timeouts = 0; longest = 0; first_cycle = None
    for i0 in starts:
        i = int(i0); seen = {i}; plies = 0; ok = None
        while True:
            cand = [j for j in succ[i] if val[j] == -1]
            if restrict == 'nonincreasing':
                cand = [j for j in cand if dtm[j] - 1 <= dtm[i]]
            elif restrict == 'optimal':
                cand = [j for j in cand if dtm[j] + 1 == dtm[i]]
            best = max(cand, key=lambda j: (score(j, kind), -int(keys[j])))
            plies += 1
            if not succ[best]:
                ok = True; break
            r = black_reply(best); plies += 1
            if r in seen:
                ok = False
                if first_cycle is None:
                    first_cycle = int(i0)
                break
            seen.add(r); i = r
            if plies > max_plies:
                ok = None; break
        if ok is True:
            mated += 1; longest = max(longest, plies)
        elif ok is False:
            cycled += 1
        else:
            timeouts += 1
    lab = f'{kind} among {restrict or "win-preserving"} moves'
    ex = f'; first cycling start {board_of(first_cycle).fen()}' if first_cycle is not None else ''
    print(f'[{which}] {lab}: starts {len(starts)}, mated {mated}, cycled {cycled}, timeouts {timeouts}, longest win {longest} plies ({time.time()-t0:.0f}s){ex}', flush=True)


runs = {'O1': ('O1', None), 'O2': ('O2', None), 'nonincreasing': ('O2', 'nonincreasing'), 'optimal': ('O2', 'optimal')}
for k, (kind, restrict) in runs.items():
    if variant == 'all' or variant == k:
        run(kind, restrict)
