"""
traps.py -- certified fortress regions as greatest fixed points (traps) inside candidate predicates.

Experiments
  1. KPvK: (a) the maximal trap of the WHOLE class must equal the exact not-lost set of the retrograde solve
            (independent-path check of the two algorithms);
           (b) the rook-pawn corner region {WP on the a-file, BK in {a8,b8,a7,b7}}: maximal trap inside it,
            its width (number of safe defender options), and comparison with the exact drawn set of the region;
           (c) other structured candidate regions (BK in front of the pawn; BK on the queening square).
  2. KBPvK, dark bishop + a-pawn ("wrong bishop"): the maximal trap inside {BK in {a8,b8,a7,b7}} with the
     boundary oracle = KPvK/KBvK tables + a one-ply escape rule after promotions (certified) or None.
     Also the exact not-lost set of the whole dark-bishop/a-pawn sub-class under the ASSUMPTION that a
     promotion the defender cannot immediately capture into a dead position is lost (flagged; sampled
     against Syzygy separately).
  3. KRvKB, dark bishop: maximal trap with the BK within Chebyshev distance 2 of the light corner a8
     ("safe corner") versus the dark corner a1 ("wrong corner"), boundaries = KRvK table + dead positions.
Run: python3 traps.py <scratch_dir_with_tables>
"""
import sys, pickle, time, collections
import numpy as np
import chess
import chesslib as cl

scr = sys.argv[1]
tables = {}
for name in ('KQvK', 'KRvK', 'KBvK', 'KNvK', 'KPvK'):
    with open(f'{scr}/table_{name}.pkl', 'rb') as f:
        T = pickle.load(f)
    U = cl.Universe(T['pieces'])
    index = {int(k): i for i, k in enumerate(T['keys'])}
    tables[name] = (U, index, T['val'])
    if name == 'KPvK':
        KP = T

A8, B8, A7, B7 = chess.A8, chess.B8, chess.A7, chess.B7
CORNER = {A8, B8, A7, B7}


def describe(U, keys, alive, defender, safe, label):
    turn_white = (keys & 1).astype(bool)
    is_def = (turn_white == (defender == chess.WHITE))
    print(f'  [{label}] trap size {int(alive.sum())}: defender-to-move {int((alive & is_def).sum())}, attacker-to-move {int((alive & ~is_def).sum())}')
    bk = np.array([U.unkey(int(k))[0][U.bk] for k in keys])
    wk = np.array([U.unkey(int(k))[0][U.wk] for k in keys])
    c = collections.Counter(chess.square_name(int(s)) for s in bk[alive])
    print('   defender king squares in the trap:', dict(sorted(c.items())))
    if is_def.any():
        w = safe[alive & is_def]
        hist = collections.Counter(int(x) for x in w)
        print('   width (safe defender options) histogram:', dict(sorted(hist.items())), ' min width', int(w.min()) if len(w) else None)
    return bk, wk


# ---------------- 1. KPvK ----------------
print('=== 1. KPvK ===', flush=True)
U = cl.Universe(['K', 'P', 'k'])
keys = KP['keys']; val = KP['val']
index = {int(k): i for i, k in enumerate(keys)}
G = dict(indptr=KP['indptr'], indices=KP['indices'], term=KP['term'], bnd_best=KP['bnd_best'],
         bnd_unc=np.zeros(len(keys), dtype=np.int32), esc=np.zeros(len(keys), dtype=np.int32), nmoves=KP['nmoves'])
alive, safe = cl.extract_trap(U, keys, G, chess.BLACK, label='whole class')
turn_white = (keys & 1).astype(bool)
not_lost = np.where(turn_white, val != 1, val != -1)
print('  maximal trap of the whole class == exact not-lost set of the solve:', bool((alive == not_lost).all()),
      f'(trap {int(alive.sum())}, not-lost {int(not_lost.sum())})')

# (b) rook-pawn corner region: enumerate the region, build its graph with escapes, extract
def region_corner(sq):
    return cl.FILE[sq[1]] == 0 and sq[2] in CORNER
keys_r, index_r = cl.enumerate_positions(U, square_filter=region_corner)
G_r = cl.build_graph(U, keys_r, index_r, cl.oracle_from_tables(tables), in_region=True)
alive_r, safe_r = cl.extract_trap(U, keys_r, G_r, chess.BLACK, label='corner region')
describe(U, keys_r, alive_r, chess.BLACK, safe_r, 'KPvK a-pawn corner')
# exact values inside the region
vals_r = np.array([val[index[int(k)]] for k in keys_r])
tw_r = (keys_r & 1).astype(bool)
nl_r = np.where(tw_r, vals_r != 1, vals_r != -1)
print(f'  exact not-lost positions in the region: {int(nl_r.sum())} of {len(keys_r)}; trap captures {int((alive_r & nl_r).sum())} of them; '
      f'trap subset of not-lost: {bool((alive_r <= nl_r).all())}')
# which region positions are drawn but outside the trap (need to leave the region)?
out = np.nonzero(nl_r & ~alive_r)[0]
print(f'  drawn-but-outside-trap: {len(out)}; examples:', [U.fen(int(keys_r[i])) for i in out[:4]])
# the positions of the region that are lost for Black
lost = np.nonzero(~nl_r)[0]
print(f'  lost positions in the region: {len(lost)}; examples:', [U.fen(int(keys_r[i])) for i in lost[:4]])
# (b2) the same with both rook files and both colours by symmetry is immediate; now the h-file as a control
def region_corner_h(sq):
    return cl.FILE[sq[1]] == 7 and sq[2] in {chess.H8, chess.G8, chess.H7, chess.G7}
keys_h, index_h = cl.enumerate_positions(U, square_filter=region_corner_h, verbose=False)
G_h = cl.build_graph(U, keys_h, index_h, cl.oracle_from_tables(tables), in_region=True, verbose=False)
alive_h, _ = cl.extract_trap(U, keys_h, G_h, chess.BLACK, label='h-file mirror control', verbose=True)
print('  mirror control equal size:', int(alive_h.sum()) == int(alive_r.sum()))

# (c) other candidate regions
def region_front(sq):   # BK directly in front of the pawn (same file, one rank above)
    p, k = sq[1], sq[2]
    return cl.FILE[k] == cl.FILE[p] and cl.RANK[k] == cl.RANK[p] + 1
def region_front2(sq):  # BK on the pawn's file, in front, at distance 1 or 2
    p, k = sq[1], sq[2]
    return cl.FILE[k] == cl.FILE[p] and 1 <= cl.RANK[k] - cl.RANK[p] <= 2
def region_queening(sq):  # BK on the queening square of the pawn's file
    p, k = sq[1], sq[2]
    return cl.FILE[k] == cl.FILE[p] and cl.RANK[k] == 7
def region_box(sq):  # BK within the 'square of the pawn' (rule of the square), pawn not on the 7th
    p, k = sq[1], sq[2]
    d = 7 - cl.RANK[p]
    return max(abs(cl.FILE[k] - cl.FILE[p]), 7 - cl.RANK[k]) <= d and cl.RANK[k] > cl.RANK[p]
for lab, pred in (('BK directly in front', region_front), ('BK in front, distance <= 2', region_front2),
                  ('BK on queening square', region_queening), ('BK inside the square of the pawn, ahead of it', region_box)):
    k2, i2 = cl.enumerate_positions(U, square_filter=pred, verbose=False)
    G2 = cl.build_graph(U, k2, i2, cl.oracle_from_tables(tables), in_region=True, verbose=False)
    a2, s2 = cl.extract_trap(U, k2, G2, chess.BLACK, label=lab, verbose=False)
    v2 = np.array([val[index[int(k)]] for k in k2]); t2 = (k2 & 1).astype(bool)
    nl2 = np.where(t2, v2 != 1, v2 != -1)
    print(f'  region "{lab}": {len(k2)} positions, exact not-lost {int(nl2.sum())}, maximal trap inside {int(a2.sum())}, '
          f'defender-to-move in trap {int((a2 & ~t2).sum())}, min width {int(s2[a2 & ~t2].min()) if (a2 & ~t2).any() else None}')

# ---------------- 2. KBPvK wrong bishop ----------------
print('\n=== 2. KBPvK, dark-squared bishop, a-pawn ===', flush=True)
U2 = cl.Universe(['K', 'B', 'P', 'k'])

def escape_oracle(board):
    """Black (bare king) to move after a promotion: 0 if a capture reaches a dead position or stalemate; -1 if mated; else None."""
    if board.turn != chess.BLACK:
        return None
    moves = list(board.legal_moves)
    if not moves:
        return -1 if board.is_check() else 0
    for m in moves:
        if board.is_capture(m):
            board.push(m)
            dead = board.is_insufficient_material()
            board.pop()
            if dead:
                return 0
    return None

def region_wb(sq):
    k_, b_, p_, bk_ = sq
    return cl.IS_DARK[b_] and cl.FILE[p_] == 0 and bk_ in CORNER
keys2, index2 = cl.enumerate_positions(U2, square_filter=region_wb)
G2 = cl.build_graph(U2, keys2, index2, cl.oracle_from_tables(tables, fallback=escape_oracle), in_region=True)
alive2, safe2 = cl.extract_trap(U2, keys2, G2, chess.BLACK, label='wrong-bishop corner region')
describe(U2, keys2, alive2, chess.BLACK, safe2, 'KBPvK corner')
t2w = (keys2 & 1).astype(bool)
dead2 = np.nonzero(~alive2)[0]
print(f'  positions of the region outside the trap: {len(dead2)}; examples:', [U2.fen(int(keys2[i])) for i in dead2[:6]])
# classify the dead ones: is the defender to move? in check? what else
cnt = collections.Counter()
for i in dead2:
    sq, t = U2.unkey(int(keys2[i]))
    b = U2.board(sq, t)
    cnt[('BTM' if t == chess.BLACK else 'WTM', 'check' if b.is_check() else 'nocheck', 'BK ' + chess.square_name(sq[3]), 'P ' + chess.square_name(sq[2]))] += 1
print('  dead by (side, check, BK, pawn):', dict(cnt.most_common(12)))
# exact not-lost set of the whole dark-bishop/a-pawn sub-class under the promotion assumption
def assume_oracle(board):
    v = escape_oracle(board)
    return -1 if v is None else v
def region_sub(sq):
    k_, b_, p_, bk_ = sq
    return cl.IS_DARK[b_] and cl.FILE[p_] == 0
keys3, index3 = cl.enumerate_positions(U2, square_filter=region_sub)
G3 = cl.build_graph(U2, keys3, index3, cl.oracle_from_tables(tables, fallback=assume_oracle), in_region=True)
print('  escapes (should be 0 since the sub-class is closed):', int(G3['esc'].sum()))
val3, dtm3 = cl.solve(U2, keys3, G3)
t3 = (keys3 & 1).astype(bool)
nl3 = np.where(t3, val3 != 1, val3 != -1)
bk3 = np.array([U2.unkey(int(k))[0][3] for k in keys3])
print(f'  sub-class: {len(keys3)} positions; not-lost for Black {int(nl3.sum())} ({100*nl3.mean():.2f}%)')
c3 = collections.Counter(chess.square_name(int(s)) for s in bk3[nl3])
print('   not-lost by BK square (top 12):', c3.most_common(12))
# the corner region inside the exact solve
in_corner = np.array([s in CORNER for s in bk3])
print(f'   corner-region positions {int(in_corner.sum())}, not-lost among them {int((nl3 & in_corner).sum())}, trap found {int(alive2.sum())}')
# boundary positions used by the assumption: count how many White-to-move positions have an uncertified-type promotion
with open(f'{scr}/kbpvk_sub.pkl', 'wb') as f:
    pickle.dump(dict(keys=keys3, val=val3, bk=bk3), f, protocol=4)

# ---------------- 3. KRvKB corners ----------------
print('\n=== 3. KRvKB, dark-squared black bishop ===', flush=True)
U3 = cl.Universe(['K', 'R', 'k', 'b'])
def near(sq, corner, d=2):
    return max(abs(cl.FILE[sq] - cl.FILE[corner]), abs(cl.RANK[sq] - cl.RANK[corner])) <= d
for corner, lab in ((chess.A8, 'light corner a8 (safe)'), (chess.A1, 'dark corner a1 (wrong)')):
    def region(sq, corner=corner):
        return cl.IS_DARK[sq[3]] and near(sq[2], corner)
    kk, ii = cl.enumerate_positions(U3, square_filter=region)
    GG = cl.build_graph(U3, kk, ii, cl.oracle_from_tables(tables), in_region=True)
    aa, ss = cl.extract_trap(U3, kk, GG, chess.BLACK, label=lab)
    bk, wk = describe(U3, kk, aa, chess.BLACK, ss, 'KRvKB ' + lab)
    tw = (kk & 1).astype(bool)
    with open(f'{scr}/krvkb_{corner}.pkl', 'wb') as f:
        pickle.dump(dict(keys=kk, alive=aa, safe=ss), f, protocol=4)
    # distance-1 sub-region count inside the trap
    d1 = np.array([near(int(s), corner, 1) for s in bk])
    print(f'   trap members with BK within distance 1 of the corner: {int((aa & d1).sum())} of {int(d1.sum())}')
print('done')
