"""
chesslib.py -- shared machinery for the chess-as-a-proof-system note (opus-2026-10-12-S23).

Objects
-------
* Universe: a material class, an ordered tuple of piece symbols (white upper case), e.g. ('K','R','k','b').
  Positions are encoded as integers: key = turn + 2 * sum(sq_i * 64**i).  Identical pieces are kept
  in sorted square order so the key is canonical.
* enumerate_positions(U, square_filter): all placements that python-chess accepts as legal in the
  tablebase sense (kings not adjacent, side not to move not in check, no pawn on rank 1/8), both
  sides to move, optionally restricted by a predicate on the square tuple and one on the board.
* build_graph(U, keys, index, boundary): CSR successor lists inside the class plus, for each
  position, the best value obtainable by a move that leaves the class (capture / promotion),
  evaluated by the caller's boundary oracle; the oracle returns the value for the side to move
  AFTER the move (+1 win, 0 draw, -1 loss, None = uncertified).
* solve(...): exact retrograde WDL+DTM inside the class given a complete boundary oracle
  (least fixed point for wins; everything unresolved is a draw).  Values are from the viewpoint of
  the side to move: +1 / 0 / -1.  DTM in plies to mate (0 = mated).
* extract_trap(...): the greatest fixed point of the defender's safety operator inside a candidate
  region: the largest F subset of the region such that
      defender to move:  stalemated, or some move stays in F, or some move reaches a certified boundary;
      attacker to move:  every move stays in F or reaches a certified boundary (mated attacker is fine).
  Only the move generator is used; no tablebase.  Survivors are PROVED >= draw for the defender
  (Knaster-Tarski / Zielonka trap argument), relative to the certified boundary.
"""
import itertools, sys, time, collections
import numpy as np
import chess

SQ = list(range(64))
FILE = [chess.square_file(s) for s in SQ]
RANK = [chess.square_rank(s) for s in SQ]
KING_ADJ = [[t for t in SQ if t != s and max(abs(FILE[s]-FILE[t]), abs(RANK[s]-RANK[t])) <= 1] for s in SQ]
IS_DARK = [((FILE[s] + RANK[s]) % 2 == 0) for s in SQ]   # a1 dark
ORDER = 'KQRBNP'


def sqname(s):
    return chess.square_name(s)


def class_name(pieces):
    w = ''.join(sorted((c for c in pieces if c.isupper()), key=lambda c: ORDER.index(c)))
    b = ''.join(sorted((c.upper() for c in pieces if c.islower()), key=lambda c: ORDER.index(c)))
    return w + 'v' + b


class Universe:
    def __init__(self, pieces):
        self.pieces = tuple(pieces)
        self.n = len(self.pieces)
        self.pobj = [chess.Piece.from_symbol(c) for c in self.pieces]
        groups = collections.defaultdict(list)
        for i, c in enumerate(self.pieces):
            groups[c].append(i)
        self.groups = [g for g in groups.values() if len(g) > 1]
        self.name = class_name(self.pieces)
        self.wk = self.pieces.index('K')
        self.bk = self.pieces.index('k')
        self.pawn_idx = [i for i, c in enumerate(self.pieces) if c in 'Pp']

    def canon(self, squares):
        sq = list(squares)
        for g in self.groups:
            vals = sorted(sq[i] for i in g)
            for i, v in zip(g, vals):
                sq[i] = v
        return tuple(sq)

    def key(self, squares, turn):
        k = 0
        for s in reversed(self.canon(squares)):
            k = k * 64 + s
        return k * 2 + (1 if turn == chess.WHITE else 0)

    def unkey(self, key):
        turn = chess.WHITE if (key & 1) else chess.BLACK
        key >>= 1
        sq = []
        for _ in range(self.n):
            sq.append(key % 64)
            key //= 64
        return tuple(sq), turn

    def board(self, squares, turn):
        b = chess.Board(None)
        for p, s in zip(self.pobj, squares):
            b.set_piece_at(s, p)
        b.turn = turn
        return b

    def fen(self, key):
        sq, t = self.unkey(key)
        return self.board(sq, t).fen()


def enumerate_positions(U, square_filter=None, board_filter=None, turns=(chess.WHITE, chess.BLACK), verbose=True):
    """Return (keys ndarray int64 sorted, index dict key->i)."""
    t0 = time.time()
    keys = []
    ranges = []
    for c in U.pieces:
        if c in 'Pp':
            ranges.append([s for s in SQ if 1 <= RANK[s] <= 6])
        else:
            ranges.append(SQ)
    b = chess.Board(None)
    for squares in itertools.product(*ranges):
        if len(set(squares)) < U.n:
            continue
        if squares[U.bk] in KING_ADJ[squares[U.wk]]:
            continue
        sq = U.canon(squares)
        if sq != squares:
            continue
        if square_filter is not None and not square_filter(sq):
            continue
        b.clear_board()
        for p, s in zip(U.pobj, sq):
            b.set_piece_at(s, p)
        for turn in turns:
            b.turn = turn
            if b.was_into_check():
                continue
            if board_filter is not None and not board_filter(b):
                continue
            keys.append(U.key(sq, turn))
    keys = np.array(sorted(keys), dtype=np.int64)
    index = {int(k): i for i, k in enumerate(keys)}
    if verbose:
        print(f"[{U.name}] enumerated {len(keys)} positions in {time.time()-t0:.1f}s", flush=True)
    return keys, index


def build_graph(U, keys, index, boundary, verbose=True, in_region=False):
    """
    For every position: internal successors (CSR), terminal status, boundary summary.
    boundary(board_after_move, move) -> value for the side to move after the move (+1/0/-1) or None.
    in_region: if True, internal successors that are not enumerated count as escapes (esc[i]).
    Arrays: indptr, indices; term (0 normal, 1 mated, 2 stalemated); bnd_best (max over boundary moves
    of -value, i.e. best for the mover; -2 if no certified boundary move); bnd_unc (number of
    uncertified boundary moves); esc; nmoves.
    """
    t0 = time.time()
    N = len(keys)
    indptr = np.zeros(N + 1, dtype=np.int64)
    indices_list = []
    term = np.zeros(N, dtype=np.int8)
    bnd_best = np.full(N, -2, dtype=np.int8)
    bnd_unc = np.zeros(N, dtype=np.int32)
    esc = np.zeros(N, dtype=np.int32)
    nmoves = np.zeros(N, dtype=np.int32)
    b = chess.Board(None)
    for i in range(N):
        sq, turn = U.unkey(int(keys[i]))
        b.clear_board()
        for p, s in zip(U.pobj, sq):
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
                v = boundary(b, m)
                b.pop()
                if v is None:
                    bnd_unc[i] += 1
                else:
                    best = max(best, -v)
                continue
            mover = sq.index(m.from_square)
            nsq = list(sq)
            nsq[mover] = m.to_square
            k2 = U.key(nsq, not turn)
            j = index.get(k2)
            if j is None:
                if in_region:
                    esc[i] += 1
                    continue
                raise RuntimeError(f"successor not enumerated: {U.fen(int(keys[i]))} move {m}")
            succ.append(j)
        bnd_best[i] = best
        indices_list.append(np.array(succ, dtype=np.int32))
        indptr[i + 1] = indptr[i] + len(succ)
        if verbose and i % 250000 == 0 and i:
            print(f"  ... {i}/{N} ({time.time()-t0:.0f}s)", flush=True)
    indices = np.concatenate(indices_list) if indices_list else np.zeros(0, dtype=np.int32)
    if verbose:
        print(f"[{U.name}] graph: {N} positions, {len(indices)} internal edges, "
              f"{int(bnd_unc.sum())} uncertified boundary moves, {int(esc.sum())} escapes; {time.time()-t0:.1f}s", flush=True)
    return dict(indptr=indptr, indices=indices, term=term, bnd_best=bnd_best, bnd_unc=bnd_unc, esc=esc, nmoves=nmoves)


def predecessors(G, N):
    indptr, indices = G['indptr'], G['indices']
    src = np.repeat(np.arange(N, dtype=np.int32), np.diff(indptr))
    order = np.argsort(indices, kind='stable')
    pred_indices = src[order]
    counts = np.bincount(indices, minlength=N)
    pred_indptr = np.zeros(N + 1, dtype=np.int64)
    pred_indptr[1:] = np.cumsum(counts)
    return pred_indptr, pred_indices


def solve(U, keys, G, verbose=True):
    """Exact retrograde WDL + DTM (plies) from the side-to-move viewpoint. Requires bnd_unc == 0 everywhere.
    Boundary wins are recorded with DTM 1 (the boundary DTM is not tracked); WDL is exact."""
    t0 = time.time()
    N = len(keys)
    if G['bnd_unc'].any():
        raise ValueError("uncertified boundary moves present; solve needs a complete boundary oracle")
    indptr, indices, term, bnd_best = G['indptr'], G['indices'], G['term'], G['bnd_best']
    val = np.zeros(N, dtype=np.int8)
    dtm = np.full(N, -1, dtype=np.int32)
    resolved = np.zeros(N, dtype=bool)
    pred_indptr, pred_indices = predecessors(G, N)
    n_unres = np.diff(indptr).astype(np.int64)
    has_nonlosing = (bnd_best >= 0)
    frontier = []
    for i in np.nonzero(term == 1)[0]:
        val[i] = -1; dtm[i] = 0; resolved[i] = True; frontier.append(int(i))
    for i in np.nonzero(term == 2)[0]:
        val[i] = 0; dtm[i] = 0; resolved[i] = True
    for i in np.nonzero((term == 0) & (bnd_best == 1))[0]:
        val[i] = 1; dtm[i] = 1; resolved[i] = True; frontier.append(int(i))
    while frontier:
        nxt = []
        for j in frontier:
            vj = val[j]
            for p in pred_indices[pred_indptr[j]:pred_indptr[j+1]]:
                if resolved[p]:
                    continue
                if vj == -1:
                    val[p] = 1; dtm[p] = dtm[j] + 1; resolved[p] = True; nxt.append(int(p))
                else:
                    n_unres[p] -= 1
                    if n_unres[p] == 0 and not has_nonlosing[p]:
                        val[p] = -1; dtm[p] = dtm[j] + 1; resolved[p] = True; nxt.append(int(p))
        frontier = nxt
    if verbose:
        nw = int((val == 1).sum()); nl = int((val == -1).sum()); nd = N - nw - nl
        mx = int(dtm[val != 0].max()) if (val != 0).any() else -1
        print(f"[{U.name}] solved: win {nw}, draw {nd}, loss {nl}; max DTM {mx} plies; {time.time()-t0:.1f}s", flush=True)
    return val, dtm


def extract_trap(U, keys, G, defender, verbose=True, label=''):
    """
    Greatest fixed point of the defender's safety operator inside the enumerated region.
    Returns (alive, safe): alive = membership in the maximal trap; safe = number of safe defender
    options (alive successors plus one if a certified boundary move exists) at defender-to-move nodes.
    Certification of boundary moves: the oracle value is for the side to move after the move, and
    bnd_best = max over boundary moves of (-value) = best outcome for the mover.
      defender to move: a boundary move is safe iff bnd_best >= 0 (the mover gets >= draw).
      attacker to move: closure needs every boundary move certified with mover value <= 0,
        i.e. bnd_unc == 0 and bnd_best <= 0.
    """
    t0 = time.time()
    N = len(keys)
    indptr, indices, term, bnd_best, bnd_unc, esc = G['indptr'], G['indices'], G['term'], G['bnd_best'], G['bnd_unc'], G['esc']
    turn_white = (keys & 1).astype(bool)
    is_def = (turn_white == (defender == chess.WHITE))
    alive = np.ones(N, dtype=bool)
    pred_indptr, pred_indices = predecessors(G, N)
    nsucc = np.diff(indptr)
    safe = nsucc.astype(np.int64).copy()
    safe += (bnd_best >= 0).astype(np.int64)
    dead_list = []
    for i in range(N):
        if is_def[i]:
            if term[i] == 1:
                alive[i] = False; dead_list.append(i)
            elif term[i] == 2:
                pass
            elif safe[i] == 0:
                alive[i] = False; dead_list.append(i)
        else:
            if term[i] != 0:
                pass
            elif esc[i] > 0 or bnd_unc[i] > 0 or bnd_best[i] > 0:
                alive[i] = False; dead_list.append(i)
    while dead_list:
        nxt = []
        for j in dead_list:
            for p in pred_indices[pred_indptr[j]:pred_indptr[j+1]]:
                if not alive[p]:
                    continue
                if is_def[p]:
                    safe[p] -= 1
                    if safe[p] == 0:
                        alive[p] = False; nxt.append(int(p))
                else:
                    alive[p] = False; nxt.append(int(p))
        dead_list = nxt
    if verbose:
        print(f"[{U.name}{(' '+label) if label else ''}] trap: {int(alive.sum())} of {N} positions survive ({time.time()-t0:.1f}s)", flush=True)
    return alive, safe


def verify_trap(U, keys, G, defender, member):
    """Closure check of a GIVEN region (member: bool array). Returns violating indices (empty = certified)."""
    N = len(keys)
    indptr, indices, term, bnd_best, bnd_unc, esc = G['indptr'], G['indices'], G['term'], G['bnd_best'], G['bnd_unc'], G['esc']
    turn_white = (keys & 1).astype(bool)
    is_def = (turn_white == (defender == chess.WHITE))
    bad = []
    for i in np.nonzero(member)[0]:
        succ = indices[indptr[i]:indptr[i+1]]
        if is_def[i]:
            if term[i] == 1:
                bad.append(int(i)); continue
            if term[i] == 2:
                continue
            if bnd_best[i] >= 0:
                continue
            if not member[succ].any():
                bad.append(int(i))
        else:
            if term[i] != 0:
                continue
            if esc[i] > 0 or bnd_unc[i] > 0 or bnd_best[i] > 0 or not member[succ].all():
                bad.append(int(i))
    return bad


def oracle_from_tables(tables, fallback=None):
    """tables: dict class-name -> (U, index, val).  Boundary oracle: insufficient material -> 0,
    solved table -> its value for the side to move, else fallback(board) or None."""
    def oracle(board, move):
        if board.is_insufficient_material():
            return 0
        pcs, sqs = [], []
        for s, p in board.piece_map().items():
            pcs.append(p.symbol()); sqs.append(s)
        name = class_name(pcs)
        ent = tables.get(name)
        if ent is None:
            return fallback(board) if fallback else None
        U2, index2, val2 = ent
        order, used = [], set()
        for c in U2.pieces:
            for t, (c2, s2) in enumerate(zip(pcs, sqs)):
                if c2 == c and t not in used:
                    used.add(t); order.append(s2); break
        j = index2.get(U2.key(order, board.turn))
        if j is None:
            return fallback(board) if fallback else None
        return int(val2[j])
    return oracle


def solve_class(pieces, tables, verbose=True, fallback=None):
    """Enumerate, build and solve a class; register in tables; return (U, keys, index, G, val, dtm)."""
    U = Universe(pieces)
    keys, index = enumerate_positions(U, verbose=verbose)
    G = build_graph(U, keys, index, oracle_from_tables(tables, fallback), verbose=verbose)
    val, dtm = solve(U, keys, G, verbose=verbose)
    tables[U.name] = (U, index, val)
    return U, keys, index, G, val, dtm
