"""audit_A_misc.py -- audit A: (1) DTM of the cited cycling-start positions from MY tables; (2) independent recount of
the visibility census (largest sealed region with p pawns on ranks 2-7, p = 1..4) with own component code.
Usage: python3 audit_A_misc.py <my_tables.pkl>"""
import sys, pickle, itertools, chess
with open(sys.argv[1], 'rb') as f:
    MY = pickle.load(f)
def lookup(name, fen):
    ent = MY[name]; idx = {st: i for i, st in enumerate(ent['states'])}
    b = chess.Board(fen); pm = b.piece_map(); sq = []
    for c in ent['pieces']:
        for s, p in pm.items():
            if p.symbol() == c: sq.append(s); break
    i = idx[(tuple(sq), b.turn)]
    return int(ent['val'][i]), int(ent['dtm'][i])
for name, fen in (('KRvK', '8/8/8/8/8/8/8/k1KR4 w - - 0 1'), ('KRvK', '8/8/8/8/8/8/8/k3KR2 w - - 0 1'), ('KRvK', '8/8/8/8/8/3K4/1R6/k7 w - - 0 1'), ('KQvK', '8/8/8/8/8/8/8/k2KQ3 w - - 0 1'), ('KQvK', '8/8/8/8/8/8/8/k3KQ2 w - - 0 1')):
    v, d = lookup(name, fen)
    print(f'{name} {fen}: value {v}, DTM {d} plies (mate in {(d+1)//2})')
# visibility census
F = [chess.square_file(s) for s in range(64)]; R = [chess.square_rank(s) for s in range(64)]
DARK = [((F[s] + R[s]) % 2 == 0) for s in range(64)]
def comps(P, dirs, verts):
    verts = [s for s in verts if s not in P]; vs = set(verts); seen = set(); out = []
    for s in verts:
        if s in seen: continue
        stack = [s]; seen.add(s); c = 0
        while stack:
            u = stack.pop(); c += 1
            for df, dr in dirs:
                f, r = F[u] + df, R[u] + dr
                if 0 <= f < 8 and 0 <= r < 8:
                    t = chess.square(f, r)
                    if t in vs and t not in seen:
                        seen.add(t); stack.append(t)
        out.append(c)
    return sorted(out, reverse=True)
ROOK = [(1,0),(-1,0),(0,1),(0,-1)]; BISH = [(1,1),(1,-1),(-1,1),(-1,-1)]
pawn_sq = [s for s in range(64) if 1 <= R[s] <= 6]
for kind, dirs, verts in (('rook', ROOK, list(range(64))), ('dark bishop', BISH, [s for s in range(64) if DARK[s]]), ('queen', ROOK + BISH, list(range(64)))):
    row = []
    for p in range(1, 5):
        best = 0; wit = None
        for P in itertools.combinations(pawn_sq, p):
            cs = comps(set(P), dirs, verts)
            second = cs[1] if len(cs) > 1 else 0
            if second > best: best, wit = second, P
        row.append(f'p={p}: {best}' + (f' (e.g. {" ".join(chess.square_name(s) for s in wit)})' if wit else ''))
    print(kind, '|', ' | '.join(row), flush=True)
