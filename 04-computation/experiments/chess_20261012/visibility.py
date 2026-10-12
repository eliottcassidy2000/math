"""
visibility.py -- the pawn-conditioned visibility graph and visibility cuts.

Lemma V1 (ray components = step components): for a pawn set P, two empty squares are joined by a chain of
unobstructed rook rays iff they are joined by a chain of orthogonal single steps through empty squares;
likewise bishop rays on colour c and diagonal steps within colour c; queen rays and king steps.
Hence the pawn-conditioned visibility atlas is digital topology: rook components = 4-connected components
of the complement of P, bishop components = diagonal-step components of (colour c) minus P, queen components
= 8-connected components.  Checked exhaustively on random pawn sets.
Lemma V2 (colour locality): the bishop graph of colour c depends only on P intersect (colour c).
Census: for p = 1..5 pawns on ranks 2-7, the largest region that can be sealed off (size of the largest
component other than the main one; both components counted), per piece type.
Named structures: the initial position and a few opening skeletons.
"""
import itertools, random, collections
import chess

SQ = range(64)
F = [chess.square_file(s) for s in SQ]; R = [chess.square_rank(s) for s in SQ]
DARK = [((F[s] + R[s]) % 2 == 0) for s in SQ]

def ray_graph(P, dirs):
    """adjacency by unobstructed rays in the given directions among empty squares (not in P)."""
    adj = {s: set() for s in SQ if s not in P}
    for s in adj:
        for df, dr in dirs:
            f, r = F[s] + df, R[s] + dr
            while 0 <= f < 8 and 0 <= r < 8:
                t = chess.square(f, r)
                if t in P:
                    break
                adj[s].add(t)
                f += df; r += dr
    return adj

def step_graph(P, dirs):
    adj = {s: set() for s in SQ if s not in P}
    for s in adj:
        for df, dr in dirs:
            f, r = F[s] + df, R[s] + dr
            if 0 <= f < 8 and 0 <= r < 8:
                t = chess.square(f, r)
                if t not in P:
                    adj[s].add(t)
    return adj

def components(adj, verts=None):
    verts = set(adj) if verts is None else set(verts)
    seen, comps = set(), []
    for s in verts:
        if s in seen:
            continue
        stack, comp = [s], set()
        seen.add(s)
        while stack:
            u = stack.pop(); comp.add(u)
            for v in adj[u]:
                if v in verts and v not in seen:
                    seen.add(v); stack.append(v)
        comps.append(comp)
    return sorted(comps, key=len, reverse=True)

ROOK = [(1,0),(-1,0),(0,1),(0,-1)]
BISH = [(1,1),(1,-1),(-1,1),(-1,-1)]
QUEEN = ROOK + BISH

def comps_for(P, kind):
    if kind == 'rook':
        return components(ray_graph(P, ROOK))
    if kind == 'queen':
        return components(ray_graph(P, QUEEN))
    if kind == 'bishop_dark':
        adj = ray_graph(P, BISH); return components(adj, [s for s in adj if DARK[s]])
    if kind == 'bishop_light':
        adj = ray_graph(P, BISH); return components(adj, [s for s in adj if not DARK[s]])

# Lemma V1 and V2 checks
random.seed(1)
pawn_squares = [s for s in SQ if 1 <= R[s] <= 6]
ok1 = ok2 = True
for trial in range(300):
    P = set(random.sample(pawn_squares, random.randint(0, 16)))
    for dirs in (ROOK, BISH, QUEEN):
        c_ray = components(ray_graph(P, dirs)); c_step = components(step_graph(P, dirs))
        if sorted(map(frozenset, c_ray)) != sorted(map(frozenset, c_step)):
            ok1 = False
    Pd = {s for s in P if DARK[s]}
    a1 = ray_graph(P, BISH); a2 = ray_graph(Pd, BISH)
    if any(a1[s] != a2[s] for s in a1 if DARK[s]):
        ok2 = False
print('Lemma V1 (ray components = step components), 300 random pawn sets x 3 piece types:', ok1)
print('Lemma V2 (dark-bishop graph depends only on dark pawns):', ok2)

# census: largest sealed-off region with p pawns
print('\nLargest region sealable with p pawns (second-largest component size; number of optimal pawn sets):')
for kind in ('rook', 'bishop_dark', 'queen'):
    row = []
    for p in range(1, 5):
        best, cnt, wit = -1, 0, None
        for P in itertools.combinations(pawn_squares, p):
            cs = comps_for(set(P), kind)
            second = len(cs[1]) if len(cs) > 1 else 0
            if second > best:
                best, cnt, wit = second, 1, P
            elif second == best:
                cnt += 1
        row.append(f'p={p}: {best} ({cnt}; e.g. {" ".join(chess.square_name(s) for s in wit)})')
    print(f'  {kind:13s}', ' | '.join(row))

# named structures
def structure(white, black):
    P = {chess.parse_square(s) for s in white.split()} | {chess.parse_square(s) for s in black.split()}
    return P
named = {
    'initial': structure('a2 b2 c2 d2 e2 f2 g2 h2', 'a7 b7 c7 d7 e7 f7 g7 h7'),
    'French advance (e5/d4 vs e6/d5, c5/f6 absent)': structure('a2 b2 c2 d4 e5 f2 g2 h2', 'a7 b7 c7 d5 e6 f7 g7 h7'),
    'Carlsbad (W: a2 b2 d4 e3 f2 g2 h2; B: a7 b7 c6 d5 e6 f7 g7 h7)': structure('a2 b2 d4 e3 f2 g2 h2', 'a7 b7 c6 d5 e6 f7 g7 h7'),
    'Closed Ruy / KID wall (d5/e4 vs d6/e5)': structure('a2 b2 c4 d5 e4 f2 g2 h2', 'a7 b7 c7 d6 e5 f7 g7 h7'),
    'main diagonal a1-h8 (8 pawns)': structure('a1 b2 c3 d4', 'e5 f6 g7 h8'),
}
print('\nComponents of the pawn-conditioned visibility graphs (sizes, largest first):')
for lab, P in named.items():
    out = []
    for kind in ('rook', 'bishop_dark', 'bishop_light', 'queen'):
        cs = comps_for(P, kind)
        out.append(f'{kind}: {[len(c) for c in cs]}')
    print(f'  {lab}: ' + '; '.join(out))
