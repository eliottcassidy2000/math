"""audit_A_strata.py -- audit A: KPvK stratum census from MY independent table (audit_A_solver.py), no chesslib.
Usage: python3 audit_A_strata.py <my_tables.pkl>"""
import sys, pickle, collections
import chess
with open(sys.argv[1], 'rb') as f:
    MY = pickle.load(f)
ent = MY['KPvK']; states, val = ent['states'], ent['val']
index = {st: i for i, st in enumerate(states)}
size = collections.Counter(); wtm = collections.Counter(); btm = collections.Counter()
for (sq, turn) in states:
    size[sq[1]] += 1
    (wtm if turn == chess.WHITE else btm)[sq[1]] += 1
sizes = sorted(size.values())
print('strata (pawn squares):', len(size), '; stratum sizes min/max:', sizes[0], sizes[-1], '; distinct sizes:', sorted(set(sizes)))
print('a2 stratum: white-to-move', wtm[chess.A2], 'black-to-move', btm[chess.A2], 'total', size[chess.A2])
king_moves = pawn_pushes = 0; nonadvancing = 0
pieces = ['K', 'P', 'k']
for (sq, turn) in states:
    b = chess.Board(None)
    for c, s in zip(pieces, sq):
        b.set_piece_at(s, chess.Piece.from_symbol(c))
    b.turn = turn
    for m in b.legal_moves:
        if b.is_capture(m) or m.promotion:
            continue
        if b.piece_type_at(m.from_square) == chess.PAWN:
            pawn_pushes += 1
            if chess.square_rank(m.to_square) <= chess.square_rank(m.from_square):
                nonadvancing += 1
        else:
            king_moves += 1
print('internal moves:', king_moves + pawn_pushes, '; king moves (stay in stratum):', king_moves, '; pawn pushes (leave):', pawn_pushes, '; non-advancing pushes:', nonadvancing)
