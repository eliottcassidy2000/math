"""audit_A_symmetry.py -- audit A: Proposition 2.1 checked computationally with python-chess.
sigma = rank reflection + colour swap + side swap + castling-right swap + ep reflection = Board.mirror().
For random playouts from the initial position (castling rights and en-passant squares occur), check that the
legal moves of sigma(s) are exactly the sigma-images of the legal moves of s, and that terminal status agrees.
Then show that the 180-degree rotation with colour swap is NOT a symmetry on states with castling rights."""
import random, chess
random.seed(3)
def mirror_move(m):
    return chess.Move(chess.square_mirror(m.from_square), chess.square_mirror(m.to_square), promotion=m.promotion)
bad = 0; checked = 0; ep_seen = 0; castle_seen = 0
for g in range(400):
    b = chess.Board()
    for ply in range(random.randint(0, 60)):
        moves = list(b.legal_moves)
        if not moves: break
        m = random.choice(moves)
        b.push(m)
        s = b.mirror()
        assert s.mirror() == b, 'not an involution'
        L1 = {mirror_move(x) for x in b.legal_moves}
        L2 = set(s.legal_moves)
        checked += 1
        if b.ep_square is not None and any(b.is_en_passant(x) for x in b.legal_moves): ep_seen += 1
        if any(b.is_castling(x) for x in b.legal_moves): castle_seen += 1
        if L1 != L2 or b.is_checkmate() != s.is_checkmate() or b.is_stalemate() != s.is_stalemate() or b.is_insufficient_material() != s.is_insufficient_material():
            bad += 1
            if bad <= 3: print('MISMATCH', b.fen(), s.fen())
print(f'sigma (Board.mirror) move-set correspondence on {checked} states from random playouts: mismatches {bad}; states with an en-passant capture available {ep_seen}; with a castling move available {castle_seen}')
# rotation with colour swap on the initial position
b = chess.Board()
r = chess.Board(None)
for sq, p in b.piece_map().items():
    r.set_piece_at(chess.square(7 - chess.square_file(sq), 7 - chess.square_rank(sq)), chess.Piece(p.piece_type, not p.color))
r.turn = chess.BLACK
print('rotated+colour-swapped initial position:', r.fen(), '; kings on', chess.square_name(r.king(chess.WHITE)), chess.square_name(r.king(chess.BLACK)), '; equals initial board?', r.board_fen() == b.board_fen())
print('rank-reflected+colour-swapped initial position equals initial board?', b.mirror().board_fen() == b.board_fen())
# the file reflection preserves everything except castling: a king with rights on d1 is not a legal castling state
f = chess.Board(None)
for sq, p in b.piece_map().items():
    f.set_piece_at(chess.square(7 - chess.square_file(sq), chess.square_rank(sq)), p)
f.castling_rights = b.castling_rights
print('file-reflected initial position with the original castling rights: python-chess status', f.status(), '(BAD_CASTLING expected)')
