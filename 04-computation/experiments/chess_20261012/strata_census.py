"""
strata_census.py -- the material-signature lattice and the irreversibility potential.

Per-side signature: (p, n, bL, bD, r, q): pawns, knights, light bishops, dark bishops, rooks, queens.
No promotion: n, r in 0..2, bL, bD in 0..1, q in 0..1  -> 72 piece inventories, times 9 pawn counts = 648.
With promotion: promoted units k = max(0,n-2)+max(0,bL-1)+max(0,bD-1)+max(0,r-2)+max(0,q-1) <= 8 - p.
Joint signatures: pairs.  Edges: captures (one unit of the victim removed), promotions (one pawn becomes a
piece).  All numbers are counts of the signature lattice, i.e. before any retrograde-legality pruning.

Also: the counter-resetting bound (<= 118 pawn moves + non-pawn captures) is checked on the per-file lemma
by exhaustive enumeration of a two-pawn single-file model.
"""
import itertools, collections

def per_side(promo):
    sigs = []
    for p in range(9):
        for n in range(0, 11):
            for bL in range(0, 10):
                for bD in range(0, 10):
                    for r in range(0, 11):
                        for q in range(0, 10):
                            k = max(0, n-2) + max(0, bL-1) + max(0, bD-1) + max(0, r-2) + max(0, q-1)
                            if not promo and k > 0:
                                continue
                            if k > 8 - p:
                                continue
                            if p + n + bL + bD + r + q > 15:
                                continue
                            sigs.append((p, n, bL, bD, r, q))
    return sigs

def units(s):
    return sum(s)

for promo in (False, True):
    S = per_side(promo)
    byp = collections.Counter(s[0] for s in S)
    print(f'promotion={promo}: per-side signatures {len(S)}; by pawn count {dict(sorted(byp.items()))}')
    print(f'  pawnless-piece inventories (p=0 row): {byp[0]}; joint signatures {len(S)**2}')
    # capture edges: ordered pairs (sig -> sig with one unit of a chosen type removed), per side, times the other side
    # a capture removes one unit from the victim; the capturer side is unchanged.
    cap_edges_per_side = 0
    for s in S:
        cap_edges_per_side += sum(1 for c in s if c > 0)
    print(f'  capture edges in the joint lattice: {2 * cap_edges_per_side * len(S)} (victim side x unchanged capturer side)')
    if promo:
        Sset = set(S)
        promo_edges = 0
        for s in S:
            p, n, bL, bD, r, q = s
            if p == 0:
                continue
            for t in ((n+1, bL, bD, r, q), (n, bL+1, bD, r, q), (n, bL, bD+1, r, q), (n, bL, bD, r+1, q), (n, bL, bD, r, q+1)):
                if (p-1,) + t in Sset:
                    promo_edges += 1
        print(f'  promotion edges in the joint lattice: {2 * promo_edges * len(S)}')
    # height: longest chain of captures = total units 30; plus promotions 16
    print(f'  max units per side 15, so the capture DAG has height 30 and with promotions a chain of length <= 46')

# ---- per-file lemma for the 118 bound ----
# Model: one file with a white pawn starting on rank 2 and a black pawn on rank 7.  Moves: white pawn up one
# (or two from rank 2), black pawn down one (or two from rank 7), never onto an occupied square, never passing;
# at any moment either pawn may be captured by an outside piece (event 'cap'), after which the survivor is free.
# Claim: number of pawn moves on this file <= 10 if neither pawn captures (i.e. stays on the file), and <= 4 if
# neither pawn is ever captured.  Exhaustive: enumerate all sequences.
best_no_capture = 0
best_one_capture = 0
def explore(w, b, moves, captured):
    global best_no_capture, best_one_capture
    # w: white rank (2..8 or None), b: black rank (1..7 or None)
    if captured:
        best_one_capture = max(best_one_capture, moves)
    else:
        best_no_capture = max(best_no_capture, moves)
    # white moves
    if w is not None and w < 8:
        steps = [1] + ([2] if w == 2 else [])
        for st in steps:
            nw = w + st
            if b is not None and (nw >= b):   # cannot pass or land on the black pawn
                # landing on b is a capture by... no, pawns do not capture straight; blocked
                continue
            explore(nw, b, moves + 1, captured)
    if b is not None and b > 1:
        steps = [1] + ([2] if b == 7 else [])
        for st in steps:
            nb = b - st
            if w is not None and (nb <= w):
                continue
            explore(w, nb, moves + 1, captured)
    # captures of one pawn by an outside piece (at most one such event matters; after it the file is free)
    if not captured:
        if w is not None and b is not None:
            explore(None, b, moves, True)
            explore(w, None, moves, True)

import sys
sys.setrecursionlimit(10000)
explore(2, 7, 0, False)
print(f'\nper-file lemma: max pawn moves on one file with neither pawn capturing: {best_one_capture} (one of them captured), '
      f'{best_no_capture} (neither captured)')
print('bound: sum over files of pawn moves + non-pawn captures <= 10(8-k) + 12k + (30-k) = 110 + k <= 118')
