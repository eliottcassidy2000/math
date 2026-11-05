#!/usr/bin/env python3
"""witness47.py -- explicit span-3 ribbon words for the strong realisers of H = 47 (order 7), with the
letters, re-verified by the transfer engine and by brute force."""
import sys
from width import gentourng, width
from ribbon import ham_count_bruteforce, is_strong, tournament_from_word, count_with_matrices, transfer_system, word_str
from functools import lru_cache
def ordering_with_span(adj, w):
    n = len(adj); full = (1 << n) - 1
    @lru_cache(maxsize=None)
    def go(S, W):
        if S == full: return ()
        outside = S
        for u in W: outside &= ~(1 << u)
        for v in range(n):
            if S >> v & 1: continue
            o = outside; ok = True
            while o:
                u = (o & -o).bit_length() - 1; o &= o - 1
                if not (adj[u] >> v & 1): ok = False; break
            if not ok: continue
            res = go(S | (1 << v), (W + (v,))[-w:] if w > 0 else ())
            if res is not None: return (v,) + res
        return None
    return go(0, ())
sysw3 = transfer_system(3)
for adj in gentourng(7):
    if ham_count_bruteforce(adj) == 47 and is_strong(adj):
        order = ordering_with_span(adj, 3)
        pos = {v: i for i, v in enumerate(order)}
        word = []
        for q, v in enumerate(order):
            letter = frozenset(q - pos[u] for u in range(7) if adj[v] >> u & 1 and pos[u] < q)
            word.append(letter)
        adj2 = tournament_from_word(3, word)
        print("H = 47 strong order-7 realiser: span-3 word", word_str(word, 3), "| transfer H =", count_with_matrices(3, word, sysw3), "| brute force on the word's tournament =", ham_count_bruteforce(adj2), "| strong:", is_strong(adj2))
