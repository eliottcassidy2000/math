#!/usr/bin/env python3
"""probe47.py -- (1) is 47 attained by any strong span-2 ribbon word (all letters, n <= 9)?
(2) per-length value sets of the 3-letter strong alphabet: min, max, solid interval, misses, tiling.
(3) numerical check of the end-twin lemma on random tournaments."""
import itertools, random
from ribbon import transfer_system, letters, count_with_matrices, tournament_from_word, ham_count_bruteforce, is_strong, word_str
w = 2
sysw = transfer_system(w)
F = frozenset
L = letters(w)
# (1)
found = {}
for n in range(2, 10):
    for word in itertools.product(L, repeat=n):
        if not word[1] <= F([1]):   # second letter cannot have span 2 (no vertex 0); normalise
            continue
        if word[0]:
            continue
        h = count_with_matrices(w, list(word), sysw)
        if h in (47, 7, 21, 67, 11, 13):
            adj = tournament_from_word(w, list(word))
            if is_strong(adj):
                found.setdefault(h, []).append(word_str(list(word), w))
print("(1) strong words (all letters, n<=9) with H in {7,11,13,21,47,67}:", {h: (len(v), v[:2]) for h, v in found.items()})
# (2)
alph = [F([1]), F([2]), F([1, 2])]
prev_d = None
for n in range(3, 13):
    vals = set()
    for word in itertools.product(alph, repeat=n):
        vals.add(count_with_matrices(w, list(word), sysw))
    mn, mx = min(vals), max(vals)
    miss = [h for h in range(mn, mx + 1, 2) if h not in vals]
    # largest solid odd interval
    best = (0, 0, 0)
    start = None
    for h in range(mn, mx + 3, 2):
        if h in vals and start is None:
            start = h
        if h not in vals and start is not None:
            if (h - 2 - start) > best[0]:
                best = (h - 2 - start, start, h - 2)
            start = None
    print(f"n={n}: {len(vals)} values in [{mn},{mx}], misses inside {len(miss)} (first {miss[:12]}), largest solid interval [{best[1]},{best[2]}] (length {best[0]//2+1})", "overlap" if prev_d is not None and best[1] <= prev_d else "")
    prev_d = best[2]
# (3) end-twin lemma
rng = random.Random(7)
bad = 0; tests = 0
for trial in range(300):
    n = rng.randint(3, 8)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
    # make vertex b = n-1 a sole exit to a = 0: b -> 0, all others -> b
    b = n - 1; a = 0
    for i in range(n - 1):
        adj[i] &= ~(1 << b); adj[b] &= ~(1 << i)
    adj[b] |= 1 << a
    for i in range(1, n - 1):
        adj[i] |= 1 << b
    # add c = n: version 1: c -> a, b -> c, others -> c ; version 2: c -> a, c -> b, others -> c
    def with_c(cb):
        A = adj[:] + [0]
        c = n
        A[c] |= 1 << a
        for i in range(1, n - 1):
            A[i] |= 1 << c
        if cb:
            A[c] |= 1 << b
        else:
            A[b] |= 1 << c
        return A
    h1 = ham_count_bruteforce(with_c(False)); h2 = ham_count_bruteforce(with_c(True))
    tests += 1
    if h1 != h2: bad += 1
print(f"(3) end-twin lemma: {tests} random tests, {bad} failures")
