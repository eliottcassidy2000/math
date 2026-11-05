#!/usr/bin/env python3
"""ribbon.py -- bounded-span ("ribbon") tournaments and the profile transfer calculus for
Hamiltonian-path counts.

A ribbon word of span w is a sequence of letters a_1 ... a_n, each letter a subset of {1..w}:
vertex q (1-based) has the back arcs q -> q - s for s in a_q with s < q; every other arc of the
transitive base order points forward.  H(word) is the number of directed Hamiltonian paths.

Two independent counts:
  (i)  brute force: subset DP over (vertex set, last vertex), O(2^n n^2);
  (ii) profile transfer: H = final . M_{a_n} ... M_{a_1} . init with nonnegative integer matrices
       indexed by the profile states below.

Profile state after processing vertices 1..q: the Hamiltonian path restricted to {1..q} is a set
of vertex-disjoint fragments.  A fragment is (head, tail): head = age of its first vertex
(age j = vertex q - j, 0 <= j < w) or -1 if the head has left the window unfed (then it must be
the global start; two such fragments kill the state); tail = age of its last vertex or -1 if it
has left the window (then every arc to a future vertex is forward: "free").
A singleton fragment has head == tail.  The state is the sorted tuple of fragments.
"""
import sys, random, itertools

def tournament_from_word(w, word):
    """adjacency as bitmask rows: adj[p] has bit q set iff p -> q (0-based vertices)."""
    n = len(word)
    adj = [0] * n
    for q in range(n):
        for p in range(q):
            s = q - p
            if s <= w and s in word[q]:
                adj[q] |= 1 << p
            else:
                adj[p] |= 1 << q
    return adj

def ham_count_bruteforce(adj):
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0:
                continue
            out = adj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1
                out &= out - 1
                dp[mask | (1 << u)][u] += c
    return sum(dp[full])

def is_strong(adj):
    n = len(adj)
    if n == 1:
        return True
    # reachability from 0 and to 0
    def reach(adjm):
        seen = 1
        frontier = 1
        while frontier:
            v = (frontier & -frontier).bit_length() - 1
            frontier &= frontier - 1
            new = adjm[v] & ~seen
            seen |= new
            frontier |= new
        return seen
    full = (1 << n) - 1
    if reach(adj) != full:
        return False
    radj = [0] * n
    for p in range(n):
        out = adj[p]
        while out:
            q = (out & -out).bit_length() - 1
            out &= out - 1
            radj[q] |= 1 << p
    return reach(radj) == full

def letters(w):
    return [frozenset(c) for r in range(w + 1) for c in itertools.combinations(range(1, w + 1), r)]

NEW = -10  # sentinel for the endpoint at the new vertex

def step(state, letter, w):
    """all successor states (with multiplicity, as a list) of a profile state under a letter."""
    frags = list(state)
    results = []
    def finish(newfrags):
        out = []
        initial = 0
        for (h, t) in newfrags:
            if h == NEW:
                h2 = 0
            elif h == -1:
                h2 = -1
            else:
                h2 = h + 1
                if h2 >= w:
                    h2 = -1
            if t == NEW:
                t2 = 0
            elif t == -1:
                t2 = -1
            else:
                t2 = t + 1
                if t2 >= w:
                    t2 = -1
            if h2 == -1:
                initial += 1
            out.append((h2, t2))
        if initial >= 2:
            return None
        return tuple(sorted(out))
    # option A: new singleton
    s = finish(frags + [(NEW, NEW)])
    if s is not None:
        results.append(s)
    tail_ok = [(t == -1) or ((t + 1) not in letter) for (h, t) in frags]
    head_ok = [(h != -1) and ((h + 1) in letter) for (h, t) in frags]
    m = len(frags)
    # option B: append after fragment i
    for i in range(m):
        if tail_ok[i]:
            nf = frags[:i] + [(frags[i][0], NEW)] + frags[i + 1:]
            s = finish(nf)
            if s is not None:
                results.append(s)
    # option C: prepend before fragment j
    for j in range(m):
        if head_ok[j]:
            nf = frags[:j] + [(NEW, frags[j][1])] + frags[j + 1:]
            s = finish(nf)
            if s is not None:
                results.append(s)
    # option D: join i -> new -> j
    for i in range(m):
        if not tail_ok[i]:
            continue
        for j in range(m):
            if i == j or not head_ok[j]:
                continue
            nf = [frags[k] for k in range(m) if k != i and k != j] + [(frags[i][0], frags[j][1])]
            s = finish(nf)
            if s is not None:
                results.append(s)
    return results

def profile_count(w, word):
    cur = {(): 1}
    for a in word:
        nxt = {}
        for st, c in cur.items():
            for s2 in step(st, a, w):
                nxt[s2] = nxt.get(s2, 0) + c
        cur = nxt
    return sum(c for st, c in cur.items() if len(st) == 1)

def transfer_system(w):
    """reachable profile states, their index, the matrices M_a (dict (i,j)->mult meaning j -> i),
    the initial index and the final (one-fragment) indices."""
    L = letters(w)
    states = {(): 0}
    order = [()]
    mats = {a: {} for a in L}
    frontier = [()]
    while frontier:
        nf = []
        for st in frontier:
            i = states[st]
            for a in L:
                for s2 in step(st, a, w):
                    if s2 not in states:
                        states[s2] = len(order)
                        order.append(s2)
                        nf.append(s2)
                    j = states[s2]
                    mats[a][(j, i)] = mats[a].get((j, i), 0) + 1
        frontier = nf
    final = [states[st] for st in order if len(st) == 1]
    return order, states, mats, final

def count_with_matrices(w, word, sysw=None):
    if sysw is None:
        sysw = transfer_system(w)
    order, states, mats, final = sysw
    vec = {0: 1}
    for a in word:
        M = mats[a]
        nv = {}
        for (j, i), mult in M.items():
            if i in vec:
                nv[j] = nv.get(j, 0) + mult * vec[i]
        vec = nv
    return sum(vec.get(j, 0) for j in final)

def random_word(w, n, rng):
    return [frozenset(s for s in range(1, w + 1) if rng.random() < 0.5) for _ in range(n)]

def word_str(word, w):
    return ".".join("".join(str(s) if s in a else "-" for s in range(1, w + 1)) for a in word)

if __name__ == "__main__":
    rng = random.Random(20261104)
    for w in (1, 2, 3):
        order, states, mats, final = transfer_system(w)
        print(f"span w={w}: {len(order)} profile states, {len(letters(w))} letters, {len(final)} one-fragment states")
        checks = 0
        for n in range(1, 10):
            for _ in range(40 if n <= 8 else 12):
                word = random_word(w, n, rng)
                adj = tournament_from_word(w, word)
                h1 = ham_count_bruteforce(adj)
                h2 = profile_count(w, word)
                h3 = count_with_matrices(w, word, (order, states, mats, final))
                assert h1 == h2 == h3, (w, n, word_str(word, w), h1, h2, h3)
                assert h1 % 2 == 1
                checks += 1
        print(f"  verified brute force == profile DP == matrix product on {checks} random words (n <= 9)")
    print("ALL CHECKS PASSED")
