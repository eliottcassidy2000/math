#!/usr/bin/env python3
"""ribbon_packing.py -- the odd-cycle packing polynomial I(Omega(T), x) = sum_k alpha_k(T) x^k of a ribbon word
(back arcs of span <= w on an ordered vertex set, THM-4880) as a finite-state series.

Frontier automaton. Process the vertices in order. A packing of vertex-disjoint odd cycles restricted to the
processed prefix {1..q} is a set of completed cycles (each weighted x) and a set of vertex-disjoint directed path
fragments of unfinished cycles. The head of a fragment (its first vertex) will be entered from a later vertex by a
back arc, so it lies in the window {q-w+1..q}; the tail (its last vertex) will be left by a forward arc to a later
vertex, which exists unless that later vertex has a back arc to the tail, so a tail outside the window is FREE. A
fragment = (head age in [0,w), tail age in [0,w) or FREE, parity of its length). The new vertex q+1 may be skipped,
start a new fragment, be appended after a tail (forward arc), be prepended before a head (back arc, span in the
letter), or both: two distinct fragments merge, or one fragment closes into a cycle, allowed iff its length + 1 is
odd (weight x). Ageing: heads reaching age w kill the state; tails reaching age w become FREE. Initial and final
state: no fragments. Hence I(Omega(T_u), x) = <empty| M_{a_n}(x) ... M_{a_1}(x) |empty>, and x = 2 gives the
Hamiltonian-path count H (Grinberg-Stanley / OCF), a test of this automaton against THM-4880's.

Usage: ribbon_packing.py selftest        (automaton vs brute force on random words; x = 2 vs ribbon.py's H)
       import: packing_system(w) -> (order, states, mats, init_index, final_index); mats[a][(j,i)] = (c0, c1)
       meaning the entry c0 + c1 x of M_a at (j, i) (j -> i is the transition from state i to state j)."""
import sys, os, random, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ribbon_spectrum_20261104'))
from ribbon import letters, tournament_from_word, ham_count_bruteforce, count_with_matrices, transfer_system, word_str

FREE = None

def step(state, letter, w):
    """yield (new_state, closed) for the new vertex with back-arc spans `letter`; closed in {0, 1}."""
    frags = list(state)
    # which fragment has a head/tail at age j
    options = []
    # 0. skip the vertex
    options.append((frags, None, None))
    # 1. new singleton fragment (head = tail = new vertex, parity 1)
    options.append((frags + [('NEW',)], None, None))
    # 2. append after tail s (forward arc s -> new): tail free, or at age j with span j+1 not in letter
    for i, (h, t, p) in enumerate(frags):
        if t is FREE or (t + 1) not in letter:
            options.append((frags, i, None))
    # 3. prepend before head h (back arc new -> h): span h+1 in letter
    for i, (h, t, p) in enumerate(frags):
        if (h + 1) in letter:
            options.append((frags, None, i))
    # 4. both
    for i, (h, t, p) in enumerate(frags):
        if not (t is FREE or (t + 1) not in letter): continue
        for k, (h2, t2, p2) in enumerate(frags):
            if (h2 + 1) in letter:
                options.append((frags, i, k))
    for fr, ai, pi in options:
        fr = [f for f in fr]
        closed = 0
        if ai is None and pi is None:
            if fr and fr[-1] == ('NEW',):
                fr = fr[:-1] + [(-1, -1, 1)]   # head age -1 (becomes 0 after ageing), tail -1, parity 1
            new = fr
        elif ai is not None and pi is None:
            h, t, p = fr[ai]
            new = fr[:ai] + fr[ai + 1:] + [(h, -1, 1 - p)]
        elif ai is None and pi is not None:
            h, t, p = fr[pi]
            new = fr[:pi] + fr[pi + 1:] + [(-1, t, 1 - p)]
        else:
            if ai == pi:
                h, t, p = fr[ai]
                if (p + 1) % 2 == 1:   # odd cycle closes
                    new = fr[:ai] + fr[ai + 1:]
                    closed = 1
                else:
                    continue
            else:
                h1, t1, p1 = fr[ai]   # tail used
                h2, t2, p2 = fr[pi]   # head used
                merged = (h2, t1, (p1 + p2 + 1) % 2)   # path: fragment pi ... then new vertex ... no: new vertex goes after tail t1 and before head h2:  frag(ai) -> new -> frag(pi); head = h1, tail = t2
                merged = (h1, t2, (p1 + p2 + 1) % 2)
                new = [f for idx, f in enumerate(fr) if idx not in (ai, pi)] + [merged]
        # ageing
        aged = []
        dead = False
        for h, t, p in new:
            h2 = h + 1
            if h2 >= w: dead = True; break
            t2 = FREE if t is FREE else (t + 1 if t + 1 < w else FREE)
            aged.append((h2, t2, p))
        if dead: continue
        yield (tuple(sorted(aged, key=lambda f: (f[0], -1 if f[1] is FREE else f[1], f[2]))), closed)

def packing_system(w):
    L = letters(w)
    states = {(): 0}; order = [()]
    mats = {a: {} for a in L}
    frontier = [()]
    while frontier:
        nf = []
        for st in frontier:
            i = states[st]
            for a in L:
                for s2, closed in step(st, a, w):
                    if s2 not in states:
                        states[s2] = len(order); order.append(s2); nf.append(s2)
                    j = states[s2]
                    c0, c1 = mats[a].get((j, i), (0, 0))
                    mats[a][(j, i)] = (c0 + (closed == 0), c1 + (closed == 1))
        frontier = nf
    return order, states, mats, 0, 0

def poly_with_matrices(w, word, sysw=None):
    """I(Omega(T_word), x) as a list of coefficients, by the automaton."""
    if sysw is None: sysw = packing_system(w)
    order, states, mats, init, final = sysw
    d = len(order)
    vec = {init: [1]}
    for a in word:
        nv = {}
        for (j, i), (c0, c1) in mats[a].items():
            if i in vec:
                v = vec[i]
                acc = nv.setdefault(j, [])
                # add c0 * v and c1 * x * v
                need = len(v) + 1
                while len(acc) < need: acc.append(0)
                for k, c in enumerate(v):
                    if c0: acc[k] += c0 * c
                    if c1: acc[k + 1] += c1 * c
        vec = nv
    out = vec.get(final, [0])
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def value_with_matrices(w, word, x, sysw=None):
    if sysw is None: sysw = packing_system(w)
    order, states, mats, init, final = sysw
    vec = {init: 1}
    for a in word:
        nv = {}
        for (j, i), (c0, c1) in mats[a].items():
            if i in vec: nv[j] = nv.get(j, 0) + (c0 + c1 * x) * vec[i]
        vec = nv
    return vec.get(final, 0)

def packing_poly_bruteforce(adj):
    """I(Omega(T), x) by the set-partition DP: cyc[S] = directed Hamiltonian cycles of T[S] (|S| odd >= 3)."""
    n = len(adj); full = 1 << n
    cyc = [0] * full
    for S in range(full):
        pc = bin(S).count('1')
        if pc < 3 or pc % 2 == 0: continue
        m = (S & -S).bit_length() - 1
        # paths from m inside S: dp[mask][v]
        dp = {(1 << m, m): 1}
        # iterate masks in increasing order
        masks = [mask for mask in range(full) if (mask & S) == mask and (mask >> m) & 1]
        masks.sort(key=lambda z: bin(z).count('1'))
        table = {}
        table[(1 << m, m)] = 1
        for mask in masks:
            for v in range(n):
                c = table.get((mask, v), 0)
                if not c: continue
                out = adj[v] & S & ~mask
                while out:
                    u = (out & -out).bit_length() - 1; out &= out - 1
                    table[(mask | (1 << u), u)] = table.get((mask | (1 << u), u), 0) + c
        total = 0
        for v in range(n):
            if (S >> v) & 1 and v != m and (adj[v] >> m) & 1:
                total += table.get((S, v), 0)
        cyc[S] = total
    p = {0: [1]}
    for S in range(1, full):
        m = (S & -S).bit_length() - 1
        acc = [0]
        sub = S
        while sub:
            if (sub >> m) & 1 and cyc[sub]:
                rest = p[S & ~sub]
                need = len(rest) + 1
                while len(acc) < need: acc.append(0)
                for k, c in enumerate(rest):
                    acc[k + 1] += cyc[sub] * c
            sub = (sub - 1) & S
        p[S] = acc
    total = [0]
    for S in range(full):
        v = p[S]
        while len(total) < len(v): total.append(0)
        for k, c in enumerate(v): total[k] += c
    while len(total) > 1 and total[-1] == 0: total.pop()
    return total

def evalpoly(c, x): return sum(ci * x**i for i, ci in enumerate(c))

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'selftest':
        rng = random.Random(7)
        for w in (1, 2, 3):
            sysw = packing_system(w)
            order = sysw[0]
            L = letters(w)
            path_sys = transfer_system(w)
            bad = 0; tested = 0
            for trial in range(120):
                n = rng.randint(1, 9)
                word = [rng.choice(L) for _ in range(n)]
                adj = tournament_from_word(w, word)
                bf = packing_poly_bruteforce(adj)
                au = poly_with_matrices(w, word, sysw)
                H = ham_count_bruteforce(adj)
                Hpath = count_with_matrices(w, word, path_sys)
                ok = (bf == au) and (evalpoly(au, 2) == H == Hpath)
                tested += 1
                if not ok:
                    bad += 1
                    if bad <= 3: print("MISMATCH w", w, word_str(word, w), "bf", bf, "aut", au, "H", H, Hpath)
            print(f"span {w}: {len(order)} reachable packing states; {tested} random words (n <= 9): {tested - bad} agree with brute force and with H at x = 2")
