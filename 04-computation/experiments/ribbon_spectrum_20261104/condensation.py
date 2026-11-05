#!/usr/bin/env python3
"""condensation.py -- a finite automaton deciding strong connectivity of span-w ribbon words.
State: the strong components of the prefix, in condensation order, from the first component that
contains a window vertex onward; a component is (frozenset of window ages, has_old).  A component
lying entirely before the first window-containing component can never be merged again, so the word
can never become strong: such prefixes are DEAD.  Verified against brute-force reachability."""
import itertools, random
from ribbon import letters, tournament_from_word, is_strong, word_str

DEAD = None

def cond_step(state, a, w):
    """state: tuple of (frozenset ages, has_old) or DEAD; a: letter (set of spans). Returns new state."""
    if state is DEAD:
        return DEAD
    comps = list(state)
    if not comps:                           # first vertex
        new = [(frozenset([0]), False)]
        return tuple(new)
    allages = frozenset().union(*[ages for ages, old in comps])
    A = frozenset(s - 1 for s in a) & allages   # current ages of the new vertex's out-neighbours (existing vertices only)
    if A:
        i0 = min(i for i, (ages, old) in enumerate(comps) if ages & A)
        m = len(comps)
        # suffix of components entirely inside A
        j1 = m - 1
        while j1 >= i0 and (not comps[j1][1]) and comps[j1][0] <= A:
            j1 -= 1
        # components i0..j1 merge with the new vertex; if j1 < i0 then only the new vertex... but then
        # comps[i0] is inside A, contradiction with i0's definition unless all are inside A: then the
        # new vertex forms its own component placed before them.
        merged_ages = set(); merged_old = False
        for i in range(i0, j1 + 1):
            merged_ages |= comps[i][0]; merged_old = merged_old or comps[i][1]
        newcomp = (frozenset(merged_ages) | frozenset(['new']), merged_old)
        comps = comps[:i0] + [newcomp] + comps[j1 + 1:]
    else:
        comps = comps + [(frozenset(['new']), False)]
    # age shift
    out = []
    for ages, old in comps:
        na = set(); no = old
        for x in ages:
            if x == 'new':
                na.add(0)
            elif x + 1 < w:
                na.add(x + 1)
            else:
                no = True
        out.append((frozenset(na), no))
    # drop leading components without window ages: if any exists and it is not the whole thing -> dead
    k = 0
    while k < len(out) and not out[k][0]:
        k += 1
    if k > 0 and k < len(out):
        return DEAD
    if k == len(out):
        # every component has no window vertex: impossible (the new vertex has age 0)
        raise RuntimeError("no window vertex")
    return tuple(out)

def is_strong_word(w, word):
    st = ()
    for a in word:
        st = cond_step(st, a, w)
        if st is DEAD:
            return False
    return len(st) == 1

def automaton(w):
    """enumerate reachable condensation states; transitions; accepting = single component."""
    L = letters(w)
    start = ()
    states = {start: 0}; order = [start]; trans = {}
    frontier = [start]
    while frontier:
        nf = []
        for st in frontier:
            for a in L:
                s2 = cond_step(st, a, w)
                if s2 is DEAD:
                    trans[(states[st], a)] = -1
                    continue
                if s2 not in states:
                    states[s2] = len(order); order.append(s2); nf.append(s2)
                trans[(states[st], a)] = states[s2]
        frontier = nf
    accept = [states[st] for st in order if len(st) == 1]
    return order, states, trans, accept

if __name__ == "__main__":
    rng = random.Random(3)
    for w in (1, 2, 3):
        L = letters(w)
        bad = 0; tests = 0
        for n in range(1, 8):
            words = itertools.product(L, repeat=n) if len(L) ** n <= 20000 else (tuple(rng.choice(L) for _ in range(n)) for _ in range(20000))
            for word in words:
                word = list(word)
                adj = tournament_from_word(w, word)
                s1 = is_strong(adj); s2 = is_strong_word(w, word)
                tests += 1
                if s1 != s2:
                    bad += 1
                    if bad <= 3: print("MISMATCH", w, word_str(word, w), s1, s2)
        order, states, trans, accept = automaton(w)
        print(f"w={w}: condensation automaton has {len(order)} states ({len(accept)} accepting); brute-force agreement on {tests} words: {bad} mismatches")
