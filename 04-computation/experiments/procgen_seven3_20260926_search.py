#!/usr/bin/env python3
"""
procgen_seven3_20260926_search.py -- systematic searches over itinerary-coded sign rules of 7n +- 1
(lane seven3, session collatz-procgen-20260922, 2026-09-26).  Exploratory: every rule reported is re-evaluated
exactly by the runner (Markov C engine + Python Markov engine + uniform seven2 engine where the depth allows).

(ii)  PATTERN SETS: complete normalized itinerary trees of symbol length L over an alphabet; every subset of the
      leaves is a flip set (negation-symmetric by construction).  Exhaustive for small L.
(i)   AUTOMATA: a start state reads the first valuation class; transient states read (relative sign, valuation
      class) letters; targets are transient states or FLIP / KEEP; valuation classes >= V are terminal.  The rule is
      truncated at bit depth D (classes still undecided get a default label).  Exhaustive for tiny sizes, random
      sampling above.
(iii) COUNTEREXAMPLE-GUIDED GROWTH: procgen_seven3_20260926_itree.local_search (toggle / refine the leaves on the
      densest cycle), optionally restricted to Lemma F-gaining patterns first.
"""
import os
import sys
import json
import time
import itertools
import random
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_seven3_20260926_lib import Q, mh, word_class, rho_markov_c
from procgen_seven3_20260926_itree import to_rule, signed_words


# ------------------------------------------------------------------------------------------------ (ii) pattern sets
def complete_tree(L, first=(2, ('ge', 3)), later=(2, ('ge', 3))):
    """leaves of the complete normalized tree of symbol length L: first symbol from `first`, later symbols
    (rel, v) with v from `later`; ('ge', V) tokens are terminal"""
    leaves = []

    def rec(word, n):
        if n == L:
            leaves.append(word)
            return
        opts = first if not word else [(r, v) for r in ('=', '!') for v in later]
        for tok in opts:
            v = tok if not word else tok[1]
            w = word + (tok,)
            if isinstance(v, tuple):
                leaves.append(w)
            else:
                rec(w, n + 1)
    rec((), 0)
    return leaves


def eval_flipset(leaves, mask):
    tree = {w: bool((mask >> i) & 1) for i, w in enumerate(leaves)}
    val, nl, x0, p, crit = rho_markov_c(to_rule(tree), check=False, want_crit=True)
    return val, crit[0], tree


def exhaustive_patterns(L, first, later, log=None):
    leaves = complete_tree(L, first, later)
    n = len(leaves)
    best = []
    t0 = time.time()
    for mask in range(1 << n):
        val, crit, tree = eval_flipset(leaves, mask)
        best.append((val, crit, mask))
    best.sort()
    return leaves, best[:20], time.time() - t0


# ------------------------------------------------------------------------------------------------ (i) automata
def automaton_rule(A, D, default=False, vclasses=(2, ('ge', 3)), q=Q):
    """A = (start, trans): start[v-index] -> target, trans[state][(rel, v-index)] -> target; targets are ints (state)
    or 'F' / 'K'.  Terminal classes ('ge', V) must map to 'F'/'K'.  Enumerate words breadth first until the class
    depth exceeds D; undecided words get `default`.  Returns the tree {word: flip}."""
    start, trans = A
    tree = {}
    stack = []
    for i, v in enumerate(vclasses):
        stack.append(((v,), start[i]))
    while stack:
        word, tgt = stack.pop()
        c, d = word_class(signed_words(word)[0], q)
        if tgt in ('F', 'K'):
            tree[word] = (tgt == 'F')
            continue
        if d >= D:
            tree[word] = default
            continue
        for rel in ('=', '!'):
            for j, v in enumerate(vclasses):
                t2 = trans[tgt][(rel, j)]
                stack.append((word + ((rel, v),), t2))
    return tree


def enum_automata(m, vclasses=(2, ('ge', 3))):
    """all automata with m transient states (plus the start state) over the given valuation classes; the last class is
    terminal"""
    nv = len(vclasses)
    targets = list(range(m)) + ['F', 'K']
    term = ['F', 'K']
    start_opts = [targets if not isinstance(vclasses[i], tuple) else term for i in range(nv)]
    letters = [(r, j) for r in ('=', '!') for j in range(nv)]
    state_opts = [targets if not isinstance(vclasses[j], tuple) else term for (r, j) in letters]
    for st in itertools.product(*start_opts):
        for rows in itertools.product(itertools.product(*state_opts), repeat=m):
            trans = [dict(zip(letters, row)) for row in rows]
            yield (list(st), trans)


def random_automaton(m, rng, vclasses=(2, ('ge', 3))):
    nv = len(vclasses)
    targets = list(range(m)) + ['F', 'K']
    term = ['F', 'K']
    st = [rng.choice(targets if not isinstance(vclasses[i], tuple) else term) for i in range(nv)]
    letters = [(r, j) for r in ('=', '!') for j in range(nv)]
    trans = [{l: rng.choice(targets if not isinstance(vclasses[l[1]], tuple) else term) for l in letters}
             for _ in range(m)]
    return (st, trans)


def eval_automaton(A, D, default, vclasses=(2, ('ge', 3))):
    tree = automaton_rule(A, D, default, vclasses)
    if not any(tree.values()):
        return Fr(1, 2), 0, tree
    val, nl, x0, p, crit = rho_markov_c(to_rule(tree), check=False, want_crit=True)
    return val, crit[0], tree


def tree_to_json(tree):
    return [[list(map(lambda t: t if not isinstance(t, tuple) else list(t), w)), f] for w, f in tree.items()]


def tree_from_json(js):
    def tok(t, first):
        if first:
            return t if not isinstance(t, list) else tuple(t)
        rel, v = t
        return (rel, v if not isinstance(v, list) else tuple(v))
    out = {}
    for w, f in js:
        out[tuple(tok(t, i == 0) for i, t in enumerate(w))] = f
    return out
