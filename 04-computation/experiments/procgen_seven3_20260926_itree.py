#!/usr/bin/env python3
"""
procgen_seven3_20260926_itree.py -- itinerary trees (lane seven3, session collatz-procgen-20260922, 2026-09-26).

An ITINERARY TREE is a negation-symmetric sign rule of q n +- 1 (q = 7) written in MH-itinerary coordinates.
Normalized words: the first token is a valuation (the first sign is normalized away by the symmetry x -> -x);
later tokens are (rel, v) with rel '=' (same sign as the previous symbol) or '!' (opposite sign).  A valuation is
an int >= 2 or ('ge', V) (valuation >= V; such a token is terminal).  A tree is a complete prefix code of such
words; each leaf carries a label (True = flip the MH sign on that cylinder).  to_rule() turns it into a bit rule
(a partition of the odd 2-adic integers into classes with signs; each cylinder is one residue class, Lemma MH(ii)).

Local search (search only, results re-evaluated exactly): counterexample-guided -- evaluate the rule (Markov C
engine), take its witness periodic orbit, and try toggling / refining-then-toggling the leaves that contain the odd
points of the orbit; keep the best candidate.
"""
import os
import sys
import random
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_seven3_20260926_lib import (Q, mh, v2, word_class, check_word_class, rho_markov_c, check_partition,
                                          flip_rule)

DEFAULT_V = (2, 3, ('ge', 4))


def children(word, split_v=DEFAULT_V):
    """itinerary children of a leaf word"""
    if not word:
        return [((v,) if not isinstance(v, tuple) else (v,)) for v in split_v]
    last = word[-1]
    lv = last if not isinstance(last, tuple) or last[0] == 'ge' else last[1]
    if isinstance(lv, tuple) and lv[0] == 'ge':
        V = lv[1]
        base = word[:-1]
        if len(word) == 1:
            return [base + (V,), base + (('ge', V + 1),)]
        rel = last[0]
        return [base + ((rel, V),), base + ((rel, ('ge', V + 1)),)]
    out = []
    for rel in ('=', '!'):
        for v in split_v:
            out.append(word + ((rel, v),))
    return out


def tok_v(tok, first):
    return tok if first else tok[1]


def signed_words(word):
    """the two signed MH words of a normalized word (first sign + and -)"""
    res = []
    for s0 in (1, -1):
        s, w = s0, []
        for i, tok in enumerate(word):
            if i == 0:
                v = tok
            else:
                rel, v = tok
                s = s if rel == '=' else -s
            w.append((s, v))
        res.append(w)
    return res


def to_rule(tree, q=Q):
    """bit rule {(c, d): sign} of a tree {word: flip}; the root (empty word) is not allowed as a leaf with a flip"""
    rule = {}
    for word, flip in tree.items():
        if not word:
            # whole odd set: MH rule at depth 2
            for c in (1, 3):
                rule[(c, 2)] = -mh(c, q) if flip else mh(c, q)
            continue
        for w in signed_words(word):
            c, d = word_class(w, q)
            s = mh(c, q)
            rule[(c, d)] = -s if flip else s
    return rule


def mh_tree():
    return {(2,): False, (3,): False, (('ge', 4),): False}


def norm_itin(x_num, x_den, nsym, q=Q):
    """normalized MH itinerary tokens of the rational x (odd numerator and denominator), nsym symbols"""
    x = Fr(x_num, x_den)
    toks, prev = [], None
    for i in range(nsym):
        D = x.denominator
        r = x.numerator * pow(D, -1, 1 << 64) % (1 << 64)
        s = mh(r, q)
        t = q * x + s
        v = 0
        while t.numerator % 2 == 0:
            t /= 2
            v += 1
            if v > 60:
                break
        toks.append(v if prev is None else ('=' if s == prev else '!', v))
        prev = s
        x = t
        if v > 60:
            break
    return toks


def class_map(tree, q=Q):
    """{(c, d): word} for the leaves of the tree (both signed versions)"""
    m = {}
    for word in tree:
        if not word:
            continue
        for w in signed_words(word):
            m[word_class(w, q)] = word
    return m


def leaf_of(tree, x, cmap=None, D=64):
    """the leaf word of the tree whose cylinder contains the odd 2-adic rational x"""
    cmap = cmap or class_map(tree)
    r = x.numerator * pow(x.denominator, -1, 1 << D) % (1 << D)
    for d in range(1, D + 1):
        key = (r % (1 << d), d)
        if key in cmap:
            return cmap[key]
    return None


def evaluate(tree, want_crit=False):
    rule = to_rule(tree)
    if want_crit:
        return rho_markov_c(rule, want_crit=True)
    val, nl, x0, p = rho_markov_c(rule)
    return val, nl, x0, p


def orbit_odd_points(x0, tree, maxlen=10 ** 5):
    """the odd points of the periodic orbit of x0 under the rule of the tree (exact rationals)"""
    rule = to_rule(tree)
    dmax = max(d for (_, d) in rule)
    from procgen_seven3_20260926_lib import rule_sign
    pts, x = [], x0
    for _ in range(maxlen):
        if x.numerator % 2 == 0:
            x = x / 2
        else:
            r = x.numerator * pow(x.denominator, -1, 1 << (dmax + 2)) % (1 << (dmax + 2))
            pts.append(x)
            s = rule_sign(rule, r, dmax)
            x = (Q * x + s) / 2
        if x == x0:
            break
    return pts


def local_search(tree, iters=200, maxdepth_sym=12, verbose=True, seed=1, tabu_len=50, on_best=None):
    random.seed(seed)
    best = dict(tree)
    bval, bnl, bx, bp = evaluate(best)
    cur, cval, cx = dict(best), bval, bx
    tabu = []
    hist = []
    for it in range(iters):
        pts = orbit_odd_points(cx, cur)
        cands = []
        cm = class_map(cur)
        for x in pts:
            lw = leaf_of(cur, x, cm)
            if lw is None:
                continue
            # toggle
            t1 = dict(cur)
            t1[lw] = not t1[lw]
            cands.append(('toggle', lw, t1))
            # refine then toggle the child containing x
            if len(lw) < maxdepth_sym:
                t2 = dict(cur)
                lab = t2.pop(lw)
                for ch in children(lw):
                    t2[ch] = lab
                cw = leaf_of(t2, x)
                if cw is not None:
                    t2[cw] = not t2[cw]
                    cands.append(('refine', cw, t2))
        scored = []
        for kind, w, t in cands:
            key = (kind, w)
            if key in tabu:
                continue
            try:
                v, nl, x0, p, crit = evaluate(t, want_crit=True)
            except Exception as ex:
                continue
            scored.append((v, crit[0], len(t), kind, w, t, x0))
        if not scored:
            break
        scored.sort(key=lambda z: (z[0], z[1], z[2]))
        v, cr, n, kind, w, t, x0 = scored[0]
        cur, cval, cx = t, v, x0
        tabu.append((kind, w))
        tabu = tabu[-tabu_len:]
        hist.append((it, v, n))
        if v < bval or (v == bval and n < len(best)):
            if on_best is not None and v < bval:
                on_best(dict(t), v, it)
            best, bval = dict(t), v
        if verbose and it % 10 == 0:
            print('  it', it, 'cur', v, float(v), 'leaves', n, 'best', bval, flush=True)
    return best, bval, hist
