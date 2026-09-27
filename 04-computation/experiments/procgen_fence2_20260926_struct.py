"""procgen_fence2_20260926_struct.py -- exact checks of the per-fence structure lemmas (Part A of
lane "fence2", session collatz-procgen-20260922, 2026-09-26) on analysed configurations.

check_structure(R) verifies, for the output R of procgen_fence2_20260926_torus.analyze:
 A1  every non-straight corner has an end ray (o + iota >= 1)            [Sector lemma / P1]
 A2  every reflex corner is R-type with o = iota = 1, at a junction with t = 0
 A3  every straight corner is S0 (o=iota=0, through) or S2 (o=iota=1, straight join)
 A4  per walk: sum(e_i - 1) = #E + #R                                      [P2, generalized]
 A5  side lengths: e=2 -> exactly j+1; e=1 -> in (j, j+1); e=0 -> in (max(0,j-1), j+1)
 A6  #O = #I = Lambda (number of landings = junction sides of a through fence with >= 1 stem)
 A7  a walk without whole side is pure O or pure I, convex (no R), turning +2pi (outer walk)
 A8  every O corner is the extreme sector next to a through ray on a stem side, and the piece of
     that fence side ending there is the incoming side; the next piece starts at the I corner
 A9  pinwheel pairing: O-walk pure O and I-walk pure I at a landing (no straight joins on the two
     sides) => that fence side has exactly one landing and the two sides add up to the fence
"""
from fractions import Fraction as Fr

import procgen_fencelim_20260926_geom as G
from procgen_fencelim_20260926_geom import sgn, sub, cross, dot
from procgen_fence2_20260926_torus import corner_type, walk_sides


def _sqrt_between(len2, lo, hi):
    """lo < sqrt(len2) < hi (lo, hi integers >= 0), exact."""
    ok_lo = (lo <= 0 and sgn(len2) > 0) or sgn(len2 - lo * lo) > 0
    ok_hi = sgn(len2 - hi * hi) < 0
    return ok_lo and ok_hi


def landings(R):
    """list of landings: (vertex key, through fence i, side +1 left / -1 right, #stems)."""
    out = []
    for kk, lst in R['order'].items():
        J = R['junctions'][kk]
        if J['t'] != 1:
            continue
        # find the through fence: a fence i with both rays (i,k,+1) and (i,k-1,-1) at this vertex
        hs = [h for (_, h) in lst]
        thr = None
        for (i, k, s) in hs:
            if s == 1 and (i, k - 1, -1) in hs:
                thr = (i, k)
        assert thr is not None
        i, k = thr
        idx_f = hs.index((i, k, 1))
        idx_b = hs.index((i, k - 1, -1))
        d = len(hs)
        left = (idx_b - idx_f - 1) % d
        right = (idx_f - idx_b - 1) % d
        if left:
            out.append((kk, i, 1, left))
        if right:
            out.append((kk, i, -1, right))
    return out


def check_structure(R, verbose=False):
    walks = R['walks']
    L = R['lattice']
    res = dict(A1=True, A2=True, A3=True, A4=True, A5=True, A6=True, A7=True, A8=True, A9=True)
    cnt = dict(O=0, I=0, E=0, R=0, S0=0, S2=0)
    npin = 0
    nwhole = 0
    walk_info = []
    for w in walks:
        cs = w['corners']
        types = [corner_type(c) for c in cs]
        for c, t in zip(cs, types):
            if t in cnt:
                cnt[t] += 1
            if c['cls'] != 0 and c['o'] + c['iota'] < 1:
                res['A1'] = False
            if c['cls'] == 1:
                if not (c['o'] == 1 and c['iota'] == 1 and R['junctions'][c['vkey']]['t'] == 0):
                    res['A2'] = False
            if c['cls'] == 0 and t not in ('S0', 'S2'):
                res['A3'] = False
        sides = walk_sides(w)
        nE = sum(1 for t in types if t == 'E')
        nR = sum(1 for t in types if t == 'R')
        if sum(s['e'] - 1 for s in sides) != nE + nR:
            res['A4'] = False
        for s in sides:
            e, j, l2 = s['e'], s['j'], s['len2']
            if e == 2:
                ok = sgn(l2 - (j + 1) ** 2) == 0
            elif e == 1:
                ok = _sqrt_between(l2, j, j + 1)
            else:
                ok = _sqrt_between(l2, max(0, j - 1), j + 1)
            if not ok:
                res['A5'] = False
            if e == 2:
                nwhole += 1
        whole = [s for s in sides if s['e'] == 2]
        ns = [t for t in types if t not in ('S0', 'S2')]
        pure = None
        if not whole:
            if set(ns) == {'O'}:
                pure = 'O'
            elif set(ns) == {'I'}:
                pure = 'I'
            else:
                res['A7'] = False
            if sgn(w['A2']) <= 0:
                res['A7'] = False
            npin += 1
        walk_info.append(dict(types=types, sides=sides, pure=pure))
    lands = landings(R)
    Lam = len(lands)
    if not (cnt['O'] == cnt['I'] == Lam):
        res['A6'] = False
    # A8/A9: locate, for each landing, the O corner and I corner and their walks
    loc = {}
    for wi, w in enumerate(walks):
        for ci, c in enumerate(w['corners']):
            loc[(c['h_in'], c['h_out'])] = (wi, ci)
    fp = R['fence_pts']
    n9 = 0
    for (kk, i, side, m) in lands:
        lst = R['order'][kk]
        hs = [h for (_, h) in lst]
        # the landing point is p_k of fence i for the k with (i,k,+1) at kk
        k = [h[1] for h in hs if h[0] == i and h[2] == 1 and (i, h[1] - 1, -1) in hs][0]
        d = len(hs)
        idx_f = hs.index((i, k, 1))
        idx_b = hs.index((i, k - 1, -1))
        if side == 1:
            # left side of i (ccw from forward to backward). Walk along the left side goes forward
            # (p_k-1 -> p_k -> p_k+1): arrives by (i,k-1,+1); O corner: arrive (i,k-1,+1), leave first
            # stem clockwise from backward ray, i.e. hs[idx_b-1]; I corner: arrive rev(hs[idx_f+1]),
            # leave (i,k,+1)
            h_in_O, h_out_O = (i, k - 1, 1), hs[(idx_b - 1) % d]
            st = hs[(idx_f + 1) % d]
            h_in_I, h_out_I = (st[0], st[1], -st[2]), (i, k, 1)
        else:
            h_in_O, h_out_O = (i, k, -1), hs[(idx_f - 1) % d]
            st = hs[(idx_b + 1) % d]
            h_in_I, h_out_I = (st[0], st[1], -st[2]), (i, k - 1, -1)
        if (h_in_O, h_out_O) not in loc or (h_in_I, h_out_I) not in loc:
            res['A8'] = False
            continue
        wO, cO = loc[(h_in_O, h_out_O)]
        wI, cI = loc[(h_in_I, h_out_I)]
        if corner_type(walks[wO]['corners'][cO]) != 'O' or corner_type(walks[wI]['corners'][cI]) != 'I':
            res['A8'] = False
            continue
        # A9
        if walk_info[wO]['pure'] == 'O' and walk_info[wI]['pure'] == 'I':
            sO = [s for s in walk_info[wO]['sides'] if s['b'] == cO][0]
            sI = [s for s in walk_info[wI]['sides'] if s['a'] == cI][0]
            if sO['j'] == 0 and sI['j'] == 0:
                n9 += 1
                # number of landings on this side of fence i
                nl = sum(1 for (kk2, i2, sd2, m2) in lands if i2 == i and sd2 == side)
                fv = sub(R['fences'][i][1], R['fences'][i][0])
                tot = (sO['vec'][0] + sI['vec'][0], sO['vec'][1] + sI['vec'][1])
                same = (sgn(tot[0] - fv[0]) == 0 and sgn(tot[1] - fv[1]) == 0) or \
                       (sgn(tot[0] + fv[0]) == 0 and sgn(tot[1] + fv[1]) == 0)
                if nl != 1 or not same:
                    res['A9'] = False
    return dict(res=res, cnt=cnt, Lambda=Lam, pinwheels=npin, whole_sides=nwhole, pair_checks=n9,
                walk_info=walk_info)
