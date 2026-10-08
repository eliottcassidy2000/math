#!/usr/bin/env python3
"""Audit E, task C: THM-4606 Proof 3 (growth lemma) and the small-start escape paths.

 (C1) Growth path from (0, e), e integer != 0, run in exact rationals with the table, for every odd p in [5, 63] and
      1 <= |e| <= 5000: closed form e* = (p/2)^(j+1) e_1 + A_j with A_j = ((p-1)/2)((p/2)^j - 1)/(p-2) + 1/2;
      0 < A_j < (p/2)^(j+1)/(p-2); the two displayed sign-wise bounds; |e*| >= (p/4)|e| - 1; e* a nonzero integer;
      no visit to (0,0); and the run length at k = -1 equals v_2((p-2) p e_1 + p - 1) (finite: (p-1)/(p-2) is not an integer).
 (C2) The 'stays at k = -1' alternative: it occurs iff p e_1 = -(p-1)/(p-2), i.e. from (0, -1/(p-2)); then |f| is
      CONSTANT (not growing). Shown for p = 5, 7, 9 (non-integer starts only).
 (C3) Own BFS (integer coordinates, exact) from every (0, e), 1 <= |e| <= 4, for p = 5 (x* = 4) and p = 7 (x* = 4/3):
      shortest coin words reaching some (0, e') with |e'| > x*, never visiting (0,0). Also confirms x* = 4/(p-4) is the
      fixed point of x -> (p/4)x - 1 and that for p >= 9 every |e| >= 1 already exceeds x*.
 (C4) Sheet-blindness side check: escape paths from non-integer starts (0, a/r) (translations of px+r not divisible by r)."""
from fractions import Fraction as Fr

def par(e):
    return e.numerator % 2      # odd denominators only

def step(p, k, e, b):
    s = par(e); P = Fr(p)
    if s == 0: return (k, e / 2) if b == 0 else (k, (p * e + 1 - P ** k) / 2)
    return (k + 1, (p * e + 1) / 2) if b == 0 else (k - 1, (e - P ** (k - 1)) / 2)

def v2(n):
    n = abs(n); c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c

def growth(p, e0, maxsteps=100000):
    k, e = 0, Fr(e0); steps = 0; visited00 = False
    halv = 0
    while par(e) == 0:                 # run at k = 0, beta = 1: e -> p e / 2
        k, e = step(p, k, e, 1); steps += 1; halv += 1
        assert k == 0 and e == Fr(p, 2) ** halv * e0
    ea = e
    k, e = step(p, k, e, 1); steps += 1    # departure down
    assert k == -1 and e == (ea - Fr(1, p)) / 2
    e1 = e
    j = 0
    while par(e) == 0:                 # run at k = -1 with beta = 1
        k, e = step(p, k, e, 1); steps += 1; j += 1
        assert k == -1
        if steps > maxsteps: return None
    k, e = step(p, k, e, 0); steps += 1    # return up
    assert k == 0
    return dict(ea=ea, e1=e1, j=j, estar=e, halv=halv, steps=steps)

def c1():
    bad = 0; n = 0; worst_ratio = None
    for p in range(5, 65, 2):
        P = Fr(p)
        for e0 in list(range(1, 5001)) + list(range(-5000, 0)):
            g = growth(p, e0); n += 1
            if g is None: bad += 1; print("RUN NEVER ENDED", p, e0); continue
            j, e1, es = g['j'], g['e1'], g['estar']
            Aj = (P - 1) / 2 * ((P / 2) ** j - 1) / (P - 2) + Fr(1, 2)
            ok = es == (P / 2) ** (j + 1) * e1 + Aj
            ok &= 0 < Aj < (P / 2) ** (j + 1) / (P - 2)
            ok &= es.denominator == 1 and es != 0
            ok &= abs(es) >= P / 4 * abs(e0) - 1
            if e0 > 0: ok &= es >= P * e0 / 4 - Fr(1, 4)
            else: ok &= abs(es) >= P / 2 * (abs(Fr(e0)) / 2 + 1 / (2 * P) - 1 / (P - 2))
            # run length at k = -1 equals v2((p-2) E1 + p - 1), E1 = p e1 (an integer)
            E1 = p * e1
            ok &= E1.denominator == 1 and j == v2((p - 2) * E1.numerator + p - 1)
            if not ok:
                bad += 1
                if bad < 10: print("C1 FAIL", p, e0, g)
    print(f"(C1) growth lemma, odd p in [5,63], 1<=|e|<=5000 ({n} starts): failures {bad}"
          f" (closed form, A_j bounds, sign-wise bounds, |e*| >= (p/4)|e|-1, e* nonzero integer, run length = v2((p-2)pe1+p-1))")
    return bad == 0

def c2():
    for p in (5, 7, 9):
        P = Fr(p)
        e0 = Fr(-1, p - 2)
        k, e = step(p, 0, e0, 1)
        traj = [(k, e)]
        for _ in range(12):
            assert par(e) == 0
            k, e = step(p, k, e, 1); traj.append((k, e))
        print(f"(C2) p={p}: from (0, {e0}) the growth path departs to (-1, {traj[0][1]}) and then stays there: "
              f"{'constant' if all(t == traj[0] for t in traj) else 'not constant'} for 12 run steps (|f| = {float(abs(traj[0][1])):.4f}, not growing)")

def bfs_escape(p, e0, thresh, depth=30):
    """shortest coin word from (0,e0) to (0,e') with |e'| > thresh, never visiting (0,0); exact; returns word, e'."""
    frontier = {(0, Fr(e0)): ''}
    seen = set(frontier)
    for d in range(depth):
        new = {}
        for (k, e), w in frontier.items():
            for b in (0, 1):
                k2, e2 = step(p, k, e, b)
                if (k2, e2) == (0, 0): continue
                if k2 == 0 and abs(e2) > thresh: return w + str(b), e2
                if (k2, e2) not in seen:
                    seen.add((k2, e2)); new[(k2, e2)] = w + str(b)
        frontier = new
    return None, None

def c3():
    ok = True
    for p in range(5, 33, 2):
        xs = Fr(4, p - 4)
        assert Fr(p, 4) * xs - 1 == xs
        if xs < 1:
            continue
        for e0 in (1, 2, 3, 4, -1, -2, -3, -4):
            if abs(e0) > xs:
                print(f"(C3) p={p} start (0,{e0:2d}): already above x*={xs}")
                continue
            w, e2 = bfs_escape(p, e0, xs)
            print(f"(C3) p={p} start (0,{e0:2d}): {'reaches (0,%s) via coins %s (length %d)' % (e2, w, len(w)) if w else 'NO PATH'}")
            ok &= w is not None
    print(f"(C3) x* = 4/(p-4) >= 1 only for p = 5, 7; escape paths found for all small starts: {'PASS' if ok else 'FAIL'}")
    return ok

def c4():
    # translations of px + r by integers a not divisible by r correspond to px+1 translations a/r
    for p, r in ((5, 3), (5, -1), (7, 3), (5, 7)):
        xs = Fr(4, p - 4)
        res = []
        for a in (1, -1, 2, -2):
            e0 = Fr(a, r)
            if e0.denominator == 1: continue
            w, e2 = bfs_escape(p, e0, max(xs, Fr(2)), depth=26)
            res.append(f"e={e0}: {'to (0,%s) via %s' % (e2, w) if w else 'none'}")
        if res: print(f"(C4) px+r with p={p}, r={r}: " + "; ".join(res))

if __name__ == '__main__':
    ok = c1()
    c2()
    ok &= c3()
    c4()
    print("ALL C CHECKS PASS" if ok else "SOME C CHECK FAILED")
