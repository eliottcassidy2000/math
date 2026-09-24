#!/usr/bin/env python3
"""Part 2 (G'): how many pairs must be flipped to make the Collatz pairing PROVABLY a tree by an L-step
descent certificate?

Pairing map (offset 0): pair i = {2i-1, 2i}; n moves up by i iff (n odd) XOR eps_i, else down by i.
eps = 0 is the Collatz shortcut T.  Class P_L: every n >= 3 has F^t(n) < n for some t <= L (then every orbit
reaches the root: PROVED by strong induction).  Collatz is in no P_L (n = 2^L - 1 rises L times).

(1) exact window optimum for P_2 (= landing => down) and window designs for P_4, P_5 (CP-SAT, 2 workers,
    time and memory capped; for L >= 4 the solver returns FEASIBLE designs, not proven optima);
(2) a sequential greedy that certifies n = 3, 4, 5, ... in turn, flipping the fewest unfrozen pairs and freezing
    the pairs its certificate uses -- an infinite construction if it never fails (checked to N).
The window problems count flips on pairs <= X and impose the descent condition on n <= 2X only, so (1) gives
exact lower bounds only for L = 2 and indicative values otherwise.  Runtime ~ 4 minutes, memory < 400 MB.
"""
import sys, time
try:
    from ortools.sat.python import cp_model
except Exception:
    cp_model = None

def pair(v):
    return (v + 1) // 2

def paths(n, L):
    """move words from n whose value first drops below n exactly at their last step (length <= L)."""
    out = []
    def rec(v, depth, acc):
        if depth == L:
            return
        i = pair(v)
        for up in (1, 0):
            w = v + i if up else v - i
            acc.append((v, up))
            if w < n:
                out.append(list(acc))
            else:
                rec(w, depth + 1, acc)
            acc.pop()
    rec(n, 0, [])
    return out

def window(X, L, tl):
    m = cp_model.CpModel(); var = {}
    def eps(i):
        if i not in var:
            var[i] = m.NewBoolVar(f"e{i}")
        return var[i]
    m.Add(eps(1) == 0)
    for n in range(3, 2 * X + 1):
        zs = []
        for p in paths(n, L):
            lits = []
            for v, up in p:
                e = eps(pair(v))
                want = (v % 2 == 0) if up else (v % 2 == 1)      # up(v) = odd(v) XOR eps
                lits.append(e if want else e.Not())
            if len(lits) == 1:
                zs.append(lits[0])
            else:
                z = m.NewBoolVar("")
                for l in lits:
                    m.AddImplication(z, l)
                zs.append(z)
        m.AddBoolOr(zs)
    m.Minimize(sum(eps(i) for i in range(1, X + 1)))
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = tl
    s.parameters.max_memory_in_mb = 600
    st = s.Solve(m)
    return s.StatusName(st), s.ObjectiveValue(), s.BestObjectiveBound()

def greedy(L, N):
    P = int(1.6 ** L * N) + 16
    val = bytearray(P + 1); fixed = bytearray(P + 1); fixed[1] = 1
    fails = 0
    for n in range(3, N + 1):
        best = None
        stack = [(n, 0, 0, ())]
        while stack:
            v, d, c, asg = stack.pop()
            if best is not None and c >= best[0]:
                continue
            i = pair(v); a = dict(asg)
            bits = [a[i] if i in a else val[i]] if (fixed[i] or i in a) else [0, 1]
            for b in bits:
                up = (v & 1) ^ b
                w = v + i if up else v - i
                free = not fixed[i] and i not in a
                nc = c + (1 if (b == 1 and free) else 0)
                nasg = asg + ((i, b),) if free else asg
                if w < n:
                    if best is None or nc < best[0]:
                        best = (nc, nasg)
                elif d + 1 < L:
                    stack.append((w, d + 1, nc, nasg))
        if best is None:
            fails += 1
            continue
        for i, b in best[1]:
            val[i] = b; fixed[i] = 1
    half = N // 2
    # independent re-check: under the final eps (unfrozen bits = 0), every certified n descends within L steps
    bad = 0
    for n in range(3, N + 1):
        v, ok = n, False
        for _ in range(L):
            i = pair(v); up = (v & 1) ^ val[i]
            v = v + i if up else v - i
            if v < n:
                ok = True; break
        bad += not ok
    return sum(val[1:half + 1]) / half, fails, bad

print("=" * 100)
print("PART 2 (G'). Flip densities needed for L-step provability (the provable rung of the pairing ladder)")
print("=" * 100)
print("  descent words: first descents happen at L = 1, 2, 4, 5, 7, 8, 10, 12, 13, ... (floor(1 + a log2 3)); other lengths")
print("  (6, 9, 11) occur only for n <= 19 (paths of length <= 13, n <= 20000), so P_6 = P_5 away from those small n")
if cp_model is not None:
    for X, L, tl in [(3000, 2, 60), (3000, 4, 90), (3000, 5, 90)]:
        t0 = time.time()
        st, ob, bd = window(X, L, tl)
        print(f"  CP-SAT window X = {X}, L = {L}: {st}, flips {ob:.0f} (density {ob/X:.4f}), proven lower bound {bd:.0f} ({bd/X:.4f})"
              f"  [{time.time()-t0:.0f}s]", flush=True)
else:
    print("  (ortools not available)")
for L, N in [(2, 20000), (4, 20000), (5, 20000), (7, 20000), (8, 20000), (4, 100000), (5, 100000)]:
    t0 = time.time()
    d, f, bad = greedy(L, N)
    print(f"  sequential greedy L = {L}, n <= {N}: flip density on pairs <= N/2 = {d:.4f}, uncertifiable n: {f},"
          f" re-check: n <= N without descent in L steps under the final eps: {bad}  [{time.time()-t0:.1f}s]", flush=True)
print("  Reading: Collatz violates P_2 on half the pairs (every even i); the cheapest P_2 member flips ~29% (exact window),")
print("  while longer certificates need far fewer flips (window designs 14.5% at L = 4, 5.9% at L = 5).  Collatz itself")
print("  is the L = infinity end of this rung, with flip density 0.")
