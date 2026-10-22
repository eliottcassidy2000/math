"""Audit S50, second pass - b2: explicit patches for the configurations that the written proof of THM-4745 (5)
does not reach (found in b1), each checked by building the stated witness path.

Patch M (MID_A): v is the middle vertex of a C_3 block of A, A with >= 2 blocks.  Then A - v is strong
  (THM-4745 (2) case (a)), so Init(A - v) = Term(A - v) = A'.  The local block is {x_(t-1), v, x_(t+1)}.
  * local arcs forward: x_(t+1) is an out-neighbour of v and beats b in B_loc'  ->  (F1) witness P Q
    (P a Hamiltonian path of A' ending at x_(t+1), Q one of B starting at b);
  * local arcs backward: b in B_loc' beats x_(t-1), an in-neighbour of v  ->  (B1) witness Q P.
  MID_B is symmetric with (F2) / (B2).
Patch E (END_A, forward): v = x_0, the first block of A is C_3 = {x_0, x_1, x_2}, A has >= 2 blocks, local arcs
  forward.  Witness  x_1 Q' v W :  Q' a Hamiltonian path of B' starting in B_loc' (beaten by x_1) and ending at an
  in-neighbour of v in B;  W a Hamiltonian path of A - {v, x_1} starting at an out-neighbour of v
  (e.g. x_3 ... x_(l_3) x_2 x_(l_3 + 1) ... x_N).  Its restriction to A starts x_1 v, and v -> x_1.
  END_A backward (v = x_N, last block C_3): the converse picture, witness  W~ v Q~ x_(N-1).
Patch E' (END_B, forward): v = y_M, the last block of B is C_3, B with >= 2 blocks, local arcs forward.
  Witness  Z v Q'' y_(M-1) :  Z a Hamiltonian path of B - {v, y_(M-1)} ending at an in-neighbour of v;  Q'' a
  Hamiltonian path of A' starting at an out-neighbour of v in A and ending in A_loc' (which beats y_(M-1)).
  Its restriction to B ends v y_(M-1), and y_(M-1) -> v.  END_B backward: the converse picture, y_1 R v S.
Every witness uses only arcs of A, of B and LOCAL cross arcs, so it is a Hamiltonian path of every completion.
Checked for every uncovered case with both chains of N <= 8 (parts in {2,3,4})."""
import time
from collections import Counter
import alib as L
import b1_lib as B1

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def hp(adj, verts, starts=None, ends=None):
    """a Hamiltonian path of the digraph adj restricted to verts, first vertex in starts, last in ends (DFS)."""
    verts = list(verts)
    vs = 0
    for x in verts:
        vs |= 1 << x
    res = [None]

    def dfs(path, used):
        if res[0] is not None:
            return
        if used == vs:
            if ends is None or path[-1] in ends:
                res[0] = list(path)
            return
        nxt = adj[path[-1]] & vs & ~used
        while nxt and res[0] is None:
            low = nxt & -nxt
            w = low.bit_length() - 1
            path.append(w)
            dfs(path, used | low)
            path.pop()
            nxt ^= low

    for s in (starts if starts is not None else verts):
        if s in verts:
            dfs([s], 1 << s)
            if res[0] is not None:
                break
    return res[0]


def is_path(adj, seq):
    return all((adj[seq[t]] >> seq[t + 1]) & 1 for t in range(len(seq) - 1))


def violates(adj, path, n, Aset, Bset):
    if sorted(path) != list(range(n)) or not is_path(adj, path):
        return False
    pa = [x for x in path if x in Aset]; pb = [x for x in path if x in Bset]
    return not is_path(adj, pa) or not is_path(adj, pb)


used = Counter(); bad = []
for (pa, pb, i, j, o, types) in B1.uncovered_cases(8):
    base, free, n, Aset, Bset, target, info = B1.glue(pa, pb, i, j, o)
    v = info['v']; mA = info['mA']; mB = info['mB']; locA = info['locA']; locB = info['locB']
    NA, NB = sum(pa), sum(pb)
    Ap = sorted(Aset - {v}); Bp = sorted(Bset - {v})
    inA = {x for x in Ap if (base[x] >> v) & 1}; outA = set(Ap) - inA
    inB = {y for y in Bp if (base[y] >> v) & 1}; outB = set(Bp) - inB
    Aadj = [base[x] & sum(1 << z for z in Aset) for x in range(n)]          # A's own arcs
    Badj = [base[x] & sum(1 << z for z in Bset) for x in range(n)]
    W = None
    if 'MID_A' in types:
        tA = i
        lo, hi = mA[tA - 1], mA[tA + 1]          # block start x_(t-1), block end x_(t+1)
        b = locB[0]
        if o == 'fwd':
            P = hp(Aadj, Ap, ends={hi}); Q = hp(Badj, Bset, starts={b})
            W = (P + Q) if P and Q else None
        else:
            Q = hp(Badj, Bset, ends={b}); P = hp(Aadj, Ap, starts={lo})
            W = (Q + P) if P and Q else None
        kind = 'M'
    elif 'MID_B' in types:
        lo, hi = mB[j - 1], mB[j + 1]
        a = locA[0]
        if o == 'fwd':
            P = hp(Aadj, Aset, ends={a}); Q = hp(Badj, Bp, starts={lo})       # (F2): a -> y_(t-1)
            W = (P + Q) if P and Q else None
        else:
            Q = hp(Badj, Bp, ends={hi}); P = hp(Aadj, Aset, starts={a})       # (B2): y_(t+1) -> a
            W = (Q + P) if P and Q else None
        kind = 'M'
    elif 'END_A' in types and o == 'fwd':          # v = x_0, first block C_3
        x1 = mA[1]
        Qp = hp(base, Bp, starts=set(locB), ends=inB)
        Wp = hp(Aadj, [x for x in Ap if x != x1], starts=outA)
        W = ([x1] + Qp + [v] + Wp) if Qp and Wp else None
        kind = 'E'
    elif 'END_A' in types and o == 'bwd':          # v = x_N, last block C_3: converse picture
        xN1 = mA[NA - 1]
        Qt = hp(base, Bp, starts=outB, ends=set(locB))
        Wt = hp(Aadj, [x for x in Ap if x != xN1], ends=inA)
        W = (Wt + [v] + Qt + [xN1]) if Qt and Wt else None
        kind = 'E'
    elif 'END_B' in types and o == 'fwd':          # v = y_M, last block C_3
        yM1 = mB[NB - 1]
        Z = hp(Badj, [y for y in Bp if y != yM1], ends=inB)
        Qpp = hp(base, Ap, starts=outA, ends=set(locA))
        W = (Z + [v] + Qpp + [yM1]) if Z and Qpp else None
        kind = "E'"
    elif 'END_B' in types and o == 'bwd':          # v = y_0, first block C_3: converse picture
        y1 = mB[1]
        R = hp(base, Ap, starts=set(locA), ends=inA)
        S = hp(Badj, [y for y in Bp if y != y1], starts=outB)
        W = ([y1] + R + [v] + S) if R and S else None
        kind = "E'"
    used[(kind, o)] += 1
    if W is None or not violates(base, W, n, Aset, Bset):
        bad.append((pa, pb, i, j, o, types))
print(f'     patch usage over {sum(used.values())} uncovered cases (N <= 8): {dict(used)}')
check(not bad, 'every uncovered case is closed by its explicit patch witness (a Hamiltonian path of every completion, '
      'built only from arcs of A, B and the local cross arcs, whose restriction to A or B is not a path)'
      + (f'; FAILED {bad[:5]}' if bad else ''))
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
