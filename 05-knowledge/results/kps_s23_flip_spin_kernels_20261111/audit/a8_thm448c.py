"""a8_thm448c.py (AUDIT) - CD_5 versus the THM-448 tower T_31: finer invariants and automorphism counts."""
import sys, os, collections, itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import a8_thm448 as m
out = []
def P(s):
    print(s); out.append(s); sys.stdout.flush()
def vertex_profile(T):
    n = len(T)
    prof = []
    for v in range(n):
        c = collections.Counter()
        for a, b in itertools.combinations([u for u in range(n) if u != v], 2):
            c[len(T[v] & T[a] & T[b])] += 1
        prof.append(tuple(sorted(c.items())))
    return prof
def aut_count(T, prof):
    n = len(T)
    img = [-1] * n; used = [False] * n; cnt = 0
    order = sorted(range(n), key=lambda v: (prof.count(prof[v]), prof[v]))
    def rec(k):
        nonlocal cnt
        if k == n:
            cnt += 1; return
        a = order[k]
        for v in range(n):
            if used[v] or prof[v] != prof[a]: continue
            if all((a in T[order[i]]) == (v in T[img[order[i]]]) and (order[i] in T[a]) == (img[order[i]] in T[v]) for i in range(k)):
                img[a] = v; used[v] = True
                rec(k + 1)
                used[v] = False
    rec(0)
    return cnt
def iso(A, pa, B, pb):
    n = len(A)
    if sorted(pa) != sorted(pb): return False
    img = [-1] * n; used = [False] * n
    order = sorted(range(n), key=lambda v: (pa.count(pa[v]), pa[v]))
    def rec(k):
        if k == n: return True
        a = order[k]
        for v in range(n):
            if used[v] or pb[v] != pa[a]: continue
            if all((a in A[order[i]]) == (v in B[img[order[i]]]) for i in range(k)):
                img[a] = v; used[v] = True
                if rec(k + 1): return True
                used[v] = False
        return False
    return rec(0)
for k in (4, 5):
    Tc = m.cd_tournament(k)
    Tt = m.core(m.tower_matrix(k))
    pc, pt = vertex_profile(Tc), vertex_profile(Tt)
    P("n = %d: classes of vertices by triple-intersection profile: CD_n %s ; tower %s" % (k, sorted(collections.Counter(pc).values()), sorted(collections.Counter(pt).values())))
    P("        multisets of vertex profiles equal: %s" % (sorted(pc) == sorted(pt)))
    P("        |Aut(CD_n)| = %d, |Aut(tower)| = %d" % (aut_count(Tc, pc), aut_count(Tt, pt)))
    Tv = m.converse(Tt); pv = vertex_profile(Tv)
    P("        CD_n isomorphic to the tower: %s ; to its converse: %s" % (iso(Tc, pc, Tt, pt), iso(Tc, pc, Tv, pv)))
open(os.path.join(m.HERE, "a8_thm448c.out"), "w").write("\n".join(out) + "\n")
