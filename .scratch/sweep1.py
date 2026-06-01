from lrc_coloring_core import *
import collections, itertools, random

def analyze(speeds):
    n = len(speeds) + 1
    crit = critical_times(speeds, n)
    mids = cell_midpoints(crit)
    chidist = collections.Counter(); omegadist = collections.Counter()
    diff_cells = []
    Ndist = collections.Counter()
    lonely_cells = 0
    for t in mids:
        k,_ = chromatic_number(speeds, n, t)
        w = clique_number_arc(speeds, n, t)
        chidist[k]+=1; omegadist[w]+=1
        if k != w:
            diff_cells.append((t,k,w))
        Nt = observer_danger_count(speeds, n, t)
        Ndist[Nt]+=1
        if Nt == 0: lonely_cells += 1
    return dict(n=n, cells=len(mids), chi=dict(sorted(chidist.items())),
               omega=dict(sorted(omegadist.items())),
               maxchi=max(chidist), minchi=min(chidist),
               diff=diff_cells, Ndist=dict(sorted(Ndist.items())),
               lonely_cells=lonely_cells)

# Curated speed sets across m=4,5,6
sets = [
    (1,2,3,4),
    (1,2,3,5),
    (1,3,4,7),
    (2,3,5,7),
    (1,2,4,8),
    (1,2,3,4,5),
    (1,2,3,4,6),
    (1,3,4,5,9),
    (2,3,4,5,6),
    (1,2,3,4,5,6),
    (1,2,3,4,5,7),
    (1,3,5,7,9,11),
    (1,2,4,8,16,32),
]
for s in sets:
    r = analyze(s)
    print("speeds", s, "n=",r['n'],"cells=",r['cells'])
    print("   chi  :", r['chi'], " max=",r['maxchi']," min=",r['minchi'])
    print("   omega:", r['omega'])
    print("   chi!=omega cells:", len(r['diff']), r['diff'][:3])
    print("   N(t) dist:", r['Ndist'], " lonely(N=0) cells:", r['lonely_cells'])
