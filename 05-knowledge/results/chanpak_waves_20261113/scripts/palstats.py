import sys
from collections import defaultdict
def longest_run(vals):
    vs=sorted(set(vals)); best=1; cur=1
    for i in range(1,len(vs)):
        if vs[i]==vs[i-1]+1: cur+=1; best=max(best,cur)
        else: cur=1
    return best if vs else 0
for m in sys.argv[1:]:
    G=defaultdict(set)
    for line in open(f"pal_{m}.txt"):
        a,q,s,e,t,f,n=map(int,line.split()[:7])
        G[(a,q,s,e)].add((t,n,f))
    # per group: number of distinct t, longest consecutive t run, ratio q/a
    rows=[]
    for k,v in G.items():
        ts=[x[0] for x in v]
        rows.append((longest_run(ts),len(set(ts)),k))
    rows.sort(reverse=True)
    ng=len(G); multi=sum(1 for r in rows if r[1]>1)
    print(f"order {m}: groups (a,q,s,e) {ng}, with >1 distinct t: {multi}")
    for r in rows[:8]:
        a,q,s,e=r[2]; print(f"   run {r[0]} distinct-t {r[1]}  a={a} q={q} q/a={q/a:.2f} s={s} e={e}")
