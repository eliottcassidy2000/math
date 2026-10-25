#!/usr/bin/env python3
"""Exact small graph/tournament quotients; no assertions, no third party modules."""
import itertools,json,math
from collections import Counter,defaultdict
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def pairs(n): return list(itertools.combinations(range(n),2))
def relabel(x,n,p,tournament=False):
    ps=pairs(n); idx={e:k for k,e in enumerate(ps)}; y=0
    for k,(i,j) in enumerate(ps):
        a,b=p[i],p[j]; bit=(x>>k)&1
        if tournament and a>b: bit^=1
        if bit:y|=1<<idx[tuple(sorted((a,b)))]
    return y

def main():
    rows=[]; witnesses={}; distributions={}
    for n in range(1,6):
        N=n*(n-1)//2; full=(1<<N)-1; perms=list(itertools.permutations(range(n)))
        gc={};tc={}
        for x in range(1<<N):
            if x not in gc:
                orb={relabel(x,n,p) for p in perms}; r=min(orb);gc.update(dict.fromkeys(orb,r))
            if x not in tc:
                orb={relabel(x,n,p,True) for p in perms};r=min(orb);tc.update(dict.fromkeys(orb,r))
        gpair=lambda x:min(gc[x],gc[x^full])
        tpair=lambda x:min(tc[x],tc[x^full])
        gr=set(gc.values());tr=set(tc.values())
        rows.append(dict(n=n,labelled=1<<N,graphs=len(gr),graph_complement_pairs=len({gpair(x) for x in gr}),tournaments=len(tr),self_converse=sum(tc[x^full]==x for x in tr),converse_pairs=len({tpair(x) for x in tr})))
        # Independent ordered-arc relabelling path and complement/reversal controls.
        for x in range(1<<N):
            for p in perms:
                arcs={(p[j],p[i]) if x>>k&1 else (p[i],p[j]) for k,(i,j) in enumerate(pairs(n))}
                direct=sum(1<<k for k,(i,j) in enumerate(pairs(n)) if (j,i) in arcs)
                require(direct==relabel(x,n,p,True),'independent arc relabelling')
                require(relabel(x^full,n,p,True)==(direct^full),'complement intertwines')
        # Same graph orbit, distinct tournament even modulo converse.
        if 'same_graph_distinct_tournaments' not in witnesses:
            for x in range(1<<N):
                for p in perms:
                    y=relabel(x,n,p)
                    if tpair(x)!=tpair(y):
                        witnesses['same_graph_distinct_tournaments']=dict(n=n,x=x,y=y,edges_x=[e for k,e in enumerate(pairs(n)) if x>>k&1],edges_y=[e for k,e in enumerate(pairs(n)) if y>>k&1]);break
                if 'same_graph_distinct_tournaments' in witnesses:break
        # Converse-compatible, order-averaged distribution of graph pair types.
        profiles={}
        for t in sorted({tpair(x) for x in tr}):
            hist=Counter(gpair(relabel(t,n,p,True)) for p in perms)
            require(sum(hist.values())==math.factorial(n),'profile mass')
            h2=Counter(gpair(relabel(t^full,n,p,True)) for p in perms)
            require(hist==h2,'converse average')
            require(n<2 or all(v%2==0 for v in hist.values()),'order reversal pairs')
            profiles[t]=tuple(sorted(hist.items()))
        distributions[n]={'distinct_profiles':len(set(profiles.values())),'converse_classes':len(profiles),'profiles':profiles}
        # Twisted Burnside: every even-order permutation fixes no tournament.
        fixed_sum=0
        for p in perms:
            actual=sum(relabel(x,n,p,True)==x for x in range(1<<N))
            seen=set();lengths=[]
            for i in range(n):
                if i not in seen:
                    j=i;l=0
                    while j not in seen:seen.add(j);l+=1;j=p[j]
                    lengths.append(l)
            odd=all(l%2 for l in lengths)
            seenpairs=set();orbs=0
            for e in pairs(n):
                if e not in seenpairs:
                    orbs+=1;f=e
                    while f not in seenpairs:seenpairs.add(f);f=tuple(sorted((p[f[0]],p[f[1]])))
            require(actual==(2**orbs if odd else 0),'Burnside fixed formula')
            fixed_sum+=actual
        require(fixed_sum==math.factorial(n)*len(tr),'Burnside total')
        if n==4:
            # Hamilton chart: path 0->1->2->3; backward chord bits.
            chords=[e for e in pairs(n) if e[1]>e[0]+1]; ix={e:k for k,e in enumerate(pairs(n))};hm=Counter()
            chart=[]
            for b in range(1<<len(chords)):
                x=sum((1<<ix[e]) for k,e in enumerate(chords) if b>>k&1)
                hm[tc[x]]+=1;chart.append([b,x,tc[x]])
            require(sorted(hm.values())==[1,1,1,5],'marked path fibres')
            require(chart[2][2]==chart[3][2] and chart[5][2]!=chart[4][2],'K depends on marked path')
            witnesses['marked_n4']={'fibres':dict(hm),'chart':chart,'chords':chords}
    require([r['tournaments'] for r in rows]==[1,1,2,4,12],'small tournament counts')
    require(witnesses['same_graph_distinct_tournaments']['n']==3,'minimal obstruction')
    result={'counts':rows,'witnesses':witnesses,'order_average':distributions}
    out=Path(__file__).with_name('certificate.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
