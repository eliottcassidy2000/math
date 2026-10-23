"""Independent adjacency-matrix audit of complement packet, orders1..5."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
import json


def run():
    checks=0
    rows=[]
    def ensure(ok, message):
        nonlocal checks
        checks+=1
        if not ok:
            raise RuntimeError(message)
    for n in range(1,6):
        pairs=list(combinations(range(n),2))
        orders=list(permutations(range(n)))
        def matrix(code, tournament=False):
            a=[[0]*n for _ in range(n)]
            for (i,j),value in zip(pairs,code):
                a[i][j]=value
                a[j][i]=1-value if tournament else value
            return a
        def code(a,p):
            return tuple(a[p[i]][p[j]] for i,j in pairs)
        def canonical(a):
            return min(code(a,p) for p in orders)
        allcodes=[tuple((mask>>i)&1 for i in range(len(pairs))) for mask in range(1<<len(pairs))]
        graphcanonical={g:canonical(matrix(g)) for g in allcodes}
        tournamentcanonical={t:canonical(matrix(t,True)) for t in allcodes}
        greps=sorted(set(graphcanonical.values()))
        treps=sorted(set(tournamentcanonical.values()))
        incidence={}
        triangle_hist={}
        single=[]
        for g in greps:
            a=matrix(g)
            ga=sum(code(a,p)==g for p in orders)
            hist=Counter()
            triangles=Counter()
            for order in orders:
                pos={v:i for i,v in enumerate(order)}
                t=[[0]*n for _ in range(n)]
                for i,j in pairs:
                    t[i][j]=int((pos[i]<pos[j]) != bool(a[i][j]))
                    t[j][i]=1-t[i][j]
                tc=code(t,tuple(range(n)))
                hist[tournamentcanonical[tc]]+=1
                cyc=sum((t[i][j] and t[j][k] and t[k][i]) or
                        (t[j][i] and t[k][j] and t[i][k])
                        for i,j,k in combinations(range(n),3))
                triangles[cyc]+=1
            triangle_hist[g]=tuple(sorted(triangles.items()))
            incidence[g]={t:v//ga for t,v in hist.items()}
            ensure(all(v%ga==0 for v in hist.values()),"graph integrality")
            degree=[sum(row) for row in a]
            ensure(Fraction(sum(k*v for k,v in triangles.items()),len(orders)) ==
                   Fraction(sum(d*(n-1-d) for d in degree),6),"expectation")
            if len(hist)==1:
                single.append(g)
            for u,v in pairs:
                e=a[u][v]
                before=[w for w in range(n) if w not in (u,v) and a[u][w]!=a[v][w] and a[u][w]!=e]
                after=[w for w in range(n) if w not in before and w not in (u,v)]
                cycle_values=[]
                for order in (before+[u,v]+after,before+[v,u]+after):
                    pos={v:i for i,v in enumerate(order)}
                    t=[[int(i!=j and ((pos[i]<pos[j]) != bool(a[i][j]))) for j in range(n)] for i in range(n)]
                    cycle_values.append(sum((t[i][j] and t[j][k] and t[k][i]) or
                                            (t[j][i] and t[k][j] and t[i][k])
                                            for i,j,k in combinations(range(n),3)))
                ensure(cycle_values[1]-cycle_values[0] == sum(a[u][w]!=a[v][w] for w in range(n) if w not in (u,v)),"adjacent swap")
        ensure(set(single)=={tuple(0 for _ in pairs),tuple(1 for _ in pairs)},"independent types")
        for t in treps:
            a=matrix(t,True)
            ta=sum(code(a,p)==t for p in orders)
            hist=Counter()
            for order in orders:
                pos={v:i for i,v in enumerate(order)}
                g=tuple(int(bool(a[i][j]) != (pos[i]<pos[j])) for i,j in pairs)
                hist[graphcanonical[g]]+=1
            for g in greps:
                ensure(hist[g]%ta==0 and hist[g]//ta==incidence[g].get(t,0),"independent incidence")
        cpairs={min(g,graphcanonical[tuple(1-x for x in g)]) for g in greps}
        distinct=len({triangle_hist[g] for g in cpairs})
        ensure(distinct==len(cpairs),"triangle histogram reconstruction in bounded universe")
        rows.append({"n":n,"graph_types":len(greps),"tournament_types":len(treps),
                     "complement_pairs":len(cpairs),"distinct_triangle_histograms":distinct})
    return {"status":"PASS","checks":checks,"universe":"all graph/tournament matrices and all total orders, n<=5","census":rows}


if __name__=="__main__":
    result=run()
    print(json.dumps(result,sort_keys=True))
