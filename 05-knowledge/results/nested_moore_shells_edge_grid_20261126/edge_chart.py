import argparse, json
from collections import Counter
from itertools import combinations
from pathlib import Path
from nested_core_check import graph, orbit_chart

def main():
    C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
    A=graph(C,S)
    for h,y in json.loads(Path("rigid_51_contacts.json").read_text()):
        v=58+56*h+y
        A[115].add(v);A[v].add(115)
    v0,v1=0,57
    assert v1 in A[v0]
    X=list(sorted(set(range(3250))-{v0,v1}-A[v0]-A[v1]))
    assert len(X)==3136
    xset=set(X)
    row={x:(x-58)//56 for x in X}
    col={x:(x-58)%56 for x in X}
    anchor=58+56*56
    for x in X:
        y,z=1+row[x],anchor+col[x]
        assert A[x]-xset=={y,z}
    r,c=Counter(),Counter()
    byrow=Counter(); bycol=Counter()
    adjacent={x:A[x]&xset for x in X}
    for x,N in adjacent.items():
        assert len(N)==len({row[y] for y in N})==len({col[y] for y in N})
        assert all(row[x]!=row[y] and col[x]!=col[y] for y in N)
        if len(N)==55:byrow[row[x]]+=1;bycol[col[x]]+=1
        for y in N:
            if y>x:
                r[tuple(sorted((row[x],row[y])))]+=1;c[tuple(sorted((col[x],col[y])))]+=1
    profiles={}
    for a in [0,28]:
        modes=[]; funcs=[]
        for b in range(56):
            if b==a:
                funcs.append(tuple(range(56)))
                continue
            f=[]
            for j in range(56):
                vs=[y for y in adjacent[58+56*a+j] if row[y]==b]
                assert len(vs)==1
                f.append(col[vs[0]])
            funcs.append(tuple(f))
            dif=Counter((z-j)%56 for j,z in enumerate(f))
            modes.append((b,max(dif.values())))
        cv=set(funcs)
        closed=all(tuple(t[j] for j in u) in cv for t in cv for u in cv)
        quotients={}
        for m in [2,4,7,8,14,28]:
            quotients[m]=sum(len({(y-x)%m for x,y in enumerate(t)})>1 for t in cv)
        profiles[a]={"projection_defects":quotients,"literal_translations":sum(n==56 for b,n in modes),
                    "translation_disagreements":sum(56-n for b,n in modes),
                    "best_shift_agreement_hist":dict(sorted(Counter(n for b,n in modes).items())),
                    "forms_group":closed}
    assert len(r)<=1540 and len(c)<=1540
    assert max(r.values())<=56 and max(c.values())<=56
    print(json.dumps({"row_cyclic_profiles":profiles,"residual_vertices":len(X),"residual_edges":sum(r.values()),
         "degree_hist":dict(sorted(Counter(map(len,adjacent.values())).items())),
         "rows_all_56_sheets_full":sorted(i for i,n in byrow.items() if n==56),
         "column_full_sheets_hist":dict(sorted(Counter(bycol[j] for j in range(56)).items())),
         "full_row_matchings":sum(v==56 for v in r.values()),
         "full_column_matchings":sum(v==56 for v in c.values()),
         "row_pair_occupancy":dict(sorted(Counter(r.get(k,0) for k in combinations(range(56),2)).items())),
         "column_pair_occupancy":dict(sorted(Counter(c.get(k,0) for k in combinations(range(56),2)).items()))},indent=2))
if __name__=="__main__":main()
