#!/usr/bin/env python3
"""Independent literal full tournament cube, not the producer orbit generator."""
import itertools, json, hashlib

def require(ok, message):
    if not ok: raise RuntimeError(message)

def check(n):
    p=n//2; q=list(reversed(range(n))); pairs=list(itertools.combinations(range(n),2))
    anti=balanced=path_total=cycle_total=0; minimum=None; digest=hashlib.sha256()
    perms=list(itertools.permutations(range(n)))
    for mask in range(1<<len(pairs)):
        a=[[0]*n for _ in range(n)]
        for k,(u,v) in enumerate(pairs):
            a[u][v]=(mask>>k)&1; a[v][u]=1-a[u][v]
        if not all(a[u][v]==a[q[v]][q[u]] for u,v in pairs): continue
        anti+=1; d=list(map(sum,a))
        if sorted(d)!=[p-1]*p+[p]*p: continue
        balanced+=1
        b=[row+[int(d[u]==p-1)] for u,row in enumerate(a)]+[[int(x==p) for x in d]+[0]]
        qq=q+[n]
        require(all(sum(row)==p for row in b),'extension not regular')
        require(all(b[u][v]==b[qq[v]][qq[u]] for u in range(n+1) for v in range(n+1)),'extension not anti')
        paths=set(); cycles=set(); centers=set()
        for word in perms:
            if not all(a[u][v] for u,v in zip(word,word[1:])): continue
            # Test fixed-path C degree criterion for every HP, independently of q.
            c=[[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1,n):
                    c[i][j]=1 if j==i+1 else 1-a[word[i]][word[j]]
                    c[j][i]=1-c[i][j]
            internal=sorted(map(sum,c))==[p-1]*p+[p]*p
            require(internal==(d[word[0]]==p and d[word[-1]]==p-1),'internal endpoint criterion')
            symmetric=all(q[word[i]]==word[-1-i] for i in range(n))
            if symmetric and internal:
                paths.add(word); centers.add(tuple(sorted(word[p-1:p+1])))
            cycle=(n,)+word
            arcs=set(zip(cycle,cycle[1:]+cycle[:1]))
            if all(b[u][v] for u,v in arcs) and {(qq[v],qq[u]) for u,v in arcs}==arcs:
                cycles.add(word)
                fixed=[(u,v) for u,v in arcs if (qq[v],qq[u])==(u,v)]
                require(fixed==[(word[p-1],word[p])],'unique middle arc')
        require(paths==cycles,'rooted cycle/path bijection')
        require(centers=={tuple(sorted((u,q[u]))) for u in range(n)},'missing center')
        minimum=len(paths) if minimum is None else min(minimum,len(paths))
        path_total+=len(paths); cycle_total+=len(cycles)
        digest.update(json.dumps([mask,sorted(paths)],separators=(',',':')).encode())
    require((anti,balanced,minimum)=={2:(2,2,1),4:(16,8,2),6:(512,112,5)}[n],'finite counts')
    return dict(n=n,full_matrices=1<<len(pairs),anti_invariant_matrices=anti,balanced_matrices=balanced,minimum_symmetric_downhill_paths=minimum,total_paths=path_total,total_cycles=cycle_total,every_paired_arc_realized=True,universe_paths_sha256=digest.hexdigest())
if __name__=='__main__':
    print(json.dumps({'status':'PASS','scope':'FINITE-EXACT; no all-order proof','rows':[check(n) for n in (2,4,6)]},sort_keys=True,indent=2))
