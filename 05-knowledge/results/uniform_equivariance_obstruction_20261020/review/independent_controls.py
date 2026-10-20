#!/usr/bin/env python3
"""Independent tiny exact controls; no imports from producer or repository code."""
from itertools import combinations, permutations
from math import factorial
import json


def parity(a): return a.bit_count() & 1

def valuation(a, p):
    s=0
    while a % p == 0:
        s += 1
        a //=p
    return s


def determinant(a):
    a=[r[:] for r in a]
    sign=1; prev=1; n=len(a)
    for k in range(n-1):
        if a[k][k]==0:
            q=next(i for i in range(k+1,n) if a[i][k])
            a[k],a[q]=a[q],a[k]; sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=a[i][j]*pivot-a[i][k]*a[k][j]
                assert z % prev == 0
                a[i][j]=z//prev
            a[i][k]=0
        prev=pivot
    return sign*a[-1][-1] if n else 1


def graph_data(n):
    edges=list(combinations(range(n),2)); m=len(edges); size=1<<m
    ei={e:i for i,e in enumerate(edges)}
    ps=list(permutations(range(n)))
    images=[]; invs=[]; affines=[]
    for p in ps:
        dest=[ei[tuple(sorted((p[i],p[j])))] for i,j in edges]
        inv=sum(1<<k for k,(i,j) in enumerate(edges) if p[i]>p[j])
        row=[0]*size
        for x in range(1,size):
            bit=x & -x; k=bit.bit_length()-1
            row[x]=row[x^bit] ^ (1<<dest[k])
        images.append(row); invs.append(inv)
        affines.append([row[x] ^ row[inv] for x in range(size)])
    gens=[]
    for i in range(n-1):
        p=list(range(n)); p[i],p[i+1]=p[i+1],p[i]
        gens.append(ps.index(tuple(p)))
    return edges,ps,images,invs,affines,gens


def invariant_subspaces(size,images,gens):
    # Breadth-first enumerate all vector subspaces, independent of invariance.
    spaces={frozenset([0])}; frontier=list(spaces)
    while frontier:
        u=frontier.pop()
        for x in range(1,size):
            if x not in u:
                v=frozenset(u | {x^a for a in u})
                if v not in spaces:
                    spaces.add(v); frontier.append(v)
    return sorted((u for u in spaces if all(all(images[g][x] in u for x in u) for g in gens)),key=lambda u:(len(u),sorted(u))),len(spaces)


def cut_space(n,edges):
    return frozenset(sum(1<<k for k,(i,j) in enumerate(edges) if ((s>>i)^(s>>j))&1) for s in range(1<<n))


def analyze(n,u,data,matching=False):
    edges,ps,images,invs,affines,gens=data
    size=1<<len(edges); c=size//len(u); gorder=factorial(n)
    canon=[min(t^x for x in u) for t in range(size)]
    classes=sorted(set(canon)); unseen=set(classes); a=[]; astab=[]
    while unseen:
        t=min(unseen)
        orb={canon[row[t]] for row in affines}; unseen-=orb
        stabilizer={g for g,row in enumerate(affines) if canon[row[t]]==t}
        a.append(t); astab.append(stabilizer)
        assert len(orb)*len(stabilizer)==gorder
    dual=[x for x in range(size) if all(parity(x&y)==0 for y in u)]
    assert len(dual)==c
    unseen=set(dual); b=[]; bstab=[]; coeffs=[]
    while unseen:
        x=min(unseen)
        orb={row[x] for row in images}; unseen-=orb
        st={g for g,row in enumerate(images) if row[x]==x}
        assert len(orb)*len(st)==gorder
        if any(parity(x&invs[g]) for g in st): continue
        coeff={}
        for g,row in enumerate(images):
            y=row[x]; sign=(-1)**parity(x&invs[g])
            assert coeff.get(y,sign)==sign
            coeff[y]=sign
        b.append(x); bstab.append(st); coeffs.append(coeff)
    assert len(a)==len(b)
    M=[[sum(s*((-1)**parity(x&t)) for x,s in co.items()) for co in coeffs] for t in a]
    N=len(a)
    # Both exact orthogonality identities, with integral scalings.
    for i in range(N):
        for j in range(N):
            left=sum((gorder//len(astab[k]))*M[k][i]*M[k][j] for k in range(N))
            right=c*(gorder//len(bstab[i]))*(i==j)
            assert left==right,(n,len(u),'first gram',i,j,left,right)
            left=sum(len(bstab[k])*M[i][k]*M[j][k] for k in range(N))
            right=c*len(astab[i])*(i==j)
            assert left==right,(n,len(u),'second gram',i,j,left,right)
    det=determinant(M); assert det and abs(det)&(abs(det)-1)==0
    prod_a=prod_b=1
    for st in astab: prod_a*=len(st)
    for st in bstab: prod_b*=len(st)
    assert det*det*prod_b==c**N*prod_a
    vals={str(p):sum(valuation(len(st),p) for st in astab) for p in (3,5)}
    assert all(vals[str(p)]==sum(valuation(len(st),p) for st in bstab) for p in (3,5))
    result={'n':n,'dim_U':len(u).bit_length()-1,'classes':c,'types':N,'determinant':det,'total_odd_valuations':vals}
    if matching:
        # Any relabelled target may represent a target type. Test every image,
        # so containment is in the required source <= target direction.
        adj=[]
        for st in astab:
            possible=[]
            for j,x in enumerate(b):
                if any(all(images[h][y]==y for h in st) for y in {row[x] for row in images}): possible.append(j)
            adj.append(possible)
        assigned={}
        def augment(i,seen):
            for j in adj[i]:
                if j in seen: continue
                seen.add(j)
                if j not in assigned or augment(assigned[j],seen):
                    assigned[j]=i; return True
            return False
        maxmatch=sum(augment(i,set()) for i in range(N))
        result['maximum_containment_matching']=maxmatch
        assert (maxmatch==N)==(n<=4)
        result['source_stabilizer_orders']=[len(st) for st in astab]
        result['target_stabilizer_orders']=[len(st) for st in bstab]
    return result


def main():
    out={'all_invariant_quotients':[], 'tournament_and_switching_controls':[]}
    for n in range(5):
        data=graph_data(n); size=1<<len(data[0])
        spaces,total=invariant_subspaces(size,data[2],data[5])
        checks=[analyze(n,u,data) for u in spaces]
        out['all_invariant_quotients'].append({'n':n,'all_subspaces':total,'invariant_subspaces':len(spaces),'checks':checks})
    for n in range(6):
        data=graph_data(n)
        for kind,u in [('tournaments',frozenset([0])),('switching',cut_space(n,data[0]))]:
            check=analyze(n,u,data,matching=True); check['kind']=kind
            out['tournament_and_switching_controls'].append(check)
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
