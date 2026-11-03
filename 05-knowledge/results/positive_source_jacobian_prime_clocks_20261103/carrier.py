"""Exact signed-forest quotient and all one-birth pullbacks at old order six."""
from itertools import combinations, permutations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction as Q
import json


def need(ok,msg):
    if not ok: raise RuntimeError(msg)


def dimension(n):
    total=1
    for j in range(1,n//2+1):
        match=factorial(2*j)//(2**j*factorial(j))
        for m in range(n-2*j+1):
            parity=sum(comb(j,w) for w in range(min(j,m)+1) if (m-w)%2==0)
            total+=comb(n,2*j)*match*comb(n-2*j,m)*parity
    return total


def forest_groups(n):
    es=list(combinations(range(n),2)); ix={e:i for i,e in enumerate(es)}
    masks=set()
    # Every undirected linear forest has an ordering of its components and
    # an orientation of each path, hence occurs as a cut permutation.
    for p in permutations(range(n)):
        bits=[1<<ix[min(a,b),max(a,b)] for a,b in zip(p,p[1:])]
        for take in range(1<<(n-1)):
            masks.add(sum(b for j,b in enumerate(bits) if take>>j&1))
    groups=defaultdict(list); finer=set()
    for mask in sorted(masks):
        nb=[[] for _ in range(n)]
        for k,(a,b) in enumerate(es):
            if mask>>k&1: nb[a].append(b); nb[b].append(a)
        used=set(); pairs=[]; old=[]; sign=1
        for a in range(n):
            if a in used or len(nb[a])!=1:continue
            p=[a]; prev=-1; cur=a
            while True:
                nxt=[b for b in nb[cur] if b!=prev]
                if not nxt:break
                prev,cur=cur,nxt[0]; p.append(cur)
            used.update(p)
            if p[0]>p[-1]:p.reverse()
            sign*=(-1)**sum(a>b for a,b in zip(p,p[1:]))
            pairs.append((p[0],p[-1],(len(p)-2)%2))
            old.append((tuple(sorted(p)),p[0],p[-1]))
        isolated=tuple(i for i in range(n) if not nb[i])
        key=(isolated,tuple(sorted(pairs)))
        groups[key].append((mask,sign))
        finer.add(tuple(sorted(old)))
    return es,dict(groups),len(finer)


def main():
    es,g6,b6=forest_groups(6); fs,g7,b7=forest_groups(7)
    need(len(g6)==dimension(6)==571,'six quotient dimension')
    need(b6==661,'inherited six bound')
    need(len(g7)==dimension(7),'seven dimension')
    l=[1]
    for n in range(1,8):
        l.append(sum(comb(n-1,k-1)*(1 if k==1 else factorial(k)//2)*l[n-k] for k in range(1,n+1)))
    need(sum(map(len,g6.values()))==l[6] and sum(map(len,g7.values()))==l[7],'complete forests')
    old_to_new=[fs.index(e) for e in es]
    new_edges=[fs.index((v,6)) for v in range(6)]
    decode={}
    for k,terms in g6.items():
        for mask,sg in terms:decode[mask]=(k,sg)
    translated={}
    for k,terms in g7.items():
        translated[k]=[(sum(((mask>>i)&1)<<j for j,i in enumerate(old_to_new)),
                        sum(((mask>>i)&1)<<j for j,i in enumerate(new_edges)),sg)
                       for mask,sg in terms]
    tested=0
    for cut in range(64):
        for key,terms in translated.items():
            poly=defaultdict(int)
            for mask,zmask,sg in terms:
                # z_(v,6)=+1 iff old vertex v beats newborn, i.e. v not in cut.
                poly[mask]+=sg*(-1)**((cut&zmask).bit_count())
            groups_seen={}; sizes=defaultdict(int)
            for mask,c in poly.items():
                if not c:continue
                need(mask in decode,'old forest support')
                group,sg=decode[mask]; value=c*sg
                if group in groups_seen:need(groups_seen[group]==value,'constant signed coefficients')
                else:groups_seen[group]=value
                sizes[group]+=1
            for group,size in sizes.items():need(size==len(g6[group]),'no missing terms in nonzero old group')
            tested+=1
    ix={e:i for i,e in enumerate(es)}
    edge_mask=lambda paths:sum(1<<ix[min(a,b),max(a,b)] for p in paths for a,b in zip(p,p[1:]))
    f1=edge_mask([(0,4,1),(2,5,3)]); f2=edge_mask([(0,5,1),(2,4,3)])
    need(f1!=f2 and decode[f1]==decode[f2],'new source-law collision')
    # Independent coefficientwise expansion of the displayed formal EGF.
    N=20
    def product(a,b):
        return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(N+1)]
    ex=[Q(1,factorial(i)) for i in range(N+1)]
    ea=[Q(1,2**(i//2)*factorial(i//2)) if i%2==0 else Q(0) for i in range(N+1)]
    co=[ex[i] if i%2==0 else Q(0) for i in range(N+1)]
    double=[Q(0)]*(N+1)
    for w in range(N//3+1):
        for r in range((N-3*w)//2+1):
            double[3*w+2*r]+=Q(1,2**w*factorial(w)*factorial(w+2*r))
    core=product(ea,double)
    core=[v-co[i]+int(i==0) for i,v in enumerate(core)]
    egf=product(ex,core)
    need(all(egf[n]*factorial(n)==dimension(n) for n in range(N+1)),'EGF coefficient agreement')
    print(json.dumps(dict(dimensions_0_to_20=[dimension(n) for n in range(21)],
                         n6=dict(old_bound=b6,new_bound=len(g6),forest_terms=l[6]),
                         n7=dict(old_bound=b7,new_bound=len(g7),forest_terms=l[7]),
                         all_birth_cuts=64,exact_basis_pullbacks=tested,
                         hostile_masks=[f1,f2],hostile_signs=[decode[f1][1],decode[f2][1]],
                         egf_coefficients_verified_through=N,
                         scope='closed rational upper bound; all-order sharpness unproved'),indent=2))


if __name__=='__main__':main()
