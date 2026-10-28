#!/usr/bin/env python3
"""Exact, standard-library, no assert guards. All 128 realizations of fixed H_1."""
from itertools import combinations, product
from collections import Counter
from fractions import Fraction
import json

def require(ok, why):
    if not ok: raise RuntimeError(why)

def build(bits):
    root=bits[0]; a=[[0]*12 for _ in range(12)]
    for i,j in combinations(range(12),2):
        b,u=divmod(i,4); c,v=divmod(j,4)
        if b!=c: arc=int((c-b)%3==1)^root
        elif u==0: arc=1^bits[1+b]
        else: arc=int((v-u)%3==1)^bits[4+b]
        a[i][j]=arc; a[j][i]=1-arc
    return a

def triangles(a):
    return frozenset(t for t in combinations(range(len(a)),3) if a[t[0]][t[1]]==a[t[1]][t[2]]==a[t[2]][t[0]])

def encode(a):return tuple(a[i][j] for i,j in combinations(range(len(a)),2))

def move(a,p):
    b=[[0]*len(a) for _ in a]
    for i in range(len(a)):
        for j in range(len(a)): b[p[i]][p[j]]=a[i][j]
    return encode(b)

def hypergraph_auts(n,H):
    # Global permutations with only incidence-based pruning: no tree assumption.
    deg=[sum(i in e for e in H) for i in range(n)]
    pair=[[sum(i in e and j in e for e in H) for j in range(n)] for i in range(n)]
    order=sorted(range(n),key=lambda i:deg.count(deg[i])); p=[-1]*n; used=set(); out=[]
    def visit(k):
        if k==n: out.append(tuple(p)); return
        i=order[k]
        for j in range(n):
            if j in used or deg[i]!=deg[j]:continue
            if any(pair[i][u]!=pair[j][p[u]] for u in order[:k]):continue
            if any((tuple(sorted((i,u,v))) in H)!=(tuple(sorted((j,p[u],p[v]))) in H) for u,v in combinations(order[:k],2)):continue
            p[i]=j; used.add(j); visit(k+1); used.remove(j)
    visit(0); return out

def ordinary_aut(a):
    n=len(a); deg=list(map(sum,a)); order=sorted(range(n),key=lambda i:deg.count(deg[i])); p=[-1]*n; used=set(); count=0
    def visit(k):
        nonlocal count
        if k==n:count+=1; return
        i=order[k]
        for j in range(n):
            if j in used or deg[i]!=deg[j]:continue
            if any(a[i][u]!=a[j][p[u]] for u in order[:k]):continue
            p[i]=j; used.add(j); visit(k+1); used.remove(j)
    visit(0); return count

def main():
    choices={bits:build(bits) for bits in product(range(2),repeat=7)}
    bank={encode(a):bits for bits,a in choices.items()}; require(len(bank)==128,'distinct labelled realizations')
    H=triangles(choices[(0,)*7]); require(all(triangles(a)==H for a in choices.values()),'common literal H')
    G=hypergraph_auts(12,H); require(len(G)==1296,'independent hypergraph group')
    # Signed node-action fixed counts, independently compared with actual arc action.
    def parity(p): return sum(p[i]>p[j] for i,j in combinations(range(len(p)),2))%2
    hist=Counter(); kernel=0
    for g in G:
        pi=tuple(g[4*i]//4 for i in range(3))
        signs=tuple(parity(tuple(g[4*i+u]%4-1 for u in (1,2,3))) for i in range(3))
        cycles=[]; left=set(range(3))
        while left:
            start=min(left); cyc=[]; j=start
            while j not in cyc:cyc.append(j); left.remove(j); j=pi[j]
            cycles.append(cyc)
        predicted=0 if parity(pi) or any(sum(signs[i] for i in c)%2 for c in cycles) else 2*4**len(cycles)
        observed=sum(move(a,g)==encode(a) for a in choices.values())
        require(predicted==observed,'signed node orbit fixed-count law')
        hist[observed]+=1; kernel+=int(observed==128)
    require(hist==Counter({0:1053,8:216,128:27}),'fixed-point histogram')
    require(sum(k*v for k,v in hist.items())==4*len(G),'Burnside full fiber')
    unseen=set(bank); rows=[]; distribution=Counter()
    while unseen:
        key=min(unseen); bits=bank[key]; a=choices[bits]; orbit={move(a,g) for g in G}
        require(orbit<=set(bank),'fiber action closure'); unseen-=orbit
        aut=ordinary_aut(a); stab=sum(move(a,g)==key for g in G)
        require(aut==stab and len(orbit)*stab==len(G),'direct ordinary group = stabilizer')
        rows.append({'representative_bits':bits,'orbit_size':len(orbit),'ordinary_aut_order':aut})
        distribution[aut]+=len(orbit)
    require(sorted((r['orbit_size'],r['ordinary_aut_order']) for r in rows)==[(16,81),(16,81),(48,27),(48,27)],'full fiber census')
    # Every one of the 128 direct arc-based automorphism counts, not just representatives.
    require(Counter(ordinary_aut(a) for a in choices.values())==distribution,'independent all-realization group distribution')
    mass=sum((Fraction(1,r['ordinary_aut_order']) for r in rows),Fraction())
    require(mass==Fraction(128,1296)==Fraction(8,81),'orbit mass')
    # Nontrivial hostile controls: hypergraph symmetries need not preserve a realization.
    require(any(move(choices[(0,)*7],g)!=encode(choices[(0,)*7]) for g in G),'not every H automorphism lifts')
    # All-height mass formula checked against independent polynomial recurrence.
    poly=Counter({0:2}); height_rows=[]
    for h in range(5):
        L=3**h; P=(3**(h+1)-1)//2
        mass_h=sum((Fraction(c,3**(L+e)) for e,c in poly.items()),Fraction())
        require(mass_h==Fraction(2**L,3**P),'all-height mass specialization')
        height_rows.append({'height':h,'isotypes':sum(poly.values()),'labelled_fiber':2**(L+P),'hypergraph_aut_order':6**P,'mass':str(mass_h)})
        cube=Counter()
        for i,a in poly.items():
            for j,b in poly.items():
                for k,c in poly.items():cube[i+j+k]+=a*b*c
        for i,a in poly.items():cube[3*i]-=a
        require(all(c%3==0 for c in cube.values()),'polynomial integrality')
        q=Counter({i:c//3 for i,c in cube.items() if c})
        for i,a in poly.items():q[3*i+1]+=a
        poly=q
    print(json.dumps({'all_height_mass_checks':height_rows,'labelled_fiber':len(bank),'hyperedges':len(H),'hypergraph_aut_order':len(G),'isotypes':rows,'labelled_aut_order_distribution':dict(sorted(distribution.items())),'orbit_mass':str(mass),'fixed_count_histogram':dict(sorted(hist.items())),'action_kernel_order':kernel,'all_checks':'PASS'},indent=2))
if __name__=='__main__': main()
