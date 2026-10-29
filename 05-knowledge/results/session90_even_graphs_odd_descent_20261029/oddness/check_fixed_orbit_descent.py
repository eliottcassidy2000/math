#!/usr/bin/env python3
"""Exact coprime fixed-orbit controls and folded path observer. Python stdlib only.
Independent formulas; no repository imports or network. All tests stay active under -O.
"""
from itertools import permutations, combinations
from math import gcd, factorial
import json

CHECKS=0
def need(x,msg='check failed'):
    global CHECKS
    CHECKS+=1
    if not x: raise AssertionError(msg)
def comp(p,q): return tuple(p[q[i]] for i in range(len(p)))
def inv(p): return tuple(p.index(i) for i in range(len(p)))
def square_root(p):
    out=[None]*len(p); seen=set()
    for i in range(len(p)):
        if i in seen: continue
        cyc=[]; j=i
        while j not in seen: cyc.append(j);seen.add(j);j=p[j]
        need(len(cyc)%2==1,'odd-order root required')
        for k,v in enumerate(cyc): out[v]=cyc[(k+(len(cyc)+1)//2)%len(cyc)]
    out=tuple(out); need(comp(out,out)==p); return out

def lift(m,n):
    rows=[0]*(2*n+1)
    def edge(a,b,bit): rows[a if bit else b] |= 1 << (b if bit else a)
    for i in range(n):
        p,q=1+i,1+n+i; d=(m>>(i*n+i))&1
        edge(p,q,1); edge(0,p,d); edge(0,q,1-d)
    for i,j in combinations(range(n),2):
        x=(m>>(i*n+j))&1; y=(m>>(j*n+i))&1
        edge(1+i,1+j,x);edge(1+n+i,1+n+j,1-x)
        edge(1+i,1+n+j,x^y);edge(1+j,1+n+i,x^y)
    return tuple(rows),(0,)+tuple(range(n+1,2*n+1))+tuple(range(1,n+1))
def arc(rows,a,b): return (rows[a]>>b)&1
def is_aut(rows,p): return all(arc(rows,i,j)==arc(rows,p[i],p[j]) for i,j in combinations(range(len(rows)),2))
def folded_paths(rows,tau):
    root=next(i for i,x in enumerate(tau) if i==x)
    pairs=[(i,tau[i]) for i in range(len(rows)) if i<tau[i]]
    n=len(pairs)
    if not n: return 1
    vertices=sum((list(pair) for pair in pairs),[])
    idx={v:i for i,pair in enumerate(pairs) for v in pair}
    dp=[{} for _ in range(1<<n)]
    for v in vertices: dp[1<<idx[v]][v]=1
    for mask in range(1,1<<n):
        for last,count in dp[mask].items():
            for nxt in vertices:
                bit=1<<idx[nxt]
                if not mask&bit and arc(rows,last,nxt):
                    target=dp[mask|bit];target[nxt]=target.get(nxt,0)+count
    return sum(count for last,count in dp[-1].items() if arc(rows,last,root))
def literal_paths(rows):
    return [p for p in permutations(range(len(rows))) if all(arc(rows,a,b) for a,b in zip(p,p[1:]))]
def actmask(p,m): return sum(1<<p[i] for i in range(len(p)) if m>>i&1)
def orbits(items,group,action):
    unseen=set(items); result=[]
    while unseen:
        a=min(unseen); orbit={action(g,a) for g in group}; need(orbit<=unseen)
        unseen-=orbit; result.append(orbit)
    return result

def cyclic_descent():
    actions=0; repairs=0
    for m in range(1,32,2):
        for d in range(1,m+1):
            if m%d: continue
            for u in range(m):
                if gcd(u,m)!=1 or (u*u-1)%m: continue
                C=[g for g in range(m) if (u*g-g)%m==0]
                for b in range(d):
                    if (u+1)*b%d: continue
                    fixed=[x for x in range(d) if (u*x+b-x)%d==0]
                    need(bool(fixed)); need({(fixed[0]+g)%d for g in C}==set(fixed))
                    for x in fixed:
                        for g in range(m):
                            y=(x+g)%d
                            if y not in fixed: continue
                            h=(u*g-g)%m; a=(h*((m+1)//2))%m; c=(g+a)%m
                            need(a%d==0);need(c in C);need((x+c)%d==y);repairs+=1
                    actions+=1
    # Even-group failures: C2 translation prevents existence; C4 inversion loses uniqueness.
    need([x for x in range(2) if (x+1)%2==x]==[])
    need({g for g in range(4) if (-g)%4==g}=={0,2})
    need(len(orbits(range(2),(0,2),lambda g,x:(g+x)%2))==2)
    return {'odd_cyclic_actions':actions,'explicit_transporter_repairs':repairs,
            'even_hostiles':['C2 regular action, commuting translation: no fixed representative',
                             'C4 on C4/C2, inversion: two centralizer orbits in one fixed class']}

def paley_control():
    n=7; rows=tuple(sum(1<<j for j in range(n) if (j-i)%7 in {1,2,4}) for i in range(n))
    tau=tuple(-i%7 for i in range(n)); G=[p for p in permutations(range(n)) if is_aut(rows,p)]
    C=[g for g in G if comp(g,tau)==comp(tau,g)]
    need(len(G)==21);need(len(C)==3)
    subset_orbits=orbits(range(128),G,actmask)
    fixed=[m for m in range(128) if actmask(tau,m)==m]
    stable=[o for o in subset_orbits if {actmask(tau,m) for m in o}==o]
    fixed_orbits=orbits(fixed,C,actmask)
    need(len(fixed_orbits)==len(stable))
    repairs=0
    for x in fixed:
        for g in G:
            y=actmask(g,x)
            if y not in fixed: continue
            ag=comp(tau,comp(g,tau));h=comp(inv(g),ag);a=square_root(h);c=comp(g,a)
            need(actmask(a,x)==x);need(c in C);need(actmask(c,x)==y);repairs+=1
    paths=literal_paths(rows);twist=lambda p:tuple(tau[v] for v in reversed(p))
    fixedpaths=[p for p in paths if twist(p)==p]; h=folded_paths(rows,tau)
    need(h==len(fixedpaths)==9)
    pathaction=lambda g,p:tuple(g[v] for v in p)
    po=orbits(paths,G,pathaction)
    stablepo=[o for o in po if {twist(p) for p in o}==o]
    fpo=orbits(fixedpaths,C,pathaction)
    need(len(paths)==189);need(len(po)==9);need(len(stablepo)==len(fpo)==3)
    need(sum(len(o) for o in stablepo)==63)
    return {'H':len(paths),'Aut':len(G),'centralizer':len(C),'fixed_paths':h,
            'ordinary_path_orbits':len(po),'blue_path_orbits':len(fpo),
            'subset_orbits':len(subset_orbits),'stable_subset_orbits':len(stable),
            'fixed_subsets':len(fixed),'transporter_repairs':repairs}

def lift_controls():
    cases=0; internal_flips=0; distributions={}
    for n in range(4):
        hist={}
        for m in range(1<<(n*n)):
            rows,tau=lift(m,n); paths=literal_paths(rows)
            fixed=sum(tuple(tau[v] for v in reversed(p))==p for p in paths)
            dp=folded_paths(rows,tau);need(dp==fixed);need(dp%2==1)
            hist[dp]=hist.get(dp,0)+1
            for bits in range(1<<n):
                rr=list(rows)
                for i in range(n):
                    if bits>>i&1:
                        p,q=1+i,1+n+i;rr[p]^=1<<q;rr[q]^=1<<p
                need(folded_paths(rr,tau)==dp);internal_flips+=1
            cases+=1
        need(sum(k*v for k,v in hist.items())==factorial(n)*(1<<(n*n)))
        distributions[n]=hist
    return {'all_matrix_lifts_through_n':3,'matrices':cases,'within_pair_flip_controls':internal_flips,
            'labelled_fixed_path_histograms':distributions}

def joint_and_path_hostiles():
    # Each marginal C3 torsor has a fixed quotient class; only one of three joint classes does.
    pairs=[(x,y) for x in range(3) for y in range(3)]
    joint=orbits(pairs,range(3),lambda g,z:((z[0]+g)%3,(z[1]+g)%3))
    twist=lambda z:(-z[0]%3,-z[1]%3)
    stable=[o for o in joint if {twist(z) for z in o}==o]
    need(len(joint)==3 and len(stable)==1)
    rows=tuple(sum(1<<j for j in range(7) if (j-i)%7 in {1,2,4}) for i in range(7))
    tau=tuple(-i%7 for i in range(7)); pairlist=[(1,6),(2,5),(3,4)]
    data=[]
    for mask in range(8):
        rr=list(rows)
        for i,(u,v) in enumerate(pairlist):
            if mask>>i&1: rr[u]^=1<<v;rr[v]^=1<<u
        G=[p for p in permutations(range(7)) if is_aut(rr,p)]
        C=[g for g in G if comp(g,tau)==comp(tau,g)]
        h=folded_paths(rr,tau); H=len(literal_paths(rr))
        need(h==9);need(h%len(C)==0)
        data.append([mask,H,len(G),len(C),h,h//len(C)])
    need(data==[[0,189,21,3,9,3],[1,159,1,1,9,9],[2,159,1,1,9,9],
                [3,129,1,1,9,9],[4,159,1,1,9,9],[5,129,1,1,9,9],
                [6,129,1,1,9,9],[7,99,9,3,9,3]])
    return {'C3_marginal_classes_each':1,'C3_joint_classes':3,'C3_stable_joint_classes':1,
            'paley_pair_flip_columns':['mask','H','Aut','centralizer','fixed_paths','blue_orbits'],
            'paley_pair_flip_rows':data}

def run():
    result={'cyclic':cyclic_descent(),'paley7':paley_control(),'lifts':lift_controls(),'hostiles':joint_and_path_hostiles()}
    result['active_checks']=CHECKS
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':run()
