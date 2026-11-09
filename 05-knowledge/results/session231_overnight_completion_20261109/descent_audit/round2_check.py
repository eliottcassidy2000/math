#!/usr/bin/env python3
"""Independent exhaustive DFS path and SCC/topology audit. No producer imports."""
import importlib.util
import itertools as it
import json
from pathlib import Path
s=importlib.util.spec_from_file_location('own_round1',Path(__file__).with_name('check.py'))
M=importlib.util.module_from_spec(s);s.loader.exec_module(M)
need=M.need

def profile(a,root=0):
    n=len(a);full=(1<<n)-1
    totals=[0]*(full+1);starts=[0]*(full+1);ends=0
    totals[0]=1
    def walk(v,mask,first):
        nonlocal ends
        totals[mask]+=1
        starts[mask]+=first==root
        if mask==full and v==root:ends+=1
        for w in range(n):
            if not mask>>w&1 and a[v][w]:walk(w,mask|(1<<w),first)
    for v in range(n):walk(v,1<<v,v)
    return (totals[full],ends,sum(starts[m]*totals[full^m] for m in range(full+1)),totals[full^(1<<root)]),sum(totals[m]*totals[full^m] for m in range(full+1))

def induced(a,vs):return tuple(tuple(a[v][w] for w in vs) for v in vs)
def decomposition(a,w):
    # P contains precisely vertices that can reach w, including its entire SCC.
    P=[v for v in range(len(a)) if w in M.reach(a,v)]
    C=[v for v in range(len(a)) if v not in P]
    need(all(a[v][z] for v in P for z in C),'ordinal SCC split')
    return induced(a,P),P.index(w),induced(a,C)

def extremal(n,root):
    a=[[i<j for j in range(n)] for i in range(n)]
    a[0][-1]=False;a[-1][0]=True
    return induced(a,[root]+[i for i in range(n) if i!=root])

def root_iso(a,b):
    if len(a)!=len(b):return False
    for p in it.permutations(range(1,len(a))):
        perm=(0,)+p
        if all(a[i][j]==b[perm[i]][perm[j]] for i in range(len(a)) for j in range(len(a))):return True
    return False

def check(a):
    x,Q=profile(a);h,e,u,d=x;delta=u-h;r=sum(a[0])
    need(delta>=r*(r-1),'inherited degree lower bound')
    if r==1:
        w=next(i for i in range(len(a)) if a[0][i])-1
        B=induced(a,list(range(1,len(a))))
        (dd,q,R,_),_=profile(B,w)
        need((e,h,u,delta)==(dd-q,dd-q+R,dd+R,q),'sole exit identities')
    if r and delta==0:
        need(r==1,'zero defect has exactly one exit')
        P,wp,C=decomposition(B,w)
        need(len(C)>0,'zero endpoint implies nonempty suffix')
        (p,ep,z,dp),_=profile(P,wp)
        (c,ec,uc,dc),Qc=profile(C)
        need((h,e,u)==(p*c+z*Qc,p*c,p*c+z*Qc),'complete zero face decomposition')
        need(e%2==1 and 3*e<=h,'odd endpoint and sharp zero face')
        need((3*e==h)==(z==p and len(C)==1),'zero face equality iff')
        if e==1:
            normal,_=M.normalize(a)
            need(root_iso(normal,extremal(len(C)+2,len(C)+1)),'e1 full topology')
            need(h==2**len(C)+1,'e1 exact scalar')
    if r:
        L=2*h+delta-e
        need(2*L>=3*h+1,'old singleton floor')
        need((2*L==3*h+1)==(h==3 and e==1 and delta==0),'old floor equality iff')
        if 2*L==3*h+1:
            normal,_=M.normalize(a)
            need(len(normal)==3 and all(sum(row)==1 for row in normal),'old equality topology')
        if h>=5:
            need(2*L>=3*h+3,'new floor')
            eq=2*L==3*h+3
            expected=(h,e,u) in {(5,1,5),(9,3,9)} or (delta==1 and h>=5 and e==(h-1)//2 and (h-1)&(h-2)==0)
            need(eq==expected,'complete equality state classification')
            if eq:
                need(L%3==0 and L>3,'composite boundary')
                if delta==1:
                    normal,_=M.normalize(a)
                    need(root_iso(normal,extremal(len(normal),len(normal)-2)),'delta1 complete normalized topology')
                elif h==9:
                    need(len(C)==1 and p==z==3,'h9 exceptional physical interface')
        return x,L
    return x,None

def main():
    graphs=zero=one_exit=old=0;equal=set();bank=set();topology=0
    for n in range(1,7):
      for a in M.tournaments(n):
        x,L=check(a);h,e,u,d=x;graphs+=1
        if n<=5:need(x==M.profile(a),'DFS versus literal permutation profile')
        one_exit+=sum(a[0])==1
        zero+=bool(any(a[0])) and h==u
        bank.add(x[:3])
        if L is not None and 2*L==3*h+1:old+=1
        if L is not None and h>=5 and 2*L==3*h+3:equal.add((h,e,u,L));topology+=1
    transfer=trans_eq=0
    for x in bank:
      h,e,u=x
      for y in bank:
        k,f,v=y;N=M.value(x,y);L=h+u-e
        need(N-k*L==u*(v-k)+e*(k-f)>=0,'exact arbitrary-input transfer')
        if e<h and h>=5:
            transfer+=1;need(2*N>=(3*h+3)*k,'transferred sharpened floor')
            if 2*N==(3*h+3)*k:
                trans_eq+=1;need(f==v==k and N%3==0 and N>3,'sink-only composite equality')
                need((h,e,u,L) in equal,'transfer source classified')
    fam=[]
    for k in range(2,10):
        a=extremal(k+2,k);x,L=check(a)
        need(x[:3]==(2**k+1,2**(k-1),2**k+2),'explicit family state')
        fam.append(list(x[:3]) + [L])
    return dict(status='PASS',labelled_rooted_graphs=graphs,permutation_crosschecks=1099,sole_exit_cases=one_exit,zero_defect_nonsinks=zero,old_floor_topology_checks=old,improved_floor_topology_checks=topology,improved_floor_equalities=sorted(equal),physical_triples=len(bank),two_input_transfer_pairs=len(bank)**2,sharpened_transfer_pairs=transfer,transferred_equality_pairs=trans_eq,literal_family_parameters=fam,scope='independent finite exact controls; all-order classification audited separately')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
