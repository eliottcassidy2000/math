#!/usr/bin/env python3
"""Exact independent controls. Run from any working directory."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json, hashlib, sys
from audit_core import *
D=Path(__file__).resolve().parent
PACKET=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else D.parent

def words(n,I,opened,chord):
    a,b,x,y,c,z=0,1,2,n-3,n-2,n-1;W=list(range(x,c));H=list(range(3,c))
    if opened:
        if y not in I:return 'open_release', [[b]+W+[c,z,a],[z]+W+[c,b,a],[a,z]+W+[c,b]]
        if c not in I:return 'open_penultimate',[[a,b,c,z,x]+H,[b,c,a,z,x]+H,[b,x]+H+[c,a,z]]
        h=next(iter(I-{y,c}));t=y-int(h==y-1);L=list(range(x,t));R=list(range(t,c))
        return 'open_endrun',[[a,b,c,z]+L+R,[b,c,z]+L+[a]+R,[b]+L+[a]+R+[c,z]]
    if c in I:
        if y not in I:
            J={n-1-i for i in I};t=min(J-{1});L=list(range(1,t));R=list(range(t,c))
            return 'closed_last_movable',[[a]+L+[z]+R+[c],L+[z]+R+[a,c],L+R+[a,c,z]]
        h=next(iter(I-{y,c}))
        if chord and h not in {x,n-4}:return 'closed_endpair_guarded',[[a,x]+H+[c,z,b],list(reversed([a,b,x,c]+list(reversed(H))+[z])),[b,z]+H+[c,x,a]]
        return None
    if not chord:return 'closed_zero',[[c]+W+[b,z,a],[z,c]+W+[b,a],[a,z,b]+W+[c]]
    if x not in I:return 'closed_one',[[c,b,z,a]+W,[z,c,b,a]+W,[a,z,b]+W+[c]]
    if y not in I:t=min({n-1-i for i in I if n-1-i>x});name='closed_second_movable'
    elif 3 not in I:t=y;name='closed_second_fixed'
    else:t=n-4;name='closed_second_tail4'
    L=list(range(x,t));R=list(range(t,c))
    return name,[[a,b]+L+[z]+R+[c],[a]+L+[z]+R+[b,c],[a]+L+R+[b,c,z]]

def boundary(n,branch):
    A=raw_matrix(n);eq=constraints(n,A);a,b,x,y,c,z=0,1,2,n-3,n-2,n-1;W=list(range(x,c));H=list(range(3,c))
    if branch=='open':
        guards=[(b,c,1),(a,y,0),(a,c,1)];q1=[a,b,c,z,x]+H;q2=[b,c,a,z,x]+H
    else:
        guards=[(b,c,0),(x,c,1),(a,x,1),(a,c,1)];q1=[c,b,z,a]+W;q2=[z,c,b,a]+W
    eq += [A[u][v]^bit for u,v,bit in guards]
    eq += blue_equations(k_move(A,q1),q2)
    basis=eliminate(eq)
    if basis is None:return {'n':n,'branch':branch,'inconsistent':True}
    forms=[reduce_form(A[0][v],basis) for v in range(2,n-1)]
    groups={};cost=0;assignment=0
    for v,f in zip(range(2,n-1),forms):
        need((f>>1).bit_count()<=1,'non-component form after elimination')
        if f<2:cost+=f==0
        else:groups.setdefault(f&~1,[]).append((v,f&1))
    for free,terms in groups.items():
        zeros=sum(sign==0 for v,sign in terms);ones=len(terms)-zeros
        cost+=min(zeros,ones)
        if zeros>ones:assignment|=free
    for pivot,equation in sorted(basis.items()):
        if ((equation&1)^((equation&~(1<<pivot)&assignment).bit_count()&1)):assignment|=1<<pivot
    need(all(((e&1)^((e&assignment).bit_count()&1))==0 for e in eq),'minimizer fails equations')
    actual=sum(not (A[0][v]&assignment).bit_count()%2 for v in range(2,n-1))
    need(actual==cost and cost>3,'boundary minimum')
    return {'n':n,'branch':branch,'inconsistent':False,'minimum_middle_inports':cost,'minimizing_inports':[v for v in range(2,n-1) if (A[0][v]&assignment).bit_count()%2==0],'projected_free_components':[[v for v,s in ts] for ts in groups.values()]}

def chain_and_orbit_controls():
    totals=Counter()
    for n in range(7,261):
        F=[0]+list(range(2,n-2))+[1,n-2,n-1];r=reversal(F);g={i:r[n-1-i] for i in range(n)}
        formula={0:0,n-1:n-1,1:2,n-4:1,n-3:n-2,n-2:3};formula.update({i:i+2 for i in range(2,n-4)})
        need(g==formula,('printed four-port permutation',n))
        lengths=sorted(map(len,permutation_cycles(g)))
        need(lengths==sorted([1,1,n-2] if n%2 else [1,1,(n-2)//2,(n-2)//2]),'orbit lengths')
        middle=[k for k in lengths if k>1]
        possible={0}
        for k in middle: possible|={s+k for s in list(possible)}
        need(any(k%2==0 for k in middle) or 4 not in possible,'four-port obstruction')
        totals['four_port_orders']+=1
        if n<27:continue
        for branch in ['open','closed']:
            A=raw_matrix(n);H=list(range(3,n-2));W=list(range(2,n-2))
            if branch=='open':
                q1=[0,1,n-2,n-1,2]+H;q2=[1,n-2,0,n-1,2]+H
                chains=[[(0,4),(n-5,n-1),(0,n-6)],[(0,8),(n-9,n-1),(4,n-6),(5,n-5),(0,n-10)],[(0,12),(n-13,n-1),(8,n-6),(5,n-9),(4,n-10),(9,n-5),(0,n-14)]]
            else:
                q1=[n-2,1,n-1,0]+W;q2=[n-1,n-2,1,0]+W
                chains=[[(0,5),(n-6,n-1),(0,n-3)],[(0,9),(n-10,n-1),(5,n-3),(2,n-6),(0,n-7)],[(0,13),(n-14,n-1),(9,n-3),(2,n-10),(5,n-7),(6,n-6),(0,n-11)]]
            r=reversal(q2);E=path_edges(q1);ends=[]
            for chain in chains:
                parity=[]
                for j,((u,v),target) in enumerate(zip(chain,chain[1:])):
                    need(u<v and target[0]<target[1],'unordered printed chain pair')
                    if j%2==0:need(target==(n-1-v,n-1-u),'P chain link');parity.append(0)
                    else:
                        u1,v1=r[v],r[u];need(pair(u1,v1)==target,'Q2 chain link')
                        parity.append(int(u1>v1)^int((u,v) in E)^int(target in E))
                need(parity==[int(j==1) for j in range(len(parity))],('printed link parities',n,branch,parity))
                ends += [chain[0][1],chain[-1][1]];totals['chain_controls']+=1
            forced=n-3 if branch=='open' else 4;need(len(set(ends+[forced]))==7,'seven physical ports')
            if branch=='closed':
                S=k_move(A,q1);rr=reversal(q2)
                need((rr[1],rr[n-3])==(n-5,n-1),'forced extra port reflection')
    return dict(totals)

def main():
    text=(PACKET/'README.md').read_text()
    for formula in ['Q1=(a,b,L,z,R,c)','Q2=(a,L,z,R,b,c)','Q3=(a,L,R,b,c,z)','Q1=(a,b,c,z,L,R)','Q2=(b,c,z,L,a,R)','Q3=(b,L,a,R,c,z)','(a,x,D,c,z,b)','reverse(a,b,x,c,reverse(D),z)','(b,z,D,c,x,a)']:
        need(formula in text,('printed formula missing',formula))
    counts=Counter();full_orders=[15,16,17]
    selected_orders=list(range(18,41))+[47,60,81,104]
    for n in full_orders+selected_orders:
        if n in full_orders:ports=list(combinations(range(2,n-1),3))
        else:
            candidates=sorted({2,3,4,5,6,n//2,n-6,n-5,n-4,n-3,n-2});ports=list(combinations(candidates,3))
        for triple in ports:
            I=set(triple)
            for opened,chord in [(True,False),(True,True),(False,False),(False,True)]:
                result=words(n,I,opened,chord)
                if result is None:counts['unclaimed_cells']+=1;continue
                name,qs=result;A=symbolic_source(n,I,[(1,n-2,int(opened)),(2,n-2,int(chord))])
                try:replay_family(A,qs)
                except RuntimeError as e:raise RuntimeError((n,triple,opened,chord,name,str(e)))
                counts[name]+=1
    finite=[boundary(n,branch) for n in range(15,27) for branch in ('open','closed')]
    expected=json.loads((PACKET/'receipt-normal.json').read_text())
    for row,ref in zip(finite,expected['finite_boundary_affine_systems']):need(all(row[k]==v for k,v in ref.items()),'producer boundary disagreement')
    hostiles=[]
    for name,hostile in expected['hostile_controls'].items():
        n=hostile['n'];I=set(hostile['ports']);A=literal_rows(hostile['rows']);u,v=hostile['illegal_join']
        if name=='old_endpair_insertion':
            oldq1=[0,1,n-2,n-1]+list(range(2,n-2));oldq2=[1,n-2,n-1]+list(range(2,n-3))+[0,n-3]
            need(all(A[u][v] for u,v in zip(oldq1,oldq1[1:])),'old first word')
            need(k_move(A,oldq1)[u][v]==0 and (u,v) in list(zip(oldq2,oldq2[1:])),'wrong illegal join')
        else:
            oldq1=list(range(n-4))+[n-1,n-4,n-3,n-2]
            need(A[u][v]==0 and (u,v) in list(zip(oldq1,oldq1[1:])),'wrong fixed cut hostile')
        label,qs=words(n,I,bool(A[1][n-2]),bool(A[2][n-2]));flags=check_literal(A,I,qs,[False]*3)
        hostiles.append({'name':name,'n':n,'ports':sorted(I),'repair':label,'blue_steps':flags})
    receipt={'status':'PASS','method':'independent GF(2) row elimination; no producer imports','source_cells':dict(counts),'all_port_orders':full_orders,'strategic_port_orders':selected_orders,'finite_boundary':finite,'hostile_replays':hostiles,'chain_and_orbit_controls':chain_and_orbit_controls(),'small_order_scope':'Only the two stored fixed-template failures are replayed here; no full n<15 theorem is claimed.'}
    print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
