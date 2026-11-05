"""Standalone exact verification of the balanced-source normal form."""
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import hashlib
import json
from core import inventories, constants, gadget, glue, adjacency_map, cycles

HERE=Path(__file__).resolve().parent
CAP=50001
GAPS=(7,11,13,21,23,41,43,47,55,67,77,83,149,163,173,231,233,247,251,
277,299,317,391,421,433,443,449,503,519,533,569,583,617,637,821,
911,1067,1339,1861,1903,2647,3187)
P={
'N':((6,4,10,9),(4,3,7,6),(10,7,15,14),(9,6,14,13)),
'M':((11,8,20,18),(8,6,14,13),(20,14,34,31),(18,13,31,29)),
'B':((14,10,26,23),(11,8,20,18),(20,15,35,32),(19,14,34,31)),
}

def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def new(s):
    D,E,U,d=s
    return D-E-d,U-D+E-d,E,d-E


def old(x):
    a,f,e,p=x
    return a+2*e+p,e,a+f+2*e+2*p,e+p


def successors(x):
    a,f,e,p=x
    return ((a+e,a+2*f+e+2*p,0,e+p),
            (2*a+f+e+p,f,e+p,0),
            (f+p,0,a+e+p,e))


def append(T,c,label):
    n=len(T);out=list(T)+[0]
    for v in range(n):
        if label in ('J','Z') and v==c:out[n]|=1<<v
        else:out[v]|=1<<n
    return out,n if label=='J' else c


def profile(T,c):
    q=inventories(T)
    return q['H'],q['e'][c],q['R'][c],q['allH'][((1<<len(T))-1)^(1<<c)]


def tournaments(n):
    pairs=[(i,j)for i in range(n)for j in range(i+1,n)]
    for mask in range(1<<len(pairs)):
        T=[0]*n
        for k,(a,b)in enumerate(pairs):
            if mask>>k&1:T[a]|=1<<b
            else:T[b]|=1<<a
        yield T


def induced(T,vertices):
    return [sum(((T[u]>>v)&1)<<j for j,v in enumerate(vertices))for u in vertices]


def canonical_inverse(T,r):
    require(len(T)>1,'singleton not in balanced domain')
    q=inventories(T)
    require(q['s'][r]==1 and q['R'][r]==q['H'],'promised balanced source')
    cands=[v for v in range(len(T))if T[r]>>v&1]
    require(len(cands)==1,'unique root out-neighbor')
    c=cands[0];mask=((1<<len(T))-1)^(1<<r)
    sinks=[v for v in range(len(T))if v!=r and not T[v]&mask]
    require(len(sinks)==1 and sinks[0]!=c,'canonical deleted sink')
    z=sinks[0];vertices=[v for v in range(len(T))if v not in(r,z)]
    W=induced(T,vertices);root=vertices.index(c)
    require(inventories(W)['s'][root]==1,'recovered physical source')
    X,u=append(*append(W,root,'O'),'J')
    require(adjacency_map(X,T,vertices+[z,r]),'canonical reconstruction every edge')
    require(u==len(X)-1,'root reconstruction')
    D,E,U,d=profile(W,root)
    require(profile(T,r)==(D+2*U,D,D+2*U,D),'inverse response')
    return W,root


def literal(T,c):
    n=len(T);full=(1<<n)-1;paths={}
    for M in range(1<<n):
        paths[M]=[p for p in permutations([v for v in range(n)if M>>v&1])
                  if all(T[a]>>b&1 for a,b in zip(p,p[1:]))]
    H=len(paths[full]);E=sum(p[-1]==c for p in paths[full])
    S=sum(p[0]==c for p in paths[full])
    U=sum(len(paths[M])*sum(bool(p)and p[0]==c for p in paths[full^M])for M in paths)
    return (H,E,U,len(paths[full^(1<<c)])),S


def complete_four(cap):
    base=(0,0,0,1);seen={base};queue=[base]
    projected={(1,1,1)};values={3}
    for x in queue:
        D,E,U,d=old(x);G=D+2*U
        require(G<=cap,'four-coordinate cap')
        projected.add((D,E,U));values.add(G)
        for z in successors(x):
            A,F,E,P=z;H=3*A+2*F+6*E+5*P
            require(2*H>=3*G,'one-letter G growth')
            if H<=cap and z not in seen:seen.add(z);queue.append(z)
    return len(seen)+1,projected,values


@lru_cache(None)
def coeff(L):return constants(L)


def complete_three(cap):
    base=(1,1,1);seen={base};queue=[base];values=set()
    for D,E,U in queue:
        G=D+2*U;values.add(G)
        cand=[(D,0,2*U)]
        for L in range(1,cap//(3*U)+1):
            g,Q,B=coeff(L)
            require(B>=L and g*D-(L%2)*E>=0,'indexed bound')
            cand.append((g*D-(L%2)*E+L*U,D-int(L==1)*E,Q*D+B*U))
        for z in cand:
            H=z[0]+2*z[2]
            require(H>G,'indexed strict G growth')
            if H<=cap and z not in seen:seen.add(z);queue.append(z)
    return seen,values


def indexed_recipe(recipe):
    s=(1,1,1,1);T=[0];c=0
    for kind,L in recipe:
        D,E,U,d=s
        if kind=='O':
            s=(D,0,2*U,d);T,c=append(T,c,'O')
        else:
            require(kind=='M' and isinstance(L,int)and L>=1,'valid source recipe')
            g,Q,B=coeff(L)
            s=(g*D-(L%2)*E+L*U,D-int(L==1)*E,Q*D+B*U,D)
            T,c=glue(T,[c],gadget(L))
    require(inventories(T)['s'][c]==1,'repair physical source')
    require(profile(T,c)==s,'repair literal graph profile')
    return s,T,c


def main():
    receipt={'status':'PASS','cap':CAP}
    count4,p4,v4=complete_four(CAP)
    p3,v3=complete_three(CAP)
    require(p4==p3,'entire projected physical source banks agree')
    require(v4==v3,'entire scalar banks agree')
    gaps=tuple(n for n in range(3,CAP+1,2)if n not in v4)
    require(gaps==GAPS,'exact balanced gaps')
    require(count4==908226 and len(p3)==788690 and len(v4)==24958,'frozen counts')
    receipt.update(source_profiles_including_singleton=count4,
        projected_source_profiles=len(p3),balanced_values=len(v4),balanced_gaps=gaps)
    del p4,p3,v4,v3
    graph_sources=balanced=literal_controls=0
    for n in range(1,6):
        for T in tournaments(n):
            inv=inventories(T)
            for c in range(n):
                if inv['s'][c]!=1:continue
                graph_sources+=1;s=profile(T,c);D,E,U,d=s
                require(literal(T,c)==(s,1),'literal path and cover counts')
                literal_controls+=1
                X,r=append(*append(T,c,'O'),'J')
                require(profile(X,r)==(D+2*U,D,D+2*U,D),'forward balanced profile')
                W,b=canonical_inverse(X,r)
                require(W==T and b==c,'left inverse on every labelled source')
                require(cycles(X)==1,'forward graph C1')
                if n>1 and D==U:
                    balanced+=1;canonical_inverse(T,c)
                if n>1:
                    a,f,e,p=new(s)
                    if f==0:require(T[c].bit_count()==1,'zero f root degree lemma')
                    if E==0:require(any(row==0 for row in T),'weak source universal sink')
    receipt.update(labelled_source_controls_through_order_5=graph_sources,
        labelled_balanced_controls_through_order_5=balanced,literal_controls=literal_controls)
    # All-order rebase, checked on every source graph through order four and k<=4.
    rebase_controls=0
    for n in range(1,5):
        for W in tournaments(n):
            inv=inventories(W)
            for c in range(n):
                if inv['s'][c]!=1:continue
                D,E,U,d=profile(W,c)
                for k in range(5):
                    A,r=append(W,c,'J')
                    for _ in range(k):A,r=append(A,r,'Z')
                    A,last=append(A,r,'J')
                    B,b=list(W),c
                    for _ in range(k+1):B,b=append(B,b,'O')
                    B,b=append(B,b,'J')
                    mapping=list(range(n))+list(range(n+1,n+k+2))+[n]
                    require(adjacency_map(B,A,mapping),'all-order graph rebase')
                    require(profile(A,r)==(D+2**(k+1)*U,D,D+2**(k+1)*U,D),
                            'balanced reroot response')
                    require(inventories(A)['H']==D+2**(k+1)*U,'rebase scalar')
                    rebase_controls+=1
    receipt['labelled_rebases_k_0_through_4']=rebase_controls
    residual_controls=0
    for n in range(1,5):
        for W in tournaments(n):
            inv=inventories(W)
            for c in range(n):
                if inv['s'][c]!=1:continue
                D,E,U,d=profile(W,c)
                for k in range(1,5):
                    T,r=append(W,c,'O')
                    for _ in range(k):T,r=append(T,r,'Z')
                    T,r=append(T,r,'J');q=inventories(T)
                    value=D+2*(2**(k+1)-1)*U-(2**(k+1)-2)*d
                    require(q['H']==value and cycles(T)==1,'residual scalar and C1')
                    balanced_roots=[a for a in range(len(T))if q['s'][a]==1 and q['R'][a]==q['H']]
                    if n==1:
                        require(c in balanced_roots,'singleton residual is balanced')
                        canonical_inverse(T,c)
                    else:
                        require(not balanced_roots,'nontrivial residual has no balanced root')
                        require([a for a,row in enumerate(T)if row.bit_count()==1]==[r],
                                'residual unique degree-one vertex')
                        mask=((1<<len(T))-1)^(1<<r)
                        require(all(T[a]&mask for a in range(len(T))if a!=r),'no residual deletion sink')
                    residual_controls+=1
    receipt['labelled_residual_controls_k_1_through_4']=residual_controls
    # Small return menu is not complete; a physical balanced 19 is a hostile control.
    menu={(3,1)};queue=list(menu)
    for h,e in queue:
        for z in ((2*h-e,e),(6*h-e,2*h-e),(3*h,h)):
            if z[0]<=19 and z not in menu:menu.add(z);queue.append(z)
    require(19 not in {h for h,e in menu},'restricted balanced menu omits 19')
    T,c=[2,0],0
    for letter in 'OZ':T,c=append(T,c,letter)
    require(profile(T,c)==(5,1,7,1),'physical precursor for balanced 19')
    T,c=append(*append(T,c,'O'),'J')
    require(profile(T,c)==(19,5,19,5),'physical balanced 19')
    repairs=json.loads((HERE/'REPAIRS.json').read_text())
    require(set(map(int,repairs))==set(GAPS)-{7,21},'every finite exception repaired')
    repair_cycles={1:0,2:0,3:0}
    for target,entry in repairs.items():
        r=entry['record'];kind=r['kind']
        if kind=='C1':
            s,T,c=indexed_recipe(r['source_recipe'])
            require(tuple(r['source'])==s,'C1 repair response')
            D,E,U,d=s;value=D-E+U
            X,b=append(T,c,'J')
            require(inventories(X)['H']==value and cycles(X)==1,'literal C1 repair')
            require(entry['cycles']==1,'C1 label')
        else:
            left,_,_=indexed_recipe(r['left_recipe']);right,_,_=indexed_recipe(r['right_recipe'])
            require(tuple(r['left'])==left and tuple(r['right'])==right,'pair repair responses')
            x,y=new(left),new(right)
            require(min(x+y)>=0,'nontrivial physical pair sources')
            value=sum(x[i]*P[kind][i][j]*y[j]for i in range(4)for j in range(4))
            require(entry['cycles']==(2 if kind in ('N','M')else 3),'pair family cycle count')
        require(value==int(target),'finite gap repair scalar')
        repair_cycles[entry['cycles']]+=1
    receipt['finite_gap_repairs_by_cycle_count']=repair_cycles
    receipt['files_sha256']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        for name in ('PROOF.md','verify.py','core.py','REPAIRS.json')}
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=='__main__':main()
