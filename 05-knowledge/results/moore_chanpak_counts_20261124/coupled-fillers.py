"""Check asymmetric fillers on the coupled pentagon comparator."""
from itertools import combinations, permutations, product
from bounded_template_profile import profile

def add(T,m):
    n=len(T)
    return tuple(row|((1-((m>>i)&1))<<n) for i,row in enumerate(T))+(m,)

def canon(T):
    # Relabel preserving score classes. Permutation signature distinguishes
    # ordinary isomorphism, not just graph-converse pairs.
    classes={}
    for u,r in enumerate(T):classes.setdefault(r.bit_count(),[]).append(u)
    cells=[classes[q] for q in sorted(classes)]
    n=len(T); E=list(combinations(range(n),2))
    return min(sum((T[p[i]]>>p[j]&1)<<e for e,(i,j) in enumerate(E))
               for parts in product(*(tuple(permutations(c)) for c in cells))
               for p in [sum(parts,())])

def fixtures(nmax=7):
    L={(0,)}
    by={}
    sizes=[]
    for n in range(1,nmax+1):
        sizes.append(len(L))
        for T in sorted(L):
            scores=tuple(profile(T,k) for k in (1,2,3))
            by.setdefault(scores,T)
        if n==nmax:break
        R={}
        for T in sorted(L):
            for m in range(1<<n):
                Q=add(T,m)
                R.setdefault(canon(Q),Q)
        L=set(R.values())
    assert sizes==[1,1,2,4,12,56,456][:nmax]
    return by

def delta(A,B):
    h,q,r=A; H,Q,R=B
    return 2*h*H+5*(q*H-h*Q)+4*(r*Q-q*R)

if __name__=='__main__':
    F=fixtures()
    both=[]; best=None
    for A,TA in F.items():
        for B,TB in F.items():
            if A[0]==B[0]==1:continue
            d=delta(A,B)
            if best is None or abs(d)<abs(best[0]):best=(d,A,B)
            if abs(d)==2:both.append((d,A,B,TA,TB))
    crossing=sum(q*H>Q*h and r*Q<R*q for h,q,r in F for H,Q,R in F)
    assert len(F)==195 and not both and crossing==1365
    print('profile pairs with opposing normalized ratio signs',crossing)
    print('distinct profiles through order seven',len(F),'closest',best,'nontransitive hits',both)
