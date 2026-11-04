"""Exact joint endpoint entry and its unique-path threshold grammar."""
from collections import deque,Counter
from itertools import permutations,combinations,product
from math import factorial
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
def need(p,msg):
    if not p:raise ValueError(msg)
def tables(T):
    n=len(T);f=(1<<n)-1;E=[[0]*n for _ in range(f+1)];H=[0]*(f+1);H[0]=1
    for i in range(n):E[1<<i][i]=1
    for s in range(1,f+1):
        for i,a in enumerate(E[s]):
            if not a:continue
            H[s]+=a;todo=T[i]&~s
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;E[s|bit][j]+=a
    J=[[0]*n for _ in range(n)];starts=[]
    for v in range(n):
        D=[[0]*n for _ in range(f+1)];D[1<<v][v]=1
        for s in range(1,f+1):
            for i,a in enumerate(D[s]):
                if not a:continue
                J[v][i]+=a*H[f^s];todo=T[i]&~s
                while todo:
                    bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;D[s|bit][j]+=a
        starts.append(sum(D[f]))
    return H,J,starts
def transitive_order(T,V):
    V=sorted(V,key=lambda i:-sum(T[i]>>j&1 for j in V))
    return V if all(T[a]>>b&1 for k,a in enumerate(V) for b in V[k+1:]) else None
def grammar(T,v,w):
    if v==w:
        A=transitive_order(T,[i for i in range(len(T)) if i!=w])
        return None if A is None else dict(path=[w],tail=A,thresholds=None)
    paths=deque([[v]]);seen={v};P=None
    while paths:
        p=paths.popleft()
        if p[-1]==w:P=p;break
        for i in range(len(T)):
            if T[p[-1]]>>i&1 and i not in seen:seen.add(i);paths.append(p+[i])
    if P is None:return None
    if any(bool(T[a]>>b&1)!=(j==i+1) for i,a in enumerate(P) for j,b in enumerate(P) if i<j):return None
    A=transitive_order(T,[i for i in range(len(T)) if i not in P])
    if A is None:return None
    ks=[]
    for a in A:
        bits=[T[a]>>p&1 for p in P];k=sum(bits)-1
        if bits!=[int(i<=k) for i in range(len(P))]:return None
        ks.append(k)
    if any(b-a>=2 for i,a in enumerate(ks) for b in ks[i+1:]):return None
    return dict(path=P,tail=A,thresholds=ks)
def main():
    atlas=json.loads((HERE.parent/'strong_atom_continuation_vitali_20261026/strong-atlas.json').read_text());checks=Counter();entries=[];root_checks=0;two_start=[];two_start_checks=0
    for z in atlas:
        T=z['source'];n=len(T);H,J,S=tables(T);f=(1<<n)-1;need(H[f]==z['H'],'inherited atlas H')
        for v in range(n):
            for w in range(n):
                G=grammar(T,v,w);need((J[v][w]==1)==(G is not None),'exact threshold recognition');checks[n]+=1
                if J[v][w]==1:
                    row=dict(n=n,H=H[f],rows=T,ports=[v,w],start_w=S[w],q=H[f^(1<<w)],grammar=G);entries.append(row)
        for w in range(n):
            if S[w]==1:
                need(any(J[v][w]==1 for v in range(n))==(H[f^(1<<w)]==1),'unique-start reentry iff transitive deletion')
                root_checks+=1
            if S[w]==2:
                accepted=[v for v in range(n) if J[v][w]==1]
                need(bool(accepted)==(n==4),'two-start reentry only in the strong four-vertex source')
                two_start_checks+=1
                if accepted:two_start.append(dict(rows=T,mark=w,ports=accepted,H=H[f],F=sum(row[w] for row in J)))
    labelled=Counter();labelled_unit=Counter();literal_checks=0
    for n in range(1,6):
        pairs=list(combinations(range(n),2))
        for code in range(1<<len(pairs)):
            T=[0]*n
            for bit,(a,b) in enumerate(pairs):
                u,v=(a,b) if code>>bit&1 else (b,a);T[u]|=1<<v
            H,J,S=tables(T)
            for v in range(n):
                for w in range(n):
                    need((J[v][w]==1)==(grammar(T,v,w) is not None),'all-labelled threshold grammar');labelled[n]+=1
                    labelled_unit[n]+=int(J[v][w]==1)
            if n<=4:
                direct=[[0]*n for _ in range(n)]
                for p in permutations(range(n)):
                    for k in range(1,n+1):
                        if all(T[a]>>b&1 for a,b in zip(p[:k],p[1:k])) and all(T[a]>>b&1 for a,b in zip(p[k:],p[k+1:])):
                            direct[p[0]][p[k-1]]+=1
                need(direct==J,'independent permutation/cut joint count');literal_checks+=n*n
    sequence=[0,1]
    for n in range(2,21):sequence.append(3*sequence[-1]-sequence[-2])
    for n in range(1,6):need(labelled_unit[n]==factorial(n)*sequence[n],'rigid bipointed class enumeration')
    threshold_counts=[]
    for L in range(2,7):
        for m in range(5):
            direct=sum(all(b-a<=1 for i,a in enumerate(word) for b in word[i+1:]) for word in product(range(L+1),repeat=m))
            # Coefficients of (1-x)^(L-1)/(1-2x)^L, exact binomial arithmetic.
            from math import comb
            formula=sum((-1)**j*comb(L-1,j)*comb(L+m-j-1,m-j)*2**(m-j) for j in range(min(m,L-1)+1))
            need(direct==formula,'literal threshold language versus GF');threshold_counts.append([L,m,direct])
    out=dict(status='FINITE_EXACT',class_count=len(atlas),strong_class_checks=dict(checks),entry_count=len(entries),
             unique_start_roots_checked=root_checks,all_labelled_checks=dict(labelled),all_labelled_unit_entries=dict(labelled_unit),
             two_start_roots_checked=two_start_checks,two_start_reentry_sources=two_start,
             class_sequence=sequence[1:],threshold_counts=threshold_counts,independent_literal_checks=literal_checks,entries=entries)
    (HERE/'entry.json').write_text(json.dumps(out,indent=2)+'\n');print('strong checks',dict(checks),'entries',len(entries),'unique-start roots',root_checks,'labelled',dict(labelled),'literal',literal_checks)
if __name__=='__main__':main()
