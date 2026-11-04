"""Small exact closure probe, complete labelled strong B through order five."""
from pathlib import Path
from itertools import combinations,permutations
from functools import lru_cache
import json
from validate import need,strong
HERE=Path(__file__).resolve().parent
@lru_cache(None)
def profile(T):
    n=len(T);f=(1<<n)-1;e=[[0]*n for _ in range(f+1)];h=[0]*(f+1);h[0]=1
    for v in range(n):e[1<<v][v]=1
    for m in range(1,f+1):
        for v in range(n):
            c=e[m][v]
            if c:
                h[m]+=c;q=T[v]&~m
                while q:
                    b=q&-q;q-=b;e[m|b][b.bit_length()-1]+=c
    F=[sum(e[m][v]*h[f^m] for m in range(f+1)) for v in range(n)]
    return h[f],tuple(F)
def substitute(T,w):
    n=len(T);U=list(T)+[T[w]|(1<<(n+1)),T[w]|(1<<w)]
    for a in range(n):
        if T[a]>>w&1:U[a]|=(1<<n)|(1<<(n+1))
    U[w]|=1<<n
    return tuple(U)
def clone(T,w):
    n=len(T);U=list(T)+[T[w]]
    for a in range(n):
        if T[a]>>w&1:U[a]|=1<<n
    U[w]|=1<<n
    return tuple(U)
def brute(T,w):
    h=F=0;n=len(T)
    for p in permutations(range(n)):
        good=[bool(T[p[i]]>>p[i+1]&1) for i in range(n-1)]
        h+=all(good)
        # First ordered slot ends at w; it is nonempty. Second may be empty.
        k=p.index(w)+1
        F+=all(good[:k-1]) and all(good[k:])
    return h,F
def run_partition(U,w):
    n=len(U);inside={w,n-2,n-1};counts=[0,0,0]
    for p in permutations(range(n)):
        if not all(U[a]>>b&1 for a,b in zip(p,p[1:])):continue
        runs=sum(p[j] in inside and (j==0 or p[j-1] not in inside) for j in range(n))
        counts[runs-1]+=1
    need(counts[0]%3==counts[1]%6==counts[2]%6==0,'C3 run decoration factors')
    return dict(path_counts_by_run_number=counts,quotient_word_coefficients=[counts[0]//3,counts[1]//6,counts[2]//6],C3_ordered_cover_profile=[3,6,6])
def main():
    seen={};witness=None;counts={};aliases=0;split=0
    for n in range(1,6):
        edges=list(combinations(range(n),2));ct=0
        for mask in range(1<<len(edges)):
            T=[0]*n
            for j,(a,b) in enumerate(edges):
                if mask>>j&1:T[a]|=1<<b
                else:T[b]|=1<<a
            T=tuple(T)
            if not strong(T):continue
            ct+=1;h,Fs=profile(T)
            for w,F in enumerate(Fs):
                U=substitute(T,w);hh,FF=profile(U);state=(hh,FF[w]);key=(n,h,F)
                if key in seen:
                    aliases+=1;old=seen[key]
                    if old['after']!=state:
                        split+=1
                        if witness is None:witness=dict(n=n,before=[h,F],left=old,right=dict(rows=T,root=w,after_rows=U,after=state))
                else:seen[key]=dict(rows=T,root=w,after_rows=U,after=state)
        counts[n]=ct
    need(witness is not None,'marked two-coordinate quotient has a hostile substitution')
    for side in ('left','right'):
        z=witness[side]
        need(brute(z['rows'],z['root'])==tuple(witness['before']),'independent permutation/cut before')
        need(brute(z['after_rows'],z['root'])==tuple(z['after']),'independent permutation/cut after')
        need(strong(z['after_rows']),'strong substitution')
        z['run_partition']=run_partition(z['after_rows'],z['root'])
        V=clone(z['rows'],z['root']);hh,ff=profile(V)
        need(strong(V) and brute(V,z['root'])==(hh,ff[z['root']]),'independent transitive-clone control')
        need(hh==witness['before'][0]+2*z['run_partition']['quotient_word_coefficients'][1],'first missing run coordinate already seen by TT2 clone')
        z['TT2_clone']=dict(rows=V,after=[hh,ff[z['root']]])
    result=dict(status='FINITE_EXACT_COUNTEREXAMPLE_WITH_INDEPENDENT_PERMUTATION_CUT_REPLAY',universe='Every labelled strong tournament B of order1 through5, every marked vertex; C3 substitution at that mark, keeping old name as first clone',strong_counts=counts,marked_aliases=aliases,split_aliases=split,witness=witness,
                consequence='No operation on (H(B),E_w(B,2)) alone can compute this marked C3 substitution on all strong B. This does not refute other operation-closed carriers or the gap conjecture.')
    (HERE/'substitution_probe.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: all labelled strong orders1..5; same-order marked (h,F) collision splits under C3 substitution; both witnesses independently replayed')
    print(json.dumps(witness,sort_keys=True))
if __name__=='__main__':main()
