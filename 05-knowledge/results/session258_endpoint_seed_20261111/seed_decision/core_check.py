"""Independent literal cover enumeration for strong-core seed classification."""
from itertools import permutations
from functools import lru_cache
from pathlib import Path
import json
import check as p

def require(x,m):
    if not x:raise RuntimeError(m)

@lru_cache(None)
def counts(T):
    n=len(T);e=[0]*n;s=[0]*n;K=[[0]*n for _ in T]
    for order in permutations(range(n)):
        good=[bool(T[order[j]]>>order[j+1]&1) for j in range(n-1)]
        if all(good):e[order[-1]]+=1;s[order[0]]+=1
        for cut in range(1,n):
            if all(good[j] for j in range(n-1) if j!=cut-1):K[order[cut-1]][order[cut]]+=1
    return sum(e),e,s,K

def stats(T,O):
    h,e,s,K=counts(T);selected=[i for i in range(len(T)) if O>>i&1]
    return h,sum(K[i][j] for i in selected for j in selected),sum(s[j]+sum(K[i][j] for i in range(len(T))) for j in selected)

def components(T):
    n=len(T);reach=[T[i]|(1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    remain=set(range(n));blocks=[]
    while remain:
        a=min(remain);C={b for b in remain if reach[a]>>b&1 and reach[b]>>a&1};remain-=C;blocks.append(sorted(C))
    return sorted(blocks,key=lambda C:-sum(len(D) for D in blocks if D!=C and T[C[0]]>>D[0]&1))

def induced(T,C):return tuple(sum(1<<j for j,b in enumerate(C) if T[a]>>b&1) for a in C)

def main():
    bases=sets=cross=seeds=0;core_seeds=[];witnesses={}
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            T=p.tournament(n,code);h,e,s,K=counts(T);Cs=components(T);which={v:i for i,C in enumerate(Cs) for v in C};bases+=1
            for a in range(n):
                for b in range(a+1,n):
                    if which[a]!=which[b]:require(K[a][b]+K[b][a]>h,('cross_bound',T,a,b));cross+=1
            for O in range(1<<n):
                hh,w,u=stats(T,O);sets+=1
                require(u>=O.bit_count()*h,('split_bound',T,O))
                if h-w!=1 or not u:continue
                seeds+=1;selected=[i for i in range(n) if O>>i&1]
                require(len({which[i] for i in selected})==1,'one_selected_component')
                z=which[selected[0]];C=Cs[z]
                require(all(len(D)==1 for j,D in enumerate(Cs) if j!=z),'singleton_exterior')
                S=induced(T,C);Ocore=sum(1<<i for i,v in enumerate(C) if O>>v&1);hs,ws,us=stats(S,Ocore)
                require(hs-ws==1 and u==(1<<(len(Cs)-1-z))*us,'padding_law')
                if len(Cs)==1 and n>1:
                    core_seeds.append((n,h,u,O.bit_count()))
                    witnesses.setdefault(u,{'rows':list(T),'selected_mask':O,'H':h,'W':w,'U':u})
    receipt=json.loads((Path(__file__).parent/'receipt.json').read_text());classified=[]
    for q in range(1,30):
        trivial=q&(q-1)==0
        qualifying=[(n,h,u,k) for n,h,u,k in core_seeds if h<=q//2 and q%u==0 and (q//u)&((q//u)-1)==0]
        exists=trivial or bool(qualifying)
        if q<=13:require(exists==(receipt['results'][q-1]['status']=='SEED'),('core_classifier',q))
        # For this range, every permitted strong core is among orders<=5.
        require(p.bound(q//2)<=5,'finite_core_bound')
        classified.append({'q':q,'status':'SEED' if exists else 'NO_SEED'})
    require([r['q'] for r in classified if r['status']=='SEED']==[1,2,4,6,8,11,12,16,20,22,24,29],'q29_spectrum')
    examples=[]
    for slope,q in ((5,20),(29,29)):
        witness=witnesses[q].copy();B=tuple(witness['rows']);O=witness['selected_mask'];tail_counts=[]
        for t in (1,2):
            T=B
            for _ in range(t):T=p.join(T,(0,))
            out=list(T);v=len(T)
            for a in range(v):
                if not(O>>a&1):out[a]|=1<<v
            out.append(O);h=counts(tuple(out))[0]
            require(h==1+(1<<t)*q,'literal_tail_witness')
            tail_counts.append(h)
        witness.update(slope=slope,tail_counts=tail_counts);examples.append(witness)
    print(json.dumps({'status':'PASS','bases':bases,'selected_sets':sets,'cross_component_pairs':cross,'positive_unit_seeds':seeds,'strong_unit_core_profiles':sorted(set(core_seeds)),'q1_to_29_core_classification':classified,'examples':examples,'scope':'all labelled bases orders1..5; all-order fixed-q1..29 via cited core floor'},indent=2))
if __name__=='__main__':main()
