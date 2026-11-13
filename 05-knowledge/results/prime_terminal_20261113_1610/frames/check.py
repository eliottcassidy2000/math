"""Independent literal controls and audit of the pinned order-nine search."""
from collections import Counter
from itertools import permutations
from math import comb, isqrt
from pathlib import Path
import gzip
import hashlib
import json

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
Q=(6,12,8,1)
SOURCE=((16,5,9,35,46,7),0)
SHIFTED=((258,4,265,51,39,135,191,287,122),1)
BEZOUT=(((2,4,25,19,3),1),((2,12,73,81,103,79,3),5),((2,12,9,17,103,79,15),5))

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def tt(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def valid(T):
    n=len(T)
    return all(0<=T[i]<(1<<n) and not (T[i]>>i&1) for i in range(n)) and all(
        (T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(n) for j in range(i))

def hp(T):
    # End-vertex deletion recurrence, independent of the C++ forward DP.
    n=len(T); full=(1<<n)-1
    need(valid(T),'tournament validity including every internal module arc')
    incoming=[full^(1<<i)^T[i] for i in range(n)]
    d=[[0]*n for _ in range(1<<n)]
    for s in range(1,1<<n):
        z=s
        while z:
            bit=z&-z;v=bit.bit_length()-1;z-=bit
            rest=s^bit
            if not rest:d[s][v]=1;continue
            pred=rest&incoming[v]
            while pred:
                b=pred&-pred;pred-=b;d[s][v]+=d[rest][b.bit_length()-1]
    return sum(d[-1])

def strong(T):
    full=(1<<len(T))-1
    for r in range(len(T)):
        old=-1;seen=1<<r
        while old!=seen:
            old=seen
            for i in range(len(T)):
                if seen>>i&1:seen|=T[i]
        if seen!=full:return False
    return True

def inflate(T,r,M):
    labels=[(i,j) for i in range(len(T)) for j in range(len(M) if i==r else 1)]
    return tuple(sum(1<<b for b,(w,k) in enumerate(labels) if
        (v!=w and T[v]>>w&1) or (v==w==r and M[j]>>k&1)) for v,j in labels)

def two_run(T,r):
    edges=[(o,i) for o in range(len(T)) for i in range(len(T))
           if T[r]>>o&1 and T[o]>>i&1 and T[i]>>r&1]
    return not any(o!=v and i!=j for o,i in edges for v,j in edges)

def skeleton(T,r):
    """Exterior permutations and forced/optional r insertion gaps."""
    out=Counter()
    for word in permutations(v for v in range(len(T)) if v!=r):
        forced=0
        optional=int(bool(T[r]>>word[0]&1))+int(bool(T[word[-1]]>>r&1))
        for x,y in zip(word,word[1:]):
            direct=T[x]>>y&1;via=(T[x]>>r&1) and (T[r]>>y&1)
            if not direct:
                if not via:break
                forced+=1
            elif via:optional+=1
        else:
            for k in range(optional+1):out[forced+k]+=comb(optional,k)
    return dict(sorted((j,a) for j,a in out.items() if j))

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def main():
    out={'status':'FINITE-EXACT controls for the proved response statements'}
    need(hp(Q)==5,'Q count')
    framed=[]
    for item,want in zip(BEZOUT,((9,7),(45,17),(41,18))):
        T,r=item
        coeff=skeleton(T,r)
        need(coeff=={1:want[0],2:want[1]},'complete Bezout frame response')
        framed.append({'rows':T,'root':r,'ordered_runs':coeff,'Q_inflation_H':hp(inflate(T,r,Q))})
    vals=[x['Q_inflation_H'] for x in framed]
    need(vals==[199,599,601] and vals[1:]==[3*vals[0]+2,3*vals[0]+4],'literal both-digit Bezout success')
    need(all(prime(h) for h in vals),'named prime controls')
    out['bezout']=framed
    for (T,r),want in ((SOURCE,{1:37,2:27}),(SHIFTED,{1:599,2:324})):
        need(strong(T) and two_run(T,r),'strong two-run named frame')
        need(skeleton(T,r)==want,'independent all-run skeleton count')
    controls=[]
    for k in range(1,6):
        s=hp(inflate(*SOURCE,tt(k)))
        need(s==27*2**k-17,'source whole count')
        row={'k':k,'source_H':s,'source_order':k+5}
        if k>=3:
            y=hp(inflate(*SHIFTED,tt(k-2)))
            need(y==81*2**k-49==3*s+2,'shifted same-target formula')
            row.update(target_H=y,target_order=k+6)
        controls.append(row)
    out['shifted_controls']=controls
    search=json.loads((HERE/'search_receipt.json').read_text())
    need(search['ordinary_parents']==6880 and search['extensions']==1736448,'declared search universe')
    audited=0
    for T in search['audit_targets']:
        h=hp(T);need(h in (95,113,115) and strong(T),'independent low target recount')
        for r in range(9):
            if two_run(T,r):
                a=(hp(inflate(T,r,tt(2)))-h)//2
                need(a!=({95:11,113:81,115:81}[h]),'negative target root audit')
                audited+=1
    for hit in search['hits']:
        T,r=hit['rows'],hit['root'];h=hp(T)
        need(h==599 and strong(T) and two_run(T,r),'hit identity and structure')
        need(hp(inflate(T,r,tt(2)))==1247,'literal shifted positive hit')
    prior=RESULTS/'prime_ternary_lift_20261113/frames.json.gz'
    bank=json.loads(gzip.decompress(prior.read_bytes()))
    blocked={(11,73),(81,-49),(81,-47),(162,-49),(162,-47),(324,-49),(324,-47)}
    need(not any((r['alpha'],r['beta']) in blocked for r in bank['profiles']),'inherited all-order<=8 profiles')
    out['order9_audit']={'target_presentations':len(search['audit_targets']),
        'two_run_root_recounts':audited,'positive_presentations':len(search['hits']),
        'negative_unshifted_profiles':[[81,-49],[81,-47]],
        'larger_shift_absence_scope':'order<=9 only; no all-order claim for shifted negative charts'}
    pins=[prior,RESULTS/'multicore_metagraph_20261113/rooted_growth.json.gz',HERE/'rooted9.cpp',HERE/'search.py']
    out['source_sha256']={str(p.relative_to(RESULTS)).replace('\\','/'):
        hashlib.sha256(p.read_bytes() if p.suffix=='.gz' else p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in pins}
    suffix='' if __debug__ else '.optimized'
    (HERE/('receipt'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'bezout':vals,'shifted':controls,'audit':out['order9_audit']}))

if __name__=='__main__':main()
