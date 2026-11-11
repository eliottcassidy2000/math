"""Independent permutation/split controls for all selected-set seeds in the SCC family."""
from itertools import permutations
from math import factorial
from pathlib import Path
from hashlib import sha256
import json

HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def words(n):
    if n==0:
        yield ()
    if n>=1:
        for rest in words(n-1):yield (1,)+rest
    if n>=3:
        for rest in words(n-3):yield (3,)+rest

def graph(word):
    n=sum(word);rows=[];at=0
    for size in word:
        later=((1<<n)-1)^((1<<(at+size))-1)
        block=(0,) if size==1 else (2,4,1)
        rows.extend((r<<at)|later for r in block);at+=size
    return rows

def literal(rows):
    """Every ordered two-cover is uniquely a permutation with one split."""
    n=len(rows);h=0;starts=[0]*n;K=[[0]*n for _ in range(n)]
    for order in permutations(range(n)):
        bad=[]
        for j in range(n-1):
            if not(rows[order[j]]>>order[j+1]&1):
                bad.append(j)
                if len(bad)==2:break
        if len(bad)>1:continue
        if not bad:
            h+=1
            if n:starts[order[0]]+=1
            cuts=range(n-1)
        else:cuts=bad
        for cut in cuts:K[order[cut]][order[cut+1]]+=1
    U=[starts[j]+sum(K[i][j] for i in range(n)) for j in range(n)]
    return h,K,U

def tail(rows,selected,t):
    n=len(rows)
    enlarged=[row|sum(1<<j for j in range(n,n+t)) for row in rows]
    enlarged.extend(sum(1<<j for j in range(i+1,n+t)) for i in range(n,n+t))
    return [row|((1<<(n+t)) if not(selected>>i&1) else 0) for i,row in enumerate(enlarged)]+[selected]

def main():
    sources=subsets=positives=0;census=[];small_profiles=[]
    for n in range(1,9):
        classes=set();labelled=0
        for word in words(n):
            rows=graph(word);h,K,U=literal(rows);sources+=1;r=word.count(3)
            need(h==3**r,('source count',word,h))
            for mask in range(1<<n):
                subsets+=1;selected=[v for v in range(n) if mask>>v&1]
                w=sum(K[i][j] for i in selected for j in selected);u=sum(U[j] for j in selected)
                counts=[];at=0
                for size in word:
                    counts.append(sum((mask>>(at+j))&1 for j in range(size)));at+=size
                predicted=(r==0 and len(selected)==1) or (
                    r==1 and sum(counts)==2 and counts[word.index(3)]==2)
                actual=h-w==1 and u>0
                need(actual==predicted,('seed classification',word,mask,h,w,u))
                if not actual:continue
                positives+=1
                selected_core=next(i for i,value in enumerate(counts) if value)
                b=len(word)-selected_core-1
                expected=(1 if r==0 else 6)*(1<<b)
                need(u==expected,('sink multiplier',word,mask,u,expected))
                classes.add(tuple(zip(word,counts)))
                labelled+=factorial(n)//(3**r)
                if n<=3:small_profiles.append(dict(word=word,selected_mask=mask,H=h,W=w,U=u))
        expected_classes=n+max(n-2,0)
        need(len(classes)==expected_classes,('ordinary coefficient',n,len(classes)))
        need(labelled==factorial(n)*expected_classes,('labelled mass',n,labelled))
        census.append(dict(n=n,ordinary_marked=expected_classes,labelled_marked=labelled,
                           EGF_coefficient=labelled//factorial(n)))
    tails=[]
    for word,mask in [((1,),1),((3,),3),((1,3,1),6)]:
        base=graph(word);h,K,U=literal(base);u=sum(U[v] for v in range(len(base)) if mask>>v&1)
        for t in (1,2):
            actual=literal(tail(base,mask,t))[0]
            need(actual==1+(1<<t)*u,('literal tail',word,mask,t,actual))
            tails.append(dict(word=word,selected_mask=mask,t=t,U=u,H=actual))
    empty=literal([])
    need(empty[0]==1 and empty[2]==[],'empty source convention')
    inherited=HERE.parents[1]/'session258_endpoint_seed_20261111/seed_decision/CORE_CLASSIFICATION.md'
    report=dict(status='PROVED family specialization; FINITE-EXACT independent permutation/split controls',
                universe='Every ordinary K1/C3 SCC word of orders1..8, every literal selected subset, including empty subsets',
                sources=sources,selected_subsets=subsets,positive_literal_selected_presentations=positives,
                ordinary_and_labelled=census,small_profiles=small_profiles,tail_controls=tails,
                positive_exact_multipliers='2^b or6*2^b, b>=0',odd_slopes=[1,3],
                ordinary_GF='(z+z^3)/(1-z)^2',labelled_EGF='(z+z^3)/(1-z)^2',
                inherited_core_classification_sha256=sha256(inherited.read_bytes()).hexdigest())
    HERE.joinpath('TRIANGLE_seeds.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(result='PASS',sources=sources,subsets=subsets,positive_presentations=positives,census=census)))

if __name__=='__main__':main()
