"""Complete ordinary triangle-chain census and internal metagraph controls."""
from collections import Counter
from math import comb,factorial
from pathlib import Path
import json
from extension_content import C3,join,extend,need
from gap_transport import flip

def words(n):
    if n==0:yield ()
    if n>=1:
        for w in words(n-1):yield (1,)+w
    if n>=3:
        for w in words(n-3):yield (3,)+w

def matching(w):
    answer=run=0
    for x in w+(1,):
        if x==3:run+=1
        else:answer+=run//2;run=0
    return answer

def classify(adj):
    n=len(adj);reach=[row|(1<<i) for i,row in enumerate(adj)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    left=set(range(n));components=[]
    while left:
        i=min(left);part={j for j in left if reach[i]>>j&1 and reach[j]>>i&1}
        if len(part) not in (1,3):return None
        components.append(part);left-=part
    components=sorted(components,key=lambda c:-sum(1 for d in components if c is not d and adj[min(c)]>>min(d)&1))
    return tuple(len(c) for c in components)

def gf_lifetime(n,l):
    if l==0:return 1
    result=0
    for a in range(l):
        for b in (0,3):
            rest=n-(4*l-1+2*a+b)
            if rest>=0:result+=comb(l-1,a)*comb(rest+l,l)
    return result

def main():
    max_n=24; recurrence=[];census=[]
    for n in range(max_n+1):
        row=Counter()
        if n==0:row[(0,0)]+=1
        if n==3:row[(1,0)]+=1
        if n>=1:row.update(recurrence[n-1])
        if n>=4:row.update({(r+1,v):count for (r,v),count in recurrence[n-4].items()})
        if n>=6:row.update({(r+2,v+1):count for (r,v),count in recurrence[n-6].items()})
        recurrence.append(row)
        exact=Counter((w.count(3),matching(w)) for w in words(n))
        need(row==exact,'independent joint rational-GF coefficient census')
        life=Counter();labelled=0
        for (r,v),count in row.items():
            life[r-v]+=count
            need(factorial(n)%3**r==0,'ordinary-to-labelled integer mass')
            labelled+=count*(factorial(n)//3**r)
        need(all(life[l]==gf_lifetime(n,l) for l in range(n//3+1)),
             'exact lifetime coefficient formula')
        census.append({'n':n,'ordinary':sum(row.values()),'labelled':labelled,
                       'lifetimes':dict(sorted(life.items()))})
    birth_cases=arc_cases=0
    for n in range(10):
        for word in words(n):
            graph=join([C3 if x==3 else (0,) for x in word])
            need(classify(graph)==word,'canonical SCC word')
            predicted=Counter()
            for p in range(len(word)+1):predicted[word[:p]+(1,)+word[p:]]+=1
            for p in range(len(word)-1):
                if word[p:p+2]==(1,1):predicted[word[:p]+(3,)+word[p+2:]]+=1
            actual=Counter()
            for mask in range(1<<n):
                result=classify(extend(graph,mask))
                if result is not None:actual[result]+=1
                birth_cases+=1
            need(actual==predicted,'complete within-family birth multiplicities')
            predicted_flips=Counter()
            for p,x in enumerate(word):
                if x==3:predicted_flips[word[:p]+(1,1,1)+word[p+1:]]+=3
            for p in range(len(word)-1):
                if word[p:p+2]==(1,1):predicted_flips[word]+=1
            for p in range(len(word)-2):
                if word[p:p+3]==(1,1,1):predicted_flips[word[:p]+(3,)+word[p+3:]]+=1
            actual_flips=Counter()
            for i in range(n):
                for j in range(i+1,n):
                    result=classify(flip(graph,i,j))
                    if result is not None:actual_flips[result]+=1
                    arc_cases+=1
            need(actual_flips==predicted_flips,'complete within-family arc-flip multiplicities')
    out={'status':'FINITE-EXACT independent controls for proved joint enumeration and internal transitions',
         'joint_GF':'(1+u*z^3)/(1-z-u*z^4-v*u^2*z^6)',
         'lifetime_GF':'(1+w*z^3)/(1-z-w*z^4-w*z^6)',
         'census':census,'complete_birth_masks':birth_cases,'complete_arc_flips':arc_cases}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'result':'PASS','orders':max_n,'birth_masks':birth_cases,'arc_flips':arc_cases}))

if __name__=='__main__':main()
