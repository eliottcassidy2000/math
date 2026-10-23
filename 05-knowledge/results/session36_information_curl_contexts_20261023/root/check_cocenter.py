#!/usr/bin/env python3
import itertools,json
from collections import Counter,defaultdict

def need(x,msg):
    if not x: raise RuntimeError(msg)
def cyc(w):return min(w[i:]+w[:i] for i in range(len(w))) if w else ()
def product(f,g):
    h=Counter()
    for u,a in f.items():
        for v,b in g.items():h[u+v]+=a*b
    return {w:c for w,c in h.items() if c}
def bracket(f,g):
    h=Counter(product(f,g))
    for w,c in product(g,f).items():h[w]-=c
    return {w:c for w,c in h.items() if c}
def trace(f):
    h=Counter()
    for w,c in f.items():h[cyc(w)]+=c
    return {w:c for w,c in h.items() if c}

def main():
    A=[1,1,1,2,4,12,56,456,6880];S=[0]*len(A)
    for n in range(1,len(A)):S[n]=A[n]-sum(S[j]*A[n-j] for j in range(1,n))
    words=[[()]]+[[] for _ in range(8)]
    for n in range(1,9):
        for k in range(1,n+1):
            for idx in range(S[k]):
                for tail in words[n-k]:words[n].append(((k,idx),)+tail)
        need(len(words[n])==A[n],'word census')
    first=None;tested=0
    for total in range(1,9):
        failures=0
        for n in range(1,total):
            groups=defaultdict(list)
            for w in words[n]:groups[cyc(w)].append(w)
            for group in groups.values():
                if len(group)<2:continue
                for tail in words[total-n]:
                    images={cyc(w+tail) for w in group};tested+=1
                    if len(images)>1:
                        failures+=1
                        if first is None:first=dict(order=total,left=group[0],right=next(w for w in group if cyc(w+tail)!=cyc(group[0]+tail)),tail=tail)
        if failures:break
    need(first is not None and first['order']==8,'first join failure')
    a=((1,0),);b=((3,0),)
    need(cyc(a+b)==cyc(b+a),'initial equality')
    need(cyc(a+b+a+b)!=cyc(b+a+a+b),'specified hostile')
    checks=0
    for length in range(2,7):
        for word in itertools.product(range(3),repeat=length):
            f={(word[0],):1}
            for x in word[1:]:f=bracket(f,{(x,):1})
            need(not trace(f),'bracket survives cocenter');checks+=1
    report=dict(status='PASS',strong_alphabet_counts=S[1:],word_counts=A[1:],first_join_failure=first,context_tests=tested,nested_bracket_tests=checks,scope='formal ordered-SCC type algebra; exact type counts through order8')
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
