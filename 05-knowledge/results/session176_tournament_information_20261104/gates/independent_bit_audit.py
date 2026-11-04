#!/usr/bin/env python3
import itertools,json

def require(p,s):
    if not p:raise RuntimeError(s)
def bits(x,w):return [(x>>i)&1 for i in range(w)]
def value(a):return sum(v<<i for i,v in enumerate(a))
def add(a,b):
    out=[];carry=0;gates=0
    for x,y in zip(a,b):
        p=x^y;g=x&y;t=p&carry;new=g|t;out.append(p^carry);carry=new;gates+=5
    require(carry==0,'overflow')
    return out,gates

def bit_circuit(word):
    n=len(word);width=n+1;P=[0]*width;L=0;h=bits(1,width);gates=0
    for c in word:
        s,g=add(P,bits(L,width));gates+=g
        nc=1-c;gates+=1
        term=[nc&v for v in s];gates+=width
        h,g=add(h,term);gates+=g
        require(P[-1]==0,'shift overflow');P=[c]+P[:-1];L=c
    return value(h),gates

def direct_formula(b):
    n=len(b)
    if not n:return 1
    ans=b[0]+1-b[-1]
    for i in range(n-1):ans+=(1-b[i])*b[i+1]
    for j in range(n):
        for i in range(j+1,n):ans+=(1<<(i-j-1))*b[j]*(1-b[i])
    return ans

def main():
    tested=0;largest=0
    for n in range(11):
        for w in itertools.product((0,1),repeat=n):
            h,g=bit_circuit(w);require(h==direct_formula(w),'bit expansion differs');require(h<=1<<n,'count bit bound');require(g==n*(11*(n+1)+1),'gate inventory');tested+=1;largest=max(largest,g)
    # Algebraic state separation, checked without calling the streaming recurrence.
    probe=0
    for k in range(13):
        pairs=set()
        for w in itertools.product((0,1),repeat=k):
            a=direct_formula(w+(0,0));b=direct_formula(w+(0,1));p=sum(c<<(k-i-1) for i,c in enumerate(w))
            require(a-b==2*p,'probe reconstruction');pairs.add((a,b));probe+=1
        require(len(pairs)==1<<k,'response collision')
    print(json.dumps({'status':'PASS','bit_level_words_through_length_10':tested,'exact_gate_formula':'n*(11*(n+1)+1)','largest_test_gate_count':largest,'direct_formula_two_suffix_probes_through_k12':probe,'semantics':'multioutput binary Hamilton count; fixed n circuit; constants and wires supplied'},indent=2,sort_keys=True))
if __name__=='__main__':main()
