from itertools import product
from functools import reduce
from math import gcd
from replay_ordinal_digits import BANK, profile, state, append, strong, external, ordinal

ORDER=(0,1,3,2,5,6,4)
A=[(profile(a,q),profile(a,r)) for a,q,r,*_ in BANK]
J,W,K0,E0=state(A,ORDER,(0,)*7)
STEPS=tuple((A[i][1][-1]-A[i][0][-1])*(J//A[i][0][0]) for i in ORDER)
PALETTE=[(sum(x*y for x,y in zip(word,STEPS)),word)
         for word in product(range(2),repeat=7)]
def need(x):
    if not x: raise AssertionError("joint residue verification failed")

def compile_data(m,reps):
    mod=14*m
    if gcd(J,mod)!=1: raise ValueError("modulus meets the deletion counts")
    layer={0:()}
    local={d%mod:(d,word) for d,word in PALETTE}
    baseline=(1,0,1,0)
    for _ in range(reps):
        target={}
        for a,words in layer.items():
            for b,(_,v) in local.items():
                u=(a+b)%mod
                if u not in target:target[u]=words+(v,)
        layer=target
        for i in ORDER:baseline=append(baseline,A[i][0])
    need(len(layer)==mod//2)
    return layer,baseline

def compile_residue(m,reps,c,k,data):
    layer,(Jr,Wr,Kr,Er)=data
    mod=14*m
    need(c%2 and k>=1 and Wr%7==0)
    power=J**(reps-1)
    residue=(7*c-Wr*2**k-Kr)*pow(power,-1,mod)%mod
    words=layer[residue]
    drift=sum(sum(a*b for a,b in zip(STEPS,v)) for v in words)
    total=Wr*2**k+Kr+power*drift
    need(total%7==0 and (total//7-c)%(2*m)==0 and total>0)
    return words,total

def main():
    need(W%7==0 and reduce(gcd,STEPS)==2)
    # All odd classes modulo 144: 3 copies of the common-boundary stack.
    m=72; reps=3; data=compile_data(m,reps)
    for k in (1,4,7):
        for c in range(1,2*m,2):compile_residue(m,reps,c,k,data)
    words,total=compile_residue(m,reps,97,1,data)
    rows=();q=0
    for v in words:
        for i,b in zip(ORDER,v):
            q|=BANK[i][1+b]<<len(rows)
            rows=ordinal(rows,BANK[i][0])
    need(strong(external(rows,q,1)))
    # A small generic modulus, using the proved upper bound on copies.
    m=5; reps=7*m-1; data=compile_data(m,reps)
    for c in range(1,2*m,2):compile_residue(m,reps,c,1,data)
    print("Three-copy joint residue classes",72,"example quotient",total//7)
    print("Coprime modulus construction and physical strongness: PASS")
if __name__=="__main__":
    main()
