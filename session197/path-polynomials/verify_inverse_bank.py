#!/usr/bin/env python3
from verify_extension import *
from math import isqrt
from functools import lru_cache

OPS=[('D',2,D),('E',3,E),('F',4,F),('G',5,G)]
def primes(n):
    sieve=bytearray(b'\1')*(n+1);sieve[:2]=b'\0\0'
    for i in range(2,isqrt(n)+1):
        if sieve[i]:sieve[i*i:n+1:i]=b'\0'*(((n-i*i)//i)+1)
    return [i for i in range(3,n+1,2) if sieve[i]]

def make_profiles(N,bound):
    levels=[{} for _ in range(N+1)];levels[1][(0,1)]=('v',)
    for n in range(2,N+1):
        target=levels[n]
        for name,size,operator in OPS:
            if n>size:
                for p in levels[n-size]:
                    q=tuple(operator(p))
                    if value(q,1)<=bound:target.setdefault(q,(name,p))
        for i in range(1,n//2+1):
            j=n-i
            for a in levels[i]:
                ha=value(a,1)
                for b in levels[j]:
                    if i==j and a>b:continue
                    if ha*value(b,1)>bound:continue
                    q=tuple(mul(a,b));target.setdefault(q,('*',a,b))
    return levels

def run():
    bound=8191;N=bound.bit_length()-1+5;require(N==17,'prime order bound')
    levels=make_profiles(N,bound)
    @lru_cache(None)
    def graph(p):
        rec=levels[len(p)-1][p];op=rec[0]
        if op=='v':return (0,)
        a=graph(rec[1])
        if op=='*':return join(a,graph(rec[2]))
        if op=='D':return wrapper(a)
        if op=='E':return edge_wrapper(a)
        if op=='F':return port_wrapper((4,1,2),(0,1),a)
        if op=='G':return port_wrapper((0,1,3,7),(0,2),a)
        raise RuntimeError(op)
    reps={}
    for level in levels:
        for p in level:reps.setdefault(value(p,1),p)
    ps=primes(bound);hits=[p for p in ps if p in reps];missing=[p for p in ps if p not in reps]
    require(len(ps)==1027 and len(hits)==77,'prime census')
    require(sum(map(len,levels))==12553,'profile census')
    witnesses=[]
    for p in hits:
        a=graph(reps[p]);validate(a);require(len(a)==len(reps[p])-1,'replay degree')
        require(h_subsets(a)[-1]==p,('physical prime replay',p,a))
        witnesses.append({'H':p,'rows':a})
    require(all(p in missing for p in (29,31,41,43)),'missing targets')
    forms=[('Dodd',3,6),('Deven',5,12),('Eodd',5,18),('Eeven',11,36),('Fodd',13,36),('Feven',25,72),('Godd',19,60),('Geven',39,120)]
    sols={}
    for p in (41,43):
        sols[p]=[(name,h,(p-a*h)//b) for name,a,b in forms for h in range(1,p+1,2) if p>=a*h and (p-a*h)%b==0]
    require(sols[41]==[('Deven',1,3),('Eodd',1,2)],('41 reduction',sols[41]))
    require(sols[43]==[('Eodd',5,1)],('43 reduction',sols[43]))
    require(40&(40-1)!=0 and 14&(14-1)!=0,'H1 power-of-two hostile')
    require({h_subsets(tournament(n,m))[-1] for n in (1,3) for m in range(1<<(n*(n-1)//2))}=={1,3},'small odd count hostile')
    a197=port_wrapper((4,1,2),(0,1),wrapper(transitive(2)))
    a211=wrapper(join((0,),edge_wrapper(transitive(2))))
    require(h_subsets(a197)[-1]==197 and h_subsets(a211)[-1]==211,'explicit mixed-generator outputs')
    print(json.dumps({'prime_bound':bound,'complete_order_bound':N,'profiles_per_order':[len(x) for x in levels[1:]],'total_profiles':sum(map(len,levels)),'odd_primes_checked':len(ps),'attainable_primes':hits,'omitted_prime_count':len(missing),'first_omitted_primes':missing[:40],'ternary_child_integer_candidates':sols,'explicit_mixed_outputs':[{'H':197,'rows':a197},{'H':211,'rows':a211}],'all_prime_witnesses':witnesses},sort_keys=True,indent=2))
if __name__=='__main__':run()
