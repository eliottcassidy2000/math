"""Exact finite controls; all acceptance gates survive python -O."""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from functools import lru_cache
from math import comb, gcd
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
def need(b,msg):
    if not b: raise RuntimeError(msg)
def pairs(n): return list(combinations(range(n),2))
def strong(E,V):
    V=set(V)
    for s in V:
        seen={s}
        for _ in V: seen|={b for a,b in E if a in seen and b in V}
        if seen!=V:return False
    return bool(V)
def completions(E,n):
    free=[e for e in pairs(n) if e not in E and e[::-1] not in E]
    for bits in product((0,1),repeat=len(free)):
        yield E|{e if bit else e[::-1] for e,bit in zip(free,bits)}
def direct(E,n):return sum(strong(F,range(n)) for F in completions(E,n))
def recurrence(E,n):
    @lru_cache(None)
    def C(mask):
        V={v for v in range(n) if mask>>v&1}
        return 1<<(comb(len(V),2)-sum(a in V and b in V for a,b in E))
    @lru_cache(None)
    def S(mask):
        if mask==0:return 0
        out=C(mask);A=(mask-1)&mask
        while A:
            B=mask^A
            if not any((B>>a&1) and (A>>b&1) for a,b in E):out-=S(A)*C(B)
            A=(A-1)&mask
        return out
    return S((1<<n)-1)
def canon(E,n):
    return min(tuple(int((p[i],p[j]) in E)-int((p[j],p[i]) in E)
                     for i,j in pairs(n)) for p in permutations(range(n)))
FREE=[(0,2),(0,3),(1,3)]
def tiling(b):
    return {(0,1),(1,2),(2,3)}|{e if b>>k&1 else e[::-1] for k,e in enumerate(FREE)}
def isoword(a,b):
    A,B=tiling(a),tiling(b)
    for p in permutations(range(4)):
        if all(((p[i],p[j]) in A)==((i,j) in B) for i,j in pairs(4)):return list(p)
    raise RuntimeError('illegal re-mark')
def quadratic(a):return (a*a-29)//4
def primes(n):
    return [x for x in range(2,n+1) if all(x%d for d in range(2,int(x**.5)+1))]
def prime(n):return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
def firsthit(m):
    p=m-1
    while not (prime(p) and p%m in (1,m-1)):p+=1
    return p
def full_pair_counts(P,cap):
    ans=Counter()
    for i,p in enumerate(P):
        for q in P[i:]:
            if p+q<=cap:ans[p+q]+=1
    return ans
def streaming(P,cap):
    seen=[];out=Counter()
    for p in P:
        for q in seen:
            if p+q<=cap:out[p+q]+=1
        if 2*p<=cap:out[2*p]+=1
        seen.append(p)
    return out
def formal_series(cap):
    U=[1<<(k*(k-1)//2) for k in range(cap+1)]
    T=[1]+[1<<comb(k-1,2) for k in range(1,cap+1)]
    s=[0]*(cap+1);J=[0]*(cap+1)
    for k in range(1,cap+1):
        s[k]=U[k]-sum(comb(k,j)*s[j]*U[k-j] for j in range(1,k))
        J[k]=T[k]-sum(J[j]*T[k-j] for j in range(1,k))
    C=[1<<(m*(m+3)//2) for m in range(cap+1)]
    star=[C[m]-sum(comb(m,k)*s[k]*C[m-k] for k in range(1,m+1)) for m in range(cap+1)]
    chain=[C[m]-sum(J[k]*C[m-k] for k in range(1,m+1)) for m in range(cap+1)]
    return U,T,s,J,C,star,chain
def main():
    B=[-7,-5,-3,-1,1,3,5,7]
    need(set(map(quadratic,B))=={-7,-5,-1,5},'quadratic image')
    need(set(Counter(map(quadratic,B)).values())=={2},'even fibres')
    magnitude={a:abs(quadratic(a)) for a in (1,3,5,7)}
    E={(0,3),(3,2),(2,0),(1,2)} # vertex labels [1,3,5,7]
    need(direct(E,4)==3,'three strong completions')
    D=defaultdict(list);A=tiling(0)
    for removed in combinations(sorted(A),2):D[canon(A-set(removed),4)].append(removed)
    need(len(D)==9 and len(D[canon(E,4)])==3,'deletion catalogue')
    catalogue=[]
    for key,rs in sorted(D.items()):
        F=A-set(rs[0]);catalogue.append({'code':key,'multiplicity':len(rs),'removed':rs,
            'outdegrees':sorted(sum(a==v for a,b in F) for v in range(4)),
            'strong_completions':direct(F,4)})
    # Whole universe of partial orientations, not just selected trees.
    partial_cases=0
    for n in range(1,5):
        for word in product((-1,0,1),repeat=comb(n,2)):
            F={e if x==1 else e[::-1] for e,x in zip(pairs(n),word) if x}
            need(direct(F,n)==recurrence(F,n),'source-SCC recurrence')
            partial_cases+=1
    U,T,s,J,C,star,chain=formal_series(12)
    for m in range(5):
        Fs={(0,1),(1,2),(2,0)}|{(i,0) for i in range(3,m+3)}
        Fc={(0,1),(1,2),(2,0)}|{(i,0 if i==3 else i-1) for i in range(3,m+3)}
        need(recurrence(Fs,m+3)==star[m],'star formula')
        need(recurrence(Fc,m+3)==chain[m],'chain formula')
        if m<=3:
            need(direct(Fs,m+3)==star[m] and direct(Fc,m+3)==chain[m],'literal family')
    forest_cases=0
    for n in range(3,7):
        for parents in product(*(range(j) for j in range(3,n))):
            F={(0,1),(1,2),(2,0)}|{(j,p) for j,p in zip(range(3,n),parents)}
            need(direct(F,n)==recurrence(F,n),'all increasing forests')
            forest_cases+=1
    fibres=defaultdict(list)
    for b in range(8):fibres[canon(tiling(b),4)].append(b)
    need(sorted(map(len,fibres.values()))==[1,1,1,5],'odd tiling fibres')
    need(sum((b^7) in (0,1,2,4,5) for b in (0,1,2,4,5))==2,'cycle placement guard')
    code={-7:3,-5:0,-3:6,-1:1,1:2,3:7,5:4,7:5}
    before=[2,4,4,3,2,5,6,7];after=[0,0,4,3,4,1,6,7]
    witnesses=[]
    for a in B:
        b=code[a];r=before[b];c=r^7;t=after[c]
        need(t==code[quadratic(a)],'native R C R')
        witnesses.append({'a':a,'b':b,'before':r,'complement':c,'after':t,
                          'first_path':isoword(b,r),'second_path':isoword(c,t)})
    # Source compiler: one control prime is legal, not arbitrary products.
    state_primes=[2,3,5,7,11,13,17,19]
    codes=dict(zip(B,state_primes));rules=[(codes[quadratic(a)],codes[a]) for a in B]
    def fracstep(x):
        for u,v in rules:
            if x%v==0:return x//v*u
        return None
    for a in B:
        for scale in (1,23,23*29):
            need(fracstep(scale*codes[a])==scale*codes[quadratic(a)],'FRACTRAN')
    # A strict radius-zero CA; shift check on a periodic test ring.
    for w in product(B,repeat=3):
        shift=lambda z:z[1:]+z[:1]
        evolve=lambda z:tuple(map(quadratic,z))
        need(evolve(shift(w))==shift(evolve(w)),'CA shift')
    clocks=[]
    for m in range(4,402,2):
        hit=firsthit(m);P=[p for p in primes(hit) if p>2]
        need(hit>=m-1,'first-hit bound')
        target=full_pair_counts([p for p in primes(m-3) if p>2],m)
        need(streaming(P,m)==target,'first-hit triangle closure')
        if m in (10,20,40,78,290):clocks.append({'modulus':m,'first_hit':hit,'pairs_at_modulus':target[m]})
    need(firsthit(26)==53 and prime(47) and not prime(21) and 47%26==21,'residue alias hostile')
    # Every formal Boolean prefix, independent of primality realization.
    for n in range(9):
        for bits in product((0,1),repeat=n):
            P=[2*i+3 for i,b in enumerate(bits) if b]
            need(streaming(P,4*n+6)==full_pair_counts(P,4*n+6),'formal pair source')
    # All CRT masks and their exact cyclic periods at these squarefree odd moduli.
    wheel=[]
    for R in (3,5,7,15,21,35,105):
        pp=[p for p in primes(R) if R%p==0]
        for h in range(R):
            word=[int(gcd(c-h,R)==gcd(c+h,R)==1) for c in range(R)]
            size=1
            for p in pp:size*=p-(1 if h%p==0 else 2)
            need(sum(word)==size,'CRT row density')
            period=next(d for d in range(1,R+1) if all(word[i]==word[(i+d)%R] for i in range(R)))
            need(period==R,'minimal row period')
        wheel.append(R)
    # The inherited clique dimension has a finite boundary impulse.
    def dim(n,s):return 1+2*comb(n,2)+sum(comb(n,k) for k in range(3,min(n,s-1)+1))+2*(n>=3)+(n>=s)
    for ss in (5,6,7,11):
        for n in range(3,60):
            delta=2*n+sum(comb(n,k) for k in range(2,min(n,ss-2)+1))+(n+1==ss)
            need(dim(n+1,ss)-dim(n,ss)==delta,'dimension birth law')
    out={'status':'FINITE-EXACT controls for scoped proofs; same-agent audit',
         'bounded_quadratic':{str(a):quadratic(a) for a in B},'magnitude_map':magnitude,
         'deletion_catalogue':catalogue,'partial_orientation_universe':partial_cases,
         'increasing_forest_universe':forest_cases,
         'star_counts':star,'chain_counts':chain,'strong_labelled_counts':s,
         'strong_fixed_path_counts':J,'tiling_fibres':list(fibres.values()),
         'native_code':code,'native_before':before,'native_after':after,'native_witnesses':witnesses,
         'fractran_rules':rules,'prime_clocks':clocks,'wheel_moduli':wheel,
         'clique_dimensions_335':[dim(n,5) for n in range(3,11)]}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'partial_cases':partial_cases,'forest_cases':forest_cases,'deletion_types':len(D),
                      'star_first':star[:5],'chain_first':chain[:5],'clocks':clocks},sort_keys=True))
if __name__=='__main__':main()
