"""Independent permutation path-bin reader, odd-cycle reader, signed inverse controls."""
from itertools import permutations
from math import comb
import json,gzip
from repair23 import hp2,near,transitive
from fibre_probe import HERE,save,need

def independent(T):
    n=len(T);covers=[0]*5;cycle=[0]*(n+1)
    for p in permutations(range(n)):
        bad=sum(not(T[a]>>b&1) for a,b in zip(p,p[1:]))
        for z in range(1,5):
            if z>=bad+1:covers[z]+=comb(n+z-1-bad,n)
        seen=set();cycles=nontrivial=0;ok=True
        for i in range(n):
            if i in seen:continue
            v=i;length=0
            while v not in seen:
                seen.add(v);length+=1
                if p[v]!=v and not(T[v]>>p[v]&1):ok=False
                v=p[v]
            cycles+=1;nontrivial+=length>1
            if length%2==0:ok=False
        if ok:cycle[cycles]+=1<<nontrivial
    need(all(covers[z]==sum(c*z**k for k,c in enumerate(cycle)) for z in range(1,5)),'independent cycle/path-bin identity')
    need(covers[2]%3==((-1)**n*covers[1])%3,'mod3 phase')
    need((covers[1],covers[2])==hp2(T),'independent reference for numeric engine')
    return covers[1:],cycle

def run():
    rows=json.loads(gzip.decompress((HERE/'fibre_probe.json.gz').read_bytes()))['records'];ordinary={}
    for row in rows:
        if row['n']<=6:ordinary.setdefault((row['n'],row['code']),tuple(row['rows']))
    inputs=list(ordinary.values())+[transitive(n) for n in range(1,7)]
    records=[]
    for T in inputs:
        covers,cycles=independent(T);records.append(dict(rows=T,covers=covers,cycle_coefficients=cycles))
    for n in range(3,11):
        m=n-2;need(hp2(near(n))==(2**m+1,4*2**m+2*3**m-2),'near-transitive formula')
    sequence=[3*n//2+1 for n in range(1,101)]
    need(sequence==[v for v in range(2,sequence[-1]+1) if v%3],'unit sequence')
    inverse=0
    for y in range(1,200,2):
        for epsilon in (-1,1):
            for k in range(1,21):
                numerator=(1<<k)*y-epsilon
                integral=numerator%3==0
                need(integral==(y%3!=0 and ((k%2==0)==(y%3==epsilon%3))),'inverse guard')
                if integral:
                    x=numerator//3
                    need(x>0 and x%2==1 and 3*x+epsilon==(1<<k)*y,'exact inverse depth')
                    need((((1<<(k+2))*y-epsilon)//3)==4*x+epsilon,'same-sheet radix4')
                    need((((1<<(k+3))*y+epsilon)//3)==8*x+3*epsilon,'cross-sheet radix8')
                    inverse+=1
    v2=lambda x:(x&-x).bit_length()-1
    phases={1:(2,1),3:(1,3),5:(3,1),7:(1,2)}
    for x in range(1,1000,2):
        for observed,expected in zip((v2(3*x+1),v2(3*x-1)),phases[x%8]):
            need((observed>=3) if expected==3 else (observed==expected),'mod8 valuation')
    save('arithmetic_probe',dict(status='PASS',scope='All43 ordinary strong tournaments through6 plus six transitive controls, independent permutation/bin and cycle readers at z1..4; R_n n3..10; all odd y<200, signs+-1, k1..20; all odd x<1000 mod8 phases',graphs=len(inputs),records=records,accepted_inverse_checks=inverse,sequence=sequence))
    print('PASS',len(inputs),'independent cycle/path graphs;',inverse,'signed inverse candidates')

if __name__=='__main__':run()
