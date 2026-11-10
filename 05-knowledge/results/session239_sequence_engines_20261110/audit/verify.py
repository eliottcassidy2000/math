#!/usr/bin/env python3
"""Independent audit; every check remains active under python -O."""
from pathlib import Path
import importlib.util, itertools, math, json, time
from collections import Counter
BASE=Path(__file__).resolve().parents[1]
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
E=load('engines',BASE/'burnside/engines.py')
B=load('bridge',BASE/'proof_bridges/verify.py')
checks=Counter()
def check(ok,label):
    checks[label]+=1
    if not ok: raise RuntimeError(label)
def vp(x): return (x & -x).bit_length()-1
def profiles(n,d=1):
    if not n: yield {}; return
    if d>n:return
    for k in range(n//d+1):
        for tail in profiles(n-k*d,d+2):
            yield ({d:k}|tail) if k else tail

def phi(d): return sum(math.gcd(i,d)==1 for i in range(1,d+1))
def independent(n):
    """Davis via divisor-square identity, grouped multiplicities and integer class sizes."""
    total=0; fact=math.factorial(n)
    for m in profiles(n):
        r=sum(m.values())
        energy=sum(phi(d)*sum(k for j,k in m.items() if j%d==0)**2 for d in range(1,n+1))
        check((energy-r)%2==0,'divisor exponent integral')
        z=math.prod(j**k*math.factorial(k) for j,k in m.items())
        check(fact%z==0,'class size integral')
        total+=(fact//z) << ((energy-r)//2)
        if n%2:
            delta=(energy-r)//2-(n-1)//2-vp(z)
            lower=(r-1)**2//2-vp(math.factorial(r))
            check(delta>=lower,'precision bound')
    check(total%fact==0,'Davis terminal integral')
    return total//fact

def main():
    start=time.perf_counter(); exact,_=E.integer_response(120)
    for n in range(46): check(independent(n)==exact[n],'independent Davis 0..45')
    ref=BASE.parent/'tournament_future_boundary_20261020/a000568-reference.txt'
    rows=[tuple(map(int,line.split())) for line in ref.read_text().splitlines() if line.strip() and not line.startswith('#')]
    for n,a in rows: check(exact[n]==a,'external reference')
    for n in [0,1,2,3,4,7,12,31,64,120]:
        check(E.integer_response(n)[0]==exact[:n+1],'horizon consistency')
        for mod in [1,2,4,6,9,72,1001]:
            check(E.integer_response(n,mod)[0]==[a%mod for a in exact[:n+1]],'composite residue')
        for bits in [1,2,3,7,32]:
            check(E.integer_response(n,binary_precision=bits)[0]==[a%(1<<bits) for a in exact[:n+1]],'binary residue')
    for args,kwargs in [((-1,),{}),((5,0),{}),((5,-1),{}),((5,),{'binary_precision':0}),((5,2),{'binary_precision':2})]:
        try:E.integer_response(*args,**kwargs)
        except AssertionError:check(True,'invalid input rejected')
        else:check(False,'invalid input rejected')
    for n in range(3,120,2):
        u=exact[n]>>((n-1)//2)
        check(u%4==(-n)%4,'mod4')
        if n>=5:
            c=sum(math.gcd(a,b)==math.gcd(a,c)==math.gcd(b,c)==1 for a in range(1,n,2) for b in range(a+2,n-a,2) for c in [n-a-b] if c>b)
            check(u%8==(pow(n,-1,8)+2*pow(n-2,-1,8)+4*c)%8,'mod8')
    check((exact[3]//2)%8!= (pow(3,-1,8)+2*pow(1,-1,8))%8,'n3 exception')
    # Independently encode every labelled tournament <=5, count paths by permutations,
    # and recognize strong connectivity by repeated set closure, not worker bit DP.
    fixed={():Counter({1:1})}; strong={():Counter()}
    for n in range(1,6):
        edges=list(itertools.combinations(range(n),2)); facts=[]
        for signs in itertools.product([0,1],repeat=len(edges)):
            arcs={(a,b) if s else (b,a) for (a,b),s in zip(edges,signs)}
            h=sum(all((p[i],p[i+1]) in arcs for i in range(n-1)) for p in itertools.permutations(range(n)))
            reach={0}
            for _ in range(n): reach|={v for u,v in arcs if u in reach}
            back={0}
            for _ in range(n): back|={u for u,v in arcs if v in back}
            facts.append((arcs,h,len(reach)==len(back)==n))
        for part in B.partitions(n):
            p=B.permutation(part)
            accepted=[(h,s) for arcs,h,s in facts if {(p[a],p[b]) for a,b in arcs}==arcs]
            fixed[part]=Counter(h for h,s in accepted)
            strong[part]=Counter(h for h,s in accepted if s)
            got=fixed[part].copy()
            for sub,rest,mult in B.subprofiles(part):
                if sub and rest:
                    for x,c in strong[sub].items():
                        for y,d in fixed[rest].items(): got[x*y]-=mult*c*d
            check(Counter({h:c for h,c in got.items() if c})==strong[part],'independent SCC recurrence')
            check(all(h%math.lcm(*part)==0 for h,s in accepted),'independent free action')
    print(json.dumps({'status':'PASS','checks':dict(checks),'reference_max':max(n for n,a in rows),'seconds':time.perf_counter()-start},indent=2))
if __name__=='__main__':main()
