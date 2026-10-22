from math import comb,gcd,prod
from itertools import permutations,product
from collections import Counter
import subprocess,re,argparse
from pathlib import Path

parser=argparse.ArgumentParser(description="Independent Session25 congruence audit; retains checks under python -O.")
parser.add_argument('--repo', type=Path, help='Path to the math repository containing checkpoint f002fc425336')
args=parser.parse_args()
def discover_repo():
    if args.repo is not None:
        return args.repo.resolve()
    for start in (Path(__file__).resolve().parent, Path.cwd().resolve()):
        for parent in (start, *start.parents):
            for candidate in (parent, parent/'math'):
                if (candidate/'01-canon').is_dir() and (candidate/'.git').exists():
                    return candidate
    parser.error('Cannot discover math repository; pass --repo /path/to/math')
repo=discover_repo()
def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
def strong_counts(N,mod):
    s=[0]*(N+1); T=[pow(2,n*(n-1)//2,mod) for n in range(N+1)]
    for n in range(1,N+1):
        b=1; acc=0
        for k in range(1,n):
            b=b*(n-k+1)//k; acc+=b*s[k]*T[n-k]
        s[n]=(T[n]-acc)%mod
    return s
primes=[p for p in range(3,100) if prime(p)]
M=prod(primes)*9
s=strong_counts(1100,M)
c=lambda n,mod: pow(2,(n-1)//2,mod) if n%2 else 0
checks=0
for p in primes:
    eps=pow(2,(p-1)//2,p)
    for m in range(1,1001):
        require((s[m+p-1]-eps*s[m])%p==0, f"quasi-periodicity p={p}, m={m}"); checks+=1
print('quasi-periodicity independent recurrence:',checks,'checks')
expected_residues={9:(7,7),21:(16,16),63:(2,2),341:(23,1),1023:(2,2)}
for n in expected_residues:
    sn=strong_counts(n,n)[n]; require((sn,c(n,n))==expected_residues[n], f'example residue n={n}'); print('n,s(n)%n,c(n)%n,D(n)%n:',n,sn,c(n,n),(sn-c(n,n))%n)
print('D(11)%31=',(s[11]-c(11,31))%31)
print('s(1)..s(7):',strong_counts(7,10**20)[1:])
# All connection sets, independent subset DP and reflected-path enumeration.
checks=0
for n in [3,5,7,9,11]:
    for bits in range(1<<((n-1)//2)):
        S={j if bits>>(j-1)&1 else n-j for j in range(1,(n+1)//2)}
        mu=sum({a*x%n for x in S}==S for a in range(1,n) if gcd(a,n)==1)
        N=1<<n; dp=[{} for _ in range(N)]; dp[1][0]=1
        for mask in range(N):
            for v,num in dp[mask].items():
                for d in S:
                    w=(v+d)%n
                    if not mask>>w&1:
                        nxt=dp[mask|1<<w]; nxt[w]=nxt.get(w,0)+num
        H=n*sum(dp[-1].values()); C=sum(num for v,num in dp[-1].items() if (-v)%n in S)
        rho=0; k=(n-1)//2
        for signs in product([1,-1],repeat=k):
            for first in permutations([j*sgn%n for j,sgn in enumerate(signs,1)]):
                path=first+(0,)+tuple(-x%n for x in first[::-1])
                rho+=all((b-a)%n in S for a,b in zip(path,path[1:]))
        require(H%n==0 and rho%2, f"path/reflection parity n={n}, S={S}")
        if prime(n):
            require(H%(n*mu)==0 and (C-k)%(n*mu)==0 and (rho-mu)%(2*mu)==0, f"prime circulant congruences n={n}, S={S}")
        checks+=1
print('all connection sets n=3,5,7,9,11:',checks,'DP/reflection checks')
# Audit advertised census denominator directly from immutable checkpoint outputs.
base='04-computation/experiments/tournament_pseudoprimes_20261021/'
files=['circ_small.out','circ_13_17.out','circ_19.out','circ_21.out','circ_23.out','paley19.out','paley23.out']
rows=[]
for f in files:
    out=subprocess.check_output(['git','-C',str(repo),'show','f002fc425336:'+base+f],text=True)
    rows += [line for line in out.splitlines() if line.startswith('n=')]
print('checkpoint rows:',len(rows),'unique exact rows:',len(set(rows)))
print('unique row orders:',dict(sorted(Counter(int(re.match(r'n=(\d+)',r)[1]) for r in set(rows)).items())))
keys=[]
for row in rows:
    mt=re.match(r'n=(\d+) S=\{([\d,]+)\}',row); n=int(mt[1]); S=set(map(int,mt[2].split(',')))
    keys.append((n,min(tuple(sorted(a*x%n for x in S)) for a in range(1,n) if gcd(a,n)==1)))
print('distinct multiplier classes:',len(set(keys)))
print('multiplier-class orders:',dict(sorted(Counter(n for n,S in set(keys)).items())))
for p in [3,5,7,11,13,17,19,23,29,31]:
    sq=p*p; vals=strong_counts(sq,sq)
    a=(vals[p]-c(p,sq))%sq; b=(vals[sq]-c(sq,sq))%sq
    print('prime-square check p,D(p),D(p^2) modulo p^2:',p,a,b)
    require((a==0)==(b==0), f"prime-square vanishing p={p}")
    eps=1 if pow(2,(p-1)//2,p)==1 else -1
    require((b-eps*a)%sq==0, f"prime-square lift p={p}")
# General p-power lifting: independent exact modular recurrence, not orbit data.
for p,max_r,max_m in [(3,5,4),(5,4,2),(7,3,3),(11,2,5)]:
    for r in range(1,max_r+1):
        mod=p**r; vals=strong_counts(max_m*mod,mod)
        eps=1 if pow(2,(p-1)//2,p)==1 else -1
        for m in range(1,max_m+1):
            n=m*mod; prev=n//p
            require((vals[n]-eps**m*vals[prev])%mod==0, f"s lift p={p}, r={r}, m={m}")
            require(((vals[n]-c(n,mod))-eps**m*(vals[prev]-c(prev,mod)))%mod==0, f"D lift p={p}, r={r}, m={m}")
        print('general prime-power lifting checked:',p,r,'m=1..'+str(max_m),'modulus',mod)

require(len(rows)==264 and len(set(keys))==262, "checkpoint census denominator")
# Incoming Lemma C: independent connected-graph recurrence and sharp boundary.
def connected_counts(N,mod):
    total=[pow(2,n*(n-1)//2,mod) for n in range(N+1)]
    connected=[0]*(N+1)
    for n in range(1,N+1):
        connected[n]=(total[n]-sum(comb(n-1,k-1)*connected[k]*total[n-k] for k in range(1,n)))%mod
    return connected
connected_checks=0
for p in [3,5,7,11,13,17,19,23,29,31]:
    counts=connected_counts(100+p,p)
    eps=pow(2,(p-1)//2,p)
    require((counts[p]-eps*counts[1])%p==p-1, f"connected initial obstruction p={p}")
    for m in range(2,101):
        require((counts[m+p-1]-eps*counts[m])%p==0, f"connected tail p={p}, m={m}")
        connected_checks+=1
print('incoming connected-graph tail independently checked:',connected_checks,'cases, plus 10 initial obstructions')
print("PASS: all runtime checks completed")
