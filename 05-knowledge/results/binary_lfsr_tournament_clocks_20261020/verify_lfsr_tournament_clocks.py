#!/usr/bin/env python3
"""Exact, dependency-free controls for binary module graphs and source towers.
All verification gates remain active under python -O."""
from collections import Counter, defaultdict
from itertools import permutations
from math import gcd, lcm
from pathlib import Path
import json


CHECKS = 0

def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)

def deg(p): return p.bit_length() - 1

def rem(a, b):
    while a and deg(a) >= deg(b): a ^= b << (deg(a) - deg(b))
    return a

def divexact(a, b):
    q = 0
    while a and deg(a) >= deg(b):
        s = deg(a) - deg(b); q ^= 1 << s; a ^= b << s
    require(not a, 'exact check at original line 20')
    return q

def pgcd(a, b):
    while b: a, b = b, rem(a, b)
    return a

def mul(a, b):
    out = 0
    while b:
        if b & 1: out ^= a
        a <<= 1; b >>= 1
    return out

def ppow(a, e):
    out = 1
    while e:
        if e & 1: out = mul(out, a)
        a = mul(a, a); e >>= 1
    return out

def powmod(a, e, f):
    out = 1
    while e:
        if e & 1: out = rem(mul(out, a), f)
        a = rem(mul(a, a), f); e >>= 1
    return rem(out, f)

def irreducibles(n):
    out = []
    for d in range(1, n + 1):
        for p in range(1 << d, 1 << (d + 1)):
            if not any(deg(q) <= d // 2 and rem(p, q) == 0 for q in out): out.append(p)
    return out

IRR = irreducibles(10)

def factor(f):
    ans = []
    for p in IRR:
        if deg(p) > deg(f): break
        e = 0
        while rem(f, p) == 0:
            f = divexact(f, p); e += 1
        if e: ans.append((p,e))
    require(f == 1, 'exact check at original line 65')
    return ans

def root_order(p):
    require(p != 2, 'exact check at original line 69')
    x = 1
    for r in range(1, (1 << deg(p))):
        x = rem(x << 1, p)
        if x == 1: return r
    raise AssertionError(p)

ORD = {p:root_order(p) for p in IRR if p != 2}

def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]

def factor_clock(fac):
    out = 1
    for p,e in fac:
        if p != 2: out = lcm(out, ORD[p] * (1 << ((e - 1).bit_length())))
    return out

def fixed_degree(fac,m):
    a = (m & -m).bit_length() - 1; u = m >> a
    return sum(deg(p)*min(e,1 << a) for p,e in fac if p != 2 and u % ORD[p] == 0)

def cycle_formula(fac):
    L = factor_clock(fac); exact = {}; cycles = {}
    for m in divisors(L):
        exact[m] = (1 << fixed_degree(fac,m)) - sum(v for d,v in exact.items() if m%d == 0)
        require(exact[m] >= 0 and exact[m] % m == 0, 'exact check at original line 94')
        if exact[m]: cycles[m] = exact[m] // m
    return cycles

def graph_stats(successors):
    N = len(successors); done = set(); cycles = Counter(); tails = [None]*N
    for start in range(N):
        if start in done: continue
        path=[]; local={}; x=start
        while x not in done and x not in local:
            local[x]=len(path); path.append(x); x=successors[x]
        if x in local:
            j=local[x]; cycles[len(path)-j]+=1
            for v in path[j:]: tails[v]=0
            path=path[:j]
        for v in reversed(path): tails[v]=tails[successors[v]]+1
        done.update(local)
    require(all(t is not None for t in tails), 'exact check at original line 111')
    return dict(sorted(cycles.items())), dict(sorted(Counter(tails).items()))

def transition(f,s):
    n=deg(f); s <<= 1
    return s ^ f if s & (1 << n) else s

def graph_signature(fac):
    h=next((e for p,e in fac if p==2),0)
    profile=defaultdict(list)
    for p,e in fac:
        if p!=2: profile[ORD[p]].append((deg(p),e))
    # Retain enough dyadic levels to reach saturation, plus the nilpotent depth.
    K=max((e for p,e in fac if p!=2),default=1).bit_length()
    vals={r:tuple(sum(d*min(e,1 << a) for d,e in es) for a in range(K+1)) for r,es in profile.items()}
    # Canonical removal of equal terminal levels makes different degree bounds comparable.
    for r,v in vals.items():
        while len(v)>1 and v[-1]==v[-2]: v=v[:-1]
        vals[r]=v
    return h,tuple(sorted(vals.items()))


def verify_all_polynomials(maxdeg=10):
    checked=states=fixchecks=0; graph_to_profile={}; profile_to_graph={}
    for n in range(1,maxdeg+1):
        for f in range(1 << n,1 << (n+1)):
            fac=factor(f); predicted=cycle_formula(fac)
            succ=[transition(f,s) for s in range(1 << n)]
            cyc,tails=graph_stats(succ)
            require(cyc==predicted, (f,cyc,predicted))
            h=next((e for p,e in fac if p==2),0); d=n-h
            wanted={0:1 << d}; wanted.update({t:1 << (d+t-1) for t in range(1,h+1)})
            require(tails==wanted, (f,tails,wanted))
            for m in divisors(factor_clock(fac)):
                require(deg(pgcd(f,powmod(2,m,f)^1))==fixed_degree(fac,m), 'exact check at original line 145')
                fixchecks+=1
            sig=graph_signature(fac); gs=(h,tuple(cyc.items()))
            require(graph_to_profile.setdefault(gs,sig)==sig, 'exact check at original line 148')
            require(profile_to_graph.setdefault(sig,gs)==gs, 'exact check at original line 149')
            checked+=1; states+=1 << n
    return dict(max_degree=maxdeg,polynomials=checked,states_enumerated=states,gcd_fixed_point_checks=fixchecks,graph_classes=len(graph_to_profile))


def tournament_reps(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    def canon(mask):
        arcs={(i,j) if (mask>>k)&1 else (j,i) for k,(i,j) in enumerate(pairs)}
        return min(sum((1 << k) for k,(i,j) in enumerate(pairs) if (p[i],p[j]) in arcs) for p in permutations(range(n)))
    return sorted({canon(mask) for mask in range(1 << len(pairs))}),canon

def source_extend(n,mask,canon):
    # New vertex n beats all older vertices, hence every new ordered pair (i,n) has bit zero.
    oldpairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    arcs={pair for k,pair in enumerate(oldpairs) if (mask>>k)&1}
    newpairs=[(i,j) for i in range(n+1) for j in range(i+1,n+1)]
    newmask=sum((1 << k) for k,pair in enumerate(newpairs) if pair in arcs)
    return canon(newmask)

def gf2_rank(rows):
    piv={}
    for x in rows:
        while x:
            j=deg(x)
            if j in piv: x^=piv[j]
            else: piv[j]=x;break
    return len(piv)

def matapply(columns,s):
    out=0
    while s:
        low=s&-s; out^=columns[low.bit_length()-1];s^=low
    return out

def matmul(A,B): return [matapply(A,b) for b in B]

def verify_tournament_towers(Nmax=5):
    reps=[]; canons=[]
    for n in range(Nmax+1):
        r,c=tournament_reps(n); reps.append(r);canons.append(c)
    t=[len(r) for r in reps]; require(t==[1,1,1,2,4,12], 'exact check at original line 190')
    summary=[]
    for N in range(Nmax+1):
        labels=[(n,m) for n in range(N+1) for m in reps[n]]; index={z:i for i,z in enumerate(labels)}; D=len(labels)
        S=[0]*D
        for (n,m),i in index.items():
            if n<N: S[i]=1 << index[n+1,source_extend(n,m,canons[n+1])]
        U=[(1 << i)^S[i] for i in range(D)]
        k=[t[r]-(t[r-1] if r else 0) for r in range(N+1)]
        fixed={}; power=U
        amax=N.bit_length()
        for a in range(amax+1):
            pred=sum(k[r]*min(1 << a,N-r+1) for r in range(N+1))
            require(pred==sum(t[max(0,N-(1 << a)+1):N+1]), 'exact check at original line 203')
            rank=gf2_rank([power[i]^(1 << i) for i in range(D)])
            require(D-rank==pred, (N,a,D-rank,pred))
            fixed[1 << a]=1 << pred
            power=matmul(power,power)
        cycles={1:fixed[1]}
        for a in range(1,amax+1):
            cycles[1 << a]=(fixed[1 << a]-fixed[1 << (a-1)])//(1 << a)
        cycles={l:c for l,c in cycles.items() if c}
        require(sum(l*c for l,c in cycles.items())==1 << D, 'exact check at original line 212')
        if N<=4:
            direct,tail=graph_stats([matapply(U,s) for s in range(1 << D)])
            require(direct==cycles and tail=={0:1 << D}, 'exact check at original line 215')
        # Distinguished empty tournament: coefficients equal Pascal row on source chain.
        cur=1 << index[0,0]
        for time in range(2*(1 << amax)):
            wanted=0;chain=0
            for j in range(N+1):
                if j & time == j: wanted ^= 1 << index[j,chain]
                if j<N: chain=source_extend(j,chain,canons[j+1])
            require(cur==wanted, (N,time,cur,wanted))
            cur=matapply(U,cur)
        summary.append(dict(N=N,dimension=D,source_free_counts=k,order=1 << amax,cycles=cycles,direct_state_enumeration=N<=4))
    return dict(tournament_counts=t,towers=summary)


def hostile_partition():
    p,q=0b1011,0b1101
    f=mul(ppow(p,2),ppow(q,4));g=mul(ppow(p,3),ppow(q,3))
    fa=[(p,2),(q,4)];ga=[(p,3),(q,3)]
    require(f!=g and deg(f)==deg(g)==18, 'exact check at original line 233')
    require(graph_signature(fa)==graph_signature(ga), 'exact check at original line 234')
    require(cycle_formula(fa)==cycle_formula(ga), 'exact check at original line 235')
    for m in range(1,113):
        require(deg(pgcd(f,powmod(2,m,f)^1))==deg(pgcd(g,powmod(2,m,g)^1)), 'exact check at original line 237')
    # First two-entry positive partition collision, mod permutation, at fixed root order.
    seen={};first=None
    for total in range(2,13):
        for a in range(1,total//2+1):
            b=total-a;profile=tuple(min(a,1 << j)+min(b,1 << j) for j in range(5))
            if profile in seen and seen[profile]!=(a,b):
                first=dict(total=total,partitions=[seen[profile],(a,b)]);break
            seen[profile]=(a,b)
        if first:break
    require(first==dict(total=6,partitions=[(2,4),(3,3)]), 'exact check at original line 247')
    direct_counts = []
    for polynomial in (f,g):
        cyc,tails=graph_stats([transition(polynomial,s) for s in range(1 << 18)])
        require(cyc==cycle_formula(fa) and tails=={0:1 << 18}, 'exact check at original line 251')
        direct_counts.append(cyc)
    # The report separates the two-entry enumeration from the proof of global minimality.
    return dict(p=p,q=q,f=f,g=g,degree=18,direct_states_enumerated=2*(1 << 18),cycles=cycle_formula(fa),profile=graph_signature(fa),first_two_entry_collision=first)


def verify_noncyclic_hostiles():
    I=[1,2];J=[1,3]
    require(graph_stats([matapply(I,s) for s in range(4)])[0]=={1:4}, 'exact check at original line 259')
    require(graph_stats([matapply(J,s) for s in range(4)])[0]=={1:2,2:1}, 'exact check at original line 260')
    def jordan(parts):
        columns=[];offset=0
        for e in parts:
            columns += [(1 << (offset+j)) ^ ((1 << (offset+j-1)) if j else 0) for j in range(e)]
            offset += e
        return columns
    A=jordan([2,4]);B=jordan([3,3])
    ca=graph_stats([matapply(A,s) for s in range(64)])[0]
    cb=graph_stats([matapply(B,s) for s in range(64)])[0]
    require(ca==cb=={1:4,2:6,4:12}, 'exact check at original line 270')
    return dict(identity_2={1:4},jordan_2={1:2,2:1},jordan_2_4_vs_3_3=ca)


def main():
    report={"status":"FINITE-EXACT","binary_companion_census":verify_all_polynomials(),"tournament_source_towers":verify_tournament_towers(),"primary_partition_hostile":hostile_partition(),"noncyclic_same_characteristic_hostile":verify_noncyclic_hostiles()}
    report['checks'] = CHECKS
    out=Path(__file__).with_name('verification.json');out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
