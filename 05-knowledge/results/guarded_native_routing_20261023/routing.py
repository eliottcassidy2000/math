#!/usr/bin/env python3
"""Exact shared-phase splice routes, primitive expansion and congestion.

The phase modulus m may encode an additional authenticated count residue;
native cells have odd order q dividing m. No electrical guard is substituted.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import json

CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def walks(r, m, beta):
    for k in range(beta, r+1, m):
        for ones in combinations(range(r), k):
            bits = [0]*r
            for j in ones: bits[j] = 1
            yield from_bits(bits, m)


def from_bits(bits, m, alpha=0):
    a = [alpha]
    for d in bits: a.append((a[-1]+d) % m)
    return tuple(a)


def bits(a, m):
    d = tuple((v-u) % m for u,v in zip(a,a[1:]))
    need(all(x in (0,1) for x in d), 'actual binary phase increments')
    return d


def bridge(u, v, m):
    L = m-1; delta = (v-u) % m
    return from_bits((1,)*delta+(0,)*(L-delta), m, u)


def splice(x, y, i, m):
    L = m-1
    return y[:i+1]+bridge(y[i],x[i+L],m)[1:-1]+x[i+L:]


def route(x, y, m, audit=True):
    """Actual phase paths; stages are retained even when a step is idle."""
    r=len(x)-1;L=m-1
    if r <= m:
        return [x,y]
    zz=[x]+[splice(x,y,i,m) for i in range(r-L+1)]+[y]
    if audit:
        for i,z in enumerate(zz[1:-1]):
            w=splice(y,x,i,m)
            dx=bits(x,m)[i:i+L];dy=bits(y,m)[i:i+L]
            rx=w[:i+1]+from_bits(dx,m,w[i])[1:]+z[i+L+1:]
            ry=z[:i+1]+from_bits(dy,m,z[i])[1:]+w[i+L+1:]
            need(rx==x and ry==y, 'complementary-splice encoder reconstructs both actual endpoints')
            need(w[0]==x[0] and w[-1]==x[-1], 'encoding remains in same endpoint carrier')
            bits(w,m)
    return zz


def sort_segment(a, start, length, m):
    """Sort the named difference interval, retaining the whole state."""
    out=[a]; d=list(bits(a,m))
    while True:
        j=next((j for j in range(start,start+length-1) if d[j:j+2]==[0,1]),None)
        if j is None: return out
        d[j],d[j+1]=1,0
        out.append(from_bits(d,m,a[0]))


def expand_edge(a, b, m):
    """Expand one block update into adjacent swaps or an m-run flip."""
    if a==b:return [a]
    r=len(a)-1
    changed=[i for i in range(1,r) if a[i]!=b[i]]
    need(changed and max(changed)-min(changed)+1<=m-1, 'contiguous repair support')
    if r<=m:
        start=0;length=r
    else:
        start=min(min(changed)-1,r-m);length=m
    da,db=bits(a,m),bits(b,m)
    need(da[:start]==db[:start] and da[start+length:]==db[start+length:], 'block exterior unchanged')
    ka=sum(da[start:start+length]);kb=sum(db[start:start+length])
    if ka!=kb:
        need(length==m and {ka,kb}=={0,m}, 'only whole-turn change inside one repair block')
        return [a,b]
    left=sort_segment(a,start,length,m);right=sort_segment(b,start,length,m)
    need(left[-1]==right[-1], 'same integer weight gives common local sorted state')
    return left+list(reversed(right[:-1]))


def primitive_options(a, m):
    d=list(bits(a,m));r=len(d)
    out=[]
    for j in range(r-1):
        e=d[:]
        if e[j]!=e[j+1]:e[j],e[j+1]=e[j+1],e[j]
        out.append(from_bits(e,m,a[0]))
    for j in range(r-m+1):
        e=d[:]
        if sum(e[j:j+m]) in (0,m):e[j:j+m]=[1-v for v in e[j:j+m]]
        out.append(from_bits(e,m,a[0]))
    return out


def primitive_matrix(states,m):
    """Lazy uniform proposal among the actual swap/run involutions."""
    pos={a:i for i,a in enumerate(states)};P={}
    for i,a in enumerate(states):
        opts=primitive_options(a,m)
        counts=Counter(pos[b] for b in opts)
        if not opts:P[i,i]=F(1);continue
        for j,c in counts.items():P[i,j]=F(c,2*len(opts))
        P[i,i]=P.get((i,i),F(0))+F(1,2)
    need(all(P.get((j,i),0)==v for (i,j),v in P.items()), 'primitive detailed balance')
    need(all(sum(v for (u,j),v in P.items() if u==i)==1 for i in range(len(states))), 'primitive probability rows')
    return P


def block_matrix(states,m):
    r=len(states[0])-1;L=m-1;N=len(states)
    if r<=m:return {(i,j):F(1,N) for i in range(N) for j in range(N)}
    B=r-L; P=defaultdict(F)
    for start in range(B):
        groups=defaultdict(list)
        for i,a in enumerate(states):groups[a[:start+1]+a[start+L+1:]].append(i)
        for group in groups.values():
            for i,j in product(group,repeat=2):P[i,j]+=F(1,B*len(group))
    return dict(P)


def bounds(r,m):
    if r<m:
        J=2*(r*r//4)
        return max(1,2*(r-1)*(2**r)*J)
    L=m-1;C=max(2,comb(m,m//2));J=2*(m*m//4)
    # For r=m the complete block has a direct route; the larger formula is valid.
    return 2*(2*r-m)*(r-m+3)*m*(m+1)*J*C*C*4**L


def host_arc(q,t,x,y):
    i,u=divmod(x,q);j,v=divmod(y,q)
    if i==j:return (v-u)%q in range(1,(q+1)//2)
    if (j-i)%t in range(1,(t+1)//2):return (v-u)%q in (1,2)
    return (u-v)%q not in (1,2)


def native_check(a,q,m):
    r=len(a)-1;t=r+1
    need(t%2==1 and q%2==1 and m%q==0,'native regular host parameters')
    p=tuple(i*q+(v+j)%q for i,v in enumerate(a) for j in range(q))
    need(len(set(p))==len(p),'all actual owners retained exactly once')
    need(all(host_arc(q,t,u,v) for u,v in zip(p,p[1:])),'all native Hamiltonian arcs')
    # Converse is checked on actual reversed owners, at reversed integer clocks.
    rev=p[::-1]
    need(all(host_arc(q,t,v,u) for u,v in zip(rev,rev[1:])),'literal converse paths')
    return p


def one_case(q,r,m,beta,full=True):
    states=list(walks(r,m,beta));pos={a:i for i,a in enumerate(states)};N=len(states)
    need(N>0,'nonempty fixed-source endpoint family')
    actual={a:native_check(a,q,m) for a in states}
    P=primitive_matrix(states,m);PB=block_matrix(states,m)
    for (i,j),prob in P.items():
        a,b=states[i],states[j]
        changed=sum((u-v)%q!=0 for u,v in zip(a,b))
        need(sum(u!=v for u,v in zip(actual[a],actual[b]))==q*changed,
             'literal changed-owner work equals q times changed-cell work')
        need(actual[a][:q]==actual[b][:q] and actual[a][-q:]==actual[b][-q:],
             'every clock in both pinned endpoint cells is retained')
    for a in (states[0],states[-1]):
        for alpha in range(q):
            lifted=native_check(tuple((v+alpha)%m for v in a),q,m)
            expected=tuple((u//q)*q+(u%q+alpha)%q for u in actual[a])
            need(lifted==expected,'actual common-phase labels and kernel lifts retained')
    load=Counter();blockload=Counter();max_events=0;max_work=0
    pair_ids=product(range(N),repeat=2) if full else [(0,N-1),(N-1,0),(N//2,0),(0,N//2)]
    pairs=0
    for ix,iy in pair_ids:
        x,y=states[ix],states[iy];rr=route(x,y,m)
        expanded=[x]
        for a,b in zip(rr,rr[1:]):
            need(PB.get((pos[a],pos[b]),0)>0,'every repair is an actual block heat-bath edge')
            ee=expand_edge(a,b,m)
            expanded.extend(ee[1:])
        need(expanded[0]==x and expanded[-1]==y,'route endpoints exactly retained')
        for a,b in zip(expanded,expanded[1:]):
            need(P.get((pos[a],pos[b]),0)>0,'every expanded edge is an inherited legal primitive')
            need(a[0]==b[0] and a[-1]==b[-1],'actual endpoint phase pins preserved')
            weight_delta=sum(bits(b,m))-sum(bits(a,m))
            need(weight_delta in (-m,0,m),'integer displacement updated explicitly')
            changed=sum((u-v)%q!=0 for u,v in zip(a,b))
            need(changed<=max(1,(q-1)*(m//q)),'sharp primitive actual-cell budget')
        events=len(expanded)-1
        work=sum(sum((u-v)%q!=0 for u,v in zip(a,b)) for a,b in zip(expanded,expanded[1:]))
        max_events=max(max_events,events);max_work=max(max_work,work)
        for a,b in zip(expanded,expanded[1:]):
            if a!=b:load[pos[a],pos[b]]+=events
        for a,b in zip(rr,rr[1:]):
            if a!=b:blockload[pos[a],pos[b]]+=len(rr)-1
        pairs+=1
    congestion=max((F(c,N)/P[e] for e,c in load.items()),default=F(0))
    if full:need(congestion<=bounds(r,m),'all-pair native congestion theorem')
    J=2*(m*m//4)
    event_bound=2*(r*r//4) if r<m else (r-m+3)*J
    need(max_events<=event_bound,'integer event-clock bound')
    out=dict(q=q,modulus=m,r=r,beta=beta,states=N,pairs=pairs,
             full_pair_universe=full,maximum_events=max_events,maximum_phase_work=max_work,
             maximum_changed_owner_work=q*max_work,
             certified_gap_denominator=bounds(r,m))
    if full:out['exact_routed_congestion']=str(congestion)
    return out


def transfer_counts(r,m,beta):
    table=[[0]*m for _ in range(r+1)];table[r][beta]=1
    for i in range(r-1,-1,-1):
        table[i]=[table[i+1][a]+table[i+1][(a+1)%m] for a in range(m)]
    return table


def unrank(rank,table,m):
    a=0;out=[a]
    for row in table[1:]:
        zero=row[a]
        if rank>=zero:rank-=zero;a=(a+1)%m
        out.append(a)
    need(rank==0,'exact uniform transfer rank exhausted')
    return tuple(out)


def rank_path(a,table,m):
    rank=0
    for i,(u,v) in enumerate(zip(a,a[1:])):
        if v==(u+1)%m:rank+=table[i+1][u]
    return rank


def main():
    cases=[]
    for q,rs in [(3,(2,4,6,8)),(5,(4,6,8)),(7,(6,8))]:
        for r in rs:
            for beta in range(min(q,r+1)):
                cases.append(one_case(q,r,q,beta))
    cases.append(one_case(9,10,9,0))
    # The actual native cells are ternary. The larger state retains count mod 2.
    for beta in range(6):cases.append(one_case(3,8,6,beta))
    # Mod-eight count guard: modulus lcm(3,8)=24, not an electrical surrogate.
    cases.append(one_case(3,26,24,0,full=False))
    transfer=[]
    for r,m,beta in [(8,3,0),(8,6,2),(128,3,1),(1024,3,0),(1024,24,0)]:
        table=transfer_counts(r,m,beta);N=table[0][0]
        need(N==sum(comb(r,k) for k in range(beta,r+1,m)),'transfer counts agree with binomial residue sum')
        ranks=range(N) if r<=8 else sorted({0,1,N//2,N-2,N-1})
        for rank in ranks:
            a=unrank(rank,table,m)
            need(rank_path(a,table,m)==rank,'lossless rank/unrank with complete integer clock')
            native_check(a,3,m)
        transfer.append(dict(r=r,modulus=m,beta=beta,count_bits=N.bit_length(),ranks_checked=len(ranks)))
    # Constant-width guard with no bounded connector: only the two constant words.
    hostiles=[]
    for r in (6,12,30,60):
        x=from_bits([0]*r,3);y=from_bits([1]*r,3)
        changed=sum(a!=b for a,b in zip(x,y))
        need(changed==2*r//3,'two-state accepted language forces unbounded native support')
        hostiles.append(dict(r=r,accepted_words=2,minimum_phase_budget=changed,
                             endpoint_phases=[x[0],x[-1]],one_event_clock_impossible=changed>2))
    return dict(status='FINITE-EXACT PASS',cases=cases,transfer=transfer,hostiles=hostiles,checks=CHECKS)


if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
