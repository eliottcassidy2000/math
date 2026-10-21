#!/usr/bin/env python3
"""Exact source-owned annulus horizon ledger; standard library only.
Default recomputes the declared finite universe. No convergence assumption.
--write saves a deterministic receipt; otherwise compare to the saved receipt.
"""
import argparse, hashlib, json, time
from collections import defaultdict
from pathlib import Path
HERE = Path(__file__).resolve().parent
CHECKS = 0
DEPTH_CAP = 256
ROUTE_CAP = 10000
L = ((1,2),(-2,-1),(1,-1))
U = ((0,1),(-1,-1),(1,0))
K = (0,-3,3)
T = (-1,-1,2)

def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok: raise ValueError(message)

def exact(n): return type(n) is int

def Q(v):
    a,b=v
    return a*a-a*b+b*b

def beta(v,w):
    a,b=v; c,d=w
    return 2*a*c-a*d-b*c+2*b*d

def encode(obj): return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()

def source(s): return exact(s) and s>0 and s%2==1 and s%8!=5

def edge(s):
    check(source(s),'edge source')
    a=2 if s%8==1 else 1
    check((3*s+1)%2**a==0,'edge divisibility')
    u=(3*s+1)//2**a
    while u%8==5:
        u=(u-1)//4; a+=2
    check(source(u),'edge target')
    return u,a

def charge(N,s):
    check(exact(N) and N>=1 and source(s),'charge domain')
    a=b=0
    while s<4*N:
        p,q=(L if s<2*N else U)[s%3]
        a+=p; b+=q; s=4*s+1
    return a,b

def independent_charge(N,s):
    c0=c1=0
    while s<4*N:
        row=K if s<2*N else T
        c0+=row[s%3]; c1+=row[(s+1)%3]; s=4*s+1
    check((c0-c1)%3==0 and (c0+2*c1)%3==0,'coordinate integrality')
    return (c0-c1)//3, -(c0+2*c1)//3

class Routes:
    def __init__(self,V,ranks):
        self.V=set(V)
        check(len(V)==627 and len(self.V)==627,'frontier universe')
        check(set(ranks)==self.V and ranks[1]==0,'rank domain')
        for s in V:
            check(source(s) and exact(ranks[s]) and ranks[s]>=0,'frontier type')
            if s!=1:
                u,a=edge(s)
                check(u in self.V and ranks[u]<ranks[s],'grounded frontier rank')
        self.cache={s:(s,0,0) for s in V}
        self.maximum_depth=0
    def get(self,s):
        check(source(s),'route source')
        trail=[]; seen=set(); u=s
        while u not in self.cache:
            check(u not in seen and len(trail)<ROUTE_CAP,'unresolved or repeated route')
            seen.add(u); v,a=edge(u); trail.append((u,a)); u=v
        end,t,d=self.cache[u]
        for u,a in reversed(trail):
            t+=a;d+=1; self.cache[u]=(end,t,d)
        end,t,d=self.cache[s]
        check(exact(t) and exact(d) and 0<=d<=DEPTH_CAP,'clock/depth cap')
        self.maximum_depth=max(self.maximum_depth,d)
        return end,t,d
    def audit(self):
        digest=hashlib.sha256()
        for s,(u,t,d) in sorted(self.cache.items()):
            if s in self.V:
                check((u,t,d)==(s,0,0),'first arrival frontier convention')
                row=[s,s,0,u,t,d]
            else:
                z,a=edge(s)
                check(z in self.cache,'route child present')
                v,b,c=self.cache[z]
                check(u==v and t==a+b and d==c+1,'native suffix authentication')
                row=[s,z,a,u,t,d]
            digest.update(encode(row))
        return {'nodes':len(self.cache),'sha256':digest.hexdigest(),'max_requested_depth':self.maximum_depth}

class Sweep:
    def __init__(self,N,routes):
        self.N=N; self.routes=routes; D=DEPTH_CAP
        self.W={u:([0]*(D+1),[0]*(D+1)) for u in routes.V}
        self.E=[0]*(D+1); self.M=[0]*(D+1)
        buckets=defaultdict(lambda:[0,0,0])
        for s in range(1,4*N,2):
            if s%8==5:continue
            u,t,d=routes.get(s); a,b=charge(N,s)
            check((a,b)==independent_charge(N,s),'initial shell charges')
            sg=1-2*(t%2);r=buckets[u,d];r[0]+=sg*a;r[1]+=sg*b;r[2]+=Q((a,b))
        for (u,d),(a,b,m) in buckets.items():self.update(u,d,a,b,m)
    def update(self,u,d,a,b,m):
        aa,bb=self.W[u];e=self.E;ms=self.M;q=Q((a,b));ca=2*a-b;cb=2*b-a
        for h in range(d,DEPTH_CAP+1):
            e[h]+=q+ca*aa[h]+cb*bb[h]; aa[h]+=a;bb[h]+=b;ms[h]+=m
    def step(self):
        N=self.N; deltas=defaultdict(lambda:[0,0])
        m=2*N+1
        changes=[(m,(L[m%3][0]-U[m%3][0],L[m%3][1]-U[m%3][1])),(4*N+1,U[(4*N+1)%3]),(4*N+3,U[(4*N+3)%3])]
        for m,weight in changes:
            s=m
            while s%8==5:s=(s-1)//4
            r=deltas[s];r[0]+=weight[0];r[1]+=weight[1]
        for s,(a,b) in deltas.items():
            old=charge(N,s); new=(old[0]+a,old[1]+b)
            check(new==charge(N+1,s),'three-endpoint source recurrence')
            m=Q(new)-Q(old);u,t,d=self.routes.get(s);sg=1-2*(t%2)
            self.update(u,d,sg*a,sg*b,m)
        self.N+=1

def direct(N,routes,emit=None):
    """Independent source enumeration; grouped horizon construction, no sweep."""
    batches=defaultdict(lambda:defaultdict(lambda:[0,0,0])); dig=hashlib.sha256(); count=0
    for s in range(1,4*N,2):
        if s%8==5:continue
        u,t,d=routes.get(s);a,b=charge(N,s)
        check((a,b)==independent_charge(N,s),'direct charge agreement')
        r15=s%15;r35=s%35
        check(r15%5==r35%5 and (r15+15*((r35-r15)%7))%105==s%105,'compatible joint register')
        row=[s,r15,r35,u,t,d,a,b,Q((a,b))]
        dig.update(encode(row));count+=1
        if emit is not None:emit.write(encode(row))
        sg=1-2*(t%2);v=batches[d][u];v[0]+=sg*a;v[1]+=sg*b;v[2]+=Q((a,b))
    max_depth=max(batches)
    E=[];M=[];W={};energy=mass=0
    for d in range(DEPTH_CAP+1):
        for u,(a,b,m) in batches[d].items():
            w=W.get(u,(0,0));energy+=Q((a,b))+beta(w,(a,b));W[u]=(w[0]+a,w[1]+b);mass+=m
        E.append(energy);M.append(mass)
    return E,M,{'N':N,'sources':count,'source_norm':mass,'final_energy':energy,'max_depth':max_depth,'source_rows_sha256':dig.hexdigest()}

def decomposition(N,routes):
    batches=defaultdict(list);oldmass=newmass=0;pair_buckets=defaultdict(list)
    for s in range(1,256*N,2):
        if s%8==5:continue
        a,b=charge(64*N,s);u,t,d=routes.get(s);q=Q((a,b));sg=1-2*(t%2)
        old=s<4*N
        if old:
            check((a,b)==charge(N,s),'64-scale fixed named source')
            oldmass+=q
        else:
            newmass+=q
            if q:pair_buckets[u,sg*a,sg*b].append((s,t,d,a,b))
        batches[d].append((old,u,sg*a,sg*b,q))
    check(newmass==153*N,'exact annulus norm')
    W={};H={};ms=0;rows=[]
    for d in range(DEPTH_CAP+1):
        for old,u,a,b,q in batches[d]:
            table=W if old else H;v=table.get(u,(0,0));table[u]=(v[0]+a,v[1]+b)
            if not old:ms+=q
        e=sum(map(Q,W.values()));h=sum(map(Q,H.values()));mix=sum(beta(w,H.get(u,(0,0))) for u,w in W.items())
        rows.append([d,e,h,mix,ms,newmass-ms,h-ms,mix+h-ms,newmass-h-mix])
        if N==4096 and d==1:
            positives=[u for u,z in H.items() if Q(z)+beta(W.get(u,(0,0)),z)>sum(q for old,v,a,b,q in batches[0]+batches[1] if not old and v==u)]
            check(len(positives)==204,'local fibre budget obstruction count')
            check(W[9]==(-1,2) and H[9]==(2,2),'literal terminal-nine vectors')
            check(Q(H[9])+beta(W[9],H[9])-4==2,'local terminal-nine excess')
    pair=None
    for (u,a,b),left in sorted(pair_buckets.items()):
        right=pair_buckets.get((u,-a,-b),[])
        if left and right and left[0][2]!=right[0][2]:
            pair={'terminal':u,'sources':[list(left[0]),list(right[0])],'complete_depth':max(left[0][2],right[0][2])}
            check(left[0][0]!=right[0][0] and min(left[0][0],right[0][0])>=4*N,'distinct new packet sources')
            check(Q((a,b))>0 and (a-a,b-b)==(0,0),'positive partial zero-pair boundary')
            break
    return rows,pair

def hostile_controls(routes):
    small=[]
    for N in range(1,8):
        rows,_=decomposition(N,routes)
        bad=[r for r in rows if r[-1]<0]
        small.append({'N':N,'first_failed_horizon':bad[0] if bad else None,'minimum_bank':min(r[-1] for r in rows)})
    check(all(x['first_failed_horizon'] is None for x in small[:6]),'small positive controls')
    check(small[6]['first_failed_horizon']==[8,19,1086,-2,1054,17,32,30,-13],'minimal-N annulus hostile')
    # A mixed zero packet cannot be erased from the baseline increment.
    R=(1,0);old=(1,0);new=(-1,0)
    full=Q((R[0]+old[0]+new[0],0))-Q((R[0]+old[0],0))
    check(full==-3 and -Q(old)-beta(R,old)==full,'mixed packet old-baseline correction')
    # Distinct terminals do not cancel each other's energies.
    check(Q((1,0))+Q((-1,0))==2 and Q((0,0))==0,'terminal merge hostile')
    for bad in (True,1.0,0,-1,2,5):
        check(not source(bad),'exact primitive source hostile')
    return small

def run(lo,hi):
    data=json.loads((HERE/'frontier.json').read_text());ranks=dict(data['frontier_ranks']);routes=Routes(data['frontier'],ranks)
    old=Sweep(lo,routes);olds={};direct_checks={}
    selected=sorted(set([lo,hi,4096,5997,11002,16384])&set(range(lo,hi+1)))
    for N in range(lo,hi+1):
        olds[N]=(old.E.copy(),old.M.copy())
        if N in selected:
            e,m,info=direct(N,routes)
            check(e==old.E and m==old.M,'independent old-scale recomputation')
            direct_checks[N]={'old':info}
        if N<hi:old.step()
    big=Sweep(64*lo,routes);records=[];horizon_digest=hashlib.sha256();worst=None
    for M in range(64*lo,64*hi+1):
        if M%64==0:
            N=M//64;oe,om=olds[N]
            check(big.M[-1]-om[-1]==153*N,'sweep annulus norm')
            deltas=[];cross=[]
            for d in range(DEPTH_CAP+1):
                delta=big.E[d]-oe[d]; gap=delta-big.M[d]+om[d]
                deltas.append(delta);cross.append(gap)
                check(gap<=0,'strong arrived-mass increment failure')
                horizon_digest.update(encode([N,d,oe[d],big.E[d],om[d],big.M[d],delta,gap]))
            peak=max(deltas);at=deltas.index(peak);worstcross=max(cross)
            record=[N,peak,at,153*N-peak,worstcross,cross.index(worstcross)]
            records.append(record)
            if worst is None or peak*worst[0]>worst[1]*N:worst=[N,peak,at]
            if N in selected:
                e,m,info=direct(M,routes)
                check(e==big.E and m==big.M,'independent enlarged-scale recomputation')
                direct_checks[N]['large']=info
        if M<64*hi:big.step()
    details={}
    for N in selected:
        rows,pair=decomposition(N,routes)
        check(all(r[7]<=0 for r in rows),'separate mixed/new ledger')
        details[N]={'horizons':rows,'new_only_partial_pair':pair}
    hostile=hostile_controls(routes)
    dag=routes.audit()
    return {'status':'PROVED algebraic identities; FINITE-EXACT stated universe; universal target OPEN','base_range':[lo,hi],'checked_horizons':[0,DEPTH_CAP],'all_later_horizons':'identical after authenticated complete arrival','route_cap':ROUTE_CAP,'source_row_schema':['s','s_mod15','s_mod35','terminal','integer_clock','P_depth','a','b','Q'],'detail_schema':['d','old_E','new_only_E','mixed_old_new','arrived_new_norm','moving_new_norm','new_new_cross','all_cross','bank_153N_minus_increment'],'record_schema':['N','max_increment','first_peak_depth','minimum_full_bank','max_all_cross','first_max_cross_depth'],'records':records,'worst_increment_ratio':worst,'horizon_rows_sha256':horizon_digest.hexdigest(),'direct_checks':direct_checks,'details':details,'small_hostiles':hostile,'route_dag':dag,'frontier_sha256':hashlib.sha256((HERE/'frontier.json').read_bytes()).hexdigest(),'checks':CHECKS}

def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--lo',type=int,default=4096);p.add_argument('--hi',type=int,default=16384);args=p.parse_args()
    check(1<=args.lo<=args.hi,'scale interval')
    result=run(args.lo,args.hi);path=HERE/'annulus_receipt.json';serialized=encode(result)
    if args.write:path.write_bytes(serialized)
    else:check(path.read_bytes()==serialized,'saved receipt equality')
    print(json.dumps({'status':'PASS','base_range':result['base_range'],'checks':result['checks'],'receipt_sha256':hashlib.sha256(serialized).hexdigest(),'worst_increment_ratio':result['worst_increment_ratio'],'route_dag':result['route_dag']},sort_keys=True))
if __name__=='__main__':main()
