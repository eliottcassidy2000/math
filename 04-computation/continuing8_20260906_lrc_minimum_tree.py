#!/usr/bin/env python3
"""Complete minimum-tree grid consumer. Exact integer NumPy operations only."""
from pathlib import Path
from math import gcd
from itertools import combinations,combinations_with_replacement
from hashlib import sha256
import json,sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
ROOT=HERE.parents[1] if HERE.parent.name=='04-computation' else Path.cwd()
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
PROFILE_SHA='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
BASELINE_SHA='c67f5e98eff3fe406a4416057c6063095290330a50f039e731ced0fc2ca4657a'
GATES=0;INF=10**10
def need(ok,why):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(why)
def canonical(v):return json.dumps(v,separators=(',',':'),sort_keys=True).encode()
def atlas_sum(s):
    p=2
    while p*p<=s:
        k=0
        while s%p==0:s//=p;k+=1
        if k and (p%3!=2 or k>2):return False
        p+=1
    return s==1 or s%3==2
def geometry(p,q):
    L=14*p*q
    A=[(max(0,14*k*q-q),min(L,14*k*q+q)) for k in range(p+1)]
    B=[(max(0,14*k*p-p),min(L,14*k*p+p)) for k in range(q+1)]
    i=j=0;out=[]
    while i<len(A) and j<len(B):
        a,b=max(A[i][0],B[j][0]),min(A[i][1],B[j][1])
        if a<b:out.append((a,b))
        if A[i][1]<B[j][1]:i+=1
        elif A[i][1]>B[j][1]:j+=1
        else:i+=1;j+=1
    return L,[(out[-1][0]-L,out[0][1])]+out[1:-1]
def capacity(n,L,intervals):
    base=0;events={};cur=0
    for a,b in intervals:
        h,r=divmod(n*(b-a),L);base+=h
        A,B=n*a%L,n*b%L
        if not r:events.setdefault(A,[0,0,0])[2]+=1
        else:
            events.setdefault(A,[0,0,0])[0]+=1
            events.setdefault(B,[0,0,0])[1]+=1
            cur+=A>B
    best=INF
    for wall,(starts,ends,holes) in sorted(events.items()):
        best=min(best,cur-ends-holes)
        cur+=starts-ends;best=min(best,cur)
    return base+best
def literal_capacity(n,L,intervals):
    walls=sorted({n*x%L for ab in intervals for x in ab})
    phases=[2*x for x in walls]+[(walls[i]+(walls[(i+1)%len(walls)]+(L if i+1==len(walls) else 0)))%(2*L) for i in range(len(walls))]
    return min(sum(-((-(2*n*b-r))//(2*L))-(2*n*a-r)//(2*L)-1 for a,b in intervals) for r in phases)
def kruskal(S,weights):
    root=list(range(7));total=0;edges=[]
    def find(i):
        while root[i]!=i:i=root[i]
        return i
    for w,i,j in sorted((weights[tuple(sorted((S[i],S[j])))][0],i,j) for i,j in combinations(range(7),2)):
        a,b=find(i),find(j)
        if a!=b:root[a]=b;total+=w;edges.append([i,j,w])
    return total,edges
def main():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    need(sha256(raw).hexdigest()==PROFILE_SHA,'complete inherited profile raw pin')
    J=json.loads(raw);P={int(k):{(c,tuple(w)) for c,w in row['profiles']} for k,row in J['levels'].items()}
    vals=sorted({c for c,w in P[6]});need(vals==J['levels']['6']['gcds'],'complete sorted sheet alphabet')
    baseline=json.loads((ROOT/'05-knowledge/results/continuing7_20260906_lrc_residue_domains_certificate.json').read_bytes())['new_scales']
    need(sha256(canonical(baseline)).hexdigest()==BASELINE_SHA,'inherited complete necessary scale array')
    need(len(baseline)==8195 and max(baseline)==14868,'declared baseline scope')
    upper=sorted((t for t in baseline if t>=12000),reverse=True)
    targets=upper+[7200];need(len(upper)==549 and 7200 in baseline,'whole upper range plus hostile')
    sub={}
    for k,rows in P.items():
        for c,w in rows:
            for n in range(k+1):sub.setdefault((k,c,n),set()).update(tuple(w[i] for i in I) for I in combinations(range(k),n))
    indexsets={n:[(I,tuple(j for j in range(n) if j not in I)) for k in range(1,min(n,6)+1) for I in combinations(range(n),k)] for n in range(1,8)}
    def valid(S):
        n=len(S)
        for I,K in indexsets[n]:
            c=gcd(*(S[i] for i in I));w=tuple(sorted(gcd(c,S[j]) for j in K))
            if w not in sub.get((7-len(I),c,len(K)),()):return False
        return n<7 or gcd(*S)==1
    domains={};domain_cert=[]
    def words_for(D):
        if D in domains:return domains[D]
        words=[];layers=[0]*8
        def visit(S,j):
            layers[len(S)]+=1
            if len(S)==7:words.append(S);return
            for k in range(j,len(D)):
                T=S+(D[k],)
                if valid(T):visit(T,k)
        visit((),0)
        need(bool(words),'nonempty full projected domain')
        # Literal unpruned completeness controls for small alphabets.
        if len(D)<=9:
            brute=[]
            for S in combinations_with_replacement(D,7):
                if gcd(*S)!=1:continue
                ok=True
                for I,K in indexsets[7]:
                    c=gcd(*(S[i] for i in I))
                    if (c,tuple(sorted(gcd(c,S[j]) for j in K))) not in P[7-len(I)]:ok=False;break
                if ok:brute.append(S)
            need(brute==words,'unpruned small-alphabet completeness')
        domains[D]=np.asarray(words,dtype=np.int64)
        domain_cert.append(dict(domain=list(D),word_count=len(words),layers=layers,words_sha256=sha256(canonical(words)).hexdigest()))
        return domains[D]
    AT=[(p,s-p,*geometry(p,s-p)) for s in range(3,357) if atlas_sum(s) for p in range(1,(s+1)//2) if gcd(p,s-p)==1]
    need(len(AT)==5855,'complete coprime strict atlas')
    # Independent all-wall literal counts, including resonance and the near origin wrap.
    controls=[]
    for p,q in [(1,1),(1,13),(1,112),(1,355),(11,263),(23,323),(33,278),(66,289)]:
        L,ints=geometry(p,q)
        for n in [1,2,7,13,14,28,29,720,968,2485]:
            c=capacity(n,L,ints);literal=literal_capacity(n,L,ints)
            need(c==literal,'strict grid sweep equals literal wall and chamber count')
            need(0<=c<=n,'intersection count bounds')
            controls.append([p,q,n,c])
    records=[];removed=[]
    for t in targets:
        D=tuple(d for d in vals if t%d==0);S=words_for(D)
        weights={pair:(INF,) for pair in combinations_with_replacement(D,2)};events=0;compatible=0
        for e in sorted({gcd(*pair) for pair in weights}):
            n=t//e
            for p,q,L,ints in AT:
                pair=tuple(sorted((e*gcd(n,p),e*gcd(n,q))))
                if pair not in weights:continue
                compatible+=1;sep=e*sum((n*(b-a)-1)//L for a,b in ints)
                # A lower bound cannot beat the incumbent minimum if already above it.
                if sep>=weights[pair][0]:continue
                c=e*capacity(n,L,ints);events+=1
                need(sep<=c<=t,'component bound and actual all-translate count')
                if c<weights[pair][0]:weights[pair]=(c,e,p,q,sep)
        W=np.full((91,91),INF,dtype=np.int64)
        for (a,b),v in weights.items():W[a,b]=W[b,a]=v[0]
        costs=W[S[:,:,None],S[:,None,:]]
        used=np.zeros_like(S,dtype=bool);dist=np.full_like(S,INF);dist[:,0]=0
        rows=np.arange(len(S));total=np.zeros(len(S),dtype=np.int64)
        for step in range(7):
            j=np.argmin(np.where(used,10*INF,dist),axis=1)
            total+=dist[rows,j];used[rows,j]=True;dist=np.minimum(dist,costs[rows,j,:])
        E=(S*((t+7*S-1)//(7*S))).sum(axis=1)-t
        need(bool(np.all(E>=0)),'seven-tail marginal excess is nonnegative')
        need(bool(np.all(E==(S*((-(t//S))%7)).sum(axis=1)//7)),'native ceiling equals exact residue cost')
        margin=total-E;worst=int(np.argmin(margin));bad=np.flatnonzero(margin<=0)
        owner=S[worst].tolist();tree,edges=kruskal(owner,weights)
        need(tree==int(total[worst]),'independent Kruskal versus vector Prim owner')
        need(valid(tuple(owner)),'minimum owner passes every profile')
        record=dict(t=t,domain=list(D),word_count=len(S),event_evaluations=events,compatible_atlas_edges=compatible,
          minimum_margin=int(margin[worst]),owner=owner,owner_E=int(E[worst]),owner_tree=tree,owner_edges=edges,
          survivor_count=len(bad),survivors=[[S[j].tolist(),int(E[j]),int(total[j])] for j in bad],
          weights=[[list(k),list(v)] for k,v in sorted(weights.items())])
        records.append(record)
        if len(bad)==0:removed.append(t)
        print('CLOCK',t,'WORDS',len(S),'MIN_MARGIN',int(margin[worst]),'SURVIVORS',len(bad),flush=True)
    need(removed==upper,'whole declared upper range closes and hostile7200 survives')
    need(records[-1]['survivor_count']==15 and records[-1]['minimum_margin']==-74,'exact hostile boundary, not unsafe physical rows')
    new=[t for t in baseline if t not in set(removed)]
    cert=dict(status='PROVED minimum-tree consumer; FINITE-EXACT selected complete domain',profile_sha256=PROFILE_SHA,baseline_semantic_sha256=BASELINE_SHA,
      atlas_count=len(AT),domains=sorted(domain_cert,key=lambda d:d['domain']),clocks=records,controls=controls,removed_scales=removed,new_scales=new,
      new_scales_sha256=sha256(canonical(new)).hexdigest(),scope='Primitive thirteen distinct speeds; selected six-component gcd t; actual strict graph on complementary seven labels connected; all inherited126profiles. Weak clearance only.')
    data=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
    print('COMPLETE',len(upper),'upper clocks excluded;',len(new),'necessary clocks remain; maximum',max(new))
    print('DOMAIN_TOTAL',len(domains),'WORDS',sum(len(v) for v in domains.values()),'MIN_UPPER_MARGIN',min(r['minimum_margin'] for r in records[:-1]))
    print('CERTIFICATE_SHA256',sha256(data).hexdigest())
    print('NEW_SCALES_SHA256',sha256(canonical(new)).hexdigest())
    print('PASS',GATES,'always-active exact gates; actual LF')
if __name__=='__main__':main()
