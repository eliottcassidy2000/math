"""Exact source-owned finite-horizon observer and sunflower ownership controls."""
from pathlib import Path
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product, combinations
from collections import defaultdict
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
FRONTIER=HERE.parent/'session18_native_annulus'/'frontier.json'
CHECKS=0

def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise ValueError(msg)

def integer(x): return type(x) is int

def Q(v):
    a,b=v
    return a*a-a*b+b*b

def add(v,w): return tuple(a+b for a,b in zip(v,w))

def scale(c,v): return tuple(c*a for a in v)

def valid_source(s): return integer(s) and s>0 and s%2==1 and s%8!=5

def edge(s):
    need(valid_source(s),'native source required')
    a=2 if s%8==1 else 1
    x=(3*s+1)//(1<<a)
    while x%8==5: x=(x-1)//4; a+=2
    need(valid_source(x),'native target required')
    return x,a

def charge(N,s):
    need(integer(N) and N>=1 and valid_source(s),'charge domain')
    out=(0,0)
    while s<4*N:
        out=add(out,(((1,2),(-2,-1),(1,-1)) if s<2*N else
                     ((0,1),(-1,-1),(1,0)))[s%3])
        s=4*s+1
    return out

def charge_check(N,s):
    c0=c1=0
    while s<4*N:
        r=(0,-3,3) if s<2*N else (-1,-1,2)
        c0+=r[s%3]; c1+=r[(s+1)%3]; s=4*s+1
    need((c0-c1)%3==0 and (c0+2*c1)%3==0,'independent charge integrality')
    return ((c0-c1)//3,-(c0+2*c1)//3)

@dataclass(frozen=True)
class Row:
    source: int
    old: bool
    terminal: int
    depth: int
    clock: int
    charge: tuple
    route: tuple

    def signed(self):
        return scale((-1)**self.clock,self.charge)

def row(N,s,V,cap=1000):
    need(integer(N) and N>=1 and valid_source(s) and s<256*N,'scale/source')
    route=[s]; clock=0; seen={s}
    while route[-1] not in V:
        need(len(route)<=cap,'finite route cap reached')
        x,a=edge(route[-1])
        need(x not in seen,'repeated unresolved route')
        seen.add(x);route.append(x);clock+=a
    r=Row(s,s<4*N,route[-1],len(route)-1,clock,charge(64*N,s),tuple(route))
    audit_row(N,r,V)
    return r

def audit_row(N,r,V):
    need(type(r) is Row and integer(N) and N>=1,'typed source row')
    need(valid_source(r.source) and r.source<256*N and type(r.old) is bool,'source mark')
    need(r.old==(r.source<4*N),'old/new source ownership')
    need(integer(r.terminal) and integer(r.depth) and integer(r.clock),'exact clocks')
    need(type(r.route) is tuple and len(r.route)>=1 and
         all(valid_source(x) for x in r.route),'exact route')
    need(r.route[0]==r.source and r.route[-1]==r.terminal and
         r.depth==len(r.route)-1,'source-bound endpoints')
    need(r.terminal in V and all(x not in V for x in r.route[:-1]),'strict first arrival')
    total=0
    for a,b in zip(r.route,r.route[1:]):
        x,t=edge(a);need(x==b,'actual native edge');total+=t
    need(total==r.clock,'actual integer clock')
    need(type(r.charge) is tuple and len(r.charge)==2 and
         all(integer(x) for x in r.charge),'exact charge')
    need(r.charge==charge(64*N,r.source)==charge_check(64*N,r.source),'original signed charge')
    return r

def rank(vectors):
    if not vectors: return 0
    a=[list(map(F,v)) for v in vectors];p=0
    for c in range(len(a[0])):
        j=next((j for j in range(p,len(a)) if a[j][c]),None)
        if j is None:continue
        a[p],a[j]=a[j],a[p]
        s=a[p][c];a[p]=[x/s for x in a[p]]
        for j in range(len(a)):
            if j!=p and a[j][c]:
                s=a[j][c];a[j]=[x-s*y for x,y in zip(a[j],a[p])]
        p+=1
        if p==len(a):break
    return p

def coordinates(basis,v):
    if not basis:
        need(not any(v),'zero-span vector')
        return ()
    k=len(basis)
    a=[[F(b[i]) for b in basis]+[F(v[i])] for i in range(len(v))]
    p=0
    for c in range(k):
        j=next(j for j in range(p,len(a)) if a[j][c])
        a[p],a[j]=a[j],a[p];s=a[p][c];a[p]=[x/s for x in a[p]]
        for j in range(len(a)):
            if j!=p:
                s=a[j][c];a[j]=[x-s*y for x,y in zip(a[j],a[p])]
        p+=1
    need(all(not any(r[:k]) and not r[-1] for r in a[k:]),'outside common feature span')
    return tuple(a[i][-1] for i in range(k))

def compress(N,rows,D,V,weights=None):
    need(type(rows) is tuple and all(type(r) is Row for r in rows),'typed row family')
    need(len({r.source for r in rows})==len(rows),'source cannot be duplicated')
    need(integer(D) and D>=0,'finite horizon')
    for r in rows: audit_row(N,r,V)
    if weights is None:weights=(1,)*len(rows)
    need(type(weights) is tuple and len(weights)==len(rows) and
         all(type(w) in (int,F) for w in weights),'exact coefficient vector')
    bins=defaultdict(list)
    for r,w in zip(rows,weights):
        if r.depth>D:
            if not r.old: bins[('moving',)].append((r.source,(Q(r.charge),),w))
        else:
            v=r.signed()
            bins[('old' if r.old else 'new',r.terminal,r.depth)].append(
                (r.source,v if r.old else v+(Q(v),),w))
    result=[]
    for key,items in sorted(bins.items()):
        basis=[]; owners=[];total=(0,)*len(items[0][1])
        for source,v,w in items:
            if rank(basis+[v])>len(basis):basis.append(v);owners.append(source)
            total=add(total,scale(w,v))
        c=coordinates(basis,total)
        recovered=tuple(sum(c[j]*basis[j][i] for j in range(len(basis)))
                        for i in range(len(total)))
        need(recovered==total,'common-space reconstruction')
        result.append({'key':key,'basis_sources':tuple(owners),'basis':tuple(basis),
                       'coordinates':c,'total':total})
    return tuple(result)

def ledger(rows,d,weights=None):
    if weights is None:weights=(1,)*len(rows)
    W=defaultdict(lambda:(0,0));H=defaultdict(lambda:(0,0));arrived=0;moving=0
    for r,w in zip(rows,weights):
        if r.depth<=d:
            dst=W if r.old else H
            dst[r.terminal]=add(dst[r.terminal],scale(w,r.signed()))
            if not r.old:arrived+=w*Q(r.charge)
        elif not r.old:moving+=w*Q(r.charge)
    delta=sum(Q(add(W[u],H[u]))-Q(W[u]) for u in set(W)|set(H))
    return delta,arrived,moving,delta-arrived

def decoded_ledger(packet,d):
    W=defaultdict(lambda:(0,0));H=defaultdict(lambda:(0,0));arrived=0;moving=0
    for b in packet:
        key=b['key'];v=b['total']
        if key[0]=='moving':moving+=v[0];continue
        kind,u,t=key
        if t<=d:
            dst=W if kind=='old' else H;dst[u]=add(dst[u],v[:2])
            if kind=='new':arrived+=v[2]
        elif kind=='new':moving+=v[2]
    delta=sum(Q(add(W[u],H[u]))-Q(W[u]) for u in set(W)|set(H))
    return delta,arrived,moving,delta-arrived

def sunflower(sets):
    need(type(sets) is tuple and len(sets)>=3 and
         all(type(s) is frozenset and all(integer(x) for x in s) for s in sets),
         'typed labelled set family')
    need(len(set(sets))==len(sets),'distinct petals')
    core=sets[0]&sets[1]
    need(all(a&b==core for a,b in combinations(sets,2)),'not a sunflower')
    return core,frozenset().union(*sets)

def interpolation(values):
    need(type(values) is tuple and len(values)>=3 and all(integer(x) for x in values),
         'exact interpolation alphabet tuple')
    a0=values[0];R=F(1);terms=[]
    for h in range(1,len(values)):
        prefix=values[:h]; coeff=[F(1)]
        if len(set(prefix))<h: coeff=[F(0)]
        else:
            for ai in prefix[1:]:
                out=[F(0)]*(len(coeff)+1)
                for j,c in enumerate(coeff):out[j]-=ai*c/(a0-ai);out[j+1]+=c/(a0-ai)
                coeff=out
        delta=int(all(x==a0 for x in values[h+1:]))
        terms.append(tuple(R*c*delta for c in coeff))
        ph=sum(c*values[h]**j for j,c in enumerate(coeff))
        R*=int(values[h]==a0)-ph
    q=F(int(len(set(values))==1))-R
    reconstructed=sum(sum(c*values[h]**j for j,c in enumerate(terms[h-1]))
                      for h in range(1,len(values)))
    need(q==reconstructed,'interpolated decomposition')
    return q,tuple(terms)

def serial(obj):
    if type(obj) is F:return [obj.numerator,obj.denominator]
    if type(obj) is Row:return {k:serial(getattr(obj,k)) for k in obj.__dataclass_fields__}
    if type(obj) is dict:return {k:serial(v) for k,v in obj.items()}
    if type(obj) in (tuple,list):return [serial(v) for v in obj]
    return obj

def verify():
    f=json.loads(FRONTIER.read_text());V=set(f['frontier']);ranks=dict(f['frontier_ranks'])
    need(len(V)==627 and ranks[1]==0,'inherited frontier')
    for s in V:
        if s!=1:
            u,t=edge(s);need(u in V and ranks[u]<ranks[s],'grounded frontier')
    # Complete small native scale, independently generated charges and clocks.
    N=7;rows=tuple(row(N,s,V) for s in range(1,256*N,2) if valid_source(s))
    need(sum(Q(r.charge) for r in rows if not r.old)==153*N,'exact source norm')
    D=max(r.depth for r in rows)
    packet=compress(N,rows,D,V);rank_all=sum(len(b['basis']) for b in packet)
    for d in range(D+1):
        need(ledger(rows,d)==decoded_ledger(packet,d),'all-cut native ledger')
    need(ledger(rows,8)==(1084,1054,17,30),'inherited N7 hostile')
    horizon_spaces={}
    for horizon in (0,1,8,D):
        pp=compress(N,rows,horizon,V)
        horizon_spaces[str(horizon)]=sum(len(b['basis']) for b in pp)
        for d in range(horizon+1):need(ledger(rows,d)==decoded_ledger(pp,d),'horizon prefix observer')
    # Signed rational coefficients are an algebraic observer control, not a new native source measure.
    weights=tuple(F((i%7)-3,3) for i in range(len(rows)))
    pp=compress(N,rows,D,V,weights)
    for d in range(D+1):
        need(ledger(rows,d,weights)==decoded_ledger(pp,d),'signed coefficient observer')
    ids=(33435,32955,33891,33897,33921)
    flower=tuple(row(4096,s,V) for s in ids)
    need(all(not r.old and r.terminal==955 and r.charge==(0,1) for r in flower),'native flower')
    need(flower[0].signed()==(0,1) and all(r.signed()==(0,-1) for r in flower[1:]),'actual clock signs')
    sets=tuple(frozenset((ids[0],s)) for s in ids[1:])
    core,union=sunflower(sets)
    need(core==frozenset((ids[0],)) and union==frozenset(ids),'literal sunflower supports')
    for d in range(max(r.depth for r in flower)+1):
        vec=lambda S:tuple(sum(r.signed()[j] for r in flower if r.source in S and r.depth<=d) for j in range(2))
        lhs=vec(union);rhs=add(tuple(sum(vec(S)[j] for S in sets) for j in range(2)),scale(-3,vec(core)))
        need(lhs==rhs,'one-copy shared-core identity')
    end=max(r.depth for r in flower)
    flower_packet=compress(4096,flower,end,V)
    need(sum(len(b['basis']) for b in flower_packet)==3,'native flower common rank')
    for d in range(end+1):
        need(ledger(flower,d)==decoded_ledger(flower_packet,d),'flower compressed cuts')
    for S in sets:
        need(ledger(tuple(r for r in flower if r.source in S),end)[0]==0,'completed pair zero')
    need(ledger(flower,end)==(9,5,0,4),'four-petal local norm failure')
    need(ledger(flower,10)==(9,3,2,6),'actual early boundary')
    # A single owned pair really cancels, but only at its completion cut.
    single=(flower[0],flower[1])
    need(ledger(single,22)[0]==1 and ledger(single,34)[0]==0,'completion clock retained')
    ninterp=0
    for k,D0 in ((3,4),(4,5)):
        for v in product(range(1,D0+1),repeat=k):
            q,terms=interpolation(v)
            allowed=len(set(v)) in (1,k)
            need(bool(q)==allowed and (len(set(v))!=1 or q==1),'tester support')
            need(all(len(c)<=h for h,c in enumerate(terms,1)),'low-degree common factor')
            ninterp+=1
    need(interpolation((1,2,3))[0]==1 and interpolation((1,3,2))[0]==F(-1,2),
         'amplitude support versus signed mass')
    # Common-space kernel/diagonal boundary: identical observer columns cannot label both sources.
    need(rank([(1,0,1),(1,0,1)])==1,'two source labels one observer direction')
    need(Q((1,0))+Q((-1,0))==2 and Q((0,0))==0,'norm sidecar cannot be inferred from total charge')
    # Arbitrary coefficients preserve the linear bill, not scaled-atom pair ownership.
    duplicate_columns=((1,0,1),(1,0,1))
    linear=lambda cs:tuple(sum(c*v[i] for c,v in zip(cs,duplicate_columns)) for i in range(3))
    need(linear((1,-1))==linear((0,0))==(0,0,0),'signed kernel retains only linear norm bill')
    need(sum(c*c for c in (1,-1))==2 and sum(c*c-c for c in (1,-1))==2,
         'scaled diagonal correction is not an owned cross term')
    def rejects(f):
        try:f()
        except (ValueError,TypeError):return
        need(False,'bad alias accepted')
    rejects(lambda: row(4096,True,V))
    rejects(lambda: audit_row(4096,Row(float(ids[0]),False,955,22,40,(0,1),flower[0].route),V))
    rejects(lambda: compress(4096,(flower[0],flower[0]),34,V))
    rejects(lambda: sunflower((sets[0],sets[0],sets[1])))
    rejects(lambda: interpolation((1,2.,3)))
    return {'status':'PROVED observer and ownership identities; FINITE-EXACT native controls; uniform153N OPEN',
            'checks':CHECKS,'paper_interpolation_tuples':ninterp,
            'small_native':{'N':N,'source_count':len(rows),'maximum_depth':D,
                'common_space_dimension':rank_all,'bin_count':len(packet),
                'horizon_common_dimensions':horizon_spaces,
                'horizon8':ledger(rows,8),'source_rows_sha256':hashlib.sha256(json.dumps(serial(rows),sort_keys=True).encode()).hexdigest()},
            'native_sunflower':{'N':4096,'rows':flower,'supports':[sorted(s) for s in sets],
                'core':sorted(core),'union':sorted(union),'final_depth':end,
                'common_space_dimension':3,'compressed_bins':flower_packet,
                'final_ledger':ledger(flower,end),'depth10_ledger':ledger(flower,10)},
            'amplitude_hostile':{'123':interpolation((1,2,3))[0],'132':interpolation((1,3,2))[0]},
            'frontier_sha256':hashlib.sha256(FRONTIER.read_bytes().replace(b'\r\n',b'\n')).hexdigest()}

if __name__=='__main__':
    output=json.dumps(serial(verify()),indent=2,sort_keys=True)+'\n'
    if '--write' in sys.argv: (HERE/'receipt.json').write_text(output,encoding='utf-8',newline='\n')
    else:
        need((HERE/'receipt.json').read_text(encoding='utf-8')==output,'saved receipt mismatch')
    print(json.dumps({'checks':json.loads(output)['checks'],'LF_SHA256':hashlib.sha256(output.encode()).hexdigest()}))
