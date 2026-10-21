"""All-future certificates from stopped finite prefixes; never reads a suffix.

Source carrier (N,D,s) is lossless: decode() recovers its literal prefix,
registers, stopped node, depth and complete integer clock. Aggregation is
only for the energy bound. No unresolved source is assigned a terminal.
"""
from collections import Counter,defaultdict
from hashlib import sha256
from importlib.util import spec_from_file_location,module_from_spec
from itertools import product
from pathlib import Path
import argparse,json,sys

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('future_native',HERE.parent/'session18_native_annulus/verify_annulus.py')
native=module_from_spec(spec);sys.modules[spec.name]=native;spec.loader.exec_module(native)
CHECKS=0

def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def integer(n,lo=0):
    if type(n) is not int or n<lo:raise ValueError('literal integer domain')
    return n

def vector(v):
    if type(v) is not tuple or len(v)!=4 or any(type(x) is not int for x in v):
        raise ValueError('four literal integer charge coordinates')
    return v

def Q(v):return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def beta(v,w):return 2*v[0]*w[0]-v[0]*w[1]-v[1]*w[0]+2*v[1]*w[1]
def delta(v):return Q(v[2:])+beta(v[:2],v[2:])
def cross(v,w):return beta(v[2:],w[2:])+beta(v[:2],w[2:])+beta(w[:2],v[2:])
def add(v,w):return tuple(a+b for a,b in zip(v,w))
def scaled(v,sg):return tuple(sg*x for x in v)
def enc(x):return (json.dumps(x,separators=(',',':'),sort_keys=True)+'\n').encode()

def envelope(arrived,pending):
    """Exact formula; repeated vector types only accelerate scalar summation."""
    for v in arrived+pending:vector(v)
    current=sum(delta(v) for v in arrived)
    types=list(Counter(max(v,scaled(v,-1)) for v in pending).items())
    ports=set(arrived)
    local=sum(n*max(0,delta(v)+max([0]+[abs(cross(w,v)) for w in ports])) for v,n in types)
    pairs=sum(n*(n-1)//2*abs(cross(v,v)) for v,n in types)
    pairs+=sum(n*m*abs(cross(v,w)) for i,(v,n) in enumerate(types) for w,m in types[:i])
    return current,local,pairs,current+local+pairs

def physical(ports):
    return sum(Q((v[0]+v[2],v[1]+v[3]))-Q(v[:2]) for v in ports)

def frontier():
    raw=(HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n',b'\n')
    need(sha256(raw).hexdigest()=='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f','fixed grounded frontier')
    obj=json.loads(raw);r=native.Routes(obj['frontier'],dict(obj['frontier_ranks']))
    return r.V

def decode(N,D,s,V):
    integer(N,1);integer(D)
    if not native.source(s) or s>=256*N:raise ValueError('literal source in declared carrier')
    u=s;t=0;word=[]
    for _ in range(D):
        if u in V:break
        # Independent literal implementation of the primitive edge.
        a=2 if u%8==1 else 1;v=(3*u+1)//(1<<a)
        while v%8==5:v=(v-1)//4;a+=2
        word.append(a);t+=a;u=v
    return {'source':s,'mod15':s%15,'mod35':s%35,'node':u,
            'P_depth':len(word),'integer_clock':t,'edge_clocks':word,
            'arrived':u in V,'old':s<4*N,'charge':list(native.charge(64*N,s))}

def prefix(N,D,V,verbose=False):
    integer(N,1);integer(D)
    arrived={};pending={};source_hash=sha256();edge_hash=sha256();drop_hash=sha256()
    initial_zero=initial_nonzero=annulus_norm=0
    # Mutable group rows: four coordinates,count,minclock,maxclock,source witness,its clock.
    for s in range(1,256*N,2):
        if not native.source(s):continue
        q=native.charge(64*N,s);old=s<4*N
        if old:need(q==native.charge(N,s),'old named charge stable under scale64')
        else:annulus_norm+=Q(q)
        need(s%15%5==s%35%5,'compatible source registers')
        source_hash.update(enc([s,s%15,s%35,int(old),*q]))
        if q==(0,0):initial_zero+=1;continue
        initial_nonzero+=1;v=q+(0,0) if old else (0,0)+q
        if s in V:arrived[s]=v
        else:pending[s]=[*v,1,0,0,s,0]
    need(annulus_norm==153*N,'unchanged exact source budget')
    history=[];cache={};drops=0;dropped_sources=0
    for d in range(D+1):
        current=physical(list(arrived.values()))
        history.append([d,current,len(pending),sum(x[4] for x in pending.values())])
        if d==D:break
        nxt={}
        for u,v in pending.items():
            if u not in cache:cache[u]=native.edge(u)
            z,c=cache[u];sg=1-2*(c%2);w=scaled(tuple(v[:4]),sg)
            lo,hi,ct=v[5]+c,v[6]+c,v[8]+c
            edge_hash.update(enc([d,u,z,c,*v]))
            if z in V:
                arrived[z]=add(arrived.get(z,(0,0,0,0)),w)
            elif z in nxt:
                r=nxt[z]
                for k in range(4):r[k]+=w[k]
                r[4]+=v[4];r[5]=min(r[5],lo);r[6]=max(r[6],hi)
                if v[7]<r[7]:r[7]=v[7];r[8]=ct
            else:nxt[z]=[*w,v[4],lo,hi,v[7],ct]
        pending={}
        for u,v in nxt.items():
            if any(v[:4]):pending[u]=v
            else:
                drops+=1;dropped_sources+=v[4];drop_hash.update(enc([d+1,u,*v]))
        if verbose and (d+1)%20==0:print(json.dumps({'N':N,'prefix_depth':d+1,'pending':len(pending)}),flush=True)
    E=envelope(list(arrived.values()),[tuple(v[:4]) for v in pending.values()])
    need(E[0]==history[-1][1],'quadratic identity agrees with physical difference')
    witnesses=[]
    for u,v in sorted(pending.items()):
        w=decode(N,D,v[7],V)
        need(not w['arrived'] and w['P_depth']==D and w['node']==u,'literal unresolved witness at cutoff')
        need(w['integer_clock']==v[8] and v[5]<=v[8]<=v[6],'retained full clock witness and range')
        witnesses.append(w)
    rows=[[u,*v] for u,v in sorted(pending.items())]
    result={'N':N,'D':D,'budget':153*N,'source_carrier':['primitive s',1,256*N,'mod15','mod35','literal prefix at depth D'],
            'source_rows_sha256':source_hash.hexdigest(),'prefix_edges_sha256':edge_hash.hexdigest(),
            'zero_drops_sha256':drop_hash.hexdigest(),'zero_groups':drops,'zero_group_sources':dropped_sources,
            'initial_zero_sources':initial_zero,'initial_nonzero_sources':initial_nonzero,
            'history':history,'arrived_ports':[[u,*v] for u,v in sorted(arrived.items())],
            'pending_row_schema':['node','W0','W1','H0','H1','source_count','min_clock','max_clock','witness','witness_clock'],
            'pending':rows,'witnesses':witnesses,'envelope':list(E),
            'early_cuts_pass':all(row[1]<=153*N for row in history),'future_tail_pass':E[-1]<=153*N}
    return result

def rejects(fn):
    try:fn()
    except (ValueError,TypeError):return True
    return False

def controls():
    ports=[(1,0,0,1),(0,1,-1,0)]
    vectors=list(product(range(-1,2),repeat=4))
    choices=[None,(0,1),(0,-1),(1,1),(1,-1)]
    for a in vectors:
        for b in vectors:
            e=envelope(ports,[a,b])[-1]
            direct=sum(max(0,delta(v)+max(abs(cross(v,p)) for p in ports)) for v in (a,b))+abs(cross(a,b))+physical(ports)
            need(e==direct,'type aggregation preserves original pair sum')
            for ca,cb in product(choices,repeat=2):
                pp=ports.copy()
                for v,c in ((a,ca),(b,cb)):
                    if c is not None:pp[c[0]]=add(pp[c[0]],scaled(v,c[1]))
                need(physical(pp)<=e,'all optional signed terminal assignments')
            need(envelope(ports,[add(a,b)])[-1]<=e,'coalescence monotonicity')
            need(envelope([add(ports[0],a),ports[1]],[b])[-1]<=e,'arrival monotonicity')
    # Opposite pending values at different states may arrive at different cuts.
    p=(0,0,1,0);m=(0,0,-1,0)
    need(physical([p])==1 and envelope([(0,0,0,0)],[p,m])[-1]>=1,'different-node zero sum cannot be erased')
    mixed=(1,0,-1,0)
    need(delta(mixed)==-1 and cross(mixed,(1,0,0,0))==-2,'mixed total zero retains old baseline and interactions')
    need(envelope(ports,[(0,0,0,0)])==envelope(ports,[]),'both-channel zero block safely erased')
    need(rejects(lambda:envelope(ports,[(False,0,0,0)])),'Boolean vector alias rejected')
    need(rejects(lambda:integer(1.0)),'floating depth rejected')

def run(cases,verbose=False):
    controls();V=frontier();out=[]
    for N,D in cases:
        r=prefix(N,D,V,verbose);out.append(r)
        if N>=4096:
            need(r['early_cuts_pass'] and r['future_tail_pass'],'all-horizon certificate consequence')
            need(len(r['pending'])>0,'nonvacuous unresolved continuation')
        elif N==7:
            need(not r['early_cuts_pass'] and not r['future_tail_pass'],'N7 actual failure cannot be erased')
            need(next(x[0] for x in r['history'] if x[1]>153*N)==8,'canonical first hostile cut')
    return {'status':'PROVED conditional all-future envelope and monotonicity; FINITE-EXACT declared scales; uniformN OPEN',
            'cases':out,'checks':CHECKS}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--verbose',action='store_true');a=p.parse_args()
    cases=[(7,12),(4096,105),(8192,110),(32768,130)]
    result=run(cases,a.verbose);data=enc(result);path=HERE/'receipt.json'
    if a.write:path.write_bytes(data)
    else:need(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved immutable receipt comparison')
    print(json.dumps({'checks':CHECKS,'sha256_lf':sha256(data).hexdigest(),
                      'cases':[[r['N'],r['D'],len(r['pending']),r['envelope'],r['budget']] for r in result['cases']]}))
