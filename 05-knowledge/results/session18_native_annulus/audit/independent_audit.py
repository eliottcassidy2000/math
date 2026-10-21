#!/usr/bin/env python3
"""Independent, non-importing audit of the frozen native-annulus packet.
Charges are rebuilt by assigning every original odd integer to its primitive
ancestor, not by the producer's per-source B-chain summation or scale sweep.
Only the explicitly listed scales below are independently recomputed.
"""
import argparse, hashlib, json
from array import array
from pathlib import Path
from itertools import product

CHECKS=0
H=257
L=((1,2),(-2,-1),(1,-1)); U=((0,1),(-1,-1),(1,0))
BASES=(1,2,3,4,5,6,7,8,448,512,4096,5997,11002,16384,28672,32768)
MANIFEST='8423a629603147141ef1a5f5570c7e74aecc6f55505f201c0d598d9aac6ad00a'

def require(test, message):
    global CHECKS
    CHECKS+=1
    if not test: raise RuntimeError(message)

def enc(x): return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def sh(b): return hashlib.sha256(b).hexdigest()
def rowenc(r): return ('['+','.join(map(str,r))+']\n').encode()
def q(a,b): return a*a-a*b+b*b
def bet(a,b,c,d): return 2*a*c-a*d-b*c+2*b*d

def primitive(n): return type(n) is int and n>0 and n%8 in (1,3,7)
def ancestor(n):
    while n%8==5: n=(n-1)//4
    return n

def edge(n):
    if not primitive(n): raise RuntimeError('invalid native edge source')
    clock=2 if n%8==1 else 1
    z=(3*n+1)//(1<<clock)
    while z%8==5:
        z=(z-1)//4; clock+=2
    if not primitive(z): raise RuntimeError('invalid native edge target')
    return z,clock

class Paths:
    def __init__(self,V,ranks,limit):
        self.V=V; self.index={s:i for i,s in enumerate(V)}; self.limit=limit
        self.terminal=array('h',[-1])*(limit//2)
        self.clock=array('I',[0])*(limit//2)
        self.depth=array('H',[0])*(limit//2)
        self.outer={};self.requested_depth=0
        require(len(V)==627 and len(set(V))==627,'627 distinct frontier seeds')
        require(set(V)==set(ranks) and ranks[1]==0,'frontier rank domain')
        for i,s in enumerate(V):
            require(primitive(s) and type(ranks[s]) is int and ranks[s]>=0,'frontier primitive rank')
            if s!=1:
                z,c=edge(s);require(z in ranks and ranks[z]<ranks[s],'frontier descending rank')
            self.put(s,(i,0,0))
    def lookup(self,s):
        if s<self.limit:
            i=s//2;t=self.terminal[i]
            return None if t<0 else (t,self.clock[i],self.depth[i])
        return self.outer.get(s)
    def put(self,s,r):
        if s<self.limit:
            i=s//2;self.terminal[i],self.clock[i],self.depth[i]=r
        else:self.outer[s]=r
    def get(self,s):
        trail=[];seen=set();z=s;r=self.lookup(z)
        while r is None:
            require(z not in seen and len(trail)<10000,'route repeats or exceeds cap')
            seen.add(z);n,c=edge(z);trail.append((z,c));z=n;r=self.lookup(z)
        t,c,d=r
        for z,dc in reversed(trail):
            c+=dc;d+=1;self.put(z,(t,c,d))
        r=self.lookup(s);require(0<=r[2]<H,'requested arrival cap')
        self.requested_depth=max(self.requested_depth,r[2]);return r
    def audit(self):
        dig=hashlib.sha256();count=0
        def verify(s,r):
            nonlocal count
            u,c,d=r;end=self.V[u]
            if s in self.index:
                require((end,c,d)==(s,0,0),'frontier first-hit stopping');row=[s,s,0,end,c,d]
            else:
                z,dc=edge(s);child=self.lookup(z)
                require(child is not None and (u,c,d)==(child[0],child[1]+dc,child[2]+1),'independent full-clock suffix')
                row=[s,z,dc,end,c,d]
            dig.update(rowenc(row));count+=1
        for i,t in enumerate(self.terminal):
            if t>=0:verify(2*i+1,(t,self.clock[i],self.depth[i]))
        for s,r in sorted(self.outer.items()):verify(s,r)
        return {'nodes':count,'sha256':dig.hexdigest(),'max_requested_depth':self.requested_depth}

def direct(M,paths,keep_rows=False):
    # Scan each original odd m once; put its shell charge at rho(m).
    ca=array('i',[0])*(2*M);cb=array('i',[0])*(2*M)
    for m in range(1,4*M,2):
        i=ancestor(m)//2;a,b=(L if m<2*M else U)[m%3];ca[i]+=a;cb[i]+=b
    nV=len(paths.V);size=H*nV
    wa=array('q',[0])*size;wb=array('q',[0])*size;mass=[0]*H
    digest=hashlib.sha256();count=0;maximum=0;kept=[]
    for s in range(1,4*M,2):
        if not primitive(s):continue
        u,t,d=paths.get(s);a,b=ca[s//2],cb[s//2];norm=q(a,b)
        r15,r35=s%15,s%35
        require(r15%5==r35%5 and next(x for x in range(r15,105,15) if x%35==r35)==s%105,'joint 15/35 register')
        row=[s,r15,r35,paths.V[u],t,d,a,b,norm];digest.update(rowenc(row));count+=1;maximum=max(maximum,d)
        if keep_rows:kept.append(row)
        k=d*nV+u;sign=-1 if t%2 else 1;wa[k]+=sign*a;wb[k]+=sign*b;mass[d]+=norm
    energy=[]
    for d in range(H):
        if d:
            mass[d]+=mass[d-1]
            for u in range(nV):
                k=d*nV+u;wa[k]+=wa[k-nV];wb[k]+=wb[k-nV]
        energy.append(sum(q(wa[d*nV+u],wb[d*nV+u]) for u in range(nV)))
    info={'N':M,'sources':count,'source_norm':mass[-1],'final_energy':energy[-1],
          'max_depth':maximum,'source_rows_sha256':digest.hexdigest()}
    return {'E':energy,'M':mass,'a':wa,'b':wb,'info':info,'rows':kept}

def compare(N,old,big,nV):
    rows=[]
    for d in range(H):
        h=mix=0
        for u in range(nV):
            k=d*nV+u;oa,ob=old['a'][k],old['b'][k];ha,hb=big['a'][k]-oa,big['b'][k]-ob
            h+=q(ha,hb);mix+=bet(oa,ob,ha,hb)
        arrived=big['M'][d]-old['M'][d];delta=big['E'][d]-old['E'][d]
        require(delta==h+mix,'source-owned arrived energy decomposition')
        rows.append([d,old['E'][d],h,mix,arrived,153*N-arrived,h-arrived,h+mix-arrived,153*N-delta])
    require(big['M'][-1]-old['M'][-1]==153*N,'annulus exact source norm')
    return rows

def explicit_pairs(N,rows):
    # Independent O(number of same-terminal pairs) event ownership at N=7.
    vs=[]
    for s,r15,r35,u,t,d,a,b,norm in rows:
        sg=(-1)**t;vs.append((s<4*N,u,d,sg*a,sg*b,norm))
    debit=[0]*H;diag=[0]*H;cross=[0]*H
    for i,(old,u,d,a,b,norm) in enumerate(vs):
        if not old:diag[d]+=norm
        for old2,v,e,c,f,norm2 in vs[:i]:
            if u==v and not (old and old2):cross[max(d,e)]+=bet(a,b,c,f)
    for d in range(H):
        debit[d]=diag[d]+cross[d]+(debit[d-1] if d else 0)
    return debit

def main():
    p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--original-root',type=Path);p.add_argument('--out',type=Path,required=True);args=p.parse_args();root=args.package
    require(sh((root/'MANIFEST.json').read_bytes())==MANIFEST,'frozen producer manifest hash')
    manifest=json.loads((root/'MANIFEST.json').read_text());provenance=json.loads((root/'source_manifest.json').read_text());original_checks=[]
    for f in manifest['files']:
        b=(root/f['path']).read_bytes();require(len(b)==f['bytes'] and sh(b)==f['sha256'],'packet member '+f['path'])
    for normal,opt,receipt in [('normal.log','optimized.log','annulus_receipt.json'),('rays-normal.log','rays-optimized.log','ray_receipt.json'),('recurrence-normal.log','recurrence-optimized.log','recurrence_receipt.json')]:
        require((root/normal).read_bytes()==(root/opt).read_bytes(),'producer normal/optimized equality')
        printed=json.loads((root/normal).read_text());require(printed['status']=='PASS' and printed['receipt_sha256']==sh((root/receipt).read_bytes()),'parse printed receipt witness')
    fd=json.loads((root/'frontier.json').read_text());original_rows={}
    if args.original_root:
        for f in provenance['source_files']:
            b=(args.original_root/f['path']).read_bytes();require(len(b)==f['bytes'] and sh(b)==f['sha256'],'original source authentication '+f['path']);original_checks.append(f['path'])
        source=args.original_root/'05-knowledge/results/native_cancellation_105_20261018g.json';b=source.read_bytes();j=json.loads(b)
        require(sh(b)==fd['input_sha256'] and fd['frontier']==j['frontier'] and fd['frontier_ranks']==j['frontier_ranks'],'original frontier provenance')
        for item in j['scales']:
            rows=item['source_rows'];n=item['statistics']['N'];digest=sh(b''.join(rowenc(r) for r in rows));original_rows[n]=rows
            require(digest==provenance['inherited_canonical_jsonl_source_rows'][str(n)]['sha256'],'original canonical row digest')
    paths=Paths(fd['frontier'],dict(fd['frontier_ranks']),4*64*max(BASES))
    ann=json.loads((root/'annulus_receipt.json').read_text());ray=json.loads((root/'ray_receipt.json').read_text());expected={}
    for item in ann['direct_checks'].values():
        for name in ('old','large'):expected[item[name]['N']]=item[name]
    for item in ray['results']:
        for name in ('old','large'):expected[item[name]['N']]=item[name]
    needed=sorted(set(BASES)|{64*n for n in BASES});ledgers={};infos=[]
    for M in needed:
        data=direct(M,paths,keep_rows=M in set(original_rows)|{448,262144});ledgers[M]=data;infos.append(data['info'])
        if M in expected:require(data['info']==expected[M],'independent source rows at '+str(M))
        if M in original_rows:require(data['rows']==original_rows[M],'all inherited source rows at '+str(M))
        print(json.dumps({'rebuilt':M,'sources':data['info']['sources'],'max_depth':data['info']['max_depth']},sort_keys=True),flush=True)
    details={};witnesses=[]
    for N in BASES:
        old,big=ledgers[N],ledgers[64*N];rows=compare(N,old,big,len(paths.V));details[N]=rows
        if str(N) in ann['details']:require(rows==ann['details'][str(N)]['horizons'],'all selected annulus horizon rows')
        if N<=7:
            bad=next((r for r in rows if r[-1]<0),None);require(bad==ann['small_hostiles'][N-1]['first_failed_horizon'],'small hostile horizon')
        if N in ray['scales']:
            saved=next(x for x in ray['results'] if x['N']==N)
            rr=[[d,old['E'][d],big['E'][d],old['M'][d],big['M'][d],big['E'][d]-old['E'][d],rows[d][7],rows[d][8]] for d in range(H)]
            require(rr==saved['horizons'],'all independent ray horizons')
        peak=max(big['E'][d]-old['E'][d] for d in range(H));at=next(d for d in range(H) if big['E'][d]-old['E'][d]==peak)
        witnesses.append({'N':N,'peak':peak,'depth':at,'budget':153*N,'min_bank':min(r[-1] for r in rows),'max_cross':max(r[7] for r in rows)})
    require(explicit_pairs(7,ledgers[448]['rows'])==[r[2]+r[3] for r in details[7]],'all explicit diagonal/pair event ownership N7')
    # Literal native partial-packet route and boundary verification.
    pair=ann['details']['4096']['new_only_partial_pair'];literal=[]
    for s,t,d,a,b in pair['sources']:
        row=next(r for r in ledgers[262144]['rows'] if r[0]==s);require(row[3:]==[pair['terminal'],t,d,a,b,q(a,b)],'literal partial packet source');literal.append(row)
    require(pair['sources']==[[276145,24,4,-2,-1],[92051,24,5,2,1]],'README native partial pair')
    # Actual positive fibre count and terminal-nine source IDs at N4096,d1.
    old,big=ledgers[4096],ledgers[262144];nV=len(paths.V);byterminal={};ids={}
    for r in big['rows']:
        s,_,_,u,t,d,a,b,norm=r
        if d<=1:
            ids.setdefault((u,s<16384),[]).append(s)
            if s>=16384:byterminal[u]=byterminal.get(u,0)+norm
    positive=[]
    for i,u in enumerate(paths.V):
        k=nV+i;oa,ob=old['a'][k],old['b'][k];ha,hb=big['a'][k]-oa,big['b'][k]-ob
        if q(ha,hb)+bet(oa,ob,ha,hb)>byterminal.get(u,0):positive.append(u)
    require(len(positive)==204 and ids[(9,True)]==[9,6371] and ids[(9,False)]==[203889,407779],'literal 204 positive fibres and terminal-nine source IDs')
    require(details[4096][1][2:8]==[2064,-1129,1860,624828,204,-925],'N4096 mixed debt')
    require(details[448][1][7]==47 and details[448][1][5]==66707 and details[512][1][7]==60 and details[512][1][5]==76477,'moving collateral witnesses')
    # Algebraic partial/mixed-packet boundary identity, including cross-packet terms.
    algebra_cases=0
    for vals in product((-1,0,1),repeat=8):
        ra,rb,ta,tb,oa,ob,ja,jb=vals;ba,bb=oa+ja,ob+jb
        actual=q(ra+ta+ba,rb+tb+bb)-q(ra+oa,rb+ob)
        corrected=q(ra+ta,rb+tb)-q(ra,rb)+q(ba,bb)-q(oa,ob)+bet(ra+ta,rb+tb,ba,bb)-bet(ra,rb,oa,ob)
        require(actual==corrected,'packet boundary identity D');algebra_cases+=1
    require(q(1,0)-q(2,0)==-3 and -q(1,0)-bet(1,0,1,0)==-3,'mixed zero packet old-baseline hostile')
    dag=paths.audit();largest=next(x for x in ray['results'] if x['N']==32768)['route_dag'];require(dag==largest,'complete independently rebuilt largest native route DAG')
    # Recheck frozen producer bytes after work.
    for f in manifest['files']:require(sh((root/f['path']).read_bytes())==f['sha256'],'producer preservation '+f['path'])
    result={'status':'PASS; selected complete universes independently replayed; no universal inference','producer_manifest_sha256':MANIFEST,'base_scales':BASES,'source_scales':needed,'horizons':[0,256],
      'source_rows_rebuilt':sum(i['sources'] for i in infos),'original_files_authenticated':original_checks,'original_rows_compared':sum(map(len,original_rows.values())),
      'independent_method':'Original odd-integer ancestor accumulation; own array-backed native first-hit DAG; full clock and joint-register rows; explicit N7 unordered pair ownership',
      'source_universes':infos,'witnesses':witnesses,'N7_first_failure':details[7][8],'N4096_d1':details[4096][1],
      'N448_d1':details[448][1],'N512_d1':details[512][1],'partial_pair':pair,'positive_fibres_N4096_d1':positive,
      'packet_identity_cases':algebra_cases,'largest_route_dag':dag,'checks':CHECKS,
      'not_replayed':'Contiguous production sweep 4096..16384 and its global maximum were inspected and authenticated, not exhaustively independently repeated by this audit. No unbounded-N or full-state moving-Gram assertion.'}
    args.out.write_bytes(enc(result));print(json.dumps({'status':'PASS','checks':CHECKS,'audit_sha256':sh(enc(result)),'source_rows_rebuilt':result['source_rows_rebuilt'],'route_dag':dag},sort_keys=True),flush=True)
if __name__=='__main__':main()
