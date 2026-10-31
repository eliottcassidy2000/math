#!/usr/bin/env python3
"""Owned Q25 head-cut carry certificates; no external Gerbicz C dependency."""
import argparse,gzip,hashlib,json,math,struct,sys
from functools import lru_cache
from pathlib import Path
sys.dont_write_bytecode=True
CATALOG_PIN='51af1a4ec0e4e8758b2ddbcb9220859c1b568c4031ee6d2a0aa642f6c093d2ab'
MASKS=[1585793692427,1599483269979,1602705117979,3372531545867,3389445068555,3526867243787]
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def E(p):return {tuple(sorted(e)) for e in zip(p,p[1:])}
def delta(p,q):return len(E(p)-E(q))
def checkpair(n,pair):
 phases=[]
 for ex,p in enumerate(pair):
  need(len(p)==n+ex and sorted(p)==list(range(1,n+ex+1)),'coverage')
  need(list(p[:2])==[1,3] and p[-1]==(10 if len(p)%2 else 8),'weak B ports')
  need(all(math.isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])),'square edge')
  need(all(not(set(edge)<={1,3,8,10}) for edge in E(p)-{(1,3)}),'non-port edge exclusion')
  phases.append({v:(j+1)%2 for j,v in enumerate(p)})
 need(all(phases[0][v]==phases[1][v] for v in range(1,n+1)),'common phase')
def boundary(word,parent_parity):
 ends=[];edges=set()
 for atom in word:
  kind,c=atom[:2]
  if kind=='single':ends.append((c,c));continue
  if kind=='head':ends.append((25-c,25-c));continue
  e,o=atom[2:];p=(parent_parity+e)%2;b=10 if p else 8
  a=25-c if kind=='block' else 75+c;z=25*b+(-1)**p*c
  ends.append((a,z) if o==1 else (z,a))
  if kind=='block':edges.add(tuple(sorted((25-c,75+c))))
 edges.update(tuple(sorted((a[1],b[0]))) for a,b in zip(ends,ends[1:]))
 need(len(edges)==61,'boundary must have61 edges')
 need(all({(u+12)//25,(v+12)//25}<={0,1,3,8,10} for u,v in edges),'boundary decoded ports')
 return edges
class Owned:
 def __init__(self,root,base_dir,small):
  self.root=Path(root);self.base_dir=Path(base_dir)
  deps=json.loads((self.root/'DEPENDENCIES.json').read_text());raw=Path(small).read_bytes();need(sha(raw)==deps['small_base_sha256'],'small base pin');self.small=json.loads(raw)
  self.base_deps=json.loads((self.root/'BASE_DEPENDENCIES.json').read_text());r=(self.base_dir/self.base_deps['manifest_file']).read_bytes();need(sha(r)==self.base_deps['manifest_sha256'],'base manifest pin');need(json.loads(r)==self.base_deps['chunks'],'base manifest content')
  self.cert_pins=json.loads((self.root/'CERTIFICATES.json').read_text());self.used=set()
 @lru_cache(None)
 def cert(self,s,p,d):
  name=f'{s}/paired/B_p{p}_d{d}.json';raw=(self.root/name).read_bytes();need(sha(raw)==self.cert_pins[name],'certificate pin');self.used.add(name);return json.loads(raw)
 @lru_cache(None)
 def chunk(self,name):
  m=next(c for c in self.base_deps['chunks'] if c['file']==name);z=(self.base_dir/name).read_bytes();need(sha(z)==m['sha256'] and len(z)==m['bytes'],'compressed base pin');raw=gzip.decompress(z);need(sha(raw)==m['decoded_sha256'] and len(raw)==m['decoded_bytes'],'decoded base pin');need(raw[:8]==b'SSBAS1\0\0' and struct.unpack('<II',raw[8:16])==(m['first'],m['last']),'base header');return raw
 @lru_cache(None)
 def base(self,n):
  need(42<=n<=1061,'base interval')
  if n<53:p=self.small[str(n)]
  else:
   m=next(c for c in self.base_deps['chunks'] if c['first']<=n<=c['last']);raw=self.chunk(m['file']);lo=m['first'];i=16+(n-lo)*(n+lo-2);p=[]
   for L in (n,n+1):
    w=[1]
    for r in raw[i:i+L-1]:w.append(r*r-w[-1])
    i+=L-1;p.append(w)
  checkpair(n,p);return tuple(map(tuple,p))
 @lru_cache(maxsize=32)
 def pair(self,s,n):
  if n<=1061:return self.base(n)
  m=(n-12)//25;d=n-25*m;source=self.pair(s,m);out=[]
  for word in self.cert(s,m%2,d)['words']:
   w=[]
   for atom in word:
    kind,c=atom[:2]
    if kind=='single':w.append(c)
    elif kind=='head':w.append(25-c)
    else:
     e,o=atom[2:];ix=range(int(kind=='tail'),len(source[e]));ix=ix if o==1 else reversed(ix)
     w.extend(25*source[e][j]+(-1)**(j+1)*c for j in ix)
   out.append(tuple(w))
  checkpair(n,out);return tuple(out)
def phase_obstruction(catalog):
 raw=Path(catalog).read_bytes();need(sha(raw)==CATALOG_PIN,'weak B catalogue pin');data=json.loads(raw);onset=data['pairs'][-1]
 # The finite onset catalogue is separately independently proved; no larger census.
 for w in onset['witnesses']:
  checkpair(42,[w['shorter'],w['longer']])
  need(sum(1<<(v-1) for v in w['shorter'][::2])==w['common_odd_mask'],'positive phase witness')
 masks=sorted(w['common_odd_mask'] for w in onset['witnesses']);need(onset['n']==42 and masks==MASKS,'six-mask onset catalogue')
 out=[]
 for m in MASKS:
  n=46 if m in [MASKS[0],MASKS[1],MASKS[3]] else 47
  ph=lambda v:(m>>(v-1))&1 if v<=42 else v%2
  ns=[u for u in range(1,n) if math.isqrt(u+n)**2==u+n and ph(u)!=ph(n)]
  expected=[35] if m==MASKS[1] else ([18] if n==46 else [34])
  need(ns==expected,'new-label degree obstruction');need(n not in [1,8 if n%2==0 else 10],'new label internal')
  out.append(dict(mask=m,failing_order=n,internal_vertex=n,compatible_neighbors=ns))
 return out
def run(a):
 o=Owned(a.operator_dir,a.base_dir,a.small_bases)
 base_stats=[(n,delta(*o.base(n))) for n in range(42,1062)]
 gap={n:delta(o.base(n-1)[1],o.base(n)[0]) for n in range(43,47)}
 families={}
 for s in ['A','B']:
  state=[];noncarry=[]
  for p in [0,1]:
   for d in range(12,37):
    cert=o.cert(s,p,d);w0,w1=cert['words'];members=[]
    for word in [w0,w1]:
     members.append({x[1]:x[2] for x in word if x[0] in ['block','tail']})
    changed=[c for c in members[0] if members[0][c]!=members[1][c]]
    need(changed==[(-1)**(p+1)*(d+1-25)],'single changed input member')
    B0=boundary(w0,p);B1=boundary(w1,p);k=len(B0-B1)
    m=42 if p==0 else 43;n=25*m+d;pair=o.pair(s,n)
    need(delta(*pair)==delta(*o.pair(s,m))+k,'exact within-pair recurrence')
    state.append(dict(parity=p,digit=d,boundary_deletions=k))
    if d<36:
     nxt=o.cert(s,p,d+1)['words'][0];noncarry.append(len(B1-boundary(nxt,p)))
  carry={}
  for p in [0,1]:
   old=o.cert(s,1-p,36)['words'][1];new=o.cert(s,p,12)['words'][0]
   carry[p]=len(boundary(old,1-p)-boundary(new,p))
  rays=[]
  for seed in [44,45]:
   r=seed;prev=None
   for depth in range(3):
    left=o.pair(s,r-1);right=o.pair(s,r);h=delta(left[1],right[0]);canonical=delta(left[0],right[0]);within=delta(*left)
    if prev is not None:need(h==25*prev+carry[seed%2],'exact full carry recurrence')
    need(h-within-1<=canonical<=h+within,'canonical triangle')
    if depth:
     hm=delta(o.pair(s,r-2)[1],left[0]);hp=delta(right[1],o.pair(s,r+1)[0]);need(max(hm,hp)<=61,'noncarry bank neighbors')
     for X in [o.pair(s,r-2)[1],left[0]]:
      need(delta(X,right[0])>=h-within-62,'short middle incoming')
     for Y in [right[1],o.pair(s,r+1)[0]]:
      need(delta(left[1],Y)>=h-delta(*right)-62,'long middle outgoing')
    rays.append(dict(seed=seed,depth=depth,order=r,gap=h,canonical_deletions=canonical,within_pair_deletions=within));prev=h;r=25*r+12
  families[s]=dict(within_pair_states=state,within_pair_K_range=[min(x['boundary_deletions'] for x in state),max(x['boundary_deletions'] for x in state)],noncarry_gap_max=max(noncarry),carry_K=carry,ray_samples=rays)
 negative={}
 for name,bad in [('duplicate',[[1,1],[1,3,10]]),('wrong_prefix',[[3,1],[1,3,10]])]:
  try:checkpair(2,bad);negative[name]=False
  except RuntimeError:negative[name]=True
 need(all(negative.values()) and math.isqrt(11)**2!=11,'hostile controls / odd nonclosure')
 return dict(status='PASS',weak_B_catalogue_sha256=CATALOG_PIN,operator_certificate_manifest_sha256=sha((Path(a.operator_dir)/'CERTIFICATES.json').read_bytes()),owned_base_manifest_sha256=o.base_deps['manifest_sha256'],small_base_sha256=sha(Path(a.small_bases).read_bytes()),hostile_controls=negative,odd_weak_B_path_cannot_close_to_1=True,scope='weak B path banks, not direct cycle closures',external_C_required=False,base_range=[42,1061],base_pair_count=len(base_stats),base_within_pair_max=max(base_stats,key=lambda x:x[1]),base_junction_gaps=gap,phase_obstructions=phase_obstruction(a.weak_b_catalogue),universal_bank_path_limsup_lower_bound='4/93',families=families,used_paired_certificates=len(o.used))
RESULTS=Path(__file__).resolve().parents[2]
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--operator-dir',default=str(RESULTS/'session124_mixed_radix_square_rays_20261031'));p.add_argument('--base-dir',default=str(RESULTS/'square_sum_B_base_interval_20261019'));p.add_argument('--small-bases',default=str(RESULTS/'session124_square_prefix_threshold42_20261031'/'SMALL_BASES.json'));p.add_argument('--weak-b-catalogue',default=str(RESULTS/'session125_origin32_numeration_20261031'/'pair_taxonomy'/'weak_b_certificate.json'));p.add_argument('--output',required=True);a=p.parse_args();out=run(a);Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('PASS',a.output)
