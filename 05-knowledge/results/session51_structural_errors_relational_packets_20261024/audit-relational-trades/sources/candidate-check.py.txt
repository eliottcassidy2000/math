#!/usr/bin/env python3
"""Exact shared-record-port restriction on the authentic nineteen-slot receiver."""
from pathlib import Path
from collections import defaultdict, Counter, deque
from hashlib import sha256
from math import comb
import argparse, importlib.util, json
P=Path(__file__).resolve().parent
CHECKS=0
Q=set(range(2048)); FULL=(1<<17)-1; TBIT=1<<16

def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(why)
def canon(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def add(x,y):return (((x//16+y//16)%128)<<4)|((x^y)&15)
def sub(x,y):return (((x//16-y//16)%128)<<4)|((x^y)&15)
def total(xs):
 q=0
 for x in xs:q=add(q,x)
 return q
def cost(d):
 t=d>>16;b=(d&65535).bit_count();return (b+3*t,2*b+8*t)
def joint(c):return [2*(c//4),c]
def load():
 for name,row in json.loads((P/'SOURCEPINS.json').read_text()).items():
  data=(P/name).read_bytes();need(len(data)==row['bytes'] and sha256(data).hexdigest()==row['sha256'],'source pin '+name)
 spec=importlib.util.spec_from_file_location('inherited_receiver',P/'inherited/receiver-thresholds/check.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 old.source_pins();n=old.native();return old,n

def enumerate_states(n):
 groups=[[s] for s in n['base']]+[[1457,1841,3015]]
 labels=[total(n['labels'][s] for s in g) for g in groups];grades=[len(g) for g in groups]
 gd=[0]*(FULL+1);ch=[0]*(FULL+1);phy=[None]*(FULL+1);phy[0]=tuple(0 for m in n['moduli']);fib=defaultdict(list);rows={}
 delta=[tuple(sum(2*(f[s]-f[n['U']]) for s in g)%m for f,m in zip(n['potentials'],n['moduli'])) for g in groups]
 for x in range(FULL+1):
  if x:
   b=x&-x;j=b.bit_length()-1;y=x^b;gd[x]=gd[y]+grades[j];ch[x]=add(ch[y],labels[j]);phy[x]=tuple((a+b)%m for a,b,m in zip(phy[y],delta[j],n['moduli']))
  key=gd[x],ch[x];fib[key].append(x)
  if key not in rows:rows[key]=phy[x]
  need(rows[key]==phy[x],'original physical source and target rows constant on full original grade/charge fibre')
 return dict(sorted(fib.items())),gd,ch,groups

def restricted(fib,j):return {k:w for k,v in fib.items() if (w:=[x for x in v if ((x>>j)&1)==(x>>16)])}

def solve(v):
 """Floyd ordinary minimax and interval DP, independently of inherited MST solver."""
 size=len(v);d=[[cost(x^y)[1] for y in v] for x in v];ordinary=[r[:] for r in d];strict=[r[:] for r in d];prim=[]
 for k in range(size):
  for i in range(size):
   for j in range(size):ordinary[i][j]=min(ordinary[i][j],max(ordinary[i][k],ordinary[k][j]))
 pairs=sorted((cost(v[i]^v[j])[0],i,j) for i in range(size) for j in range(i+1,size))
 for _,i,j in pairs:
  middle=[k for k in range(size) if k not in (i,j) and ((v[k]^v[i])&~(v[i]^v[j]))==0]
  if not middle:prim.append((i,j))
  strict[i][j]=strict[j][i]=min([d[i][j]]+[max(strict[i][k],strict[k][j]) for k in middle])
  need(ordinary[i][j]<=strict[i][j]<=d[i][j],'pair bottleneck inequality')
 need(max((d[i][j] for i,j in prim),default=0)==max(map(max,strict),default=0),'primitive exact uniform criterion')
 return d,ordinary,strict,prim

def endpoint_dag(v,x,y):
 interval=sorted((z for z in v if ((z^x)&~(x^y))==0),key=lambda z:((z^x).bit_count(),z))
 best={x:0}
 for z in interval[1:]:
  candidates=[max(value,cost(w^z)[1]) for w,value in best.items() if ((w^x)&~(z^x))==0]
  need(bool(candidates),'endpoint DAG has predecessor')
  best[z]=min(candidates)
 return best[y]

def component(v,bound):
 left=set(v);out=[]
 while left:
  x=min(left);left.remove(x);todo=[x];part=[x]
  while todo:
   x=todo.pop();nxt=[y for y in left if cost(x^y)[1]<=bound]
   for y in nxt:left.remove(y);todo.append(y);part.append(y)
  out.append(sorted(part))
 return out

def route(v,x,y,bound):
 prev={x:None};todo=deque([x])
 while todo:
  a=todo.popleft()
  if a==y:break
  for b in v:
   if b not in prev and cost(a^b)[1]<=bound:prev[b]=a;todo.append(b)
 need(y in prev,'ordinary route exists');out=[y]
 while out[-1]!=x:out.append(prev[out[-1]])
 return out[::-1]

def main():
 old,n=load();fib,gd,ch,groups=enumerate_states(n);base=n['base'];selected=1394;need(base[0]==selected,'selected bit index')
 # Literal owner-clock relation on the two cells' existing record ports.
 pairs=[]
 for bit in (0,1):pairs.append([n['literal'][1457][bit]['global_path'][1],n['literal'][selected][bit]['global_path'][2]])
 need(pairs==[[194,162],[191,159]],'literal two-entry shared-port table')
 for a in (0,1):
  for b in (0,1):
   ports=[n['literal'][1457][a]['global_path'][1],n['literal'][selected][b]['global_path'][2]]
   need((ports in pairs)==(a==b),'literal relation iff target bits agree')
   # K has the same path; converse reverses the actual eight-entry path.
   reflected=[n['literal'][1457][a]['global_path'][::-1][6],n['literal'][selected][b]['global_path'][::-1][5]]
   need(reflected==ports,'same physical port values on converse at reflected clocks')
 # Core private frames must not be silently identified.
 need(not set(n['literal'][1457][0]['global_path'])&set(n['literal'][selected][0]['global_path']),'no physical owner identification')
 rf=restricted(fib,0);summary={k:{'states':0,'charges':0,'pairs':0,'primitive_pairs':0,'ordinary_clocks':0,'strict_clocks':0} for k in range(20)}
 digest=sha256();catalog={};phist=Counter();ohist=Counter();shist=Counter();critical=Counter();lower={};per_fibre=[];new_witness=None;independent_pair_count=0
 for (k,q),v in rf.items():
  d,o,s,pr=solve(v);od,oo,ss,pp,_=old.solve_fibre(v)
  need(d==od and o==oo and s==ss and set(pr)==pp,'independent Floyd/interval and inherited MST/interval agree on every pair')
  row=summary[k];row['states']+=len(v);row['charges']+=1;row['pairs']+=len(v)*(len(v)-1)//2;row['primitive_pairs']+=len(pr)
  ob=max(map(max,o));sb=max(map(max,s));row['ordinary_clocks']=max(row['ordinary_clocks'],ob);row['strict_clocks']=max(row['strict_clocks'],sb)
  per_fibre.append([k,q,len(v),ob,sb,len(pr)])
  for i in range(len(v)):
   for j in range(i+1,len(v)):
    independent_pair_count+=1;ohist[o[i][j]]+=1;shist[s[i][j]]+=1
    need(endpoint_dag(v,v[i],v[j])==s[i][j],'independent directed endpoint-subset DAG optimum')
    digest.update(canon([k,q,v[i],v[j],joint(o[i][j]),joint(s[i][j])]))
    if s[i][j]==34:critical[k]+=1
  for i,j in pr:
   x,y=v[i],v[j];common=x&y;a,b=sorted((x^common,y^common));c=d[i][j];phist[c]+=1
   if (a,b) not in catalog:
    original=fib[(gd[a],ch[a])];middle=[z for z in original if z not in (a,b) and ((z^a)&~(a^b))==0]
    catalog[(a,b)]={'masks':[a,b],'cost':list(cost(a^b)),'originally_primitive':not middle,'original_interval_interior':middle}
    if middle and (new_witness is None or c>new_witness['cost'][1]):new_witness=catalog[(a,b)]
  if k in (8,9,10,11) and ob==34 and k not in lower:
   # A complete lower fibre and its actual bottleneck cut, not just a crossing minimum.
   cut=component(v,30);need(len(cut)>1,'sharp grade lower cut')
   x=cut[0][0];y=cut[1][0];ordinary_old=old.solve_fibre(fib[(k,q)])[1];vi=fib[(k,q)];old_bound=ordinary_old[vi.index(x)][vi.index(y)]
   lower[k]={'grade':k,'charge':q,'complete_restricted_fibre':v,'components_below34':cut,'endpoints':[x,y],'old_endpoint_optimum':joint(old_bound),'old_optimal_route':route(vi,x,y,old_bound),'new_optimal_route':route(v,x,y,34)}
 # Independent branch-conditioned polynomial and exact holes.
 B=[s for s in base if s!=selected];labels=[n['labels'][s] for s in B];dp={(0,0):1}
 for label in labels:
  nxt=dict(dp)
  for (k,q),count in dp.items():key=(k+1,add(q,label));nxt[key]=nxt.get(key,0)+count
  dp=nxt
 q4=total(n['labels'][s] for s in [selected,1457,1841,3015]);need(q4==1009,'native four-block charge')
 support=[set() for _ in range(20)];branch=[];dpfull=Counter()
 for a in (0,1):
  states=[x for v in rf.values() for x in v if x>>16==a]
  for k in range(20):
   vals={add(q,q4 if a else 0) for (r,q) in dp if r+4*a==k};support[k]|=vals
   actual={ch[x] for x in states if gd[x]==k};need(vals==actual,'every literal branch-conditioned charge support')
   branch.append({'branch':a,'grade':k,'support_count':len(vals),'holes':sorted(Q-vals),'support_sha256':sha256(sum(1<<q for q in vals).to_bytes(256,'little')).hexdigest()})
  for (k,q),count in dp.items():dpfull[(k+4*a,add(q,q4 if a else 0))]+=count
 need(dpfull==Counter({key:len(v) for key,v in rf.items()}),'weighted B plus four-block polynomial exactly counts every fibre')
 holes={str(k):sorted(Q-support[k]) for k in range(20)}
 for k in range(20):
  hh=[set(row['holes']) for row in branch if row['grade']==k];need(set(holes[str(k)])==hh[0]&hh[1],'conditional holes are intersection, not unconditional old supports')
 # New primitive generating function independently counts contexts.
 model=Counter();normalized_hist=Counter();new_hist=Counter()
 for (a,b),record in catalog.items():
  d=a^b;t=d>>16;bs=(d&65535).bit_count()-t;s,c=cost(d);normalized_hist[c]+=1
  if not record['originally_primitive']:new_hist[c]+=1
  need((a&1)==a>>16 and (b&1)==b>>16,'normalization retains shared relation')
  for j in range(16-bs):
   for z in range(2-t):model[(s//2+j+4*z,c)]+=comb(15-bs,j)
 observed=Counter()
 for (k,q),v in rf.items():
  for i,x in enumerate(v):
   for y in v[i+1:]:
    if not any(z not in (x,y) and ((z^x)&~(x^y))==0 for z in v):observed[(k,cost(x^y)[1])]+=1
 need(model==observed,'normalized primitive-context polynomial equals direct pair census')
 nmax=normalized_hist[34];need(nmax==100 and dict(critical)=={8:100,9:300,10:300,11:100},'critical polynomial 100 z8 (1+z)^3')
 # All sixteen edge choices: same topology, unequal charge holes.
 census=[]
 for idx,source in enumerate(base):
  ff=restricted(fib,idx);hist={k:{'charges':0,'ordinary_clocks':0,'strict_clocks':0} for k in range(20)}
  for (k,q),v in ff.items():
   _,oo,ss,_=solve(v);r=hist[k];r['charges']+=1;r['ordinary_clocks']=max(r['ordinary_clocks'],max(map(max,oo)));r['strict_clocks']=max(r['strict_clocks'],max(map(max,ss)))
  need(all(hist[k]['ordinary_clocks']==hist[k]['strict_clocks']==34 for k in range(8,12)),'every one-edge link central threshold')
  need(not any(r['charges']==2048 for r in hist.values()),'no linked source retains a complete grade')
  census.append({'source':source,'grades':hist})
 # Literal whole-packet arrivals for every lower witness and the newly primitive trade.
 literal={}
 for k,w in lower.items():
  literal[str(k)]={'before_restriction':old.literal_route(n,groups,w['old_optimal_route']),'after_restriction':old.literal_route(n,groups,w['new_optimal_route'])}
  need(any(((x&1)!=(x>>16)) for x in w['old_optimal_route']) if k in (9,10) else True,'new relation removes at least one cheap old-route arrival')
 if new_witness:
  literal['new_primitive_direct']=old.literal_route(n,groups,new_witness['masks'])
  a,b=new_witness['masks'];need([a,b]==[53028,77843],'frozen newly primitive witness')
  original=fib[(gd[a],ch[a])];di,oo,ss,pp,_=old.solve_fibre(original);i,j=original.index(a),original.index(b)
  need(oo[i][j]==ss[i][j]==18,'new primitive old strict optimum8/18')
  path=[a,67124,b];need(not ((a^67124)&(67124^b)),'old two-event route has disjoint target-directed switches')
  literal['new_primitive_before_restriction']=old.literal_route(n,groups,path)
  new_witness.update({'grade':gd[a],'charge':ch[a],'old_ordinary_optimum':joint(oo[i][j]),'old_strict_optimum':joint(ss[i][j]),'old_strict_route':path})
 # Equality tree, redundant cycle, and inconsistent signed cycle controls.
 cycle_counts=Counter()
 for x in range(16):
  a,b,c,d=[(x>>j)&1 for j in range(4)];tree=a==b and b==c and a==d;redundant=tree and c==d;twisted=tree and c!=d
  need(tree==redundant,'redundant equality cycle changes no assignment');cycle_counts['tree']+=tree;cycle_counts['redundant_cycle']+=redundant;cycle_counts['twisted_cycle']+=twisted
 need(dict(cycle_counts)=={'tree':2,'redundant_cycle':2,'twisted_cycle':0},'signed cycle control counts')
 # Treewidth-one equal-margin hostile: alternating physical-slot path.
 abstract=[]
 for m in range(1,9):
  v=[x for x in range(1<<(2*m)) if all(((x>>j)^(x>>(j+1)))&1 for j in range(2*m-1)) and x.bit_count()==m]
  need(len(v)==2 and (v[0]^v[1]).bit_count()==2*m,'equal-margin path has exactly two antipodal states')
  abstract.append({'pairs':m,'states':v,'fixed_population':m,'forced_slots':2*m,'clock_cost_with_two_per_slot':4*m})
 return {'status':'FINITE-EXACT authenticated nineteen-slot shared-record-port restriction; abstract and conditional extensions scoped separately','checks':CHECKS,'literal_relation':{'cell1':[1457,1841,3015],'cell2':[1457,1394],'shared_existing_source':1457,'ports':[[1457,1],[1394,2]],'allowed_global_owner_pairs':pairs,'converse_clocks':[6,5],'four_block_charge':q4,'four_block_clock_cost':10},'state_count':sum(len(v) for v in rf.values()),'nonempty_fibres':len(rf),'endpoint_pairs':independent_pair_count,'largest_fibre':max(map(len,rf.values())),'summary':summary,'pair_digest':digest.hexdigest(),'ordinary_histogram':dict(sorted(ohist.items())),'strict_histogram':dict(sorted(shist.items())),'primitive_endpoint_histogram':dict(sorted(phist.items())),'normalized_primitive_histogram':dict(sorted(normalized_hist.items())),'newly_primitive_histogram':dict(sorted(new_hist.items())),'primitive_catalogue':list(catalog.values()),'new_primitive_witness':new_witness,'critical_by_grade':dict(critical),'exact_holes':holes,'branch_conditioned_holes':branch,'all_link_census':census,'sharp_grade_witnesses':lower,'literal_full_packet_routes':literal,'cycle_controls':dict(cycle_counts),'abstract_equal_margin_controls':abstract,'fibre_summaries':per_fibre}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(main(),sort_keys=True));path=P/'certificate.json'
 if args.write:path.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 else:need(out==json.loads(path.read_text()),'complete frozen certificate equality')
 print(json.dumps({'status':'PASS','states':out['state_count'],'fibres':out['nonempty_fibres'],'pairs':out['endpoint_pairs'],'normalized_primitives':len(out['primitive_catalogue']),'newly_primitive':sum(out['newly_primitive_histogram'].values()),'pair_digest':out['pair_digest'],'checks':out['checks']},sort_keys=True))
