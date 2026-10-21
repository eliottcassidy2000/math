#!/usr/bin/env python3
"""Exact source-owned fixed-frontier credits and weighted interval handoffs.
Standard library only. Source inputs authenticated under --repo-root.
No convergence or unbounded-annulus claim.
"""
import argparse,json
from pathlib import Path
from collections import Counter,defaultdict
from hashlib import sha256
HERE=Path(__file__).resolve().parent
CHECKS=0
L=((1,2),(-2,-1),(1,-1)); U=((0,1),(-1,-1),(1,0))
PINS={
 'session18_native_annulus/frontier.json':'391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f',
 'annulus_residual_flow_20261020k/receipt.json':'35974096cbf893fb209ce5c47d915a9837b4a048ac7579b1d0806b46e106ed4d',
 'residual_external_credit_20261020k/receipt.json':'d702d85fd51e0f9778ac1fb9a784a1426b3af76ea4aa0cb1d5a32392711eef82',
 'annulus_block_completion_20261020j/receipt.json':'b85cfd8a95962c827d9e21d1418a8fe7f4086c6cca6203f042f547d32ad3b609',
 'sibling_annulus_guards_20261020i/README.md':'1f4964733fd8d0f85d3b17b33628d739972f03d74687affb1973c695ed65adef',
 'session18_native_phase_lifts/README.md':'20d2cf37ed149345e259b1107509b7a9480de709461bad326df57b396f1b61ff'}
def ck(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(msg)
def integer(x,lo=0):
 if type(x) is not int or x<lo:raise ValueError('literal integer domain')
 return x
def enc(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def source(s):return type(s) is int and s>0 and s%2==1 and s%8!=5
def Q(v):return v[0]**2-v[0]*v[1]+v[1]**2
def beta(v,w):return 2*v[0]*w[0]-v[0]*w[1]-v[1]*w[0]+2*v[1]*w[1]
def rotate(v,k):
 for _ in range(k%3):v=(-v[1],v[0]-v[1])
 return v
def edge(s):
 if not source(s):raise ValueError('primitive source')
 a=2 if s%8==1 else 1;t=a;u=(3*s+1)//2**a
 while u%8==5:u=(u-1)//4;t+=2
 ck(source(u),'primitive native edge endpoint')
 return u,t

def charge(N,s):
 integer(N,1)
 if not source(s):raise ValueError('primitive charge source')
 a=b=0
 while s<4*N:
  x,y=(L if s<2*N else U)[s%3];a+=x;b+=y;s=4*s+1
 return a,b

def independent_charge(N,s):
 a=b=0;K=(0,-3,3);T=(-1,-1,2)
 while s<4*N:
  row=K if s<2*N else T;a+=row[s%3];b+=row[(s+1)%3];s=4*s+1
 ck((a-b)%3==0 and (a+2*b)%3==0,'three-row divisibility')
 return (a-b)//3,-(a+2*b)//3

def inverse(u,a,j):
 integer(u,1);integer(a,1);integer(j)
 if not source(u) or a not in (1,2):raise ValueError('inverse parent and branch')
 z=2**a*(3*u+1)*4**j-(2**a+3)
 if z%9:raise ValueError('inverse integrality')
 s=z//9
 if not source(s) or (2 if s%8==1 else 1)!=a:raise ValueError('inverse literal first branch')
 return s

def lifted_mod(u,a,j,k,m):
 integer(k);integer(m,1);inverse(u,a,j)
 z=(2**a*(3*u+1)*pow(4,j+3*k,9*m)-(2**a+3))%(9*m)
 ck(z%9==0,'full-denominator register')
 return z//9

class Routes:
 def __init__(self,data):
  self.V=set(data['frontier']);ranks=dict(data['frontier_ranks'])
  ck(len(self.V)==627 and max(self.V)==63665,'frozen frontier size and cutoff')
  ck(set(ranks)==self.V and ranks[1]==0,'frontier rank universe')
  for s in sorted(self.V):
   if s!=1:
    u,t=edge(s);ck(u in self.V and ranks[u]<ranks[s],'grounded native rank')
  self.cache={s:(s,0,0)for s in self.V}
 def get(self,s):
  if not source(s):raise ValueError('route source')
  seen=set();trail=[];u=s
  while u not in self.cache:
   ck(u not in seen and len(trail)<10000,'finite route not silently omitted')
   seen.add(u);v,t=edge(u);trail.append((u,t));u=v
  terminal,t,d=self.cache[u]
  for u,c in reversed(trail):t+=c;d+=1;self.cache[u]=(terminal,t,d)
  result=self.cache[s];ck(result[2]<=256,'declared finite horizon cap');return result
 def row(self,N,s):
  u,t,d=self.get(s);v=charge(64*N,s)
  ck(v==independent_charge(64*N,s),'independent source charge')
  ck(s%15%5==s%35%5 and (s%15+15*((s%35-s%15)%7))%105==s%105,'joint prime register')
  return [s,s%15,s%35,u,t,d,*v,Q(v)]

def catalogue(N,routes):
 integer(N,1);out=[]
 for u in sorted(routes.V):
  if u>=4*N:continue
  for a in (1,2):
   j=(3-a-u)%3;s=inverse(u,a,j)
   while s<4*N:s=64*s+7*(2**a+3);j+=3
   ck(s<256*N,'each old parent branch has one annular candidate')
   ck(s==inverse(u,a,j) and edge(s)==(u,a+2*j),'literal candidate edge and clock')
   if s in routes.V:continue
   old=routes.row(N,u);new=routes.row(N,s)
   ck(old[3:6]==[u,0,0] and new[3:6]==[u,a+2*j,1],'first-arrival credit cutoff')
   ck(tuple(old[6:8])==charge(N,u),'fixed old source charge')
   sg=1-2*(a%2);w=tuple(sg*x for x in new[6:8]);b=beta(tuple(old[6:8]),w)
   out.append({'old':u,'new':s,'branch':a,'strippings':j,'beta':b,'old_row':old,'new_row':new})
 ck(len({(r['old'],r['new'])for r in out})==len(out),'unique original pair addresses')
 ck(len(out)<=1254,'fixed-frontier address cap')
 return out

def source_bank(N,routes):
 blocks=[];rows={};used=set()
 for band in range(5):
  lo,hi=(4<<band)*N,(8<<band)*N;runs=[];run=[]
  for x in range(lo+(1-lo)%8,hi,8):
   if x in routes.V or 2*x+1 in routes.V:
    if run:runs.append(run);run=[]
   else:run.append(x)
  if run:runs.append(run)
  for run in runs:
   for j in range(0,len(run)-2,3):
    terms=[];ids=[]
    for x in run[j:j+3]:
     y=2*x+1;p=routes.row(N,x);q=routes.row(N,y);rows[x]=p;rows[y]=q
     ck(p[3]==q[3] and p[5]==q[5] and q[4]==p[4]+1,'actual sibling route and clock')
     terms.append((-beta(tuple(p[6:8]),tuple(q[6:8])),p[5],x,y));ids.extend((x,y))
    ck(not used.intersection(ids),'sibling source disjointness');used.update(ids)
    pos=[d for w,d,x,y in terms if w==3];neg=[d for w,d,x,y in terms if w==-3]
    if not pos and not neg:sign=start=end=0
    else:
     ck(len(pos)==len(neg)==1,'native interval alphabet');start,end=sorted((pos[0],neg[0]));sign=0 if start==end else (1 if pos[0]<neg[0] else -1)
    blocks.append({'id':len(blocks),'band':band,'terms':tuple(terms),'sign':sign,'start':start,'end':end})
 return blocks,rows

def select_exact(rows,target):
 negative=sorted((r for r in rows if r['beta']<0),key=lambda r:(r['beta'],r['old'],r['branch']))
 selected=[];K=0
 for r in negative:
  w=-r['beta']
  if K+w<=target:selected.append(r);K+=w
  if K==target:break
 ck(K==target,'literal exact permanent credit budget')
 largest=sum(-r['beta']for r in negative[:len(selected)-1])
 ck(largest<target,'minimum number of owners in complete catalogue')
 return selected,largest

def allocate(blocks,credits,D=257):
 """Unit indices partition real weights, not copies of a pair term."""
 jobs={(b['id'],k):(b['start'],b['end'])for b in blocks if b['sign']==1 for k in range(3)}
 owners={(0,b['id'],k):(b['start'],b['end'])for b in blocks if b['sign']==-1 for k in range(3)}
 owners.update({(1,i,k):(1,D)for i,r in enumerate(credits)for k in range(-r['beta'])})
 mapping={};began={};segments=[]
 for d in range(D+1):
  active={p for p,(a,z)in jobs.items()if a<=d<z}
  available={n for n,(a,z)in owners.items()if a<=d<z}
  for p,n in list(mapping.items()):
   if p not in active or n not in available:
    segments.append([*p,*n,began.pop(p),d]);del mapping[p]
  pending=sorted(active-set(mapping));free=sorted(available-set(mapping.values()))
  ck(len(pending)<=len(free),'actual source-owned cut capacity')
  for p,n in zip(pending,free):mapping[p]=n;began[p]=d
 ck(not mapping,'all positive intervals finished')
 segments.sort();audit_segments(jobs,owners,segments)
 return jobs,owners,segments

def audit_segments(jobs,owners,segments):
 pp=defaultdict(list);nn=defaultdict(list)
 for seg in segments:
  ck(type(seg) is list and len(seg)==7 and all(type(x)is int for x in seg),'literal handoff segment')
  p=tuple(seg[:2]);n=tuple(seg[2:5]);a,z=seg[5:]
  ck(p in jobs and n in owners,'source-addressed unit owner')
  ck(jobs[p][0]<=a<z<=jobs[p][1] and owners[n][0]<=a<z<=owners[n][1],'both actual birth/death guards')
  pp[p].append((a,z));nn[n].append((a,z))
 ck(set(pp)==set(jobs),'every positive unit owned')
 for p,parts in pp.items():
  parts.sort();ck(parts[0][0]==jobs[p][0] and parts[-1][1]==jobs[p][1],'whole interval endpoints')
  ck(all(parts[i][1]==parts[i+1][0]for i in range(len(parts)-1)),'whole positive interval exact partition')
 for parts in nn.values():
  parts.sort();ck(all(parts[i][1]<=parts[i+1][0]for i in range(len(parts)-1)),'no simultaneous negative-unit reuse')

def reject(fn):
 try:fn()
 except ValueError:return True
 return False

def run(root):
 inputs={}
 for rel,h in PINS.items():
  b=(root/'05-knowledge/results'/rel).read_bytes().replace(b'\r\n',b'\n');ck(sha256(b).hexdigest()==h,'authenticated source '+rel)
  if rel.endswith('.json'):inputs[rel]=json.loads(b)
 routes=Routes(inputs['session18_native_annulus/frontier.json'])
 flow=inputs['annulus_residual_flow_20261020k/receipt.json'];cases=[]
 for N,target,number in ((4096,273,79),(8192,474,154)):
  prior=next(c for c in flow['native_cases']if c['N']==N);rows=catalogue(N,routes);selected,largest=select_exact(rows,target)
  ck(len(selected)==number,'minimal pair-owner count')
  blocks,srows=source_bank(N,routes)
  ck(sha256(enc(blocks)).hexdigest()==prior['block_sha256'],'exact frozen sibling bank')
  ck(sha256(enc([srows[s]for s in sorted(srows)])).hexdigest()==prior['source_rows_sha256'],'complete actual sibling source rows')
  J=[]
  for d in range(257):
   j=sum(3*b['sign']for b in blocks if b['start']<=d<b['end']);J.append(j)
  ck(J==prior['original_bank_debit'] and max(J)==target and J[0]==0,'sharp all-cut permanent capacity')
  jobs,owners,segs=allocate(blocks,selected)
  ck(reject(lambda:audit_segments(jobs,owners,segs+[segs[0]])),'duplicate handoff unit hostile')
  premature=[list(segs[0])];premature[0][5]=0
  ck(reject(lambda:audit_segments(jobs,owners,premature)),'premature incomplete handoff hostile')
  record={'N':N,'catalogue_count':len(rows),'negative_count':sum(r['beta']<0 for r in rows),'available_credit':sum(-r['beta']for r in rows if r['beta']<0),
   'weight_histogram':dict(sorted(Counter(r['beta']for r in rows).items())),'extracted_credit':target,'selected_count':len(selected),'largest_one_fewer_capacity':largest,
   'selected_credits':selected,'catalogue':rows,'block_sha256':prior['block_sha256'],'source_rows_sha256':prior['source_rows_sha256'],'completed_positive_intervals':len(jobs)//3,
   'peak_deficit_horizons':[d for d,j in enumerate(J)if j==target],'original_bank_debit':J,'extended_bank_debit':[j-(target if d else 0)for d,j in enumerate(J)],
   'handoff_segments':segs,'handoff_schema':['positive_block','positive_unit','owner_kind_0_sibling_1_external','owner_id','owner_unit','start','end'],
   'handoff_count':len(segs)}
  if N==4096:
   old=inputs['residual_external_credit_20261020k/receipt.json'];negative={(r['old'],r['new'])for r in rows if r['beta']<0}
   ck(all(tuple(c[:2])in negative for c in old['credits']),'old94 contained not stacked')
   ledger=inputs['annulus_block_completion_20261020j/receipt.json']['full_native_control']['full_ledger'];out=[]
   for d,r in enumerate(ledger):
    C=-target if d else 0;Rp=r[6]-C;I=r[5];M=r[2]
    ck(I==r[6]+J[d]==Rp+J[d]+C,'unchanged original interaction ledger')
    ck(J[d]+C<=0,'complete sibling payment with no reused external credit')
    out.append([d,C,r[6],Rp,I,M,Rp-M])
   record['remainder_ledger']=out;record['maximum_Rprime_minus_M']=max(r[-1]for r in out)
  cases.append(record)
 # Structural controls independent of a large native orbit census.
 phases=[];phase_rows=[]
 for k in (1,2,3,4):
  N=4096*64**k;rows=catalogue(N,routes);phase_rows.append(rows)
  ck(len(rows)==1254,'stable complete frontier branches')
  phases.append({'k':k,'N':N,'negative_count':sum(r['beta']<0 for r in rows),'available_credit':sum(-r['beta']for r in rows if r['beta']<0),'weight_histogram':dict(sorted(Counter(r['beta']for r in rows).items()))})
 for before,after in zip(phase_rows,phase_rows[1:]):
  for r,t in zip(before,after):
   c=7*(2**r['branch']+3)
   ck(t['old']==r['old'] and t['branch']==r['branch'] and t['new']==64*r['new']+c,'exact complete-bank lift')
   ck(t['strippings']==r['strippings']+3 and t['new_row'][4]==r['new_row'][4]+6,'actual P-depth and unbounded clock')
   ck(tuple(t['new_row'][6:8])==rotate(tuple(r['new_row'][6:8]),c),'annular band phase rotation')
   ck(t['old_row'][6:8]==r['old_row'][6:8],'fixed parent charge')
 ck([r['beta']for r in phase_rows[0]]==[r['beta']for r in phase_rows[3]],'exact period-three scalar weights')
 ck(cases[0]['catalogue_count']!=len(phase_rows[2]),'small-scale eligibility hostile to premature periodicity')
 # Exact mod105 stratification, full-denominator source reading, and six-step return.
 register_periods=Counter()
 for r in phase_rows[0]:
  u,a,j,s=r['old'],r['branch'],r['strippings'],r['new'];c=7*(2**a+3)
  period=3 if s%5==(0 if a==1 else 2)else 6;register_periods[period]+=1;x=s
  for k in range(7):
   for m in (15,35,105):ck(lifted_mod(u,a,j,k,m)==x%m,'guard register via full denominator')
   if 0<k<period:ck(x%105!=s%105,'exact register period has no early return')
   if k==period:ck(x%105==s%105,'exact register return')
   x=64*x+c
 for N in (1,7,15916,15917,10**40+7):
  rows=catalogue(N,routes)
  if N>=15917:ck(len(rows)==1254,'large arbitrary-scale branch coverage')
 # Literal native N7 at actual P-depth8 remains hostile after mere ledger repartition.
 W=defaultdict(lambda:[0,0]);H=defaultdict(lambda:[0,0]);arrived=0;total=0
 for s in range(1,256*7,2):
  if not source(s):continue
  u,t,d=routes.get(s);v=charge(64*7,s);sg=1-2*(t%2)
  if s>=4*7:
   total+=Q(v)
   if d<=8:arrived+=Q(v)
  if d<=8:
   z=(W if s<4*7 else H)[u];z[0]+=sg*v[0];z[1]+=sg*v[1]
 oldE=sum(Q(v)for v in W.values());newE=sum(Q(v)for v in H.values());mixed=sum(beta(v,H.get(u,(0,0)))for u,v in W.items())
 hostile=[oldE,newE,mixed,arrived,total-arrived,newE+mixed-arrived,newE+mixed,total]
 ck(hostile==[19,1086,-2,1054,17,30,1084,1071],'native N7 hostile unchanged')
 # Bounded direct source scan independently verifies the inverse catalogue, not just its formulas.
 for N in range(1,33):
  direct=[]
  for s in range(4*N+1,256*N,2):
   if not source(s) or s in routes.V:continue
   u,t=edge(s)
   if u in routes.V and u<4*N:direct.append((u,s,t))
  symbolic=catalogue(N,routes)
  ck(sorted(direct)==sorted((r['old'],r['new'],r['new_row'][4])for r in symbolic),'independent bounded inverse completeness')
 for bad in (True,1.0,-1,0):ck(reject(lambda bad=bad:catalogue(bad,routes)),'scale type hostile')
 ck(reject(lambda:inverse(87,True,1)) and reject(lambda:inverse(87,1,1.0)),'inverse type hostiles')
 return {'scope':'source-owned fixed-frontier credits; minimum finite bank payment; uniform153N OPEN','source_pins':PINS,'cases':cases,'stable_ray':phases,'register_period_histogram':dict(register_periods),'native_N7_depth8':hostile,'checks':CHECKS}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--write',action='store_true');p.add_argument('--verify',action='store_true');args=p.parse_args()
 result=run(args.repo_root.resolve());data=enc(result);out=HERE/'receipt.json'
 if args.write:out.write_bytes(data)
 else:ck(out.read_bytes().replace(b'\r\n',b'\n')==data,'saved receipt exact equality')
 print(json.dumps({'checks':CHECKS,'sha256':sha256(data).hexdigest(),'cases':[{k:c[k]for k in ['N','available_credit','extracted_credit','selected_count','handoff_count','completed_positive_intervals']}for c in result['cases']],'stable_ray':result['stable_ray']}))
