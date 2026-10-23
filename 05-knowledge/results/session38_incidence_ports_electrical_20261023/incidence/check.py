"""Finite incidence menu and transport audit. No solver; no removable assertions.
Run with --classes classes_n8.tsv --ports ports.tsv.gz --capacities endpoint_capacities.tsv.gz.
--write deliberately regenerates receipt; default compares frozen receipt.
"""
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
import argparse,csv,gzip,hashlib,json,math,heapq
HERE=Path(__file__).resolve().parent
checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise AssertionError(why)
def matrix(code,n):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  A[i][j]=(code>>(n*(n-1)//2-1-k))&1;A[j][i]=1-A[i][j]
 return A
def stats(A):
 n=len(A);c=[0]*n;E=[[0]*n for _ in range(n)];v4=[0]*n
 for i,j,k in combinations(range(n),3):
  if A[i][j]==A[j][k]==A[k][i]:
   for u in (i,j,k):c[u]+=1
   for u,v in ((i,j),(i,k),(j,k)):E[u][v]+=1;E[v][u]+=1
 for B in combinations(range(n),4):
  if sorted(sum(A[u][v] for v in B) for u in B)==[1,1,2,2]:
   for u in B:v4[u]+=1
 out=[[v for v in range(n) if A[u][v]] for u in range(n)]
 for u in range(n):
  need(c[u]==sum(E[u][v] for v in out[u]),'outgoing cyclic incidence balance')
  need(c[u]==sum(E[u][v] for v in range(n) if A[v][u]),'incoming cyclic incidence balance')
 return c,E,v4,out
def esig(c,E,out,p=None,q=None):
 return tuple(sorted(tuple(sorted((c[v] if p is None else c[v]%p,E[u][v] if q is None else E[u][v]%q) for v in out[u])) for u in range(len(c))))
def get_signatures(A):
 c,E,v4,out=stats(A);n=len(A);s=list(map(len,out))
 def N(tag):return tuple(sorted((tag[u],tuple(sorted(tag[v] for v in out[u]))) for u in range(n)))
 J=N(list(zip(s,c)));J4=N(list(zip(s,c,v4)));q=[sum(E[u][v]*(E[u][v]-1)//2 for v in out[u]) for u in range(n)]
 M=sum(q[u]*q[v]*q[v] for u in range(n) for v in out[u]);eh=tuple(sorted(E[u][v] for u in range(n) for v in out[u]));M2=sum(e*(e-1)//2 for e in eh)
 D={'J':J,'J_M2':(J,M2),'J4':J4,'J4_edgehist':(J4,eh),'J4_Q':(J4,M),'E':esig(c,E,out),'E43':esig(c,E,out,4,3),'E42':esig(c,E,out,4,2),'E_neighbor2':esig(c,E,out,2,None),'E_neighbor3':esig(c,E,out,3,None),'E_edge2':esig(c,E,out,None,2),'E_no_neighbor':tuple(sorted(tuple(sorted(E[u][v] for v in out[u])) for u in range(n))),'E_no_edge':tuple(sorted(tuple(sorted(c[v] for v in out[u])) for u in range(n))),'E_unpaired':tuple(sorted((tuple(sorted(c[v] for v in out[u])),tuple(sorted(E[u][v] for v in out[u]))) for u in range(n))),'E_global':tuple(sorted((c[u],c[v],E[u][v]) for u in range(n) for v in out[u]))}
 return D,(c,E,v4,out,M)
def peeling(rowcap,colcap,allowed,actual):
 r=dict(rowcap);c=dict(colcap);rn=defaultdict(set);cn=defaultdict(set);queue=[];fixed={};trace=[]
 for a,b in allowed:rn[a].add(b);cn[b].add(a)
 for a in r:
  if len(rn[a])==1:heapq.heappush(queue,(0,a))
 for b in c:
  if len(cn[b])==1:heapq.heappush(queue,(1,b))
 def remove(a,b):
  rn[a].discard(b);cn[b].discard(a)
  if r[a] and len(rn[a])==1:heapq.heappush(queue,(0,a))
  if c[b] and len(cn[b])==1:heapq.heappush(queue,(1,b))
 while queue:
  side,v=heapq.heappop(queue)
  if side==0:
   if not r[v] or len(rn[v])!=1:continue
   a=v;b=next(iter(rn[a]));z=r[a]
  else:
   if not c[v] or len(cn[v])!=1:continue
   b=v;a=next(iter(cn[b]));z=c[b]
  need(z<=r[a] and z<=c[b],'capacity conflict');need(actual.get((a,b),0)==z,'forced cell mismatch')
  fixed[a,b]=z;r[a]-=z;c[b]-=z;trace.append((side,a,b,z))
  if not r[a]:
   for bb in list(rn[a]):remove(a,bb)
  if not c[b]:
   for aa in list(cn[b]):remove(aa,b)
 rr={a:z for a,z in r.items() if z};cc={b:z for b,z in c.items() if z};edges=sorted((a,b) for a in rr for b in rn[a])
 need(sum(rr.values())==sum(cc.values()),'residual mass mismatch')
 return fixed,trace,rr,cc,edges

def main():
 p=argparse.ArgumentParser();p.add_argument('--classes',type=Path,default=HERE/'inherited/classes_n8.tsv');p.add_argument('--ports',type=Path,default=HERE.parent/'session38-wiring/ports.tsv.gz');p.add_argument('--capacities',type=Path,default=HERE.parent/'session38-wiring/endpoint_capacities.tsv.gz');p.add_argument('--write',action='store_true');args=p.parse_args()
 R=sorted(csv.DictReader(args.classes.open(),delimiter='\t'),key=lambda r:int(r['canon']));need(len(R)==6880,'class universe size');need(len({r['canon'] for r in R})==6880,'duplicate class canon');need(sum(math.factorial(8)//int(r['Aut']) for r in R)==1<<28,'declared orbit mass')
 groups=defaultdict(lambda:defaultdict(list));modgroups={(p,q):set() for p in range(1,9) for q in range(1,8)};sigby=[];hostiles=[]
 for cid,row in enumerate(R):
  A=matrix(int(row['canon']),8);D,(c,E,v4,out,M)=get_signatures(A);sigby.append(D)
  for k,v in D.items():groups[k][v].append(cid)
  for (p,q),S in modgroups.items():
   # Exact self-delimiting byte serialization: row length then sorted pair-code c*q+e.
   S.add(b''.join(sorted(bytes([len(out[u])]+sorted((c[v]%p)*q+(E[u][v]%q) for v in out[u])) for u in range(8))))
  if int(row['canon']) in (27467784,27467793,6427781,6427778):hostiles.append({'class_id':cid,'canon':int(row['canon']),'rows':[sum(x<<i for i,x in enumerate(r)) for r in A],'H':int(row['H']),'Q':M,'edge_histogram':dict(sorted(Counter(E[u][v] for u in range(8) for v in out[u]).items()))})
 O={'scope':'FINITE-EXACT all6880 ordinary n8 classes, plus complete supplied fixed-path ledger; class exhaustiveness independently audited elsewhere','inputs':{role:hashlib.sha256(path.read_bytes()).hexdigest() for role,path in [('classes',args.classes),('ports',args.ports),('capacities',args.capacities)]},'class_menu':{},'modulus_menu':[[p,q,len(S)] for (p,q),S in modgroups.items()],'J4_hostiles':hostiles}
 for k,G in groups.items():
  cg=[[int(R[i]['canon']) for i in ids] for ids in G.values() if len(ids)>1]
  O['class_menu'][k]={'distinct':len(G),'collision_groups':cg,'largest_fiber':max(map(len,G.values()))}
 need(len(groups['E43'])==6880,'E43 failed');need(len(groups['E'])==6880,'E failed');need(len(groups['J4_Q'])==6880,'J4Q failed');need(len(groups['J4'])==6878,'inherited J4 census changed');need(len(groups['J4_edgehist'])==6878,'edge histogram unexpected')
 # Signature IDs are canonical-order first-occurrence IDs, not polynomial hashes.
 gids={};colliders={};members={}
 for name,G in groups.items():
  gids[name]={};members[name]={}
  for gid,ids in enumerate(G.values()):
   members[name][gid]=ids
   for cid in ids:gids[name][cid]=gid
  colliders[name]={cid for ids in G.values() if len(ids)>1 for cid in ids}
 with (HERE/'signature_ids.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['class_id','canon']+list(groups))
  for cid,r in enumerate(R):w.writerow([cid,r['canon']]+[gids[k][cid] for k in groups])
 models=('J','J_M2','J4','J4_Q','E','E43','E_unpaired','E_edge2','E_neighbor2','E_neighbor3')
 union=set().union(*(colliders[nm] for nm in models));bytarget=defaultdict(list);target_entries=Counter();target_mass=Counter();allmass=0;allentries=0;replayed_open=Counter()
 with gzip.open(args.ports,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   i,j,a,b,blue,z=(int(row[k]) for k in ('source_class','target_class','start_score','end_score','blue','multiplicity'));allmass+=z;allentries+=1;target_entries[j]+=1;target_mass[j]+=z;replayed_open[j,8-a,6-b,blue]+=z
   if j in union:bytarget[j].append((i,j,a,b,blue,z))
 need(allmass==1<<20,'closed mass wrong');need(allentries==956339,'native entry count changed')
 caps={};declared_open={}
 with gzip.open(args.capacities,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   j,a,b,z=(int(row[k]) for k in ('class_id','start_score','end_score','open'));blue=int(row['blue'] in ('1','True','true'))
   if z:declared_open[j,a,b,blue]=z
   if j in union and z:caps[j,a,b,blue]=z
 need(dict(replayed_open)==declared_open,'full open endpoint capacity ledger differs from replay')
 O['transport']={}
 for name in models:
  coll=colliders[name];rows={};rowkeys=[];rowcap=Counter();act0=Counter()
  for j in sorted(coll):
   for i,j,a,b,blue,z in bytarget[j]:
    key=(i,a,b,blue,gids[name][j])
    if key not in rows:rows[key]=len(rowkeys);rowkeys.append(key)
    ri=rows[key];rowcap[ri]+=z;act0[ri,(j,8-a,6-b,blue)]+=z
  colkeys=sorted(k for k in caps if k[0] in coll);colids={k:i for i,k in enumerate(colkeys)};colcap={i:caps[k] for i,k in enumerate(colkeys)};actual={(ri,colids[k]):z for (ri,k),z in act0.items()};allow=[]
  for ri,(i,a,b,blue,gid) in enumerate(rowkeys):
   for j in members[name][gid]:
    ck=(j,8-a,6-b,blue)
    if ck in colids:allow.append((ri,colids[ck]))
  need(set(actual)<=set(allow),'actual entry not allowed');need(sum(rowcap.values())==sum(colcap.values()),'column capacities differ from replay');need(sum(actual.values())==sum(rowcap.values()),'actual mass differs')
  fixed,tr,rr,cc,edges=peeling(rowcap,colcap,allow,actual);unq_entries=sum(target_entries[j] for j in range(6880) if j not in coll)
  O['transport'][name]={'class_signatures':len(groups[name]),'row_count':unq_entries+len(rows),'candidate_entries':unq_entries+len(allow),'forced_entries':unq_entries+len(fixed),'unresolved_mass':sum(rr.values()),'unresolved_rows':len(rr),'unresolved_columns':len(cc),'residual_allowed_entries':len(edges),'ambiguous_target_classes':len(coll),'peeling_sha256':hashlib.sha256(json.dumps(tr,separators=(',',':')).encode()).hexdigest()}
  if rr:
   residual={'signature':name,'row_keys':[[a,list(rowkeys[a]),z] for a,z in sorted(rr.items())],'column_keys':[[b,list(colkeys[b]),z] for b,z in sorted(cc.items())],'allowed_entries':[[a,b,actual.get((a,b),0)] for a,b in edges],'class_canon':{str(j):int(R[j]['canon']) for j in sorted(coll)},'note':'Actual values retain original mass on residual cells; all removed incident cells are zero or forced. Row key=source,start,end,blue,signatureID; column=target,start,end,blue.'}
   (HERE/f'residual_{name}.json').write_text(json.dumps(residual,separators=(',',':'))+'\n')
  print(name,O['transport'][name],flush=True)
 O['checks']=checks
 O=json.loads(json.dumps(O));out=HERE/'receipt.json'
 if args.write:out.write_text(json.dumps(O,indent=2)+'\n')
 else:need(json.loads(out.read_text())==O,'frozen receipt differs')
 print('checks',checks,'receipt', 'wrote' if args.write else 'matched',flush=True)
if __name__=='__main__':main()
