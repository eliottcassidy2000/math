#!/usr/bin/env python3
"""Independent relative-owner cycle and weighted complete-menu audit."""
import argparse,csv,hashlib,itertools,json,math
from collections import Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();checks=0

def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
R=a.repo/'05-knowledge/results';nodefile=R/'metagraph_order8_20261023/merged_nodes.tsv';electrical=R/'metagraph_order8_electrical_20261023/receipt.json';nodes=list(csv.DictReader(nodefile.open(),delimiter='\t'));E=json.loads(electrical.read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0]
options={r['source']:r for r in json.loads((a.packet/'pool.json').read_text())['options']};case=json.loads((a.packet/'candidates.json').read_text())['found']['9'];ss=case['plus']+case['minus'];sign=[1]*9+[-1]*9;cost=[];reports=[]
for s,t in zip(ss,sign):
 old=options[s]['before'if t==1 else'after']['path'];new=options[s]['after'if t==1 else'before']['path'];rho=dict(zip(old,new));cycles=[];seen=set()
 for v in range(8):
  if v in seen:continue
  cyc=[];j=v
  while j not in seen:seen.add(j);cyc.append(j);j=rho[j]
  need(j==v,'relative permutation closes each cycle')
  if len(cyc)>1:cycles.append(cyc)
 distance=sum(x!=y for x,y in zip(old,new));cost.append(distance);need(sum(map(len,cycles))==distance,'cycle moved-owner support equals indexed clock distance')
 code=int(nodes[s]['canon_a']);M=[[0]*8 for _ in range(8)]
 for bit,(i,j)in zip(range(27,-1,-1),itertools.combinations(range(8),2)):M[i][j]=(code>>bit)&1;M[j][i]=1-M[i][j]
 need(any(M[i][j]!=M[rho[i]][rho[j]]for i in range(8)for j in range(8)),'relative nonidentity permutation is not source automorphism')
 entry={'source':s,'old_path':old,'new_path':new,'relative_nontrivial_cycles':cycles,'owner_clock_cost':distance}
 if s==2835:
  need(cycles==[[0,4,6]],'literal exceptional three-cycle')
  third=[rho[v]for v in new];need(third==[2,1,7,5,4,0,6,3],'literal second ambient iterate')
  failures=[[u,v]for u,v in zip(third,third[1:])if M[u][v]!=1];need(bool(failures),'second ambient iterate leaves native Hamiltonian domain')
  entry.update({'third_ambient_path':third,'failed_directed_arcs':failures})
 else:need(len(cycles)==1 and len(cycles[0])==2,'every other selected permutation is transposition')
 reports.append(entry)
need(Counter(cost)=={2:17,3:1}and sum(cost)==37,'literal cost profile')
need(sum(X[s]for s in case['plus'])==sum(X[s]for s in case['minus'])==1161,'claimed unsigned primitive potential sums')
allpoly=Counter();guardpoly=Counter();accepted=[];levels={m:Counter()for m in [2,4,8,16,32,64,128,256]};odd_index=ss.index(2835)
for mask in range(1<<18):
 chosen=[j for j in range(18)if mask>>j&1];key=(len(chosen),sum(cost[j]for j in chosen));allpoly[key]+=1
 need(key[1]%2==(mask>>odd_index)&1,'every assignment clock parity equals source2835 bit')
 if sum(sign[j]for j in chosen):continue
 par=any(sum(sign[j]*Y[h][ss[j]]for j in chosen)%2 for h in range(4));mass=sum(sign[j]*X[ss[j]]for j in chosen)
 for m,poly in levels.items():
  if m==2 or(not par and mass%(m//2)==0):poly[key]+=1
 if not par and mass%128==0:guardpoly[key]+=1;accepted.append(mask)
formula=Counter()
for k in range(18):formula[(k,2*k)]+=math.comb(17,k);formula[(k+1,2*k+3)]+=math.comb(17,k)
for key in set(allpoly)|set(formula):need(allpoly[key]==formula[key],'every unguarded bivariate polynomial coefficient')
need(sum(allpoly.values())==262144 and guardpoly=={(0,0):1,(18,37):1}and accepted==[0,262143],'complete guarded weighted enumerator')
producer=json.loads((a.packet/'cycles-normal.json').read_text());tower=[]
for m,poly in levels.items():
 saved=next(z for z in producer['dyadic_polynomial_tower']if z['modulus']==m);expected={(z['changed_slots'],z['changed_clocks']):z['coefficient']for z in saved['coefficients']}
 need(dict(poly)==expected,'every original-coordinate tower polynomial coefficient matches producer')
 need(all(poly.get((18-k,37-c),0)==n for(k,c),n in poly.items()),'full bivariate complement reciprocity')
 odd=sum(n for(k,c),n in poly.items()if c%2);need(2*odd==sum(poly.values()),'half of accepted states have odd clock cost')
 need(sum((-1)**c*n for(k,c),n in poly.items())==0,'univariate clock polynomial P(1,t) has factor1+t')
 tower.append({'modulus':m,'accepted_states':sum(poly.values()),'odd_clock_cost_states':odd,'coefficients':[{'slots':k,'clocks':c,'coefficient':n}for(k,c),n in sorted(poly.items())]})
report={'status':'PASS','checks':checks,'source_sha256':{'producer_cycles':hashlib.sha256((a.packet/'cycles-normal.json').read_bytes()).hexdigest(),'nodes':hashlib.sha256(nodefile.read_bytes()).hexdigest(),'electrical':hashlib.sha256(electrical.read_bytes()).hexdigest(),'candidates':hashlib.sha256((a.packet/'candidates.json').read_bytes()).hexdigest(),'pool':hashlib.sha256((a.packet/'pool.json').read_bytes()).hexdigest()},'sources':reports,'dyadic_weighted_polynomial_tower':tower,'unguarded_polynomial':[{'slots':k,'clocks':c,'coefficient':n}for (k,c),n in sorted(allpoly.items())],'guarded_polynomial':[{'slots':k,'clocks':c,'coefficient':n}for (k,c),n in sorted(guardpoly.items())],'scope':'Relative permutations of each source owner frame, not source automorphisms, a legal S3 action on the two-state native domain, or a permutation of a common host. Full-bank interpretation uses independently audited complete original paired-charge coordinates and direct native current audit.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'guarded_polynomial':report['guarded_polynomial'],'exceptional_source':next(z for z in reports if z['source']==2835)},sort_keys=True))
