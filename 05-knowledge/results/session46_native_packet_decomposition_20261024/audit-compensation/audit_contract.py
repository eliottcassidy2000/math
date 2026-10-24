#!/usr/bin/env python3
"""Independent exact control for the guard/reuse interpretation contract."""
import argparse,csv,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--contract',type=Path,required=True);p.add_argument('--foundation',type=Path,required=True);p.add_argument('--routes',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
checks=0
def require(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
e=a.repo/'05-knowledge/results/metagraph_order8_20261023/merged_edges.tsv.gz';f=json.loads(a.foundation.read_text());r=json.loads(a.routes.read_text());U,V=1965,2234;N=3528
require(f['input_sha256']['edges']==sha(e),'same literal original graph as domain foundation')
chi=[0]*N;chi[U]=4;image=[0]*N;odd=[];degree=0;incidences=0
for row in csv.DictReader(gzip.open(e,'rt'),delimiter='\t'):
 i,j,b,k=(int(row[t])for t in('source','target','blue','black'));w=b+k;z=w*(chi[i]-chi[j]);image[i]+=z;image[j]-=z
 if i!=j and U in(i,j):
  degree+=w;incidences+=1
  if w%2:odd.append((i,j,w))
require(not odd,'all nonloop U incidences have even conductance')
for value in image:require(value%8==0,'actual sparse character harmonic at every original physical row')
for row in f['sources']:
 s=row['source'];require(s not in(U,V),'every supplied source outside targets')
 require((chi[s]-chi[U])%8==4 and(chi[s]-chi[V])%8==0,'every U/V switch changes this separately marked ordinary phase')
 require(2*((chi[s]-chi[U])-(chi[s]-chi[V]))%8==0,'actual black-pair multiplicity annihilates this coordinate in aggregate current')
base=[z for z in r['routes']if z['name'].startswith('no_buffer')];compensated=[z for z in r['routes']if not z['name'].startswith('no_buffer')]
require(len(base)==3 and len(compensated)==5,'three primitive endpoint pairs and five advertised compensated routes')
for z in base:
 require(len(z['accepted_states'])==2,'original complete pinned fibre primitive')
 require((z['advertised_path'][0]^z['advertised_path'][-1])==(1<<len(z['sources']))-1,'all original binary slots change endpoint choice')
 require(z['resource_reentry_count']==0,'single primitive event is temporally nonreturning')
for z in compensated:
 require(z['resource_reentry_count']>0 and not z['satisfies_strict_temporal_nonreturn'],'no positive compensated route is mislabeled strict temporal nonreturn')
report={'status':'PASS','checks':checks,'input_sha256':{'contract':sha(a.contract),'edges':sha(e),'foundation':sha(a.foundation),'routes':sha(a.routes)},'ordinary_phase_obstruction':{'target_U':U,'target_V':V,'modulus':8,'potential_nonzero':[[U,4]],'U_nonloop_incidences':incidences,'U_nonloop_weighted_degree':degree,'odd_U_incidences':odd,'full_physical_rows_checked':N,'U_port_phase':4,'V_port_phase':0},'strict_temporal_nonreturn_proof':'Every restored binary buffer that changes must reacquire its initial resource, so it must instead remain fixed. Every original binary slot can change at most once. An intermediate therefore projects to a conformal original state; each authenticated original fibre has only its two endpoints. The first nontrivial move must be the entire primitive16/18 atom.','scope':'Independent mathematical review and exact sparse-character test of the separately hashed interpretation contract. The contradiction applies only to demanding every ordinary individual-port phase, not to full nodewise dyadic current. The no-return obstruction applies to the exact binary menus, fixed margins and restored buffers.'}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'U_degree':degree,'U_incidences':incidences,'contract_sha256':sha(a.contract)},sort_keys=True))
