#!/usr/bin/env python3
"""Independent singleton/profile/family/highest-node synthesis replay."""
from pathlib import Path
from itertools import combinations,product
import hashlib,json,subprocess,tempfile,sys
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;SRC=BASE/'cyclic_ambient_synthesis'
sys.path.insert(0,str(BASE/'even_mod16_audit'))
# Reuse only independently written physical graph/profile helpers, never producer code.
import audit_ambient as I
need=I.need;op=I.op;blow=I.blow;profile=I.profile;modules=I.modules;reverse=I.reverse;augment=I.augment;endpoints=I.endpoints;build=I.build
PINS={'PROOF.md':'396a8876d98a7900a7b76797ea369d211583b744c43f01d80b498f7a09829f3e','verify.py':'b3c341501408c3946b6113fbac94e1654f6a8d7cea8a85a841f6be3f705c2e92','physical.cpp':'4d980994dfafa599ddb065768ce34fdf54abacc936a9fb89e177834e38b3152c','normal.json':'e6116adaf3f31444245f24abfcb4107b28af59ee4c0983b0aa41038c9bb84d16','optimized.json':'e6116adaf3f31444245f24abfcb4107b28af59ee4c0983b0aa41038c9bb84d16'}
def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen source '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source status')
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((BASE/f).read_bytes()).hexdigest()==h,'dependency '+f)
 data=json.loads((BASE/'five_profile_types'/'normal.json').read_text());table={tuple(r['types']):r for r in data['canonical_records']};single=[r for r in data['canonical_records'] if 'S' in r['types']];need(len(single)==61==d['singleton_canonical_triples'],'complete singleton domain');canon={'S':[0],'P':[2,0],'E':[14,12,8,0],'O':[6,4,0],'C':[2,4,1]};keys=list(product(range(8),repeat=2));jobs={};picked={}
 def retain(T,a,c,values,n=None):
  G=augment(T,a,c)
  if n is not None:
   M=((1<<n)-1)|(1<<len(T))|(1<<(len(T)+1));need(all(G[v]&M in(0,M) for v in range(n,len(T))),'same actual ambient core module')
  k=(G,len(values));need(k not in jobs or jobs[k]==values,'exact repeated graph consistency');jobs[k]=values
 need([r['types'] for r in d['singleton_records']]==[r['types'] for r in single],'entire singleton list')
 for r in d['singleton_records']:
  label=tuple(r['types']);c0=table[label];res=[v%4 for v in c0['contrasts_mod8']];need(res==r['complete_mod4_bank'],'full mod4 canonical row');index=next(i for i,v in enumerate(res) if v);query=keys[index];need(res[index]==2 and list(query)==r['block_query'],'first separator');T,U,a,c,n=build([canon[t] for t in label],[0],query);need(T==r['source_old_rows'] and U==r['target_old_rows'] and [a,c]==r['physical_query'],'canonical physical labels');vals=[r['source_first_three'],r['target_first_three']];need([v[0] for v in vals]==r['frozen_H_values']==c0['values'][index] and [(a-b)%4 for a,b in zip(*vals)]==[2,0,0],'all first-three cover residues');picked[label]=query
  for V,v in zip((T,U),vals):retain(V,a,c,v)
 for r in d['larger_child_ambient_controls']:
  B=r['child_rows'];label=tuple(r['types']);query=picked[label];need(list(query)==r['block_query'],'type-selected query')
  for actual,t in zip(B,label):
   q=profile(tuple(actual));p=profile(tuple(canon[t]));unit=q[1]*pow(p[1],-1,8)%8;need(all(q[k]%8==unit*(p[k] if k<len(p) else 0)%8 for k in range(len(q))),'complete larger-child profile type')
  T,U,a,c,n=build(B,r['outer_rows'],query);need(T==r['source_old_rows'] and U==r['target_old_rows'] and [a,c]==r['physical_query'],'larger fixed-ambient graph');need((r['values'][0]-r['values'][1])%4==2,'actual singleton ambient residue')
  for V,v in zip((T,U),r['values']):retain(V,a,c,[v],n)
 B=[[0],[6,4,0],[6,4,0]];mm=[modules(b) for b in B];selections=list(product(*[[[M for j,M in enumerate(z) if bits>>j&1] for bits in range(1<<len(z))] for z in mm]));z=picked['S','O','O'];dual=(7^z[1],7^z[0]);need([list(z),list(dual)]==d['root_family_block_queries'],'only core duality');strict=0
 recs=d['root_family_records'];need(len(recs)==128==d['root_family_contexts'],'complete family count')
 for index,r in enumerate(recs):
  K=[0] if index<64 else[2,4,1];selection=selections[index%64];need(r['outer_rows']==K and tuple(r['selected_original_child_modules'])==selection,'all original-module selections');target=[op(reverse(b,sel)) for b,sel in zip(B,selection)];need(target==r['target_child_rows'],'original XOR plus root');qp=[profile(tuple(b)) for b in target];old=[profile(tuple(b)) for b in B];need(qp==r['target_profiles'] and all(y>=x for a,b in zip(old,qp) for x,y in zip(a,b)),'complete profile monotonicity');need((qp!=old)==r['profile_strict'],'strict flag');strict+=qp!=old
  base=r['baseline_values'];vals=r['target_values'];diff=[s-t for s,t in base];need(diff[0]==-diff[1] and diff[0]%4==2,'opposite ambient baseline signs');choice=0 if diff[0]<0 else 1;need(choice==r['negative_query_index'] and vals[choice][0]-vals[choice][1]<=diff[choice]<0,'one negative member for every target')
  for query,bv,tv in zip((z,dual),base,vals):
   for tgt,vv in((None,bv),(target,tv)):
    T,U,a,c,n=build(B,K,query,tgt)
    for V,v in zip((T,U),vv):retain(V,a,c,[v],n)
 need(strict==d['profile_strict_family_contexts'],'strict family census')
 h=d['multiple_highest_control'];M,_=blow([2,4,1],[[2,0],[2,4,1],[0]]);V,_=blow([4,1,2],[[2,0],[2,4,1],[0]]);T,_=blow([2,4,1],[M,[2,4,1],[0]]);B0,_=blow([2,4,1],[V,[2,4,1],[0]]);U=reverse(T,[63,448]);need([T,B0,U]==[h['source_rows'],h['baseline_rows'],h['target_rows']],'two-highest physical targets');allmods=[1<<i for i in range(len(T))]+modules(T);strong=[m for m in allmods if not any((m&x) and(m&~x) and(x&~m) for x in allmods)];nodes=[];changed=[]
 for m in sorted(strong):
  children=[x for x in sorted(strong) if x!=m and x&m==x and not any(x!=y and y!=m and x&y==x and y&m==y for y in strong)]
  if len(children)<2:continue
  reps=[(x&-x).bit_length()-1 for x in children];flag=any(bool(T[i]>>j&1)!=bool(U[i]>>j&1) for i in reps for j in reps if i!=j);nodes.append(dict(mask=m,children=children,changed_quotient=flag))
  if flag:changed.append(m)
 highest=[m for m in changed if not any(m!=x and m&x==m for x in changed)];need(nodes==h['original_strong_nodes'] and highest==h['highest_changed_masks'] and sorted(highest)==[63,448],'genuine highest strong nodes');need(sum(endpoints(T)[-1])==sum(endpoints(U)[-1]),'scalar equality');query=picked['P','C','S'];a,c=[sum(((1<<sz)-1)<<o for i,(sz,o) in enumerate(zip([2,3,1],[0,2,5])) if mask>>i&1)|(1<<9) for mask in query];need([a,c]==h['physical_query'] and h['values'][1]==h['values'][2] and (h['values'][0]-h['values'][1])%4==2,'exact other-highest restoration')
 for V,v in zip((T,B0,U),h['values']):retain(V,a,c,[v],6)
 primes=[]
 for code in range(64):
  Q=[0]*4
  for k,(i,j) in enumerate(combinations(range(4),2)):
   v,w=(i,j) if code>>k&1 else(j,i);Q[v]|=1<<w
  if not any(M.bit_count()<4 for M in modules(Q)):primes.append(code)
 need(primes==d['prime_four_codes']==[],'no four-name prime quotient')
 with tempfile.TemporaryDirectory(prefix='cyclic-synthesis-audit-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);items=list(jobs);wire=str(len(items))+'\n'+'\n'.join(str(len(G))+' '+str(K)+' '+' '.join(map(str,G)) for G,K in items)+'\n';actual=[list(map(int,line.split())) for line in subprocess.check_output([str(exe)],input=wire.encode()).decode().splitlines()];need(actual==[jobs[k] for k in items],'all independent physical profile values')
 helpers=['audit_ambient.py','cover_automaton.py','cover_physical.py','automaton.py']
 return dict(status='PASS',source_sha256=PINS,dependencies=d['dependency_pins'],independent_helper_sha256={f:hashlib.sha256((BASE/'even_mod16_audit'/f).read_bytes()).hexdigest() for f in helpers},singleton_types=61,root_families=128,strict_families=strict,highest_changed=highest,actual_unique_profile_graphs=len(jobs),complete_actual_values=[[list(G),K,v] for (G,K),v in zip(items,actual)])
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
