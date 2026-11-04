"""Complementary physical/scoping audit; full kernel banks audited separately."""
from pathlib import Path
from itertools import permutations,combinations
import hashlib,json,math,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
P=ROOT.parent
NAMES=['prime_five_all_profiles','prime_five_all_ambient']
def need(x,m):
 if not x:raise RuntimeError(m)
def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
CAN={'S':tt(1),'P':tt(2),'E':tt(4),'O':tt(3),'C':[2,4,1]}
def op(Q):return [((1<<len(Q))-1)^r^(1<<i)for i,r in enumerate(Q)]
def blow(Q,blocks):
 off=[];n=0
 for B in blocks:off.append(n);n+=len(B)
 out=[0]*n
 for i,B in enumerate(blocks):
  for u in range(len(B)):
   for j,C in enumerate(blocks):
    for v in range(len(C)):
     if(i==j and B[u]>>v&1)or(i!=j and Q[i]>>j&1):out[off[i]+u]|=1<<(off[j]+v)
 return out,off
def masks(q,blocks,off):return sum(1<<(off[i]+j)for i,B in enumerate(blocks)for j in range(len(B))if q>>i&1)
def aug(R,a,c):
 n=len(R);out=list(R)+[1<<(n+1),0]
 for i in range(n):
  for j,M in [(n,a),(n+1,c)]:
   if M>>i&1:out[i]|=1<<j
   else:out[j]|=1<<i
 return out
def profiles(B):
 n=len(B);q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(1<<i for i,(u,v)in enumerate(zip(w,w[1:]))if not B[u]>>v&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad==bad:q[cuts.bit_count()+1]+=1
 return q
def rows(n,code):
 Q=[0]*n
 for k,(i,j)in enumerate(combinations(range(n),2)):
  if code>>k&1:Q[i]|=1<<j
  else:Q[j]|=1<<i
 return Q
def prime(Q):
 n=len(Q)
 return not any(M.bit_count()>1 and all(Q[v]&M in(0,M)for v in range(n)if not M>>v&1)for M in range(1,(1<<n)-1))
def rename(Q,p):
 R=[0]*len(Q)
 for i,j in permutations(range(len(Q)),2):
  if Q[i]>>j&1:R[p[i]]|=1<<p[j]
 return tuple(R)
def main():
 local=json.loads((P/NAMES[0]/'normal.json').read_text());ambient=json.loads((P/NAMES[1]/'normal.json').read_text())
 for name,data in zip(NAMES,[local,ambient]):
  need((P/name/'normal.json').read_bytes()==(P/name/'optimized.json').read_bytes(),'paired producer bytes')
  need(all(hashlib.sha256((P/name/f).read_bytes()).hexdigest()==h for f,h in data['source_pins'].items()),'frozen source pins')
 reps={r['code']:r['rows']for r in local['records']};orbits={code:{rename(Q,p)for p in permutations(range(5))}for code,Q in reps.items()};allprime={tuple(rows(5,c))for c in range(1024)if prime(rows(5,c))};need(len(allprime)==264 and set.union(*orbits.values())==allprime,'complete physical quotient relabelings');need({c:len(o)for c,o in orbits.items()}=={10:120,41:120,76:24},'three disjoint prime orbits')
 qp={name:profiles(B)for name,B in CAN.items()};normalized={name:[(z*pow(q[1],-1,8))%8 for z in(q+[0]*4)[:4]]for name,q in qp.items()};need(normalized=={'S':[0,1,0,0],'P':[0,1,2,0],'E':[0,1,6,4],'O':[0,1,6,6],'C':[0,1,2,2]},'all actual canonical profiles and odd scalars')
 graphs=[];ctrl=[]
 def add(R):graphs.append(R);return len(graphs)-1
 for r in local['literal_physical_controls']:
  blocks=[CAN[t]for t in r['types']];Q=reps[r['quotient_code']];Q=Q if r['side']==0 else op(Q);R,off=blow(Q,blocks);a,c=[masks(q,blocks,off)for q in r['block_query']];unit=math.prod(qp[t][1]for t in r['types'])%8
  need(R==r['old_rows']and[a,c]==r['physical_query'],'local physical labels and masks');need(unit==r['odd_normalization_unit'],'literal scalar')
  bank=next(v for rr in local['records']if rr['code']==r['quotient_code']for v in rr['full_run_bank']if v[:2]==r['block_query']);channel=bank[2+r['side']];weights=[1]+[normalized[t][2]for t in r['types']]+[2*int(t in 'OC')for t in r['types']]+[4*int(r['types'][i]!='S')*int(r['types'][j]!='S')for i,j in combinations(range(5),2)]
  need(unit*sum(x*y for x,y in zip(channel,weights))%8==r['predicted_mod8'],'all21-channel scalar contraction')
  ctrl.append({'kind':'local','source_record':r,'id':add(aug(R,a,c))})
 for r in ambient['physical_ambient_controls']:
  blocks=[CAN[t]for t in r['types']];Q=reps[r['quotient_code']];Q=Q if r['side']==0 else op(Q);core,off=blow(Q,blocks);a,c=[masks(q,blocks,off)for q in r['block_query']];outer,z=blow([2,4,1],[core,[0],[0]]);incoming=1<<z[2];unit=math.prod(qp[t][1]for t in r['types'])%8
  need(outer==r['old_rows']and[a|incoming,c|incoming]==r['query'],'ambient names and both root-incidence masks');need(r['outer_rows']==[2,4,1]and unit==r['odd_normalization_unit'],'actual outer quotient and scalar')
  ctrl.append({'kind':'ambient','source_record':r,'id':add(aug(outer,a|incoming,c|incoming)),'core_id':add(aug(core,a,c))})
 singleton=[]
 for rr in local['records']:
  Q=rr['rows']
  for M,(a,c,gap)in rr['singleton_pattern_mod4']:
   blocks=[CAN['S'if M>>i&1 else'P']for i in range(5)];ids=[]
   for QQ in [Q,op(Q)]:R,off=blow(QQ,blocks);ids.append(add(aug(R,masks(a,blocks,off),masks(c,blocks,off))))
   singleton.append({'quotient_code':rr['code'],'singleton_mask':M,'block_query':[a,c],'gap_mod4':gap,'ids':ids})
 with tempfile.TemporaryDirectory()as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);data=str(len(graphs))+'\n'+''.join(str(len(R))+' '+' '.join(map(str,R))+'\n'for R in graphs);values=[list(map(int,line.split()))for line in subprocess.check_output([str(exe)],input=data.encode()).decode().splitlines()]
 need(len(values)==len(graphs),'all independent full integer profiles')
 for x in ctrl:
  v=values[x.pop('id')];r=x.pop('source_record');expect=r['actual_mod16']if x['kind']=='local'else r['actual_ambient_mod16'];need(v[0]%16==expect,'literal exact H versus producer residue');x.update({k:r[k]for k in ['quotient_code','types','side','block_query','odd_normalization_unit']});x['full_physical_profile']=v
  if x['kind']=='ambient':
   core=values[x.pop('core_id')];need([z%8 for z in core]==r['core_first_three_mod8'],'literal core normalization for all three coordinates');need(v[0]==3*core[0]+core[1],'exact strong-ambient identity');need(v[0]%8==r['predicted_ambient_mod8'],'actual ambient residue');x['full_core_profile']=core
  else:need(v[0]%8==r['predicted_mod8'],'actual local normalized prediction')
 for x in singleton:
  v=[values[i]for i in x.pop('ids')];need((v[0][0]-v[1][0])%4==x['gap_mod4']==2,'every actual singleton-position separator');need([(a-b)%4 for a,b in zip(*v)]==[2,0,0],'complete singleton mod4 profile lift');x['full_profiles']=v
 need(len(ctrl)==22 and len(singleton)==93,'complete physical scopes')
 return {'status':'PASS','source_pins':{name:data['source_pins']for name,data in zip(NAMES,[local,ambient])},'auditor_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in ['audit.py','count.cpp']},'prime_quotient_count':len(allprime),'orbit_sizes':{str(k):len(v)for k,v in orbits.items()},'literal_canonical_profiles':qp,'normalized_canonical_profiles':normalized,'all22_source_physical_controls':ctrl,'all93_singleton_pattern_controls':singleton,'physical_profile_occurrences':len(graphs),'distinct_physical_graphs':len({tuple(R)for R in graphs}),'scope':'Complementary scope/normalization/physical-lift audit; full finite run banks are independently reviewed in separate packets.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
