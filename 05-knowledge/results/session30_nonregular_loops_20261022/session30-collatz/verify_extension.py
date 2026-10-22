#!/usr/bin/env python3
"""Extra guarded relation families in the same rank-one contraction cone."""
import argparse,hashlib,json
from pathlib import Path
import verify as v

def extended(N):
 for tag,u,s,c in v.bank(N):yield tag,u,((s,c),)
 cut=(3*N+1)//2
 for tag,rs in [('C',(3,39,57)),('D',(23,41,59))]:
  for r in rs:
   for u in range(cut+(r-cut)%72,2*N,72):
    cs=(((32*u+9)//3,5),((16*u+3)//3,4)) if tag=='C' else (((128*u+41)//3,7),)
    yield tag,u,cs

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args();V=set(json.loads((args.repo_root/'05-knowledge/results/session18_native_annulus/frontier.json').read_text())['frontier'])
 rows=0;digest=hashlib.sha256();scales=[];routes=0
 for N in list(range(1,1001))+[4096,8192,32768,63666,65536,100000]:
  owned=set();counts=dict.fromkeys('ABCD',0)
  for tag,u,cs in extended(N):
   counts[tag]+=1;ids=[u]+[s for s,c in cs];v.ck(len(ids)==len(set(ids)) and not owned.intersection(ids),'source ownership across all four families');owned.update(ids)
   w=v.charge(N,u);h=v.Z
   for s,c in cs:
    v.ck(v.edge(s)==(u,c),'literal edge and exact clock');h=v.add(h,v.mul((-1)**c,v.charge(64*N,s)))
    v.ck(4*N<=s<256*N,'same original annulus')
   target=v.A if tag in 'AC' else v.B
   v.ck((w,h)==target,'same rank-one cone after signed synchronous child sum')
   z=v.add(w,h);v.ck(z==(0,sum(w)//3) and sum(w)%3==0,'integer-lattice contraction')
   if N>=7959:v.ck(all(s not in V for s,c in cs),'all children outside frontier')
   if N==63666:
    ru=v.route(u,V)
    v.ck(ru is not None,'finite parent arrival test only')
    for s,c in cs:v.ck(v.route(s,V)==(ru[0],ru[1]+c,ru[2]+1),'literal synchronized first-arrival clock and depth');routes+=1
   digest.update(json.dumps([N,tag,u,u%15,u%35,[[s,c,s%15,s%35] for s,c in cs]],separators=(',',':')).encode());rows+=1
  cut=(3*N+1)//2
  for tag,lo,hi,rs in [('A',cut,2*N,(9,27,63)),('B',N,cut,(17,35,71)),('C',cut,2*N,(3,39,57)),('D',cut,2*N,(23,41,59))]:v.ck(counts[tag]==v.count(lo,hi,rs),'exact floor count')
  if N%144==0:v.ck(set(counts.values())=={N//48},'equal aligned-scale density')
  if N>=4096:scales.append({'N':N,'families':counts,'distinct_original_sources':len(owned)})
 # Sharp operator identity 4Q(w)-9Q(Sw)=3(w1-w2)^2.
 for a in range(-30,31):
  for b in range(-30,31):
   v.ck(4*v.q((a,b))-(a+b)**2==3*(a-b)**2,'sharp 4/9 norm contraction')
 out={'checks':v.CHECKS,'literal_clusters':rows,'finite_child_arrival_checks':routes,'source_clock_register_hash':digest.hexdigest(),'scales':scales,'uniform153N':'OPEN','outside_bank_cross_terms':'retained by REPORT section4, with new bank partition','completed_contraction':'4/9, sharp over real span','frozen_base_verifier_sha256':hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest()};args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
