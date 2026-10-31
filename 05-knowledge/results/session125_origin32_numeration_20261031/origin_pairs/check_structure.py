#!/usr/bin/env python3
"""Exact structural consequences and hostile controls; standard library only."""
import argparse,itertools,json,math,pathlib
from check_independent import read_reference,verify_path

def require(claim,message):
 if not claim:raise RuntimeError(message)
def square(x):return math.isqrt(x)**2==x
def cycle_edges(p):return {tuple(sorted(e)) for e in zip(p,p[1:]+p[:1])}
def insertions(v):
 ns=[a for a in range(1,v) if square(v+a)]
 return [(a,b) for a,b in itertools.combinations(ns,2) if square(a+b)]
def projections(masks,labels):return sorted({sum(((m>>(v-1))&1)<<j for j,v in enumerate(labels)) for m in masks})
def rejected(fn):
 try:fn()
 except (ValueError,RuntimeError):return True
 return False

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',default='structure_certificate.json');args=parser.parse_args()
 cert=json.loads(read_reference())
 D={r['n']:r for r in cert['path_orders']};C={r['n']:r for r in cert['cycles']}
 require(C[31]['component_counts']=={'14+17':1,'11+20':1},'G31 obstruction')
 require(len(C[32]['cycles'])==len(C[33]['cycles'])==1,'unique cycles')
 e32=cycle_edges(C[32]['cycles'][0]);e33=cycle_edges(C[33]['cycles'][0])
 require(e32-e33=={(5,31),(16,20)},'deleted edges')
 require(e33-e32=={(5,20),(16,33),(31,33)},'added edges')
 require(all(not insertions(v) for v in range(33,44)),'no insertion 33..43')
 require(insertions(44)==[(5,20)],'first arithmetic opportunity')
 c43=[1,15,21,43,38,26,10,6,19,30,34,2,7,18,31,5,20,29,35,14,11,25,39,42,22,3,33,16,9,27,37,12,4,32,17,8,28,36,13,23,41,40,24]
 k=c43.index(5)+1;c44=c43[:k]+[44]+c43[k:]
 for n,p in ((43,c43),(44,c44)):
  require(sorted(p)==list(range(1,n+1)),'first insertable permutation')
  require(all(square(a+b) for a,b in cycle_edges(p)),'first insertable cycle edges')
 require(cycle_edges(c43)-cycle_edges(c44)=={(5,20)},'pure insertion deletion')
 require(cycle_edges(c44)-cycle_edges(c43)=={(5,44),(20,44)},'pure insertion additions')
 # Independent root parametrization for all insertion triangles through largest vertex100.
 param=[]
 for r in range(1,15):
  for s in range(r+1,15):
   for t in range(s+1,15):
    if t*t>=r*r+s*s or (r*r+s*s+t*t)%2:continue
    a=(r*r+s*s-t*t)//2;b=(r*r+t*t-s*s)//2;v=(s*s+t*t-r*r)//2
    if v<=100:param.append((v,a,b))
 brute=[(v,a,b) for v in range(1,101) for a,b in insertions(v)]
 require(sorted(param)==brute,'triangle parametrization')
 # Nice-pair first onset and complete normalized-pair multiplicities.
 require(all(x['common_mask_count']==0 for x in cert['nice_pairs'][:-1]),'early compatible mask')
 onset=cert['nice_pairs'][-1]
 require((onset['n'],onset['common_mask_count'],onset['path_pair_count'])==(41,3,6),'first pair count')
 labels=(6,7,30,33)
 low=projections(D[40]['masks'],labels);high=projections(D[41]['masks'],labels)
 require(low==[1,2,7,8,13] and high==[3,4,6,9,12],'four-label spectra')
 require(set(low).isdisjoint(high),'four-label separation')
 # The last obstruction genuinely requires at least four coordinates in a pure projection certificate.
 no_three=True
 for labels3 in itertools.combinations(range(1,41),3):
  if set(projections(D[40]['masks'],labels3)).isdisjoint(projections(D[41]['masks'],labels3)):no_three=False;break
 require(no_three,'a smaller coordinate separator exists')
 first=onset['witnesses'][0];p=first['shorter'];q=first['longer']
 require(verify_path(p,41)==verify_path(q,42),'shared positional parity')
 require(q.index(42)+1==12,'new vertex need not be terminal')
 require(p.index(2)%2==0,'integer parity differs from position parity')
 neg={}
 bad=p[:];bad[1]=bad[0];neg['duplicate_label']=rejected(lambda:verify_path(bad,41))
 bad2=p[:];bad2[1],bad2[2]=bad2[2],bad2[1];neg['nonsquare_edge']=rejected(lambda:verify_path(bad2,41))
 neg['false_odd_mask']=rejected(lambda:verify_path(p,41,first['odd_mask_common']^2))
 neg['wrong_endpoint']=rejected(lambda:verify_path(p[::-1],41))
 # Valid paths are insufficient: select the first exact G40 and G41 witnesses.
 r=D[40]['witnesses'][str(D[40]['masks'][0])]['path'];s=D[41]['witnesses'][str(D[41]['masks'][0])]['path']
 neg['incompatible_valid_paths']=verify_path(r,40)!=(verify_path(s,41)&((1<<40)-1))
 require(all(neg.values()),'hostile control failure')
 out={'status':'PASS','first_nice_pair_n':41,'nice_pair_count':6,'parity_mask_count':3,'last_obstruction_labels':labels,'G40_projection':low,'G41_projection':high,'no_separator_on_three_labels':no_three,'minimum_C32_C33_old_edge_deletions':2,'insertion_free_arrivals':list(range(33,44)),'first_arithmetic_insertion_after_origin':{'new_vertex':44,'old_edge':[5,20]},'first_pure_insertion_cycles':{'C43':c43,'C44':c44},'insertion_triangles_through_100':brute,'hostile_controls':neg}
 pathlib.Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
