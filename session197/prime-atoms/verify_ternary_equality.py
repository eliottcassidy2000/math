#!/usr/bin/env python3
from verify_minimum_ternary import hp,rank,clone,strong,canonical,require
from pathlib import Path
from collections import Counter
import hashlib,itertools,json,os,re,subprocess,tempfile,math
P=Path(__file__).resolve().parent
EXPECTED=[([2,20,41,35,45,195,31,95],3),([4,5,56,67,43,75,151,63],4)]

def main():
 with tempfile.TemporaryDirectory(prefix='ternary-equality-',dir=P) as tmp:
  exe=Path(tmp)/'search';bankfile=Path(tmp)/'bank.jsonl';mode=os.environ.get('TERNARY_CPP_MODE','optimized')
  flags=['-O3'] if mode=='optimized' else ['-O0','-D_GLIBCXX_ASSERTIONS']
  subprocess.run([os.environ.get('CXX','g++'),'-std=c++17',*flags,str(P/'verify_ternary_equality.cpp'),'-o',str(exe)],check=True)
  result=subprocess.run([str(exe),str(bankfile)],text=True,capture_output=True,check=True)
  require(result.stdout.splitlines()[-1]=='DONE tested=152 hits=2','complete equality search outcome')
  entries=[json.loads(line) for line in bankfile.read_text().splitlines()];counts={3:1,4:1,5:6,6:35,7:352,8:736,9:380,10:84}
  require(dict(Counter(e['n'] for e in entries))==counts,'complete equality bank counts')
  checked=first=0;hits=[];seen=set()
  for e in entries:
   n=e['n'];rows=e['rows'];H=e['H'];key=(n,tuple(rows))
   require(key not in seen and canonical(rows)==tuple(rows),'exact distinct canonical cores');seen.add(key)
   require(strong(rows) and hp(rows)==H<=179,'independent core counts')
   if H!=179:continue
   for r in range(n):
    if rank(rows,r)!=2:continue
    checked+=1;H2=hp(clone(rows,r,2))
    if H2!=539:continue
    first+=1;H3=hp(clone(rows,r,3))
    if H3==1619:hits.append((rows,r))
  require(checked==152 and hits==EXPECTED,'complete exact equality cases')
  A,a=hits[0];B,b=hits[1];p=list(reversed(range(8)))
  require(p[a]==b and all(((A[i]>>j)&1)==((B[p[j]]>>p[i])&1) for i in range(8) for j in range(8)),'rooted converse certificate')
  require(A!=B,'distinct unrooted isomorphism classes')
  original=[2,12,145,149,99,143,175,19];q=[0,1,2,4,5,6,7,3]
  require(q[7]==a and all(((original[i]>>j)&1)==((A[q[i]]>>q[j])&1) for i in range(8) for j in range(8)),'original attaining graph certificate')
  auto=[]
  for rows,r in hits:
   ds=[v.bit_count() for v in rows];G=[]
   for g in itertools.permutations(range(8)):
    if any(ds[i]!=ds[g[i]] for i in range(8)):continue
    if all(((rows[i]>>j)&1)==((rows[g[i]]>>g[j])&1) for i in range(8) for j in range(i+1,8)):G.append(g)
   require(G==[tuple(range(8))],'trivial full automorphism group');auto.append(len(G))
  require(5**10>179**3,'all-order equality core bound')
  generation=[]
  for line in result.stderr.splitlines():
   m=re.fullmatch(r'order (\d+) candidates (\d+) accepted (\d+) classes (\d+)',line)
   if m:generation.append(dict(zip(['order','extension_candidates','accepted_before_isomorphism','classes'],map(int,m.groups()))))
  receipt={'status':'PASS','strong_H_cap':179,'retained_classes_by_order':counts,'generation':generation,'retained_bank_sha256':hashlib.sha256(bankfile.read_bytes()).hexdigest(),'independent_graph_checks':len(entries),'candidate_roots':checked,'first_equation_survivors':first,'complete_rooted_strong_cores':[{'rows':rows,'root':r,'profile_runs':[179,180,60]} for rows,r in hits],'converse_map':p,'attaining_graph_isomorphism':q,'automorphism_orders':auto,'source_class_counts':[{'n':n,'unlabelled':2*(n-7),'labelled':2*(n-7)*math.factorial(n)} for n in range(8,13)],'scope':'All amplitude60 fixed-root transitive-clone sources, with arbitrary neutral ordered singleton padding.'}
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
