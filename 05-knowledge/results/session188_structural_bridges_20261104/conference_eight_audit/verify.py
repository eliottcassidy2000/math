"""Independent degree-constrained conference enumeration and full bank audit."""
from pathlib import Path
from itertools import combinations,permutations
from collections import Counter
import struct,json,hashlib,tempfile,subprocess
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'conference_eight';BANK=ROOT.parent/'eight_pair_quotient'/'blind8.bin'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def encode(a):return sum(((a[i]>>j)&1)<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))
def rows(n,c):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if c>>k&1 else(j,i);a[u]|=1<<v
 return a
def regulars():
 n=7;a=[0]*n;wins=[0]*n
 def walk(i):
  if i==n:
   if all(x==3 for x in wins):yield a[:]
   return
  need_wins=3-wins[i];rem=list(range(i+1,n))
  if not 0<=need_wins<=len(rem):return
  for chosen in combinations(rem,need_wins):
   chosen=set(chosen);changes=[];wins[i]+=need_wins
   for j in rem:
    if j in chosen:a[i]|=1<<j;changes.append((i,j))
    else:a[j]|=1<<i;wins[j]+=1;changes.append((j,i))
   if all(wins[j]<=3<=wins[j]+n-1-i for j in rem):yield from walk(i+1)
   for u,v in changes:a[u]^=1<<v
   for j in rem:
    if j not in chosen:wins[j]-=1
   wins[i]-=need_wins
 yield from walk(0)
def switched(a,mask):
 n=len(a);b=[0]*n
 for i,j in combinations(range(n),2):
  forward=((a[i]>>j)&1)^((mask>>i)&1)^((mask>>j)&1);u,v=(i,j) if forward else(j,i);b[u]|=1<<v
 return b
def conference(a):
 full=(1<<len(a))-1
 return all((((a[i]^a[j])&(full^(1<<i)^(1<<j))).bit_count()==3) for i,j in combinations(range(8),2))
def isprime(a):
 n=len(a);full=(1<<n)-1
 return all(any(a[v]&M not in(0,M) for v in range(n) if not(M>>v&1)) for M in range(1,full) if M.bit_count()>1)
def main():
 regular=0;bases=[]
 for a in regulars():
  regular+=1
  if all((a[i]&a[j]).bit_count()==1 for i,j in combinations(range(7),2)):bases.append(encode(a))
 bases.sort();need(regular==2640 and len(bases)==240,'independent constrained enumeration counts')
 src=json.loads((SRC/'census.json').read_text());need(bases==src['normalized_seven_codes'],'complete normalized code list')
 raw=BANK.read_bytes();need(len(raw)==30720*4,'exact bank size');bank=list(struct.unpack('<30720I',raw));need(bank==sorted(set(bank)),'bank distinct sorted')
 generated=set();normalizations=Counter()
 for code in bank:
  a=rows(8,code);need(conference(a),'all bank Gram identities');mask=sum(1<<j for j in range(1,8) if not(a[0]>>j&1));b=switched(a,mask);need(b[0]==254,'root normalization');base=encode([r>>1 for r in b[1:]]);normalizations[base]+=1
 need(set(normalizations)==set(bases) and set(normalizations.values())=={128},'every unique normalized class represented')
 for c in bases:
  a=rows(8,(c<<7)|127)
  for mask in range(128):generated.add(encode(switched(a,mask<<1)))
 need(sorted(generated)==bank,'complete independently generated conference set equals h-zero bank')
 orbit_records=[];allorbits=[];covered=set();reps=[600626,600653,602967,1472106]
 for rep,size in zip(reps,[1920,13440,13440,1920]):
  a=rows(8,rep);orbit=set()
  for p in permutations(range(8)):orbit.add(encode([sum(1<<j for j,v in enumerate(p) if a[u]>>v&1) for u in p]))
  need(len(orbit)==size and not(orbit&covered),'complete relabeling orbit');covered|=orbit;allorbits.append(orbit);orbit_records.append({'minimum_code':min(orbit),'size':len(orbit),'prime':isprime(a),'universal_vertex':any(r.bit_count() in(0,7) for r in a),'converse_code':rep^((1<<28)-1)})
 need(all((reps[i]^((1<<28)-1)) in allorbits[j] for i,j in enumerate([3,2,1,0])),'converse orbit pairing');need(covered==generated,'orbit partition of whole class');need([r['prime'] for r in orbit_records]==[False,True,True,False],'prime orbit classification')
 # Independent subset start/end parity check, without t7 lookup tables.
 with tempfile.TemporaryDirectory() as t:
  exe=Path(t)/'endpoint';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'endpoint.cpp'),'-o',str(exe)],check=True);checked=int(subprocess.check_output([str(exe)],input=(str(len(bank))+'\n'+'\n'.join(map(str,bank))+'\n').encode()))
 need(checked==30720,'all conference shared cubic and endpoint checks')
 # At order4 only two-element even supports contribute to a shared cubic;
 # each anchor has two partners after excluding the other named vertex.
 need((2%2)==0,'TT4 shared cubic boundary');T=rows(4,63);dot=sum((1 if T[0]>>k&1 else -1)*(1 if T[1]>>k&1 else -1) for k in(2,3));need(dot==2,'TT4 nonconference Gram boundary')
 return {'status':'PASS','regular_degree_constrained_count':regular,'all_normalized_codes':bases,'conference_class_size':len(generated),'exact_complete_hzero_set_equality':True,'shared_cubic_zero_and_endpoint_ones_checks':checked,'normalization_multiplicity':128,'orbit_records':orbit_records,'TT4_boundary_inner_product':dot,'source_pins':{f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','census.cpp','census.json','normal.json','optimized.json']},'bank_sha256':hashlib.sha256(raw).hexdigest(),'scope':'Exact order-eight equivalence only; no general-order equivalence or response collision.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
