"""Exact prime endpoint-parity obstruction and complete physical response controls."""
from pathlib import Path
from itertools import permutations,product
import hashlib,importlib.util,json,random,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def subsets(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M

def main():
 Q=[238,4,240,118,195,146,34,74];n=8;need(all(not(Q[i]>>i&1)for i in range(n))and all(bool(Q[i]>>j&1)!=bool(Q[j]>>i&1)for i in range(n)for j in range(i+1,n)),'valid physical tournament rows');P=u.paths(Q);E=u.paths(u.op(Q));Svec=P[-1];Evec=E[-1];need(Svec==[87,13,13,42,31,19,18,26]and Evec==[13,65,29,22,13,29,52,26],'exact endpoint vectors');need(sum(Svec)==sum(Evec)==249 and all((a+b)%2==0 for a,b in zip(Svec,Evec)),'balanced endpoint parity')
 paths=[w for w in permutations(range(n))if all(Q[a]>>b&1 for a,b in zip(w,w[1:]))];sv=[sum(w[0]==i for w in paths)for i in range(n)];ev=[sum(w[-1]==i for w in paths)for i in range(n)];need(sv==Svec and ev==Evec,'independent complete permutation enumeration');mods=[M for M in range(1,255)if M.bit_count()>1 and all(Q[v]&M in(0,M)for v in range(n)if not M>>v&1)];need(not mods,'no proper nontrivial module')
 h=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(n):
   if i==j:continue
   h[i][j]=sum(P[M][i]+E[M][i]for M in range(1,256)if M>>i&1 and not M>>j&1 and M.bit_count()%2==0)%2
 i,j=next((i,j)for i in range(n)for j in range(n)if h[i][j]);blocks=[u.tt(2)]*8;T,o=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);cq=[]
 for bits in range(8):
  a=(3<<(2*i)if bits&1 else 0)|(3<<(2*j)if bits&2 else 0);c=3<<(2*j)if bits&4 else 0;cq.append((a,c))
 half=[m for m in range(256)if m.bit_count()<=2];bank=list(product(half,repeat=2));rng=random.Random(188808);extra=[(rng.randrange(256),rng.randrange(256))for _ in range(32)];cases=[(Q,u.op(Q),bank+extra),(T,U,cq)];lines=[str(len(cases))]
 for A,B,qs in cases:lines += [f'{len(A)} {len(qs)}',' '.join(map(str,A)),' '.join(map(str,B))]+[f'{a} {c}'for a,c in qs]
 with tempfile.TemporaryDirectory()as td:
  lit=Path(td)/'literal';ce=Path(td)/'census';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(lit)],check=True);subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'census.cpp'),'-o',str(ce)],check=True);raw=list(map(int,subprocess.check_output([str(lit)],input=('\n'.join(lines)+'\n').encode()).split()));census=[list(map(int,line.split()))for line in subprocess.check_output([str(ce)]).decode().splitlines()]
 need(all(row[2]==0 and sum(row[3:])==row[1]for row in census)and[row[:2]for row in census]==[[2,2],[4,64],[6,32768]],'complete smaller-even census');values=[raw[2*k:2*k+2]for k in range(len(bank)+len(extra))];pv=[raw[2*(len(bank)+len(extra)+k):2*(len(bank)+len(extra)+k)+2]for k in range(8)];need(len(raw)==2*(len(bank)+len(extra)+8),'all physical outputs');coef=sum((-1)**(3-b.bit_count())*(a-c)for b,(a,c)in enumerate(pv));need(coef%8==4,'actual shared-cubic pair-inflation separator')
 vdict=dict(zip(bank,values));coeff={};records=[];hist={k:0 for k in range(5)}
 for a,c in bank:
  z=[0,0]
  for x in subsets(a):
   for y in subsets(c):
    sign=(-1)**(a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count())
    for k,v in enumerate(vdict[x,y]):z[k]+=sign*v
  coeff[a,c]=z;records.append([a,c,*vdict[a,c],*z]);hist[a.bit_count()+c.bit_count()]+=z[0]!=z[1]
 full=[]
 for k,(a,c)in enumerate(extra):
  z=[sum(v[t]for(x,y),v in coeff.items()if x&a==x and y&c==y)for t in range(2)];need(z==values[len(bank)+k],'full named Boolean reconstruction');full.append([a,c,*z])
 first=next(row[:4]for row in records if row[2]!=row[3]);need(hist[4]==0,'whole-converse quartic symmetry retained as actual zero data')
 deps=['parity_profile_transfer/PROOF.md','eight_pair_quotient/PROOF.md','uniform_block_starts/verify.py']
 return{'status':'PASS','quotient_rows':Q,'Hamilton_start_counts':Svec,'Hamilton_end_counts':Evec,'Hamilton_paths':paths,'common_endpoint_parity_mask':55,'proper_nontrivial_modules':mods,'smaller_even_census_columns':['order','labeled_tournaments','balanced_count','histogram_by_parity_mask_weight...'],'smaller_even_census':census,'complete_named_bank':{'columns':['a','c','F_source','F_target','coefficient_source','coefficient_target'],'records':records,'difference_histogram':hist,'full_mask_reconstructions':full,'first_actual_separator':first},'pair_inflation':{'source_old_rows':T,'target_old_rows':U,'shared_cubic_matrix':h,'coefficient_indices':[i,j,j],'physical_queries':cq,'physical_values':pv,'exact_contrast_coefficient':coef,'contrast_mod8':coef%8},'physical_counts':len(raw),'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in['PROOF.md','verify.py','literal.cpp','census.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()for f in deps},'scope':'Smallest even-order endpoint-parity criterion failure, witnessed by a prime eight-name quotient; no exact or full named-response collision.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
