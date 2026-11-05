#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,tempfile,itertools as it,math
from functools import lru_cache
P=Path(__file__).resolve().parent
if hashlib.sha256((P/'verify_c3_complete.cpp').read_bytes()).hexdigest()!='e49f0b9ff6e29f7b7260d6516ccb0fafe03ea4bbcf06dccc0664b7fddead806f':raise RuntimeError('C3 dependency changed')
from verify_start_two import words,relabel,labelled,need

def coefficients(numerator,N=101):
 a=[0]*N
 for n in range(N):a[n]=numerator.get(n,0)+sum(c*a[n-j]for j,c in ((1,6),(2,-11),(3,6),(4,-1))if n>=j)
 return a

def main():
 sym=coefficients({5:1,6:-2});root=coefficients({5:2,6:-4,7:1});exception=[0]*101
 for n in range(101):exception[n]=(n==5)+(3*exception[n-1]if n else 0)-(exception[n-2]if n>1 else 0)
 counts={n:[0,0,0,0]for n in range(4,13)}
 with tempfile.TemporaryDirectory(prefix='c3-aut-')as d:
  d=Path(d);exe=d/'verifier';csv=d/'records.csv'
  subprocess.run(['g++','-O2','-std=c++17',str(P/'verify_c3_complete.cpp'),'-o',str(exe)],check=True)
  subprocess.run([str(exe),str(csv)],check=True,stdout=subprocess.DEVNULL)
  for line in csv.read_text().splitlines():
   n,c,*orig=map(int,line.split(','))
   if c!=3:continue
   counts[n][0]+=1;r=orig;roots=[u for u in range(n)if r[u].bit_count()==1]
   if not roots:
    r=[sum(1<<v for v in range(n)if orig[v]>>u&1)for u in range(n)];roots=[u for u in range(n)if r[u].bit_count()==1]
   need(roots,'missing fixed degree-one root');w=words(r,roots[0],4);need(len(w)==3,'cycle count at degree-one root')
   forms={relabel(r,p)for p in w};need(len(forms)in(1,3),'automorphism divisibility')
   if len(forms)==3:continue
   counts[n][1]+=1;p=[0]*n
   for u,v in zip(w[0],w[1]):p[u]=v
   need(p!=list(range(n))and all(p[p[p[u]]]==u for u in range(n)),'nontrivial order-three generator')
   need(all(bool(orig[u]>>v&1)==bool(orig[p[u]]>>p[v]&1)for u in range(n)for v in range(n)),'automorphism arcs')
   full=(1<<n)-1
   @lru_cache(None)
   def paths(mask,u):
    if mask==full:return 1
    total=0;todo=orig[u]&~mask
    while todo:
     b=todo&-todo;todo-=b;total+=paths(mask|b,b.bit_length()-1)
     if total>=4:return 4
    return total
   for u in range(n):
    if p[u]==u and paths(1<<u,u)==3:
     counts[n][2]+=1
     if not any(orig[v].bit_count()==n-2 for v in range(n)if v!=u):counts[n][3]+=1
 for n,(total,s,r,e)in counts.items():need((s,r,e)==(sym[n],root[n],exception[n]),('symmetric inventories',n,(s,r,e),(sym[n],root[n],exception[n])))
 literal=[]
 for n in range(1,7):
  count=0
  if n>=3:
   for r in labelled(n):
    C=sum(all(r[u]>>v&1 for u,v in zip((0,)+p,p+(0,)))for p in it.permutations(range(1,n)))
    count+=C==3
  expected=math.factorial(n)*(3*counts[n][0]-2*counts[n][1])//3 if n>=4 else 0
  need(count==expected,'complete labelled orbit mass');literal.append(count)
 result={'status':'PASS','inventory':[{'order':n,'C3_classes':a,'automorphism_order_three':b,'rooted_start_three_stabilizer_three':c,'exceptional_rooted_stabilizer_three':d,'labelled_C3':math.factorial(n)*(3*a-2*b)//3}for n,(a,b,c,d)in counts.items()],'complete_labelled_C3_counts_through6':literal,'nonrigid_generating_function':'z^5(1-2z)/(1-3z+z^2)^2'}
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'c3_automorphisms_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
