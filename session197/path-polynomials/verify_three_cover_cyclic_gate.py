#!/usr/bin/env python3
from verify import *
from math import factorial,isqrt
import hashlib,argparse
from pathlib import Path
K=[[5,6,0],[6,28,24],[0,24,72]]
M=[[15,-11,2],[-11,13,-4],[2,-4,2]]
def cov3(a):return (cover_polynomial(a)+[0,0,0])[1:4]
def form(v,w):return sum(v[i]*K[i][j]*w[j] for i in range(3) for j in range(3))
def phi(n):return 2*9**n-8*6**n+13*4**n+4*3**n-22*2**n+15
def converse(a):return tuple(row^((1<<len(a))-1)^(1<<i) for i,row in enumerate(a))
def strong(a):
 for v in range(len(a)):
  seen=1<<v
  while True:
   new=seen
   for u in range(len(a)):
    if seen>>u&1:new|=a[u]
   if new==seen:break
   seen=new
  if seen!=(1<<len(a))-1:return False
 return True

def inverse(target):
 n=1
 while phi(n)<=target:
  for code in range(1<<(n*(n-1)//2)):
   t=tournament(n,code);v=cov3(t)
   if form(v,v)==target:return {'input_rows':list(t),'cover_vector':v,'output_rows':list(cyclic(t,converse(t),transitive(2)))}
  n+=1
 return None

def run():
 got=[[0]*3 for _ in range(3)];ns=[]
 for k in (1,2):
  N=[[0]*3 for _ in range(3)]
  for i,j in product(range(1,4),repeat=2):
   for start in range(3):
    word=tuple((start+d)%3 for d in range(i+j+k))
    if tuple(word.count(c) for c in range(3))==(i,j,k):N[i-1][j-1]+=1
   got[i-1][j-1]+=factorial(i)*factorial(j)*factorial(k)*N[i-1][j-1]
  ns.append(N)
 require(got==K and ns==[[[3,1,0],[1,1,0],[0,0,0]],[[1,1,0],[1,3,1],[0,1,1]]],'literal cyclic-run kernel')
 inputs=[tournament(n,code)for n in range(1,4)for code in range(1<<(n*(n-1)//2))];physical=0;profiles=0
 for t,u in product(inputs,repeat=2):
  v=cov3(t);w=cov3(u);a=cyclic(t,u,transitive(2));h=h_subsets(a)[-1]
  require(h==form(v,w),'actual two-input cyclic gate')
  tv=[value(bin_polynomial(t),x) for x in (1,2,3)];uv=[value(bin_polynomial(u),x)for x in (1,2,3)]
  require(h==sum(tv[i]*M[i][j]*uv[j] for i in range(3)for j in range(3)),'bin-evaluation coordinate transform');physical+=1
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   t=tournament(n,code);v=cov3(t);h,b,c=[value(bin_polynomial(t),x)for x in (1,2,3)];out=form(v,v)
   require((-4*c+13*b-11*h)%2==1,'uncancellable odd binary coefficient')
   require(out>=phi(n) and (out==phi(n))==(h==1),'sharp paired order bound/equality')
   require(out%8==5 and out%3==(0 if n%2==0 else (-h*h)%3),'paired congruences');profiles+=1
 require([phi(n)for n in range(1,9)]==[5,45,509,6069,69485,744405,7545389,73529109],'exact order thresholds')
 possible=sorted({form(cov3(t),cov3(t))for t in inputs});require(possible==[5,45,509,621],'complete paired values through2000')
 inverse_controls=0
 for target in range(1,2001):
  ans=inverse(target);require((ans is not None)==(target in possible),'complete restricted inverse')
  if ans:require(h_subsets(tuple(ans['output_rows']))[-1]==target,'inverse literal witness')
  inverse_controls+=1
 t=(78,124,120,112,97,65,0);v=cov3(t);require(v==[23,273,611],'witness input covers')
 a=cyclic(t,converse(t),transitive(2));validate(a);n=len(t);tau=[i+n for i in range(n)]+list(range(n))+[2*n+1,2*n]
 require(all(tau[tau[i]]==i and tau[i]!=i for i in range(len(a))),'fixed-point-free involution')
 require(all(((a[i]>>j)&1)==((a[tau[j]]>>tau[i])&1)for i,j in combinations(range(len(a)),2)),'literal anti-involution')
 require(strong(a),'actual strong graph')
 H=h_subsets(a)[-1];require(H==form(v,v)==37050461,'direct sixteen-vertex count')
 p=5292923;bound=isqrt(p);require(H==7*p and bound==2300,'factor and exact primality bound')
 rem=[p%d for d in range(2,bound+1)];require(all(rem),'complete integer trial division')
 digest=hashlib.sha256(json.dumps(rem,separators=(',',':')).encode()).hexdigest()
 witness={'input_rows':t,'input_cover_vector':v,'rows':a,'anti_involution':tau,'order':len(a),'H':H,'prime_quotient':p,'trial_division_max':bound,'number_of_nonzero_remainders':len(rem),'remainder_sha256':digest,'strong':True,'H_source':'direct vertex-subset Hamilton DP and independent exact fragment kernel'}
 Path(__file__).with_name('three_cover_cyclic_witness.json').write_text(json.dumps(witness,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':'PASS','kernel':K,'run_color_tables':ns,'physical_small_pairs':physical,'profile_controls':profiles,'complete_inverse_targets':inverse_controls,'paired_values_through2000':possible,'witness':witness},sort_keys=True,indent=2))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--target',type=int);args=parser.parse_args()
 if args.target is not None:print(json.dumps({'target':args.target,'paired_witness':inverse(args.target)},sort_keys=True,indent=2))
 else:run()
