"""All-order argument controls; complete finite data are not extrapolated."""
from pathlib import Path
from itertools import combinations,product
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
from itertools import permutations

def op(R):return [((1<<len(R))-1)^r^(1<<i)for i,r in enumerate(R)]
def starts(R,twice=None):
 n=len(R);N=1<<n;dp={}
 for i in range(n):dp[1<<i,0,i]=[int(i==j)for j in range(n)]
 for M in range(1,N):
  for rep in range(2):
   for last in range(n):
    z=dp.get((M,rep,last))
    if z is None:continue
    for nxt in range(n):
     if not R[last]>>nxt&1:continue
     if M>>nxt&1:
      if nxt!=twice or rep:continue
      key=M,1,nxt
     else:key=M|1<<nxt,rep,nxt
     y=dp.setdefault(key,[0]*n)
     for i in range(n):y[i]+=z[i]
 return [sum(dp.get((N-1,int(twice is not None),last),[0]*n)[i]for last in range(n))for i in range(n)]
def refine(Q):
 U=op(Q);S=starts(Q);E=starts(U);rS=[starts(Q,j)for j in range(len(Q))];rE=[starts(U,j)for j in range(len(Q))];sigma=[(S[i]-E[i]+2*sum(x[i]-y[i]for x,y in zip(rS,rE)))%4 for i in range(len(Q))]
 return S,E,sigma
def rows(n,code):
 R=[0]*n
 for k,(i,j)in enumerate(combinations(range(n),2)):
  if code>>k&1:R[i]|=1<<j
  else:R[j]|=1<<i
 return R
def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 return [sum(((1<<len(C))-1)<<off[j]for j,C in enumerate(B)if Q[i]>>j&1)|(r<<off[i])for i,b in enumerate(B)for r in b],off
def aug(R,a,c=0):
 n=len(R);out=list(R)+[1<<(n+1),0]
 for i in range(n):
  for v,m in [(n,a),(n+1,c)]:
   if m>>i&1:out[i]|=1<<v
   else:out[v]|=1<<i
 return out
def physical(Q,B,i,K=None):
 R,off=blow(Q,B);a=((1<<len(B[i]))-1)<<off[i];c=0
 if K is not None:
  R,z=blow(K,[R]+[tt(1)]*(len(K)-1));outside=sum(1<<z[j]for j in range(1,len(K))if K[j]&1);a|=outside;c|=outside
 return aug(R,a,c)
def multiset_words(n,j,copies):
 for positions in combinations(range(n+copies-1),copies):
  if any(v==u+1 for u,v in zip(positions,positions[1:])):continue
  for other in permutations([i for i in range(n)if i!=j]):
   it=iter(other);yield tuple(j if k in positions else next(it)for k in range(n+copies-1))
def swap_returns(w,j):
 a,b,c=[k for k,v in enumerate(w)if v==j]
 return w[:a+1]+w[b+1:c+1]+w[a+1:b+1]+w[c+1:]
def module_masks(Q):
 n=len(Q);return [M for M in range(1,(1<<n)-1)if M.bit_count()>1 and all(Q[v]&M in(0,M)for v in range(n)if not M>>v&1)]
def main():
 table=[];literal_r=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   Q=rows(n,code);S,E,sigma=refine(Q);p=sum(((S[i]+E[i])&1)<<i for i in range(n));ref=sum((sigma[i]==2 and not(p>>i&1))<<i for i in range(n));need(all(z%2==(p>>i&1)for i,z in enumerate(sigma)),'sigma parity')
   if n<=4:
    literal=[w for w in permutations(range(n))if all(Q[u]>>v&1 for u,v in zip(w,w[1:]))];need(S==[sum(w[0]==i for w in literal)for i in range(n)]and E==[sum(w[-1]==i for w in literal)for i in range(n)],'literal endpoint recurrence')
    for j in range(n):
     words=[w for w in multiset_words(n,j,2)if all(Q[u]>>v&1 for u,v in zip(w,w[1:]))];need(starts(Q,j)==[sum(w[0]==i for w in words)for i in range(n)],'literal indistinguishable two-run recurrence');literal_r+=1
   table.append([n,code,p,ref,*sigma])
 involution=[]
 for n in range(3,7):
  count=0
  for j in range(n):
   for w in multiset_words(n,j,3):
    v=swap_returns(w,j);need(v!=w and swap_returns(v,j)==w,'free return-segment involution');need(w[0]==v[0]and w[-1]==v[-1]and sorted(zip(w,w[1:]))==sorted(zip(v,v[1:])),'same quotient arc conditions and endpoints');count+=1
  involution.append([n,count])
 C3=[2,4,1];types=[tt(2),tt(4),tt(3),C3];graphs=[];records=[]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   Q=rows(n,code);S,E,sigma=refine(Q)
   for typ in product(range(4),repeat=n):
    B=[types[t]for t in typ]
    for i in range(n):
     if(S[i]+E[i])%2:continue
     idx=len(graphs);graphs.extend([physical(Q,B,i),physical(op(Q),B,i)]);records.append({'order':n,'quotient_code':code,'types':typ,'anchor':i,'sigma':sigma[i],'indices':[idx,idx+1]})
 Q5=[20,1,2,7,14];Q8=[238,4,240,118,195,146,34,74];regular=[24,17,3,6,12];special=[]
 for name,Q,B,i,K in [('W5_pairs',Q5,[tt(2)]*5,1,None),('W5_C3',Q5,[C3]*5,2,None),('W5_mixed',Q5,[C3,tt(4),tt(2),tt(3),tt(2)],3,None),('W5_mixed_ambient',Q5,[C3,tt(4),tt(2),tt(3),tt(2)],4,C3),('prime8_pairs',Q8,[tt(2)]*8,0,None)]:
  S,E,sigma=refine(Q);need(not module_masks(Q),'prime example');idx=len(graphs);graphs.extend([physical(Q,B,i,K),physical(op(Q),B,i,K)]);special.append({'name':name,'quotient':Q,'blocks':B,'anchor':i,'ambient':K,'start':S,'end':E,'sigma':sigma,'indices':[idx,idx+1]})
 need(refine(Q5)[2]==[0,2,2,2,2]and refine(Q8)[2]==[2,0,0,0,0,2,2,2],'displayed refined vectors');need(refine(C3)[2]==[0]*3 and refine(regular)[2]==[0]*5,'genuine blind boundaries')
 need(all(Q5[u]>>v&1 for k,u in enumerate([4,3,2,1])for v in [4,3,2,1][k+1:])and Q5[0]==(1<<4)|(1<<2),'alternating W5 label map')
 idx=len(graphs);graphs.extend([aug(Q8,1),aug(op(Q8),1)]);singleton=[idx,idx+1]
 with tempfile.TemporaryDirectory()as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True)
  data=str(len(graphs))+'\n'+''.join(str(len(R))+' '+' '.join(map(str,R))+'\n'for R in graphs);values=[list(map(int,s.split()))for s in subprocess.check_output([str(exe)],input=data.encode()).decode().splitlines()]
 need(len(values)==len(graphs),'all physical outputs')
 for rec in records:
  v=[values[k]for k in rec.pop('indices')];need([(x-y)%8 for x,y in zip(*v)]==[2*rec['sigma']%8,0,0],'universal complete physical profile control');rec['values']=v
 for rec in special:
  v=[values[k]for k in rec.pop('indices')];need([(x-y)%8 for x,y in zip(*v)]==[4,0,0],'special prime/ambient full profile');rec['values']=v
 sing=[values[k]for k in singleton];need(sing[0][0]==771 and sing[1][0]==2403 and(sing[0][0]-sing[1][0])%8==0,'singleton boundary')
 deps=['uniform_block_starts/PROOF.md','parity_profile_transfer/PROOF.md','mod8_transfer/PROOF.md','module_decomposition/PROOF.md','monotone_dual_queries/PROOF.md','highest_quotient_reduction/PROOF.md','prime_eight_uniform_anchor/PROOF.md']
 return {'status':'PASS','quotient_table_columns':['order','edge_code','endpoint_parity_mask','refined_anchor_mask','sigma_vector...'],'complete_quotient_table':table,'literal_two_run_controls':literal_r,'involution_word_counts':involution,'complete_small_physical_controls':records,'special_controls':special,'singleton_boundary_values':sing,'physical_graph_occurrences':len(graphs),'distinct_physical_graphs':len({tuple(R)for R in graphs}),'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in ['PROOF.md','verify.py','physical.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()for f in deps},'scope':'All-order sufficient second-endpoint transfer for arbitrary nontrivial children and fixed ambient; no completeness or collision claim.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
