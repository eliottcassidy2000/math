"""Independent literal word census and predecessor-DP physical audit."""
from pathlib import Path
from itertools import permutations,combinations
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'even_endpoint_obstruction'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def op(R):return [((1<<len(R))-1)^r^(1<<i) for i,r in enumerate(R)]
def aug(R,a,c):
 n=len(R);out=list(R)+[1<<(n+1),0]
 for v in range(n):
  for z,m in [(n,a),(n+1,c)]:
   if m>>v&1:out[v]|=1<<z
   else:out[z]|=1<<v
 return out
def pair(R):
 return [sum((3<<(2*j)) for j in range(len(R))if R[i]>>j&1)+(1<<(2*i+1)if k==0 else 0)for i in range(len(R))for k in range(2)]
def census(n):
 edges=list(combinations(range(n),2));at={e:k for k,e in enumerate(edges)};N=1<<len(edges);s=bytearray(N);e=bytearray(N)
 for w in permutations(range(n)):
  forced=0;one=0
  for u,v in zip(w,w[1:]):
   b=1<<at[min(u,v),max(u,v)];forced|=b
   if u<v:one|=b
  free=(N-1)^forced;M=free
  while True:
   code=M|one;s[code]^=1<<w[0];e[code]^=1<<w[-1]
   if M==0:break
   M=(M-1)&free
 hist=[0]*(n+1)
 for a,b in zip(s,e):hist[(a^b).bit_count()]+=1
 return [n,N,hist[0],*hist]
# No producer helper imports.
def masks(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M
def main():
 r=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired producer bytes')
 Q=[238,4,240,118,195,146,34,74];need(Q==r['quotient_rows'],'literal rows');n=8
 need(all(not(Q[i]>>i&1)for i in range(n))and all(bool(Q[i]>>j&1)!=bool(Q[j]>>i&1)for i in range(n)for j in range(i+1,n)),'tournament')
 counts={};hp=[]
 for k in range(1,n+1):
  for w in permutations(range(n),k):
   if not all(Q[u]>>v&1 for u,v in zip(w,w[1:])):continue
   M=sum(1<<v for v in w)
   if M not in counts:counts[M]=[[0]*n,[0]*n]
   counts[M][0][w[0]]+=1;counts[M][1][w[-1]]+=1
   if k==n:hp.append(list(w))
 S,E=counts[255];need(hp==r['Hamilton_paths'] and S==r['Hamilton_start_counts']and E==r['Hamilton_end_counts'],'complete literal path words')
 need(sum(S)==249 and all((x+y)%2==0 for x,y in zip(S,E)),'endpoint obstruction')
 bad=[]
 for M in range(1,255):
  if M.bit_count()<2:continue
  members=[u for u in range(n)if M>>u&1]
  if all(len({bool(Q[v]>>u&1)for u in members})==1 for v in range(n)if not M>>v&1):bad.append(M)
 need(not bad and bad==r['proper_nontrivial_modules'],'prime module test')
 h=[[sum(counts[M][0][i]+counts[M][1][i]for M in counts if M.bit_count()%2==0 and M>>i&1 and not M>>j&1)%2 if i!=j else 0 for j in range(n)]for i in range(n)]
 need(h==r['pair_inflation']['shared_cubic_matrix'],'literal even-subset endpoint matrix')
 cs=[census(k)for k in (2,4,6)];need(cs==r['smaller_even_census'],'independent word-completion census')
 bank=r['complete_named_bank'];queries=[row[:2]for row in bank['records']]+[row[:2]for row in bank['full_mask_reconstructions']];half=[m for m in range(256)if m.bit_count()<=2];need(queries[:1369]==[[a,c]for a in half for c in half],'complete bidegree bank')
 pairdata=r['pair_inflation'];T=pair(Q);U=pair(op(Q));need(T==pairdata['source_old_rows']and U==pairdata['target_old_rows'],'pair orientations/physical labels')
 graphs=[];expected=[]
 for (a,c),row in zip(queries,bank['records']+bank['full_mask_reconstructions']):
  for R,z in zip((Q,op(Q)),row[2:4]):graphs.append(aug(R,a,c));expected.append(z)
 for (a,c),row in zip(pairdata['physical_queries'],pairdata['physical_values']):
  for R,z in zip((T,U),row):graphs.append(aug(R,a,c));expected.append(z)
 with tempfile.TemporaryDirectory()as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True)
  data=str(len(graphs))+'\n'+''.join(str(len(R))+' '+' '.join(map(str,R))+'\n'for R in graphs)
  actual=list(map(int,subprocess.check_output([str(exe)],input=data.encode()).split()))
 need(actual==expected and len(actual)==r['physical_counts']==2818,'every exact physical response')
 values={(row[0],row[1]):row[2:4]for row in bank['records']};coeff={};hist=[0]*5
 for row in bank['records']:
  a,c=row[:2];z=[sum((-1)**((a.bit_count()+c.bit_count())-(x.bit_count()+y.bit_count()))*values[x,y][t]for x in masks(a)for y in masks(c))for t in range(2)]
  need(z==row[4:],'full finite difference coefficient');coeff[a,c]=z;hist[a.bit_count()+c.bit_count()]+=z[0]!=z[1]
 need(hist==[int(bank['difference_histogram'][str(k)])for k in range(5)],'full difference histogram')
 for row in bank['full_mask_reconstructions']:
  a,c=row[:2];need(row[2:]==[sum(z[t]for(x,y),z in coeff.items()if x&a==x and y&c==y)for t in range(2)],'full mask reconstruction')
 i,j,k=pairdata['coefficient_indices'];need(i!=j and j==k,'shared cubic names')
 need(pairdata['physical_queries']==[[(3<<(2*i)if b&1 else 0)|(3<<(2*j)if b&2 else 0),3<<(2*j)if b&4 else 0]for b in range(8)],'physical alias battery')
 delta=sum((-1)**(3-b.bit_count())*(x-y)for b,(x,y)in enumerate(pairdata['physical_values']));need(delta==pairdata['exact_contrast_coefficient']==-56897620 and delta%8==4 and h[i][j]==1,'cubic contrast')
 pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest()for f in ['PROOF.md','verify.py','literal.cpp','census.cpp']};need(pins==r['source_pins'],'frozen producer pins')
 return {'status':'PASS','source_pins':pins,'auditor_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in ['audit.py','count.cpp']},'literal_Hamilton_paths':len(hp),'literal_subsets':len(counts),'start':S,'end':E,'prime':True,'census':cs,'physical_counts':len(actual),'physical_values_sha256':hashlib.sha256(json.dumps(actual,separators=(',',':')).encode()).hexdigest(),'complete_named_coordinates':len(coeff),'difference_histogram':hist,'full_mask_reconstructions':len(bank['full_mask_reconstructions']),'pair_shared_cubic':delta,'scope':'Endpoint-parity obstruction, not a response collision.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
