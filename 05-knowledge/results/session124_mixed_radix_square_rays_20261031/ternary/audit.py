#!/usr/bin/env python3
"""Independent ternary composition, cycle, decoding and frontier controls."""
import argparse,hashlib,itertools,json,math
from pathlib import Path
OPS={'A':(25,12,(1,3)),'B':(25,8,(1,3,6,19)),'C':(49,4,(1,3,6,19))};checks=0

def need(x,msg):
 global checks
 checks+=1
 if not x:raise ValueError(msg)
def square(x):return x>0 and math.isqrt(x)**2==x
def product(w):return math.prod(OPS[x][0] for x in w)
def threshold(w):return (85*product(w)-1)//2
def prefix(w):
 R=1;C=0;out=[]
 for letter in w:
  Q,c,H=OPS[letter];out.extend(R*v+(-C if j%2==0 else C) for j,v in enumerate(H));C+=R*c;R*=Q
 out.extend([R-C,3*R+C]);return out

def decode(P):
 result=[]
 while len(P)>2:
  if P[2]==13:letter='A'
  elif P[2]==6 and len(P)>=6 and P[4] in [17,45]:letter='B' if P[4]==17 else 'C'
  else:raise ValueError('uncodable prefix')
  Q,c,H=OPS[letter];need(P[:len(H)]==list(H),'decoded literal prefix');tail=P[len(H):];following=[]
  for j,v in enumerate(tail):
   x=v+(c if j%2==0 else -c);need(x%Q==0,'exact inverse affine decoding');following.append(x//Q)
  result.append(letter);P=following
 need(P==[1,3],'decoding terminates at base');return ''.join(result)

def validate_pair(pair,pref):
 n=len(pair[0]);need(len(pair[1])==n+1,'mixed consecutive intervals');color={}
 for i,P in enumerate(pair):
  need(P[:len(pref)]==pref and len(P)==n+i and set(P)==set(range(1,n+i+1)),'mixed input coverage/prefix')
  need(P[-1]==(8 if len(P)%2==0 else 10),'mixed input endpoint')
  for j,v in enumerate(P):
   if not i:color[v]=j%2
   elif v<=n:need(color[v]==j%2,'mixed input compatibility')
   if j:need(square(P[j-1]+v),'mixed input square edge')
def stream(Q,pair,atoms):
 for atom in atoms:
  kind,c=atom[:2]
  if kind=='single':yield c
  elif kind=='head':yield Q-c
  else:
   _,c,i,o=atom;P=pair[i];lo=int(kind=='tail');positions=range(lo,len(P)) if o==1 else range(len(P)-1,lo-1,-1)
   for j in positions:yield Q*P[j]+(c if j%2 else -c)

def leaves(N,complete_ternary=False):
 todo=[('',1)];answer=[]
 while todo:
  w,R=todo.pop();can_expand=85*(49 if complete_ternary else 25)*R<=2*N+1
  if can_expand:
   todo.extend([(w+'A',25*R),(w+'B',25*R)])
   if 85*49*R<=2*N+1:todo.append((w+'C',49*R))
  else:answer.append((w,R))
 return answer

def count_dynamic(N):
 # Bottom-up on distinct products, independently of word-frontier traversal.
 active={1};todo=[1]
 while todo:
  R=todo.pop()
  for Q in [25,49]:
   nxt=R*Q
   if 85*nxt<=2*N+1 and nxt not in active:active.add(nxt);todo.append(nxt)
 count={}
 for R in sorted(active,reverse=True):
  count[R]=2*count[25*R]+(count[49*R] if 49*R in count else 0) if 25*R in count else 1
 return count[1],len(count)

def main():
 ap=argparse.ArgumentParser()
 for x in ['a-paired','a-cycles','b-paired','b-cycles','c-paired','c-cycles','small-base-file']:ap.add_argument('--'+x,type=Path,required=True)
 args=ap.parse_args();small=json.loads(args.small_base_file.read_text());hashes={}
 def cert(letter,p,d,terminal=False):
  folder=getattr(args,letter.lower()+('_cycles' if terminal else '_paired'));stem=('cycle' if letter=='C' else'C') if terminal else'B';name=f'{stem}_p{p}_d{d}.json';raw=(folder/name).read_bytes();hashes[letter+'/'+('cycle/' if terminal else'pair/')+name]=hashlib.sha256(raw).hexdigest();return json.loads(raw)
 def build(n,w):
  need(n>=threshold(w),'mixed prefix domain')
  if not w:
   need(str(n) in small,'finite control reaches certified small base');pair=small[str(n)]
  else:
   letter=w[0];Q=OPS[letter][0];h=(Q-1)//2;m,d=divmod(n-h,Q);d+=h
   source=build(m,w[1:]);witness=cert(letter,m%2,d);pair=[list(stream(Q,source,x)) for x in witness['words']]
  need(len(pair[0])==n,'mixed paired order');validate_pair(pair,prefix(w));return pair
 word_controls=list(itertools.product('ABC',repeat=2))+list(itertools.permutations('ABC'))+[tuple(x*3) for x in 'ABC'];requests=[]
 for w in word_controls:
  for extra in ([0,1] if len(w)==2 else[0]):requests.append((''.join(w),threshold(w)+extra))
 for w in ['AB','BC','CA']:requests.append((w,threshold(w)+OPS[w[0]][0]))
 controls=[]
 for w,N in requests:
  Q=OPS[w[0]][0];h=(Q-1)//2;m,d=divmod(N-h,Q);d+=h;pair=build(m,w[1:]);witness=cert(w[0],m%2,d,True);atoms=witness['word'] if w[0]=='C' else witness['words'][0]
  seen=bytearray(N+1);first=last=None;count=0;wanted=prefix(w);initial=[];H=hashlib.sha256()
  for v in stream(Q,pair,atoms):
   need(1<=v<=N and not seen[v],'mixed cycle ownership');seen[v]=1;count+=1;H.update((str(v)+'\n').encode())
   if len(initial)<len(wanted):initial.append(v)
   if first is None:first=v
   if last is not None:need(square(last+v),'mixed cycle square')
   last=v
  need(count==N and sum(seen)==N and initial==wanted,'mixed cycle complete prefix');need(square(last+first),'mixed cycle closure')
  controls.append({'word':w,'N':N,'prefix':wanted,'cycle_decimal_lines_sha256':H.hexdigest()})
 decoding_cases=0
 for l in range(1,8):
  distinct=set()
  for w in itertools.product('ABC',repeat=l):
   w=''.join(w);P=prefix(w);need(decode(P)==w,'ternary unique word recovery');need(len(set(P))==len(P) and min(P)>0,'mixed positive simple prefix');need(all(square(x+y) for x,y in zip(P,P[1:])),'mixed closed-form square sums');need(len(P)==2+2*w.count('A')+4*(w.count('B')+w.count('C')),'mixed length formula');distinct.add(tuple(P));decoding_cases+=1
  need(len(distinct)==3**l,'distinct coded prefixes')
 frontier_controls=[]
 for N in [1062,2082,26562,52062,102042,664062,16601562,10**9,10**12]:
  full=leaves(N);words={w for w,R in full};dynamic,states=count_dynamic(N);need(len(full)==dynamic,'exact frontier enumeration versus bottom-up count')
  for w,R in full:
   need(w and 85*R<=2*N+1 and 85*25*R>2*N+1,'exact frontier leaf conditions')
   for k in range(len(w)):need(w[:k] not in words,'frontier prefix antichain')
  if N>=2082:
   weighted=leaves(N,True);polynomial={}
   for w,R in weighted:
    need(85*R<=2*N+1<85*49*R,'weighted frontier budget band');u=len(w)-w.count('C');v=w.count('C')
    for k in range(v+1):polynomial[u+k]=polynomial.get(u+k,0)+math.comb(v,k)*(-2)**k
   need({k:v for k,v in polynomial.items() if v}=={0:1},'exact weighted mass after v=1-2u substitution');need(len(full)>=len(weighted),'full frontier at least weighted frontier')
  frontier_controls.append({'N':N,'maximum_coded_antichain':dynamic,'distinct_product_states':states})
 expected={1062:2,2082:3,26562:5,102042:9,16601562:35,10**9:145}
 for row in frontier_controls:
  if row['N'] in expected:need(row['maximum_coded_antichain']==expected[row['N']],'published frontier example')
 neighbors={6,prefix('B')[4],prefix('C')[4]};need(neighbors=={6,17,45},'first branch coexistence would violate cycle degree')
 print(json.dumps({'status':'passed','checks':checks,'mixed_cycle_controls':controls,'decoding_words_checked':decoding_cases,'decoding_max_depth':7,'frontier_controls':frontier_controls,'certificate_sha256':hashes,'producer_imports':0},indent=2,sort_keys=True))
if __name__=='__main__':main()
