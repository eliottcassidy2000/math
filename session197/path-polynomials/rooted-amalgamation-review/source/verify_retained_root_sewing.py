#!/usr/bin/env python3
from verify import *
from collections import deque
from random import Random
import hashlib
import argparse
def state(A,a):
 n=len(A);M=1<<n;dp=[[0]*n for _ in range(M)];sp=[[0]*n for _ in range(M)]
 for v in range(n):dp[1<<v][v]=1
 sp[1<<a][a]=1
 for S in range(1,M):
  for v in range(n):
   U=A[v]&~S
   while U:
    b=U&-U;U-=b;w=b.bit_length()-1;dp[S|b][w]+=dp[S][v];sp[S|b][w]+=sp[S][v]
 h=[sum(z)for z in dp];h[0]=1;F=M-1;D=F^(1<<a)
 Rs=sum(sum(sp[S])*h[F^S]for S in range(M));Re=sum(dp[S][a]*h[F^S]for S in range(M))
 b=sum(h[S]*h[D^S]for S in range(M)if S&D==S);bfull=sum(h[S]*h[F^S]for S in range(M))
 return (h[F],dp[F][a],sum(sp[F]),Rs,Re,h[D],b,int(n==1)),bfull
def sew(A,a,B,b):
 C=list(join(A,B));C[a]^=1<<(len(A)+b);C[len(A)+b]^=1<<a;validate(C);return tuple(C)

def four(A,a):
 z,_=state(A,a);h,e,_,R,_,d,_,_=z
 return (h,R-d,d-e,d)

def signature(A,a):
 z,b=state(A,a)
 return (z[0],z[1],z[2],z[3],z[4],b)

def literal_signatures(A):
 n=len(A);h=0;E=[0]*n;S=[0]*n;Q=[0]*n;L=[0]*n;C=0
 for word in permutations(range(n)):
  legal=[bool(A[word[j]]>>word[j+1]&1)for j in range(n-1)]
  if all(legal):h+=1;E[word[-1]]+=1;S[word[0]]+=1
  for cut in range(n+1):
   if all(legal[:max(0,cut-1)])and all(legal[cut:]):
    C+=1
    if cut<n:Q[word[cut]]+=1
    if cut>0:L[word[cut-1]]+=1
 return [(h,E[a],S[a],Q[a],L[a],C)for a in range(n)]

def positive_transfer(st,sig):
 h,x,y,d=st;H,E,S,Q,L,C=sig
 return(H*h+L*x+S*y+(L-S)*d,C*x+(C-Q-H)*d,(H-E)*d,H*d)

def make_bank(nmax):
 out={};graphs=roots=0
 for n in range(1,nmax+1):
  for c in range(1<<(n*(n-1)//2)):
   A=tournament(n,c);literal=literal_signatures(A);graphs+=1
   for a in range(n):
    sig=signature(A,a);require(sig==literal[a],'independent literal gadget signature');roots+=1
    H,E,S,Q,L,C=sig
    require(C>=2*Q>=2*H and L>=H>=S and H>=E and C>=(n+1)*H,'physical positivity inequalities')
    require((C-Q-H==0)==(n==1),'exact singleton activation boundary')
    out.setdefault(sig,(A,a))
 return list(out.items()),graphs,roots

def inverse(M,letters):
 zero=(1,0,0,1);seen={zero:None};queue=deque([zero]);hits={1:zero};depth={zero:0}
 while queue:
  st=queue.popleft()
  if st[0]+st[1]>M:continue
  for k,(sig,_)in enumerate(letters):
   z=positive_transfer(st,sig)
   if z[0]>M or z in seen:continue
   h,x,y,d=z
   require(1<=d<=h and 0<=y<=d and x>=0 and h%2==1,'all reached states physical inequalities')
   seen[z]=(st,k);depth[z]=depth[st]+1
   require(depth[z]<=h.bit_length(),'exact normalized depth bound')
   queue.append(z);hits.setdefault(h,z)
 digest=hashlib.sha256(json.dumps(sorted(seen),separators=(',',':')).encode()).hexdigest()
 return seen,hits,{'bound':M,'signatures':len(letters),'reachable_states':len(seen),'max_representative_depth':max(depth.values()),'attained_odds':len(hits),'missing_odds':[h for h in range(1,M+1,2)if h not in hits],'state_sha256':digest}

def graph_word(letters,word):
 A=(0,)
 for k in word:
  B,b=letters[k][1];A=sew(A,0,B,b)
 return A

def witness(target,letters,seen,hits):
 z=hits[target];word=[]
 while seen[z] is not None:z,k=seen[z];word.append(k)
 word.reverse();A=graph_word(letters,word);st=four(A,0)
 require(st==hits[target]and h_subsets(A)[-1]==target,'literal witness and complete four-state')
 return {'target':target,'word':word,'rows':A,'root':0,'four_state':st}

def run(bound):
 letters,graphs,roots=make_bank(5)
 require((graphs,roots,len(letters))==(1099,5405,65),'complete order-five labelled inventory')
 letters4=[v for v in letters if len(v[1][0])<=4];require(len(letters4)==19,'complete order-four inventory')
 small=[(tournament(n,c),r)for n in range(1,4)for c in range(1<<(n*(n-1)//2))for r in range(n)]
 pairs=0
 for A,a in small:
  for B,b in small:
   require(four(sew(A,a,B,b),a)==positive_transfer(four(A,a),signature(B,b)),'actual tiny pair full-state transfer');pairs+=1
 extras=0;rng=Random(1971504)
 for turn in range(150):
  n,m=rng.randrange(1,6),rng.randrange(1,6);A=tournament(n,rng.randrange(1<<(n*(n-1)//2)));B=tournament(m,rng.randrange(1<<(m*(m-1)//2)));a,b=rng.randrange(n),rng.randrange(m)
  require(four(sew(A,a,B,b),a)==positive_transfer(four(A,a),signature(B,b)),'larger actual retained-root transfer');extras+=1
 # Identity and the precise two rooted edge letters, not arbitrary formal states.
 F=(transitive(2),0);K=(transitive(2),1);single=((0,),0)
 inputs=[(0,),sew((0,),0,*F),sew((0,),0,*K),sew(sew((0,),0,*F),0,*F)]
 contexts=[[],[single],[single,single],[K]];matrix=[]
 for word in contexts:
  row=[]
  for A in inputs:
   for B,b in word:A=sew(A,0,B,b)
   row.append(h_subsets(A)[-1])
  matrix.append(row)
 require(matrix==[[1,1,3,3],[1,3,5,9],[1,5,9,19],[3,5,9,15]]and det(matrix)==16,'physical rank-four continuation minor')
 # Full-state equality of initial singleton padding, and arbitrary subsequent bank words.
 prefix=(0,)
 for k in range(1,7):
  prefix=sew(prefix,0,(0,),0);require(four(prefix,0)==(1,0,0,1),'exact initial padding transport')
 seen4,hits4,receipt4=inverse(2000,letters4)
 require(receipt4['reachable_states']==20389 and receipt4['max_representative_depth']==10 and receipt4['missing_odds']==[7,21,173],'complete order-four bank spectrum')
 seen,hits,receipt=inverse(bound,letters)
 if bound==5000:
  require(receipt['reachable_states']==152044 and receipt['max_representative_depth']==12 and receipt['missing_odds']==[7,21],'complete order-five bank spectrum')
 targets=[67,163,173,179,181,247,511,539,541,821,4999]
 witnesses=[witness(t,letters,seen,hits)for t in targets if t<=bound and t in hits]
 inventory=[{'index':k,'signature':sig,'rows':A,'root':a}for k,(sig,(A,a))in enumerate(letters)]
 print(json.dumps({'status':'PASS','literal_gadget_graphs':graphs,'literal_rooted_signatures':roots,'tiny_actual_pairs':pairs,'larger_actual_pairs':extras,'continuation_matrix':matrix,'continuation_determinant':16,'order_four_census':receipt4,'order_five_census':receipt,'inventory':inventory,'witnesses':witnesses},indent=2,sort_keys=True))

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--bound',type=int,default=5000);args=parser.parse_args();require(args.bound>=1,'positive bound');run(args.bound)
