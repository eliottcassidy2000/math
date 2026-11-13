"""Bounded full-root-state controls of the all-core neutral gap law."""
from pathlib import Path
import json
from check import tables,opposite,need
HERE=Path(__file__).resolve().parent

def shell(x):
 z=1
 while x>1:
  q=x//2
  if not x&1:z+=q+(q&1)
  x=q
 return z

def moments(C):
 H,E=tables(C);_,S=tables(opposite(C));full=len(H)-1;n=len(C)
 F=[sum(E[m][i]*H[full^m] for m in range(full+1)) for i in range(n)]
 G=[sum(S[m][i]*H[full^m] for m in range(full+1)) for i in range(n)]
 return H[-1],sum(H[m]*H[full^m] for m in range(full+1)),S[-1],E[-1],F,G

def construct(C,J,epsilon,word,before):
 n=len(C);t=len(word);u=n;w=n+t+1;tail=((1<<w)-1)^((1<<(n+1))-1)
 P=[row|tail|(0 if before else 1<<u) for row in C]
 P.append(tail|(((1<<n)-1) if before else 0))
 P.extend(((1<<w)-1)^((1<<(n+2+i))-1) for i in range(t))
 X=sum(1<<j for j in J)|(epsilon<<u)|sum(1<<(n+1+i) for i,b in enumerate(word) if b==0)
 return tuple(row|((1<<w) if not(X>>i&1) else 0) for i,row in enumerate(P))+(X,)

def profile(T,r=None):
 H,E=tables(T);_,S=tables(opposite(T));full=len(H)-1
 if r is None:r=len(T)-1
 return H[-1],S[-1][r],H[full^(1<<r)],sum(E[m][r]*H[full^m] for m in range(full+1))

def main():
 cores=[(0,),(6,4,0),(2,4,1),(2,4,9,3),(22,8,2,5,14)]
 words=[(),(0,),(1,),(0,1),(1,0),(0,0,1),(1,0,1)]
 checks=0;examples=[]
 for C in cores:
  h,p2,S,E,F,G=moments(C);n=len(C)
  _,end_rows=tables(C);_,start_rows=tables(opposite(C));full=(1<<n)-1
  subsets={(),tuple(range(n)),(0,),tuple(range(1,n)),tuple(range(0,n,2))}
  for J in sorted(subsets):
   K=tuple(i for i in range(n) if i not in J)
   joint=sum(sum(end_rows[m][i] for i in J)*sum(start_rows[full^m][i] for i in J) for m in range(full+1))
   for eps in (0,1):
    for word in words:
     x=sum(b<<i for i,b in enumerate(word));beta=x&1
     post=construct(C,J,eps,word,False);pre=construct(C,J,eps,word,True)
     A=profile(post);B=profile(pre)
     sumj=lambda v:sum(v[i] for i in J)
     sumk=lambda v:sum(v[i] for i in K)
     delta=sumk(F)+sumk(S)+x*sumk(G)-beta*sumk(E) if eps else -(x+1)*sumj(G)-(1-beta)*sumj(E)
     change=(delta,sumk(S),0,sumk(F)) if eps else (delta,-sumj(S),0,-sumj(F))
     need(tuple(b-a for a,b in zip(A,B))==change,'full marked relocation law')
     need(A[1]==sumj(S) and B[1]==eps*h and A[2]==B[2]==h,'absolute start/delete coordinates')
     need(A[3]==(2*x+2-eps)*p2-sumj(F),'absolute ending-cover coordinate')
     need(B[3]==2*(x+1)*p2-2*sumj(F)-eps*h,'prepared ending-cover coordinate')
     need((delta>0)==bool(eps and K) and (delta<0)==bool(not eps and J),'strict sign and equality')
     need(delta%2==0,'even count difference')
     need(post[-1]==pre[-1],'same named newborn cut')
     handoff=None
     if eps:
      handoff=profile(pre,n);qu=(x+1)*sumj(G)+shell(x)*h-(1-beta)*sumj(E)-joint
      need(handoff==(B[0],qu-sumj(S),qu,B[0]+sumj(S)),'retained-u handoff fullstate')
      need([i for i,row in enumerate(pre) if row>>n&1]==[len(pre)-1],'retained-u unique incoming port')
     checks+=1
     if C==(2,4,1) and J==(0,) and eps==1 and word==(1,):
      examples.append(dict(core=C,J=J,epsilon=eps,word=word,post_rows=post,pre_rows=pre,post_profile=A,pre_profile=B,change=change,retained_u_profile=handoff))
 # Exact inverse examples: each parity determines at most one nonnegative suffix address.
 inverses=[]
 C=(2,4,1);h,p2,S,E,F,G=moments(C);K=(1,2);f=sum(F[i] for i in K);s=sum(S[i] for i in K);g=sum(G[i] for i in K);e=sum(E[i] for i in K)
 for D in (8,12,20,24,48):
  hits=[]
  for beta in(0,1):
   offset=f+s+beta*(g-e)
   if D>=offset and (D-offset)%(2*g)==0:hits.append(2*((D-offset)//(2*g))+beta)
  need(all(f+s+x*g-(x&1)*e==D for x in hits),'inverse difference gate')
  inverses.append(dict(desired_difference=D,addresses=hits))
 out=dict(status='PASS',scope='Five specified nonempty cores, selected entry sets, epsilon0/1, seven specified suffixes; literal full-state controls',
  cores=cores,suffixes=words,full_state_controls=checks,examples=examples,difference_inverse_examples=inverses)
 (HERE/'neutral_gap_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','full_state_controls':checks,'examples':examples,'inverse':inverses},indent=2))

if __name__=='__main__':main()
