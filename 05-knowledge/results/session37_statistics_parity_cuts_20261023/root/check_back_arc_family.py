#!/usr/bin/env python3
import itertools,json,math

def need(x,msg):
 if not x:raise RuntimeError(msg)
def make(m):
 n=3*m;A=[[0]*n for _ in range(n)]
 for i,j in itertools.combinations(range(n),2):
  w=i if i//3!=j//3 or (j-i)%3==1 else j;l=j if w==i else i;A[w][l]=1
 A[0][-1]=0;A[-1][0]=1
 return A

def H(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 outs=[sum(A[i][j]<<j for j in range(n)) for i in range(n)]
 for mask in range(1,1<<n):
  for i in range(n):
   v=D[mask][i]
   if not v:continue
   rest=outs[i]&~mask
   while rest:
    b=rest&-rest;rest-=b;D[mask|b][b.bit_length()-1]+=v
 return sum(D[-1])

def cycle_poly(A):
 n=len(A);C={};N=1<<n
 for root in range(n):
  D={(1<<root,root):1}
  for mask in range(1<<root,N):
   if mask&((1<<(root+1))-1)!=(1<<root):continue
   count=0
   for i in range(root,n):
    val=D.get((mask,i),0)
    if not val:continue
    if A[i][root]:count+=val
    for j in range(root+1,n):
     if not mask>>j&1 and A[i][j]:D[mask|1<<j,j]=D.get((mask|1<<j,j),0)+val
   if mask.bit_count()>=3 and mask.bit_count()%2 and count:C[mask]=count
 P=[None]*N;P[0]=[1]
 for mask in range(1,N):
  v=mask&-mask;p=P[mask^v][:];sub=mask
  while sub:
   if sub&v and sub in C:
    q=P[mask^sub]
    if len(p)<len(q)+1:p.extend([0]*(len(q)+1-len(p)))
    for k,a in enumerate(q):p[k+1]+=C[sub]*a
   sub=(sub-1)&mask
  P[mask]=p
 return P[-1]
def formula(m):
 p=[math.comb(m,k) for k in range(m+1)]
 for j in range(m-1):
  num=math.comb(m-2,j)*(9*10**(m-2-j)-(-2)**(m-2-j))
  need(num%2==0,'nonintegral coefficient');p[j+1]+=num//2
 return p
def val(p,z):return sum(c*z**k for k,c in enumerate(p))
def main():
 rows=[]
 for m in range(2,6):
  A=make(m);p=formula(m);h=H(A);need(h==val(p,2),'H mismatch')
  direct=cycle_poly(A) if m<=4 else None
  if direct is not None:need(direct==p,'cycle packing mismatch')
  rows.append(dict(m=m,n=3*m,H=h,polynomial=p,direct_cycle_check=m<=4))
 tests=0
 for m in range(3,31):
  p=formula(m);need((m-1)*p[-2]**2<2*m*p[-3],'Newton failure absent')
  if m>=5:
   for z in range(-5,6):
    rhs=(3*z+9)*val(formula(m-1),z)-(3*z*z+18*z-12)*val(formula(m-2),z)+(z**3+9*z*z-12*z-20)*val(formula(m-3),z)
    need(val(p,z)==rhs,'recurrence');tests+=1
 need(val(formula(4),2)!=15*val(formula(3),2)-36*val(formula(2),2),'exception lost')
 print(json.dumps(dict(status='PASS',direct_rows=rows,recurrence_controls=tests,newton_controls=28,scope='Direct H n6,9,12,15; direct odd cycle packings n6,9,12; exact recurrence controls only thereafter'),indent=2,sort_keys=True))
if __name__=='__main__':main()
