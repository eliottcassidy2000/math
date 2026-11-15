from itertools import product
from math import gcd

B=(44,121,10,0,13,92,29)
A=(0,)
K=3315840

def table(rows):
 n=len(rows); size=1<<n
 f=[[0]*n for _ in range(size)]
 r=[[0]*n for _ in range(size)]
 for x in range(n): f[1<<x][x]=r[1<<x][x]=1
 for m in range(1,size):
  for x in range(n):
   for arr,mask in ((f,rows[x]),(r,~rows[x])):
    p=mask&(size-1)&~m
    while p:
     b=p&-p;p-=b
     arr[m|b][b.bit_length()-1]+=arr[m][x]
 h=[1]+[sum(f[m]) for m in range(1,size)]
 return h,r,f

def literal(board):
 tables=[]
 for rows,mask in board:
  h,r,f=table(rows); all=len(h)-1
  starts=[sum(z for x,z in enumerate(row) if (mask>>x)&1) for row in r]
  ends=[sum(z for x,z in enumerate(row) if not ((mask>>x)&1)) for row in f]
  tables.append((h,all,starts,ends))
 ans=0; n=len(board)
 # i is the last block with a nonempty left side, or -1;
 # j is the first with a nonempty right side, or n.
 for i,j in product(range(-1,n),range(n+1)):
  out=1
  for b,(h,all,starts,ends) in enumerate(tables):
   total=0
   for s in range(all+1):
    t=all^s
    if (b>i and s) or (b==i and not s): continue
    if (b<j and t) or (b==j and not t): continue
    total+=(ends[s] if b==i else h[s])*(starts[t] if b==j else h[t])
   out*=total
  ans+=out
 return ans

def local(mask):
 h,r,f=table(B); all=len(h)-1
 q=sum(h[m]*h[all^m] for m in range(all+1))
 vals=[0]*4
 for m in range(1,all+1):
  s=sum(z for x,z in enumerate(r[m]) if (mask>>x)&1)*h[all^m]
  t=sum(z for x,z in enumerate(f[m]) if not ((mask>>x)&1))*h[all^m]
  sg=1 if m.bit_count()&1 else -1
  for i,x in enumerate((s,sg*s,t,sg*t)): vals[i]+=x
 return h[-1],q,literal([(B,mask)]),*vals

def transfer(masks):
 H,N,S,d=(1,1,1,1)
 for mask in masks:
  T,e={2:(630,32),32:(624,30)}[mask]
  H,N,S,d=(37*H,37*N+118*H+T*S+e*d,704*S+118*H,0)
 return H,N+S+d,2*S,0

def board(bits):
 masks=[2]
 for x in bits: masks.extend((32,2) if x else (2,32))
 return masks

def main():
 assert local(2)==(37,704,155,118,0,630,32)
 assert local(32)==(37,704,155,118,0,624,30)
 for bits in product(range(2),repeat=1):
  masks=board(bits)
  graph=[(A,1)]+[(B,x) for x in masks]+[(A,0)]
  assert literal(graph)==transfer(masks)[1]
 m=8; reps=set(); N0=transfer(board([0]*m))
 for bits in product(range(2),repeat=m):
  H,N,S,d=transfer(board(bits)); reps.add(N)
  assert (H,S,d)==(N0[0],N0[2],N0[3])
  assert N-N0[1]==K*sum(b*37**(2*(m-1-j))*704**(2*j) for j,b in enumerate(bits))
  assert N%37==(16 if bits[-1] else 5)
 assert len(reps)==2**m
 assert gcd(*(x-N0[1] for x in reps))==K
 for m in range(1,160):
  P,N,_,_=transfer(board([0]*m))
  for p,t in [(3,1+2*m),(5,2),(11,4*m),(157,(69+77*m)*pow(37,-1,157))]:
   assert (N-P*t)%p==0
 P,N,_,_=transfer(board([0,0]))
 assert gcd(N,37*K)==1
 print("Literal partition replay, positional digits, and primitive residue formulas passed")
if __name__=='__main__': main()
