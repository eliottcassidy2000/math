from graphs import SEED,birth
PORTS=[0,1,4,5,6]
SEED_B=[2,12,145,149,355,399,175,19,207]
def coefficient(t):return 40*3**t+160*2**t-20,2**(t+1)-5
def value(m,t):
 A,C=coefficient(t);return A*3**m+C

def canonical(m,t):
 rows=SEED[:];ports=PORTS[:]
 for _ in range(m-1):rows,ports=birth(rows,ports,'R',7,5)
 rows,ports=birth(rows,ports,'L',12,5)
 for _ in range(t):rows,ports=birth(rows,ports,'L',24,5)
 return rows

def physical(word):
 rows=SEED_B[:];ports=PORTS[:];roots=[7];leaves=[]
 for side in word:
  v=len(rows);rows,ports=birth(rows,ports,side,7 if side=='R'else 24,5);(roots if side=='R'else leaves).append(v)
 return rows,ports,list(range(7))+roots+[8]+leaves

def inverse(N):
 ans=[];t=0
 while True:
  A,C=coefficient(t)
  if 3*A+C>N:return ans
  z,rem=divmod(N-C,A)
  if rem==0:
   m=0
   while z>0 and z%3==0:z//=3;m+=1
   if z==1 and m>=1:ans.append((m,t))
  t+=1
