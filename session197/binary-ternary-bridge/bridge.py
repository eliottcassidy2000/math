from graphs import SEED,birth
PORTS=[0,1,4,5,6]
def graph(m,t):
 if m<1 or t<0:raise ValueError('m>=1,t>=0 required')
 n=7+m+t;rows=[x&127 for x in SEED[:7]]+[0]*(m+t)
 for r in range(7,7+m):
  rows[r]=19
  for v in range(7):
   if not 19>>v&1:rows[v]|=1<<r
  for z in range(r+1,7+m):rows[r]|=1<<z
 for a in range(7+m,n):
  rows[a]=31|sum(1<<r for r in range(7,7+m))
  rows[5]|=1<<a;rows[6]|=1<<a
  for b in range(a+1,n):rows[a]|=1<<b
 return rows

def physical(word):
 rows=SEED[:];ports=PORTS[:];rs=[7];ls=[]
 for side in word:
  label=len(rows);rows,ports=birth(rows,ports,side,7 if side=='R'else 24,5)
  (rs if side=='R'else ls).append(label)
 return rows,ports,list(range(7))+rs+ls[::-1]

def value(m,t):return 60*2**t*3**m+(1 if t else -1)
def inverse(N):
 if N<=1:return None
 z=N-1;a=b=0
 while z%2==0:a+=1;z//=2
 while z%3==0:b+=1;z//=3
 return (b-1,a-2)if z==5 and a>=3 and b>=2 else None
