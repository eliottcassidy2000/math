"""Three physical symmetric two-integer families of strong C2 tournaments."""
from core import birth,constants,require

def g(u):return constants(u)[0]
def scalar(kind,u,v):
 require(u>=1 and v>=1,'positive parameters');x,y=g(u),g(v);p,t=u%2,v%2
 if kind=='N':
  require(u>=2 and v>=2,'N needs nonempty prefix and suffix')
  return x*y+(u-p)*y+(v-t)*x+3*u*v-u-v-p-t+p*t
 if kind=='M':
  require(u>=2 and v>=2,'M release uses nonempty prefix and suffix')
  return x*y+2*u*y+2*v*x+6*u*v-u-v-u*t-v*p-p-t
 require(kind=='P','unknown family')
 return x*y+4*v*x+4*u*y+6*u*v-u*t-v*p

def word(kind,u,v):
 require(u>=1 and v>=1,'positive parameters')
 w=[] if u==1 else [int(c) for c in bin(u-1)[2:]]
 z=[] if v==1 else [int(c) for c in bin(v-1)[2:][::-1]]
 if kind=='N':
  require(u>=2 and v>=2,'N boundary')
  return [1 if b else 3 for b in w]+[2]+z
 if kind=='M':
  require(u>=2 and v>=2,'M boundary')
  return [1 if b else 3 for b in w]+[0,3]+z
 require(kind=='P','unknown family')
 return [2 if b else 3 for b in w]+[0,2,3]+[2*b for b in z]

def graph(kind,u,v):
 r=[2,0];ports=[0,1]
 for S in word(kind,u,v):r,ports=birth(r,ports,S,2)
 return r

def prefix_graph(u):
 require(u>=2,'strong one-apex prefix')
 r=[2,0];ports=[0,1]
 for c in bin(u-1)[2:]:r,ports=birth(r,ports,1 if c=='1' else 3,2)
 return r

def inverse(kind,h):
 """All parameter pairs u<=v in the proved finite search region."""
 require(h>=1,'positive target');out=[]
 lower=1 if kind=='P' else 2
 u=lower
 while True:
  if kind=='N':bound=(h-1)//3
  elif kind=='M':bound=(h-5)//6
  elif kind=='P':bound=(h-7)//6
  else:raise RuntimeError('unknown family')
  if u*u>bound:break
  for v in range(u,bound//u+1):
   if scalar(kind,u,v)==h:out.append((u,v))
  u+=1
 return out
