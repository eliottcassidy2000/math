import itertools,json,math
import sympy as s
class Field:
 def __init__(self,e):
  self.e=e; self.q=3**e; x=s.Symbol('x')
  for a in itertools.product(range(3),repeat=e):
   if a[0] and s.Poly(x**e+sum(a[i]*x**i for i in range(e)),x, modulus=3).is_irreducible:
    self.mod=list(a);break
  self.digits=[[(v//3**i)%3 for i in range(e)] for v in range(self.q)]
  squares={self.mul(x,x) for x in range(1,self.q)}
  self.chi=[0]+[1 if x in squares else -1 for x in range(1,self.q)]
  assert len(squares)==(self.q-1)//2
 def sub(self,x,y):return sum(((a-b)%3)*3**i for i,(a,b) in enumerate(zip(self.digits[x],self.digits[y])))
 def mul(self,x,y):
  e=self.e;c=[0]*(2*e-1)
  for i,a in enumerate(self.digits[x]):
   for j,b in enumerate(self.digits[y]):c[i+j]=(c[i+j]+a*b)%3
  for k in range(2*e-2,e-1,-1):
   for j,a in enumerate(self.mod):c[k-e+j]=(c[k-e+j]-c[k]*a)%3
  return sum(c[i]*3**i for i in range(e))
 def edge(self,x,y):return self.chi[self.sub(y,x)]==1
 def embed(self,A):
  xs=[]
  for v in range(len(A)):
   for x in range(self.q):
    if x not in xs and all(self.edge(xs[i],x)==bool(A[i]>>v&1) for i in range(v)):
     xs.append(x);break
   else:raise AssertionError('no extension')
  return xs

def H(A):
 n=len(A);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1<<n):
  for v in range(n):
   for w in range(n):
    if not mask>>w&1 and A[v]>>w&1:d[mask|1<<w][w]+=d[mask][v]
 return sum(d[-1])
def matrices(n):
 pairs=list(itertools.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  A=[0]*n
  for k,(i,j) in enumerate(pairs):
   u,v=(i,j) if code>>k&1 else (j,i);A[u]|=1<<v
  yield A
f=Field(7); tests=0;cert={}
for n in range(1,5):
 for A in matrices(n):
  xs=f.embed(A); h=H(A)
  raw=sum(math.prod(1+f.chi[f.sub(xs[p[j+1]],xs[p[j]])] for j in range(n-1)) for p in itertools.permutations(range(n)))
  assert raw==2**(n-1)*h
  tests+=1;cert.setdefault(h,{'matrix':A,'field_vertices':xs})
g=Field(3); patterns=0; minimum=100
for a,b in itertools.combinations(range(g.q),2):
 counts=[0]*4
 for x in range(g.q):
  if x not in (a,b):counts[(g.chi[g.sub(x,a)]==1)+2*(g.chi[g.sub(x,b)]==1)]+=1
 assert min(counts)>0
 minimum=min(minimum,min(counts));patterns+=4
print(json.dumps({'field_q':f.q,'irreducible_coefficients_low_to_high':f.mod+[1],'literal_embeddings_and_integer_character_equations':tests,'certificates':cert,'F27_two_anchor_patterns':patterns,'minimum_pattern_population':minimum,'note':'Finite checks support identities, not universal target coverage.'},indent=2))
