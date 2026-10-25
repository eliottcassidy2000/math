class Field:
 def __init__(self,n,f):
  self.n=n;self.f=f;self.q=3**n
  self.d=[tuple((i//3**k)%3 for k in range(n)) for i in range(self.q)]
 def enc(self,t):return sum(a%3*3**k for k,a in enumerate(t))
 def add(self,a,b):return self.enc(x+y for x,y in zip(self.d[a],self.d[b]))
 def neg(self,a):return self.enc(-x for x in self.d[a])
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):
  z=[0]*(2*self.n-1)
  for i,x in enumerate(self.d[a]):
   for j,y in enumerate(self.d[b]):z[i+j]+=x*y
  for k in range(len(z)-1,self.n-1,-1):
   for j in range(self.n):z[k-self.n+j]-=z[k]*self.f[j]
  return self.enc(z[:self.n])
