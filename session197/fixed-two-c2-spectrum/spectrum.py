"""Complete arithmetic inverse for C2 tournaments with feedback number <=2."""
from families import g,scalar,graph,word,inverse as family_inverse
from core import birth,require

def Q(u,v):
 require(u>=2 and v>=1,'Q parameter domain');D,p=g(u),u%2
 return g(v)*(D+3*u-p-2)-(v%2)*(u-1)+v*(2*D+4*u-p-1)
def q_word(u,v):
 require(u>=2 and v>=1,'Q parameter domain')
 return [2 if c=='1' else 3 for c in bin(u-1)[2:]]+[1,3]+([]if v==1 else[2*int(c)for c in bin(v-1)[2:][::-1]])
def from_word(w):
 T=[2,0]
 for S in w:T,_=birth(T,[0,1],S,2)
 return T
def inverse(h):
 require(h>=1,'positive target');out=[]
 if h==15:out.append(dict(family='regular5',word=[0,1,3]))
 for f in ('N','M','P'):
  for u,v in family_inverse(f,h):out.append(dict(family=f,u=u,v=v,word=word(f,u,v)))
 for u in range(2,(h-1)//6+1):
  for v in range(1,(h-2*u)//(4*u+1)+1):
   if Q(u,v)==h:out.append(dict(family='Q',u=u,v=v,word=q_word(u,v)))
 return out
