"""Independent nilpotent/commutative polynomial audit of Sections 1--3."""
import itertools
NV=13
zero=(0,)*NV
def constant(c):return {(0,zero):c} if c else {}
def variable(i):
 e=list(zero);e[i]=1;return {(0,tuple(e)):1}
def z(i):return {(1<<i,zero):1}
def add(*polys):
 out={}
 for p in polys:
  for e,c in p.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def scale(p,c):return {e:c*v for e,v in p.items() if c*v}
def mul(*polys):
 out=constant(1)
 for p in polys:
  new={}
  for (m,e),c in out.items():
   for (n,f),d in p.items():
    if m&n:continue
    key=(m|n,tuple(x+y for x,y in zip(e,f)));new[key]=new.get(key,0)+c*d
  out={e:c for e,c in new.items() if c}
 return out
def coeff(p,mask=7):return {(0,e):c for (m,e),c in p.items() if m==mask}
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
a,b,t,Ei,Ej,Si,Sj,Gij,Gji,li,lj,ri,rj=[variable(i) for i in range(NV)]
one=constant(1)
M=[[{} for _ in range(3)] for _ in range(3)]
for length in range(1,4):
 for word in itertools.permutations(range(3),length):
  term=mul(*[z(i) for i in word])
  for u,v in zip(word,word[1:]):term=mul(term,add(t,one if (v-u)%3==1 else {}))
  M[word[0]][word[-1]]=add(M[word[0]][word[-1]],term)
C=[add(*[M[v][u] for v in range(3)]) for u in range(3)];R=[add(*M[u]) for u in range(3)];H=add(*C)
need(add(M[0][1],scale(M[1][0],-1))==mul(z(0),z(1),add(one,scale(mul(z(2),add(one,scale(t,2))),-1))),'M skew')
need(not coeff(add(mul(C[0],C[1]),scale(mul(R[0],R[1]),-1))),'CC/RR coefficient')
need(coeff(add(mul(C[0],R[1]),scale(mul(R[0],C[1]),-1)))==constant(2),'CR/RC coefficient')
need(coeff(mul(H,add(M[0][1],scale(M[1][0],-1))))==one,'H times M skew')
def full(transpose=False):
 m=[[M[v][u] if transpose else M[u][v] for v in range(3)] for u in range(3)]
 c=[add(*[m[v][u] for v in range(3)]) for u in range(3)];r=[add(*m[u]) for u in range(3)]
 E={u:mul(a,c[u]) for u in range(3)};S={u:mul(b,r[u]) for u in range(3)};g={(u,v):m[u][v] for u in range(3) for v in range(3)}
 for u,e,s,l,rho in [(3,Ei,Si,li,ri),(4,Ej,Sj,lj,rj)]:
  E[u]=add(e,mul(a,H,rho));S[u]=add(s,mul(b,l,H))
  for v in range(3):g[u,v]=mul(l,c[v]);g[v,u]=mul(r[v],rho)
 g[3,4]=add(Gij,mul(li,H,rj));g[4,3]=add(Gji,mul(lj,H,ri))
 return E,S,g
def quartic(first,second,transpose=False):
 E,S,g=full(transpose);terms=[]
 for x,y in [first,first[::-1]]:
  for u,v in [second,second[::-1]]:terms.extend([mul(E[x],g[y,u],S[v]),mul(E[u],g[v,x],S[y])])
 return add(*terms)
actual=coeff(add(quartic((0,3),(1,4)),scale(quartic((0,3),(1,4),True),-1)))
expected=add(mul(scale(mul(a,b),2),add(Gij,scale(Gji,-1))),scale(mul(add(one,scale(t,2)),add(mul(Ei,Sj),scale(mul(Ej,Si),-1))),-1),mul(a,add(mul(ri,Sj),scale(mul(rj,Si),-1))),mul(b,add(mul(Ei,lj),scale(mul(Ej,li),-1))))
need(actual==expected,'complete formal disjoint kernel')
need(not add(quartic((0,1),(3,4)),scale(quartic((0,1),(3,4),True),-1)),'same-port formal zero before A extraction')
RESULT={'status':'PASS','formal_M_terms':sum(len(p) for row in M for p in row),'kernel_monomials':len(actual),'exact_CR_identity_checks':4,'same_port_zero':True,'exterior_variables':['alpha','beta','tau','E_i','E_j','S_i','S_j','G_ij','G_ji','ell_i','ell_j','rho_i','rho_j']}
if __name__=='__main__':
 import json;print(json.dumps(RESULT,sort_keys=True,indent=2))
