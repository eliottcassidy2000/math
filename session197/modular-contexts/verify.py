#!/usr/bin/env python3
"""Exact finite-context congruence controls. Writes only stdout."""
from core import *
from collections import deque
from itertools import product
from functools import lru_cache
import hashlib

def Epoly(f):return mul([0,1],add(mul([-2,0,1],f),mul([1,-2],shift(f,-1)),mul([1,2],shift(f,1))))
def weight(k,odd):return math.factorial(k)*gamma(k,odd)
def cutoff(M,odd):
 k=0
 while weight(k,odd)%M:k+=1
 return k
def Rcut(M):
 r=1
 while math.factorial(r-1)*math.factorial(r)%M:r+=1
 return r
def matrix(M,odd,op):
 outodd=odd if op=='D'else not odd;n,m=cutoff(M,odd),cutoff(M,outodd);A=[[0]*n for _ in range(m)]
 for k in range(n):
  f=(D if op=='D'else Epoly)(basis(k,odd));u=coords(f)
  for j in range(m):
   z=weight(j,outodd)*(u[j]if j<len(u)else 0);need(z%weight(k,odd)==0,'integral central transition');A[j][k]=z//weight(k,odd)%M
 return A

def wcoords(f,R):
 f=newton(f);return tuple((math.factorial(i-1)*f[i]if i<len(f)else 0)for i in range(1,R))
def acoords(f,M):
 odd=(len(f)-1)%2;u=coords(f);return(bool(odd),tuple(weight(k,odd)*(u[k]if k<len(u)else 0)%M for k in range(cutoff(M,odd))))
def substitute(q,r,t):
 vs=[(v,0)for v in range(len(q))if v!=r]+[(r,i)for i in range(len(t))];a=[0]*len(vs)
 for i,(v,k)in enumerate(vs):
  for j,(z,l)in enumerate(vs):
   if i!=j and(t[k]>>l&1 if v==z else q[v]>>z&1):a[i]|=1<<j
 return tuple(a)
@lru_cache(None)
def literal_profile(a):return profile(a)

def compile_newton(q,root,R):
 # operator on positive monomials, computed from actual transitive substitutions
 images={s:literal_profile(substitute(q,root,trans(s)))for s in range(1,R)}
 A=[[0]*(R-1)for _ in range(R-1)];fall=[1]
 for r in range(1,R):
  fall=mul(fall,[-(r-1),1]);out=add(*(scale(images[s],c)for s,c in enumerate(fall)if c));ws=wcoords(out,R);den=math.factorial(r-1)*math.factorial(r)
  for i,z in enumerate(ws):need(z%den==0,'integer Newton transition');A[i][r-1]=z//den
 return A

def finite_graph(M):
 n=cutoff(M,True);need(n==cutoff(M,False),'equal prime cutoffs');AA={(o,c):matrix(M,o,c)for o in(True,False)for c in'DE'}
 states=[(o,v)for o in(True,False)for v in product(range(M),repeat=n)]
 def step(s,c):
  o,v=s;return(o if c=='D'else not o,tuple(sum(x*y for x,y in zip(row,v))%M for row in AA[o,c]))
 start=(True,(1,)+(0,)*(n-1));words={start:''};queue=deque([start])
 while queue:
  s=queue.popleft()
  for c in'DE':
   z=step(s,c)
   if z not in words:words[z]=words[s]+c;queue.append(z)
 target=48 if M==5 else 686;need(len(words)==target,'complete prime orbit')
 rev={s:[]for s in states}
 for s in states:
  for c in'DE':rev[step(s,c)].append(s)
 back={start:0};queue=deque([start])
 while queue:
  s=queue.popleft()
  for t in rev[s]:
   if t not in back:back[t]=back[s]+1;queue.append(t)
 need(len(back)==2*(M**n-1)and all(any(s[1])for s in back),'entire nonzero communicating part')
 if M==5:
  for c in'DE':need(len({step(s,c)for s in states})==len(states),'mod5 permutations')
 for o in(True,False):need(all(step((o,(0,)*n),c)[1]==(0,)*n for c in'DE'),'zero absorbing')
 maxlen,maxcost=(4,10)if M==5 else(5,13);catalog=[]
 for length in range(1,maxlen+1):
  for word in product('DE',repeat=length):
   cost=sum(2 if c=='D'else 3 for c in word)
   if cost<=maxcost:catalog.append((''.join(word),cost))
 checked=0;worst_steps=worst_cost=0;witness=None
 for s in states:
  if not any(s[1]):continue
  got={}
  for word,cost in catalog:
   z=s
   for c in word:z=step(z,c)
   got.setdefault(z[1][0],[]).append((len(word),cost,word))
  need(len(got)==M,'joint positive word/order steering bound')
  for r,ws in got.items():
   mn=min(x[0]for x in ws);mc=min(x[1]for x in ws);worst_steps=max(worst_steps,mn)
   if mc>worst_cost:worst_cost=mc;witness=[s,r,min(ws,key=lambda x:(x[1],x[0],x[2]))]
   checked+=1
 # Prime-function descriptions are injective within parity and meet only at zero.
 funcs=[]
 for odd in(True,False):
  fs=set()
  for v in product(range(M),repeat=n):
   f=add(*(scale(basis(k,odd),a*pow(weight(k,odd),-1,M)%M)for k,a in enumerate(v)))
   fs.add(tuple(val(f,x)%M for x in range(M)))
  need(len(fs)==M**n,'same parity injective function map');funcs.append(fs)
 need(len(funcs[0]&funcs[1])==1 and len(funcs[0]|funcs[1])==2*M**n-1,'global exact quotient count')
 # A nonempty return word yields a strong representative of the initial state.
 returnword=None
 for c in'DE':
  z=step(start,c);q=deque([(z,c)]);seen={z}
  while q and returnword is None:
   z,w=q.popleft()
   if z==start:returnword=w;break
   for t in'DE':
    y=step(z,t)
    if y not in seen:seen.add(y);q.append((y,w+t))
  if returnword:break
 need(returnword is not None,'nonempty initial-state return')
 digest=hashlib.sha256(json.dumps(sorted((int(o),v,w)for(o,v),w in words.items()),separators=(',',':')).encode()).hexdigest()
 return {'M':M,'matrices':{('odd'if o else'even')+c:A for(o,c),A in AA.items()},'tagged_orbit':len(words),'maximum_shortest_word':max(map(len,words.values())),'nonzero_return_states':len(back),'global_classes':len(funcs[0]|funcs[1]),'positive_steering_targets':checked,'joint_length_bound':maxlen,'joint_added_order_bound':maxcost,'observed_max_min_length':worst_steps,'observed_max_min_order':worst_cost,'worst_order_example':witness,'nonempty_return_word':returnword,'state_word_sha256':digest}

def run():
 moduli=(3,4,5,7,9);contexts=physical=closure=0
 inputs=[tournament(n,c)for n in range(1,5)for c in range(1<<(n*(n-1)//2))]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   q=tournament(n,code)
   for root in range(n):
    contexts+=1;A=compile_newton(q,root,max(Rcut(M)for M in moduli))
    for t in inputs:
     f=literal_profile(t);g=literal_profile(substitute(q,root,t));physical+=1
     for M in moduli:
      R=Rcut(M);v=wcoords(f,R);out=wcoords(g,R)
      predicted=tuple(sum(A[i][j]*v[j]for j in range(R-1))%M for i in range(R-1))
      need(predicted==tuple(x%M for x in out),'literal arbitrary context modular closure');closure+=1
 wrapper_controls=0
 for t in inputs:
  f=literal_profile(t)
  for op,graph in(('D',W(t)),('E',E(t))):
   pf=(D if op=='D'else Epoly)(f);need(literal_profile(graph)==pf and strong(graph),'literal strong wrapper polynomial')
   for M in (3,4,5,7,9,16,25,49):
    o,v=acoords(f,M);A=matrix(M,o,op);z=(o if op=='D'else not o,tuple(sum(a*b for a,b in zip(row,v))%M for row in A))
    need(z==acoords(pf,M),'composite and parity-changing central state');wrapper_controls+=1
 # Cutoffs and high-coordinate vanishing, including signed profiles of high order.
 algebra=0
 for p in(3,5,7,11,13,17,19):
  for odd in(False,True):
   need(cutoff(p,odd)==(p-1)//2,'prime cutoff');need(cutoff(p*p,odd)==(p if odd else p-1),'prime-square cutoff')
   for e in range(1,8):need(cutoff(p**e,odd)<=p*((e+2)//3),'prime-power bound');algebra+=1
 for M in(2,3,4,5,7,8,9,16,25,27,49):
  for odd in(False,True):
   for k in range(cutoff(M,odd),cutoff(M,odd)+8):need(weight(k,odd)%M==0,'high coordinate zero');algebra+=1
 # Every input order: exact matrix characteristic and explicit mod7 periodicity.
 rec=[]
 for odd in(False,True):
  A=[[0]*3 for _ in range(3)]
  for k in range(3):
   A[k][k]=k*k+6*k+(3 if odd else 5)
   if k<2:A[k][k+1]=1
   if k:A[k][k-1]=beta(k,odd)
  def mm(X,Y):return[[sum(X[i][k]*Y[k][j]for k in range(3))for j in range(3)]for i in range(3)]
  I=[[int(i==j)for j in range(3)]for i in range(3)];A2=mm(A,A);A3=mm(A2,A);co=(336,-231,32)if odd else(708,-345,38)
  need(A3==[[co[2]*A2[i][j]+co[1]*A[i][j]+co[0]*I[i][j]for j in range(3)]for i in range(3)],'exact cubic characteristic')
  if not odd:
   P=I
   for _ in range(16):P=[[z%7 for z in row]for row in mm(A,P)]
   need(P==I,'even matrix period16')
  else:need(all((A3[i][j]-4*A2[i][j])%7==0 for i in range(3)for j in range(3)),'odd order-three recurrence')
  q=[0,1]if odd else[0,0,1];hs=[]
  for j in range(34):hs.append(val(q,1));q=D(q)
  if not odd:
   expected=[1,5,2,4,6,5,0,3,0,1,6,2,2,1,1,3];need([h%7 for h in hs[:16]]==expected,'literal vacuum residues');need(all((h%7==0)==(j%16 in(6,8))for j,h in enumerate(hs)),'zero residue classes');need(hs[6]==7*89*3851 and hs[8]==7*241*367*3709,'composite early targets')
  rec.append({'odd':odd,'matrix':A,'first16_mod7':[h%7 for h in hs[:16]]})
 cyc5=tuple(sum(1<<((i+j)%5)for j in(1,2))for i in range(5));f=literal_profile(cyc5);need(strong(cyc5)and all(x%5==0 for x in wcoords(f,5)),'actual zero class')
 graphs=[finite_graph(5),finite_graph(7)]
 print(json.dumps({'status':'PASS','arbitrary_rooted_contexts':contexts,'actual_context_profiles':physical,'Newton_modular_closure_controls':closure,'central_wrapper_controls':wrapper_controls,'cutoff_controls':algebra,'zero_class_C5':{'rows':cyc5,'profile':f},'uniform_recurrences':rec,'finite_quotients':graphs},sort_keys=True,indent=2))
if __name__=='__main__':run()
