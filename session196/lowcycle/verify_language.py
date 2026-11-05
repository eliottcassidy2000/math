"""All-order finite equivalence certificate for the exact C=2 expression.
Compiles the independently specified expression with Thompson's construction,
then checks every reachable pair against the saturated-count Moore DFA.
Also checks state distinguishability and the rational word-count identity.
"""
from collections import deque
from apex_transfer import build_cycle_dfa,minimize_dfa

def A(x):return ('a',x)
def E():return ('e',)
def U(*xs):return ('u',)+xs
def C(*xs):return ('c',)+xs
def S(x):return ('s',x)
def word(s):return C(*(A(int(x)) for x in s))

P=U(E(),C(S(U(A(1),A(3))),A(1)))
Q=U(E(),C(S(U(A(2),A(3))),A(2)))
R=U(E(),C(S(A(3)),A(1),P),C(S(A(3)),A(2),Q))
plus3=C(A(3),S(A(3)))
plus0=C(A(0),S(A(0)))
L1=C(A(1),S(U(A(0),A(1))),U(C(U(word('20'),word('310'),C(word('32'),S(A(3)),A(1))),P),C(A(2),plus3,A(0))))
L2=C(A(2),S(U(A(0),A(2))),U(C(U(C(A(1),S(A(3)),A(2)),word('30')),Q),C(word('33'),S(A(3)),A(0))))
L3=C(A(3),U(C(plus0,U(C(A(2),R),C(plus3,A(2),Q))),C(word('10'),P),word('20'),C(A(2),S(A(3)),A(1),P)))
EXPRESSION=U(L1,L2,L3)

class NFA:
 def __init__(self):self.edges=[]
 def state(self):self.edges.append([]);return len(self.edges)-1
 def edge(self,a,letter,b):self.edges[a].append((letter,b))
 def build(self,r):
  start=self.state();end=self.state();kind=r[0]
  if kind=='a':self.edge(start,r[1],end)
  elif kind=='e':self.edge(start,None,end)
  elif kind=='u':
   for x in r[1:]:
    a,b=self.build(x);self.edge(start,None,a);self.edge(b,None,end)
  elif kind=='c':
   previous=start
   for x in r[1:]:
    a,b=self.build(x);self.edge(previous,None,a);previous=b
   self.edge(previous,None,end)
  elif kind=='s':
   a,b=self.build(r[1]);self.edge(start,None,end);self.edge(start,None,a);self.edge(b,None,a);self.edge(b,None,end)
  else:raise ValueError(kind)
  return start,end
 def close(self,states):
  ans=set(states);todo=list(ans)
  while todo:
   for letter,target in self.edges[todo.pop()]:
    if letter is None and target not in ans:ans.add(target);todo.append(target)
  return frozenset(ans)
 def step(self,states,letter):return self.close({target for state in states for x,target in self.edges[state] if x==letter})

def main():
 raw,rt,ro=build_cycle_dfa();partition,t,o=minimize_dfa(rt,ro)
 assert len(raw)==251 and len(t)==32
 import json
 from pathlib import Path
 stored=json.loads(Path(__file__).with_name("cycle_dfa.json").read_text())
 assert stored["initial_state"]==partition[0]
 assert stored["transitions"]==t and stored["outputs"]==o
 assert stored["raw_states"]==len(raw) and stored["states"]==len(t)
 assert all(o[partition[i]]==ro[i] for i in range(len(raw)))
 assert all(partition[rt[i][a]]==t[partition[i]][a] for i in range(len(raw)) for a in range(4))
 # Every pair of quotient states has a distinguishing continuation.
 distinguish={(i,j) for i in range(32) for j in range(i) if o[i]!=o[j]}
 while True:
  old=len(distinguish)
  for i in range(32):
   for j in range(i):
    if any(tuple(sorted((t[i][a],t[j][a]),reverse=True)) in distinguish for a in range(4)):
     distinguish.add((i,j))
  if old==len(distinguish):break
 assert len(distinguish)==32*31//2
 nfa=NFA();start,final=nfa.build(EXPRESSION)
 initial=(partition[0],nfa.close({start}));seen={initial};todo=deque([initial])
 while todo:
  ds,ns=todo.popleft()
  assert (o[ds]==2)==(final in ns),(ds,ns)
  for letter in range(4):
   nxt=(t[ds][letter],nfa.step(ns,letter))
   if nxt not in seen:seen.add(nxt);todo.append(nxt)
 print('all-order regex equivalence product states:',len(seen))
 # Rational expression, from the unambiguous first-exit decomposition.
 import sympy as sp
 z=sp.symbols('z');p=(1-z)/(1-2*z);r=1/(1-2*z)
 l1=z/(1-2*z)*((z**2+z**3+z**3/(1-z))*p+z**3/(1-z))
 l2=z/(1-2*z)*((z**2/(1-z)+z**2)*p+z**3/(1-z))
 l3=z*(z/(1-z)*(z*r+z**2/(1-z)*p)+z**2*p+z**2+z**2/(1-z)*p)
 expected=z**3*(7-14*z+8*z**2-5*z**3)/((1-z)*(1-2*z)**2)
 assert sp.cancel(l1+l2+l3-expected)==0
 partial=sp.Rational(5,4)*z**3+z**2/2+sp.Rational(47,16)*z+sp.Rational(61,16)-sp.Rational(5,32)/(1-2*z)+sp.Rational(11,32)/(1-2*z)**2-4/(1-z)
 assert sp.cancel(expected-partial)==0
 v=[0]*32;v[partition[0]]=1
 for m in range(41):
  got=sum(v[i] for i in range(32) if o[i]==2)
  if m>=4:assert got==((11*m+6)*2**m)//32-4,(m,got)
  elif m==3:assert got==7
  else:assert got==0
  w=[0]*32
  for i,n in enumerate(v):
   for j in t[i]:w[j]+=n
  v=w
 print('Moore quotient verified: 251 -> 32, all pairs distinguishable')
 print('rational word-count identity verified; exact count recurrence checked through length 40')

if __name__=='__main__':main()
