"""Verify finite spectral census and the new rank-two recurrence. Run from any cwd."""
from screen import *
import sympy as s,random
from ribbon import profile_count,ham_count_bruteforce,tournament_from_word
root=Path(__file__).parent
for w,d in [(2,9),(3,5)]:
 A,v,h=load(w);r=len(v);I=np.eye(r,dtype=np.int64);hits=[];tested=[0]
 def visit(B,word):
  if word:
   tested[0]+=1
   if singular(B-2*I):hits.append(word)
  if len(word)<d:
   for c,C in A.items():visit(C@B%1009,word+[c])
 visit(I,[])
 assert hits==json.loads((root/f'screen{w}_{d}.json').read_text())
 survivors=[]
 for word in hits:
  B=I
  for c in word:B=A[c]@B%1013
  if not singular(B-2*I,1013):continue
  Q=s.eye(r)
  for c in word:Q=s.Matrix(A[c])*Q
  if (Q-2*s.eye(r)).det()==0:
   survivors.append(word)
   assert w==2 and len(word)==8
   x=s.Symbol('lambda');assert Q.charpoly().as_expr()==s.expand(x**12*(x-265)*(x-2))
 assert len(survivors)==(16 if w==2 else 0)
 if w==2:
  bases=[(0,2,0,3,2,3,3,3),(0,2,1,3,3,2,3,2)]
  assert {tuple(t) for t in survivors}=={b[i:]+b[:i] for b in bases for i in range(8)}
 rng=random.Random(220+w)
 for _ in range(100):
  word=[rng.choice(list(A)) for _ in range(rng.randint(1,8))]
  y=v.copy()
  for c in word:y=A[c]@y
  ww=[frozenset(i+1 for i in range(w) if c>>i&1) for c in word]
  assert h@y==profile_count(w,ww)==ham_count_bruteforce(tournament_from_word(w,ww))
 print(w,'words',tested[0],'mod1009 candidates',len(hits),'exact eigen2 survivors',len(survivors),flush=True)
A,v,h=load(2)
for word in [[0,2,0,3,2,3,3,3],[0,2,1,3,3,2,3,2]]:
 B=s.eye(14)
 for c in word:B=s.Matrix(A[c])*B
 assert B.rank()==2 and B**3==267*B**2-530*B
 y=s.Matrix(v); f=s.Matrix([h.tolist()])
 for k in range(1,8):
  y=B*y;assert int((f*y)[0])==(63*265**k+200*2**k)//263
 print('recurrence checked',word,flush=True)
print('PASS')
