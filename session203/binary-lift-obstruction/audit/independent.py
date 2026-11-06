import sys,pathlib,itertools,math,json
P=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P))
from compiler import move,base,SEED
ops=('O','Z','J','K1','K2')
states={SEED};layer={SEED}
for depth in range(5):
 layer={move(p,o) for p in layer for o in ops};states|=layer
checks=0;accepted=0;small_defects=0
for p in states:
 B=base(p)
 for c in range(1,5):
  H=2*math.factorial(c)
  for o in ('J','K1','K2'):
   D=base(move(p,o))-2*B
   if D<=c:assert B<=H;small_defects+=1
  if B<=H:continue
  for length in range(4):
   for w in itertools.product(ops,repeat=length):
    q=p
    for o in w:q=move(q,o)
    checks+=1
    if abs(base(q)-2*B)<=c:
     accepted+=1;assert w and w[0]=='Z' and all(o=='O' for o in w[1:]),(p,c,w,q)
     if len(w)>1:
      assert base(move(q,'Z'))>2*base(q)+c
for c in range(1,5):
 H=2*math.factorial(c)
 for p in states:
  B=base(p)
  if B<=H:continue
  q=move(p,'Z')
  if base(q)-2*B>c:continue
  h,s,e,L,R,d,C,b,T=q;delta=s+(R-h+e-d)+(b-L-d)
  assert d==e and 1<=delta<=c
  z=q
  for n in range(1,6):
   z=move(z,'Z');assert base(z)==2**n*(base(q)+delta)-delta
print(json.dumps(dict(status='PASS',states=len(states),scope='All distinct transfer states from words of length<=5; c=1..4; all candidate lift words of length<=3; first-lift-normalized Z rays',candidate_words=checks,bounded_lifts=accepted,bounded_nonZ_defects=small_defects),indent=2))
