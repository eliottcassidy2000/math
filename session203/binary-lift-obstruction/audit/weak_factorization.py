import pathlib,sys,json,math
P=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P))
from counting import nine,labelled,relabel

def sub(a,ids):return tuple(sum(1<<j for j,v in enumerate(ids) if a[u]>>v&1) for u in ids)
def reach(a,i):
 seen=1<<i
 while True:
  q=seen
  for j in range(len(a)):
   if seen>>j&1:q|=a[j]
  if q==seen:return q
  seen=q
ct=weak=0
for n in range(1,6):
 for a in labelled(n):
  p=nine(a);h,s,e,L,R,d,C,b,T=p
  if not s:continue
  ct+=1;rs=[reach(a,i) for i in range(n)];A=[i for i in range(n) if rs[i]&1];rest=set(range(n))-set(A);q=1
  while rest:
   i=min(rest);D={j for j in rest if rs[i]>>j&1 and rs[j]>>i&1};q*=nine(sub(a,sorted(D)))[0];rest-=D
  aa=nine(sub(a,A));hh,ss,ee,ll,*_=aa
  assert h==q*hh and s==q*ss and L==q*ll and h+L-s==q*(hh+ll-ss)
  assert len(A)<=hh and q*len(A)<=h
  assert h+L-s<=2*q*math.factorial(len(A))<=2*math.factorial(h)
  if len(A)<n:assert e==0;weak+=1
print(json.dumps(dict(status='PASS',scope='All positive-start root-zero labeled tournaments through n5',positive_roots=ct,weak_roots=weak,checks='SCC h/s/L/B factorization; root-component order; factorial inequality; weak e=0'),indent=2))
