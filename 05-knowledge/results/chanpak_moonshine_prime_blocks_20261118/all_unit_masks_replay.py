from unit_gain_replay import K,L,R,forward,H,bits
from collections import Counter

def rev(x,n):return int(f'{x:0{n}b}'[::-1],2)
for r in (2,6,10,14):
 B=1<<r;n=r+1;s=(4*B-1)//5;rows=[];bins=Counter()
 for x in range(2*B-s):
  a,b=rev(x,n),rev(x+s,n)
  assert 4*B+5*(rev(a,n)-rev(b,n))==1
  for eps in (0,1):
   aa=a+2*B*eps;bb=b+2*B*eps
   C=23*aa-25*bb-5*(K(aa)-K(bb))-(aa&1)+(bb&1)
   assert C%2==0 and (C>0)==bool(x&1)
   if C<0:continue
   V=C//2;wl=aa+4*B*V;wr=bb+4*B*V
   hl=9+23*wl-5*K(wl)-(wl&1)
   hr=9+25*wr-5*K(wr)-(wr&1)
   assert hl==hr and 100*B*B<hl<2200*B*B
   rows.append(hl);bins[V]+=1
   if r==2:
    assert H(forward(L,3,bits(wl)))==hl
    assert H(forward(R,3,bits(wr)))==hr
 assert len(rows)==(6*B-4)//5
 print({'r':r,'pairs':len(rows),'distinct_heights':len(set(rows)),'max_same_V':max(bins.values())})
print('all unit masks: PASS')
