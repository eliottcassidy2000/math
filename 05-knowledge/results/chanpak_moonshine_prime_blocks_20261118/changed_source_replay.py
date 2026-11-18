from unit_gain_replay import K,state,H,forward,bits
L=(16,9,3,5,14);R=(24,17,3,6,12)
assert state(L,1)==(15,3,19,5)
assert state(R,0)==(15,3,21,5)
for j in range(1,11):
 M=16**j;r=4*j+2;B=4*M;u=(7*M-2)//5
 C=(247*M+28-10*j)//5
 assert C%2==0
 V=C//2;wl=B+4*B*V;wr=u+4*B*V
 hl=15+19*wl-5*K(wl)-3*(wl&1)
 hr=15+21*wr-5*K(wr)-3*(wr&1)
 assert hl==hr and hl%5==(1-2*j)%5
 if j<=5:print('periodic source return',j,hl,hl%5)

rev=lambda x,n:int(f'{x:0{n}b}'[::-1],2)
def bank(r,A):
 B=1<<r;s=(4*B-A)//5
 assert A%2 and A%5==1 and 5*s==4*B-A
 rows=[]
 for x in range(2*B-s):
  for bit in (0,1):
   a=rev(x,r+1)+2*B*bit;b=rev(x+s,r+1)+2*B*bit
   C=19*a-21*b-5*(K(a)-K(b))-3*(a&1)+3*(b&1)
   if C<0 or C%(2*A):continue
   V=C//(2*A);wl=a+4*B*V;wr=b+4*B*V
   h=15+19*wl-5*K(wl)-3*(wl&1)
   assert h==15+21*wr-5*K(wr)-3*(wr&1)
   rows.append(h)
   if r==2:
    assert H(forward(L,1,bits(wl)))==h
    assert H(forward(R,0,bits(wr)))==h
 return set(rows)
bank(2,1)
old=bank(6,1)
for A in (1,11,21,31,41):
 S=bank(6,A)
 print('response',A,'distinct_centers',len(S),'outside_unit_bank',len(S-old),'first',sorted(S)[:3])
print('changed physical sources and general response: PASS')
