from unit_gain_replay import K

def data(j):
 r=4*j+2;B=2**r;u=(7*16**j-2)//5
 C=(299*16**j+36-10*j)//5
 z=[-(u>>i&1) for i in range(r)]+[1]
 free=[i for i in range(r) if not (u>>i&1)]+[r+1]
 def w(i,k):i,k=sorted((i,k));return 2**(k-i-1)+(k==i+1)
 c=[2*2**i+5*sum(z[k]*w(i,k) for k in range(r+1) if i!=k) for i in free]
 return r,B,u,C,free,c

def bank(j):
 r,B,u,C,free,c=data(j)
 for mask in range(2**len(free)):
  P=sum(2**i for k,i in enumerate(free) if mask>>k&1)
  V=(C-sum(v for k,v in enumerate(c) if mask>>k&1))//2
  wl=B+P+4*B*V;wr=u+P+4*B*V
  hl=9+23*wl-5*K(wl)-(wl&1);hr=9+25*wr-5*K(wr)-(wr&1)
  assert hl==hr and 400*B*B<=hl<=1000*B*B
  yield hl,wl,wr

def select(N):
 j=1
 while 400*(4*16**j)**2<=N:
  if N<=1000*(4*16**j)**2:
   for item in bank(j):
    if item[0]==N:return item
   return None
  j+=1
 return None

for j in range(1,11):
 r,B,u,C,free,c=data(j)
 expected=[4,6,6]+[v for k in range(j-1) for v in (104*16**k-2,208*16**k-2)]+[104*16**(j-1)-2]
 assert sorted(c)==sorted(expected)
 assert sum(c)==(39*16**j-14-20*j)//5
 assert C-sum(c)==52*16**j+10+2*j
 if j<=5:rows=list(bank(j));print(j,len(rows),len(set(t[0] for t in rows)),min(t[0] for t in rows),max(t[0] for t in rows))
 if j<=4:
  for row in rows:assert select(row[0]) is not None
print('complete response family and target inverse: PASS')
