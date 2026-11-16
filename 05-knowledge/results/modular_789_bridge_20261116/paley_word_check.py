from itertools import product,combinations
N=7
Q={1,2,4}
B=lambda *r:{frozenset((a+x)%N for x in r) for a in range(N)}
S1,S2=B(0,1,3),B(0,1,5)
def compose(g,f):return tuple(g[x] for x in f)
for b in (-1,1):
  F=lambda k: tuple(pow(2,-k,N)*(3*x+b)%N for x in range(N))
  from itertools import permutations
  for k in range(1,8):
    f=F(k)
    assert {frozenset(f[x] for x in row) for row in S1}==S2
    assert all((f[y]-f[x])%N not in Q for x,y in permutations(range(7),2) if (y-x)%7 in Q)
  P={tuple(range(7))}
  while True:
    new=P|{compose(F(k),w) for w in P for k in range(1,4)}
    if new==P:break
    P=new
  assert len(P)==42
  for k in product(range(1,4),repeat=4):
    r,s=1,1
    for t in reversed(k):
      s+=t
      r=((2**t*r-b)*pow(3,-1,2**s))%(2**s)
    for x in range(7):
      n=r+2**s*((x-r)*pow(2**s,-1,7)%7+7)
      a=n%7
      for t in k:
        z=3*n+b
        assert z//2**t%2==1 and z%2**t==0
        n=z//2**t
        a=F(t)[a]
      assert n%7==a
print('mod-seven signed affine group, Fano sheet exchange, and finite native CRT chart: PASS')
