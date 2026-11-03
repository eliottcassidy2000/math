from pathlib import Path
exec((Path(__file__).resolve().parent / 'cone.py').read_text().split('for n in range(1,6):')[0])
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for code in range(1<<len(pairs)):
  rows=[0]*n
  for i,(u,v) in enumerate(pairs): rows[(u,v)[not bool(code>>i&1)]] |= 1<<(v,u)[not bool(code>>i&1)]
  h,e=dpc(rows)
  for p in range(n):
   aa=add_a(add_a(rows,p),n);haa,_=dpc(aa); xy=(haa-h)//2
   cnt=0
   for s in range(1,1<<n):
    if s>>p&1:continue
    l=[i for i in range(n) if s>>i&1]; r=[i for i in range(n) if not s>>i&1]
    left=[sum((rows[x]>>y&1)<<j for j,y in enumerate(l)) for x in l]
    right=[sum((rows[x]>>y&1)<<j for j,y in enumerate(r)) for x in r]
    rr=[sum((rows[y]>>x&1)<<j for j,y in enumerate(r)) for x in r]
    cnt+=dpc(left)[0]*dpc(rr)[1][r.index(p)]
   start=dpc([sum((rows[j]>>i&1)<<j for j in range(n)) for i in range(n)])[1][p]
   if xy!=start+cnt: raise AssertionError((n,code,p,xy,start,cnt))
 print('n',n,'verified')
