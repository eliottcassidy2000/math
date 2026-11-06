exec(open('session202-bridge-eventual/probe.py').read().split("D={'':I}")[0])
D={'':I}; mats={O:'O',Z:'Z',J:'J'};seen={I:''}; groups=[{} for _ in range(4)]
for k in range(13):
 hits=[]
 for w,m in D.items():
  a=mul(Q,m)
  for p,q in itertools.product(range(2),repeat=2):
   key=(a[0][0],a[0][2],a[0][3]+q*a[0][1],a[2][0],a[2][2],a[2][3]+q*a[2][1],a[3][0]+p*a[1][0],a[3][2]+p*a[1][2]);v=a[3][3]+p*a[1][3]+q*a[3][1]+p*q*a[1][1];g=groups[2*p+q]
   if key in g:
    w0,v0=g[key]
    if v0!=v:hits.append((p,q,w0,w,v0,v))
   else:g[key]=(w,v)
 print(k,len(D),'delta_constants',hits[:15],flush=True)
 nxt={}
 for w,m in D.items():
  for t,letter in mats.items():
   n=mul(t,m)
   if n not in seen:seen[n]=w+letter;nxt[w+letter]=n
 D=nxt
