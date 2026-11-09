import itertools,json,hashlib

def need(x,msg):
 if not x: raise RuntimeError(msg)
def tours(n):
 edges=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(edges)):
  a=[0]*n
  for k,(i,j) in enumerate(edges):
   if mask&(1<<k):a[i]|=1<<j
   else:a[j]|=1<<i
  yield a

def paths(a):
 n=len(a); size=1<<n
 ends=[[0]*n for _ in range(size)];h=[0]*size;h[0]=1
 for v in range(n):ends[1<<v][v]=1
 for s in range(1,size):
  h[s]=sum(ends[s])
  for v,c in enumerate(ends[s]):
   if c:
    bits=a[v]&~s
    while bits:
     b=bits&-bits;bits-=b
     ends[s|b][b.bit_length()-1]+=c
 return h,ends

def state(a,r=0):
 h,e=paths(a);n=len(a);full=(1<<n)-1
 # Reverse paths ending at r are exactly paths starting at r.
 rev=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)]
 _,back=paths(rev)
 return (h[full],e[full][r],sum(back[s][r]*h[full^s] for s in range(1,full+1)),h[full^(1<<r)])
def reach(a,r):
 s=1<<r
 while True:
  z=s
  for i,row in enumerate(a):
   if s>>i&1:z|=row
  if z==s:return s
  s=z

def strong(a):return all(reach(a,r)==(1<<len(a))-1 for r in range(len(a)))
def join(a,b):
 n=len(a);return [r|(((1<<len(b))-1)<<n) for r in a]+[r<<n for r in b]
def graft(a,r,b):
 n=len(a);extra=((1<<len(b))-1)<<n
 return [row|(extra if i!=r else 0) for i,row in enumerate(a)]+[(row<<n)|(1<<r) for row in b]
def brute(a,r):
 n=len(a);full=(1<<n)-1;h=[0]*(full+1);s=[0]*(full+1);e=[0]*(full+1);h[0]=1
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(a[u]>>v&1 for u,v in zip(p,p[1:])):
    mask=sum(1<<v for v in p);h[mask]+=1;s[mask]+=p[0]==r;e[mask]+=p[-1]==r
 return h[full],e[full],sum(s[m]*h[full^m] for m in range(full+1)),h[full^(1<<r)]

def main():
 starts={};atoms={};gadgets={};seen=0;brutes=0;prefix=0;small=[]
 for n in range(1,7):
  for a in tours(n):
   seen+=1;h,_=paths(a);H=h[-1];Q=sum(h[m]*h[(1<<n)-1-m] for m in range(1<<n))
   need(Q>=2*H and Q>=1<<n,'gadget floors')
   if n<=4:
    gadgets.setdefault((H,Q),a);small.append(a)
   if n<=5:
    x=state(a);need(x==brute(a,0),'independent permutation path counts');brutes+=1
    suffixmask=reach(a,0)
    suffix=[v for v in range(n) if suffixmask>>v&1];pre=[v for v in range(n) if not suffixmask>>v&1]
    sub=lambda vs:[sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs]
    if pre:
     factor=paths(sub(pre))[0][-1];y=state(sub(suffix),suffix.index(0));need(x==tuple(factor*t for t in y),'SCC root suffix factor');prefix+=1
   if not strong(a):continue
   if Q<=31:atoms.setdefault((H,Q),a)
   if H<=20:
    x=state(a)
    if x[2]<=31:starts.setdefault(x,a)
 states=dict(starts);stack=list(starts)
 while stack:
  x=stack.pop();D,E,U,d=x
  for (h,q),b in atoms.items():
   y=(D*h,0,U*q,d*h)
   if y[0]<=20 and y[2]<=31 and y not in states:
    states[y]=join(states[x],b);stack.append(y)
 for x,a in states.items():need(state(a)==x and reach(a,0)==(1<<len(a))-1,'source transfer and reach')
 out={};pairs=0
 for x,a in states.items():
  D,E,U,d=x
  if not a[0]:continue
  for (H,Q),b in gadgets.items():
   N=H*(D-E+U)+(Q-2*H)*(U-d)
   if N>31:continue
   g=graft(a,0,b);need(paths(g)[0][-1]==N and strong(g),'literal output');pairs+=1
   need(H+Q<=N and H<=N//3 and U<=N//H-1 and 3*H*D<=2*N-2*Q+3*H and D<N,'descent')
   out.setdefault(N,dict(source=a,state=x,gadget=b,H=H,Q=Q))
 # All arbitrary sources and arbitrary gadgets through order four, including sinks.
 controls=0;sinks=0
 for a in small:
  D,E,U,d=state(a)
  for b in small:
   hp,_=paths(b);H=hp[-1];Q=sum(hp[m]*hp[(1<<len(b))-1-m] for m in range(1<<len(b)))
   N=paths(graft(a,0,b))[0][-1]
   need(N==H*(D-E+U)+(Q-2*H)*(U-d),'unfiltered literal scalar law')
   if a[0]:need(H+Q<=N and 3*H*D<=2*N-2*Q+3*H and H<N and D<N,'arbitrary nonsink descent')
   else:need((D,E,U,d)==(d,d,d,d) and N==H*d,'sink boundary');sinks+=1
   controls+=1
 need(len(states)==77,'state total');need(sorted(out)==[3,5,9,11,13,15,17,19,23,25,27,29,31],'image')
 print(json.dumps(dict(status='PASS',labelled_tournaments=seen,root_convention='Root fixed at 0, all labelings',permutation_controls=brutes,SCC_prefix_controls=prefix,initial_states=len(starts),atom_pairs=sorted(atoms),gadget_pairs=sorted(gadgets),source_states=len(states),literal_low_target_pairs=pairs,unfiltered_grafts=controls,sink_controls=sinks,image=sorted(out),witnesses=out),sort_keys=True,indent=2))
if __name__=='__main__':main()
