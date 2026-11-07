"""Independently check retained recipes by the closed binary formula and literal DP controls."""
import json,pathlib,struct,random,itertools
D=pathlib.Path(__file__).resolve().parent
CORES=json.loads((D/'CORES.json').read_text())
def need(b,m):
 if not b:raise RuntimeError(m)
def tables(a):
 n=len(a);N=1<<n;h=[0]*N;h[0]=1;s=[[0]*n for _ in range(N)];e=[[0]*n for _ in range(N)]
 for v in range(n):s[1<<v][v]=e[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if not(mask>>v&1):continue
   prev=mask^(1<<v)
   if prev:
    s[mask][v]=sum(s[prev][u] for u in range(n) if a[v]>>u&1)
    e[mask][v]=sum(e[prev][u] for u in range(n) if a[u]>>v&1)
  h[mask]=sum(e[mask])
 return h,s,e
def hp(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,N):
  for v,x in enumerate(dp[m]):
   if not x:continue
   z=a[v]&~m
   while z:
    b=z&-z;dp[m|b][b.bit_length()-1]+=x;z-=b
 return sum(dp[-1])
C=[]
for core in CORES:
 n=len(core);h,s,e=tables(core);N=(1<<n)-1
 for cut in range(1<<n):
  A=sum(h[N^m]*sum(s[m][v] for v in range(n) if cut>>v&1) for m in range(1,N+1))
  G=sum(h[N^m]*sum(e[m][v] for v in range(n) if cut>>v&1) for m in range(1,N+1))
  T=sum(sum(e[m][v] for v in range(n) if cut>>v&1)*sum(s[N^m][v] for v in range(n) if cut>>v&1) for m in range(1,N))
  q=sum(h[m]*h[N^m] for m in range(N+1));S=sum(s[N][v] for v in range(n) if cut>>v&1);E=sum(e[N][v] for v in range(n) if cut>>v&1)
  C.append((h[N],q,A,G,T,S,E))
need(C==[tuple(map(int,l.split())) for l in (D/'CORE_STATES.txt').read_text().splitlines()],'core table')
def ff(p):
 # Independent pair-of-selected-positions definition, not the producer recurrence.
 pos=[j for j in range(p.bit_length()) if p>>j&1]
 return sum((1<<(j-i-1))+(j==i+1) for k,i in enumerate(pos) for j in pos[k+1:])
def decode(code):
 w=code-1;return w&127,(w>>7)&31,(w>>12)&31,(w>>17)&65535,(w>>33)&65535
def formula(c,a,b,p,r):
 h,q,A,G,T,S,E=C[c];rev=int(f'{r:0{b}b}'[::-1],2) if b else 0
 return (q*p+A)*((1<<(b+1))-rev)+h*(2*r+1-ff(r))-p*G-h*ff(p)-T-(p&1)*S-E*((r>>(b-1))&1 if b else 0)
def corecut(c):
 return (CORES[0],c) if c<16 else (CORES[1+(c-16)//32],(c-16)%32)
def graph(c,a,b,p,r):
 core,cut=corecut(c);n=len(core);l=a+n+b;v=l;x=l+1;out=[0]*(l+2)
 for i in range(l):
  for j in range(i+1,l):out[i]|=1<<j
 for i in range(n):out[a+i]=(out[a+i]&~(((1<<n)-1)<<a))|(core[i]<<a)
 S=0
 for i in range(a):S|=((p>>(a-i-1))&1)<<i
 S|=cut<<a
 for i in range(b):S|=((r>>(b-i-1))&1)<<(a+n+i)
 for i in range(l):
  out[i]|=1<<v
  if S>>i&1:out[x]|=1<<i
  else:out[i]|=1<<x
 out[v]|=1<<x
 return out
raw=(D/'RECIPES.bin').read_bytes();need(len(raw)==4000008,'recipe size')
miss=[];checked=0;small=[];max_order=0
for target in range(23,1000001,2):
 code=struct.unpack_from('<Q',raw,8*(target//2))[0]
 if not code:miss.append(target);continue
 c,a,b,p,r=decode(code);need(c<80 and a+b<=16 and p<(1<<a) and r<(1<<b),'recipe range')
 core,cut=corecut(c);need((a and p>>(a-1)&1) or (not a and cut),'strongness condition')
 need(formula(c,a,b,p,r)==target,'recipe value');checked+=1;max_order=max(max_order,len(core)+a+b+2)
 if a+b<=5:small.append((target,c,a,b,p,r))
need(miss==[47,63],'finite missing list')
rng=random.Random(216);literal=0
for target,c,a,b,p,r in rng.sample(small,min(120,len(small))):
 need(hp(graph(c,a,b,p,r))==target,'literal control');literal+=1
# All legal templates with padding<=3: the elementary inverse bounds prove this suffices for47/63.
low=set();lowchecks=0
for c in range(80):
 core,cut=corecut(c)
 for total in range(4):
  for a in range(total+1):
   b=total-a
   for p in range(1<<a):
    if (a and not(p>>(a-1)&1)) or (not a and not cut):continue
    for r in range(1<<b):
     value=hp(graph(c,a,b,p,r));need(value==formula(c,a,b,p,r),'small universe formula');low.add(value);lowchecks+=1
need(47 not in low and 63 not in low,'permanent exceptions')
exc=json.loads((D/'EXCEPTIONS.json').read_text());need(set(exc)=={'47','63'},'exception witnesses')
for target,w in exc.items():
 a=w['rows'];cy=w['cycle'];need(hp(a)==int(target),'exception count');need(sorted(cy)==list(range(len(a))) and all(a[x]>>y&1 for x,y in zip(cy,cy[1:]+cy[:1])),'exception cycle')
print(json.dumps({'status':'PASS','max_retained_order':max_order,'extra_exception_witnesses':2,'core_cut_states':len(C),'recipe_checks':checked,'padding_bound':16,'target_interval':[23,1000000],'missing':[47,63],'random_literal_controls':literal,'all_padding_le3_template_checks':lowchecks,'padding_le3_values':sorted(low)},indent=2))
