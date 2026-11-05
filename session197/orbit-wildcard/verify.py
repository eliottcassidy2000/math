"""Exact independent controls; standard library only, no removable assertions."""
import itertools as it, json, math, pathlib
P=pathlib.Path(__file__).resolve().parent

def gate(ok,msg):
 if not ok: raise RuntimeError(msg)
def validate(r):
 n=len(r)
 gate(all(not(r[i]>>i&1) and r[i]>>n==0 for i in range(n)), 'bad rows')
 gate(all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)), 'not tournament')
def graphs(n):
 pairs=list(it.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  r=[0]*n
  for k,(i,j) in enumerate(pairs):
   if mask>>k&1:r[i]|=1<<j
   else:r[j]|=1<<i
  yield r
def trans(n):return [(1<<n)-(1<<(i+1)) for i in range(n)]
def cyclic(q,ds):return [sum(1<<((i+d)%q) for d in ds) for i in range(q)]
def h_all(r):
 n=len(r);e=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
 for v in range(n):e[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1 and s!=(1<<v):e[s][v]=sum(e[s^(1<<v)][u] for u in range(n) if r[u]>>v&1)
  h[s]=sum(e[s])
 return h
def h_backward(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    d[s][v]=sum(d[t][u] for u in range(n) if r[v]>>u&1) if t else 1
 return sum(d[-1]) if n else 1
def covers(r,h=None):
 n=len(r);h=h or h_all(r);f=[[0]*(n+1) for _ in h];f[0][0]=1
 for s in range(1,1<<n):
  pivot=s&-s;t=s
  while t:
   if t&pivot:
    for k in range(1,n+1):f[s][k]+=h[t]*f[s^t][k-1]
   t=(t-1)&s
 return f[-1]
def cycles(r):
 out=[]
 for k in range(3,len(r)+1,2):
  for s in it.combinations(range(len(r)),k):
   for tail in it.permutations(s[1:]):
    w=(s[0],)+tail
    if all(r[w[i]]>>w[(i+1)%k]&1 for i in range(k)):out.append((sum(1<<i for i in w),k,w))
 return out
def compression(r):
 cs=cycles(r);out={}
 def rec(start,used,count,j):
  out[j]=out.get(j,0)+(1<<count)
  for i in range(start,len(cs)):
   s,k,w=cs[i]
   if not(s&used):rec(i+1,used|s,count+1,j-k+1)
 rec(0,0,0,len(r));return out
ST=[[0]*32 for _ in range(32)];ST[0][0]=1
for n in range(1,32):
 for k in range(1,n+1):ST[n][k]=ST[n-1][k-1]+k*ST[n-1][k]
def substitute(q,root,u):
 verts=[('u',i) for i in range(len(u))]+[('q',i) for i in range(len(q)) if i!=root];r=[0]*len(verts)
 for x,(s,i) in enumerate(verts):
  for y,(t,j) in enumerate(verts):
   if x==y:continue
   arc=(u[i]>>j&1) if s==t=='u' else ((q[i]>>j&1) if s==t=='q' else ((q[root]>>j&1) if s=='u' else (q[i]>>root&1)))
   if arc:r[x]|=1<<y
 return r
def runs(q,root):
 n=len(q);a=[0]*(n+1);ext=[i for i in range(n) if i!=root]
 for p in it.permutations(ext):
  for mask in range(1,1<<n):
   w=[]
   for i in range(n):
    if mask>>i&1:w.append(root)
    if i<n-1:w.append(p[i])
   if all(q[x]>>y&1 for x,y in zip(w,w[1:])):a[mask.bit_count()]+=1
 while not a[-1]:a.pop()
 return a
def clone_count(a,j):return sum(a[k]*math.factorial(k)*ST[j][k] for k in range(1,len(a)))
def beta(a):return [0]+[sum(a[k]*(-1)**(k-b)*math.comb(k,b) for k in range(b,len(a))) for b in range(1,len(a))]
def poly_mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def recurrence(a,e):
 D=len(a)-1;coeff=[1]
 for b in range(1,D+1):coeff=poly_mul(coeff,[-b*b,1])
 vals=[clone_count(a,e+2*t) for t in range(12)]
 for t in range(len(vals)-D):gate(sum(coeff[i]*vals[t+i] for i in range(D+1))==0,'recurrence')
 g=math.gcd(*vals[:D]);gate(g==math.gcd(*vals),'gcd continuation control')
 return g

def factor(n):
 d=2;out=[]
 while d*d<=n:
  if n%d==0:
   e=0
   while n%d==0:n//=d;e+=1
   out.append((d,e))
  d+=1
 if n>1:out.append((n,1))
 return out
def prime_cert(p):
 if p==2:return {'p':2}
 fs=factor(p-1)
 for a in range(2,p):
  if pow(a,p-1,p)==1 and all(math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q,e in fs):break
 else:raise RuntimeError('no Lucas witness')
 return {'p':p,'a':a,'factors':[{'e':e,'certificate':prime_cert(q)} for q,e in fs]}
def check_cert(c):
 p=c['p']
 if p==2:return
 gate(p>2 and p%2==1,'prime certificate order')
 prod=1;seen=set()
 for t in c['factors']:
  q=t['certificate']['p'];gate(q not in seen and q<p and t['e']>0,'factor structure');seen.add(q)
  check_cert(t['certificate']);prod*=q**t['e']
 gate(prod==p-1,'complete factorization')
 a=c['a'];gate(pow(a,p-1,p)==1,'Fermat certificate')
 for q in seen:gate(math.gcd(pow(a,(p-1)//q,p)-1,p)==1,'order certificate')
def autos(r):
 n=len(r);out=[]
 for p in it.permutations(range(n)):
  if all(((r[i]>>j)&1)==((r[p[i]]>>p[j])&1) for i in range(n) for j in range(i+1,n)):out.append(p)
 return out

def main():
 stats={'labelled_profiles':0,'prime_coefficient_checks':0,'context_graph_checks':0,'clone_graph_checks':0}
 banks=[]
 for n in range(1,6):
  for r in graphs(n):
   validate(r);h=h_all(r);f=covers(r,h);c=compression(r)
   gate(h[-1]==h_backward(r),'two HP recurrences')
   gate(sum(c.values())==h[-1],'OCF')
   gate(all(j>0 and j%2==n%2 for j in c),'positive parity')
   gate(f==[sum(w*ST[j][k] for j,w in c.items()) for k in range(n+1)],'Touchard resolution')
   if n%2:
    for ell in (3,5,7):
     if ell-1<=n:gate(f[ell-1]%ell==0,'prime coefficient');stats['prime_coefficient_checks']+=1
   stats['labelled_profiles']+=1
   if n<=4:banks.append((r,f,c))
 templates=[('C3',cyclic(3,[1]),[0,3,1],3,1),('C5',cyclic(5,[1,2]),[0,15,13,2],15,1),('C7',cyclic(7,[1,2,3]),[0,175,227,98,6],5,1)]
 table={}
 for name,q,expected,go,ge in templates:
  a=runs(q,0);gate(a==expected,'complete run inventory '+name);b=beta(a)
  gate((recurrence(a,1),recurrence(a,2))==(go,ge),'context parity gcd')
  D=len(a)-1
  table[name]={'a':a[1:],'beta':b[1:],'odd_initial':[clone_count(a,1+2*t) for t in range(D)],'even_initial':[clone_count(a,2+2*t) for t in range(D)],'gcd':[go,ge]}
  for m in range(1,2*D+1):
   r=substitute(q,0,trans(m));val=clone_count(a,m)
   gate(val==sum(b[x]*x**m for x in range(1,len(b))),'power formula')
   gate(h_all(r)[-1]==val and h_backward(r)==val,'actual clone graph');stats['clone_graph_checks']+=1
  for u,f,c in banks:
   value=sum(a[k]*math.factorial(k)*f[k] for k in range(1,min(len(a),len(f))))
   mixture=sum(w*clone_count(a,j) for j,w in c.items())
   actual=h_all(substitute(q,0,u))[-1]
   gate(value==mixture==actual,'arbitrary physical context')
   gate(actual%(go if len(u)%2 else ge)==0,'parity divisor')
   stats['context_graph_checks']+=1
 A=cyclic(7,[1,2,3]);B=cyclic(7,[1,2,4]);profiles=[]
 for q,expected,linear in [(A,[0,175,819,1113,630,168,21,1],34),(B,[0,189,819,1113,630,168,21,1],48)]:
  f=covers(q);c=compression(q)
  gate(f==expected,'seven full profile');gate(c=={7:1,5:28,3:112,1:linear},'seven cycle compression')
  profiles.append(f[1:])
 for q in graphs(3):
  for root in range(3):
   gate(h_all(substitute(q,root,B))[-1]-h_all(substitute(q,root,A))[-1]==14*h_all(q)[-1],'physical switch')
 r=[0]*14
 for offset in (0,7):
  for i in range(7):r[offset+i]|=B[i]<<offset
 for i in range(7):
  for j in range(7):
   if (j-i)%7 in (1,3):r[i]|=1<<(7+j)
   else:r[7+j]|=1<<i
 expected=[1302,2604,5208,10289,4322,8517,651,2863,5726,11325,6394,12661,8939,1495]
 gate(r==expected,'physical rows');validate(r)
 forward=h_all(r)[-1];backward=h_backward(r);gate(forward==backward==15144269,'two-orbit HP count')
 p=2163467;gate(forward==7*p,'7p factorization');cert=prime_cert(p);check_cert(cert)
 gate(all(p%d for d in range(2,math.isqrt(p)+1)),'independent trial primality')
 perm=[(i+1)%7 for i in range(7)]+[7+(i+1)%7 for i in range(7)]
 gate(all(((r[i]>>j)&1)==((r[perm[i]]>>perm[j])&1) for i in range(14) for j in range(14)),'rotation')
 degrees=[x.bit_count() for x in r];gate(degrees==[5]*7+[8]*7,'degree split')
 ab=autos(B);gate(len(ab)==21,'Paley automorphism bank');full=[]
 for x in ab:
  for y in ab:
   perm2=list(x)+[7+i for i in y]
   if all(((r[i]>>j)&1)==((r[perm2[i]]>>perm2[j])&1) for i in range(7) for j in range(7,14)):full.append(perm2)
 gate(len(full)==7,'independent complete automorphism count')
 modules=[]
 for s in range(1<<14):
  if not (2<=s.bit_count()<=13):continue
  if all((s>>v&1) or (r[v]&s in (0,s)) for v in range(14)):modules.append(s)
 gate(not modules,'all subsets module check')
 for v in range(14):
  seen=1<<v
  while True:
   new=seen
   for u in range(14):
    if seen>>u&1:new|=r[u]
   if new==seen:break
   seen=new
  gate(seen==(1<<14)-1,'strongness')
 witness={'rows':r,'H':forward,'prime_quotient':p,'rotation':perm,'all_automorphisms':full,'primality_certificate':cert}
 (P/'witness.json').write_text(json.dumps(witness,indent=2,sort_keys=True)+'\n')
 receipt={'controls':stats,'template_profiles':table,'seven_cover_profiles':profiles,'two_orbit':{'H':forward,'prime':p,'aut_order':len(full),'module_subset_universe':1<<14,'nontrivial_modules':len(modules),'trial_division_bound':math.isqrt(p)},'status':'PASS'}
 text=json.dumps(receipt,indent=2,sort_keys=True)+'\n';(P/'receipt.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
