"""Independent endpoint-subset DP. No producer imports; checks survive python -O."""
import json,itertools,math

def demand(ok,why):
 if not ok: raise ValueError(why)
def endpoints(a):
 n=len(a);demand(all(not(a[i]>>i&1) for i in range(n)),'loops')
 demand(all(((a[i]>>j&1)+(a[j]>>i&1)==1) for i in range(n) for j in range(i)),'tournament')
 d=[{} for s in range(1<<n)]
 for i in range(n):d[1<<i][i,i]=1
 for s in range(1,1<<n):
  for (first,last),count in list(d[s].items()):
   for v in range(n):
    if not(s>>v&1) and a[last]>>v&1:
     key=first,v;new=d[s|1<<v];new[key]=new.get(key,0)+count
 return d
def stats(a,m):
 d=endpoints(a);full=(1<<len(a))-1
 totals=[sum(x.values()) for x in d];totals[0]=1
 start=[sum(c for (f,l),c in x.items() if m>>f&1) for x in d]
 end=[sum(c for (f,l),c in x.items() if m>>l&1) for x in d]
 Q=A=G=T=0
 for s in range(full+1):
  r=full^s;Q+=totals[s]*totals[r];A+=totals[s]*start[r];G+=end[s]*totals[r];T+=end[s]*start[r]
 return [totals[full],Q,A,G,start[full],end[full],T]
def join(a,b):
 n=len(a);return [x|(((1<<len(b))-1)<<n) for x in a]+[x<<n for x in b]
def close(a,m):
 b=join(a,[0]);n=len(a);b.append(0)
 for i in range(n+1):
  if i<n and m>>i&1:b[n+1]|=1<<i
  else:b[i]|=1<<(n+1)
 return b
def strong(a,m):
 b=join(join([0],a),[0]); n=len(b);b.append(0)
 for i in range(n):
  if i==0 or 1<=i<=len(a) and m>>(i-1)&1:b[n]|=1<<i
  else:b[i]|=1<<n
 return b
def compose(c,d):
 h,q,A,G,S,E,T=c;k,r,B,J,V,F,U=d
 return [h*k,q*r,A*r+h*B,J*q+k*G,S*k,h*F,A*J+T*k+h*U+E*V]
B=[26,44,57,16,34,9]; M=[9,10,17,18]
expected=[[29,320,81,129,7,15,22],[29,320,81,117,7,13,16],[29,320,81,93,7,9,12],[29,320,81,81,7,7,14]]
demand([stats(B,m) for m in M]==expected,'rich stats')
weights=[stats(B,1<<i)[2] for i in range(6)]
demand(weights==[50,50,86,31,31,43],'weights')
demand([m for m in range(64) if stats(B,m)[2]==81]==M,'entire alpha81 fiber')
closures=[sum(endpoints(strong(B,m))[-1].values()) for m in M]
demand(closures==[673,691,719,729],'strong closure')
# Residue proof uses pairwise distinct G mod Q, coprime H and Q.
demand(len({c[3]%320 for c in expected})==4 and math.gcd(29,320)==1,'all-word induction hypotheses')
SW=[([64,121,35,37,13,81,28],[12,20],[[95,1276,312,270,24,18,38],[95,1276,312,270,24,18,36]]),([20,85,72,195,204,159,33,71],[47,181],[[389,5600,3369,3456,250,271,1858],[389,5600,3369,3456,250,271,1860]])]
result=[]
for a,ms,ex in SW:
 got=[stats(a,m) for m in ms];demand(got==ex,'switch stats');demand(math.gcd(*got[0][:2])==1,'coprime')
 values=[sum(endpoints(close(a,m))[-1].values()) for m in ms]
 demand(values==[2*c[2]+c[0]-c[6] for c in got],'literal star counts')
 # Check literal composition, for all binary singleton contexts up to total length 3.
 count=0
 for length in range(4):
  for bits in itertools.product((0,1),repeat=length):
   for cut in range(length+1):
    outputs=[]
    for m,c in zip(ms,got):
     aa=a[:];mm=m;cc=c[:]
     for bit in reversed(bits[:cut]):aa=join([0],aa);mm=(mm<<1)|bit;cc=compose([1,2,bit,bit,bit,bit,0],cc)
     for bit in bits[cut:]:mm|=bit<<len(aa);aa=join(aa,[0]);cc=compose(cc,[1,2,bit,bit,bit,bit,0])
     actual=stats(aa,mm);demand(actual==cc,'literal composition');outputs.append(actual)
    demand(outputs[0][:-1]==outputs[1][:-1] and outputs[0][-1]-outputs[1][-1]==got[0][-1]-got[1][-1],'context gap');count+=1
 result.append(dict(block=a,masks=ms,coordinates=got,closure=values,literal_contexts=count))
# Automorphism counts by independent literal permutation action, for all three bases.
autos=[]
for a in [B]+[s[0] for s in SW]:
 total=sum(all(bool(a[i]>>j&1)==bool(a[p[i]]>>p[j]&1) for i in range(len(a)) for j in range(i)) for p in itertools.permutations(range(len(a))))
 demand(total==1,'rigid');autos.append(total)
# Additional structural checks.
def transfer(c):
 h,q,A,G,S,E,T=c
 w={(0,1):S,(0,2):h-S,(1,0):h-E,(2,0):E,(2,1):T,(1,1):A-S-T,(2,2):G-E-T,(1,2):q-2*h-A+S-G+E+T}
 states=[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)];mat=[[0]*6 for _ in range(6)]
 for j,(p,r) in enumerate(states):
  for (pp,rr),v in w.items():
   if r==0 and rr==2:continue
   i=states.index((pp or p,int(bool(r or rr))));mat[i][j]+=v
 return mat
for a,ms,cs in SW:
 for start in range(len(a)):
  seen={start}
  while True:
   more=seen|{j for i in seen for j in range(len(a)) if a[i]>>j&1}
   if more==seen:break
   seen=more
  demand(len(seen)==len(a),'strong core')
 mat=[transfer(c) for c in cs];delta=cs[0][-1]-cs[1][-1]
 demand(all(mat[0][i][j]-mat[1][i][j]==-delta*([0,0,0,1,0,-1][i])*([1,0,1,0,1,0][j]) for i in range(6) for j in range(6)),'central defect sign')
for c in expected+[c for a,m,cs in SW for c in cs]+[[1,2,b,b,b,b,0] for b in (0,1)]:
 mat=transfer(c);w=[0,0,0,1,0,-1];rho=[1,0,1,0,1,0]
 demand([sum(mat[i][j]*w[j] for j in range(6)) for i in range(6)]==[c[0]*x for x in w],'M w')
 demand([sum(rho[i]*mat[i][j] for i in range(6)) for j in range(6)]==[c[0]*x for x in rho],'rho M')

# Closure unsaturation is not inherited from its body: explicit hostile controls.
closure_controls=[]
for (a,ms,cs),want in zip(SW, [[(681,15600,3),(683,15664,1)],[(5269,111760,11),(5267,111608,1)]]):
 got=[]
 for m in ms:
  h,q=stats(close(a,m),0)[:2];got.append((h,q,math.gcd(h,q)))
 demand(got==want,'closure saturation controls')
 closure_controls.append(dict(base=a,masks=ms,closure_H_Q_gcd=got))

print(json.dumps(dict(status='PASS',method='subset DP by both endpoints, complementary-subset products',rich=expected,rich_strong_counts=closures,switches=result,automorphism_orders=autos,closure_saturation_controls=closure_controls),indent=2))
