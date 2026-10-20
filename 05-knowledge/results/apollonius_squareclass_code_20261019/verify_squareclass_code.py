#!/usr/bin/env python3
import itertools,json,sympy as s
checks=0
def req(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def rank2(rows):
 rows=list(rows);rank=0
 for col in range(max((x.bit_length() for x in rows),default=0)-1,-1,-1):
  i=next((i for i in range(rank,len(rows)) if rows[i]>>col&1),None)
  if i is None:continue
  rows[rank],rows[i]=rows[i],rows[rank]
  for j in range(len(rows)):
   if j!=rank and rows[j]>>col&1:rows[j]^=rows[rank]
  rank+=1
 return rank
# Algebraic elimination, exact independent coefficients.
a,b,c,r,u,v=s.symbols('a b c r u v')
x0=(a*a+r*r-u*u)/(2*a);x1=(r-u)/a
y0=(b*b+c*c+r*r-v*v-2*b*x0)/(2*c);y1=(r-v-b*x1)/c
A=x1*x1+y1*y1-1;B=2*(x0*x1+y0*y1-r);C=x0*x0+y0*y0-r*r
P=(a*a-(r-u)**2)*(b*b+c*c-(r-v)**2)*((a-b)**2+c*c-(u-v)**2)
req(s.cancel(B*B-4*A*C-P/(a*a*c*c))==0,'discriminant factorization')
rows=[0b010101,0b101001,0b100110,0b011010]
req(rank2(rows)==3,'rank three')
rel=[mask for mask in range(16) if __import__('functools').reduce(int.__xor__,[rows[i] for i in range(4) if mask>>i&1],0)==0]
req(rel==[0,15],'unique relation')
columns=[sum(((rows[i]>>j)&1)<<i for i in range(4)) for j in range(6)]
req(set(columns)=={3,5,9,6,10,12},'K4 columns')
image={__import__('functools').reduce(int.__xor__,[columns[j] for j in range(6) if m>>j&1],0) for m in range(64)}
req(image=={x for x in range(16) if x.bit_count()%2==0},'even code')
hist={z:sum(2*(4-x.bit_count())==z for x in image) for z in [0,4,8]}
req(hist=={0:1,4:6,8:1},'fixed counts')
graphs=0
for n in range(1,6):
 E=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(E)):
  es=[e for j,e in enumerate(E) if mask>>j&1];parent=list(range(n))
  def find(x):
   while parent[x]!=x:x=parent[x]
   return x
  for i,j in es:parent[find(i)]=find(j)
  comps=len({find(i) for i in range(n)})
  cols=[(1<<i)|(1<<j) for i,j in es]
  req(rank2(cols)==n-comps,'graph rank');graphs+=1
# Exhaustive normalized finite-field generic configurations.
fields={}
for q in [5,7,11]:
 hist={0:0,4:0,8:0};excluded=0
 for b,c,r,u,v in itertools.product(range(q),range(1,q),range(1,q),range(1,q),range(1,q)):
  counts=[];ds=[];generic=True
  for su,sv in [(1,1),(1,-1),(-1,1),(-1,-1)]:
   U=su*u%q;V=sv*v%q;inv2=pow(2,-1,q);invc=pow(c,-1,q)
   X0=(1+r*r-U*U)*inv2%q;X1=(r-U)%q
   Y0=(b*b+c*c+r*r-V*V-2*b*X0)*inv2*invc%q;Y1=(r-V-b*X1)*invc%q
   aa=(X1*X1+Y1*Y1-1)%q;bb=2*(X0*X1+Y0*Y1-r)%q;cc=(X0*X0+Y0*Y0-r*r)%q
   d=(bb*bb-4*aa*cc)%q;prod=(1-(r-U)**2)*(b*b+c*c-(r-V)**2)*((1-b)**2+c*c-(U-V)**2)%q
   req(d*c*c%q==prod,'finite discriminant')
   if aa==0 or cc==0 or d==0:generic=False;break
   ds.append(d);counts.append(2 if pow(d,(q-1)//2,q)==1 else 0)
  if not generic:excluded+=1;continue
  total=sum(counts);req(total in hist,'generic forbidden count')
  req(pow(__import__('functools').reduce(lambda x,y:x*y%q,ds,1),(q-1)//2,q)==1,'parity relation')
  hist[total]+=1
 fields[q]={'generic_counts':hist,'excluded_degenerate':excluded}
print(json.dumps({'status':'PASS','checks':checks,'graph_count':graphs,'binary_relations':rel,'fixed_point_histogram':{0:1,4:6,8:1},'finite_fields':fields},indent=2))
