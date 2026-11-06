"""Independent literal graph/sector and conditional-arithmetic controls.
Uses private vertices first, then marked r,c, unlike the author labeling.
No author imports; no scalar-extraction constructor is assumed implemented.
"""
import functools,json

def need(b,m):
 if not b:raise RuntimeError(m)
def g(u):
 answer=1
 while u>1:
  if u%2:u=(u+1)//2;answer+=u-u%2
  else:u//=2
 return answer

def build(u):
 m=(u-1).bit_length();r=m;c=m+1;n=m+2;rows=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if j<m:winner=i
   elif j==r:winner=i
   elif i==r:winner=r
   else:winner=c if (u-1)>>(m-1-i)&1 else i
   loser=j if winner==i else i;rows[winner]|=1<<loser
 return tuple(rows),r,c

def literal(rows,r,c):
 n=len(rows);full=(1<<n)-1
 need(all(rows[i]>>i&1==0 for i in range(n)),'diagonal')
 need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)),'tournament')
 @functools.cache
 def walk(last,left,first):
  if not left:return 1,int(first==1)
  choices=rows[last]&left;H=F=0
  while choices:
   bit=choices&-choices;choices-=bit;v=bit.bit_length()-1
   mark=first or (1 if v==r else 2 if v==c else 0)
   h,f=walk(v,left^bit,mark);H+=h;F+=f
  return H,F
 def paths(S):
  if not S:return 1,0
  H=F=0;todo=S
  while todo:
   bit=todo&-todo;todo-=bit;v=bit.bit_length()-1
   h,f=walk(v,S^bit,1 if v==r else 2 if v==c else 0);H+=h;F+=f
  return H,F
 @functools.cache
 def cycle(last,left):
  if not left:return int(bool(rows[last]>>r&1))
  todo=rows[last]&left;out=0
  while todo:
   bit=todo&-todo;todo-=bit;out+=cycle(bit.bit_length()-1,left^bit)
  return out
 @functools.cache
 def ends_c(last,left):
  if not left:return int(last==c)
  todo=rows[last]&left;out=0
  while todo:
   bit=todo&-todo;todo-=bit;out+=ends_c(bit.bit_length()-1,left^bit)
  return out
 H,F=paths(full);D=paths(full^(1<<r))[0]
 E=sum(ends_c(v,full^(1<<r)^(1<<v)) for v in range(n) if v!=r)
 return H,F,H-F,cycle(r,full^(1<<r)),D,E

def extend(rows,r,c,kind):
 m=r;n=len(rows);z=m;rr=r+1;cc=c+1;out=[0]*(n+1)
 relabel=lambda i:i if i<m else i+1
 for i,row in enumerate(rows):
  for j in range(n):
   if row>>j&1:out[relabel(i)]|=1<<relabel(j)
 out[z]|=1<<rr
 if kind=='A':out[cc]|=1<<z
 else:out[z]|=1<<cc
 for i in range(m):out[i]|=1<<z
 return tuple(out),rr,cc

records=[];extensions=0
for u in range(1,130):
 rows,r,c=build(u);H,F,A,C,D,E=literal(rows,r,c)
 need((H,F,A,C,D,E)==(u+g(u)-u%2,u,g(u)-u%2,int(u>=2),g(u),u%2),'literal source profile')
 need(len(rows)==2+(u-1).bit_length(),'source order')
 records.append((u,H,F,A))
 if u>=2:
  for kind,v in [('A',2*u),('B',2*u-1)]:
   enlarged=extend(rows,r,c,kind);need(enlarged==build(v),'physical edge extension')
   h,f,a,cy,dd,ee=literal(*enlarged)
   if kind=='A':need((h,f,a)==(2*u+g(u),2*u,A+u%2),'move A arithmetic')
   else:need((h,f,a)==(3*u+g(u)-u%2-2,2*u-1,A+u-1),'move B arithmetic')
   need(cy==1,'extension cycle');extensions+=1
# A remains physically canonical at u=1; B does not, so its graph identity requires u>=2.
need(extend(*build(1),'A')==build(2),'A singleton-sector boundary')
need(extend(*build(1),'B')!=build(1),'B noncanonical boundary retained')
first13=[H for u,H,F,A in records[:13]]
need(first13==[1,3,5,5,9,9,11,9,17,15,19,15,23] and 13 not in first13,'complete H13 family hostile')
h13=(18,8,3,5,14);h13profile=literal(h13,0,1);need(h13profile[0]==13,'H13 unrestricted witness')
need(literal(*build(39))[0]==63,'E39 repair')
for k in range(1,101):
 u=2**k+1;need(g(u)==u,'nonproduct infinite-family arithmetic')
 need(u+g(u)-1==2*u-1 and 2*u+g(u)==3*u,'ratio family')
# Algebraic stress control of the explicitly CONDITIONAL reduction.
# No graph with the two hypothetical extraction counts is claimed constructed.
conditional_cases=0
for n in range(1,100002,2):
 if n in (7,21):continue
 a=0;b=n
 while b%3==0:a+=1;b//=3
 if b==1:reconstructed=3**a
 elif b==5:reconstructed=5*3**a
 elif b==7:
  need(a>=2,'unrepaired exception');reconstructed=63*3**(a-2)
 else:
  sign=1 if b%6==1 else -1
  u=(b-sign)//6
  need(u>=2 and 6*u+sign==b,'OPEN extraction parameter')
  reconstructed=(6*u+sign)*3**a
 need(reconstructed==n,'conditional target factorization');conditional_cases+=1
print(json.dumps({'status':'PASS','scope':'Literal finite controls and conditional arithmetic only; extraction and eventual coverage OPEN','indices_checked':129,'literal_extensions_checked':extensions,'largest_graph_order':len(build(258)[0]),'H13_hostile_table':first13,'unrestricted_H13':h13profile[0],'unrestricted_H13_cycles':h13profile[3],'E39_H':63,'nonproduct_arithmetic_k_max':100,'conditional_factorization_cases':conditional_cases,'conditional_factorization_cap':100001,'excluded_from_conditional_claim':[7,21]},indent=2,sort_keys=True))
