import json,pathlib
P=pathlib.Path(__file__).resolve().parent

def parts(n):
 if not n:yield()
 for a in range(2,n+1):
  for s in parts(n-a):yield(a,)+s
def f(p):return (1<<(p-1))+1
def product(c):
 x=1
 for p in c:x*=f(p)
 return x
out=[];checked=0
for n in range(9,25):
 r=n%3;start,coef={0:(9,26),1:(10,46),2:(11,78)}[r];k=(n-start)//3;minimum=coef*5**k
 actual=[];expected=[]
 for c in parts(n-1):
  for j,p in enumerate(c):
   checked+=1
   if p==2:delta=product(c[:j])*product(c[j+1:])+product(c[:max(j-1,0)])*product(c[j+2:])
   else:delta=(1<<(p-2))*product(c[:j]+c[j+1:])
   eligible=False
   if p==2 and 0<j<len(c)-1:
    if r==0:eligible=(c.count(2)==1 and c.count(3)==k+2 and len(c)==k+3)
    if r==1:eligible=(c.count(2)==1 and c.count(4)==1 and c.count(3)==k+1 and sorted((c[j-1],c[j+1]))==[3,4]and len(c)==k+3)
    if r==2:eligible=(c.count(2)==2 and c.count(3)==k+2 and c[j-1]==c[j+1]==3 and len(c)==k+4)
   if delta<minimum:raise RuntimeError('minimum failed')
   if delta==minimum:actual.append((c,j))
   if eligible:expected.append((c,j))
 if actual!=expected:raise RuntimeError('classification mismatch')
 want={0:k+1,1:2*(k+1),2:(k+1)*(k+2)}[r]
 if len(actual)!=want:raise RuntimeError('class count')
 sc=sum(c==c[::-1]and j==len(c)-1-j for c,j in actual)
 if sc!=int(r==0 and k%2==0):raise RuntimeError('rooted converse symmetry')
 out.append({'order':n,'minimum':minimum,'rooted_classes':len(actual),'self_converse_rooted_classes':sc})
receipt={'status':'PASS','marked_block_positions_checked':checked,'maximum_order':24,'results':out}
(P/'chain_equality_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
