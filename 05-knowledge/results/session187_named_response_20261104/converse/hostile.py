import json
from verify import hp,ext,setup,conv,need

def polynomial(a,b):
 n=len(a);sm=[x for x in range(1<<n) if x.bit_count()<=2];vals={(x,y):hp(ext(a,x,y))-hp(ext(b,x,y)) for x in sm for y in sm};coeff=[]
 for x,y in vals:
  v=0;i=x
  while True:
   j=y
   while True:
    v+=(-1)**((x^i).bit_count()+(y^j).bit_count())*vals[i,j]
    if not j:break
    j=(j-1)&y
   if not i:break
   i=(i-1)&x
  if v:coeff.append([x,y,v])
 return coeff
q=[0]*10
base=[8,1,3,6];q[:4]=base
for cycle in ([4,5,6],[7,8,9]):
 for u,v in zip(cycle,cycle[1:]+cycle[:1]):q[u]|=1<<v
for u in [4,5,6]:
 for v in [0,1,2,3,7,8,9]:q[u]|=1<<v
for u in range(4):
 for v in [7,8,9]:q[u]|=1<<v
A,B,f,lift,parts=setup(q,[2,0]);Ar,Br,fr,liftr,_=setup(q,[0,1]);need(f==fr==9,'factor9')
x=(1<<2)|(1<<4);y=(1<<2)|(1<<1);terms=[];total=0;i=x
while True:
 j=y
 while True:
  v=hp(ext(A,lift(i),lift(j)))-hp(ext(Ar,lift(i),lift(j)))
  sign=(-1)**((x^i).bit_count()+(y^j).bit_count());total+=sign*v
  terms.append({'active_a':i,'active_c':j,'ambient_a':lift(i),'ambient_c':lift(j),'difference':v,'sign':sign})
  if not j:break
  j=(j-1)&y
 if not i:break
 i=(i-1)&x
ambient=polynomial(B,Br);glob=polynomial(B,conv(B));need(total==18,'quartic amplification');need(max(x.bit_count()+y.bit_count() for x,y,v in ambient)==4,'internal quartic');need(max(x.bit_count()+y.bit_count() for x,y,v in glob)==3,'whole cubic')
print(json.dumps({'status':'PASS','quotient_rows':q,'ambient_rows':A,'ambient_reversed_module_rows':Ar,'active_rows':B,'active_reversed_rows':Br,'factor':f,'specialized_ambient_quartic':total,'sixteen_literal_queries':terms,'complete_active_module_reversal_polynomial':ambient,'complete_active_global_converse_polynomial':glob,'omitted_basis_coefficients':'zero in full named bidegree-(2,2) basis'},indent=2,sort_keys=True))
