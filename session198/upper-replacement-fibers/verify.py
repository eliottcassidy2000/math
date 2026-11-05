"""Exact checks for physical upper-child response fibers."""
from pathlib import Path
import json,hashlib
import positive_model as m
from positive_model import source_bank,old,new,GEN,mv,mm,BASE,build,profile,require
from core import cycles
import bilinear,c3
HERE=Path(__file__).resolve().parent

def word(y,w):
 for c in w:y=mv(GEN[c],y)
 return y

def C(s):
 D,E,U,d=s
 return 9*D-3*E-9*U+4*d

def n(s):
 D,E,U,d=s
 return 3*D+12*U-2*d

def hp(T):
 size=len(T);dp=[[0]*size for _ in range(1<<size)]
 for i in range(size):dp[1<<i][i]=1
 for mask,row in enumerate(dp):
  for i,v in enumerate(row):
   if v:
    rem=T[i]&~mask
    while rem:
     b=rem&-rem;rem-=b;dp[mask|b][b.bit_length()-1]+=v
 return sum(dp[-1])

basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
for y in basis:
 D,E,U,d=old(y)
 require(word(y,'JZ')==word(y,'OJ'),'balanced profile reassociation')
 require(old(word(y,'OOJJ'))==(D+8*U,4*U,2*D+8*U,D+4*U),'A row')
 require(old(word(y,'OZJ'))==(D+6*U-2*d,D+2*U-d,D+6*U-d,D+2*U),'B row')
 require(C(old(word(y,'J')))==D-6*E,'precursor response')
 require(C(old(word(y,'OOJJ')))==4*U-5*D,'A response')
 require(C(old(word(y,'OZJ')))==D+2*U-6*d,'B response')
 require(n(old(word(y,'OOJJ')))==25*D+112*U,'A parent row')
 require(n(old(word(y,'OZJ')))==13*D+86*U-18*d,'B parent row')
 for k in range(1,7):
  for l in range(1,7):
   actual=old(word(y,'O'*l+'Z'*k))
   expected=(D+(2**k-1)*2**l*U-(2**k-2)*d,d,2**(l+k)*U-(2**k-1)*d,d)
   require(actual==expected,'residual block')
   require(old(word(y,'O'*l+'J'))==(D+2**l*U,D,D+2**l*U,D),'balanced block')
for h in range(1,101):
 y=new((4*h-1,h,4*h-1,h));s=old(word(y,'JJ'));q=old(word(y,'JJOZJ'))
 require(s[0]+2*s[2]-6*s[3]==1,'lower to upper seed bridge')
 require(C(q)==1 and n(q)==1320*h-347,'lower to upper readout bridge')

pump_tests=[]
for functional,block,alpha,beta,base in [('S','OZ',2,3,4),('A','OOZZ',4,5,16)]:
 for y in basis:
  D,E,U,d=old(y);DD,EE,UU,dd=old(word(y,block))
  require(alpha*UU-beta*DD==alpha*U-beta*D-2*d and dd==d,'arithmetic pump row')
 for target in range(-101,102,2):
  l=0;u=2
  while alpha*u-beta<target+2:l+=1;u*=2
  k=(alpha*u-beta-target)//2
  yy=word(BASE,'O'*l+block*k);D,E,U,d=old(yy)
  require(k>=1 and E==d==1 and alpha*U-beta*D==target,'exact physical response target')
  den=3 if functional=='S' else 5
  expected=((den*u-1)*base**k+1)//den
  require(U==expected and (den*u-1)*base**k+1==den*U,'closed pump formula')
  pump_tests.append((functional,target,l,k,len('O'*l+block*k)+2))

bank=source_bank(5000);records={old(y):''.join(w)for y,w in bank.items()}
A={};B={};Bbridge={};upper=[];strong=0
for y,w in bank.items():
 D,E,U,d=old(y);a,f,e,p=y
 if E:
  strong+=1
  require(f<=a and p<=e and E<=d<=2*E and 2*D>=3*d and U<=2*(D-E),'strong cones')
  if p==0 and f>0:require(2*f>=e,'nonbalanced Z cone')
  if 4*U-5*D==1:
   q=old(word(y,'OOJJ'));require(q==(11*D+2,5*D+1,12*D+2,6*D+1),'A reduced profile')
   require(C(q)==1 and n(q)==165*D+28,'A child')
   A[q]=(old(y),''.join(w))
  if D+2*U-6*d==1:
   q=old(word(y,'OZJ'));require(q==(4*U+4*d+1,5*d+1,4*U+5*d+1,6*d+1),'B reduced profile')
   require(C(q)==1 and n(q)==60*(U+d)+13,'B child')
   B[q]=(old(y),''.join(w))
  if D==U and D==4*d-1:
   ss=old(word(y,'JJ'));Bbridge[ss]=(old(y),''.join(w))
  if D==4*E-1:require(D==U and E==d,'balanced negative-unit necessity')
for y,w in bank.items():
 s=old(y);D,E,U,d=s
 if E and D+2*U-6*d==1:require(s in Bbridge,'B bridge necessity')
 if C(s)==1:
  if D==U and E==d==1:
   j=(D-1).bit_length()-1
   require(D==2**j+1 and j>=1 and old(word(BASE,'O'*(j-1)+'J'))==s,'balanced ray')
   label='ray'
  elif s in A:label='A'
  elif s in B:label='B'
  else:raise RuntimeError('unclassified physical upper source')
  upper.append({'profile':s,'word':''.join(w),'branch':label,'parent':n(s),'child':2*n(s)+1})
for ss,(y,w)in Bbridge.items():
 if ss[2]<=5000:
  require(ss in records and ss[0]+2*ss[2]-6*ss[3]==1,'B bridge sufficiency')

literal=[]
for label,w in [('ray','J'),('A','OJOZOOJJ'),('B','ZJJOZJ')]:
 q=old(word(BASE,w));T,c=build(w);require(profile(T,c)==q,'literal source profile')
 source_parent,pc=build(w+'O');t,tc=build('')
 parent=bilinear.construct('N',t,tc,source_parent,pc)[0]
 source_child,jc=build(w+'J');child=c3.construct('J0',t,tc,source_child,jc)[0]
 actual=(hp(parent),cycles(parent),hp(child),cycles(child))
 require(actual==(n(q),2,2*n(q)+1,3),'literal parent and child')
 literal.append({'branch':label,'source_word':w,'source_profile':q,'parent_order':len(parent),'child_order':len(child),'counts':actual})

out={'status':'PASS','source_cap':5000,'source_profiles':len(bank),'strong_profiles':strong,'upper_response_one':sorted(upper,key=lambda r:r['profile'][2]),'A_seed_profiles':sum(1 for y in bank if old(y)[1]and 4*old(y)[2]-5*old(y)[0]==1),'B_seed_profiles':sum(1 for y in bank if old(y)[1]and old(y)[0]+2*old(y)[2]-6*old(y)[3]==1),'literal_graph_checks':literal,'symbolic_basis_rows':4,'residual_block_basis_tests':4*6*6,'bridge_parameter_tests':100,'exact_arithmetic_target_tests':len(pump_tests),'maximum_pump_source_order':max(r[-1]for r in pump_tests),'dependency_hashes':{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((HERE/'dependencies').glob('*'))if p.is_file()},'scope':'Complete physical upper-response-one fiber and exact construction branches, not cofinite scalar coverage.'}
text=json.dumps(out,sort_keys=True,indent=2)+'\n';(HERE/'VALIDATION.json').write_text(text);print(text,end='')
