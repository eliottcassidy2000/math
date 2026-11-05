"""Standard-library, optimization-safe verification of the new exact reduction."""
from pathlib import Path
from itertools import combinations
import json,hashlib
import c3,bilinear
from core import require,inventories,cycles,glue
from source_oracle import source_bank,build_source
from coverage import *
P=Path(__file__).resolve().parent

# Separate direct digit gadget counter; it does not call constants().
def gadget_constants(L):
 bits=[] if L==1 else [1-int(x) for x in bin(L-1)[2:][::-1]]
 # Transitive base and root; exact subset counters independently recover constants.
 T=[0]*(len(bits)+1)
 for i,b in enumerate(bits,1):
  T[0 if b else i]|=1<<(i if b else 0)
  for j in range(i+1,len(T)):T[i]|=1<<j
 p=inventories(T)
 return p['H'],p['R'][0],1<<len(bits)

def source_value(recipe):
 # The binary constants use a separate iterative low-bit recurrence.
 D=E=U=d=1
 for kind,L in recipe:
  if kind=='O':E=0;U*=2;continue
  trail=[];k=L
  while k>1:trail.append(k);k=(k+1)//2
  g=Q=B=1
  for k in reversed(trail):
   if k&1:ell=(k+1)//2;g+=ell-ell%2;Q+=B
   B*=2
  D,E,U,d=g*D-(L&1)*E+L*U,D-(E if L==1 else 0),Q*D+B*U,D
 return D,E,U,d

def hc(T):
 # Independent forward subset paths, separate from inherited end-subset counter.
 n=len(T);N=1<<n;A=[[0]*n for _ in range(N)];C=[[0]*n for _ in range(N)]
 for i in range(n):A[1<<i][i]=1
 C[1][0]=1
 for S in range(1,N):
  for u in range(n):
   if not S>>u&1:continue
   outs=T[u]&~S
   while outs:
    bit=outs&-outs;outs-=bit;v=bit.bit_length()-1
    A[S|bit][v]+=A[S][u];C[S|bit][v]+=C[S][u]
 return sum(A[-1]),sum(C[-1][u]for u in range(n)if T[u]&1)

def main():
 for L in range(1,100):
  from core import constants
  require(constants(L)==gadget_constants(L),'literal gadget constants')
 checks=0
 for kind,M in EXPECTED.items():
  num,den=LOWER[kind]
  for i,v in enumerate(BASIS):
   for j,w in enumerate(BASIS):
    require(c3.value(kind,v,w)==M[i][j],'new matrix entry')
    require(den*M[i][j]>=num*WEIGHT[i]*WEIGHT[j],'product lower matrix')
    checks+=1
 for v in BASIS:
  for w in BASIS:
   require(c3.value('T',v,w)==c3.value('B',w,v),'T transpose')
   require(c3.value('J2',v,w)==c3.value('J1',w,v),'J2 transpose')
   require(c3.value('J3',v,w)==c3.value('J0',w,v),'J3 transpose')
   for k in ('A','C'):require(c3.value(k,v,w)==c3.value(k,w,v),'symmetric form')
 for v in BASIS:
  require(c3.value('T',ONE,v)==c3.value('J0',v,ONE)==c3.value('J3',ONE,v),'singleton transpose boundary')
 require([c3.value('C',ONE,v) for v in BASIS]==[21,16,36,33],'C singleton row')
 require([c3.value('T',ONE,v) for v in BASIS]==[12,9,15,15],'T singleton row')
 for k in range(1,12):
  z=(1<<k)+1;pz=(2*z-1,z-1,2*z,z)
  require(source_value((('M',1<<k),('M',1)))==pz,'sharp physical source family')
  require(c3.value('T',ONE,pz)==15*z,'sharp uniform cutoff family')
  for n in range(1,12):
   t=(1,0,1<<n,1)
   require(c3.value('J0',pz,t)==12*z*(1<<n)+3*z,'J0 sharp family')
   require(c3.value('B',(1,0,2,1),t)==14*(1<<n)+3,'B sharp family')
   for m in range(1,7):
    u=(1,0,1<<m,1);a=1<<m;b=1<<n
    require(c3.value('A',u,t)==7*a*b+a+b-3,'A sharp family')
    require(c3.value('C',u,t)==14*a*b+3*a+3*b-1,'C sharp family')
 t1=(1,0,2,1)
 require(c3.value('J1',t1,t1)==c3.value('J2',t1,t1)==29,'J1 exact equality')
 B=5001
 C1,bank1=c1_support(B);C2,bank2,pairs2=bilinear.support(B);C3,bank3,tests=c3_support(B)
 old,oldbank,oldtests=c3.support(B)
 require(set(C3)==set(old),'five-form and eight-form full images agree')
 for v,recipe in bank1.items():
  require(source_value(recipe)==v,'independent source recurrence certificate')
  D,E,U,d=v;require(E<=d,'physical E<=d')
  if v!=ONE:require(min(ycoords(v))>=0 and U>=2*d,'physical nonnegative cone')
 for kind,support in ((1,C1),(2,C2),(3,C3)):
  for H,r in support.items():
   if kind==1:
    v=source_value(r['recipe']);require(v==r['source'] and v[0]-v[1]+v[2]==H,'C1 certificate')
   elif r['kind']=='R5':require(H==15,'R5 exception')
   else:
    v=source_value(r['left_recipe']);w=source_value(r['right_recipe'])
    require(v==tuple(r['left']) and w==tuple(r['right']),'pair certificate sources')
    f=bilinear.value if kind==2 else c3.value
    require(f(r['kind'],v,w)==H,'pair certificate value')
    if kind==3:require(max(v[2],w[2])*15<=2*H,'improved cap on every witness')
 selected=[]
 for H in (15,29,31,67,533,821,1067,1339,1861,2647):
  for kind,support,realizer in ((2,C2,bilinear.realize),(3,C3,c3.realize)):
   if H not in support:continue
   T,r=realizer(support[H]);require(hc(T)==(H,kind),'independent physical witness path/cycle count')
   selected.append(dict(H=H,C=kind,rows=T,root=r,certificate=support[H]))
 missing1=[x for x in range(3,B+1,2)if x not in C1]
 missing12=[x for x in range(3,B+1,2)if x not in C1 and x not in C2]
 missing123=[x for x in missing12 if x not in C3]
 require(missing12==[7,21,67] and missing123==[7,21],'finite support controls')
 report=dict(status='PASS',bound=B,matrix_entries=checks,source_recurrence_checks=len(bank1),C1_values=len(C1),C2_values=len(C2),C3_values=len(C3),new_C3_source_bank=len(bank3),old_C3_source_bank=len(oldbank),new_C3_tests_by_form=tests,new_C3_pair_tests=sum(tests.values()),old_C3_pair_tests=oldtests,C2_pair_tests=pairs2,C1_missing=missing1,C12_missing=missing12,C123_missing=missing123,selected_physical_witnesses=len(selected),max_selected_graph_order=max(len(r['rows']) for r in selected))
 encoded=json.dumps(report,indent=2,sort_keys=True)+'\n';(P/'VALIDATION.json').write_text(encoded)
 (P/'WITNESSES.json').write_text(json.dumps(selected,indent=2,sort_keys=True)+'\n')
 # Full finite support preserves compact source recipes for every value and class.
 (P/'SUPPORT.json').write_text(json.dumps(dict(bound=B,C1=C1,C2=C2,C3=C3),separators=(',',':'),sort_keys=True)+'\n')
 print(encoded,end='')
if __name__=='__main__':main()
