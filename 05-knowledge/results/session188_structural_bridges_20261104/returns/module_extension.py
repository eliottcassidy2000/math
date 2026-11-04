"""Exact arbitrary-module extension checks; standard library only."""
import itertools,json,hashlib
from pathlib import Path
from inherited import tables,strong,tournament
from verify import response,signature,need
from module_coeffs import module_coeffs
from carrier import base,outer,v

def replace(A):
 m=len(A);R=[A[i]+sum(1<<(m+j) for j,u in enumerate(outer) if base[v]>>u&1) for i in range(m)]
 R += [sum(1<<(m+j) for j,w in enumerate(outer) if base[u]>>w&1)+(((1<<m)-1) if base[u]>>v&1 else 0) for u in outer]
 return R,{m+j for j,u in enumerate(outer) if 42>>u&1}

def rows_from_mask(n,mask):
 R=[0]*n
 for k,(a,b) in enumerate(itertools.combinations(range(n),2)):
  if mask>>k&1:R[a]|=1<<b
  else:R[b]|=1<<a
 return R

def decorated(A):
 n=len(A);N=(1<<n)-1;h,e=tables(A);rev=[sum(1<<j for j in range(n) if A[j]>>i&1) for i in range(n)];_,s=tables(rev)
 q2=[sum(h[B]*h[S^B] for B in range(1,N+1) if B!=S and B&S==B) for S in range(N+1)]
 Z=[]
 for u in range(n):
  E1=e[N][u];S1=s[N][u]
  E2=sum(e[B][u]*h[N^B] for B in range(1,N));S2=sum(s[B][u]*h[N^B] for B in range(1,N))
  E3=sum(e[B][u]*q2[N^B] for B in range(1,N+1))
  Z.append((E1,E2,E3,S1,S2))
 return h[N],q2[N],Z

def run():
 expected={'h':{'0':5,'1':17,'2':9},'Jempty':{'0':2,'1':8,'2':6},'E':{'1':25,'2':44,'3':18},'S':{'1':1,'2':5}}
 need(module_coeffs()==expected,'nonnegative exterior skeleton table')
 counts={'labelled_modules':0,'marked_modules':0,'singleton_formulas':0,'injection_checks':0};small=[]
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):
   A=rows_from_mask(n,mask);h,q2,Z=decorated(A);B,W=replace(A);R=response(B);counts['labelled_modules']+=1
   need(tournament(B) and strong(B),'literal module substitution strong')
   need(R[0]==17*h+9*q2 and signature(R,W)[4]==8*h+6*q2,'unmarked arbitrary module formulas')
   for u,(E1,E2,E3,S1,S2) in enumerate(Z):
    pred=9*h+3*q2-25*E1-44*E2-18*E3-S1-5*S2
    z=signature(R,W|{u});need(z[0]-z[4]==pred,'singleton arbitrary-module constant')
    need(h<=E1+E2 and q2<=2*(E2+E3),'injection inequalities');counts['singleton_formulas']+=1;counts['injection_checks']+=2
    need(pred<=-16*E1-29*E2-12*E3-S1-5*S2,'paid positive terms')
    if n>=2:need(E2>=1 and S2>=1 and E3>=2**(n-1)-2 and pred<=-10-12*2**(n-1),'universal order bound')
   for ports in range(1,1<<n):
    U={u for u in range(n) if ports>>u&1};z=signature(R,W|U);c=z[0]-z[4];counts['marked_modules']+=1
    for u in U:
     one=signature(R,W|{u});need(c<=one[0]-one[4],'actual marking monotonicity')
    if n>=4:need(c<=-106,'all larger module orders below targets')
    if n>=2:
     for C in (-103,-101):
      if c>=C:need(3*c<C<c and c-min(z[2],z[3])<C,'complete small-order contextual reduction')
    if n==3 and A==[2,4,1] and len(U)==1:need(c==-115,'directed triangle hostile');small.append({'A':A,'U':sorted(U),'signature':z,'c':c})
 result={'status':'PASS','counts':counts,'coefficients':expected,'triangle_controls':small,'proof_sha256':hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest(),'source_sha256':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest() for p in ('module_extension.py','module_coeffs.py','verify.py','inherited.py','carrier.py','skeleton.py')}}
 Path(__file__).with_name('module_extension.receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
