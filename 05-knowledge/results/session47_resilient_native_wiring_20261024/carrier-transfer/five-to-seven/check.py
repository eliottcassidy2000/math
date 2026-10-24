#!/usr/bin/env python3
"""Separate next-size native suspension probe; does not modify frozen 4->6 packet."""
import argparse,itertools,json,math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import native_suspension_probe as base
CHECKS=0

def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def mv(A,x,p):return tuple(sum(a*b for a,b in zip(r,x))%p for r in A)
def echelon(A,b,p):
 A=[[a%p for a in r] for r in A];b=[a%p for a in b];pivots=[];ops=[];t=0
 for col in range(len(A[0])):
  row=next((i for i in range(t,len(A)) if A[i][col]),None)
  if row is None:continue
  if row!=t:A[row],A[t]=A[t],A[row];b[row],b[t]=b[t],b[row];ops.append(('swap',row,t))
  inv=pow(A[t][col],-1,p)
  if inv!=1:A[t]=[a*inv%p for a in A[t]];b[t]=b[t]*inv%p;ops.append(('scale',t,inv))
  for i in range(t+1,len(A)):
   q=A[i][col]
   if not q:continue
   A[i]=[(a-q*c)%p for a,c in zip(A[i],A[t])];b[i]=(b[i]-q*b[t])%p;ops.append(('add',i,t,-q))
  pivots.append(col);t+=1
  if t==len(A):break
 return A,b,pivots,ops

def kernel_basis(A,p):
 B,b,pivots,ops=echelon(A,[0]*len(A),p);free=[j for j in range(len(A[0])) if j not in pivots];out=[]
 for j in free:
  x=[0]*len(A[0]);x[j]=1
  for i,c in reversed(list(enumerate(pivots))):x[c]=-sum(a*z for a,z in zip(B[i],x))%p
  need(not any(mv(A,x,p)),'old nullbasis literal test');out.append(x)
 return out,free

def solve(A,b,p):
 B,c,pivots,ops=echelon(A,b,p);bad=next((i for i in range(len(pivots),len(A)) if c[i]),None)
 if bad is not None:
  z=[0]*len(A);z[bad]=1
  for op in reversed(ops):
   if op[0]=='swap':_,i,j=op;z[i],z[j]=z[j],z[i]
   elif op[0]=='scale':_,i,q=op;z[i]=z[i]*q%p
   else:_,i,j,q=op;z[j]=(z[j]+q*z[i])%p
  need(all(sum(z[i]*A[i][j] for i in range(len(A)))%p==0 for j in range(len(A[0]))),'left obstruction annihilates full constraint matrix')
  value=sum(a*c for a,c in zip(z,b))%p;need(value!=0,'left obstruction detects requested readout')
  return None,{'dual':z,'pairing':value,'rank':len(pivots)}
 x=[0]*len(A[0])
 for i,j in reversed(list(enumerate(pivots))):x[j]=(c[i]-sum(a*y for a,y in zip(B[i],x)))%p
 need(mv(A,x,p)==tuple(a%p for a in b),'full native extension verified')
 return x,{'rank':len(pivots)}

def integer_determinant(A):
 B=[r[:] for r in A];previous=1
 for k in range(len(A)-1):
  pivot=B[k][k];need(pivot!=0,'positive definite old leading principal pivot')
  for i in range(k+1,len(A)):
   for j in range(k+1,len(A)):
    numerator=pivot*B[i][j]-B[i][k]*B[k][j]
    need(numerator%previous==0,'Bareiss exact division')
    B[i][j]=numerator//previous
  for i in range(k+1,len(A)):B[i][k]=0
  previous=pivot
 return B[-1][-1]

def field_determinant(A,p):
 B,b,pivots,ops=echelon(A,[0]*len(A),p);u=1
 for op in ops:
  if op[0]=='swap':u=-u%p
  elif op[0]=='scale':u=u*op[2]%p
 return len(pivots),pow(u,-1,p) if len(pivots)==len(A) else 0

def main():
 old=base.build(5);new=base.build(7);idx={v:i for i,v in enumerate(new['keys'])};readouts=[];sources={}
 for w,u in enumerate(old['owner']):
  B=base.matrix(5,w);S=base.suspend(B);K=base.partner(S,(5,0,1,2,3,4,6));v=old['owner'][63^w]
  ns=idx[base.merged(S)];nt=idx[base.merged(K)];sources.setdefault(u,set()).add(ns);readouts.extend(((ns,u),(nt,v)))
  need([r[:5] for r in K[:5]]==base.partner(B,(0,1,2,3,4)),'literal n5 inner target restriction')
 readouts=sorted(set(readouts));root=min(sources[0]);coords=[i for i in range(len(new['keys'])) if i!=root]
 A=[r[1:] for r in old['L'][1:]];N=[[new['L'][i][j] for j in coords] for i in coords]
 R=[[int(j==node) for j in coords] for node,owner in readouts];S=[[int(j==owner) for j in range(1,len(old['keys']))] for node,owner in readouts]
 M=N+R;J=[[0]*len(A) for _ in N]+S;tests=[];lifts2=None;free2=None
 olddet=integer_determinant(A);need(olddet==41096==8*11*467,'old exact determinant factorization')
 odd_determinants=[]
 for prime,expected in ((11,8),(467,397)):
  rank,det=field_determinant(N,prime)
  need(rank==len(N) and det==expected,'new full Laplacian excludes old odd-primary direction')
  odd_determinants.append({'prime':prime,'new_rank':rank,'new_size':len(N),'new_determinant_mod_prime':det})
 for p in (2,11,467):
  basis,free=kernel_basis(A,p);lifts=[];failed=[]
  for h in basis:
   y,witness=solve(M,mv(J,h,p),p)
   if y is None:failed.append({'old_vector':h,**witness})
   else:lifts.append({'old_vector':h,'new_vector':y,**witness})
  tests.append({'prime':p,'old_nullity':len(basis),'old_free_columns':free,'successful_generators':lifts,'failures':failed})
  if p==2 and not failed:lifts2=[a['new_vector'] for a in lifts];free2=free
 certificate=None
 if lifts2 is not None:
  from fractions import Fraction
  # Exact rational inverse of the small original reduced matrix.
  n=len(A);aug=[[Fraction(a) for a in r]+[Fraction(int(i==j)) for j in range(n)] for i,r in enumerate(A)]
  for col in range(n):
   row=next(i for i in range(col,n) if aug[i][col]);aug[col],aug[row]=aug[row],aug[col];v=aug[col][col];aug[col]=[a/v for a in aug[col]]
   for i in range(n):
    if i==col:continue
    q=aug[i][col]
    if q:aug[i]=[a-q*b for a,b in zip(aug[i],aug[col])]
  inv=[r[n:] for r in aug];Z=[[0]*n for _ in coords]
  for f,y in zip(free2,lifts2):
   for i,a in enumerate(y):Z[i][f]=a
  defect=[[sum(M[i][k]*Z[k][j] for k in range(len(coords)))-J[i][j] for j in range(n)] for i in range(len(M))]
  rational=[[sum(Fraction(defect[i][k])*inv[k][j] for k in range(n)) for j in range(n)] for i in range(len(M))]
  d=math.lcm(*(a.denominator for r in rational for a in r));need(d%2==1,'all dyadic extension certificate has odd denominator')
  Q=[[int(d*a) for a in r] for r in rational]
  need([[sum(Q[i][k]*A[k][j] for k in range(n)) for j in range(n)] for i in range(len(M))]==[[d*a for a in r] for r in defect],'exact odd-multiplier certificate')
  certificate={'identity':'d (M Z - J) = Q A over integers','odd_multiplier':d,'Z':Z,'Q':Q}
 # Tiny exact full-Laplacian obstruction, independent of Gaussian bookkeeping.
 rowset=[10,41,196,244]
 row_sum=[sum(new['L'][i][j] for i in rowset)%2 for j in range(len(new['L']))]
 need(row_sum==[int(j in (10,41)) for j in range(len(new['L']))],'four physical rows force equality of new nodes10,41')
 h=[int(i in (3,4,6,7)) for i in range(len(old['keys']))]
 need(not any(mv(old['L'],h,2)) and h[3]!=h[9],'old full harmonic distinguishes3,9')
 source_records=[]
 for node,word,target in ((3,55,10),(9,0,41)):
  B=base.matrix(5,word);suspended=base.suspend(B)
  need(old['owner'][word]==node and idx[base.merged(suspended)]==target,'literal source suspension obstruction readout')
  source_records.append({'old_node':node,'word':word,'old_matrix':B,'suspended_matrix':suspended,'new_node':target})
 dyadic=tests[0];need(len(dyadic['successful_generators'])==2 and len(dyadic['failures'])==1,'exact two-dimensional successful subspace')
 successful=[]
 for bits in itertools.product((0,1),repeat=2):
  h=tuple(sum(b*v['old_vector'][i] for b,v in zip(bits,dyadic['successful_generators']))%2 for i in range(len(A)))
  y=tuple(sum(b*v['new_vector'][i] for b,v in zip(bits,dyadic['successful_generators']))%2 for i in range(len(coords)))
  need(mv(M,y,2)==mv(J,h,2),'all four successful old harmonics have verified extensions')
  successful.append(h)
 need(len(set(successful))==4,'successful span has dimension two')
 tiny={'new_laplacian_rows':rowset,'new_relation_support':[10,41],'old_harmonic_support':[3,4,6,7],'old_relation_support':[3,9],'source_records':source_records,'old_kernel_dimension':3,'extendible_dimension':2,'extendible_old_vectors':successful,'total_old_harmonics':8,'extendible_old_harmonics':4}
 collisions={str(node):sorted({owner for nn,owner in readouts if nn==node}) for node in {nn for nn,o in readouts} if len({owner for nn,owner in readouts if nn==node})>1}
 return {'status':'FINITE-EXACT native5-to7 suspension extension; not8-to10','old':old,'new':new,'new_root':root,'new_reduced_coordinates':coords,'readouts':readouts,'record_collisions':collisions,'tests':tests,'dyadic_integer_certificate':certificate,'tiny_obstruction':tiny,'old_determinant':olddet,'old_smith_factors':[1,1,1,1,1,1,2,2,10274],'odd_determinants':odd_determinants,'checks':CHECKS,'base_checks':base.CHECKS}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path);a=p.parse_args();result=main();text=json.dumps(result,sort_keys=True,indent=2)+'\n'
 if a.out:a.out.write_text(text)
 else:print(text,end='')
