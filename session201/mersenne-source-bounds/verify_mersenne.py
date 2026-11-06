"""Exact controls for MERSENNE.md; not a global occupancy census."""
from pathlib import Path
from functools import lru_cache
import hashlib,json
ROOT=Path(__file__).resolve().parent
MENU=json.loads((ROOT/'dependencies/TWELVE.json').read_text())
def need(c,m):
 if not c:raise RuntimeError(m)
@lru_cache(None)
def g(n):
 need(n>=1,'positive index')
 if n==1:return 1
 if n%2==0:return g(n//2)
 h=(n+1)//2;return g(h)+h-h%2

def g_bits(n):
 n-=1;ans=1
 while n:
  if n%2==0:ans+=n//2+(n//2)%2
  n//=2
 return ans

def bilinear(M,u,v,A,B):
 t=(A,u%2,u,1);s=(B,v%2,v,1)
 return sum(t[i]*M[i][j]*s[j] for i in range(4) for j in range(4))

def main():
 forms=[]
 for e in MENU:
  for transpose in (False,True):
   M=list(map(list,zip(*e['matrix']))) if transpose else e['matrix']
   if not any(M==x[1] for x in forms):forms.append((e['word']+('*' if transpose else ''),M))
 need(len(forms)==15,'transpose closure has 15 matrices')
 # Finite exact coefficient certificate, applying to every actual source.
 for w,M in forms:
  if not w:continue
  lam,c,d,b=M[0][0],M[0][2],M[2][0],M[2][2]
  e,f,h,i,j=-M[0][1],-M[1][0],-M[1][2],-M[2][1],M[1][1]
  need(lam in (1,3,5) and 0<=c<=7 and 0<=d<=7 and 0<=b<=8,'core coefficients')
  need(min(e,f,h,i)>=0 and j in (-1,0,1),'parity coefficients')
  need(lam+2*c-e>=1 and lam+2*d-f>=1,'both fixed-source g coefficients positive')
  need(all(M[3][k]==M[k][3]==0 for k in range(4)),'no nonempty affine constants')
 # Recheck the all-source physical product bound by its cone certificate.
 P=[[1,0,2,1],[0,0,1,0],[1,1,2,2],[0,0,1,1]];weights=[1,1,2,2]
 for w,M in forms:
  for i in range(4):
   for j in range(4):
    C=sum(P[k][i]*M[k][l]*P[l][j] for k in range(4) for l in range(4))
    need(2*C>=weights[i]*weights[j],'physical uv/2 coefficient bound')
 for n in range(1,8193):
  A=g(n);need(A==g_bits(n),'independent g implementations')
  need(A%2==1 and 1<=A<=n,'g elementary bounds')
  a=n
  while a%2==0:a//=2
  need(g(a)==A and a<=2*A,'odd core bound')
  if a>=3:need(2*A>=a+3,'odd core lower bound')
 for j in range(1,31):need(g((1<<j)+1)==(1<<j)+1,'binary-plus-one endpoint')
 controls=[];total=0
 for G in (1,3,5,9):
  # Under the proof, all Mersenne witnesses with this G have r,s <= floor(log2(90G^4)).
  cap=(90*G**4).bit_length()-1
  cores=[a for a in range(1,2*G+1,2) if g(a)<=G]
  indices=[(a,r,g(a),a<<r) for a in cores for r in range(cap+1) if a<<r>=2]
  vals={};checks=0
  for a,r,A,u in indices:
   for b0,s,B,v in indices:
    if r<s:continue
    for w,M in forms:
     n=bilinear(M,u,v,A,B);checks+=1
     if (n+1)&n:continue
     need(n<=18630*G**8,'Mersenne polynomial ceiling')
     vals.setdefault(n,[w,u,v,A,B])
     k=(n+1).bit_length()-1
     if r==0:continue
     need(k>r,'larger source height below target height')
     if w:
      lam,c=M[0][0],M[0][2]
      if s:
       need((lam*A*B+1)%(1<<s)==0,'small height nonempty divisor')
       need(1<<s<=6*G*G and v<=12*G**3,'small height nonempty bound')
       C=A*(lam*B+c*v)
       need(C>0 and (C+1)%(1<<r)==0,'large height nonempty divisor')
       need(1<<r<=90*G**4 and u<=180*G**5,'large height nonempty bound')
      else:
       C=A*(lam*B+c*v+M[0][1])
       need(C>0 and (C+1)%(1<<r)==0 and 1<<r<=20*G**2,'one odd nonempty bound')
     elif s:
      Z=A+B-v+2
      need(Z!=0 and Z%(1<<r)==0,'empty even/even zero boundary and divisor')
      need((A+B+2)%(1<<s)==0 and v<=8*G**2,'empty small-height bound')
      need(1<<r<=12*G**2 and u<=24*G**3,'empty large-height bound')
     else:
      Z=A+B-v+1
      need(Z!=0 and Z%(1<<r)==0 and 1<<r<=5*G,'empty even/odd zero boundary and divisor')
  total+=checks
  controls.append({'G':G,'odd_cores':cores,'maximum_height_tested':cap,'source_indices':len(indices),'matrix_pair_checks':checks,'Mersenne_values':sorted(vals),'largest_Mersenne_witness':vals[max(vals)]})
 need(31 in controls[0]['Mersenne_values'],'positive G=1 hostile to overstrong exclusion')
 # Deliberately realize the algebraic zero-C branch at a non-Mersenne value.
 # Empty u=2, v=4: A=B=1 gives Z=0 but n+1=6, which has odd factor 3.
 need(bilinear(MENU[0]['matrix'],2,4,1,1)==5,'zero constant can occur outside Mersenne targets')
 return {'status':'PASS','finite_g_controls_through':8192,'formal_transpose_matrices':len(forms),'matrix_pair_checks':total,'controls':controls,'empty_zero_constant_non_Mersenne_control':{'u':2,'v':4,'n':5},'source_menu_sha256':hashlib.sha256((ROOT/'dependencies/TWELVE.json').read_bytes()).hexdigest(),'claim_scope':'Only Mersenne targets represented with both source g coordinates <= G','uniform_occupancy':'OPEN'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
