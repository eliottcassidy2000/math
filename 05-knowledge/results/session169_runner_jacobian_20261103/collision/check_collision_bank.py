#!/usr/bin/env python3
"""Exact, standard-library controls for S169 full-query collision banks."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json

P=1009
C=((0,1,0),(1,1,1),(1,0,0))
D=((1,0,1),(0,2,0),(0,0,0))
I=((1,0,0),(0,1,0),(0,0,1))
V=((1,2,3),(2,3,6),(0,1,2))
BASIS=('', 'C','D','CC','CD','DC','DD','CCD','CDD')
BANK=('CDDD','DD','CDCD','CCD','CDC','DDCC','DCD','DCDC','CDD','CCC','CCCC','C')
checks=0

def require(b,msg):
 global checks
 checks+=1
 if not b:raise RuntimeError(msg)

def mm(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))) for i in range(len(A)))
def flat(A):return tuple(x for row in A for x in row)
def mat(word):
 M=I
 for letter in word:M=mm(C if letter=='C' else D,M)
 return M

def rankmod(A,p=P):
 A=[[int(x)%p for x in row] for row in A];r=0
 if not A:return 0
 for j in range(len(A[0])):
  i=next((i for i in range(r,len(A)) if A[i][j]),None)
  if i is None:continue
  A[r],A[i]=A[i],A[r];z=pow(A[r][j],-1,p);A[r]=[x*z%p for x in A[r]]
  for i in range(r+1,len(A)):
   z=A[i][j]
   if z:A[i]=[(a-z*b)%p for a,b in zip(A[i],A[r])]
  r+=1
  if r==len(A):break
 return r

def inverse(A):
 n=len(A);A=[[F(x) for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n) if A[i][j]);A[j],A[i]=A[i],A[j]
  z=A[j][j];A[j]=[x/z for x in A[j]]
  for i in range(n):
   if i!=j and A[i][j]:
    z=A[i][j];A[i]=[a-z*b for a,b in zip(A[i],A[j])]
 return [row[n:]for row in A]

def det(A):
 A=[[F(x)for x in row]for row in A];d=F(1)
 for j in range(len(A)):
  i=next((i for i in range(j,len(A))if A[i][j]),None)
  if i is None:return F(0)
  if i!=j:A[j],A[i]=A[i],A[j];d=-d
  z=A[j][j];d*=z;A[j]=[a/z for a in A[j]]
  for i in range(j+1,len(A)):
   z=A[i][j];A[i]=[a-z*b for a,b in zip(A[i],A[j])]
 return d

def mod(x):
 x=F(x);return x.numerator*pow(x.denominator,-1,P)%P

def coordinate_rows(t):
 n=len(t);pairs=[(i,j)for i in range(n)for j in range(n)if i!=j];ix={v:k for k,v in enumerate(pairs)}
 rows=[]
 for i,j,k in combinations(range(n),3):
  row=[F(0)]*len(pairs)
  for pair,z in [((j,k),t[i]*t[j]),((i,k),-t[i]*t[j]),((i,j),t[i]*t[k]),((k,j),-t[i]*t[k]),((k,i),t[j]*t[k]),((j,i),-t[j]*t[k])]:row[ix[pair]]=z
  rows.append(row)
 return rows

def common(n):return [[int(j==k)for i in range(n)for j in range(n)if i!=j]for k in range(n)]

def generic_tests():
 out=[]
 for n in range(3,10):
  extra=[tuple(1 for i in range(n)),tuple(i+1 for i in range(n)),tuple((i+1)**2 for i in range(n))]
  rows=[r for t in extra for r in coordinate_rows(t)]
  rk=rankmod(rows);want=n*(n-2)
  require(rk==want,'generic minimal bank')
  require(rankmod(common(n))==n,'common motion injective')
  out.append({'dimension':n,'directions':n+3,'normalized_constraint_rank':rk,'common_kernel':n})
 return out

def hp(rows,start=None,deleted=None):
 n=len(rows);mask=(1<<n)-1
 if deleted is not None:mask^=1<<deleted
 @lru_cache(None)
 def tail(m,v):
  if m==(1<<v):return 1
  rem=m^(1<<v)
  return sum(tail(rem,u)for u in range(n)if rem>>u&1 and rows[v]>>u&1)
 if start is not None:return tail(mask,start)
 return sum(tail(mask,v)for v in range(n)if mask>>v&1)

def strong(rows):
 n=len(rows)
 for root in range(n):
  seen=1<<root;old=0
  while old!=seen:
   old=seen
   for i in range(n):
    if seen>>i&1:seen|=rows[i]
  if seen!=(1<<n)-1:return False
 return True

def physical(word):
 rows=[2,4,1];q=0
 for letter in word:
  n=len(rows)
  rows=[r|((1<<n)if i==q else 0)for i,r in enumerate(rows)]+[((1<<n)-1)^(1<<q)]
  if letter=='C':q=n
 return rows,q

def query(rows,p):
 n=len(rows)
 return [r|((1<<n)if i!=p else 0)for i,r in enumerate(rows)]+[1<<p],n

def state(rows,q):
 pre=[i for i in range(len(rows))if rows[i]>>q&1]
 require(len(pre)==1,'unique in-neighbor')
 s=hp(rows,start=q);h=hp(rows)
 return s,h-s,hp(rows,start=pre[0],deleted=q)

def literal_tests():
 out=[]
 for word in sorted(set(BASIS+BANK),key=lambda x:(len(x),x)):
  rows,q=physical(word);cols=[];p=1
  for k in range(3):
   require(strong(rows),'strong literal query graph')
   cols.append(state(rows,q))
   if k<2:rows,p=query(rows,p)
  observed=tuple(zip(*cols));expected=mm(mat(word),V)
  require(observed==expected,'nine-state raw-query transfer')
  out.append({'word':word or 'empty','vertices':3+len(word),'state_row_major':flat(observed)})
 return out

def physical_bank_test(words):
 basis=[];extra=[]
 for word in words:
  t=flat(mm(mat(word),V))
  if len(basis)<9 and rankmod([v for _,v in basis]+[t])>len(basis):basis.append((word,t))
  else:extra.append((word,t))
 require(len(basis)==9,'physical spanning bank')
 U=list(zip(*[v for _,v in basis]));Ui=inverse(U)
 rows=[];relations=[]
 ordered=basis+extra
 for extra_index,(word,t) in enumerate(extra):
  a=[sum(x*y for x,y in zip(row,t))for row in Ui]
  rows.extend(coordinate_rows(a))
  relations.append([-x for x in a]+[F(j==extra_index)for j in range(len(extra))])
 period=[]
 for relation in relations:
  for i,j in combinations(range(9),2):
   row=[F(0)]*(9*len(ordered))
   for k,(_,t) in enumerate(ordered):
    row[9*k+j]=relation[k]*t[i];row[9*k+i]=-relation[k]*t[j]
   period.append([mod(x)for x in row])
 period_rank=rankmod(period)
 require(period_rank==87,'independent original wedge-period rank')
 rk=rankmod([[mod(x)for x in row]for row in rows])
 require(rk==63,'physical complete collision bank')
 return {'words':[w or 'empty'for w in words],'basis_words':[w or 'empty'for w,_ in basis],'physical_basis_determinant':str(det(U)),'normalized_constraint_rank_mod_1009':rk,'common_motion_kernel_dimension':9,'max_moves':max(map(len,words)),'direct_period_matrix_shape':[len(period),len(period[0])],'direct_period_matrix_rank_mod_1009':period_rank}

def radius3_test():
 U=list(zip(*[flat(mat(w))for w in BASIS]));Ui=inverse(U)
 require(det(U)==16,'matrix-coordinate basis determinant')
 pairs=[(i,j)for i in range(9)for j in range(9)if i!=j];ix={v:k for k,v in enumerate(pairs)}
 H=[[0]*72 for _ in range(3)]
 H[0][ix[4,5]]=H[0][ix[5,4]]=1
 for j,z in enumerate((-1,-1,-1,2,1,1,1,0,0)):
  if j!=7:H[1][ix[7,j]]=z
 for j,z in enumerate((1,1,1,0,-2,-2,-1,0,0)):
  if j!=8:H[2][ix[8,j]]=z
 words=[''.join(x)for k in range(4)for x in product('CD',repeat=k)]
 rows=[]
 for word in words:
  a=[sum(x*y for x,y in zip(row,flat(mat(word))))for row in Ui]
  rr=coordinate_rows(a);rows.extend(rr)
  for h in H:
   require(all(sum(x*y for x,y in zip(row,h))==0 for row in rr),'exact radius-three hostile')
 rk=rankmod([[mod(x)for x in row]for row in rows])
 require(rk==60,'radius-three rank')
 require(rankmod(common(9)+H)==12,'three independent extra classes')
 require(rankmod([flat(mat(w))for w in words if len(w)<=2])==7,'radius-two span')
 return {'word_count':len(words),'basis_determinant_matrix_coordinates':16,'normalized_constraint_rank_mod_1009':rk,'kernel_dimension_exact':12,'blind_quotient_dimension_exact':3,'radius_two_span':7}

def main():
 require(det(V)==-2,'seed raw-query rank')
 result={'generic_dimension_controls':generic_tests(),'physical_bank':physical_bank_test(BANK),'radius_three':radius3_test(),'literal_states':literal_tests()}
 result['checks']=checks
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
