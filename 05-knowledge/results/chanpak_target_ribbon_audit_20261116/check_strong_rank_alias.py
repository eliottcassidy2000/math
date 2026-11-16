"""Same-size strong obstruction, after the independently replayed rational quotient."""
from verify_rank_alias import product,O,B,h,v0,rank,piv,rank_pair,Matrix
from cap_check import need
import json
from pathlib import Path

def brute_hamiltonian_count(w):
 letters=(0,)+tuple(w); n=len(letters)
 need(all(a < (1 << min(i,3)) for i,a in enumerate(letters)))
 out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   backward=j-i<=3 and bool(letters[j]&(1<<(j-i-1)))
   if backward:out[j]|=1<<i
   else:out[i]|=1<<j
 reach=out.copy()
 for i in range(n):reach[i]|=1<<i
 for k in range(n):
  for i in range(n):
   if reach[i]&(1<<k):reach[i]|=reach[k]
 need(all(x==(1<<n)-1 for x in reach))
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   rest=out[i]&~mask
   while rest:
    bit=rest&-rest;rest-=bit;j=bit.bit_length()-1
    dp[mask|bit][j]+=dp[mask][i]
 return sum(dp[-1])
u=(0,2,7,0,2);v=(0,1,7,0,2);tail=1
result={'alphabet':'bit mask of backward gaps 1,2,3','u':u,'v':v,'common_suffix':u[-3:],'appended_letter':tail,'minimal_dimension':rank_pair,'ranks':{}}
for name,w in [('u',u),('v',v),('u1',u+(tail,)),('v1',v+(tail,))]:
 P=product(w)
 K=O*P*B
 actual=int((Matrix([list(h)])*P*Matrix(v0))[0])
 need(actual==brute_hamiltonian_count(w))
 result[name+'_H']=actual
 r=rank(K);result['ranks'][name]=r
 print(name,w,'rank',r,'H',actual,flush=True)
 # Retain compact explicit nonsingular minor and exact rank factorization.
 cs=piv(K);C=K[:,list(cs)];rs=piv(C.T)
 D=C[list(rs),:]
 need(D.det()!=0)
 F=D.inv()*K[list(rs),:]
 need(C*F==K)
 result[name+'_minor']={'rows':list(rs),'cols':list(cs),'det':str(D.det())}
need(result['ranks']=={'u':3,'v':3,'u1':3,'v1':2})
need(json.dumps(result,indent=2)+'\n'==Path(__file__).with_name('STRONG_RANK_ALIAS.json').read_text())
print('PASS: same last three letters and same rank do not determine next rank.',flush=True)
