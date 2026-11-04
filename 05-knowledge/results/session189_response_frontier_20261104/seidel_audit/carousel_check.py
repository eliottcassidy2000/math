from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json
D=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
def bareiss(A):
 n=len(A)
 if not n:return 1
 M=[r.copy() for r in A];sign=1;previous=1
 for k in range(n-1):
  pivot=next((i for i in range(k,n) if M[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:M[k],M[pivot]=M[pivot],M[k];sign=-sign
  p=M[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    top=p*M[i][j]-M[i][k]*M[k][j];need(top%previous==0,'exact Bareiss division');M[i][j]=top//previous
  for i in range(k+1,n):M[i][k]=0
  previous=p
 return sign*M[-1][-1]

def main():
 records=[]
 for n in range(3,16,2):
  k=(n-1)//2;Q=[sum(1<<((i+d)%n) for d in range(1,k+1)) for i in range(n)];rest=[i for i in range(n) if i not in (0,k)];total=Fraction(0);terms=[]
  for size in range(0,len(rest)+1,2):
   for chosen in combinations(rest,size):
    ids=sorted((0,k)+chosen);A=[[0 if i==j else 1 if Q[i]>>j&1 else -1 for j in ids] for i in ids];d=bareiss(A);need(d==1,'carousel even principal Seidel determinant is one');rs=[]
    for label in (0,k):
     col=ids.index(label);N=[[1 if j==col else v for j,v in enumerate(row)] for row in A];rs.append(Fraction(bareiss(N),d))
    term=rs[0]**2*rs[1];total+=term;terms.append([sum(1<<i for i in ids),str(rs[0]),str(rs[1]),str(term)])
  coefficient=total*Fraction(2**max(6-n,0),2**max(n-6,0));need(coefficient==-8*(n-4),'all-subset independent exact carousel coefficient');need(len(terms)==2**(n-3),'all even subsets containing the two fixed names');records.append(dict(order=n,rows=Q,shared_names=[0,k,0],subset_terms=terms,cubic_coefficient=str(coefficient)))
 (D/'carousel_data.json').write_text(json.dumps(records,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',orders=[r['order'] for r in records],exact_subsets=sum(len(r['subset_terms']) for r in records),coefficients=[r['cubic_coefficient'] for r in records]),sort_keys=True,indent=2))
if __name__=='__main__':main()
