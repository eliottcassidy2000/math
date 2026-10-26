"""Bounded standard-library incoming review; no external census/binary required."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import json, csv, math
ROOT=Path(__file__).resolve().parent

def c3(M):
 return [list(t) for t in combinations(range(len(M)),3) if all(sum(M[i][j] for j in t if j!=i)==1 for i in t)]
def converse(M): return [list(x) for x in zip(*M)]
def sub(Q,children):
 sizes=[len(c) for c in children]; off=[sum(sizes[:i]) for i in range(len(sizes))]; n=sum(sizes); M=[[0]*n for _ in range(n)]
 for a,A in enumerate(children):
  for i in range(len(A)):
   for b,B in enumerate(children):
    for j in range(len(B)):M[off[a]+i][off[b]+j]=A[i][j] if a==b else Q[a][b]
 return M
def isocount(M,N):
 n=len(M); dM=list(map(sum,M));dN=list(map(sum,N)); candidates={i:[j for j in range(n) if dM[i]==dN[j]] for i in range(n)}
 if sorted(dM)!=sorted(dN):return 0
 order=sorted(range(n),key=lambda i:len(candidates[i])); mp={}; used=set()
 def go(k):
  if k==n:return 1
  i=order[k]; tot=0
  for j in candidates[i]:
   if j not in used and all(M[i][x]==N[j][y] for x,y in mp.items()):
    mp[i]=j;used.add(j);tot+=go(k+1);used.remove(j);del mp[i]
  return tot
 return go(0)
def hp(M):
 n=len(M);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if d[mask][i]:
    for j in range(n):
     if not(mask>>j&1) and M[i][j]:d[mask|1<<j][j]+=d[mask][i]
 return sum(d[-1])
C=[[0,1,0],[0,0,1],[1,0,0]]; A=sub([[0,1],[0,0]],[[[0]],C]); T=sub(C,[A,A,A]);U=sub(C,[converse(A),A,A]);V=sub(C,[A,[[0]],[[0]]]);W=converse(V)
witness={'construction':'A=singleton -> C3; T=C3[A,A,A]; U=C3[A^op,A,A]','A_scores':sorted(map(sum,A)),'Aop_scores':sorted(map(sum,converse(A))),'T':T,'U':U,'same_labelled_c3':c3(T)==c3(U),'cyclic_triples':len(c3(T)),'aut_T':isocount(T,T),'aut_U':isocount(U,U),'H_T':hp(T),'H_U':hp(U),'tree_fibre_N':128,'n6_strong_types_counterexample':{'V':V,'W':W,'same_labelled_c3':c3(V)==c3(W),'scores_V':sorted(map(sum,V)),'scores_W':sorted(map(sum,W)),'isomorphisms':isocount(V,W)}}
# Independent exact coefficient evaluation from aggregate table, not producer code.
def ct_parse(s):
 if s=='-':return []
 ans=[]
 for tok in s.split():
  ct,m=tok.rsplit('x',1);ans.append(([(int(l),int(e)) for l,e in (piece.split('^') for piece in ct.split(','))],int(m)))
 return ans
rows=[]
for row in csv.reader((ROOT/'klein-packet/data/prime_cycle_types_n9.tsv').read_text().splitlines()[1:],delimiter='\t'):
 k,sc,aut,ats,anti,count=row;rows.append((int(k),int(sc),int(aut),ct_parse(ats),ct_parse(anti),int(count)))
def mul(a,b,n):return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(n+1)]
def coeff(ct,seq,n):
 p=[1]+[0]*n
 for l,e in ct:
  q=[seq[j//l] if j and j%l==0 else 0 for j in range(n+1)]
  for _ in range(e):p=mul(p,q,n)
 return p[n]
def ci(terms,seq,n):return sum(m*coeff(ct,seq,n) for ct,m in terms)
a=[0,1]+[0]*9;s=a.copy();c=a.copy();cs=a.copy()
for n in range(2,11):
 s[n]=sum(F(count,aut)*ci(ats,a,n) for k,sc,aut,ats,anti,count in rows if k<=n)+(7906161 if n==10 else 0)
 a[n]=s[n]+sum(s[j]*a[n-j] for j in range(1,n))
 cs[n]=sum(F(count,aut*(2 if sc else 1))*F(1 if sc else 1,1 if sc else 2)*ci(ats+anti,c,n) for k,sc,aut,ats,anti,count in rows if k<=n)+(F(7906161+5525,2) if n==10 else 0)
 # Euler transform via logarithmic derivative of product (1-x^j)^(-cs[j]).
 z=[1]+c[1:n]+[0]
 c[n]=F(sum(sum(d*cs[d] for d in range(1,k+1) if k%d==0)*z[n-k] for k in range(1,n+1)),n)
R5=[[int(i!=j and (j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]; reflection=[(-i)%5 for i in range(5)]; inversions=sum(reflection[i]>reflection[j] for i in range(5) for j in range(i+1,5))
witness['A5_bit_counterexample']={'R5':R5,'reflection':reflection,'inversions':inversions,'permutation_even':inversions%2==0,'maps_to_converse':all(R5[i][j]==R5[reflection[j]][reflection[i]] for i in range(5) for j in range(5)),'same_labelled_C3':c3(R5)==c3(converse(R5))}
expected_a=[1,1,2,4,12,56,456,6880,191536,9733056];expected_c=[1,1,2,3,9,32,247,3328,93612,4774417]
receipt={'status':'pass','witness':witness,'species':{'method':'new Fraction polynomial arithmetic from 23 aggregated source cycle-index rows; order10 prime totals treated as input, not recounted','A000568':[int(x) for x in a[1:]],'strong':[int(x) for x in s[1:]],'C3_classes':[int(x) for x in c[1:]],'strong_C3_classes':[int(x) for x in cs[1:]],'matches_A000568':a[1:]==expected_a,'matches_C3':c[1:]==expected_c}}
if not (witness['same_labelled_c3'] and witness['aut_T']==81 and witness['aut_U']==27 and a[1:]==expected_a and c[1:]==expected_c):raise RuntimeError('check failed')
(ROOT/'independent-review.json').write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps({**receipt,'witness':{k:v for k,v in witness.items() if k not in ('T','U')}},indent=2))
