"""Independent reviewer: no imports from author implementation."""
import itertools as I,math,json,random
from fractions import Fraction

def need(c,m):
 if not c:raise RuntimeError(m)

def tournament(n,z):
 a=[[False]*n for _ in range(n)]
 for b,(i,j) in enumerate(I.combinations(range(n),2)):a[i][j]=bool(z>>b&1);a[j][i]=not a[i][j]
 return a

def qprofile(a):
 n=len(a);q=[0]*n
 for p in I.permutations(range(n)):
  forced=sum(not a[x][y] for x,y in zip(p,p[1:]))
  for cuts in range(forced,n):q[cuts]+=math.comb(n-1-forced,cuts-forced)
 return q

def context(a,r,out):
 n=len(a);m=len(r);b=[row[:]+[bool(out>>j&1) for j in range(m)] for row in a]
 for j in range(m):b.append([not bool(out>>j&1)]*n+r[j][:])
 return b

def coeff(r,out):
 m=len(r);ans=[0]*(m//2+1)
 def edge(x,y):
  if x==-1:return bool(out>>y&1)
  if y==-1:return not bool(out>>x&1)
  return r[x][y]
 for p in I.permutations(range(m)):
  for k in range(1,len(ans)+1):
   for gaps in I.combinations(range(m+1),k):
    w=[]
    for j in range(m+1):
     if j in gaps:w.append(-1)
     if j<m:w.append(p[j])
    if all(edge(x,y) for x,y in zip(w,w[1:])):ans[k-1]+=1
 return ans

def rank(rows):
 piv={}
 for row in rows:
  row=list(map(Fraction,row))
  for col in range(len(row)):
   if not row[col]:continue
   if col in piv:
    c=row[col];row=[x-c*y for x,y in zip(row,piv[col])]
   else:
    c=row[col];piv[col]=[x/c for x in row];break
 return len(piv)

def main():
 rng=random.Random(17699);ranks=[];profile_checks=0
 for n in range(1,6):
  evals=[]
  for z in range(1<<(n*(n-1)//2)):
   q=qprofile(tournament(n,z));profile_checks+=1
   evals.append([sum(v*math.comb(t,k) for k,v in enumerate(q,1) if k<=t) for t in range(1,n+1)])
  ranks.append(rank(evals))
 need(ranks==[1,1,2,2,3],'realized ranks')
 contexts=0
 for n,m in [(2,4),(3,3),(4,3),(5,2)]:
  for _ in range(12):
   a=tournament(n,rng.randrange(1<<(n*(n-1)//2)));r=tournament(m,rng.randrange(1<<(m*(m-1)//2)));out=rng.randrange(1<<m)
   q=qprofile(a);c=coeff(r,out)
   need(qprofile(context(a,r,out))[0]==sum(x*y for x,y in zip(q,c)),'arbitrary actual context');contexts+=1
 # Reversal gives equal q; two nested actual contexts must preserve complete q.
 nested=0
 for _ in range(12):
  a=tournament(3,rng.randrange(8));b=[[a[j][i] for j in range(3)] for i in range(3)]
  need(qprofile(a)==qprofile(b),'reverse profile')
  r=tournament(2,rng.randrange(2));s=tournament(2,rng.randrange(2));o=rng.randrange(4);p=rng.randrange(4)
  aa=context(context(a,r,o),s,p);bb=context(context(b,r,o),s,p)
  need(qprofile(aa)==qprofile(bb),'nested profile congruence');nested+=1
 kernels=[]
 for r in range(5):
  n=2*r;tr=[[i<j for j in range(n)] for i in range(n)];out=sum(1<<i for i in range(0,n,2));c=coeff(tr,out);kernels.append(c)
  need(c[-1]==math.factorial(r),'probe diagonal')
 need(kernels==[[1],[3,1],[9,10,2],[31,83,50,6],[117,684,908,336,24]],'probe rows')
 print(json.dumps({'status':'PASS','method':'independent literal forced-cut profiles; actual random contexts; nested substitutions; evaluation-basis rank','labelled_profile_census':profile_checks,'realized_ranks':ranks,'arbitrary_uniform_contexts':contexts,'nested_profile_comparisons':nested,'literal_probe_rows':kernels},sort_keys=True,indent=2))
if __name__=='__main__':main()
