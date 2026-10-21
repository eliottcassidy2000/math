#!/usr/bin/env python3
"""All128 marked q2,m5 frames; exact subset-DP and permutation-search audit."""
import itertools,json,pathlib,hashlib,time,datetime
ROOT=pathlib.Path(__file__).resolve().parent
checks=0
def need(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)
def circ(m,mask):
 D={i if mask>>(i-1)&1 else m-i for i in range(1,(m+1)//2)}
 return [[int(i!=j and (j-i)%m in D) for j in range(m)] for i in range(m)]
def frame(ma,mb,bits):
 A=[[0]*10 for _ in range(10)];C=circ(5,ma);D=circ(5,mb)
 for j in range(5):
  for k in range(5):
   A[2*j][2*k]=C[j][k];A[2*j+1][2*k+1]=D[j][k]
   z=bits[(k-j)%5];A[2*j][2*k+1]=z;A[2*k+1][2*j]=1-z
 return A
def maps(A,anti=False):
 n=len(A);scores=list(map(sum,A));p=[-1]*n;out=[]
 candidates=[[j for j in range(n) if scores[j]==(n-1-scores[i] if anti else scores[i])] for i in range(n)]
 def rec(i,used):
  if i==n:out.append(tuple(p));return
  for j in candidates[i]:
   if used>>j&1:continue
   if all(A[i][k]==(A[p[k]][j] if anti else A[j][p[k]]) for k in range(i)):
    p[i]=j;rec(i+1,used|1<<j)
 rec(0,0);return out
def shore_product_maps(A,anti=False):
 # Independent completeness route: scores force preservation (aut) or
 # interchange (anti) of the two shores; exhaust all120 maps per shore.
 shores=[list(range(0,10,2)),list(range(1,10,2))];choices=[]
 for side,src in enumerate(shores):
  target=shores[1-side if anti else side];valid=[]
  for image in itertools.permutations(target):
   p=dict(zip(src,image))
   if all(A[x][y]==(A[p[y]][p[x]] if anti else A[p[x]][p[y]]) for x in src for y in src):valid.append(p)
  choices.append(valid)
 out=[]
 for left in choices[0]:
  for right in choices[1]:
   p=left|right
   if all(A[x][y]==(A[p[y]][p[x]] if anti else A[p[x]][p[y]]) for x in shores[0] for y in shores[1]):out.append(tuple(p[i] for i in range(10)))
 return out
def dp(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S in range(1,1<<n):
  for v,x in enumerate(D[S]):
   if x:
    for w in range(n):
     if not(S>>w&1) and A[v][w]:D[S|1<<w][w]+=x
 return sum(D[-1])
def words(A,antis):
 n=len(A);H=B=0;first_blue=None
 def rec(p,unused):
  nonlocal H,B,first_blue
  if not unused:
   H+=1;r=[None]*n
   for i in range(n):r[p[i]]=p[n-1-i]
   if tuple(r) in antis:B+=1;first_blue=first_blue or p[:]
   return
  for w in unused:
   if not p or A[p[-1]][w]:rec(p+[w],[v for v in unused if v!=w])
 rec([],list(range(n)));return H,B,first_blue
def blue_by_antis(A,antis):
 n=len(A);total=0
 for a in antis:
  if any(a[a[i]]!=i or a[i]==i for i in range(n)):continue
  pairs=[(i,a[i]) for i in range(n) if i<a[i]]
  for flips in itertools.product((0,1),repeat=len(pairs)):
   for half in itertools.permutations([pair[f] for pair,f in zip(pairs,flips)]):
    P=half+tuple(a[x] for x in reversed(half))
    total+=all(A[x][y] for x,y in zip(P,P[1:]))
 return total
start=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.monotonic();rows=[]
for ma in range(4):
 for mb in range(4):
  for free in itertools.product((0,1),repeat=3):
   bits=(1,)+free+(0,);A=frame(ma,mb,bits)
   need(all(A[i][(i+1)%10] for i in range(10)),'marked Hamilton cycle')
   need(all(A[i][j]==A[(i+2)%10][(j+2)%10] for i in range(10) for j in range(10)),'step2 automorphism')
   au=maps(A);anti=maps(A,True)
   need(set(au)==set(shore_product_maps(A)),'independent automorphism maps')
   need(set(anti)==set(shore_product_maps(A,True)),'independent anti maps')
   H=dp(A);H2,B,blue=words(A,set(anti))
   need(H==H2,'independent HP counts')
   need(B==blue_by_antis(A,anti),'independent blue counts')
   need(H%len(au)==0 and B%len(au)==0,'free word orbit divisibility')
   h=H//len(au);b=B//len(au)
   need(h%2==1,'odd path orbits')
   if anti:need(b%2==1 and (h-b)%2==0,'self converse parity')
   rows.append({'frame':[ma,mb,list(bits)],'source_rows':[''.join(map(str,r)) for r in A],'H':H,'aut_count':len(au),'anti_count':len(anti),'h':h,'b':b,'delta':h-b,'blue_word':blue,'automorphisms':au,'anti_automorphisms':anti})
sc=[r for r in rows if r['anti_count']];surv=[r for r in sc if r['delta']==6]
summary={'frames':len(rows),'self_converse_frames':len(sc),'deficit_six_frames':len(surv),'minimum_self_converse_deficit':min(r['delta'] for r in sc),'self_converse_deficit_histogram':{str(d):sum(r['delta']==d for r in sc) for d in sorted({r['delta'] for r in sc})}}
receipt={'status':'FINITE-EXACT complete marked q2,m5 universe, not all order10 tournaments','start_utc':start,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.monotonic()-t,'checks':checks,'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'summary':summary,'frames':rows}
(ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='frames'},indent=2))
