"""Exact cyclic equal-profile gap; standard-library symbolic and literal checks."""
import itertools,json,math,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'response'))
from verify import tournament,endpoints,extend,covers,opposite,need

def hp(rows):
 # Independent predecessor recurrence; C++ uses forward endpoint extension.
 n=len(rows)
 if not n:return 1
 incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
 dp=[[0]*n for _ in range(1<<n)]
 for mask in range(1,1<<n):
  rem=mask
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;prev=mask^bit
   if not prev:dp[mask][v]=1;continue
   src=prev&incoming[v];z=0
   while src:
    ub=src&-src;src-=ub;z+=dp[prev][ub.bit_length()-1]
   dp[mask][v]=z
 return sum(dp[-1])


def add(*ps):
 d={}
 for p in ps:
  for a,c in p.items():d[a]=d.get(a,0)+c
 return {a:c for a,c in d.items() if c}
def scale(p,c):return {a:b*c for a,b in p.items() if b*c}
def mul(p,q):
 d={}
 for a,c in p.items():
  for b,e in q.items():
   z=tuple(x+y for x,y in zip(a,b));d[z]=d.get(z,0)+c*e
 return {a:c for a,c in d.items() if c}
def prod(*ps):
 z={(0,0,0):1}
 for p in ps:z=mul(z,p)
 return z
ONE={(0,0,0):1};X={(1,0,0):1};Y={(0,1,0):1};Z={(0,0,1):1};P=prod(X,Y,Z);R=add(ONE,scale(P,-1))
G=[[X,prod(X,Y),P],[P,Y,prod(Y,Z)],[prod(X,Z),P,Z]]
def Fnum(G):
 def left(I):return add(R,*[G[u][v] for u in range(3) for v in I])
 def right(I):return add(R,*[G[u][v] for u in I for v in range(3)])
 def middle(I,J):return add(*[G[u][v] for u in I for v in J])
 return add(prod(left([0]),add(R,middle([1,2],[1])),right([0,2])),prod(left([1]),middle([0,2],[0]),right([1,2])))
def permute(p,pi):
 d={}
 for a,c in p.items():
  b=[0,0,0]
  for i,j in enumerate(pi):b[j]=a[i]
  d[tuple(b)]=c
 return d
DELTA=add(Fnum(G),scale(Fnum([list(r) for r in zip(*G)]),-1))
CNUM=mul(prod(X,Y),add(prod(add(ONE,Z,prod(Y,Z)),add(ONE,Z,prod(X,Z))),scale(prod(add(ONE,X,prod(X,Z)),add(ONE,Y,prod(Y,Z))),-1)))
N=add(prod(X,X,Y,Y),prod(X,X,Z,Z),prod(Y,Y,Z,Z),prod(X,X,Y),prod(X,X,Z),prod(X,Y,Y),prod(X,Z,Z),prod(Y,Y,Z),prod(Y,Z,Z),scale(mul(P,add(X,Y,Z)),-1),scale(P,-6))

def functional(num,power,q):
 n=len(q)-1;result=0
 for k in range(n+1):
  fac=math.comb(k+power-1,power-1)
  for e,c in num.items():
   term=fac*c
   for t in e:term*=q[k+t] if k+t<len(q) else 0
   result+=term
 return result

def gap(q):
 a=q+[0,0];return -sum((k+1)*(a[k+2]+2*a[k+1])*(a[k+1]**2-a[k]*a[k+2]) for k in range(len(q)-1))
def square_certificate(q):
 a=q+[0,0];m=len(q)-1
 first=sum(j*a[j]*(a[j-1]-a[j+1])**2 for j in range(1,m+1))
 total=a[1]**3
 for j in range(1,m+2):
  total+=j*(a[j-1]-a[j+1])**2*(a[j-1]+3*a[j]+2*a[j+1])
 for j in range(1,m+1):
  total+=(a[j]-a[j+1])**2*((j+2)*a[j]+(2*j+1)*a[j+1])
 return total

def symbolic_square_certificate(m):
 q=[{}]+[{tuple(int(i==j) for i in range(m)):1} for j in range(m)]+[{},{}]
 left={};right=mul(mul(q[1],q[1]),q[1])
 for j in range(1,m+1):
  a,b,c=q[j-1:j+2]
  left=add(left,scale(mul(add(c,scale(b,2)),add(mul(b,b),scale(mul(a,c),-1))),3*j))
  d=add(b,scale(c,-1))
  right=add(right,mul(mul(d,d),add(scale(b,j+2),scale(c,2*j+1))))
 for j in range(1,m+2):
  a,b,c=q[j-1:j+2];d=add(a,scale(c,-1))
  right=add(right,scale(mul(mul(d,d),add(a,scale(b,3),scale(c,2))),j))
 need(left==right,'coefficientwise positive-square identity, including upper boundary')
 return {'cutoff':m,'cubic_monomials':len(left)}

def lc(q):
 a=q+[0];return all(a[k]**2>=a[k-1]*a[k+1] for k in range(1,len(q)))
def profile(a):
 n=len(a);full=(1<<n)-1;H=[sum(r) for r in endpoints(a)];H[0]=1;D=[[0]*(n+1) for _ in range(1<<n)];D[0][0]=1
 for mask in range(1,full+1):
  sub=mask
  while sub:
   rest=mask^sub
   for k in range(1,n+1):D[mask][k]+=H[sub]*D[rest][k-1]
   sub=(sub-1)&mask
 return D[full]
def tt(n):return [sum(1<<v for v in range(u+1,n)) for u in range(n)]
def blowup(a,reverse=False):
 m=len(a);rows=[0]*(3*m)
 for i in range(3):
  for u in range(m):
   rows[i*m+u]|=a[u]<<(i*m)
   j=(i+(2 if reverse else 1))%3;rows[i*m+u]|=((1<<m)-1)<<(j*m)
 return rows

def main():
 permutations=list(itertools.permutations(range(3)));sym=add(*[permute(DELTA,pi) for pi in permutations]);need(sym==scale(mul(N,R),2),'six-permutation literal-query identity')
 need(add(*[permute(add(DELTA,mul(CNUM,R)),pi) for pi in permutations])=={},'query equals negative block skew after symmetry')
 certificates=[symbolic_square_certificate(m) for m in range(1,21)]
 arbitrary=0
 for m in range(1,6):
  for tail in itertools.product(range(4),repeat=m-1):
   for first in range(1,4):
    q=[0,first]+list(tail);need(square_certificate(q)==-3*gap(q)>0,'unconditional nonnegative-profile strict separator');arbitrary+=1
 profiles=[];count=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code);q=profile(a);need(lc(q),'finite q-log-concavity control');need(functional(DELTA,3,q)==gap(q)==-functional(CNUM,2,q),'exact profile functional identity');need(gap(q)<0,'finite positive-gap control');profiles.append([n,code,q]);count+=1
 for q in ([0,1,20,1],[0,100,1,100],[0,1,1,40,2],[0,3,6,6]):need(functional(DELTA,3,list(q))==gap(list(q)),'algebra does not assume log-concavity')
 transitive=[];q=[0,1]
 for m in range(1,51):
  if m>1:q=[0]+[k*((q[k] if k<len(q) else 0)+q[k-1]) for k in range(1,m+1)]
  need(lc(q),'transitive profile recurrence');need(gap(q)<0,'transitive strict separator');transitive.append({'module_size':m,'q':q,'gap':gap(q)})
 literal=[]
 for name,a in [('TT1',tt(1)),('TT2',tt(2)),('TT3',tt(3)),('TT4',tt(4)),('C3',[2,4,1]),('strong4',[8,1,3,6])]:
  q=profile(a);T=blowup(a);U=blowup(a,True);m=len(a);A=(1<<m)-1;B=A<<m;f=hp(extend(T,A,B));g=hp(extend(U,A,B));need(f-g==gap(q),'literal uniform separator')
  C=covers(T);skew=sum(C[u][v]-C[v][u] for u in range(m) for v in range(m,2*m));need(skew==-gap(q),'independent block cover skew')
  literal.append({'name':name,'module_rows':a,'q':q,'T':T,'U':U,'a_mask':A,'c_mask':B,'F_T':f,'F_U':g,'gap':f-g,'block_cover_skew':skew})
 full=json.loads((ROOT/'full_response.json').read_text());need(full['status']=='PASS' and full['basis_coordinates']==2116,'complete nine-vertex response')
 need(full['rows'][0]==blowup([2,4,1]) and full['rows'][1]==blowup([2,4,1],True),'nine-vertex physical names')
 for rows in full['rows']:
  need(all(any(((rows[w]>>u)^(rows[w]>>v))&1 for w in range(9) if w not in (u,v)) for u in range(9) for v in range(u+1,9)),'no homogeneous pair')
  for u in range(9):
   reached=1<<u
   while True:
    nxt=reached
    for v in range(9):
     if reached>>v&1:nxt|=rows[v]
    if nxt==reached:break
    reached=nxt
   need(reached==511,'strong connectivity')
 half=[a for a in range(512) if a.bit_count()<=2];need(len(full['coordinate_data'])==len(half)**2,'full coordinate count');lookup={(a,c):(f,g,o,cf,cg,co) for a,c,f,g,o,cf,cg,co in full['coordinate_data']}
 hist={};literal_extra=[]
 for a,c in itertools.product(half,repeat=2):
  row=lookup[a,c]
  for h in range(3):
   coef=0;aa=a
   while True:
    cc=c
    while True:
     coef+=(-1)**((a^aa).bit_count()+(c^cc).bit_count())*lookup[aa,cc][h]
     if cc==0:break
     cc=(cc-1)&c
    if aa==0:break
    aa=(aa-1)&a
   need(coef==row[3+h],'independent full interpolation')
  degree=a.bit_count()+c.bit_count()
  if degree==4:need(row[3]==row[5],'global quartic cancellation')
  if row[3]!=row[4]:hist[str(degree)]=hist.get(str(degree),0)+1
 for a,c in [(7,56),(341,170),(511,511),(17,258)]:
  vals=[hp(extend(rows,a,c)) for rows in full['rows']]
  for j in range(3):need(vals[j]==sum(row[3+j] for (x,y),row in lookup.items() if x&a==x and y&c==y),'independent dense literal response')
  literal_extra.append([a,c]+vals)
 need(hist.get('4',0)>0,'covered partial reversal really quartic')
 return {'status':'PASS','scope':'All-order exact equal-profile gap formula and unconditional positive-square noncancellation for every three nonempty equal-profile children. No arbitrary unequal-profile or arbitrary strong ambient transfer is claimed.','positive_square_symbolic_checks':certificates,'arbitrary_nonnegative_profiles_checked':arbitrary,'symbolic_delta_numerator':[[*a,c] for a,c in sorted(DELTA.items())],'symbolic_skew_numerator':[[*a,c] for a,c in sorted(CNUM.items())],'symbolic_symmetric_numerator':[[*a,c] for a,c in sorted(N.items())],'finite_tournament_profiles_checked':count,'finite_profile_records':profiles,'transitive_recurrence_controls':transitive,'literal_blowup_controls':literal,'pair_free_nine':{'complete_coordinates_per_response':2116,'nonzero_partial_difference_by_degree':hist,'independent_dense_literal_queries':literal_extra,'full_response_sha256':hashlib.sha256((ROOT/'full_response.json').read_bytes()).hexdigest()}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
