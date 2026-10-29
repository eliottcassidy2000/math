#!/usr/bin/env python3
import itertools,json,math
from collections import Counter
checks=0
def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)
def tournament(n,mask):
 S=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  S[i][j]=1 if mask>>k&1 else -1;S[j][i]=-S[i][j]
 return S
def rankone_cut(S,M):
 O=[i for i in range(len(S)) if i not in M];r=O[0];c=M[0]
 return all(S[i][j]*S[r][c]==S[i][c]*S[r][j] for i in O for j in M)
def module(S,M):return all(len({S[z][x] for x in M})==1 for z in range(len(S)) if z not in M)
def normalize(S,r):
 d=[S[r][i] if i!=r else 1 for i in range(len(S))]
 return [[d[i]*d[j]*S[i][j] for j in range(len(S))] for i in range(len(S))]
def sign_type(S,p):
 n=len(S);answer=[]
 for e in (1,-1):
  d=[1]+[e*S[0][i]*S[p[0]][p[i]] for i in range(1,n)]
  if all(S[p[i]][p[j]]==e*d[i]*d[j]*S[i][j] for i in range(n) for j in range(i+1,n)):answer.append(e)
 return answer
def groups(S):
 P=[];N=[]
 for p in itertools.permutations(range(len(S))):
  typ=sign_type(S,p)
  need(len(typ)<=1,'sign unique for n>=3')
  if typ:(P if typ[0]==1 else N).append(p)
 return P,N
def twopart(n):return n&-n
def paley_root(q):
 sq={i*i%q for i in range(1,q)};S=[[0]*(q+1) for _ in range(q+1)]
 for j in range(1,q+1):S[0][j]=1;S[j][0]=-1
 for i in range(1,q+1):
  for j in range(1,q+1):
   if i!=j:S[i][j]=1 if (j-i)%q in sq else -1
 return S
def projective(q,special):
 out=set()
 for a,b,c,d in itertools.product(range(q),repeat=4):
  det=(a*d-b*c)%q
  if det==0 or (special and det!=1):continue
  p=[]
  for index in range(q+1):
   if index==0:
    p.append(0 if c==0 else 1+a*pow(c,-1,q)%q)
   else:
    x=index-1;den=(c*x+d)%q
    p.append(0 if den==0 else 1+(a*x+b)*pow(den,-1,q)%q)
  out.add(tuple(p))
 return out
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def simple_certificate(G):
 rem=set(G);sizes=[]
 while rem:
  g=next(iter(rem));cl={compose(compose(h,g),inv(h)) for h in G}
  need(cl<=set(G),'conjugacy stays in group');rem-=cl;sizes.append(len(cl))
 sizes.sort();need(sizes==[1,21,24,24,42,56],'PSL2(7) classes')
 other=sizes[1:];possible=[]
 for mask in range(1<<len(other)):
  size=1+sum(x for i,x in enumerate(other) if mask>>i&1)
  if len(G)%size==0:possible.append(size)
 need(sorted(set(possible))==[1,len(G)],'no proper normal subgroup order')
 return {'conjugacy_sizes':sizes,'possible_normal_orders':sorted(set(possible))}

def main():
 cut_count=0;cases=0;hist={}
 for n in range(2,6):
  h=Counter()
  for mask in range(1<<(n*(n-1)//2)):
   S=tournament(n,mask);cases+=1
   for bits in range(1,(1<<n)-1):
    M=[i for i in range(n) if bits>>i&1];r=next(i for i in range(n) if i not in M)
    need(rankone_cut(S,M)==module(normalize(S,r),M),'root normalization bijoin');cut_count+=1
    if module(S,M):
     for dmask in (0,1,(1<<n)-1,mask% (1<<n)):
      d=[-1 if dmask>>i&1 else 1 for i in range(n)];T=[[d[i]*d[j]*S[i][j] for j in range(n)] for i in range(n)]
      need(module(T,M)==(len({d[i] for i in M})==1),'ordinary module survival')
   if n>=3:
    P,N=groups(S);h[str((len(P),len(P)+len(N)))]+=1
    need(n%twopart(len(P))==0,'positive2part divides n')
    need((2*n)%twopart(len(P)+len(N))==0,'extended2part divides2n')
    R=normalize(S,0);aut=sum(all(R[p[i]][p[j]]==R[i][j] for i in range(1,n) for j in range(1,n)) for q in itertools.permutations(range(1,n)) for p in [(0,)+q])
    need(sum(p[0]==0 for p in P)==aut,'root stabilizer exact')
  if n>=3:hist[n]=dict(sorted(h.items()))
 examples={}
 for q in (3,7):
  S=paley_root(q);P,N=groups(S);PSL=projective(q,True);PGL=projective(q,False)
  need(set(P)==PSL,'PSL exact signed automorphisms');need(set(P+N)==PGL,'PGL exact extended')
  examples[q]={'positive_order':len(P),'extended_order':len(P)+len(N),'root_stabilizer':sum(p[0]==0 for p in P),'PSL_matches':True,'PGL_matches':True}
  if q==7:examples[q]['simplicity_certificate']=simple_certificate(P)
 atoms={}
 R5=[[0 if i==j else (1 if (j-i)%5 in (1,2) else -1) for j in range(5)] for i in range(5)]
 for label,S in [('R5',R5),('root_Paley7',paley_root(7))]:
  n=len(S);mods=[];cuts=[]
  for bits in range(1,2**n-1):
   M=[i for i in range(n) if bits>>i&1]
   if 2<=len(M)<=n-1 and module(S,M):mods.append(M)
   if 2<=len(M)<=n-2 and rankone_cut(S,M):cuts.append(M)
  atoms[label]={'ordinary_nontrivial_modules':mods,'nontrivial_signed_cut_sides':cuts}
 need(not atoms['R5']['ordinary_nontrivial_modules'] and len(atoms['R5']['nontrivial_signed_cut_sides'])==10,'ordinary prime signed decomposable')
 need(len(atoms['root_Paley7']['ordinary_nontrivial_modules'])==1 and not atoms['root_Paley7']['nontrivial_signed_cut_sides'],'signed prime ordinary decomposable')
 print(json.dumps({'status':'FINITE-EXACT','checks':checks,'labelled_cases_n2_through5':cases,'proper_cuts':cut_count,'group_histogram':hist,'atomicity_hostiles':atoms,'projective_examples':examples,'scope':'Elementary proofs in companion report; classical bijoin and switching-group theory credited.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
