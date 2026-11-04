#!/usr/bin/env python3
import itertools,json,math,hashlib
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'returns'
if not SRC.is_dir(): SRC=HERE.parent/'session189-returns'

def need(x,msg):
 if not x:raise RuntimeError(msg)

def tournament(n,mask):
 rows=[0]*n
 for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if mask>>b&1:rows[i]|=1<<j
  else:rows[j]|=1<<i
 return rows

def strong(rows):
 full=(1<<len(rows))-1
 for v in range(len(rows)):
  seen=1<<v;todo=seen
  while todo:
   b=todo&-todo;todo-=b;new=rows[b.bit_length()-1]&~seen;seen|=new;todo|=new
  if seen!=full:return False
 return True

def paths(rows):
 n=len(rows);N=1<<n;H=[0]*N;H[0]=1;S=[[0]*n for _ in range(N)];E=[[0]*n for _ in range(N)];J={}
 def visit(s,v,mask):
  H[mask]+=1;S[mask][s]+=1;E[mask][v]+=1
  J[mask,s,v]=J.get((mask,s,v),0)+1
  opts=rows[v]&~mask
  while opts:
   b=opts&-opts;opts-=b;visit(s,b.bit_length()-1,mask|b)
 for s in range(n):visit(s,s,1<<s)
 return H,S,E,J

def fields(rows):
 n=len(rows);full=(1<<n)-1;H,S,E,J=paths(rows);out=[]
 for v in (0,1,2):
  w=4
  p2=sum(H[a]*H[full^a] for a in range(full+1))
  end=lambda u:sum(E[a][u]*H[full^a] for a in range(full+1))
  start=lambda u:sum(S[a][u]*H[full^a] for a in range(full+1))
  cross=lambda u,t:sum(E[a][u]*S[full^a][t] for a in range(full+1))
  out.extend([H[full],p2,E[full][v],S[full][v],H[full^(1<<v)],end(v),start(v),E[full][w],S[full][w],H[full^(1<<w)],end(w),start(w),cross(v,v),cross(v,w),cross(w,v),cross(w,w)])
 return out

def make_D(U,r,m):
 n=len(U);w=n;Z=list(range(n+1,n+1+m));rows=list(U)+[0]*(m+1)
 for i in range(n):
  if i==r:rows[i]|=1<<w
  else:rows[w]|=1<<i
 for j,z in enumerate(Z):
  rows[w]|=1<<z
  rows[z]|=(1<<n)-1
  for zz in Z[j+1:]:rows[z]|=1<<zz
 return rows,w

def birth(rows,r):
 n=len(rows);new=list(rows)+[((1<<n)-1)^(1<<r)];new[r]|=1<<n;return new

def generic_checks():
 controls=[]
 for n in (3,4):
  for mask in range(1<<(n*(n-1)//2)):
   U=tournament(n,mask)
   if not strong(U):continue
   h,s,e,j=paths(U);full=(1<<n)-1
   for r in range(n):
    f=sum(e[a][r]*h[full^a] for a in range(full+1))
    joint=[sum(count*h[full^a] for (a,start,end),count in j.items() if start==v and end==r) for v in range(n)]
    for m in (1,2):
     D,w=make_D(U,r,m);hd,sd,ed,_=paths(D);fd=(1<<len(D))-1
     ew=sum(ed[a][w]*hd[fd^a] for a in range(fd+1))
     need(hd[fd]==h[full]+2**m*f,'general D Hamilton')
     need(sd[fd][w]==h[full] and hd[fd^(1<<w)]==h[full] and ew==hd[fd],'retained root response')
     X=D
     for parent in (w,len(D),len(D)+1):X=birth(X,parent)
     hx,sx,_,_=paths(X);fx=(1<<len(X))-1
     need(hx[fx]==5*h[full]+6*2**m*f,'general X Hamilton')
     need(sx[fx][:n]==joint,'general core joint restoration')
     need(strong(D) and strong(X),'general core strongness')
     controls.append({'n':n,'mask':mask,'root':r,'m':m,'H_core':h[full],'E2_root':f,'joint':joint,'H_D':hd[fd],'H_X':hx[fx],'old_starts_X':sx[fx][:n]})
 return controls

def inverse(M):
 n=len(M);A=[[Fraction(x)for x in row]+[Fraction(i==j)for j in range(n)]for i,row in enumerate(M)]
 for j in range(n):
  pivot=next(i for i in range(j,n) if A[i][j]);A[j],A[pivot]=A[pivot],A[j];q=A[j][j];A[j]=[x/q for x in A[j]]
  for i in range(n):
   if i!=j:
    q=A[i][j]
    if q:A[i]=[x-q*y for x,y in zip(A[i],A[j])]
 return [r[n:]for r in A]

def main():
 saved=json.loads((SRC/'response.json').read_text());rec=json.loads((SRC/'receipt.json').read_text());C=rec['symbolic_coefficients_including_zero_base'];samples=saved['samples']
 direct=[]
 for sample in samples[:4]:
  calculated=fields(sample['rows']);need(calculated==sample['fields'],'independent simple-path DFS fields')
  direct.append({'m':sample['m'],'old_order':len(sample['rows']),'all_48_fields':calculated})
 inv=inverse([[b**m for b in range(11)]for m in range(11)])
 recovered=[]
 for f in range(48):
  co=[sum(inv[b][m]*samples[m]['fields'][f] for m in range(11))for b in range(11)]
  need(all(x.denominator==1 for x in co),'nonintegral inverse certificate')
  co=[int(x)for x in co];need(co==C[f],'independent Vandermonde coefficient recovery');recovered.append(co)
 # Independent row identities and image sizes; all ten positive bases retained.
 need(all(all(row[b]==0 for b in range(4,11))for row in C),'unexpected higher mode')
 mat=[row[1:4]for row in C];S=mat[3];need(S==[1,0,0],'unit first direction')
 D=[11*mat[34][b]+mat[21][b]-mat[5][b]-73*S[b] for b in range(3)]
 K=[2*mat[37][b]-mat[5][b]+45*S[b]-16*D[b]for b in range(3)]
 need(D==[0,2,0] and K==[0,0,38],'primitive lattice generators')
 need(all(row[1]%2==0 and row[2]%38==0 for row in mat),'lattice upper inclusion')
 mods=[]
 for N in range(1,13):
  seen=set();kernel=[]
  for u in itertools.product(range(N),repeat=3):
   sig=tuple(sum(a*b for a,b in zip(row,u))%N for row in mat);seen.add(sig)
   if not any(sig):
    future=(u[0],2*u[1],3*u[2]);need(not any(sum(a*b for a,b in zip(row,future))%N for row in mat),'future reveals hidden state');kernel.append(u)
  want=N*(N//math.gcd(N,2))*(N//math.gcd(N,38));need(len(seen)==want,'modular image cardinality');mods.append({'N':N,'image_size':len(seen),'kernel_size':len(kernel)})
 spacing=[]
 for m in range(1,21):
  h=25+36*2**m;f=[190*3**m+288*2**m-45,190*3**m+224*2**m-27,114*3**m+160*2**m-45]
  vals=sorted((h+2**k*ff,k,v)for k in range(1,13)for v,ff in enumerate(f));diff=min(b[0]-a[0]for a,b in zip(vals,vals[1:]));need(diff==128*2**m-36,'finite spacing minimum');spacing.append({'m':m,'clocks':[1,12],'minimum_difference':diff})
 controls=generic_checks()
 audit={'result':'PASS','independent_method':'simple-path DFS; rational inverse of 11 exact exponential evaluations; full small modular image enumeration','complete_field_count':48,'direct_field_controls':direct,'all_recovered_symbolic_coefficients':recovered,'general_core_controls':controls,'modular_image_controls':mods,'spacing_controls':spacing,'scope':'all-order proof reviewed separately; finite tests do not establish reusable two-branch closure','source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in(SRC/'PROOF.md',SRC/'verify.py',SRC/'count.cpp',SRC/'response.json',SRC/'receipt.json')}}
 (HERE/'audit_receipt.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps({'result':'PASS','general_core_cases':len(controls),'independent_full_48_field_cases':len(direct),'recovered_symbolic_fields':48,'modular_image_moduli':12,'spacing_control_orders':20}))
if __name__=='__main__':main()
