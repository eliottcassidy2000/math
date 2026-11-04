"""Independent universal eight-term Schur contraction, with free exterior symbols."""
from pathlib import Path
from itertools import permutations
import hashlib,json
NAMES=['tau','alpha','beta','ell_i','ell_j','rho_i','rho_j','E_i','E_j','S_i','S_j','G_ij','G_ji']
def need(c,m):
 if not c:raise RuntimeError(m)
def add(*Ps):
 D={}
 for P in Ps:
  for k,v in P.items():D[k]=D.get(k,0)+v
 return {k:v for k,v in D.items() if v}
def sc(c,P):return {k:c*v for k,v in P.items() if c*v}
def mul(A,B):
 D={}
 for (m,t),v in A.items():
  for (s,u),w in B.items():
   if not m&s:
    key=(m|s,tuple(sorted(t+u)));D[key]=D.get(key,0)+v*w
 return {k:v for k,v in D.items() if v}
def var(n):return {(0,(NAMES.index(n),)):1}
one={(0,()):1}
def coeff(P):return {t:v for (m,t),v in P.items() if m==7}
def matrix(T):
 M=[[{} for _ in T] for _ in T]
 for k in range(1,4):
  for p in permutations(range(3),k):
   z={(sum(1<<v for v in p),()):1}
   for u,v in zip(p,p[1:]):z=mul(z,add(sc(int(bool(T[u]>>v&1)),one),var('tau')))
   M[p[0]][p[-1]]=add(M[p[0]][p[-1]],z)
 return M
alpha,beta,li,lj,ri,rj,Ei,Ej,Si,Sj,Gij,Gji=[var(x) for x in NAMES[1:]]
def quantities(M):return [add(*(M[u][v] for u in range(3))) for v in range(3)],[add(*row) for row in M],add(*(x for row in M for x in row))
def quartic(M,p,q):
 C,R,H=quantities(M);ep=mul(alpha,C[p]);eq=mul(alpha,C[q]);sp=mul(beta,R[p]);sq=mul(beta,R[q]);ei=add(Ei,mul(mul(alpha,H),ri));ej=add(Ej,mul(mul(alpha,H),rj));si=add(Si,mul(mul(beta,li),H));sj=add(Sj,mul(mul(beta,lj),H));giq=mul(li,C[q]);gij=add(Gij,mul(mul(li,H),rj));gpj=mul(R[p],rj);gjp=mul(lj,C[p]);gji=add(Gji,mul(mul(lj,H),ri));gqi=mul(R[q],ri)
 return add(mul(ep,add(mul(giq,sj),mul(gij,sq))),mul(ei,add(mul(M[p][q],sj),mul(gpj,sq))),mul(eq,add(mul(gjp,si),mul(gji,sp))),mul(ej,add(mul(M[q][p],si),mul(gqi,sp))))
def main():
 records=[]
 for T,pairs in [([6,4,0],[(0,1),(1,2),(0,2)]),([4,1,2],[(0,2)])]:
  M=matrix(T);Mt=[list(v) for v in zip(*M)];C,R,H=quantities(M);KS=add(mul(alpha,add(mul(li,Sj),mul(lj,Si))),sc(-1,mul(beta,add(mul(Ei,rj),mul(Ej,ri)))));KD=mul(mul(alpha,beta),add(Gij,sc(-1,Gji)));B0=add(mul(Ei,Sj),sc(-1,mul(Ej,Si)));B1=add(mul(alpha,add(mul(ri,Sj),sc(-1,mul(rj,Si)))),mul(beta,add(mul(Ei,lj),sc(-1,mul(Ej,li)))))
  for p,q in pairs:
   S=add(mul(C[p],C[q]),sc(-1,mul(R[p],R[q])));D=add(mul(C[p],R[q]),sc(-1,mul(C[q],R[p])));E=add(M[p][q],sc(-1,M[q][p]));rhs=add(mul(S,KS),mul(D,KD),mul(E,B0),mul(mul(H,E),B1));lhs=add(quartic(M,p,q),sc(-1,quartic(Mt,p,q)));need(coeff(lhs)==coeff(rhs),'universal full eight-term Schur identity with free exterior symbols');records.append({'rows':T,'pair':[p,q],'full_support_polynomial':[[list(k),v] for k,v in sorted(coeff(lhs).items())]})
 maps=[dict((tuple(k),v) for k,v in r['full_support_polynomial']) for r in records];combined={}
 for c,P in zip([1,1,-1],maps[:3]):
  for k,v in P.items():combined[k]=combined.get(k,0)+c*v
 need({k:v for k,v in combined.items() if v}==maps[3],'universal circulation, not only literal finite controls')
 print(json.dumps({'status':'PASS','exterior_variables':NAMES,'records':records,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()
