"""Independent two-exterior-singleton contraction with free A path symbols."""
from pathlib import Path
import hashlib,json
names=['H','Cp','Cq','Rp','Rq','Gpq','Gqp']
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
    k=(m|s,tuple(sorted(t+u)));D[k]=D.get(k,0)+v*w
 return {k:v for k,v in D.items() if v}
def prod(*Ps):
 R={(0,()):1}
 for P in Ps:R=mul(R,P)
 return R
def coeff(P):return {t:v for (m,t),v in P.items() if m==3}
def need(c,m):
 if not c:raise RuntimeError(m)
def expr(H,Cp,Cq,Rp,Rq,Gpq,Gqp):
 one={(0,()):1};u={(1,()):1};v={(2,()):1};uv=mul(u,v)
 Ep=mul(Cp,add(one,v,mul(uv,add(one,H))));Eq=mul(Cq,add(one,v,mul(uv,add(one,H))));Sp=mul(Rp,add(one,u,mul(uv,add(one,H))));Sq=mul(Rq,add(one,u,mul(uv,add(one,H))))
 Eu=mul(u,add(one,H,mul(v,H)));Ev=mul(v,add(one,mul(u,add(one,H))));Su=mul(u,add(one,mul(v,add(one,H))));Sv=mul(v,add(one,H,mul(u,H)))
 Guq=mul(uv,Cq);Guv=uv;Gpq2=add(Gpq,prod(uv,Rp,Cq));Gpv=prod(u,v,Rp);Gvp=mul(v,Cp);Gvu=mul(uv,H);Gqp2=add(Gqp,prod(uv,Rq,Cp));Gqu=mul(u,Rq)
 return add(prod(Ep,Guq,Sv),prod(Ep,Guv,Sq),prod(Eu,Gpq2,Sv),prod(Eu,Gpv,Sq),prod(Eq,Gvp,Su),prod(Eq,Gvu,Sp),prod(Ev,Gqp2,Su),prod(Ev,Gqu,Sp))
def main():
 H,Cp,Cq,Rp,Rq,Gpq,Gqp=[{(0,(i,)):1} for i in range(7)];one={(0,()):1}
 lhs=coeff(expr(H,Cp,Cq,Rp,Rq,Gpq,Gqp));expected=add(mul(Cp,Cq),mul(Cp,Rq),Gqp,prod(add(one,H),add(one,H),Gpq),prod(H,Cq,Rp),mul(Rp,Rq));need(lhs=={t:v for (m,t),v in expected.items()},'complete source eight-term contraction')
 target=coeff(expr(H,Rp,Rq,Cp,Cq,Gqp,Gpq));contrast=dict(lhs)
 for k,v in target.items():contrast[k]=contrast.get(k,0)-v
 contrast={k:v for k,v in contrast.items() if v};D=add(mul(Cp,Rq),sc(-1,mul(Cq,Rp)));E=add(Gpq,sc(-1,Gqp));rhs=add(mul(add(one,sc(-1,H)),D),mul(add(sc(2,H),mul(H,H)),E));need(contrast=={t:v for (m,t),v in rhs.items()},'exact arbitrary-child quartic polynomial')
 print(json.dumps({'status':'PASS','variables':names,'source_full_exterior_coefficient':[[list(k),v] for k,v in sorted(lhs.items())],'contrast_polynomial':[[list(k),v] for k,v in sorted(contrast.items())],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()
