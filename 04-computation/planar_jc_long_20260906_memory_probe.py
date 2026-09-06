#!/usr/bin/env python3
"""Exact filtered response certificate: lower-weight row-fifteen memory.

Candidate until independent audit. Imports only pinned THM-4308/4315 row
operators; reconstructs the complete raw joint response, not a frozen row.
"""
from pathlib import Path
import sys
import hashlib
import json
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parent))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R
import jc2_source_normal_student_stein_row9_thm4315 as D
CHECKS=0
def check(ok,label):
 global CHECKS
 if not ok: raise AssertionError(label)
 CHECKS+=1

def pin(module,expected):
 actual=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
 check(actual==expected,'dependency '+Path(module.__file__).name)
 print('DEPENDENCY',Path(module.__file__).name,actual)

x,t=R.x,R.t
eps=s.Symbol('eps')
def ex(v): return s.cancel(s.expand(v))
def co(v,j): return s.expand(v).coeff(x,j)
def lin(v): return s.expand(v).coeff(eps,1)
def solve(eq,vs):
 M,b=s.linear_eq_to_matrix(eq,vs)
 piv=M.rref()[1]; rows=M.T.rref()[1]
 vals=M.extract(rows,piv).inv()*b.extract(rows,[0])
 mp={vs[c]:ex(vals[i]) for i,c in enumerate(piv)}
 raw=[ex(q.subs(mp)) for q in eq]
 return mp,raw,M,piv

def main():
 pin(R,'3703758f87a628583cf0f2f9e8fb1973f8ee65c875ab679acc20a9867c27e7f1')
 pin(D,'c9625d0b974c4c579388ec14bcafaacfac5bf4a525d3d02c649d1bbf510b29fa')
 baseA=[R.A0,R.A1,R.A2,R.A3.subs({R.Phi:0,R.Delta:s.Rational(896,15)})]
 baseC=[R.C0,R.C1,R.C2,R.C3.subs({R.Phi:0,R.Delta:s.Rational(896,15)})]
 # A3,C3 contain K as the already defined expression in Delta.
 da=[s.Integer(0)]*12; dc=da.copy(); theta=[]
 for n in range(12,16):
  ba=baseA+[s.Integer(0)]*(n-len(baseA));bc=baseC+[s.Integer(0)]*(n-len(baseC))
  aa=[ba[i]+eps*da[i] for i in range(n)]
  cc=[bc[i]+eps*dc[i] for i in range(n)]
  bp=lin(R.B_row(n,aa,cc))
  pa,pc=R.particular_row(n,bp)
  vs=s.symbols(f'u{n}_0:{n+1}');theta.extend(vs)
  th=sum(v*x**j for j,v in enumerate(vs))
  da.append(s.expand(pa+th*s.diff(R.A0,x)))
  dc.append(s.expand(pc+th*s.diff(R.C0,x)))
  check(s.degree(da[n],x)<=n+1 and s.degree(dc[n],x)<=n+2,'response degree caps')
 mons=[(13,b) for b in range(3,7)]+[(14,b) for b in range(5,8)]+[(15,7),(15,6)]
 rs=s.symbols('r0:'+str(len(mons)))
 source=sum(r*R.p**(v-2*b)*R.y**b for r,(v,b) in zip(rs,mons))
 eq=[]
 for n in range(13,16):
  ba=baseA+[s.Integer(0)]*(n-len(baseA));bc=baseC+[s.Integer(0)]*(n-len(baseC))
  aa=[ba[i]+eps*da[i] for i in range(n)];cc=[bc[i]+eps*dc[i] for i in range(n)]
  delta=lin(R.predicted_G(n,aa,cc))
  target=s.expand(source).coeff(t,n)
  eq.extend(co(target-delta,j) for j in range(n+1))
 check(len(eq)==45,'all bracket rows 13 through 15')
 for depth,polys in ((2,da),(3,dc)):
  coords,mat=D.depth_matrix(depth,15)
  vec=s.Matrix([co(polys[n],r) for n,r in coords])
  eq.extend(ex((left.T*vec)[0]) for left in mat.T.nullspace())
 check((len(eq),len(theta),len(rs))==(136,58,9),'complete joint universe')
 # Raw projected system, with all tangent coordinates and no deduplication.
 mp,raw,M,piv=solve(eq,tuple(theta))
 E,b=s.linear_eq_to_matrix(raw,rs)
 check(M.rank()==48 and E.rank()==3,'joint ranks')
 expected=s.zeros(3,9)
 expected[0,1]=1;expected[0,8]=-s.Rational(27945,235202)
 expected[1,3]=1;expected[1,8]=s.Rational(39123,470404)
 expected[2,5]=1;expected[2,8]=s.Rational(52578,117601)
 reduced=E.rref()[0][:3,:]
 check(reduced==expected,'complete source relation space')
 check(all(v==0 for v in b),'homogeneous response target')
 print('UNIVERSE bracket45 depth91 tangent58 source9')
 print('RANK tangent48 source3 joint51 terminal10')
 print('MONOMIALS',[(v-2*j,j,2*v-j) for v,j in mons])
 print('RELATION_RREF',reduced.tolist())
 # Every weight<=23 monomial of valuation>=13 visible through row15.
 universe={(a,j) for a in range(16) for j in range(8)
           if 13<=a+2*j<=15 and 2*a+3*j<=23}
 check(universe=={(v-2*j,j) for v,j in mons[:-1]},'complete filtered Hamiltonian universe')
 for w in (20,21,22,23):
  cols=[i for i,(v,j) in enumerate(mons[:-1]) if 2*v-j<=w]
  rank=E[:,cols].rank()
  augmented=E[:,cols].row_join(E[:,-1]).rank()
  print('WEIGHT',w,'rank',rank,'with_high',augmented)
  check((rank,augmented)==({20:(1,2),21:(1,2),22:(3,3),23:(3,3)}[w]),'sharp filtered replacement threshold')
 # Earlier joint rows: the 22-weight kernel carries a nonzero later debt.
 prefix=list(eq[:29])
 for depth,polys in ((2,da),(3,dc)):
  coords,mat=D.depth_matrix(depth,14)
  vec=s.Matrix([co(polys[n],j) for n,j in coords])
  prefix.extend(ex((left.T*vec)[0]) for left in mat.T.nullspace())
 _,prefix_raw,_,_=solve(prefix,tuple(theta[:42]))
 PE,_=s.linear_eq_to_matrix(prefix_raw,rs)
 prefix_reduced=PE.rref()[0][:PE.rank(),:]
 expected_prefix=s.zeros(2,9)
 expected_prefix[0,1]=1;expected_prefix[0,5]=s.Rational(135,508)
 expected_prefix[1,3]=1;expected_prefix[1,5]=-s.Rational(189,1016)
 check(PE.rank()==2 and prefix_reduced==expected_prefix,'complete earlier-row kernel')
 print('PREFIX_RREF',prefix_reduced.tolist())
 balanced=s.zeros(9,1);balanced[1]=1;balanced[3]=-s.Rational(7,10);balanced[5]=-s.Rational(508,135)
 check(PE*balanced==s.zeros(PE.rows,1) and E*balanced!=s.zeros(E.rows,1),'prefix-neutral later-visible packet')
 # Difference L-R, where R=p^3*y^6 and all weights of L are <=22.
 packet=s.zeros(9,1)
 packet[1]=-s.Rational(27945,235202)
 packet[3]=s.Rational(39123,470404)
 packet[5]=s.Rational(52578,117601)
 packet[8]=-1
 check(E*packet==s.zeros(E.rows,1),'lower packet equals high response')
 parameters=dict(zip(rs,packet))
 final={u:ex(mp.get(u,s.Integer(0)).subs(parameters)) for u in theta}
 outA=[ex(v.subs(final).subs(parameters)) for v in da]
 outC=[ex(v.subs(final).subs(parameters)) for v in dc]
 check(not set().union(*(v.free_symbols for v in outA+outC))-{x},'rational response rows')
 for i,q in enumerate(eq):
  check(ex(q.subs(final).subs(parameters))==0,'joint coefficient '+str(i))
 for n in range(12,16):
  print('DELTA_A',n,s.sstr(s.factor(outA[n])))
  print('DELTA_C',n,s.sstr(s.factor(outC[n])))
 # Separate literal nonlinear coefficient check, with a generic low-degree
 # prefix. Rows >=4 cannot affect a perturbation starting at row12 by row15.
 B=sum(baseA[i]*t**i for i in range(4));C=sum(baseC[i]*t**i for i in range(4))
 dA=sum(outA[n]*t**n for n in range(12,16));dC=sum(outC[n]*t**n for n in range(12,16))
 lam=s.Symbol('lambda')
 target=ex(source.subs(parameters))
 fullP=R.P(B+lam*dA,C+lam*dC)-R.P(B,C)-lam*target
 fullJ=R.jacobian(B+lam*dA,C+lam*dC)-R.jacobian(B,C)
 for n in range(16):check(ex(s.expand(fullP).coeff(t,n))==0,'literal nonlinear P row '+str(n))
 for n in range(15):check(ex(s.expand(fullJ).coeff(t,n))==0,'literal nonlinear J row '+str(n))
 check(s.expand(fullP).coeff(lam,2).coeff(t,24)!=0,'first potential nonlinear layer is genuinely nonzero')
 print('NONLINEAR first_quadratic_P_row24_nonzero')
 # Dropping the compensation coordinates destroys the simultaneous response.
 lone=s.zeros(9,1);lone[1]=packet[1];lone[8]=-1
 check(E*lone!=s.zeros(E.rows,1),'new monomial alone cannot replace high packet')
 odd=s.zeros(9,1);odd[6]=1
 check(E*odd==s.zeros(E.rows,1),'y7 odd channel is jointly neutral')
 wrong=packet.copy();wrong[5]+=1
 check(E*wrong!=s.zeros(E.rows,1),'rho coefficient off graph hostile')
 print('HOSTILES single_new_monomial_fails odd_y7_neutral rho_plus1_fails weight21_fails')
 semantic={'monomials':[(v-2*j,j) for v,j in mons],
           'rref':[[str(v) for v in row] for row in reduced.tolist()],
           'deltaA':[str(s.expand(v)) for v in outA[12:]],
           'deltaC':[str(s.expand(v)) for v in outC[12:]]}
 print('SEMANTIC_SHA256',hashlib.sha256(json.dumps(semantic,sort_keys=True,separators=(',',':')).encode()).hexdigest())
 print('SCOPE exact_row15_filtered_replacement_valuation_at_least13_base_rows0to3_fixed_not_full_source_or_termination_or_JC2')
 print('CHECKS',CHECKS,'PASS')
if __name__=='__main__': main()
