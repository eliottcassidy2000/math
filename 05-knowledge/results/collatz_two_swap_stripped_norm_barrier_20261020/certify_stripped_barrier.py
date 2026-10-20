#!/usr/bin/env python3
"""Exact controls for an infinite 2,3-free two-edit norm family.
The analytic limit is proved in the adjacent note, not inferred from this run.
"""
import json,math
import sympy as sp
Z=sp.Symbol('Z');checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(msg)
# F9=F3[t]/(t^2+1), represented by pairs. No external finite-field code.
def add(x,y):return ((x[0]+y[0])%3,(x[1]+y[1])%3)
def neg(x):return ((-x[0])%3,(-x[1])%3)
def mul(x,y):return ((x[0]*y[0]-x[1]*y[1])%3,(x[0]*y[1]+x[1]*y[0])%3)
def power(x,n):
 y=(1,0)
 while n:
  if n&1:y=mul(x,y)
  x=mul(x,x);n//=2
 return y
field=[]
for zeta in [(0,1),(0,2)]:
 zm1=add(zeta,(-1,0))
 need(power(zeta,2)==(2,0) and power(zeta,4)==(1,0),'root order4')
 need(power(zm1,2)==zeta and power(zm1,4)==(2,0) and power(zm1,8)==(1,0),'difference order8')
 for r in range(8):need(power(zeta,r)!=zm1,'root collision')
 field.append(dict(zeta=zeta,zeta_minus_one=zm1,order_zeta=4,order_difference=8))
def carry(bits):
 b=0
 for i,v in enumerate(bits):b=(3 if v else 1)*b+(v<<i)
 return b
def swap(bits,i):
 out=bits[:];j=(i+1)%len(out);need(out[i]!=out[j],'equal bits in path')
 out[i],out[j]=out[j],out[i];return out
rows=[]
for a in range(2,9):
 p=2*3**a;three=3**p;L=three.bit_length();q=L+(1-L)%6;r=q-p;n=p-r
 need(2**(L-1)<three<2**L,'critical ceil')
 need(q>=L and q-L<=5 and q%6==1,'least adjusted q')
 need(math.gcd(p,r)==1 and p<q<2*p and 2*r>p,'shape and legality hypotheses')
 need(n>=5 and n%6==5,'nonconstant term separation')
 follower=pow(r,-1,p);leader=(follower+1)%p
 sf=follower+(follower*r+p-1)//p-1
 sl=leader+(leader*r+p-1)//p-1
 sf%=q;sl%=q
 need((-follower*r)%p==p-1 and (-leader*r)%p==p-1-r,'residue labels')
 source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
 need(source[sf]==source[sl]==1 and source[(sl+1)%q]==0 and sl==(sf+1)%q,'source110')
 first=swap(source,sl);last=swap(first,sf)
 need(first!=source and last!=first and sum(last)==p,'actual swaps')
 D=2**q-three;B=carry(last);g=math.gcd(D,B)
 # Mod3 polynomial arithmetic can be reduced to Z^2+1 by Frobenius.
 norm3=1
 for zeta in [(0,1),(0,2)]:
  Q=add(add((0,0),neg(power(zeta,p-1))),add(power(zeta,n),neg(power(zeta,n-1))))
  need(Q!=(0,0),'Q common root mod3')
  # Product over the two conjugates, then the odd Frobenius multiplicity.
  norm3=mul((norm3,0),Q) if isinstance(norm3,int) else mul(norm3,Q)
 need(norm3[1]==0 and norm3[0] in[1,2],'field norm unit')
 row=dict(a=a,p=p,critical_length=L,q=q,integer_offset=q-L,r=r,n=n,
          follower_label=follower,leader_label=leader,ordered_edges=[sl,sf],
          actual_content=g,norm_mod2=1,norm_mod3=norm3[0])
 if a<=4:
  Q=3-Z**(p-1)+Z**n-Z**(n-1)
  R=int(sp.resultant(Z**p-2,Q,Z));M=abs(R)
  need(math.gcd(M,6)==1,'resultant not2,3free')
  need(R%3==norm3[0] and R%2==1,'modular resultant mismatch')
  row.update(raw_resultant=R,stripped_norm=M,three_to_p=three)
 rows.append(row)
# Do not silently apply the theorem at the two bad small indices.
small=[]
for a in [0,1]:
 p=2*3**a;L=(3**p).bit_length();q=L+(1-L)%6
 need(q>=2*p,'small-index exception not detected');small.append(dict(a=a,p=p,q=q,reason='outside k=1 band'))
print(json.dumps(dict(status='PASS_EXACT_CONTROLS',checks=checks,finitefield_rows=field,
   excluded_small_indices=small,sequence_rows=rows,
   theorem_scope='All a>=2 by the proof; finite concrete controls a2..8; exact resultants a2..4. Not exactcritical-only and not an actual-content lower bound.'),indent=2))
