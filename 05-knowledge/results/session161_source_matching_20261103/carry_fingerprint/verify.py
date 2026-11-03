from mechanisms import *
from profile import response,triangle
from inherited_recurrences import matrix, words_matrix, krylov
import sympy as s, json

def need(a,msg):
 if not a: raise AssertionError(msg)

# Local copy of the incoming gap automaton, transformed to an exact linear
# certificate in the common power parameter.
# P=2^j+6,R=2^j+8, j>=4; 3 = 11, 0 = 00, 2 = 01, 1 = 10.
dec=lambda x:(x&1,x>>1)
M=words_matrix([dec(int(c)) for c in '330213330'])
A=matrix((0,0));B=matrix((1,1));systems=[]; dimensions=[]
for reverse in (False,True):
 AA=matrix((0,0),reverse);BB=matrix((1,1),reverse)
 alpha=s.zeros(1,27);alpha[0,0]=1
 alpha=alpha*words_matrix([(1,1)],reverse)*AA
 beta=words_matrix([(0,0)],reverse)*s.Matrix([int(all(q!=1 for q in state) and (state[1]!=0 or reverse)) for state in STATES])
 mid=words_matrix([dec(int(c)) for c in '330213330'],reverse)
 left,LA=krylov(alpha,AA,True);right,RB=krylov(beta,BB,False)
 moment=s.Matrix([[(u.T*mid*v)[0] for v in right] for u in left])
 op=s.kronecker_product(LA,RB);initial=s.zeros(1,op.rows);initial[0]=1
 systems.append((op,initial,s.Matrix(list(moment))))
 dimensions.append([len(left),len(right)])
op=s.diag(*(x[0] for x in systems));initial=s.Matrix.hstack(*(x[1] for x in systems));final=s.Matrix.vstack(*(x[2] for x in systems))
x=s.Symbol('x'); poly=s.Poly(x*(x-1)*(x-2)*(x-4),x)
res=s.zeros(op.rows,1);power=final
for c in reversed(poly.all_coeffs()): res+=c*power;power=op*power
row=initial
for i in range(op.rows):need((row*res)[0]==0,'recurrence exact Krylov certificate');row=row*op
vals=[]
for j in range(4,13):
 P=2**j+6;w,u=chart(P,P+2)
 expected_template='3'+'0'*(j-3)+'330213330'+'3'*(j-4)+'0'
 need(''.join(str(a+2*b) for a,b in zip(w,u))==expected_template,'word template')
 H=gap_h(w,u,False)
 formula=s.Rational(26665,64)*4**j+s.Rational(31479,16)*2**j+1721-191*(j==4)
 need(H==formula,'closed count');need(gap_h(w,u)==H+2,'unit edge jump')
 vals.append({'j':j,'lower':H})

# Independent contraction formula checked symbolically for all equal-length pairs P,P+2 through 511.
r=s.Symbol('r');checked=0; seen={}
for P in range(1,512):
 R=P+2;k=P.bit_length()
 if R.bit_length()!=k:continue
 t=0
 while (P>>(t+1))&1:t+=1
 H=r**(t+1)+(r*r-1)*((r-1)**t-(r+1)**t)/2
 want=2*r**(2*k-t-1)*H
 w,u=chart(P,R)
 need(s.expand(response(w,u,r)-want)==0,('profile',P));checked+=1
 key=(k,t);polynomial=s.Poly(want,r);least=min(mon[0] for mon,co in polynomial.terms());coeff=polynomial.coeff_monomial(r**least)
 if t==0:need((least,coeff)==(2*k,2),'zero-carry fingerprint')
 else:
  z=(t+1)//2;need((least,coeff)==(2*k-2*z,2 if t%2 else 4*z),'carry fingerprint')
 if key in seen:need(s.expand(seen[key]-want)==0,'equal k,t response')
 else:seen[key]=want

P=22;w,u=chart(P,P+2);lower=source(w,u,False);upper=source(w,u,True)
example=[hp(lower),hp(upper),hp(triangle(lower)),hp(triangle(upper))]
need(example==[139669,139671,11830709,11829687],'subset DP literal example')
out={'status':'PASS','profile_symbolic_pair_checks':checked,'distinct_fingerprints':len(seen),'exact_recurrence_certificate':{'annihilator':str(poly.as_expr()),'two_krylov_dimensions':dimensions,'dimension':op.rows,'zero_residual_tests':op.rows},'lower_values':vals,'literal_counts':example}
print(json.dumps(out,indent=2))
