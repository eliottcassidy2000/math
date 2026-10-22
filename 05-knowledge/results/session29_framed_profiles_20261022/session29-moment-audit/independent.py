"""Independent symbolic audit. Does not import producer code."""
import json
import sympy as S

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

z,x=S.symbols('z x')
K=S.Matrix([[172,-24,408],[24,-4,68],[408,-68,706]])
J=S.diag(1,-1,1)
e=S.Matrix([1,0,0])
p=x**3-874*x**2-43344*x-50400
D=1-874*z-43344*z**2-50400*z**3
N=1-702*z+1800*z**2
check(K.charpoly(x).as_expr()==p,'characteristic polynomial')
check((S.eye(3)-z*K).det()==D,'determinant resolvent')
check(S.cancel(((S.eye(3)-z*K).inv())[0,0]-N/D)==0,'scalar resolvent')
check(S.gcd(N,D)==1,'no cancellation')
check(J*K==K.T*J,'self-adjoint relation')
B=S.Matrix.hstack(e,K*e,K**2*e)
check(B.det()==-4406400 and K.det()==50400,'Krylov and transfer determinants')
check((K**2-874*K-45144*S.eye(3))*e==S.Matrix([0,10800,0]),'negative polynomial')

# Derive the sequence independently from the rational numerator/denominator,
# then compare Hankel matrices and matrix-power controls.
h=[]
num=[1,-702,1800]
for r in range(80):
    value=num[r] if r<3 else 0
    for shift, coefficient in ((1,874),(2,43344),(3,50400)):
        if r>=shift:
            value+=coefficient*h[r-shift]
    h.append(value)
check(h[:6]==[1,172,195472,178348096,164357443072,151388576906752],'initial moments')
q=[-45144,-874,1]
negative_square=sum(q[i]*q[j]*h[i+j] for i in range(3) for j in range(3))
check(negative_square==-116640000,'moment negative square')
for s in range(76):
    H=S.Matrix(3,3,lambda i,j:h[s+i+j])
    check(H.det()==-19416360960000*50400**s,'shifted Hankel determinant')
    check(H==B.T*J*K**s*B,'shifted Gram identity')
    check(h[s]>0,'positive moment excludes negative definiteness')

# Spectral weights can be certified without radicals: residues of the
# resolvent in the field Q[x]/(p(x)).
num_x=x*x-702*x+1800
weight=num_x/S.diff(p,x)
# Rational partial fraction weights have power sums given by the trace
# in the quotient algebra. This uses polynomial modular inverses, not roots.
inv_deriv=S.invert(S.diff(p,x),p,x)
residue_poly=S.rem(num_x*inv_deriv,p,x)
power_sums=[S.Integer(3),S.Integer(874),S.Integer(874**2+2*43344)]
for r in range(3,82):
    power_sums.append(874*power_sums[-1]+43344*power_sums[-2]+50400*power_sums[-3])
coeffs=[S.expand(residue_poly).coeff(x,i) for i in range(3)]
for r in range(80):
    spectral_sum=sum(coeffs[i]*power_sums[r+i] for i in range(3))
    check(spectral_sum==h[r],'signed spectral expansion via quotient trace')
L=S.Rational('921.115395614954');U=S.Rational('921.115395614956')
check(p.subs(x,-46)<0<p.subs(x,-45),'left negative root')
check(p.subs(x,-2)>0>p.subs(x,-1),'middle negative root')
check(p.subs(x,L)<0<p.subs(x,U),'dominant positive root')
check(S.diff(p,x).subs(x,L)>0,'positive derivative bound')
check(2*L-702>0 and 6*L-1748>0,'monotonic numerator and derivative')
cl=num_x.subs(x,L)/S.diff(p,x).subs(x,U)
cu=num_x.subs(x,U)/S.diff(p,x).subs(x,L)
check(0<cl<cu,'coefficient enclosure')
print(json.dumps({'verdict':'PASS','independent_moment_count':len(h),'shifted_gram_checks':76,'negative_square':int(negative_square),'residue_polynomial_mod_p':str(residue_poly),'root_brackets':[['-46','-45'],['-2','-1'],[str(L),str(U)]],'c_interval':[str(cl),str(cu)],'c_decimal_for_orientation':[str(S.N(cl,18)),str(S.N(cu,18))],'rho_over_18_decimal_for_orientation':[str(S.N(L/18,18)),str(S.N(U/18,18))]},indent=2))
