from fractions import Fraction
from math import gcd,lcm
import sympy as s,json
pi=[0,4,25,168,1229,9592,78498,664579,5761455,50847534,455052511,4118054813]
a=[s.Rational(pi[k]-pi[k-1],9*10**(k-1)) for k in range(1,len(pi))]
rows=[]
for shift in [0,1]:
 size=6 if shift==0 else 5
 M=s.Matrix(size,size,lambda i,j:a[i+j+shift]);A=M[:-1,:-1];v=M[:-1,-1];x=-A.inv()*v;c=list(x)+[s.Integer(1)];den=lcm(*(z.q for z in c));ci=[int(z*den) for z in c];g=gcd(*ci);ci=[z//g for z in ci];q=(s.Matrix(ci).T*M*s.Matrix(ci))[0]
 assert q<0
 rows.append({'shift':shift,'size':size,'matrix':[[str(z) for z in row] for row in M.tolist()],'determinant':str(M.det()),'integer_negative_vector':ci,'quadratic_value':str(q),'last_schur_pivot':str(M[-1,-1]-(v.T*A.inv()*v)[0])})
x=s.symbols('x')
cert=[]
for shift,p in [(0,x*(x-1)**2*(x+1)*(2*x-1)),(0,x-3*x**2+5*x**4-3*x**5),(1,x*(x-1)**2*(2*x-1))]:
 value=s.factor(sum(c*a[e[0]+shift] for e,c in s.Poly(p*p,x).terms()))
 assert value<0
 cert.append({'shift':shift,'polynomial':str(s.expand(p)),'square':str(s.expand(p*p)),'negative_value':str(value)})
print(json.dumps({'pi_counts':pi,'densities':[str(z) for z in a],'negative_gram_witnesses':rows,'small_certificates':cert},indent=2))
