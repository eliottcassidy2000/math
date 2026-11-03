"""Exact formal unipotent clocks, and their obstruction as legal birth words."""
import json,math
import sympy as s
from common import HERE,need,save,data,rref

def modular(M,m):return s.Matrix([[int(x.p)*pow(int(x.q),-1,m)%m for x in row] for row in M.tolist()])
def power(A,k,m):
    R=s.eye(A.rows)
    while k:
        if k&1:R=(R*A).applyfunc(lambda x:x%m)
        A=(A*A).applyfunc(lambda x:x%m);k//=2
    return R

def main():
    a=json.loads((HERE/'algebra_exact.json').read_text());d=json.loads((HERE/'multiplication.json').read_text())
    coefficients=[1,0,1,0,0,1,1,1,2,1,0,1,2,1,2,1,0,2,0,1,0,1,2,2,1,0,2,0,0,1,0,0,0,0,2,0,0,2,1,1,2,1,0,1,2,1,0,1,0,1,1,0,0,2,0,2,0,2,2,0,2,1,0,2,2,2,0,1,2]
    v=s.Matrix([coefficients])*s.Matrix(a['trace_kernel_coordinates']);M=s.Matrix(27,27,list(v*s.Matrix(a['basis'])));I=s.eye(27)
    need(M**4==s.zeros(27) and M**3!=s.zeros(27),'length-four literal nilpotent')
    D=int(s.ilcm(*[x.q for x in M]));content=math.gcd(*(int(D*x) for x in M));need(D==2 and content==1,'odd-prime integral and primitive model')
    need(modular(M**3,3)!=s.zeros(27),'ternary cube survives')
    ranks={k:len(rref(s.Matrix(A))[1]) for k,A in a['generators'].items()};need(max(ranks.values())<27,'all nonempty legal positive words are singular')
    tests=[]
    for p in (3,5,7,11,13):
        for exponent in (1,2,3):
            m=p**exponent;U=modular(I+M,m);order=p**(exponent+int(p==3))
            need(power(U,order,m)==I and power(U,order//p,m)!=I,'exact prime-power algebraic clock')
            tests.append({'prime':p,'precision':exponent,'modulus':m,'exact_order':order})
    save('radical_clocks.json',{'status':'EXACT_FORMAL_OPERATOR_CLOCKS_WITH_OPERATIONAL_HOSTILE',
      'radical_basis_coefficients':coefficients,'nilpotent_matrix':data(M),'denominator':D,'integer_content':content,'legal_generator_ranks':ranks,
      'clocks':tests,'scope':'I+M is an algebra element, not a nonempty legal C3p/C5p/N/J word; its full rank forbids that interpretation'})
    print('PASS: formal clock p^e (p>=5), 3^(e+1); all positive generators singular',ranks,flush=True)
if __name__=='__main__':main()
