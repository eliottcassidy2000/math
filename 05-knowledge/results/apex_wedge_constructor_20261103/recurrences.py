"""All-order diagonal-run recurrence certificates, not interpolation alone."""
from pathlib import Path
from math import isqrt
import json
import sympy as s
from mechanisms import *
HERE=Path(__file__).resolve().parent

def matrix(letter,reverse=False):
    a,b=letter
    if reverse:a,b=b,a
    M=s.zeros(27)
    for i,state in enumerate(STATES):
        for out in step(state,a,b):M[i,INDEX[out]]+=1
    return M
def words_matrix(word,reverse=False):
    M=s.eye(27)
    for x in word:M=M*matrix(x,reverse)
    return M
def krylov(v,M,row):
    base=[];cur=v
    for k in range(28):
        column=cur.T if row else cur
        joined=s.Matrix.hstack(*base,column)
        if joined.rank()==len(base):
            coeff=s.Matrix.hstack(*base).gauss_jordan_solve(column)[0]
            d=len(base);T=s.zeros(d)
            for i in range(d-1):T[i,i+1]=1
            for j in range(d):T[d-1,j]=coeff[j]
            return base,T
        base.append(column);cur=cur*M if row else M*cur
    raise RuntimeError('Cayley-Hamilton bound')

# A template is prefix A^(t+a0) middle B^(t+b0) suffix, t=r-r0.
templates={
 '-2':{'r0':2,'prefix':[],'A':(1,1),'a0':1,'middle':[(0,0),(0,1),(1,0)],'B':(1,1),'b0':2,'suffix':[(0,0)],'roots':[1,2,3,4,6]},
 '0':{'r0':2,'prefix':[(1,1)],'A':(0,0),'a0':2,'middle':[(0,1),(1,0),(1,1),(0,0)],'B':(1,1),'b0':0,'suffix':[(0,0)],'roots':[0,1,2,4]},
 '3':{'r0':3,'prefix':[(1,1)],'A':(0,0),'a0':1,'middle':[(1,1),(1,1),(0,1),(1,0),(0,0),(1,1),(0,0)],'B':(1,1),'b0':0,'suffix':[(0,0)],'roots':[0,1,2,4]}}

def closed(off,r):
    if off==-2:return -7-s.Rational(45,2)*2**r+4*3**r-s.Rational(19,2)*4**r+s.Rational(46,3)*6**r
    if off==0:return s.Rational(269,4)*4**r+s.Rational(315,4)*2**r+5-45*(r==2)
    if off==3:return s.Rational(4303,16)*4**r+s.Rational(7455,8)*2**r+623-263*(r==3)
    raise RuntimeError(off)

def main():
    result={};x=s.Symbol('x')
    for name,t in templates.items():
        off=int(name);systems=[];dims=[]
        for reverse in (False,True):
            A=matrix(t['A'],reverse);B=matrix(t['B'],reverse)
            alpha=s.zeros(1,27);alpha[0,0]=1
            alpha=alpha*words_matrix(t['prefix'],reverse)*A**t['a0']
            # Lower-count graph has apex b -> a: order(a,b) rejects an empty middle.
            terminal=s.Matrix([int(all(q!=1 for q in state) and (state[1]!=0 or reverse)) for state in STATES])
            beta=B**t['b0']*words_matrix(t['suffix'],reverse)*terminal
            middle=words_matrix(t['middle'],reverse)
            left,LA=krylov(alpha,A,True);right,RB=krylov(beta,B,False)
            moment=s.Matrix([[(u.T*middle*v)[0] for v in right] for u in left])
            operator=s.kronecker_product(LA,RB)
            initial=s.zeros(1,operator.rows);initial[0,0]=1
            final=s.Matrix(list(moment))
            systems.append((operator,initial,final));dims.append([len(left),len(right)])
        operator=s.diag(*(t[0] for t in systems));initial=s.Matrix.hstack(*(t[1] for t in systems));final=s.Matrix.vstack(*(t[2] for t in systems))
        poly=s.Poly(s.prod(x-root for root in t['roots']),x)
        residue=s.zeros(operator.rows,1);power=final
        for coefficient in reversed(poly.all_coeffs()):
            residue+=coefficient*power;power=operator*power
        # Exactly dim tests, then Cayley-Hamilton proves all future residuals zero.
        checkrow=initial
        for j in range(operator.rows):
            need((checkrow*residue)[0]==0,('Krylov annihilator',name,j));checkrow=checkrow*operator
        values=[];row=initial
        for k in range(max(8,poly.degree()+2)):
            r=t['r0']+k;w,u=chart((1<<r)+off,(1<<r)+off+2)
            literal_types=list(zip(w,u))
            template=t['prefix']+[t['A']]*(k+t['a0'])+t['middle']+[t['B']]*(k+t['b0'])+t['suffix']
            need(literal_types==template,'literal diagonal word template')
            h=int((row*final)[0]);need(h==gap_h(w,u,False)==closed(off,r),'closed formula and evaluator')
            need(gap_h(w,u,True)==h+2,'unit neighbor throughout family')
            values.append(h);row=row*operator
        result[name]={'r0':t['r0'],'krylov_dimensions_by_apex_order':dims,'certificate_dimension':operator.rows,
                      'annihilator_coefficients':[int(c) for c in poly.all_coeffs()],
                      'factorization':str(s.factor(poly.as_expr())),'initial_lower_counts':values,
                      'template':t}
    # Source-preserving 7p phase: the upper count of the offset-three family.
    phases=[];residue=0
    def quotient_mod(t,mod):
        M=7*mod;z=pow(8,t,M)
        value=(275392*z*z+29820*z+625)%M
        need(value%7==0,'7-lane quotient integrality');return value//7
    for e in range(1,8):
        modulus=7**e;previous=7**(e-1)
        candidates=[residue+j*previous for j in range(7) if quotient_mod(residue+j*previous,modulus)==0]
        need(len(candidates)==1,'unique forbidden phase lift');residue=candidates[0]
        phases.append({'e':e,'t_residue':residue,'t_modulus':modulus,'r_residue':5+3*residue,'r_modulus':3*modulus})
    need(phases[0]['t_residue']==1,'mod49 phase')
    prime=43691
    need(all(prime%d for d in range(2,isqrt(prime)+1)),'exact trial division prime')
    need(int(closed(3,5))+2==7*prime,'explicit 7p source')
    output={'status':'PROVED_FINITE_RECURRENCE_CERTIFICATE','families':result,'seven_adic_phases':phases,
            'prime_control':{'r':5,'lower':305835,'upper':305837,'upper_over_7':prime,'order':16}}
    (HERE/'recurrences.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')
    print('PROVED recurrence certificates', {k:v['certificate_dimension'] for k,v in result.items()})
if __name__=='__main__':main()
