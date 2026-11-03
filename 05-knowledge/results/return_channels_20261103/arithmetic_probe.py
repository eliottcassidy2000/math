"""Exact loop valuation language and its Liouville/Mobius observation carriers."""
from collections import defaultdict
from math import prod
import json,sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from common import *
N=48
MACROS=(('B',3,1,-1),('Q',4,1,1),('L2',8,2,-1),('L6',9,6,1),('L45',19,45,1),('L13',17,13,1))
SF=(1,2,6,13,26,78)

def mul(a,b):
    c=[0]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:N+1-i]):c[i+j]+=x*y
    return c

def monomial(k,a=1):
    r=[0]*(N+1);r[k]=a;return r

def add(*polys):return [sum(p[i] for p in polys) for i in range(N+1)]

def reciprocal(d):
    need(d[0]==1,'unit formal denominator');a=[1]+[0]*N
    for n in range(1,N+1):a[n]=-sum(d[k]*a[n-k] for k in range(1,n+1))
    return a

def arith(g):
    ex=[]
    for p in (2,3,5,13):
        e=0
        while g%p==0:g//=p;e+=1
        ex.append(e)
    need(g==1,'declared prime support')
    lam=(-1)**sum(ex);mu=lam if max(ex)<=1 else 0
    return tuple(ex),lam,mu

def main():
    dp=[defaultdict(int) for _ in range(N+1)];dp[0][(1,1)]=1
    for n in range(N+1):
        for (g,phase),count in list(dp[n].items()):
            for name,cost,factor,sg in MACROS:
                if n+cost<=N:dp[n+cost][(g*factor,phase*sg)]+=count
    allg={g for d in dp for g,phase in d};sf={g for g in allg if arith(g)[2]}
    need(sf==set(SF),'exact squarefree products in tested language')
    for g in allg:
        (u,v,w,x),lam,mu=arith(g)
        need(v>=2*w and u>=v-2*w,'valuation cone')
        need(mu==int(sp.mobius(g)),'independent arithmetic squarefree function')
    spectra={}
    for s in (1,-1):
        direct_lambda=[];direct_mu=[];direct_sf=[]
        for d in dp:
            direct_lambda.append(sum(c*arith(g)[1]*(phase if s==-1 else 1) for (g,phase),c in d.items()))
            direct_mu.append(sum(c*arith(g)[2]*(phase if s==-1 else 1) for (g,phase),c in d.items()))
            direct_sf.append(sum(c*int(arith(g)[2]!=0)*(phase if s==-1 else 1) for (g,phase),c in d.items()))
        neutral=add(monomial(0),monomial(3,-s),monomial(4,-1));A=reciprocal(neutral)
        lambda_den=add(neutral,monomial(8,s),monomial(9,-1),monomial(19),monomial(17))
        need(reciprocal(lambda_den)==direct_lambda,'Liouville series')
        w=add(monomial(8,s),monomial(9));t=monomial(17)
        squarefree=add(A,mul(add(w,t),mul(A,A)),[2*x for x in mul(mul(w,t),mul(mul(A,A),A))])
        mw=add(monomial(8,-s),monomial(9));mt=monomial(17,-1)
        mobius=add(A,mul(add(mw,mt),mul(A,A)),[2*x for x in mul(mul(mw,mt),mul(mul(A,A),A))])
        need(squarefree==direct_sf and mobius==direct_mu,'squarefree/Mobius rational series')
        spectra[str(s)]={'lambda_coefficients':direct_lambda,'mu_coefficients':direct_mu,'squarefree_coefficients':direct_sf}
    # Three concrete contexts separate all seven reachable Mobius classes.
    representatives=SF+(4,);contexts=(1,2,13)
    signatures={g:tuple(arith(g*c)[2] for c in contexts) for g in representatives}
    need(len(set(signatures.values()))==7,'minimum seven Mobius response classes')
    joint_reps=SF+(4,45)
    joint={g:(arith(g)[1],)+signatures.get(g,tuple(arith(g*c)[2] for c in contexts)) for g in joint_reps}
    need(len(set(joint.values()))==8,'minimum joint arithmetic classes')
    need({(arith(g)[1],arith(g)[2],phase) for d in dp for (g,phase) in d if g in (4,45)}=={(l,0,e) for l in (-1,1) for e in (-1,1)},'dead arithmetic classes and independent phase')
    U=('r5','a','a');L2=('r24','a','a','a');contexts=((),('b',),('q2',),L2)
    columns=sp.Matrix.hstack(*[sp.Matrix(matrix(w))[:,0] for w in contexts])
    need(hermite_normal_form(columns)==sp.diag(1,2,2),'full integral observation lattice from four macro contexts')
    rows=sp.Matrix([vector(U+w) for w in contexts[:3]])
    hankel=rows*columns[:,:3];need(hankel.det()==309024,'actual rational minimum three')
    for i,left in enumerate(contexts[:3]):
        for j,right in enumerate(contexts[:3]):
            need(literal_frontier(actual_target(U+left+right))[0]==int(hankel[i,j]),'literal Hankel witness')
    result={'status':'PASS','checked_vertex_cost':N,'reachable_gain_values':len(allg),'squarefree_gains':SF,
            'minimum_response_states':{'lambda':2,'mu':7,'lambda_mu':8,'lambda_mu_orientation_phase':16},
            'separating_context_multipliers':(1,2,13),'mu_signatures':signatures,
            'joint_signatures':joint,'series':spectra,
            'macro_hamilton_observer':{'integer_lattice_basis':[[1,0,0],[0,2,0],[0,0,2]],
                'hankel_matrix':[[int(x) for x in row] for row in hankel.tolist()],
                'hankel_determinant':309024,'rational_dimension':3},
            'scope':'Arithmetic of the error multiplier in the declared six-macro language; not H or its prime factors.'}
    (HERE/'arithmetic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','gain_values':len(allg),'squarefree_gains':SF,'minimal_states':result['minimum_response_states']}))
if __name__=='__main__':main()
