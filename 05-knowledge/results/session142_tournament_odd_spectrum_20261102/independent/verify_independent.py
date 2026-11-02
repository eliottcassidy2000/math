#!/usr/bin/env python3
"""Independent exact audit; stdlib only, no assert-dependent gates.
Polynomial exponents use carry-free Kronecker indices a+32*b+1024*c.
Input coefficients are expanded by hand, not imported from producer scripts.
Tournament enumeration independently uses recursive permutation backtracking.
"""
from fractions import Fraction as Q
from itertools import product
from collections import Counter
import json

def need(c, m):
    if not c:
        raise RuntimeError(m)

def clean(p): return {k:v for k,v in p.items() if v}
def add(p,q):
    r=p.copy()
    for k,v in q.items(): r[k]=r.get(k,0)+v
    return clean(r)
def neg(p): return {k:-v for k,v in p.items()}
def sub(p,q): return add(p,neg(q))
def mul(p,q):
    r={}
    for a,u in p.items():
        for b,v in q.items():
            # Prevent any Kronecker collision via exponent carries.
            need(all(((a>>(5*i))&31)+((b>>(5*i))&31)<32 for i in range(3)), 'exponent carry')
            r[a+b]=r.get(a+b,0)+u*v
    return clean(r)
def mon(a,b,c): return a+32*b+1024*c
def poly(terms): return {mon(*m):v for m,v in terms}
def deriv(p,i):
    r={}
    for m,c in p.items():
        e=(m>>(5*i))&31
        if e: r[m-(1<<(5*i))]=e*c
    return r
def evaluate(p,xyz):
    return sum(c*xyz[0]**(m&31)*xyz[1]**((m>>5)&31)*xyz[2]**((m>>10)&31) for m,c in p.items())
def det3(j):
    return add(sub(mul(j[0][0],sub(mul(j[1][1],j[2][2]),mul(j[1][2],j[2][1]))),
                   mul(j[0][1],sub(mul(j[1][0],j[2][2]),mul(j[1][2],j[2][0])))),
               mul(j[0][2],sub(mul(j[1][0],j[2][1]),mul(j[1][1],j[2][0]))))
def formula(x,y,z):
    u=1+x*y
    return (u**3*z+y*y*u*(4+3*x*y),y+3*x*u*u*z+3*x*y*y*(4+3*x*y),2*x-3*x*x*y-x*x*x*z)
def jacobian():
    f=[poly([((0,0,1),1),((1,1,1),3),((2,2,1),3),((3,3,1),1),((0,2,0),4),((1,3,0),7),((2,4,0),3)]),
       poly([((0,1,0),1),((1,0,1),3),((2,1,1),6),((3,2,1),3),((1,2,0),12),((2,3,0),9)]),
       poly([((1,0,0),2),((2,1,0),-3),((3,0,1),-1)])]
    for p in product((Q(-2),Q(-1,2),Q(0),Q(1),Q(3,2)),repeat=3):
        need(tuple(evaluate(v,p) for v in f)==formula(*p),'manual expansion mismatch')
    j=[[deriv(v,i) for i in range(3)] for v in f]
    need(det3(j)=={0:-2},'exact Jacobian determinant')
    pts=[(Q(0),Q(0),Q(-1,4)),(Q(1),Q(-3,2),Q(13,2)),(Q(-1),Q(3,2),Q(13,2))]
    ims=[formula(*p) for p in pts]
    need(len(set(pts))==3 and len(set(ims))==1 and ims[0]==(Q(-1,4),Q(0),Q(0)),'three-point collision')
    # Hostile coefficient change must be caught by the exact determinant gate.
    g=[v.copy() for v in f]; g[0][mon(1,3,0)]+=1
    bad=det3([[deriv(v,i) for i in range(3)] for v in g])
    need(bad!={0:-2},'mutated coefficient escaped')
    need(det3([[{0:1},{},{}],[{},{0:1},{}],[{},{},{0:1}]])=={0:1},'identity determinant control')
    return {'determinant':-2,'term_counts':[len(v) for v in f], 'manual_expansion_grid_points':125,
            'points':[[str(x) for x in p] for p in pts],'images':[[str(x) for x in p] for p in ims],
            'hostile_changed_F1_xy3_coefficient':{'determinant_terms':len(bad),'rejected':True}}

def edgepairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def matrix(n,m,varied=None,t=None):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(edgepairs(n)):
        v=t if k==varied else ((m>>k)&1)
        a[i][j],a[j][i]=v,1-v
    return a
def path_count(a):
    n=len(a)
    def step(last,unused):
        if not unused: return 1
        return sum(a[last][v]*step(v,unused-{v}) for v in unused if a[last][v])
    return sum(step(v,set(range(n))-{v}) for v in range(n))
def connected(a):
    n=len(a)
    for s in range(n):
        seen={s}; todo=[s]
        while todo:
            for v in range(n):
                if a[todo[-1]][v] and v not in seen: seen.add(v); todo.insert(0,v)
            todo.pop()
        if len(seen)!=n: return False
    return True

def tournaments():
    spectra={}; cases=0; parity_coefficients={}
    for n in range(1,7):
        counts=Counter(); strong=Counter(); values=[]
        for code in range(1<<(n*(n-1)//2)):
            a=matrix(n,code); h=path_count(a); counts[h]+=1; values.append(h)
            if connected(a): strong[h]+=1
            cases+=1
        need(all(h%2 and h not in (7,21) for h in counts),'finite parity/gap census')
        if n<=5:
            coef=values.copy()
            for bit in range(n*(n-1)//2):
                for mask in range(len(coef)):
                    if mask&(1<<bit): coef[mask]-=coef[mask^(1<<bit)]
            need(coef[0]==1 and all(c%2==0 for c in coef[1:]),'multiaffine coefficient parity')
            parity_coefficients[n]={'coefficient_slots':len(coef),'nonconstant_nonzero':sum(c!=0 for c in coef[1:]),'all_nonconstant_even':True}
        spectra[n]={'labelled_count':sum(counts.values()),'strong_labelled_count':sum(strong.values()),
                    'all_values':sorted(counts),'strong_values':sorted(strong)}
    need(spectra[5]['strong_values']==[9,11,13,15],'strong n5 independent census')
    need(spectra[6]['strong_values']==[15,17,19,23,25,27,29,31,33,37,41,43,45],'strong n6 independent census')
    hostile=[]
    for n,lo,hi,h0 in [(5,4,260,5),(6,20,28,19)]:
        edge=(lo^hi).bit_length()-1
        need(lo^hi==1<<edge,'one-edge mask')
        samples=[]
        for t in (Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)):
            h=path_count(matrix(n,lo,edge,t)); need(h==h0+4*t,'weighted line')
            samples.append([str(t),str(h)])
        need(path_count(matrix(n,lo,edge,Q(1,2))) in (7,21),'midpoint target')
        hostile.append({'n':n,'masks':[lo,hi],'edge':edgepairs(n)[edge],'values':samples,'boolean_midpoint':False})
    unimodular=[]
    for n,lo,hi,h0 in [(5,4,5,5),(6,20,16,19)]:
        edge=(lo^hi).bit_length()-1
        need(lo^hi==1<<edge,'unimodular one-edge masks')
        samples=[]
        for t in (Q(-1),Q(0),Q(1,2),Q(1)):
            p=1-t if (lo>>edge)&1 else t
            h=path_count(matrix(n,lo,edge,p))
            need(h==h0-2*t,'unimodular normalized one-edge response')
            samples.append([str(t),str(h),str((h-1)/2)])
        unimodular.append({'n':n,'masks':[lo,hi],'edge':edgepairs(n)[edge],'values_t_H_Q':samples,'forbidden_target_integer_inverse':-1,'legal_inverse':False})
    return {'unimodular_hostiles':unimodular,'cases':cases,'spectra':spectra,'weighted_hostiles':hostile,'coefficient_parity':parity_coefficients,'scope':'Exact n<=6 Boolean census; interpolation is not tournament realizability.'}

def main():
    print(json.dumps({'status':'PASS','jacobian':jacobian(),'tournaments':tournaments()},indent=2,sort_keys=True))
if __name__=='__main__':main()
