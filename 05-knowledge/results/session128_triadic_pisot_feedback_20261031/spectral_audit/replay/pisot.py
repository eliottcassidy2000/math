#!/usr/bin/env python3
"""Exact, dependency-free Session128 Pisot/GIFS/feedback certificates.
All checks use explicit exceptions, so python -O retains them.
No numerical root approximation is used as proof.
"""
from fractions import Fraction as F
from itertools import permutations
from collections import Counter
import hashlib, json

def need(c,m):
    if not c: raise RuntimeError(m)
def parity(p): return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def mm(A,B): return [[sum(A[i][r]*B[r][j] for r in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def eye(n): return [[int(i==j) for j in range(n)] for i in range(n)]
def det(A): return sum(parity(p)*prod(A[i][p[i]] for i in range(len(A))) for p in permutations(range(len(A))))
def prod(xs):
    v=1
    for x in xs:v*=x
    return v

def polyadd(a,b):
    c=[0]*max(len(a),len(b))
    for i,v in enumerate(a):c[i]+=v
    for i,v in enumerate(b):c[i]+=v
    return c

def polymul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def charpoly(M):
    n=len(M); out=[0]*(n+1)
    for perm in permutations(range(n)):
        term=[parity(perm)]
        for i,j in enumerate(perm): term=polymul(term,[-M[i][j],int(i==j)])
        out=polyadd(out,term)
    return out

class Ring:
    def __init__(self,k):self.k=k;self.zero=(F(0),)*k;self.one=(F(1),)+(F(0),)*(k-1);self.t=(F(0),F(1))+(F(0),)*(k-2)
    def scalar(self,q):return (F(q),)+(F(0),)*(self.k-1)
    def add(self,a,b):return tuple(x+y for x,y in zip(a,b))
    def neg(self,a):return tuple(-x for x in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        v=[F(0)]*(2*self.k-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):v[i+j]+=x*y
        for j in range(len(v)-1,self.k-1,-1):
            for r in range(self.k):v[j-self.k+r]+=v[j]
        return tuple(v[:self.k])
    def power(self,a,n):
        r=self.one
        while n:
            if n&1:r=self.mul(r,a)
            a=self.mul(a,a);n//=2
        return r
    def smul(self,q,a):return tuple(F(q)*x for x in a)

def rmm(R,A,B):
    def dot(a,b):
        v=R.zero
        for x,y in zip(a,b):v=R.add(v,R.mul(x,y))
        return v
    return [[dot(row,col) for col in zip(*B)] for row in A]
def radd(R,A,B):return [[R.add(a,b) for a,b in zip(x,y)] for x,y in zip(A,B)]
def rvec(R,A,v):return tuple(row[0] for row in rmm(R,A,[[x] for x in v]))
def encode(v):return [str(x) for x in v]
def evalpoly(k,x):return x**k-sum(x**j for j in range(k))

def family_certificate(k):
    M=[[int(i==0 or i==j+1) for j in range(k)] for i in range(k)]
    C=[[int(i==j+1 or j==k-1) for j in range(k)] for i in range(k)]
    U=[[1 if i==j else -1 if i<j else 0 for j in range(k)] for i in range(k)]
    need(mm(U,M)==mm(C,U),f'field basis intertwiner k={k}')
    need(det(U)==1,f'unimodular basis k={k}')
    need(charpoly(M)==[-1]*k+[1],f'characteristic polynomial k={k}')
    P=eye(k)
    for _ in range(k):P=mm(P,M)
    need(all(x>0 for row in P for x in row),f'primitive kth power k={k}')
    R=Ring(k);inv=(-F(1),)*(k-1)+(F(1),)
    need(R.mul(R.t,inv)==R.one,f'inverse k={k}')
    need(R.power(inv,k)==R.sub(R.scalar(2),R.t),f'baseline inverse unit k={k}')
    return {'k':k,'M':M,'C':C,'U':U,'charpoly_low_to_high':charpoly(M),'det_U':det(U),'M_power_k_strictly_positive':True}

def graph_certificate(depth=12):
    R=Ring(3);o=R.one;z=R.zero;t=R.t
    basis=[o,(-F(1),F(1),F(0)),(-F(1),-F(1),F(1))]
    # Edge tuple (target i, source j, prefix digit); indices 0-based.
    edges=[(0,0,z),(1,0,o),(0,1,z),(2,1,o),(0,2,z)]
    words=[(0,),(1,),(2,)]
    sigma=[(0,1),(0,2),(0,)]
    paths=[[Counter({z:1}) if i==j else Counter() for j in range(3)] for i in range(3)]
    by_depth=[]
    for n in range(depth+1):
        direct=[[Counter() for j in range(3)] for i in range(3)]
        for j,w in enumerate(words):
            prefix=z
            for i in w:
                direct[i][j][prefix]+=1
                prefix=R.add(prefix,basis[i])
        need(paths==direct,f'graph/ordered-word Fourier coefficient equality depth={n}')
        sets=[set().union(*(set(paths[i][j]) for j in range(3))) for i in range(3)]
        by_depth.append({'depth':n,'word_lengths':[len(w) for w in words],'component_support_counts':[len(s) for s in sets]})
        if n==depth:break
        nxt=[[Counter() for _ in range(3)] for _ in range(3)]
        for i,j,d in edges:
            for source in range(3):
                for v,c in paths[j][source].items():nxt[i][source][R.add(R.mul(t,v),d)]+=c
        paths=nxt;words=[tuple(x for a in w for x in sigma[a]) for w in words]
    serialized=[[[[encode(v),c] for v,c in sorted(paths[i][j].items())] for j in range(3)] for i in range(3)]
    digest=hashlib.sha256(json.dumps(serialized,separators=(',',':')).encode()).hexdigest()
    return R,sets,{'depths':by_depth,'all_i_j_fourier_coefficient_dictionaries_equal':True,'final_dictionary_sha256':digest,'universe':'All occurrences in sigma^n(j), j=1,2,3, 0<=n<=12; multiplicities retained per i,j.'}

def feedback_certificate(R,sets):
    # Source triples in sets are coefficients of the ABSTRACT generator;
    # their star evaluation uses alpha. Below, scalar-ring triples instead
    # evaluate at the REAL beta. coord() is the explicit bridge, not an
    # identification alpha=beta. The two fields have the same defining f.
    o,z,t=R.one,R.zero,R.t;p=(-F(1),-F(1),F(1));s=R.sub(o,t)
    A=[[z,R.neg(p)],[o,s]]
    I=[[o,z],[z,o]];d=(o,z)
    Gram=[[o,R.smul(F(1,2),s)],[R.smul(F(1,2),s),p]]
    transpose=lambda X:list(map(list,zip(*X)))
    pGram=[[R.mul(p,v) for v in row] for row in Gram]
    need(rmm(R,rmm(R,transpose(A),Gram),A)==pGram,"exact complex Gram metric")
    def coord(a):return (R.sub(R.scalar(a[0]),R.smul(a[2],p)),R.add(R.scalar(a[1]),R.smul(a[2],s)))
    tested=[]
    for h in [F(0),F(1),F(-2),F(10),F(7,3)]:
        H=R.scalar(h);S=[[o,H],[z,o]];Si=[[o,R.neg(H)],[z,o]]
        B=rmm(R,rmm(R,S,A),Si)
        G=[[H,R.sub(R.smul(h,s),R.scalar(h*h))],[z,R.neg(H)]]
        need(B==radd(R,A,G),'closed loop shear formula')
        need(rmm(R,B,S)==rmm(R,S,A),'intertwining')
        need(rvec(R,S,d)==d,'digit fixed')
        Gh=rmm(R,rmm(R,transpose(Si),Gram),Si)
        need(rmm(R,rmm(R,transpose(B),Gh),B)==[[R.mul(p,v) for v in row] for row in Gh],'transported exact Gram metric')
        need(R.add(B[0][0],B[1][1])==s,'trace retained')
        need(R.sub(R.mul(B[0][0],B[1][1]),R.mul(B[0][1],B[1][0]))==p,'det retained')
        for a in set().union(*sets):
            need(rvec(R,A,coord(a))==coord(R.mul(t,a)),'complex-coordinate multiplication')
            need(rvec(R,B,rvec(R,S,coord(a)))==rvec(R,S,coord(R.mul(t,a))),'finite orbit conjugacy')
        tested.append(str(h))
    # A+gI has trace s+2g: g=1/10 is an exact nonconjugacy witness.
    g=F(1,10);need(R.add(s,R.scalar(2*g))!=s,'additive scalar gain hostile')
    # Independently recover the twisted digits from ordered prefixes.
    basis=[o,R.sub(t,o),p]
    def digit_graph(words):
        edges=[]
        for j,w in enumerate(words):
            prefix=z
            for i in w:
                edges.append((i,j,prefix));prefix=R.add(prefix,basis[i])
        return sorted(edges)
    standard=digit_graph([(0,1),(0,2),(0,)])
    twisted=digit_graph([(1,0),(0,2),(0,)])
    need([(i,j) for i,j,digit in standard]==[(i,j) for i,j,digit in twisted],'same incidence matrix, different word order')
    need(twisted==sorted([(0,0,R.sub(t,o)),(1,0,z),(0,1,z),(2,1,o),(0,2,z)]),'twisted ordered digits')
    # Twisted digit graph common affine conjugacy forces two incompatible scales.
    inv=p
    c_from_edge21=R.sub(inv,t);c_from_edge32=inv
    need(R.sub(c_from_edge32,c_from_edge21)==t and t!=z,'twisted same-matrix hostile')
    return {'exact_shear_gains_tested':tested,'tested_points':len(set().union(*sets)), 'all_checks_pass':True,'exact_Gram_Lyapunov_identities':True,
            'scalar_additive_gain_hostile':'g=1/10 changes trace by 1/5',
            'twisted_affine_scale_requirements':[encode(c_from_edge21),encode(c_from_edge32)],
            'twisted_scale_difference':encode(t)}

def bounds_and_units(R):
    lower,upper=F(183,100),F(184,100)
    need(evalpoly(3,lower)<0<evalpoly(3,upper),'rational isolating signs')
    need(lower>F(9,5)>F(16,9),'beta lower bounds')
    need(F(5,9)<F(9,16),'conjugate modulus squared below (3/4)^2')
    inv=(-F(1),-F(1),F(1));v=R.power(inv,16)
    need(v==(-F(103),F(56),F(0)),'Session127 inverse16 identity')
    # U^-1 maps the sparse power-basis vector to (a0+a1,a1,0).
    substitution_coordinates=(v[0]+v[1],v[1],v[2])
    need(substitution_coordinates==(-F(47),F(56),F(0)),'Lucas47 substitution-coordinate junction')
    need(F(9,5)**8>100,'inverse16 internal image exceeds100')
    return {'beta_isolating_interval':['183/100','184/100'],'polynomial_endpoint_values':[str(evalpoly(3,lower)),str(evalpoly(3,upper))],
        'certified_contraction_upper_bound':'3/4 (strict)','all_compact_components_radius_upper_bound':'4 (strict)',
        'depth12_hausdorff_error_strict_upper_bound':str(4*F(3,4)**12),
        'inverse16_power_basis':encode(v),'inverse16_substitution_basis':encode(substitution_coordinates),'inverse16_star_modulus':'beta^8 > (9/5)^8 > 100; outside all components'}

def main():
    R,sets,graph=graph_certificate()
    out={'status':'FINITE-EXACT checks; general claims require accompanying proofs','family':[family_certificate(k) for k in (3,4,5)],
      'graph':graph,'feedback':feedback_certificate(R,sets),'bounds':bounds_and_units(R),
      'not_certified':['Rauzy tiling or disjointness','pure dynamical spectrum','Amanda threshold dynamics conjugate to substitution','field embeddings between degrees3,4,5','original arithmetic module preserved by nonzero shear']}
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
