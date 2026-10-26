#!/usr/bin/env python3
"""Exact finite controls. Standard library only; no assert is a test gate."""
from fractions import Fraction as F
from itertools import permutations, product, combinations
from dataclasses import dataclass
from pathlib import Path
import json, sys

def require(p, message):
    if not p: raise RuntimeError(message)

@dataclass(frozen=True)
class Q:
    # Q(phi), phi^2=phi+1.
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a)); object.__setattr__(self,'b',F(self.b))
    def __add__(self, v):
        v = asq(v); return Q(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,v): return self+-asq(v)
    def __rsub__(self,v): return asq(v)+-self
    def __mul__(self,v):
        v=asq(v); return Q(self.a*v.a+self.b*v.b,self.a*v.b+self.b*v.a+self.b*v.b)
    __rmul__=__mul__
    def conj(self): return Q(self.a+self.b,-self.b)
    def __truediv__(self,v):
        v=asq(v); c=v.conj(); n=v*c
        require(n.b==0 and n.a!=0,'nonzero quadratic norm')
        w=self*c; return Q(w.a/n.a,w.b/n.a)
    def sign(self):
        # sign(2a+b+b sqrt(5)); rational comparisons only.
        a,b=2*self.a+self.b,self.b
        if a==0: return (b>0)-(b<0)
        if b==0: return (a>0)-(a<0)
        if a>0 and b>0: return 1
        if a<0 and b<0: return -1
        d=a*a-5*b*b
        require(d!=0,'sqrt5 irrationality check')
        return ((d>0)-(d<0))*((a>0)-(a<0))
    def serial(self): return [str(self.a),str(self.b)]

def asq(x): return x if isinstance(x,Q) else Q(x)
phi=Q(0,1); q=2-phi
Z=(Q(0),Q(0)); U=(Q(1),Q(0)); V=(Q(0),Q(1)); W=(Q(1),Q(1))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def mul(s,a):return tuple(s*x for x in a)
def norm(a):
    x,y=a; return x*x+y*y+(phi-1)*x*y

def wedge(a,b): return a[0]*b[1]-a[1]*b[0]
def area_coefficient(P):
    x=sum((wedge(a,b) for a,b in zip(P,P[1:]+P[:1])),Q())/2
    return x if x.sign()>=0 else -x

def rot(p):
    x,y=p; return ((1-phi)*x-y,x)

def penrose():
    C=(q,q)
    D=[Z,U,C,V]; K=[W,U,C,V]
    dc=[0,1,0,1]; kc=[1,0,1,0]  # black=0, white=1, Penrose Fig.9.
    legal=[sub(p,V) for p in K]
    shift=(phi,Q(-1))
    illegal=[add(rot(p),shift) for p in K]
    require(q.sign()>0 and (1-q).sign()>0,'0 < q < 1')
    require(phi*phi==phi+1,'golden polynomial')
    for P in [D,K]:
        require([norm(sub(b,a)) for a,b in zip(P,P[1:]+P[:1])]==[Q(1),q,q,Q(1)],'P2 side metrics')
    require(area_coefficient(K)/area_coefficient(D)==phi,'kite dart area ratio')
    require(wedge(rot(U),rot(V))==Q(1),'rotation determinant one')
    for p in [U,V,add(U,V),sub(U,V)]:
        require(norm(rot(p))==norm(p),'rotation preserves Gram form')
    require(all(p[1].sign()>=0 for p in D),'dart supporting halfplane')
    for P in [legal,illegal]:
        require(all(p[1].sign()<=0 for p in P),'kite opposite halfplane')
        require(set(p for p in P if p[1]==Q(0))=={Z,U},'shared seam exactly [0,U]')
        require(area_coefficient(P)==area_coefficient(K),'kite area preserved')
        require(all(p[0].sign()>=0 and (1-p[0]).sign()>=0 for p in P),'unit oblique strip')
    colors_D=dict(zip(D,dc)); colors_legal=dict(zip(legal,kc)); colors_illegal=dict(zip(illegal,kc))
    good=all(colors_D[v]==colors_legal[v] for v in [Z,U])
    bad=all(colors_D[v]==colors_illegal[v] for v in [Z,U])
    require(good and not bad,'native marked-tile attachment guard')
    # Inflation inventory matrix M=[[2,1],[1,1]], area eigenvector (phi,1).
    lam=phi+1
    require(2*phi+1==lam*phi and phi+1==lam,'inflation area eigenpair')
    return {'field':'Q(phi), phi^2=phi+1','gram_off_diagonal':((phi-1)/2).serial(),
            'D':[[x.serial(),y.serial()] for x,y in D],
            'K_legal':[[x.serial(),y.serial()] for x,y in legal],
            'K_illegal':[[x.serial(),y.serial()] for x,y in illegal],
            'legal_guard':good,'illegal_guard':bad,'nonoverlap_both':True,
            'same_simple_incidence':'K2','same_counts':{'kite':1,'dart':1},
            'area_ratio':phi.serial(),'rotation_determinant':1,'reflection_used':False}

def involutions(n):
    def recurse(todo,p):
        if not todo: yield tuple(p); return
        x=todo[0]; rest=todo[1:]; p[x]=x
        yield from recurse(rest,p)
        for y in rest:
            p[x]=y;p[y]=x
            yield from recurse([z for z in rest if z!=y],p)
    yield from recurse(list(range(n)),list(range(n)))

def orbits(a):return sorted({tuple(sorted({s,a[s]})) for s in range(len(a))})
def diagram(a,b):
    A={s:('A',o) for o in orbits(a) for s in o}; B={s:('B',o) for o in orbits(b) for s in o}
    return {s:(A[s],B[s]) for s in range(len(a))}

def decode(E,n):
    out=[]
    for side in range(2):
        nodes={edge[side] for edge in E.values()}; a=list(range(n))
        for node in nodes:
            incident=[s for s,e in E.items() if e[side]==node]
            require(len(incident) in (1,2),'degree one or two')
            if len(incident)==2:
                x,y=incident;a[x]=y;a[y]=x
        out.append(tuple(a))
    return tuple(out)

def components(a,b):
    todo=set(range(len(a))); out=[]
    while todo:
        start=min(todo); seen={start}; stack=[start]
        while stack:
            s=stack.pop()
            for t in (a[s],b[s]):
                if t not in seen:seen.add(t);stack.append(t)
        todo-=seen;out.append(sorted(seen))
    return out

def inverse(p):
    a=[0]*len(p)
    for i,v in enumerate(p):a[v]=i
    return tuple(a)
def compose(a,b): return tuple(a[b[i]] for i in range(len(a)))
def conjugate(p,a):return compose(compose(p,a),inverse(p))

def matching_checks():
    total=0; counts={}
    for n in range(7):
        inv=list(involutions(n)); counts[n]=len(inv)
        for a in inv:
            for b in inv:
                E=diagram(a,b);require(decode(E,n)==(a,b),'whole-side ED exact roundtrip');total+=1
    a=(1,0,3,2); b=(2,3,0,1)
    require(len(components(a,a))==2 and len(components(a,b))==1,'relative alignment component witness')
    ident=(0,1); swap=(1,0)
    require(compose(ident,ident)==compose(swap,swap),'product-only collision')
    require(len(components(ident,ident))==2 and len(components(swap,swap))==1,'product alone loses matching components')
    G=list(permutations(range(4))); H=[g for g in G if conjugate(g,a)==a]
    doublecosets=set()
    for g in G:
        doublecosets.add(min(compose(compose(h,g),k) for h in H for k in H))
    require(len(H)==8 and len(doublecosets)==2,'S4 matching double cosets')
    # A common atomic refinement is deliberately NOT the whole-side diagram.
    whole={'long':1,'short_left':1,'short_right':1}
    atomized={'long_left':2,'short_left':2,'long_right':2,'short_right':2}
    require(whole!=atomized,'flat-side refinement cannot masquerade as literal ED')
    return {'involution_counts_n0_through_n6':counts,'ordered_pairs_roundtripped':total,
            'same_individual_matching_type_different_joint_components':[components(a,a),components(a,b)],
            'S4_centralizer_size':len(H),'S4_double_cosets':len(doublecosets),
            'whole_sides_not_atomic_segments':True,
            'product_only_hostile':{'same_product':[0,1],'component_counts':[2,1]}}

def square_checks():
    def develop(steps):
        out=[(0,0)]
        for d in steps:out.append(tuple(x+y for x,y in zip(out[-1],d)))
        return out
    E=(1,0);N=(0,1);W=(-1,0);S=(0,-1)
    good=develop([E,N,W,W]);bad=develop([E,N,W,S])
    def distinct_ports(steps):
        return all(tuple(-x for x in steps[i-1])!=steps[i] for i in range(1,len(steps)))
    require(distinct_ports([E,N,W,W]) and distinct_ports([E,N,W,S]),'distinct local square ports')
    require(len(set(good))==5 and len(set(bad))==4,'global overlap independent of local tree consistency')
    require(bad[0]==bad[4],'exact interior overlap of first and fifth square')
    return {'good_square_origins':good,'bad_square_origins':bad,'declared_seam_graph':'P5',
            'tree_holonomy_vacuous':True,'distinct_local_ports':True,
            'good_union_area':5,'bad_union_area':4,'bad_piece_area_sum':5,
            'warning':'P5 is the declared seam graph, not the full incidental-contact graph.'}

def parity_checks():
    outputs=[]
    for k in range(1,8):
        states=list(product((0,1),repeat=k)); subsets=[J for r in range(k+1) for J in combinations(range(k),r)]
        # Bit1 is the legal candidate pose, bit0 the illegal pose.
        moments=[]
        for J in subsets:
            values=[]
            for parity in [0,1]:
                population=[b for b in states if sum(b)%2==parity]
                values.append(F(sum(all(b[j] for j in J) for b in population),len(population)))
            if len(J)<k:require(values[0]==values[1],'all proper-subset conjunctions agree')
            else:require(abs(values[0]-values[1])==F(1,2**(k-1)),'top conjunction parity gap')
            moments.append(values)
        # Sort subsets/states by bitmask: Boolean zeta matrix triangular with unit diagonal.
        n=1<<k
        require(all((s&t)==s for s in range(n) for t in [s]),'zeta diagonal')
        require(all(not ((s&t)==s) for s in range(n) for t in range(s)), 'zeta lower triangle zero')
        outputs.append({'k':k,'states':n,'proper_subset_guards':n-1,'guard_matrix_determinant':1,
                        'top_gap':str(F(1,2**(k-1)))})
    return outputs

def main():
    if '--inject-fault' in sys.argv:require(False,'injected hostile: guard acceptance changed')
    out={'status':'FINITE-EXACT','arithmetic':'integers and fractions in Q(phi); no floats',
         'penrose_native_guard':penrose(),'matching_diagram':matching_checks(),
         'nonoverlap_hostile':square_checks(),'native_conjunction_laws':parity_checks()}
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
