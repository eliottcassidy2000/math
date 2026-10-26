#!/usr/bin/env python3
"""Independent exhaustive small-universe audit. No imported repository code.
All operations exact; no assertions (the verifier remains active with python -O).
"""
import ast
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

CHECKS = 0

def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

def parity(x):
    return x.bit_count() & 1

def bareiss(matrix):
    a = [r[:] for r in matrix]
    n = len(a)
    if not n:
        return 1
    sign = 1
    prev = 1
    for k in range(n-1):
        pivot = next((i for i in range(k,n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        p = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num = a[i][j]*p-a[i][k]*a[k][j]
                check(num % prev == 0, 'Bareiss exact division')
                a[i][j] = num//prev
            a[i][k] = 0
        prev = p
    return sign*a[-1][-1]

class Universe:
    def __init__(self, n):
        self.n=n
        self.edges=list(itertools.combinations(range(n),2))
        self.m=len(self.edges)
        self.size=1<<self.m
        self.perms=list(itertools.permutations(range(n)))
        idx={e:i for i,e in enumerate(self.edges)}
        self.linear=[]
        self.flips=[]
        self.invsrc=[]
        for p in self.perms:
            eperm=[idx[tuple(sorted((p[u],p[v])))] for u,v in self.edges]
            flip=sum(1<<eperm[i] for i,(u,v) in enumerate(self.edges) if p[u]>p[v])
            invsrc=sum(1<<i for i,(u,v) in enumerate(self.edges) if p[u]>p[v])
            action=[0]*self.size
            for x in range(1,self.size):
                b=x & -x
                action[x]=action[x^b]^(1<<eperm[b.bit_length()-1])
            self.linear.append(action)
            self.flips.append(flip)
            self.invsrc.append(invsrc)
        self.gorbits,self.gcanon=self.orbits(False)
        self.torbits,self.tcanon=self.orbits(True)
        self.even=[]
        self.even_orbits={}
        for x,orb in self.gorbits.items():
            if all(not parity(x & self.invsrc[i]) for i,a in enumerate(self.linear) if a[x]==x):
                self.even.append(x)
                self.even_orbits[x]=orb
        self.euler=[x for x in range(self.size) if all(not parity(x & sum(1<<i for i,e in enumerate(self.edges) if v in e)) for v in range(n))]
    def orbits(self, affine):
        reps={}
        canon=[None]*self.size
        for x in range(self.size):
            if canon[x] is not None:
                continue
            orb={a[x]^(self.flips[i] if affine else 0) for i,a in enumerate(self.linear)}
            check(min(orb)==x,'orbit minimum')
            reps[x]=sorted(orb)
            for y in orb: canon[y]=x
        return reps,canon
    def aut_t(self,t):
        return [i for i,a in enumerate(self.linear) if a[t]^self.flips[i]==t]
    def graph_edges(self,x):
        return [e for i,e in enumerate(self.edges) if x>>i&1]
    def census(self):
        ts=list(self.torbits)
        es=self.even
        even_labelled=set(itertools.chain.from_iterable(self.even_orbits.values()))
        check(len(ts)==len(es),'Royle count')
        euler_types={self.gcanon[x] for x in self.euler}
        ee_types=euler_types.intersection(es)
        # Independently compute signed orbit polynomials and all per-g trace identities.
        coeff=[0]*(self.m+1)
        for p,a in enumerate(self.linear):
            signed=[0]*(self.m+1)
            for x in range(self.size):
                if a[x]==x: signed[x.bit_count()]+=(-1 if parity(x & self.invsrc[p]) else 1)
            tfixed=sum(a[t]^self.flips[p]==t for t in range(self.size))
            check(sum(signed)==tfixed,'signed trace identity')
            coeff=[u+v for u,v in zip(coeff,signed)]
        expected=Counter(x.bit_count() for x in es)
        check(all(c==len(self.perms)*expected.get(k,0) for k,c in enumerate(coeff)),'graded identity')
        # Walsh invariant columns, by actual transported sign and direct evaluation.
        columns=[]
        for x in es:
            orbit_sign={}
            for p,a in enumerate(self.linear):
                y=a[x]
                sign=-1 if parity(x & self.invsrc[p]) else 1
                check(y not in orbit_sign or orbit_sign[y]==sign,'even sign well-defined')
                orbit_sign[y]=sign
            vals=[sum(s*(-1 if parity(t&y) else 1) for y,s in orbit_sign.items()) for t in range(self.size)]
            for t in range(self.size):
                check(vals[t]==vals[self.tcanon[t]],'Walsh relabel invariance')
                check(sum(vals[t^(1<<j)] for j in range(self.m))==(self.m-2*x.bit_count())*vals[t],'flip eigenfunction')
            columns.append([vals[t] for t in ts])
        M=[list(r) for r in zip(*columns)]
        det=bareiss(M)
        check(det != 0 and abs(det)&(abs(det)-1)==0,'Walsh determinant power of two')
        # Weighted native single-edge flip operator and exact intertwiner.
        index={t:i for i,t in enumerate(ts)}
        A=[[0]*len(ts) for t in ts]
        for i,t in enumerate(ts):
            for j in range(self.m): A[i][index[self.tcanon[t^(1<<j)]]]+=1
        for i in range(len(ts)):
            for j,x in enumerate(es):
                check(sum(A[i][k]*M[k][j] for k in range(len(ts)))==(self.m-2*x.bit_count())*M[i][j],'quotient spectral intertwiner')
        # All equivariant orbit-bijective candidates, via stabilizer containment.
        adjacency=[]
        for t in ts:
            H=self.aut_t(t)
            adjacency.append([j for j,x in enumerate(es) if any(all(self.linear[h][y]==y for h in H) for y in self.even_orbits[x])])
        matching={}
        def augment(i,seen):
            for j in adjacency[i]:
                if j in seen: continue
                seen.add(j)
                if j not in matching or augment(matching[j],seen):
                    matching[j]=i
                    return True
            return False
        for i in range(len(ts)): augment(i,set())
        hall=None
        if len(matching)<len(ts):
            for bits in range(1,1<<len(ts)):
                src=[i for i in range(len(ts)) if bits>>i&1]
                nbr=sorted(set(itertools.chain.from_iterable(adjacency[i] for i in src)))
                if len(nbr)<len(src):
                    hall={'source_representatives':[ts[i] for i in src], 'target_representatives':[es[j] for j in nbr], 'sizes':[len(src),len(nbr)],'source_automorphism_orders':[len(self.aut_t(ts[i])) for i in src]}
                    break
            check(hall is not None,'Hall failure certificate')
        hostile=None
        for x in sorted(even_labelled):
            for y in sorted(even_labelled):
                if x^y not in even_labelled:
                    hostile={'x':x,'y':y,'xor':x^y,'edge_sets':[self.graph_edges(z) for z in (x,y,x^y)]}
                    break
            if hostile: break
        # Cocycle obstruction at full S_n; positive control C_n, n odd.
        full_fixed=[t for t in range(self.size) if all(a[t]^self.flips[p]==t for p,a in enumerate(self.linear))]
        check(bool(full_fixed)==(self.n<2),'full symmetric affine fixed-point obstruction')
        cyclic_fixed=None
        if self.n and self.n%2:
            cyc=tuple((i+1)%self.n for i in range(self.n))
            p=self.perms.index(cyc)
            cyclic_fixed=[t for t in range(self.size) if self.linear[p][t]^self.flips[p]==t]
            check(bool(cyclic_fixed),'odd cyclic fixedpoint control')
            base=cyclic_fixed[0]
            for x in range(self.size):
                check((self.linear[p][x]^base)==(self.linear[p][x^base]^self.flips[p]),'odd cyclic equivariant torsor trivialization')
        return {'n':self.n,'labelled_tournaments':self.size,'labelled_sign_even_graphs':len(even_labelled),'labelled_eulerian_graphs':len(self.euler),'tournament_types':len(ts),'sign_even_graph_types':len(es),'eulerian_types':len(euler_types),'sign_even_eulerian_types':len(ee_types),'sign_even_edges_polynomial':[expected.get(k,0) for k in range(self.m+1)],'walsh_determinant':det,'native_flip_eigenvalues':[self.m-2*x.bit_count() for x in es],'hall_maximum_matching':len(matching),'hall_hostile':hall,'xor_closure_hostile':hostile,'odd_cyclic_fixed_tournaments':len(cyclic_fixed) if cyclic_fixed is not None else None}

def main():
    rows=[Universe(n).census() for n in range(6)]
    check([r['tournament_types'] for r in rows]==[1,1,1,2,4,12],'known census control')
    check([r['labelled_sign_even_graphs'] for r in rows]==[1,1,1,4,20,396],'labelled sign-even control')
    check(rows[4]['xor_closure_hostile'] is not None and rows[3]['xor_closure_hostile'] is None,'smallest XOR closure hostile')
    check(rows[5]['hall_maximum_matching']==11 and all(r['hall_maximum_matching']==r['tournament_types'] for r in rows[:5]),'smallest directional Hall hostile')
    tree=ast.parse(Path(__file__).read_text())
    check(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'no assert gates')
    check(not any(isinstance(x,ast.Constant) and isinstance(x.value,float) for x in ast.walk(tree)),'no floating arithmetic')
    print(json.dumps({'status':'PASS','universe':'All labelled tournaments and graphs on n=0..5, all vertex permutations; no filters','checks':CHECKS,'rows':rows},indent=2,sort_keys=True))

if __name__=='__main__': main()
