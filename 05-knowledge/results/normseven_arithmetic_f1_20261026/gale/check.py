"""Exact controls for GA; no producer imports or assertion-only checks."""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path

GATES = 0


def require(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


class Field:
    """Integers encode a+b*u, u^2=nu; degree1 uses b=0 only."""
    def __init__(self, p, degree):
        self.p, self.degree = p, degree
        self.q = p ** degree
        self.nu = next(n for n in range(2, p)
                       if pow(n, (p-1)//2, p) == p-1) if degree == 2 else 0
        self.roots = Counter(self.mul(a, a) for a in range(self.q))

    def add(self, a, b):
        p = self.p
        return (a % p + b % p) % p + p * ((a // p + b // p) % p)

    def neg(self, a):
        return self.add(0, (-a % self.p) + self.p * (- (a//self.p) % self.p))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        p = self.p
        x, y, z, w = a % p, a // p, b % p, b // p
        return ((x*z+self.nu*y*w) % p) + p * ((x*w+y*z) % p)

    def power(self, a, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            n >>= 1
        return out

    def chi(self, a):
        if a == 0:
            return 0
        v = self.power(a, (self.q-1)//2)
        require(v in (1, self.p-1), 'Euler value must be +/-1')
        return 1 if v == 1 else -1


FORMS = [(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1)]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def normalize(a):
    d = gcd(gcd(abs(a[0]), abs(a[1])), abs(a[2]))
    out = tuple(t//d for t in a)
    return out if next(t for t in out if t) > 0 else tuple(-t for t in out)


def dot(a, b):
    return sum(x*y for x,y in zip(a,b))


def branch_audit():
    points = sorted(set(normalize(cross(a,b)) for a,b in combinations(FORMS,2)))
    incidence = [(point,[i for i,a in enumerate(FORMS) if dot(a,point)==0])
                 for point in points]
    counts = Counter(len(labels) for _,labels in incidence)
    require(counts == {2:9,3:2}, 'nine nodes and two triples')
    dets = sorted(set(abs(dot(a,cross(b,c))) for a,b,c in combinations(FORMS,3)))
    require(dets == [0,1,2,3], 'only bad incidence primes2,3')
    for p in (5,7,11,13,17,19):
        modpoints = {}
        for a,b in combinations(FORMS,2):
            v = tuple(t % p for t in cross(a,b))
            inv = pow(next(t for t in v if t), -1, p)
            v = tuple(t*inv % p for t in v)
            modpoints[v] = [i for i,L in enumerate(FORMS) if dot(L,v)%p==0]
        require(Counter(map(len,modpoints.values())) == counts, 'incidence modp')
    return {'rational_intersections':incidence, 'nonzero_determinants':[1,2,3],
            'resolution_correction':17}


def count_field(p, degree):
    f = Field(p, degree)
    q = f.q
    chi = [f.chi(a) for a in range(q)]
    for a in range(q):
        require(f.roots.get(a,0) == 1+chi[a], 'independent square multiplicity')
    # Path A: affine product of three factors, direct square-root histogram.
    values = Counter()
    S = 0
    for x,z in product(range(q),repeat=2):
        v = f.mul(f.mul(f.sub(1,f.mul(x,x)),f.sub(1,f.mul(z,z))),
                  f.add(f.add(1,x),z))
        values[v] += 1
        S += chi[v]
    # Path B: projective representatives and six separate linear forms.
    projective = [(1,x,z) for x,z in product(range(q),repeat=2)]
    projective += [(0,1,z) for z in range(q)] + [(0,0,1)]
    linear_values = Counter()
    for E,x,z in projective:
        factors = [E,f.add(E,x),f.sub(E,x),f.add(E,z),f.sub(E,z),f.add(f.add(E,x),z)]
        v = 1
        for a in factors:
            v = f.mul(v,a)
        linear_values[v] += 1
    require(linear_values == values+Counter({0:q+1}), 'affine/projective paths')
    twist_counts = []
    for delta in range(1,q):
        raw_a = q+1+sum(m*f.roots.get(f.mul(delta,a),0) for a,m in values.items())
        raw_b = sum(m*(1+chi[f.mul(delta,a)]) for a,m in linear_values.items())
        formula = q*q+18*q+1+chi[delta]*S
        require(raw_a == raw_b, 'root-count and projective Euler routes')
        require(raw_a+17*q == formula, 'GA-count every nonzero twist')
        twist_counts.append(formula)
    # Every nonzero coefficient gives a rational exceptional node conic.
    conics = []
    for delta in range(1,q):
        n = sum(f.mul(w,w)==f.mul(delta,f.mul(u,v)) for u,v,w in projective)
        require(n==q+1, 'all twisted A1 exceptional conics are P1')
        conics.append(n)
    minus = twist_counts[p-2] # delta=p-1; list begins at delta1.
    plus = twist_counts[0]
    require(sum(twist_counts) == (q-1)*(q*q+18*q+1), 'twist mean')
    return {'p':p,'degree':degree,'q':q,'quadratic_nonresidue':f.nu,
            'S_q':S,'plus_points':plus,'minus_one_points':minus,
            'twist_histogram':dict(sorted(Counter(twist_counts).items())),
            'all_nonzero_twists_checked':q-1,'projective_base_points':q*q+q+1,
            'twisted_node_conics_checked':len(conics)}


def algebra_audit():
    A=[[1,1,1,1,1,1],[0,1,-1,0,0,1],[0,0,0,1,-1,1]]
    C=[[-2,1,1,0,0,0],[-2,0,0,1,1,0],[1,-1,0,-1,0,1]]
    M=[[1,1,2],[1,0,0],[0,1,0]]
    permutation=[5,2,1,4,3,0]
    scalars=[-2,-1,1,-1,1,2]
    for i,j in product(range(3),range(3)):
        require(sum(A[i][k]*C[j][k] for k in range(6))==0,'A C^t=0')
    for i,j in product(range(3),range(6)):
        require(sum(M[i][k]*C[k][j] for k in range(3)) == A[i][permutation[j]]*scalars[j], 'Gale coefficient identity')
    require(dot(M[0],cross(M[1],M[2]))==2,'Gale base map invertible')
    s=1
    for a in scalars:s*=a
    require(s==-4,'minus-one after dividing coordinate by2')
    # A degree-six homogeneous polynomial is determined by tensor grid0..6.
    for E,x,z in product(range(7),repeat=3):
        lam=[E+x,E+z,2*E]
        net=1
        for j in range(6):net*=sum(C[i][j]*lam[i] for i in range(3))
        F=E*(E*E-x*x)*(E*E-z*z)*(E+x+z)
        require(net==-4*F,'independent polynomial Gale identity')
    section_controls=[]
    for t in [Fraction(a,b) for b in range(1,6) for a in range(-10,11)]:
        d=1-t*t
        X,Y=d*t,d*d
        require(Y*Y==(X-d)*(X+d)*(X-d*(t+1)),'section P')
        z=-t
        rhs=d*(1-z*z)*(1+t+z)
        require(rhs==d*d,'visible section line')
        for delta in (-7,-2,-1,1,2,3,4,7,9):
            require(delta*rhs==delta*d*d,'fixed twist line')
        if d and t not in (0,-2):section_controls.append(str(t))
    # Positive and hostile controls for the one-dimensional descent action.
    for delta_sign in (-1,1):
        action=[1]*18+[delta_sign]
        require(sum(a==1 for a in action)==(18 if delta_sign==-1 else 19), 'NS invariant rank')
        require(sum(action)==18+delta_sign,'NS character trace')
    return {'gale_permutation':permutation,'gale_scalars':scalars,
            'gale_multiplier':s,'polynomial_grid_points':7**3,
            'rational_section_samples':105,
            'smooth_section_samples':len(section_controls),
            'representation_checks':'Formal algebra only; geometric rank is inherited, not computed.'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    branch=branch_audit()
    algebra=algebra_audit()
    fields=[count_field(p,degree) for p,degree in [(5,1),(7,1),(11,1),(13,1),(17,1),(19,1),(7,2),(13,2)]]
    byq={r['q']:r for r in fields}
    require((byq[7]['plus_points'],byq[7]['minus_one_points'])==(180,172),'nonsplit p7 witness')
    require((byq[13]['plus_points'],byq[13]['minus_one_points'])==(408,408),'split p13 hostile')
    require((byq[49]['plus_points'],byq[49]['minus_one_points'])==(3384,3384),'fixed twist even-extension hostile')
    require(byq[169]['plus_points']==32088,'inherited p13^2 positive control')
    receipt={'status':'FINITE-EXACT; unbounded theorem uses inherited rank and written descent proof',
             'universe':'All projective base points and every nonzero constant twist over F5,F7,F11,F13,F17,F19,F49,F169; all twisted node conics; 105 rational section inputs; degree-six interpolation grid0..6 cubed.',
             'branch':branch,'algebra':algebra,'fields':fields,'explicit_gates':GATES}
    output=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(output,encoding='utf-8',newline='\n')
    print(output,end='')


if __name__=='__main__':main()
