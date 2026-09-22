"""Exact controls for the Collatz-edge triangle encoding and sharp angle gaps."""
from fractions import Fraction as F
from hashlib import sha256
from math import atan, degrees, gcd, isqrt, pi, sqrt
from pathlib import Path
import json

CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def step(x, b):
    z = 3*x+b
    if z <= 0:
        raise ValueError('positive numerator required')
    k = (z & -z).bit_length()-1
    return z >> k, k


def triangle(x, y):
    return x*y, abs(x*x-y*y)//2, (x*x+y*y)//2


def tangent(x, y):
    a, b, _ = triangle(x, y)
    return F(min(a,b), max(a,b))


def branch_count(X, b):
    maximum_k = (3*isqrt(2*X)+1).bit_length()-1
    total = 0
    counts = {}
    for k in range(1,maximum_k+1):
        A = 4**k+9
        D = 4**k*(2*X*A-b*b)
        endpoint = (-3*b+isqrt(D))//A
        modulus = 2**(k+1)
        residue = ((2**k-b)*pow(3,-1,modulus)) % modulus
        check(residue % 2 == 1, 'valuation residue odd')
        number = max(0,(endpoint-residue)//modulus+1)
        if number:
            first_y, first_k = step(residue,b)
            last = residue+(number-1)*modulus
            last_y,last_k = step(last,b)
            check(first_k == last_k == k, 'branch endpoints have exact valuation')
            check(triangle(last,last_y)[2] <= X, 'last residue fits height')
            nxt = last+modulus
            check(triangle(nxt,step(nxt,b)[0])[2] > X, 'next residue exceeds height')
        counts[k] = number
        total += number
    # Both parameters have fixed positive state1, whose triangle is degenerate.
    if X >= 1:
        counts[2 if b==1 else 1] -= 1
        total -= 1
    return total,counts


def enumerate_edges(X,b):
    edges=[]
    counts={}
    for x in range(1,isqrt(2*X)+1,2):
        y,k=step(x,b)
        if x!=y and triangle(x,y)[2]<=X:
            edges.append((x,y))
            counts[k]=counts.get(k,0)+1
    return edges,counts


def ppt_count_roots(X):
    return sum(gcd(u,v)==1 and (u*u+v*v)//2<=X
               for u in range(3,isqrt(2*X)+1,2) for v in range(1,u,2))


def ppt_count_euclid(X):
    return sum(gcd(m,n)==1 and (m-n)%2==1 and m*m+n*n<=X
               for m in range(2,isqrt(X)+1) for n in range(1,m))


def main():
    extrema=[]
    for parameter in (-1,1):
        top=F(0)
        witnesses=[]
        row_counts={}
        for x in range(1,100002,2):
            y,k=step(x,parameter)
            check(k>=1 and y%2==1 and y>0, 'legal accelerated edge')
            check(gcd(x,y)==1, 'unit-parameter edge coprime')
            if x==y:
                check(x==1, 'only fixed state in finite positive universe')
                continue
            a,b,c=triangle(x,y)
            check(a*a+b*b==c*c and min(a,b)>0, 'right triangle')
            check(gcd(gcd(a,b),c)==1 and a%2==1 and b%2==0, 'primitive parity')
            check(c+b==max(x,y)**2 and c-b==min(x,y)**2, 'recover both roots')
            z,_=step(y,parameter)
            check(2*(triangle(y,z)[2]-c)==z*z-x*x, 'two-step height transport')
            tan=tangent(x,y)
            check(tan<=F(65,72) if parameter==1 else tan<F(48,55), 'global angle gate')
            row=min(k,4)
            row_counts[row]=row_counts.get(row,0)+1
            if parameter==1:
                bounds={1:F(8,15),2:F(7,24),3:F(65,72),4:F(8,15)}
                check(tan<=bounds[row], 'positive branch bound')
                if k==3:
                    check(x%16==13 and x>=13, 'sharp positive guard')
            else:
                bounds={1:F(5,12),2:F(12,35),3:F(48,55),4:F(96,247)}
                check(tan<=bounds[row], 'negative branch bound')
                if k==3:
                    check(x%16==3 and x>=3, 'negative supremum guard')
            if tan>top:
                top,witnesses=tan,[[x,y,k]]
            elif tan==top:
                witnesses.append([x,y,k])
        extrema.append({'parameter':parameter,'finite_max_tangent':str(top),
                        'witnesses':witnesses,'branch_populations_k4plus_combined':row_counts})
    check(step(13,1)==(5,3) and tangent(13,5)==F(65,72), 'sharp positive witness')
    check(triangle(13,5)==(65,72,97), 'sharp witness sides')
    previous=F(0)
    for j in range(1000):
        x=3+16*j
        y,k=step(x,-1)
        t=tangent(x,y)
        check(k==3 and previous<t<F(48,55), 'negative limit family increasing')
        previous=t

    check(step(7,1)[0]==11 and triangle(7,11)==(77,36,85), 'legal child')
    check((2*7-3,7)==(11,7), 'Berggren parent roots')
    check(step(3,1)[0]!=7 and step(7,1)[0]!=3, 'parent loses positive legality')
    check(step(5,-1)[0]==7 and abs(7-2*5)==3, 'shell move on legal minus edge')
    check(step(3,-1)[0]!=7 and step(7,-1)[0]!=3, 'shell move loses minus legality')

    heights=[1,5,13,29,97,100,1000,10000,100000]
    census=[]
    for X in heights:
        p=ppt_count_roots(X)
        check(p==ppt_count_euclid(X), 'two independent PPT parametrizations')
        row={'hypotenuse_bound':X,'PPT_count':p,'edges':[]}
        for b in (-1,1):
            edges,by_start=enumerate_edges(X,b)
            count,by_residue=branch_count(X,b)
            check(count==len(edges), 'independent edge-height counts')
            check({k:v for k,v in by_residue.items() if v}==by_start, 'branch counts agree')
            unique={triangle(x,y) for x,y in edges}
            check(len(unique)<=len(edges)<=2*len(unique), 'direction-loss fibre bound')
            row['edges'].append({'parameter':b,'directed_count':count,
                                 'unmarked_count':len(unique), 'by_exponent':by_start})
        census.append(row)
    return {'status':'PASS','source_sha256_lf':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
            'checks_passed':CHECKS,'universe':{'parameters':[-1,1],
            'positive_odd_starts':[1,100001],'negative_limit_family_j':[0,999],
            'hypotenuse_bounds':heights,'exclude_fixed_edges':True},
            'extrema':extrema,'height_census':census,
            'illustrative_decimals_not_proof_gates':{
                'C_edge_first99terms':sum(1/sqrt(2*(4**k+9)) for k in range(1,100)),
                'theta_plus_degrees':degrees(atan(F(65,72))),
                'theta_minus_sup_degrees':degrees(atan(F(48,55))),
                'ambient_PPT_leading_coefficient':1/(2*pi)},
            'infinite_claims_depend_on_written_proofs':True,'collatz_convergence_proved':False}


if __name__=='__main__':
    result=main()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',
                                                 encoding='utf-8',newline='\n')
    print('PASS:',result['checks_passed'],'exact edge/triangle controls')
