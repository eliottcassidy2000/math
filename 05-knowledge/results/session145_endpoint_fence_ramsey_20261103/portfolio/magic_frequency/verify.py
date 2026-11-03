#!/usr/bin/env python3
"""Exact nine-label/Latin/Frobenius bridge. No assertions: valid under python -O."""
from fractions import Fraction as Q
from itertools import permutations, product
from collections import Counter
from pathlib import Path
import argparse, copy, hashlib, json

FREQUENCIES = (174,285,396,417,528,639,741,852,963)
TARGET = (108,148,122,140,126,112,130,104,144)
RANKS = (2,9,4,7,5,3,6,1,8)
U = (0,1,-1,-1,0,1,1,-1,0)
V = (-1,1,0,1,0,-1,0,-1,1)
T = (-3,2,4,1,0,-1,-4,-2,3)
LINES = ((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))

def require(pred, msg):
    if not pred: raise ValueError(msg)

def sums(grid):
    require(len(grid)==9, 'grid length')
    return [sum(grid[i] for i in line) for line in LINES]

def magic(grid):
    return len(set(sums(grid))) == 1

def rank(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    if not a: return 0
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        d=a[r][j]; a[r]=[x/d for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:
                d=a[i][j]; a[i]=[x-d*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r

def full_constraint_matrix():
    return [[int(i in line)-int(i in LINES[0]) for i in range(9)] for line in LINES[1:]]

def magic_grids_constructive(values):
    c = Q(sum(values),9)
    if c not in values: return set()
    out=set()
    for a,b in product(values,repeat=2):
        # a,b are centered corner parameters, not input values.
        x=a-c; y=b-c
        grid=(c+x,c-x-y,c+y,c+y-x,c,c+x-y,c-y,c+x+y,c-x)
        if Counter(grid)==Counter(values): out.add(grid)
    return out

def magic_grids_bruteforce(values):
    c=Q(sum(values),9)
    if c not in values: return set()
    rem=list(values); rem.remove(c)
    out=set()
    for p in permutations(rem):
        grid=p[:4]+(c,)+p[4:]
        if magic(grid): out.add(grid)
    return out

def node_polynomial(t):
    ans=t
    for a in range(1,5): ans*=t*t-a*a
    return ans

def u_polynomial(t):
    return Q(t*(t*t-9)*(3*t**4-57*t*t+124),560)

def v_polynomial(t):
    return -Q(t*(t*t-16)*(t**4-12*t*t+67),840)

def frequency(a):
    return 100*a+10*((a-4)%9+1)+((a-7)%9+1)

def digit_coords(a):
    return ((a-1)%3-1,(a-1)//3-1)

def order_map(a,c=126,p=4,q=18):
    u,v=digit_coords(a)
    return c+p*u+q*v

def trace_map(a):
    t=a-5
    return 126-t**3+15*t

def scalar(x):
    if isinstance(x,Q): return int(x) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
    if isinstance(x,dict): return {k:scalar(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [scalar(v) for v in x]
    return x

def verify():
    require(tuple(frequency(a) for a in range(1,10))==FREQUENCIES,'digit encoding')
    require(sum(FREQUENCIES)==4995 and Q(sum(FREQUENCIES),9)==555,'frequency mean')
    require(555 not in FREQUENCIES,'frequency no-magic center')
    require(sums(TARGET)==[378]*8,'target magic')
    require(tuple(126-t**3+15*t for t in T)==TARGET,'Frobenius evaluation')
    require(tuple(sorted(trace_map(a) for a in range(1,10)))==tuple(sorted(TARGET)),'trace bijection')
    require(tuple(order_map(a) for a in range(1,10))==tuple(sorted(TARGET)),'order bijection')
    require(tuple(order_map(a) for a in RANKS)==TARGET,'order placement')
    require(set(zip(U,V))==set(product((-1,0,1),repeat=2)),'Latin orthogonality')
    for z in (U,V):
        require(sums(z)==[0]*8,'zero-sum coordinate')
        for line in LINES[:6]: require({z[i] for i in line}=={-1,0,1},'Latin row/column')
    near=tuple(FREQUENCIES[a-1] for a in RANKS)
    require(sums(near)==[1665]*7+[1584],'seven-line near-magic')
    require(all(f==528+111*u+Q(567,2)*v+Q(81,2)*v*v for f,u,v in zip(near,U,V)),'carry curvature')
    require(all(n==126+4*u+18*v for n,u,v in zip(TARGET,U,V)),'target chart')
    require(rank(full_constraint_matrix())==6,'full magic constraints rank six')
    require(rank([list((1,u,v)) for u,v in zip(U,V)])==3,'magic basis independent')
    # Tensor monomials form a complete nine-dimensional function basis.
    require(rank([[u**i*v**j for i,j in product(range(3),repeat=2)] for u,v in zip(U,V)])==9,'tensor basis')
    projection=tuple(555+111*u+Q(567,2)*v for u,v in zip(U,V))
    residual=tuple(f-g for f,g in zip(near,projection))
    require(magic(projection),'orthogonal projection magic')
    for z in ((1,)*9,U,V): require(sum(r*c for r,c in zip(residual,z))==0,'projection orthogonality')
    require(sum(r*r for r in residual)==Q(6561,2),'projection squared error')
    require(tuple(u_polynomial(t) for t in T)==U,'degree-seven U interpolant')
    require(tuple(v_polynomial(t) for t in T)==V,'degree-seven V interpolant')
    require(all(4*u_polynomial(t)+18*v_polynomial(t)==-t**3+15*t for t in range(-6,7)),'polynomial cubic cancellation')
    require(all(node_polynomial(t)==0 for t in range(-4,5)),'node polynomial')
    degree_ranks=[]
    for degree in range(9):
        A=[[sum(T[i]**k for i in line)-sum(T[i]**k for i in LINES[0]) for k in range(degree+1)] for line in LINES[1:]]
        degree_ranks.append(rank(A))
    require(degree_ranks==[0,1,2,2,3,4,5,5,6],'fixed-placement polynomial dimensions')
    constructive=magic_grids_constructive(TARGET)
    exhaustive=magic_grids_bruteforce(TARGET)
    require(constructive==exhaustive and len(exhaustive)==8,'two-path arrangement census')
    require(not magic_grids_constructive(FREQUENCIES),'frequency arrangements empty')
    fibers={str(r):[trace_map(a) for a in range(1,10) if a%3==r] for r in range(3)}
    require([sorted(fibers[str(r)]) for r in (1,2,0)]==[[112,130,148],[108,126,144],[104,122,140]],'trace fiber APs')
    # The digit-rotation C3 actions preserve fibers, not numerical magic.
    rotate=lambda a:(a-4)%9+1
    rotated=tuple(order_map(rotate(a)) for a in RANKS)
    require(not magic(rotated),'transported digit rotation is not magic symmetry')
    for a in range(1,10):
        b=rotate(a)
        require(frequency(b)==10*(frequency(a)%100)+frequency(a)//100,'decimal cyclic rotation')
        require(a%3==b%3 and rotate(rotate(rotate(a)))==a,'rotation fiber action')
    return scalar({
      'status':'PASS','universe':'Nine fixed frequency labels and all 3x3 placements; real-valued maps on those nine states.',
      'frequencies':FREQUENCIES,'frequency_mean':555,'frequency_magic_arrangements':0,
      'target':TARGET,'target_line_sums':sums(TARGET),'target_magic_arrangements':8,
      'rank_placement':RANKS,'U':U,'V':V,'frequency_near_grid':near,'near_line_sums':sums(near),
      'trace_map':[trace_map(a) for a in range(1,10)],'order_map':[order_map(a) for a in range(1,10)],
      'trace_fibers':fibers,'order_map_rotation_line_sums':sums(rotated),
      'magic_constraint_rank':6,'magic_function_dimension':3,'polynomial_constraint_ranks_degree_0_to_8':degree_ranks,
      'nearest_magic_grid':projection,'projection_residual':residual,'projection_squared_error':Q(6561,2),
      'source_revision':'7b76044ca13e6684d1746923285b83dbf83828d2'
    })

def validate_bundle(data):
    expected=verify()
    require(data==expected,'bundle differs from independently recomputed exact result')
    return True

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    result=verify();out=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(out)
    print(out,end='')
