"""Exact finite head and symbolic chord certificates for two owner lines.

Universe: every primitive sorted distinct positive odd ternary-unit triple
with c<=51. The all-height geometry is proved in the companion report.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
from functools import reduce
import sys

CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def cross(u, v):
    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])


def carriers(w):
    # Independent x,z scan, solving the middle coordinate exactly.
    a, b, c = w
    bx, by, bz = ((3*(sum(w)-wi)-1)//14 for wi in w)
    out = set()
    for x in range(-bx, bx+1):
        if not x % 3:
            continue
        for z in range(-bz, bz+1):
            if not z % 3 or (a*x+c*z) % b:
                continue
            y = -(a*x+c*z)//b
            if abs(y) <= by and y % 3:
                out.add((x,y,z))
    return out


def owner_line(w, points):
    if not points or all(cross(u,v)==(0,0,0) for u,v in combinations(points,2)):
        return None
    A = sorted(u for u in points if w[0]*u[0] % 3 == 1)
    need(len(A)*2==len(points) and len(A)>=2, ('two owners',w))
    diff = tuple(A[1][i]-A[0][i] for i in range(3))
    g = reduce(gcd, map(abs,diff))
    r = tuple(x//g for x in diff)
    if next(x for x in r if x) < 0:
        r = tuple(-x for x in r)
    if any(cross(r,tuple(u[i]-A[0][i] for i in range(3)))!=(0,0,0) for u in A):
        return None
    need(any(x%3==0 for x in r), ('invisible primitive step',w,r))
    need(g==3, ('consecutive full owner-line spacing',w,r,g))
    X = cross(A[0],r)
    need(X==w or X==tuple(-x for x in w), ('saturated transverse basis',w,r,X))
    M = max(map(abs,r))
    need(M>=3, ('no smaller invisible relation',w,r))
    if M==3:
        need(sorted(map(abs,r))==[1,2,3], ('small direction classification',w,r))
        need(Q(len(points)) < Q(2*w[2],21)+2, ('small direction chord bound',w,r))
    need(Q(len(points)) < Q(w[2],7)+2, ('universal owner-line tail bound',w,r))
    return r


def affine_add(f,g):
    return (f[0]+g[0],f[1]+g[1])


def affine_scale(f,t):
    return (f[0]*t,f[1]*t)


def chord_rows():
    # Normalize c=1 and a=-(r2*b+r3)/r1. Recover every possible sign/order
    # case directly from 0<=a<=b<=1. Endpoints are allowed for an upper bound.
    rows=[]
    for r in product(range(-3,4),repeat=3):
        if r[0]<=0 or 0 in r or sorted(map(abs,r))!=[1,2,3]:
            continue
        a=(Q(-r[1],r[0]),Q(-r[2],r[0])); b=(Q(1),Q(0)); one=(Q(0),Q(1))
        lo,hi=Q(0),Q(1)
        for slope,intercept in (a,affine_add(b,affine_scale(a,-1)),affine_add(one,affine_scale(b,-1))):
            if slope>0:
                lo=max(lo,-intercept/slope)
            elif slope<0:
                hi=min(hi,-intercept/slope)
            elif intercept<0:
                lo=Q(2)
        if lo>=hi:
            continue
        u=(affine_scale(b,Q(-1,r[2])),affine_scale(a,Q(1,r[2])),(Q(0),Q(0)))
        ws=(a,b,one); lower=[]; upper=[]
        for i in range(3):
            roof=affine_scale(affine_add(ws[(i+1)%3],ws[(i+2)%3]),Q(3,14))
            L=affine_scale(affine_add(affine_scale(roof,-1),affine_scale(u[i],-1)),Q(1,r[i]))
            U=affine_scale(affine_add(roof,affine_scale(u[i],-1)),Q(1,r[i]))
            if r[i]<0:
                L,U=U,L
            lower.append(L); upper.append(U)
        certs=[]
        for i,j in product(range(3),repeat=2):
            slope,intercept=affine_add(upper[i],affine_scale(lower[j],-1))
            bound=max(slope*lo+intercept,slope*hi+intercept)
            certs.append((bound,i+1,j+1,slope,intercept))
        best=min(certs)
        need(best[0]<=Q(1,7), ('global chord certificate',r,best))
        rows.append((r,(str(lo),str(hi)),tuple(map(str,best))))
    need(len(rows)==7, 'complete small-direction sign/order classification')
    return rows


def main():
    sys.stdout.reconfigure(newline='\n')
    rows=chord_rows()
    eligible=hits=0; best=(Q(0),None); sizes={}
    for w in combinations([x for x in range(1,52,2) if x%3],3):
        if gcd(gcd(w[0],w[1]),w[2])!=1:
            continue
        eligible+=1
        points=carriers(w)
        r=owner_line(w,points)
        if r is None:
            continue
        hits+=1; sizes[len(points)]=sizes.get(len(points),0)+1
        ratio=Q(len(points),w[2])
        if ratio>best[0]:
            best=(ratio,w)
        need(Q(len(points))<Q(2*w[2],11), ('complete strict low count gate',w,r,len(points)))
    named=[]
    for w in ((17,23,25),(19,23,29),(23,29,37),(41,47,49),(1,137,277),(1,499,1001)):
        points=carriers(w); r=owner_line(w,points)
        named.append((w,len(points),r))
        need((r is None)==(w in ((19,23,29),(23,29,37))), ('positive and excluded controls',w))
    need(Q(53,7)+2<Q(106,11),'strict all-height tail onset')
    print('PROVED_ANALYTICALLY + FINITE_EXACT: rank-two full supports with each owner color collinear are count-safe')
    print('finite_head=all primitive sorted distinct odd ternary-unit triples c<=51')
    print('eligible='+str(eligible)+'; owner_line_rows='+str(hits)+'; N_counts='+str(sorted(sizes.items())))
    print('finite_max_N_over_c='+str(best))
    print('small_direction_chord_certificates=(r,b_interval,(max_bound,upper_i,lower_j,slope,constant))')
    for row in rows:
        print(row)
    print('named_controls='+str(named))
    print('all_height_basis=u_cross_r=+-w; owner_line_step=3r')
    print('all_height_count=N<c/7+2; c>=53 implies N<2c/11')
    print('finite_low_count_gate=PASS; general_full_support_count_gate=OPEN; LRC14=OPEN')
    print('checks='+str(CHECKS))


if __name__=='__main__':
    main()
