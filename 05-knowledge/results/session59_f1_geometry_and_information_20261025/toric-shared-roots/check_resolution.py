#!/usr/bin/env python3
"""Reproduce bounded checks supporting RESOLUTION.md."""
from itertools import product
from pathlib import Path
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
from check import Field, projective_tuples, pair_counts


def norm_projective(field, point):
    first=next(v for v in point if v)
    inverse=field.inv(first)
    return tuple(field.mul(x,inverse) for x in point)


def conic_param(field, s, t):
    s2,t2=field.square(s),field.square(t)
    middle=field.add(s2,field.neg(t2))
    twice_st=field.mul(field.constant(2),field.mul(s,t))
    return (field.add(s2,t2),field.add(middle,twice_st),
            field.add(middle,field.neg(twice_st)))


def field_geometry():
    rows=[]
    for field in [Field(p) for p in [3,5,7,11,13]]+[Field(p,2) for p in [3,5,7]]:
        q=field.q
        conic={point for point in projective_tuples(field,3)
               if field.add(field.square(point[1]),field.square(point[2])) ==
                  field.mul(field.constant(2),field.square(point[0]))}
        param={norm_projective(field,conic_param(field,*point))
               for point in projective_tuples(field,2)}
        assert param == conic and len(conic)==q+1
        boundary={point for point in conic if point[0]==0}
        eps=field.chi(field.neg(1))
        assert len(boundary)==1+eps
        basepoints=set(product(boundary,boundary))
        assert len(basepoints)==2*(1+eps)
        # Blowing up a smooth rational point replaces it by its tangent P^1.
        tangent_directions=list(projective_tuples(field,2))
        blowup_count=(q+1)**2-len(basepoints)+len(basepoints)*len(tangent_directions)
        predicted=q*q+(4+2*eps)*q+1
        assert blowup_count==predicted
        # Independent base count from the actual root equations.
        counts=pair_counts(field)
        base_count=(sum(n*n for n in counts)-1)//(q-1)
        assert blowup_count-base_count==q*2*(1+eps)
        # Verify the stated rational map on every rational non-basepoint.
        map_images=set()
        for first,second in product(conic,conic):
            c1,a,b=first; c2,d,e=second
            target=(field.mul(c1,c2),field.mul(a,c2),field.mul(b,c2),
                    field.mul(d,c1),field.mul(e,c1))
            if not any(target):
                assert (first,second) in basepoints
                continue
            c,A,B,D,E=target
            rhs=field.mul(field.constant(2),field.square(c))
            assert field.add(field.square(A),field.square(B))==rhs
            assert field.add(field.square(D),field.square(E))==rhs
            map_images.add(norm_projective(field,target))
        rows.append(dict(q=q,p=field.p,degree=field.degree,epsilon=eps,
                         conic=len(conic),basepoints=len(basepoints),
                         singular_base=base_count,resolution=blowup_count,
                         rational_map_image_before_blowup=len(map_images)))
    return rows


def fan_checks():
    rays=[(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]
    determinant=lambda u,v:u[0]*v[1]-u[1]*v[0]
    assert all(determinant(rays[j],rays[(j+1)%8])==1 for j in range(8))
    diagonal_rays=rays[1::2]
    assert all(determinant(diagonal_rays[j],diagonal_rays[(j+1)%4])==2
               for j in range(4))
    intersection=sp.zeros(8)
    self_intersections=[]
    for j in range(8):
        old,nxt=rays[(j-1)%8],rays[(j+1)%8]
        total=tuple(a+b for a,b in zip(old,nxt))
        coefficient=next(k for k in [1,2] if total==tuple(k*a for a in rays[j]))
        self_intersections.append(-coefficient)
        intersection[j,j]=-coefficient
        intersection[j,(j-1)%8]=intersection[j,(j+1)%8]=1
    assert self_intersections==[-2,-1]*4
    assert intersection.rank()==6
    assert sum(intersection)==4  # (-K)^2
    characters=sp.Matrix(rays)
    assert intersection*characters==sp.zeros(8,2)
    smith=smith_normal_form(characters,domain=ZZ)
    assert abs(smith[0,0])==abs(smith[1,1])==1
    involution=sp.zeros(8)
    for j in range(8):
        involution[(j+4)%8,j]=1
    assert involution*characters==-characters
    assert involution.T*intersection*involution==intersection
    x=sp.Symbol('x')
    ambient=sp.factor(involution.charpoly(x).as_expr())
    picard=sp.cancel(ambient/(x+1)**2)
    assert sp.expand(picard-(x-1)**4*(x+1)**2)==0
    traces=[]
    for n in range(1,9):
        trace=int((involution**n).trace())-2*(-1)**n
        assert trace==4+2*(-1)**n
        traces.append(trace)
    return dict(rays=rays,self_intersections=self_intersections,
                intersection_rank=6,anticanonical_square=4,
                picard_characteristic_polynomial=str(sp.factor(picard)),
                involution_power_traces=traces)


def symbolic_checks():
    s,t=sp.symbols('s t')
    c=s*s+t*t; a=s*s-t*t+2*s*t; b=s*s-t*t-2*s*t
    assert sp.expand(a*a+b*b-2*c*c)==0
    # In the w-chart, strict transform of wz=c^2 is Z=C^2.
    w,Z,C=sp.symbols('w Z C')
    assert sp.cancel((w*(w*Z)-(w*C)**2)/w**2)==Z-C*C
    q,tau=sp.symbols('q tau')
    # Formal power-sum coefficients of the displayed zeta denominator.
    for eps in [-1,1]:
        for n in range(1,9):
            count=q**(2*n)+(4+2*eps**n)*q**n+1
            eigen_sum=1+q**(2*n)+4*q**n+2*(eps*q)**n
            assert sp.expand(count-eigen_sum)==0
    return dict(parametrization=True,blowup_chart=True,zeta_power_sums=True)


def main():
    data=dict(status='PASS',field_rows=field_geometry(),fan=fan_checks(),
              symbolic=symbolic_checks())
    Path(__file__).with_name('resolution_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS:',len(data['field_rows']),
          'field geometry cases; smooth eight-ray fan; Picard sign decomposition;')
    print('parametrization, blowup chart, rational-map evaluations, zeta coefficients')


if __name__=='__main__':
    main()
