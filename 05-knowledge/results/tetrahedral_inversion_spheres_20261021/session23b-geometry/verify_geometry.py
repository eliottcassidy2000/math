"""Exact, reproducible checks for the seven-label tetrahedron sphere audit."""
from pathlib import Path
import random
import sympy as S

V = S.Matrix
Z = V([0, 0, 0])

def simp(x):
    return S.simplify(x)

def construction(vertices):
    v = vertices[0]
    face = vertices[1:]
    s = sum(face, Z)
    t = v + s
    g = t / 4
    n = (face[1] - face[0]).cross(face[2] - face[0])
    height = n.dot(face[0])
    if n.dot(n) == 0 or height == 0 or g.dot(g) == 0:
        raise ValueError('An excluded definition degeneracy occurred.')
    w = height * n / n.dot(n)
    q = t - 4 * v
    points = {
        'A': simp(v - 2 * v.dot(q) * q / q.dot(q)),
        'B': simp(w / w.dot(w)),
        'C': v,
        'H': simp(s - 2 * w),
        'E': simp(g / g.dot(g)),
        'D': t,
        'F': 4 * w,
    }
    assert simp(points['A'].dot(points['A']) - 1) == 0
    for a, b, c in [(face[0], face[1], face[2]),
                    (face[1], face[0], face[2]),
                    (face[2], face[0], face[1])]:
        assert simp((points['H'] - a).dot(b - c)) == 0
    for a in face:
        assert simp(points['B'].dot(a) - 1) == 0
    return s, t, w, points

counterexample = [V(x) for x in [(0,0,-1), (1,0,0), (0,1,0), (0,0,1)]]
s,t,w,points = construction(counterexample)
assert S.det(S.Matrix.hstack(*(x-counterexample[0] for x in counterexample[1:]))) != 0
print('Rational counterexample:')
for label, point in points.items():
    print(label, list(point))
assert (points['B']-points['H']).cross(points['F']-points['H']) == Z
assert len({tuple(points[k]) for k in ['B','H','F']}) == 3
print('B, H, F are three distinct collinear points: VERIFIED')
det = S.det(V([[p[0], p[2], p.dot(p), 1] for p in [points[k] for k in ['A','C','D','E']]]))
assert det == -S.Rational(56,9)
print('Coplanar circle determinant for A,C,D,E:', det)

T = S.symbols('T', nonzero=True)
c = (T-4)/3
residual = 16/T - 4*(4*(1+c))/T + c
assert simp(residual - (T-4)*(T-16)/(3*T)) == 0
print('Four-point residual factorization:', S.factor(residual))

p = S.symbols('p')
center_norm_numerator = 1 + 48*p*p + 4 + 8*p + 24*p
assert S.factor(center_norm_numerator) == (4*p+1)*(12*p+5)
print('Exceptional center norm numerator:', S.factor(center_norm_numerator))

a = S.symbols('a', real=True, nonzero=True)
root3 = S.sqrt(3)
v = V([-root3*(1-a*a)/(2*(1+a*a)), 2*a/(1+a*a), (1-a*a)/(2*(1+a*a))])
vertices = [v, V([root3/2, 0, S.Rational(1,2)]),
            V([0, root3/2, S.Rational(1,2)]),
            V([0, -root3/2, S.Rational(1,2)])]
s,t,w,points = construction(vertices)
center = V([0,a/2,1])
radius2 = 1+a*a/4
assert simp(t.dot(t)-4) == 0
assert simp(w.dot(w)-S.Rational(1,4)) == 0
assert simp(s.dot(s)-3) == 0
for label, point in points.items():
    assert simp((point-center).dot(point-center)-radius2) == 0, label
p = simp(v.dot(w))
formula_center = simp((v+4*p*s+4*w)/(2*(1+4*p)))
assert simp(formula_center-center) == Z
assert simp((5+12*p)/(4*(1+4*p))-radius2) == 0
assert simp(points['D']-points['E']) == Z
assert simp(points['B']-points['F']) == Z
print('One-parameter exceptional family: all seven exact residuals vanish.')
print('Center = (0,a/2,1), radius squared = 1+a^2/4.')

generic_vertices = [V([0,0,-1]), V([1,0,0]), V([0,1,0]), V([S.Rational(3,5),0,S.Rational(4,5)])]
s,t,w,points = construction(generic_vertices)
v=generic_vertices[0];T=t.dot(t);k=w.dot(w);delta=v.dot(s.cross(w))
u=simp((5*s.cross(w)+(T-1)*w.cross(v)+(1+4*k)*v.cross(s))/(2*delta))
assert u == V([-S.Rational(43,24),S.Rational(37,6),-S.Rational(5,2)])
assert simp(u.dot(u)-4) == S.Rational(25049,576)
for label in ['B','C','D','E','F']:
    x=points[label]
    assert simp(x.dot(x)-2*u.dot(x)+4) == 0
x=4*v
assert simp(x.dot(x)-2*u.dot(x)+4) == 0
assert simp(points['A'].dot(points['A'])-2*u.dot(points['A'])+4) == S.Rational(1178,225)
assert simp(points['H'].dot(points['H'])-2*u.dot(points['H'])+4) == S.Rational(38,5)
print('Generic six-point inversion sphere: exact zero residuals at B,C,D,E,F,4C.')
print('A and H are strictly outside: residuals 1178/225 and 38/5.')

rng = random.Random(2302)
def rational_sphere_point():
    u, v = S.Rational(rng.randint(-5,5),rng.randint(1,5)), S.Rational(rng.randint(-5,5),rng.randint(1,5))
    d=1+u*u+v*v
    return V([2*u/d, 2*v/d, (1-u*u-v*v)/d])

accepted = 0
while accepted < 100:
    vertices = [rational_sphere_point() for _ in range(4)]
    if S.det(S.Matrix.hstack(*(x-vertices[0] for x in vertices[1:]))) == 0:
        continue
    try:
        s,t,w,points = construction(vertices)
    except ValueError:
        continue
    if t.dot(t) == 4:
        continue
    rows=[];rhs=[]
    for label in ['A','C','D','E']:
        point=points[label]
        rows.append([2*point[0],2*point[1],2*point[2],-1])
        rhs.append(point.dot(point))
    M=V(rows)
    assert M.rank() < M.row_join(V(rhs)).rank()
    v=vertices[0];delta=v.dot(s.cross(w));T=t.dot(t);k=w.dot(w)
    if delta != 0:
        u=simp((5*s.cross(w)+(T-1)*w.cross(v)+(1+4*k)*v.cross(s))/(2*delta))
        for label in ['B','C','D','E','F']:
            x=points[label]
            assert simp(x.dot(x)-2*u.dot(x)+4) == 0
        A,H=points['A'],points['H']
        mu=-2*v.dot(t-4*v)/(t-4*v).dot(t-4*v)
        assert simp(A.dot(A)-2*u.dot(A)+4-mu*(16-T)) == 0
        assert simp(H.dot(H)-2*u.dot(H)+4-(6-2*v.dot(s))) == 0
    accepted += 1
print('100 nondegenerate rational unit-sphere tetrahedra: exact inconsistent sphere systems verified.')
print('All checks passed.')
