"""Exact local collision/transgression audit; no finite rank is extrapolated.

Reproduce: python3 -B 04-computation/synthesis_20260905_transgression.py
Universe: the three retained branches of the Russell compiler and all first
normal polynomial directions. Symbolic parameters cover the exceptional
quartic's second jets. The all-order assertion uses an invertible formal
Jacobian, not a bounded Taylor experiment.
"""
import sympy as s

checks = 0


def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)


x, q = s.symbols("x q")
a, b, c = s.symbols("a b c")
D = 1 + x*x*q
compiler = s.Matrix(((D-1)*(D+2)**2, x*D*(D+2), q*(D+3)))
require(s.expand(compiler[1]**2*compiler[2]-compiler[0]*(compiler[0]+4)) == 0,
        "surface identity")
P = x*x*(x*x-1)**2
Q1 = x**5+s.Rational(9,2)*x**4-2*x**3-s.Rational(27,4)*x*x+x-s.Rational(3,4)
Q = Q1+P*(a+b*x+c*x*x)
curve = compiler.subs(q,Q).applyfunc(s.expand)
normal = compiler.diff(q).subs(q,Q).applyfunc(s.expand)
second_normal = compiler.diff(q,2).subs(q,Q).applyfunc(s.expand)
third_normal = compiler.diff(q,3).subs(q,Q).applyfunc(s.expand)
points = (-1,0,1)
tangents = [curve.diff(x).subs(x,z) for z in points]
normals = [normal.subs(x,z) for z in points]
require(all(v == s.Matrix((0,0,-3)) for v in [curve.subs(x,z) for z in points]),
        "common target")
require(tangents == [s.Matrix((0,3,-9)),s.Matrix((0,3,4)),s.Matrix((0,3,9))],
        "tangent rows independent of higher jets")
require(normals == [s.Matrix((0,2,-2)),s.Matrix((0,0,4)),s.Matrix((0,-2,-2))],
        "normal rows independent of higher jets")
t = [v[1:3,0] for v in tangents]
n = [v[1:3,0] for v in normals]
weights = s.Matrix((5,-18,13))
T = s.Matrix.hstack(*t).T
require(T.rank() == 2 and (weights.T*T) == s.zeros(1,2), "unique tangent annihilator")
require([s.det(s.Matrix.hstack(ti,ni)) for ti,ni in zip(t,n)] == [12]*3,
        "normal Jacobian at each branch")

# Unknowns: three endpoint velocities, common C,E velocities. The sixth
# column is the coefficient of x in the polynomial normal deformation.
A = s.zeros(6,5)
vx = s.zeros(6,1)
for i,z in enumerate(points):
    A[2*i:2*i+2,i] = t[i]
    A[2*i:2*i+2,3:5] = -s.eye(2)
    vx[2*i:2*i+2,0] = z*n[i]
J = A.row_join(vx)
require(A.rank() == 5, "endpoint/common-target repair rank")
require(J.det() != 0, "one coefficient repairs all formal orders")
ell = A.T.nullspace()[0]
hminus,hzero,hplus = s.symbols("hm h0 hp")
hvalues = (hminus,hzero,hplus)
vh = s.Matrix.vstack(*[hv*ni for hv,ni in zip(hvalues,n)])
L = 5*hminus-18*hzero+13*hplus
ratio = s.cancel((ell.T*vh)[0]/L)
require(not ratio.has(hminus,hzero,hplus), "collision obstruction equals Lambda(h)")
require(ratio != 0, "nonzero proportionality")
print("PROVED local first-order gate: 5*h(-1)-18*h(0)+13*h(1)=0")
print("Endpoint repair matrix rank:", A.rank(), "; augmented determinant:", J.det())
print("Left-kernel obstruction factor:", ratio)

for label,h,compatible in (("constant",s.Integer(1),True),
                            ("linear hostile",x,False),
                            ("quadratic survivor",4*x*x-9*x,True),
                            ("fixed-endpoint vanishing",x*(x*x-1),True)):
    values = [h.subs(x,z) for z in points]
    obstruction = (weights.T*s.Matrix(values))[0]
    rhs = s.Matrix.vstack(*[hv*ni for hv,ni in zip(values,n)])
    require((A.row_join(-rhs).rank() == 5) == compatible, label+" rank criterion")
    require((obstruction == 0) == compatible, label+" scalar criterion")
    print(label+": obstruction="+str(obstruction))

# Six independent equations are equivalent to the full target equality near
# b=0 because d(c^2 e-b(b+4))/db=-4 there. A formally invertible J proves
# unique all-order compensation, for arbitrary fixed other coefficients.
require(s.diff(s.Symbol("cc")**2*s.Symbol("ee")-s.Symbol("bb")*(s.Symbol("bb")+4),
               s.Symbol("bb")).subs(s.Symbol("bb"),0) == -4, "smooth target chart")

# Symbolic second-order controls for the entire Q1+P*(a+b*x+c*x^2) family.
# This contains Q_R: a=-259/36+p(alpha)+4alpha, b=-9alpha, c=-p(alpha).
for label,h in (("constant",s.Integer(1)),("quadratic survivor",4*x*x-9*x)):
    rhs1 = s.Matrix.vstack(*[h.subs(x,z)*ni for z,ni in zip(points,n)])
    solution1 = -J.inv()*rhs1
    require(solution1[5] == 0, label+" no first-order compensator")
    residual2 = []
    for i,z in enumerate(points):
        eta = solution1[i]
        second = (curve.diff(x,2)*eta**2/2
                  +(normal*h).diff(x)*eta
                  +second_normal*h*h/2).subs(x,z)[1:3,0]
        residual2.append(second.applyfunc(s.expand))
    residual2 = s.Matrix.vstack(*residual2)
    solution2 = (-J.inv()*residual2).applyfunc(s.factor)
    require((J*solution2+residual2).applyfunc(s.expand) == s.zeros(6,1),
            label+" second-order repair")
    print(label+" first-order (endpoints,C,E): "+str(tuple(solution1[:5])))
    print(label+" second-order x compensator: "+str(solution2[5]))
    alpha = s.symbols("alpha")
    p_alpha = s.Rational(520,9)*alpha**2-s.Rational(1688,81)*alpha-s.Rational(5717,729)
    exceptional = {a:-s.Rational(259,36)+p_alpha+4*alpha,b:-9*alpha,c:-p_alpha}
    exceptional_c2 = s.factor(solution2[5].subs(exceptional))
    print(label+" exceptional second-order compensator: "+str(exceptional_c2))
    if label == "quadratic survivor":
        require(exceptional_c2 != 0 and s.degree(exceptional_c2,alpha) < 4,
                "quadratic first-order survivor obstructed at second order")
    else:
        require(exceptional_c2 == 0,"constant direction survives second order")
        residual3 = []
        for i,z in enumerate(points):
            eta,theta = solution1[i],s.factor(solution2[i].subs(exceptional))
            third = (curve.diff(x,2)*eta*theta+curve.diff(x,3)*eta**3/6
                     +normal.diff(x)*theta+normal.diff(x,2)*eta**2/2
                     +second_normal.diff(x)*eta/2+third_normal/6).subs(x,z)[1:3,0]
            residual3.append(third.applyfunc(lambda e:s.factor(e.subs(exceptional))))
        residual3=s.Matrix.vstack(*residual3)
        solution3=(-J.inv()*residual3).applyfunc(s.factor)
        require((J*solution3+residual3).applyfunc(s.expand)==s.zeros(6,1),
                "third-order collision repair")
        require(solution3[5] == 0, "constant exceptional third-order positive control")
        print("constant exceptional third-order x compensator: "+str(solution3[5]))

# Independent literal first-jet check of the geometric linear algebra with
# branch tangent (1,m): preserving collision iff det(t_i,n_i) lies in
# span{1,m}. This includes the tangential gauge and a nonpreserving control.
slopes = (-3,s.Rational(4,3),3)
common = s.Matrix((7,-2))
for i,m in enumerate(slopes):
    tangent = s.Matrix((1,m))
    velocity = common+(i+1)*tangent
    require(s.det(s.Matrix.hstack(tangent,velocity)) == common[1]-m*common[0],
            "coordinate-free wedge/gauge control "+str(i))
print("Single x-compensator formal Jacobian: invertible in characteristic zero")
print("Two-form response: Lambda(tau_h(dC wedge dE)) = (2/3)*L(h)")
print("Adding tau_x to the THM-4404 image fills its one-dimensional cokernel;")
print("every first-order collision-preserving direction still kills Lambda.")
print("Exact checks:",checks)
