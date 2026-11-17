"""Exact symbolic controls, not a replacement for the all-order proofs."""
import json
from pathlib import Path
import sympy as s
x,y,z,w=s.symbols('x y z w')
A=s.Matrix([[1,-2,2],[2,-1,2],[2,-2,3]])
v=s.Matrix([x,y,z]); u=A*v
assert A.det()==1
assert s.expand(u[0]**2+u[1]**2-u[2]**2)==x*x+y*y-z*z
assert A*s.Matrix([3,4,5])==s.Matrix([5,12,13])
assert A*s.Matrix([1,2,3])==s.Matrix([3,6,7])
assert 1+4+4==9 and 9+36+4==49
assert 3**3+4**3+5**3==6**3
for d in range(3,9):
    f=x**d+y**d-z**d
    h=s.det(s.hessian(f,(x,y,z)))
    assert s.expand(h+(d*(d-1))**3*(x*y*z)**(d-2))==0
    assert 3**d+4**d-5**d!=0
    g=f.subs(dict(zip((x,y,z),u)), simultaneous=True)
    assert s.expand(s.det(s.hessian(g,(x,y,z)))-A.det()**2*h.subs(dict(zip((x,y,z),u)), simultaneous=True))==0
F=s.Matrix([x+x**3,y]); J=F.jacobian([x,y]).det()
assert J==1+3*x*x
assert s.expand(J.subs(x,s.I/s.sqrt(3)))==0
for d in range(1,12):
    G=s.Matrix([x+y**d,y])
    assert G.jacobian([x,y]).det()==1
    assert s.expand(G[0].subs(x,x-y**d))==x
receipt={'status':'PASS','quadratic_matrix_controls':4,'higher_degree_hessian_controls':6,'arbitrary_degree_triangular_controls':11,'scope':'Exact controls support the symbolic proofs; no conjecture completion.'}
Path(__file__).with_name('RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
