from pathlib import Path
import sympy as s,json
x=s.symbols('t');A=s.Matrix([[2,1],[1,1]]);B=s.Matrix([[3,-1],[1,0]]);U=s.Matrix([[1,-1],[0,1]])
assert A*U==U*B and U.det()==1
rows=[]
for n in range(1,13):
 L=s.zeros(n)
 for i in range(n-1):L[i,i]+=1;L[i+1,i+1]+=1;L[i,i+1]-=1;L[i+1,i]-=1
 f=s.expand((s.eye(n)+x*L).det())
 assert f.subs(x,1)==s.fibonacci(2*n)
 if n>=3:assert s.expand(f-(1+2*x)*prev+x*x*prevprev)==0
 if n>=2:prevprev=prev
 prev=f
 rows.append({'n':n,'polynomial':str(f),'forest_count':int(f.subs(x,1)),'trace':int(s.trace(A**n))})
assert int(s.trace(A**2))==7 and int(s.fibonacci(8))==21
assert int(s.fibonacci(6)*s.fibonacci(8)*s.fibonacci(10))==9240
out={'status':'PASS exact symbolic controls','scope':'Path orders1..12, integral conjugacy, distinct trace and boundary observables','rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
