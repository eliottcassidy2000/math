import sympy as S,json
h,s,e,L,R,d,C,b,T=S.symbols('h s e L R d C b T');v=S.Matrix([h,s,e,L,R,d,C,b,T])
expr={
'O':[h,s,0,L,2*R,d,2*C,2*b,3*T],
'Z':[h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b],
'J':[h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C]}
m={k:S.Matrix(a).jacobian(v) for k,a in expr.items()}
l=S.Matrix([[1,0,0,1,0,0,0,0,0]])
words=['Z','J','OZ','OJ','ZZ','ZJ','JJ','JOZ'];rows=[]
for w in words:
 a=S.eye(9)
 for c in w:a=m[c]*a
 rows.append(list(l*a-l))
a=S.Matrix(rows);det=a[:,[0,2,3,4,5,6,7,8]].det()
unit=S.eye(9)[:,1]
if det!=32 or a.nullspace()!=[unit]:raise RuntimeError('observability certificate failed')
for x in m.values():
 if x*unit!=unit:raise RuntimeError('start-only difference not invariant')
print(json.dumps({'status':'PASS','words':words,'rows':rows,'eight_coordinate_determinant':int(det),'scope':'All-word synchronization criterion only; physical pair existence and spectrum occupancy remain open.'},default=int,indent=2))
