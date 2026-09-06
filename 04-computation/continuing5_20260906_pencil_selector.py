"""Exact quartic endpoint selectors for the anchored residue pencils.

Universal claims have separate proofs. This producer retains three strict
fibre witnesses, two nonnested-channel hostiles and a positive lower endpoint.
"""
from pathlib import Path
import sympy as S
import json,sys
sys.stdout.reconfigure(newline='\n')
v,z=S.symbols('v z'); x,y=S.symbols('x y',real=True)
gates=0
def check(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)

def moments(xx,yy,label):
    if label=='C':m=[S.Integer(1),S.Integer(1),S.Integer(3),xx/3-16,16*xx/3-373-4*yy/7]
    else:m=[S.Integer(1),S.Integer(2),S.Integer(7),7*xx/12-19,115*xx/12-632-6*yy/7]
    for j in range(5,9):
        m.append(S.expand(13*m[j-1]-55*m[j-2]+xx*m[j-3]-yy*m[j-4]+z*m[j-5]))
    return m

def source(xx,yy,label):
    F=v**5-13*v**4+55*v**3-xx*v**2+yy*v
    A=(v**4-12*v**3+45*v**2-2*xx*v/3+3*yy/7 if label=='C' else
       v**4-11*v**3+36*v**2-5*xx*v/12+yy/7)
    return F,A

def hankel(m):return S.Matrix(5,5,lambda i,j:m[i+j])
def minors(H):return [S.factor(H[:j,:j].det()) for j in range(1,6)]
def polynomial(p):return [str(t) for t in S.Poly(p,z).all_coeffs()]

for label,slopes in [('C',[1,14,130,904+4*x/3]),('D',[1,15,147,1067+19*x/12])]:
    m=moments(x,y,label);H=hankel(m);T=H.diff(z)
    for j in range(9):check(S.diff(m[j],z,2)==0,'affine moments')
    check([S.diff(t,z) for t in m[5:]]==slopes,'exact slope coefficients')
    check(T[:1,:]==S.zeros(1,5) and T[:,:1]==S.zeros(5,1),'fixed null coordinate')
    check(S.expand(T[1:,1:].det())==1,'rank-four slope determinant')
    check(T.subs({x:0})[1:,1:].det()==1,'connected nonsingular inertia homotopy')

fixtures=[('C_upper',S.Rational(155,2),S.Integer(9),S.Rational(1,200)),
          ('D_upper',S.Rational(155,2),S.Rational(37,4),S.Rational(1,5)),
          ('positive_lower',S.Rational(311,4),S.Rational(21,2),S.Rational(1,5))]
record={'status':'FINITE-EXACT controls for proved quartic selector', 'fixtures':[]}
for name,xx,yy,zz in fixtures:
    row={'name':name,'x':str(xx),'y':str(yy),'strict_z':str(zz),'channels':{}}
    ranges={}
    for label in ('C','D'):
        F,A=source(xx,yy,label);H=hankel(moments(xx,yy,label))
        R=S.Poly(S.resultant(F-z,A,v),z,domain=S.QQ)
        check(S.expand(H.det(method='domain-ge')-R.as_expr())==0,'Hankel resultant identity')
        check(R.degree()==4 and R.LC()==1,'monic quartic')
        ds=minors(H.subs(z,zz))
        for d in ds:check(d>0,'strict reference point Sylvester minor')
        intervals=R.intervals(eps=S.Rational(1,100000))
        check(len(intervals)==4 and all(k==1 for ab,k in intervals),'four simple real determinant roots')
        boxes=[ab for ab,k in intervals]
        check(boxes[1][1]<zz<boxes[2][0],'strict witness between middle roots')
        ranges[label]=boxes
        check(S.Poly(A,v).count_roots(0,S.oo)==4,'positive interlacer roots')
        bpoly=S.Poly((F-z).subs(z,zz),v)
        check(bpoly.count_roots(0,S.oo)==5,'five positive beta roots')
        for at in (boxes[0][0]-1,boxes[1][1]+(boxes[2][0]-boxes[1][1])/2,boxes[3][1]+1):
            md=minors(H.subs(z,at))
            check(all(d>0 for d in md)==(boxes[1][1]<at<boxes[2][0]),'outside determinant-positive hostile versus actual PSD')
        row['channels'][label]={'quartic_descending':polynomial(R.as_expr()),
            'isolating_intervals':[[str(a),str(b)] for a,b in boxes],
            'reference_leading_minors':[str(d) for d in ds]}
        print(name,label,'middle-root boxes',boxes[1:3],flush=True)
    if name=='C_upper':check(ranges['C'][2][1]<ranges['D'][2][0],'C strictly sets upper endpoint')
    if name=='D_upper':check(ranges['D'][2][1]<ranges['C'][2][0],'D strictly sets upper endpoint')
    if name=='positive_lower':
        check(max(ranges['C'][1][0],ranges['D'][1][0])>S.Rational(7,100),'strict positive lower domain endpoint')
        check(min(ranges['C'][2][0],ranges['D'][2][0])>S.Rational(1,3),'nonvacuous domain width')
        for label in ('C','D'):
            H=hankel(moments(xx,yy,label)).subs(z,0)
            row['channels'][label]['z0_leading_minors']=[str(d) for d in minors(H)]
        check(any(not all(d>0 for d in minors(hankel(moments(xx,yy,label)).subs(z,0)))
                  for label in ('C','D')),'zero is outside full fibre')
    record['fixtures'].append(row)

record['channel_hostiles']=[]
for label,xx,yy,zz in [('D',S.Rational(155,2),S.Integer(9),S.Rational(1,10)),
                      ('C',S.Rational(155,2),S.Rational(37,4),S.Rational(13,50))]:
    other='C' if label=='D' else 'D'
    good=minors(hankel(moments(xx,yy,label)).subs(z,zz))
    bad=minors(hankel(moments(xx,yy,other)).subs(z,zz))
    check(all(d>0 for d in good),'one actual strict channel survives')
    check(bad[-1]<0,'other full Hankel determinant fails')
    record['channel_hostiles'].append({'survives':label,'x':str(xx),'y':str(yy),'z':str(zz),
        'surviving_leading_minors':[str(d) for d in good],
        'failed_leading_minors':[str(d) for d in bad]})
    print('NONNESTED',label,'survives',(xx,yy,zz),'opposite determinant',bad[-1],flush=True)

out=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name(
    'continuing5_20260906_pencil_selector_certificate.json')
out.write_bytes((json.dumps(record,indent=2)+'\n').encode())
print('PASS',gates,'always-active exact gates')
