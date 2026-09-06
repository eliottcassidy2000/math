"""Exact degree-seven hostile and degree-eight geometry-decoder controls.

Standard library only. Numerical exploration located the rational witness;
no numerical quantity is used below. All gates remain active under python -O.
Writes an LF JSON certificate beside this source, or into ../05-knowledge/results
when filed in 04-computation. The universal theorem is proved in the report.
"""
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys

sys.stdout.reconfigure(newline='\n')

GATES = 0


def need(ok, label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def trim(a):
    a = list(map(F, a))
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0) for i in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c * t for t in a])


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def val(a, x):
    out = F(0)
    for c in reversed(a):
        out = out * x + c
    return out


def divmod_poly(a, b):
    a, b = trim(a), trim(b)
    if b == [0]:
        raise ZeroDivisionError()
    q = [F(0)] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        j, c = len(a) - len(b), a[-1] / b[-1]
        q[j] = c
        a = add(a, [F(0)] * j + scale(b, -c))
    return trim(q), a


def gcd_poly(a, b):
    a,b=trim(a),trim(b)
    while trim(b) != [0]:
        a, b = b, divmod_poly(a, b)[1]
    return scale(a, 1 / a[-1])


def sturm(a):
    a = trim(a)
    b = [i * a[i] for i in range(1, len(a))]
    out = [a, b]
    while out[-1] != [0]:
        r = scale(divmod_poly(out[-2], out[-1])[1], -1)
        if r == [0]:
            break
        # Positive scaling controls fraction growth without changing signs.
        r = scale(r, 1 / abs(r[-1]))
        out.append(r)
    return out


def variation(seq, endpoint):
    signs = []
    for a in seq:
        if endpoint == '+inf':
            t = a[-1]
        elif endpoint == '-inf':
            t = a[-1] * (-1) ** (len(a) - 1)
        else:
            t = val(a, endpoint)
        if t:
            signs.append(1 if t > 0 else -1)
    return sum(a != b for a, b in zip(signs, signs[1:]))


def count(a, left, right):
    seq = sturm(a)
    return variation(seq, left) - variation(seq, right)


def det(a):
    a = [list(map(F, row)) for row in a]
    out = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            out = -out
        t = a[i][i]
        out *= t
        for j in range(i + 1, len(a)):
            r = a[j][i] / t
            for k in range(i + 1, len(a)):
                a[j][k] -= r * a[i][k]
    return out


def solve(a, b):
    a = [list(map(F, row)) + [F(c)] for row, c in zip(a, b)]
    for i in range(len(a)):
        j = next(j for j in range(i, len(a)) if a[j][i])
        a[i], a[j] = a[j], a[i]
        t = a[i][i]
        a[i] = [v / t for v in a[i]]
        for j in range(len(a)):
            if j != i:
                t = a[j][i]
                a[j] = [u - t * v for u, v in zip(a[j], a[i])]
    return [row[-1] for row in a]


def lead(a):
    return [det([row[:k] for row in a[:k]]) for k in range(1, len(a) + 1)]


def polynomials(x, y, z):
    return ([-z, y, -x, 55, -13, 1],
            [3 * y / 7, -2 * x / 3, 45, -12, 1],
            [y / 7, -5 * x / 12, 36, -11, 1])


def quotient_moments(numerator, denominator, order):
    den, num = list(reversed(trim(denominator))), list(reversed(trim(numerator)))
    out = []
    for j in range(order + 1):
        out.append(((num[j] if j < len(num) else 0) -
                    sum(den[k] * out[j-k] for k in range(1, min(j, len(den)-1)+1))) / den[0])
    return out


def hankel(m, n, shift=0):
    return [[m[i+j+shift] for j in range(n)] for i in range(n)]


def response(x, y, z):
    # Literal original convolutions; exponents start at -1.
    from math import comb
    beta = {-1: F(1), 0: F(13), 1: F(55), 2: x, 3: y, 4: z}
    cc = {-1: F(1), 0: F(12), 1: F(45), 2: 2*x/3, 3: 3*y/7}
    dd = {-1: F(1), 0: F(11), 1: F(36), 2: 5*x/12, 3: y/7}
    def conv(a,b):
        out = {}
        for i, u in a.items():
            for j, v in b.items():
                out[i+j] = out.get(i+j,F(0)) + u*v
        return out
    odd = {j: F(comb(14,2*j+1)) for j in range(7)}
    even = {j: F(comb(14,2*j)) for j in range(8)}
    carrier = conv(odd,odd)
    for j, c in conv(even,even).items():
        carrier[j-1] = carrier.get(j-1,F(0)) + c
    raw = conv(beta,beta)
    for j,c in conv(cc,dd).items():
        raw[j+1] = raw.get(j+1,F(0)) + 2*c
    qq = {j:c*carrier.get(j,0) for j,c in raw.items()}
    pp = {j:c*odd.get(j,0) for j,c in beta.items()}
    return pp,qq


def qvalue(p, t):
    return sum(c*t**j for j,c in p.items())


# Universal common-zero algebra, in the variable a=Re(r).
h = [F(63,2), -21, 3]
imag_y = add([0,126,-84,12], mul([28,-12],h))
need(imag_y == scale(mul([-7,2],[21,-14,2]),-6), 'imaginary elimination')
need(val(h,F(7,2)) == F(-21,4), 'middle root excludes nonzero imaginary part')
need(h == scale([21,-14,2],F(3,2)), 'other two roots force zero imaginary part')
# The same x(v), y(v) annihilate C,D identically.
xv, yv = [0,108,-36,F(24,7)], [0,0,63,-28,3]
need(add(add([0,0,45,-12,1],scale(mul([0,1],xv),F(-2,3))),scale(yv,F(3,7))) == [0], 'C common-zero parameter')
need(add(add([0,0,36,-11,1],scale(mul([0,1],xv),F(-5,12))),scale(yv,F(1,7))) == [0], 'D common-zero parameter')

x,y,s = F(78071,1000),F(601,50),F(57,2)
z = 12*y/(7*s)-x/s**2+10/s**3-1/(11*s**4)
M0 = F(707,100)
B,C,D = polynomials(x,y,z)
pp,qq = response(x,y,z)
need(qvalue(pp,-s) == 0, 'original root retained')
qpositive = qvalue(qq,-s)
need(qpositive > 0, 'positive complete response')
need(qq[-1] == 28 and set(j for j,c in qq.items() if c) == set(range(-1,9)), 'all original carries retained')
need(count(B,0,'+inf') == 3 and count(B,'-inf',0) == 0, 'beta has exactly three real roots')
need(len(gcd_poly(B,[i*B[i] for i in range(1,len(B))])) == 1, 'beta simple')
P = [pp.get(i,F(0))*(-1)**i for i in range(5)]
need(count(P,0,'+inf') == 4 and len(gcd_poly(P,[i*P[i] for i in range(1,5)])) == 1, 'all original phases positive simple')
need(5*M0*M0-26*M0-67 < 0 and M0 > F(13,5), 'chosen support is below exact anchor endpoint')
need(y > F(161875,888583) and x > 75 and y > F(8750,8241)*(x-75), 'proved floor and slope survive')
need(5*z < M0*y, 'coefficient support inequality survives')
es = [F(1),F(13),F(55),x,y,z]
from math import comb
newton = [(es[k]/comb(5,k))**2-(es[k-1]/comb(5,k-1))*(es[k+1]/comb(5,k+1)) for k in range(1,5)]
for t in newton:
    need(t > 0, 'strict Newton margin')

record = {'witness': {'x':str(x),'y':str(y),'z':str(z),'s':str(s),'support_cap':str(M0)},
          'P_minus_s': [str(t) for t in P], 'Q_at_minus_s':str(qpositive),
          'Q_coefficients':{str(j):str(c) for j,c in sorted(qq.items())},
          'newton_margins':list(map(str,newton))}
quadratures=[]
for label,A in [('C',C),('D',D)]:
    m = quotient_moments(A,B,9)
    need(m[:3] == ([F(1),F(1),F(3)] if label=='C' else [F(1),F(2),F(7)]), label+' anchor moments')
    H,K = hankel(m,4),hankel(m,4,1)
    U = [[M0*H[i][j]-K[i][j] for j in range(4)] for i in range(4)]
    minors = {'ordinary4':lead(H),'shifted4':lead(K),'upper4':lead(U)}
    for name, vv in minors.items():
        for d in vv:
            need(d > 0,label+' '+name+' positive definite')
    V = [[M0*m[i+j+1]-m[i+j+2] for j in range(3)] for i in range(3)]
    minors['quadratic_localizer3']=lead(V)
    for d in minors['quadratic_localizer3']:
        need(d > 0,label+' lower-degree quadratic localizer')
    bad = det(hankel(m,5))
    need(bad < 0,label+' eighth moment rejects surrogate')
    quartic = solve(H,[-m[j+4] for j in range(4)])+[F(1)]
    num_desc=[]
    dd=list(reversed(quartic))
    for j in range(4):
        num_desc.append(sum(dd[k]*m[j-k] for k in range(j+1)))
    numerator=list(reversed(num_desc))
    qm=quotient_moments(numerator,quartic,8)
    for j in range(8):
        need(qm[j] == m[j],label+' positive quadrature reproduces moment '+str(j))
    need(qm[8] != m[8],label+' quadrature loses native eighth moment')
    need(count(quartic,0,M0) == 4,label+' quadrature has four interior nodes')
    need(len(gcd_poly(quartic,[j*quartic[j] for j in range(1,5)])) == 1,label+' quadrature simple')
    quadratures.append(quartic)
    record[label]={'moments':list(map(str,m)),
                   'leading_minors':{name:list(map(str,vv)) for name,vv in minors.items()},
                   'ordinary5_determinant':str(bad),
                   'quadrature_denominator':list(map(str,quartic)),
                   'quadrature_numerator':list(map(str,numerator)),
                   'quadrature_eighth_moment':str(qm[8])}
need(len(gcd_poly(*quadratures))==1,'the two canonical quadratures have disjoint nodes')

# Independent algebraic controls: literal beta roots and a shared canceled node.
controls=[]
for name,xx,yy,zz in [('strict',F(2127,25),F(27144,625),F(3564,625)),
                     ('common_node',F(648,7),F(54),F(27,7))]:
    BB,CC,DD=polynomials(xx,yy,zz)
    need(count(BB,0,'+inf')==5,name+' all five beta nodes positive')
    data={}
    for lab,AA in [('C',CC),('D',DD)]:
        mm=quotient_moments(AA,BB,9)
        ll=lead(hankel(mm,5))
        for d in ll[:4]:need(d>0,name+' '+lab+' first four pivots')
        need(ll[4]>0 if name=='strict' else ll[4]==0,name+' '+lab+' exact rank')
        data[lab]=list(map(str,ll))
    gg=gcd_poly(gcd_poly(BB,CC),DD)
    need(gg==([F(1)] if name=='strict' else [F(-3),F(1)]),name+' common denominator factor')
    if name=='strict':
        literal=[F(1)]
        roots=[F(1,5),F(3,5),F(9,5),F(22,5),F(6)]
        for root in roots:literal=mul(literal,[-root,1])
        need(literal==BB,'literal strict anchor roots')
        derivative=[j*BB[j] for j in range(1,6)]
        for root in roots:
            for lab,AA in [('C',CC),('D',DD)]:
                need(val(AA,root)*val(derivative,root)>0,lab+' literal positive residue')
    else:
        reduced=divmod_poly(BB,[-3,1])[0]
        bb_intervals=[(F(8,100),F(9,100)),(F(104,100),F(105,100)),
                      (F(225,100),F(226,100)),(F(662,100),F(663,100))]
        cc_intervals=[(F(5,10),F(7,10)),(F(2),F(21,10)),(F(63,10),F(64,10))]
        dd_intervals=[(F(2,10),F(3,10)),(F(16,10),F(17,10)),(F(61,10),F(62,10))]
        for lo,hi in bb_intervals:need(count(reduced,lo,hi)==1,'common control beta interval')
        for AA,intervals in [(CC,cc_intervals),(DD,dd_intervals)]:
            reduced_a=divmod_poly(AA,[-3,1])[0]
            for i,(lo,hi) in enumerate(intervals):
                need(count(reduced_a,lo,hi)==1,'common control numerator interval')
                need(bb_intervals[i][1]<lo<hi<bb_intervals[i+1][0], 'literal weak interlacing after cancellation')
    controls.append({'name':name,'x':str(xx),'y':str(yy),'z':str(zz),'leading_minors':data,'gcd':list(map(str,gg))})

# A repeated-root C-only boundary must not be silently accepted as two-interlacer.
BB,CC,DD=polynomials(F(75),F(0),F(0))
mm=quotient_moments(CC,BB,8)
need(mm == [F(1)]+[F(3)**(j-1) for j in range(1,9)],'weak repeated C measure')
need(det(hankel(quotient_moments(DD,BB,4),3)) == F(-37,16),'D excludes C-only repeated boundary')
record['controls']=controls

stem=Path(__file__).stem
outdir=Path(__file__).resolve().parent
if outdir.name=='04-computation':outdir=outdir.parent/'05-knowledge'/'results'
outdir.mkdir(parents=True,exist_ok=True)
certificate=outdir/(stem+'_certificate.json')
certificate.write_bytes((json.dumps(record,indent=2,sort_keys=True)+'\n').encode())
print('Exact original-phase hostile:', {k:str(t) for k,t in zip(('x','y','z','s'),(x,y,z,s))})
print('Full Q(-s) =',qpositive)
print('Support cap =',M0,'; below (13+6 sqrt(14))/5')
for label in ('C','D'):
    print(label,'leading minors =',record[label]['leading_minors'])
    print(label,'ordinary 5x5 determinant =',record[label]['ordinary5_determinant'])
print('Beta real-root count = 3; original positive phase count = 4; both simple.')
print('Both exact four-atom quadratures match moments 0..7 and fail the native moment 8.')
print('Strict and common-node weak controls PASS; C-only repeated boundary rejected by D.')
print('Certificate SHA256',hashlib.sha256(certificate.read_bytes()).hexdigest())
print('PASS',GATES,'always-active exact gates')
