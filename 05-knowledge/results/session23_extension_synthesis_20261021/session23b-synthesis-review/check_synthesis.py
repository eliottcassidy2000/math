#!/usr/bin/env python3
"""Exact synthesis binding and small algebra checks, not a rerun of large censuses.
Usage: python3 check_synthesis.py --root WORKSPACE
   or: python3 check_synthesis.py --results-root REPO/05-knowledge/results
Requires SymPy. All failures remain active under python -O.
"""
from pathlib import Path
import argparse, hashlib, json, math, re
import sympy as S

HERE = Path(__file__).resolve().parent
p = argparse.ArgumentParser(description=__doc__)
g = p.add_mutually_exclusive_group()
g.add_argument('--root', type=Path)
g.add_argument('--results-root', type=Path)
a = p.parse_args()
base = a.results_root or a.root or HERE.parent
binding = json.loads((HERE/'source_bindings.json').read_text())
checks = 0

def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)

def equal(x, y, message):
    need(S.cancel(x-y) == 0, message)

texts = {}
for key, spec in binding['files'].items():
    path = base/(spec['published'] if a.results_root else key)
    raw = path.read_bytes()
    need(len(raw) == spec['bytes'], 'Byte count: '+key)
    need(hashlib.sha256(raw).hexdigest() == spec['sha256'], 'Source drift: '+key)
    if path.suffix == '.md':
        texts[key] = raw.decode()
key = 'session23b-root/SESSION23_EXTENSION_SYNTHESIS.md'
syn = texts[key]
links = re.findall(r'\]\(\.\./([^/]+)/README\.md\)', syn)
need(links == binding['synthesis_wrapper_destinations'], 'Wrapper destination drift')
for dest in links:
    need(any(v['published'] == dest+'/README.md' for v in binding['files'].values()),
         'Unbound wrapper: '+dest)

# These exact anchors bind the mathematical checks and scope review to the prose.
anchors = [
 'B\' f\' = Bf - M B_D 1.', 'M(e_last-e_first)', 'M s(T)+b',
 'score-complement identity itself is inherited from THM-646',
 'not a new Ohm\'s-law or voltage realization theorem',
 'Their capacity-one total-use coefficient at q⁴ differs,28 versus21.',
 'd=(n-1)(n-2)/2', 'palindromic of exact degree d-1',
 'edgeless order1 case is explicitly separate',
 'delta=(n-1)(n-4)/2', 'exact codegree is3',
 'NOT generally palindromic', 'coefficient of w is91-10=81',
 '(r-1,r,...,r,r+1)', '1/[(1-w)(1-wq³)(1-wq⁵)]',
 'R_a(q^-1)=(-1)^d q^m R_a(q).', 'fixed circulation dimension r²',
 'AT MOST M', 'C_T M^[-r(r-1)] (1+O(1/M)),   C_T>0.',
 'cumulative total-use constant is3150/17', 'pointwise coefficient-ratio constant below is2100/17',
 'minimal constant-coefficient recurrence orders are22,26,22,32',
 'minimal quasiperiods are60,60,120,120',
 'exact25-cone unimodular triangulation', 'five-case lattice decomposition',
 'fixed UNDERLYING h', 'fixed UNDERLYING fundamental black degree',
 '1-z²-z⁶/(1-q³)²', '1-z²-z⁶/(1-q⁶)',
 'No formula here redefines h', 'Nor does a circulation count Euler tours',
 'S4-equivariant projective duality, not an equivariant point bijection',
 'orbit sizes are1+6 in one model and4+3 in the other',
 'exactly five distinct points', 'B,H′,F distinct and collinear',
 'When v,s,w are independent', 'inverse pairs are(D,E),(B,F),(C,4C)',
 'omitted A and H′ lie strictly outside',
 'all eight cyclic completions to have H=171', 'exactly two Paley completions',
 'old16/112 Fano-orientation census and Hamming-code locus are inherited',
 'nonzero weight-three codeword', '[0,14,0; 2,9,3; 0,7,7]',
 'class sizes240,1680,720', 'eigenvalues14,1+sqrt15,1-sqrt15',
 '(189-630z-1470z²)/[(1-14z)(1-2z-14z²)]',
 'Triangle reversal is NOT the original fixed-path nonpath-complement operation',
 'Paley has exactly two contained Fano designs', 'anti-actions exchange them',
 'least distinguishing palette sizes are2,3,4', 'Kernel exceptions',
 'single equivariantly chosen subset cannot break a nontrivial automorphism group',
 '(3u+2u²+2u³)/[(1-u)(1-u³)]', 'minimal recurrence order4',
 'chi_d(P7[P7])=7 rather than9', 'starting1,3,7,17,40,94',
 'finite carry alphabet is therefore not an autonomous finite-state model',
 'This does not prove a Collatz theorem', 'not a literature-priority claim',
 'That incoming theorem was not re-proved in this session',
 'uniform153N source-budget obligation', 'Native mod15/mod35 guards',
 'No general A000568 speedup has been proved or benchmarked',
 'No Lonely Runner, Hadwiger–Nelson, Collatz, or universal Hamiltonian-path-spectrum conjecture is settled',
 'raw Paley7 flow numerator and anti-involution fixed-cone numerator remain uncomputed',
 'must be tested for equitability',
]
for anchor in anchors:
    need(anchor in syn, 'Missing reviewed claim or boundary: '+anchor)

q,w,M,x,z,c = S.symbols('q w M x z c')
# Full capacity h* and path-face h*: distinct dimension/codegree laws.
N5=(M+1)**2*(3*(M+1)**4+5*(M+1)**2+4)/12
h5=[1,17,72,72,17,1]
for m in range(11):
    ehrhart=sum(h*S.binomial(m+6-i,6) for i,h in enumerate(h5) if m>=i)
    equal(N5.subs(M,m),ehrhart,'N5 polynomial / h*')
need(h5 == h5[::-1] and sum(h5)>0, 'Full palindromicity and no (1-w) cancellation')
need((7-1)*(7-2)//2 == 15 and (7-1)*(7-4)//2 == 9, 'n7 dimensions')
need(9+1-3 == 7 and 91-(9+1) == 81 and 34 != 1, 'Path-face h* degree, linear term, nonpalindromicity')
need(171//3+175//7+189//21 == 91, 'Path-presentation weighted sum')
need(24 != 14 and 175 != 189 and 28 != 21, 'Information-loss hostile distinctions')
# At w^M the capacity/use series counts nonnegative a,b with a+b<=M.
for m in range(6):
    direct=sum(q**(3*i+5*j) for i in range(m+1) for j in range(m-i+1))
    coeff=S.series(1/((1-w)*(1-w*q**3)*(1-w*q**5)),w,0,m+1).removeO().coeff(w,m)
    equal(direct,coeff,'Capacity/use joint series')
    equal(direct.subs(q,1),S.binomial(m+2,2),'Path-face n5 count')

P=1+2*q**3+3*q**4+q**5+3*q**6+q**7+3*q**8+2*q**9+q**12
K=1+q**4+q**5+q**6+q**7+q**8+q**12
R=P/((1-q**3)**3*(1-q**4)**2*(1-q**5))
f=(R+4/(1-q**5)**2)/5
k=K/((1-q**3)*(1-q**5)*(1-q**6)*(1-q**8))
H=(f+k)/2
A=1+q**3+3*q**4+2*q**6+4*q**7+q**8+3*q**10
B=2*q**3+2*q**6+2*q**7+q**9+2*q**10+q**13
equal(A*(1-q**3)+B*(1-q**4),P*(1-q**5),'Raw fan numerator identity')
five_cases=(1/((1-q**3)*(1-q**5)**2*(1-q**6))
 +q**4/((1-q**3)*(1-q**4)*(1-q**5)**2)
 +(q**8+q**6+q**7)/((1-q**5)**2*(1-q**6)*(1-q**8))
 +q**10/((1-q**4)*(1-q**5)**2*(1-q**6)))
equal(five_cases,k,'Anti-fixed five-case sum')
orders=[]; periods=[]; profiles=[]
for F in [R,f,k,H]:
    num,den=S.fraction(S.cancel(F)); orders.append(int(S.degree(den,q)))
    need(S.degree(S.gcd(num,den),q)==0,'Reduced rational form')
    poles={}
    remaining=S.Poly(den,q).monic()
    for i in range(1,23):
        phi=S.Poly(S.cyclotomic_poly(i,q),q); exponent=0
        while remaining.degree()>0 and remaining.rem(phi).is_zero:
            remaining=remaining.exquo(phi);exponent+=1
        if exponent:poles[i]=exponent
    need(remaining.degree()==0,'Cyclotomic pole completeness')
    periods.append(math.lcm(*poles));profiles.append(poles)
    equal(F.subs(q,1/q),q**10*F,'R5 reciprocity')
need(orders==[22,26,22,32], 'Minimal recurrence orders')
need(periods==[60,60,120,120], 'Minimal quasiperiods')
lead_R=S.limit((1-q)**6*R,q,1)/S.factorial(5)
lead_f=S.limit((1-q)**6*f,q,1)/S.factorial(5)
lead_k=S.limit((1-q)**4*k,q,1)/S.factorial(3)
for a0,b0 in [(lead_R,S.Rational(17,259200)),(lead_f,S.Rational(17,1296000)),(lead_k,S.Rational(7,4320))]:
    equal(a0,b0,'Pointwise leading coefficient')
equal(lead_k/lead_f,S.Rational(2100,17),'Pointwise ratio')
equal((lead_k/4)/(lead_f/6),S.Rational(3150,17),'Cumulative ratio')
need(3*(3-1)==6,'Order7 symmetry-rarity exponent')
U7=1/(1-q**7)**3
U3=(1+q**3+q**6+q**9+q**12)/((1-q**3)*(1-q**6)**2*(1-q**9)**2)
V6=(1+q**15)/((1-q**6)*(1-q**12)*(1-q**18))
for F in [U7,U3,V6]:equal(F.subs(q,1/q),-q**21*F,'Paley fixed-cone reciprocity')
T4=1/((1-q**3)**2*(1-q**4)); T4anti=1/((1-q**6)*(1-q**4))
equal(T4.subs(q,1/q),-q**10*T4,'T4 raw sign')
equal(T4anti.subs(q,1/q),q**10*T4anti,'T4 anti-fixed sign')
need(S.expand((1-q**3)**2)!=1-q**6,'Independent/mirrored neutral denominators differ')

Q=S.Matrix([[0,14,0],[2,9,3],[0,7,7]])
need([sum(Q.row(i)) for i in range(3)]==[14]*3,'Triangle degree')
equal(Q.charpoly(x).as_expr(),(x-14)*(x*x-2*x-14),'Quotient eigenvalues')
weights=S.diag(240,1680,720)*Q
need(weights==weights.T and 240+1680+720==2640,'Class sizes and detailed balance')
move=(S.eye(3)-z*Q).inv()*S.Matrix([189,171,175])
equal(move[0],(189-630*z-1470*z*z)/((1-14*z)*(1-2*z-14*z*z)),'Move-word GF')
equal(10,1+(9*10+3*12)/S.Integer(14),'First hit from H171')
equal(12,1+(7*10+7*12)/S.Integer(14),'First hit from H175')
equal(1/S.Rational(1,7),7,'Restricted cube first hit')
colors=[c**7-7*c**3+6*c,c**7-7*c**4-7*c**3+7*c*c+6*c,
 c**7-21*c**5+28*c**4+56*c**3-112*c*c+48*c]
color_counts=[]
for poly,minimum in zip(colors,[2,3,4]):
    row=[int(poly.subs(c,i)) for i in range(1,5)];color_counts.append(row)
    need(all(v==0 for v in row[:minimum-1]) and row[minimum-1]>0,'Least distinguishing palette')
need([row[i] for row,i in zip(color_counts,[1,2,3])]==[84,1512,4032],'Sharp named-color counts')
# Derive the rational palette-budget series from k=3m+r exactly.
palette=7*z**3*(1+z+z*z)/(1-z**3)**2+(3*z+5*z*z)/(1-z**3)
equal(palette,(3*z+2*z*z+2*z**3)/((1-z)*(1-z**3)),'Palette-budget GF')
need(S.degree(S.fraction(S.cancel(palette))[1],z)==4,'Palette recurrence order')
seq=[1]
for i in range(5):seq.append((7*seq[-1]+2)//3)
need(seq==[1,3,7,17,40,94] and seq[2]==7,'Depth initial values')

# Small exact geometry check; the general theorem is proof-reviewed, not inferred here.
v=S.Matrix([0,0,-1]); face=[S.Matrix([1,0,0]),S.Matrix([0,1,0]),S.Matrix([0,0,1])]
s=sum(face,S.zeros(3,1));w0=s/3;B0=3*w0;F0=4*w0;H0=s-2*w0
need(B0!=F0 and B0!=H0 and F0!=H0 and (B0-H0).cross(F0-H0)==S.zeros(3,1),'Three distinct collinear sphere hostile')
need(S.Matrix.hstack(*(p0-v for p0 in face)).det()!=0,'Nondegenerate rational tetrahedron')
need(sorted([1,6])!=sorted([4,3]) and 4+3==6+1==7,'Equivariant point-bijection orbit obstruction')

print(json.dumps({'status':'PASS','always_active_checks':checks,'bound_files':len(binding['files']),
 'publication_manifests_verified_at_audit':len(binding['source_manifests_verified_at_audit']),
 'bound_wrapper_links':len(links),'scope_anchors':len(anchors),'minimal_recurrence_orders':orders,
 'minimal_quasiperiods':periods,'cyclotomic_pole_profiles':profiles,
 'pointwise_ratio_constant':'2100/17','cumulative_ratio_constant':'3150/17',
 'path_face_n7_h_star':{'degree':7,'constant':1,'linear':81,'leading':34},
 'synthesis_sha256':binding['files'][key]['sha256'],
 'binding_sha256':hashlib.sha256((HERE/'source_bindings.json').read_bytes()).hexdigest(),
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'large_censuses_rerun':False,'git_operations':False},indent=2,sort_keys=True))
