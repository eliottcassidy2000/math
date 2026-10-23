#!/usr/bin/env python3
"""Independent affine-lattice audit: direct permutations and field norms."""
import ast,json,math,itertools,pathlib
# Only load the independent auditor's function definitions, never producer code.
parent=pathlib.Path(__file__).resolve().parent.parent
candidates=[parent/name/'audit.py' for name in ('session43-audit-cyclotomic','audit-cyclotomic')]
source=next((p for p in candidates if p.is_file()),None)
if source is None:raise RuntimeError('Missing independent audit-cyclotomic/audit.py or session43-audit-cyclotomic/audit.py sibling')
tree=ast.parse(source.read_text()); defs=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(source),'exec'))

def gens(n):
    G=[eulerian(n)]
    for s in range(1,(n-1)//2+1):
        z=[0]*n
        for i,c in enumerate(eulerian(n-2*s)):
            for j in range(2*s+1):z[i+j]+=2*c*math.comb(2*s,j)*(-1)**j
        G.append(z)
    return G

def coeffs(f,G):
    # Recover coefficients by exact lowest nonzero Taylor term at1.
    r=[a-b for a,b in zip(f,G[0])];out=[]
    for s,C in enumerate(G[1:],1):
        value=sum(math.comb(i,2*s)*r[i] for i in range(2*s,len(r)))
        den=2*math.factorial(len(r)-2*s)
        check(value%den==0,'integer excess');a=value//den;out.append(a)
        r=[v-a*w for v,w in zip(r,C)]
    check(not any(r),'full excess reconstruction')
    return out

def block(n,s):
    pairs=list(itertools.combinations(range(n),2));code=(1<<len(pairs))-1
    if s:code^=1<<pairs.index((0,2*s))
    return direct_f(matrix(n,code))

count=0;hostile=(7,38,30,38,7);fives=set()
for n in range(1,6):
    G=gens(n)
    for code in range(1<<(n*(n-1)//2)):
        f=direct_f(matrix(n,code));cs=coeffs(f,G)
        check(all(c>=0 for c in cs),'nonnegative coordinates');count+=1
        if n==5:fives.add(tuple(f))
check(hostile not in fives,'nonrealizable hostile')
check(coeffs(hostile,gens(5))==[3,0],'hostile cycle counts')
check(list(hostile)==[3*x-2*y for x,y in zip(block(5,1),eulerian(5))],'affine certificate')
family=0
for n in range(1,9):
    G=gens(n)
    for s in range(len(G)):
        actual=block(n,s);expect=G[0][:]
        for j in range(1,s+1):expect=[x+math.comb(2*s-1,2*j-1)*y for x,y in zip(expect,G[j])]
        check(actual==expect,'actual block family');family+=1
# Norm valuation table for explicit generators and their integral block mixtures.
rows=[]
for n in range(1,32):
    G=gens(n);B=[]
    for s in range(len(G)):
        f=G[0][:]
        for j in range(1,s+1):f=[x+math.comb(2*s-1,2*j-1)*y for x,y in zip(f,G[j])]
        B.append(f)
    for p in (2,3,5,7,11):
        vg=[vprime(norm(f,p),p) for f in G];vb=[vprime(norm(f,p),p) for f in B]
        finite=lambda a:min((x for x in a if x is not None),default=None)
        exact=finite(vg);check(exact==finite(vb),'norm basis minima')
        E=n-max(ds(m,p) for m in range(1,n+1))
        if p>2:check(exact==2*((E+1)//2),'tested odd-prime conjecture')
        elif n%2 and n>=3:check(exact==bound(n),'binary sharp')
        elif n%2==0:check(exact is None,'even binary infinity')
        rows.append({'n':n,'p':p,'exact_valuation':exact})
print(json.dumps({'status':'PASS','direct_tournament_coordinate_checks':count,'direct_block_family_checks':family,'norm_valuation_cells':len(rows),'hostile':hostile,'table':rows},sort_keys=True,indent=2))
