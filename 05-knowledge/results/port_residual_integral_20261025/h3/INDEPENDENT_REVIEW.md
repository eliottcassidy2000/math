# Independent review of the h=3 forward repair

2026-10-25. **REVIEW PASS:** the all-order proof in `README.md` is sound
for its stated source family and n>=12. **FINITE-EXACT:** producer replay
and the independent Gaussian/degree controls below. No global black-core
connectivity claim follows from this review.

## Source universe and quantifiers

The source is a labelled tournament T with directed P=(0,...,n-1),
closing arc (n-1)->0, actual positional reversal r_P anti, exactly
N^-(0)={3,n-3,n-2,n-1}, and 1->3. No unspecified middle edge is fixed.
These conditions are all affine equations over F2 on unordered-pair
orientation bits. The imported producer engine imposes exactly those
equations; it adds no degree, clone, isomorphism, or middle-orientation
filter. Its source-consistency gate prevents vacuous certification.

For every candidate word, directedness is checked by implication from
that source system after the actual preceding flips. Blackness is the
inconsistency of adding the actual positional anti equations at that
stage. Thus each affine cell quantifies every source completion, whereas
the seeded literal models only corroborate that exact calculation.
The independent code below replaces the imported parity-union engine
with Gaussian elimination on augmented integer bit rows.

## All-order argument audited

The three words are Hamiltonian without repeated vertices for n>=12.
Their joins and edge-set XOR hold literally. The last move's blackness
can be checked at B=K_P(T): positional reversal preserves its own path
edge set, so its anti status is the same at both ends of that K move.
The first witness maps (0,n-1) to (1,2), whose T values are0 and1.

For Q3, the displayed composition g3=r_Q3 r_P has the claimed cycles.
If it were an automorphism of B, its orbit of0 would be a transitive
induced tournament. For odd n this is all vertices and d_B^+(0)=5,
forcing n=11. For even n its size is n/2 and its only outgoing ports
from0 are1 and3, forcing n=10. The exact intersection with the port
row is essential; global regularity alone would not settle the even case.

For Q2, the comparison S=K_Q1(T) versus B flips precisely
F={12,23,13,2(n-1)}. In B these orientations give the degree changes
delta_2=-1, delta_(n-1)=+1 and zero elsewhere. If r_Q2 were anti on S,
the degree equation really is

    beta_(gv)-beta_v = -delta_(r_P v)-delta_(gv),
    g=r_Q2 r_P.

This follows by applying the S anti equation at r_P(v) and the B anti
equation at v; it does not assume g is an automorphism of B. The claimed
increments and the two odd-order cycle tables agree. Their offset sums
are0 for n=1 mod4 and -n for n=3 mod4. Since beta_0=5, the total degree
sum forces n=11 or9 respectively, outside the theorem's range.

At n=2 mod4, summing increments around the displayed cycle through0
gives -1. This already contradicts returning to the same degree.

The load-bearing case is n=4t. The displayed C has size2t, its
complement is r_P(C), and both are g-invariant. Only23 from F is internal
to C, and no F pair is internal to its complement. Therefore
S[r_P(C)]=B[r_P(C)], whereas S[C] reverses the B arc2->3. Applying the
two anti maps between these sets yields

    alpha_(gv)-alpha_v = -delta_C(gv),
    delta_C(2)=-1, delta_C(3)=+1.

This is exactly equation(10), including its signs. It is an induced-degree
identity, not an unsupported internal automorphism assertion. The actual
port row gives alpha_0=2. In the cycle order, precisely t vertices have
degree1 and t have degree2. Their sum3t must equal binom(2t,2), forcing
t=2. This contradicts n>=12. Cross-orbit arcs need not be known.

No proof gap or producer correction was found. The boundary examples
at n=8 (Q2 blue) and n=10 (Q3 blue) concern this template; they are not
counterexamples to connectivity. Exact clone identities and all middle
arcs remain in the literal tournaments throughout the construction.

## Replayed and independent controls

The producer `check.py` was run normally and with `python -O`. It checks
all affine completions for n=12..100,120,160, two literal models per order,
and its permutation/degree identities through500, retaining both small
blue hostiles and the old illegal-template hostile.

The standalone code below imports no producer. It certifies all source
completions for both new h3 and h4 families at n=12..40, checks the two
h3 blue boundary examples and h4's extra-port hostile, and verifies the
induced/global degree contradictions and h4 witness images for every
n=12..512. Gaussian consistency is finite evidence. The all-order
conclusion rests on the audited proof above.

To run the embedded control from the repository root, use the following
PowerShell command; replace `python` by `python -O` for optimized mode:

    @'
    from pathlib import Path
    p=Path('05-knowledge/results/port_residual_integral_20261025/h3/INDEPENDENT_REVIEW.md')
    code=p.read_text(encoding='utf-8').split('```python\n',1)[1].split('\n```',1)[0]
    exec(compile(code,str(p),'exec'))
    '@ | python -

```python
from itertools import combinations
import json

def need(ok, text):
    if not ok:
        raise RuntimeError(text)

def edge(a,b):
    return (a,b) if a<b else (b,a)

def path(q):
    return {edge(a,b) for a,b in zip(q,q[1:])}

def reverse(q):
    return dict(zip(q,q[::-1]))

def words(n,h):
    if h==3:
        u=tuple(range(3,n-3)); tail=tuple(range(n-3,n))
        return [(0,1)+u+tail+(2,), (1,)+u+(0,)+tail+(2,),
                (1,2)+u+(0,)+tail]
    v=tuple(range(5,n))
    return [(0,)+v+(1,2,3,4), (2,0)+v+(1,3,4), (2,0,1,3,4)+v]

class Gaussian:
    def __init__(self,n):
        self.index={p:i for i,p in enumerate(combinations(range(n),2))}
        self.constant=1<<len(self.index)
        self.variables=self.constant-1
        self.basis={}
        self.ok=True
    def copy(self):
        out=object.__new__(Gaussian)
        out.index=self.index; out.constant=self.constant
        out.variables=self.variables; out.basis=self.basis.copy(); out.ok=self.ok
        return out
    def bit(self,a,b):
        return 1<<self.index[edge(a,b)]
    def reduce(self,row):
        while row & self.variables:
            pivot=(row & self.variables)&-(row & self.variables)
            if pivot not in self.basis:
                break
            row ^= self.basis[pivot]
        return row
    def equation(self,mask,rhs):
        row=self.reduce(mask | (self.constant if rhs else 0))
        if row & self.variables:
            pivot=(row & self.variables)&-(row & self.variables)
            self.basis[pivot]=row
        elif row:
            self.ok=False
    def force(self,a,b,value):
        self.equation(self.bit(a,b),value^(a>b))
    def forced(self,a,b,value):
        row=self.bit(a,b) | (self.constant if value^(a>b) else 0)
        return self.reduce(row)==0
    def anti(self,q,flips=set()):
        r=reverse(q)
        for a,b in self.index:
            x,y=r[b],r[a]
            self.equation(self.bit(a,b)^self.bit(x,y),
                          (x>y)^((a,b) in flips)^(edge(x,y) in flips))

def source(n,h,extra=()):
    s=Gaussian(n); s.anti(tuple(range(n)))
    for v in range(n-1):
        s.force(v,v+1,1)
    incoming={h,n-3,n-2,n-1}|set(extra)
    for v in range(1,n):
        s.force(v,0,int(v in incoming))
    a,b=(1,3) if h==3 else (3,1)
    s.force(a,b,1)
    need(s.ok,('nonempty exact source',n,h,extra))
    return s

def cell(n,h):
    s=source(n,h); flips=set(); allpairs=set(s.index); feasible=[]
    for q in words(n,h):
        need(sorted(q)==list(range(n)),('Hamilton word',n,h))
        for a,b in zip(q,q[1:]):
            need(s.forced(a,b,1^(edge(a,b) in flips)),('unforced join',n,h,a,b))
        trial=s.copy(); trial.anti(q,flips); feasible.append(trial.ok)
        flips ^= allpairs^path(q)
    need(flips==allpairs^path(tuple(range(n))),('literal XOR',n,h))
    return feasible

def orbit(g,start=0):
    result=[]; at=start
    while at not in result:
        result.append(at); at=g[at]
    need(at==start,'permutation orbit closes')
    return result

def walk_degrees(g,c,start,increment):
    value=start; total=0
    for v in c:
        total+=value
        value+=increment(v)
    return total,value

def geometry(n):
    q1,q2,q3=words(n,3); ports={1,3,n-3,n-2,n-1}
    r3=reverse(q3); g3={v:r3[n-1-v] for v in range(n)}
    c3=orbit(g3)
    need(len(c3)==(n if n%2 else n//2),'Q3 orbit size')
    d0=len(set(c3)&ports)
    need(d0==(5 if n%2 else 2),'Q3 actual induced degree')
    need(len(c3)*d0!=len(c3)*(len(c3)-1)//2,'Q3 contradiction')
    r2=reverse(q2); g={v:r2[n-1-v] for v in range(n)}; c=orbit(g)
    delta={2:-1,n-1:1}
    inc=lambda v:-delta.get(n-1-v,0)-delta.get(g[v],0)
    if n%2:
        need(len(c)==n,'Q2 odd transitivity')
        total,last=walk_degrees(g,c,5,inc)
        need(last==5 and total!=n*(n-1)//2,'Q2 global degree contradiction')
    elif n%4==2:
        need(sum(inc(v) for v in c)!=0,'Q2 degree-cycle contradiction')
    else:
        C=set(c); R={n-1-v for v in C}
        need(R==set(range(n))-C and {g[v] for v in R}==R,'Q2 complementary orbits')
        F={(1,2),(2,3),(1,3),(2,n-1)}
        need({p for p in F if set(p)<=C}=={(2,3)},'one internal defect')
        need(not any(set(p)<=R for p in F),'no complementary defect')
        need(C&ports=={3,n-2},'actual induced degree two')
        internal={2:-1,3:1}
        total,last=walk_degrees(g,c,2,lambda v:-internal.get(g[v],0))
        need(last==2 and total!=len(c)*(len(c)-1)//2,'Q2 induced degree contradiction')
    a,b,c=words(n,4); ra,rb,rc=map(reverse,(a,b,c))
    need((ra[3],ra[0])==(5,4),'h4 first witness image')
    need((rb[n-1],rb[3])==(6,0),'h4 middle witness image')
    need((rc[3],rc[0])==(n-4,n-2),'h4 last witness image')
    need(path(a)^path(b)^path(c)==path(tuple(range(n))),'h4 literal XOR')

for n in range(12,41):
    for h in (3,4):
        need(cell(n,h)==[False]*3,('all-completion black certificate',n,h))
need(cell(8,3)[1], 'h3 lower-order Q2 blue hostile')
need(cell(10,3)[2], 'h3 lower-order Q3 blue hostile')
extra=source(12,4,(5,))
need(extra.forced(0,5,0),'h4 extra incoming port blocks first join')
for n in range(12,513):
    geometry(n)
print(json.dumps({'status':'PASS','independent_gaussian_cells':58,
                  'affine_orders':[12,40],'degree_and_witness_orders':[12,512],
                  'h3_blue_hostiles':[[8,2],[10,3]],'h4_extra_port_hostile':5},indent=2))
```
