# Independent review of the order-twelve three-chord selector

2026-10-25. **ACCEPTED / PROVED:** the two complete source-family repairs
and the selector in [README.md](README.md). **FINITE-EXACT:** independent
raw-arc signed-union-find verification below, including all sixteen
guard-cube cells and all sixteen literal H-blue sources.
The complete order12 endpoint-four consequence uses the already
audited base/residual partition at its stated thresholds; it does not
claim all-layer black-core connectivity.

## Source scope and inherited part

The exact source consists of P=(0,...,11) directed,11->0, P anti,
the exact incoming interior row I={2,7,9} or{6,7,9}, and both1->10
and2->10. These are precisely the two cells in
[LOW_ORDERS.md](../../port_residual_integral_20261025/LOW_ORDERS.md);
the [independent low-order review](../../port_residual_integral_20261025/LOW_ORDERS_REVIEW.md)
certified their original dimensions and failed words.
There is no hidden constraint on unspecified interior arcs, strongness,
clones, isomorphism type, or reflected-port overlap.

The old three words are uniformly directed and have the literal target
K_P(T). H2 and H3 are universally black in both18-dimensional cells.
The new checker and the independent implementation below re-establish
those facts on the full affine cells. A consistent H-blue subsystem
has dimension3; those small loci are consequences of adding symmetry,
not an initial source restriction.

## Forced chords, exact target and three witnesses

If H is blue, g=r_H r_P is an automorphism of T with the stated four
three-cycles. Reading the exact endpoint rows gives8->11,3->11,
and6->11 in both families. Two applications of g transport(4,9)
through(10,5) to(8,11); one application sends(7,5) to(3,11).
Therefore H-blue implies all three guards, with the directions stated
in the producer. The proof uses an actual automorphism of the source,
not a relation between different intermediate tournaments.

Under those guards, every new join is justified: Q1 uses the first
two guards and8->11; Q2 uses the complements of7->5 and6->11;
Q3 restores the three omitted original path arcs after two flips.
The displayed X,Y identities hold on the actual unordered pairs,
so three K moves reach the original labelled target. Shared vertices
cause no ambiguity because XOR acts on pair sets. Every original
vertex and clone adjacency is retained.

The individual blackness witnesses are correct at their specified
sources. Q1 compares T(0,9)=0 with T(10,11)=1. Q2 compares
S1(0,4)=0 with S1(5,6)=1; the first pair was complemented and the
second retained by Q1. Q3 is tested at B=K_P(T), where
B(1,10)=0 and B(2,11)=1. The exact positional reversals give those
three comparisons. Testing Q3 at B is legitimate because a word's
reversal preserves its path-edge set, so its own K move preserves
its blue/black status.

The guards are sufficient for this ledger; they are not claimed
equivalent to H-blue. If any guard is false, H-blue is impossible
by the proved implication, so the old ledger is three-black.
This proves the selector on every source. The independent check
also splits all three guard bits, obtaining eight nonempty
15-dimensional cells per source family, and checks the selected
ledger in each entire cell.

Each guard is needed for the displayed new ledger: with the other
two true, the missing guard's false direction makes the specified
Q1 or Q2 join false. The six such cells are nonempty and are
covered by the old ledger. Thus no necessary/sufficient direction
has been reversed.

## Fallback and full order12 consequence

The five-black fallback is compatible with the proof. H-blue allows
the same literal source to be remarked by H. The remarked word is
closed because9->0; its endpoint degrees are7 and4. Its interior
port positions become{3,4,9} or{3,8,9}, both outside the failed
order12 list. The remarked middle is T(1,8) and its next chord is
T(10,8). The three branches printed in FALLBACK partition that
source locus. Independent affine replay confirms dimensions2,1,1
for each original cell, and verifies their three-black replacement
of H followed by H2,H3. This auxiliary proof is not needed for the
shorter selector but has the stated exact target.

At order12 the inherited complete source partition has336 cells:
322 base cells and14 closed-endpair residual cells. Only the two
base cells reviewed here lacked a three-black ledger. The residual
sets are I={h,9,10},2<=h<=8. The second-port theorem applies to h2
from n>=11; the complementary h3/h4 repairs and the h>=5 theorem
apply from n>=12. Their paired chord cases cover every residual
completion without an extra branch restriction. These are the
precise domains linked in
[PORT_FRONTIER.md](../../port_residual_integral_20261025/PORT_FRONTIER.md),
section3. Consequently the full closed endpoint-four layer at
order12 is indeed closed by this repair plus the inherited results.

## Reproduction

The producer passed normal and optimized Python with byte-identical
outputs. The following standalone code imports no producer or
predecessor solver. It starts with all66 forward-arc bits and uses
signed union-find, independently of the producer's Gaussian elimination.
It verifies the whole guard cube, H-blue implication, actual-source
anti witnesses, exact edge targets, all six fallback cells, all six
missing-guard cells, and all16 literal H-blue sources.

Both normal and optimized executions return the same PASS output:
two source dimensions18, guard dimensions15, H-blue dimensions3;
16 selector cells,16 literal symmetric sources,6 fallback cells,
and6 omitted-guard controls.

```python
from itertools import combinations,product
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

class Signs:
    def __init__(self,n):
        self.p=list(range(n));self.s=[0]*n;self.rank=[0]*n
    def copy(self):
        d=Signs(0);d.p=self.p[:];d.s=self.s[:];d.rank=self.rank[:];return d
    def find(self,v):
        if self.p[v]!=v:
            r,t=self.find(self.p[v]);self.s[v]^=t;self.p[v]=r
        return self.p[v],self.s[v]
    def join(self,a,b,value):
        a,x=self.find(a);b,y=self.find(b);value^=x^y
        if a==b:return value==0
        if self.rank[a]<self.rank[b]:a,b=b,a
        self.p[b]=a;self.s[b]=value
        if self.rank[a]==self.rank[b]:self.rank[a]+=1
        return True
    def implies(self,a,b,value):
        r,s=self.find(a);t,u=self.find(b);return r==t and s^u==value
    def dimension(self):
        return len({self.find(i)[0]for i in range(len(self.p))})-1

n=12;P=tuple(range(n));pairs=list(combinations(P,2));index={p:i+1 for i,p in enumerate(pairs)};ALL=set(pairs)
OLD=((0,1,10,11,2,3,4,5,6,7,8,9),(1,10,0,11,2,3,4,5,6,7,8,9),
     (1,2,3,4,5,6,7,8,9,10,0,11))
NEW=((0,1,2,3,4,9,10,5,6,7,8,11),(0,1,2,3,4,9,10,5,7,8,11,6),
     (0,1,2,3,4,5,7,8,9,10,11,6))
guards=((4,9),(10,5),(7,5))
cycles=((0,6,2),(1,7,3),(4,10,8),(5,11,9))
g={c[i]:c[(i+1)%3]for c in cycles for i in range(3)}
FALLBACK=(
 ((0,1,8,9,10,11,2,3,4,5,6,7),(1,8,0,9,10,11,2,3,4,5,6,7),(1,10,11,2,3,4,5,6,7,8,0,9)),
 ((8,10,11,2,3,4,5,6,7,1,9,0),(9,8,10,11,2,3,4,5,6,7,1,0),(0,9,1,10,11,2,3,4,5,6,7,8)),
 ((8,1,9,0,10,11,2,3,4,5,6,7),(9,8,1,0,10,11,2,3,4,5,6,7),(0,9,1,10,11,2,3,4,5,6,7,8)))

def edges(q):return {tuple(sorted((a,b)))for a,b in zip(q,q[1:])}
def expr(a,b,flips=frozenset()):
    e=tuple(sorted((a,b)));return index[e],int(a>b)^int(e in flips)
def force(d,a,b,value):
    v,s=expr(a,b);return d.join(v,0,value^s)
def known(d,a,b,value,flips=frozenset()):
    v,s=expr(a,b,flips);return d.implies(v,0,value^s)
def anti(d,q,flips=frozenset()):
    r=dict(zip(q,q[::-1]))
    for a,b in pairs:
        v,s=expr(a,b,flips);w,t=expr(r[b],r[a],flips)
        if not d.join(v,w,s^t):return False
    return True
def source(I):
    d=Signs(67)
    need(anti(d,P) and all(force(d,i,i+1,1)for i in range(11))
         and force(d,11,0,1),'path and anti source')
    need(all(force(d,0,v,int(v not in I))for v in range(1,11))
         and force(d,1,10,1)and force(d,2,10,1),'exact source rows')
    return d
def replay(d,qs,target=P):
    flips=set();blacks=[]
    for q in qs:
        need(sorted(q)==list(P),'word permutation')
        need(all(known(d,a,b,1,flips)for a,b in zip(q,q[1:])),'uniform legality')
        blacks.append(not anti(d.copy(),q,flips))
        flips^=ALL^edges(q)
    need(flips==ALL^edges(target),'literal target')
    return blacks

need({v:dict(zip(OLD[0],OLD[0][::-1]))[11-v]for v in P}==g,'symmetry product')
report=[];literal_count=0;guard_cells=0;fallback_cells=0;missing_guard_cells=0
for I in ((2,7,9),(6,7,9)):
    base=source(I)
    need(base.dimension()==18,'base dimension')
    need(replay(base,OLD)==[False,True,True],'old full-cell blackness')
    for bits in product((0,1),repeat=3):
        d=base.copy()
        need(all(force(d,*edge,bit)for edge,bit in zip(guards,bits)),'nonempty guard cell')
        need(d.dimension()==15,'guard cube dimension')
        choice=NEW if all(bits)else OLD
        need(replay(d,choice)==[True]*3,'complete selector cell')
        if all(bits):
            r1=dict(zip(NEW[0],NEW[0][::-1]));r2=dict(zip(NEW[1],NEW[1][::-1]))
            r3=dict(zip(NEW[2],NEW[2][::-1]));S=ALL^edges(NEW[0]);B=ALL^edges(P)
            need((r1[9],r1[0])==(10,11) and known(d,0,9,0) and known(d,10,11,1),
                 'first literal anti witness')
            need((r2[4],r2[0])==(5,6) and known(d,0,4,0,S) and known(d,5,6,1,S),
                 'second literal anti witness')
            need((r3[10],r3[1])==(2,11) and known(d,1,10,0,B) and known(d,2,11,1,B),
                 'third literal anti witness at B')
        if sum(bits)==2:
            omitted=bits.index(0)
            flips=set() if omitted<2 else ALL^edges(NEW[0])
            a,b=guards[omitted] if omitted<2 else (5,7)
            need(known(d,a,b,0,flips),'necessary guard hostile')
            missing_guard_cells+=1
        guard_cells+=1
    H=base.copy()
    need(anti(H,OLD[0])and H.dimension()==3,'H-blue locus')
    need(all(known(H,a,b,1)for a,b in guards),'H-blue implies guards')
    for a,b in pairs:
        v,s=expr(a,b);w,t=expr(g[a],g[b])
        need(H.implies(v,w,s^t),'exact failed symmetry automorphism')
    for branch,assumptions in enumerate((((1,8,1),),((1,8,0),(10,8,0)),((1,8,0),(10,8,1)))):
        d=H.copy()
        need(all(force(d,*term)for term in assumptions),'fallback cell feasible')
        need(d.dimension()==(2,1,1)[branch],'fallback dimensions')
        need(replay(d,FALLBACK[branch],OLD[0])==[True]*3,'fallback replacement')
        need(replay(d,FALLBACK[branch]+OLD[1:])==[True]*5,'complete fallback')
        fallback_cells+=1
    # Enumerate only the eight models of the proved three-dimensional locus.
    zero,zero_sign=H.find(0)
    roots=sorted({H.find(v)[0]for v in range(67)}-{zero})
    for bits in product((0,1),repeat=len(roots)):
        values=dict(zip(roots,bits));values[zero]=zero_sign
        def bit(a,b):
            v,s=expr(a,b);r,t=H.find(v);return values[r]^t^s
        T=[[0 if a==b else bit(a,b)for b in P]for a in P]
        need(all(T[a][b]^T[b][a]==1 for a,b in pairs),'literal tournament')
        need({v for v in range(1,11)if T[v][0]}==set(I),'literal exact ports')
        def move(R,q):
            E=edges(q)
            return [[0 if a==b else R[a][b]^int(tuple(sorted((a,b)))not in E)
                     for b in P]for a in P]
        def literal(R,qs):
            target=move(R,P);flags=[]
            for q in qs:
                need(all(R[a][b]for a,b in zip(q,q[1:])),'literal joins')
                r=dict(zip(q,q[::-1]))
                flags.append(all(R[a][b]==R[r[b]][r[a]]for a,b in pairs))
                R=move(R,q)
            need(R==target,'literal target');return flags
        need(literal(T,OLD)==[True,False,False] and literal(T,NEW)==[False]*3,
             'complete symmetric literal sources')
        literal_count+=1
    report.append({'ports':I,'base_dimension':18,'guard_dimension':15,'H_blue_dimension':3})
print(json.dumps({'status':'PASS','families':report,'guard_selector_cells':guard_cells,
                  'all_literal_H_blue_sources':literal_count,'fallback_cells':fallback_cells,
                  'omitted_guard_cells':missing_guard_cells},sort_keys=True))

```

