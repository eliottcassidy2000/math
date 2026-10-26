# Independent review of the order-fourteen repair

2026-10-25. **REVIEW PASS:** the source classification, common relocation
rectangle, and complete four-cell repair in `README.md` are correct.
**FINITE-EXACT:** producer normal/optimized replay and the independent
parity-equation control below. No producer correction is needed.

## Exact source and all-completion certification

The variables are all91 unordered-pair orientations on the literal vertices
0,...,13. The source equations impose directed P,13->0, actual P-anti,
the complete incoming row at0, the closed middle12->1, and the stated
2,12 chord. Endpoint degrees(9,4) follow from that row and anti symmetry;
there is no extra degree, strongness, module, or clone filter. Each assigned
cell has affine dimension28. The producer's substitution of its Gaussian
basis leaves a free-variable affine parametrization, so inconsistency of
adding an anti system checks all completions, not selected examples.

Its120-port-set classifier explicitly keeps12 outside the port set and
leaves the2,12 chord free. This is the correct declared finite universe.
The independent control below also allows12 as a port, testing all165
three-element subsets of{2,...,12}; it finds the same four exceptional
points. The proof itself forces12 to be a nonport.

## The failed symmetry and four points

Composing the actual reversals gives exactly

    g=(0,11,9,7,5,3,1)(2,13,12,10,8,6,4).

Both factors are anti on B, so g is an actual automorphism. Its two
seven-element orbits are interchanged by r_P. The known B row at0 has
outdegree5. Closed middle gives B(1,12)=1, and g sends that pair to(0,10),
forcing10 to be a port. The induced first orbit is a transitive cyclic
action on a seven-vertex tournament, hence regular of degree3. Its paired
cyclic differences are{11,1},{9,3},{7,5}. The retained edge0->1 fixes
the first choice; the others give exactly the four claimed port sets.

The chord computation uses ordered pairs correctly:
B(2,12)=B(1,11)=B(0,9). These are nonpath pairs when transferred back
to T, so2->12 in T holds precisely when3 rather than9 is the port.
The complete row at0, g-invariance, orientation reversal, and P-anti
then determine all arcs uniquely. These steps establish at most one
point per cell. The affine and direct carrier constructions agree and
give existence for all four, with no isomorphism quotient.

The consequences T(3,5)=0 and degree8 at every vertex of the second
orbit follow exactly as stated: g^2 sends(3,5) to(0,1), while P-anti
sends degree5 at0 to degree8 at13. The internal vertices4 and12 of P
therefore have degree5 in T.

## Literal ledger and degree receipts

The new joins follow from the exact port row and T(3,5)=0. The edge
rectangle remains valid although its two relocations share vertex3.
XOR of the three path-edge sets is E(P), so three complements reach
the exact labelled B. Legality has been checked separately at each
actual preceding tournament; no Hamiltonian word is changed for free.

The degree receipt is consistent with the proof, but its two ways of
checking Q3 must be distinguished. At the three actual pre-move sources,
the endpoint degrees are respectively(9,5),(8,7),(6,4), with sums14,15,10.
The proof's sum16 is the alternative Q3 test at its other endpoint B,
where its endpoint degrees are(8,8). Every sum differs from13. Blue/black
status is unchanged across a word's own K edge because its reversal
preserves that word's path-edge set. Thus either Q3 calculation is valid;
the receipt does not contradict the proof.

The inherited first two moves are uniformly black on the four full
28-dimensional cells, and their final F has the same status at B as at
its actual source. The conditional rule therefore covers every completion:
old ledger when F is black, new ledger when F is blue. The two-pair
perturbation stays in its source cell and invalidates the new join3->5,
while the old ledger succeeds. It is a genuine missing-hypothesis hostile.

## Stronger guard extension

The separate [GUARD_EXTENSION.md](GUARD_EXTENSION.md) records a stronger
result found during this audit: T(3,5)=0 alone suffices for the new ledger
on each full cell. Its three fixed-pair witnesses avoid the degree and
symmetry assumptions. The symmetry classification still establishes that
T(3,5)=1 rules out blue F, so the two-ledger decision can use that one bit.
This extension does not modify the reviewed original producer.

## Reproduction

Producer `check.py` passes normally and with `python -O`; both outputs
agree. The independent code below imports no producer. It uses parity
components rather than the producer's Gaussian reduction, covers all165
port sets with free chord, verifies all four full source dimensions and
unique exceptions, reconstructs the literal matrices independently,
checks the actual-stage degree certificates and rectangle, and tests
the eight27-dimensional guard half-cells. All gates survive optimization.

From the repository root, run the embedded block with:

    @'
    from pathlib import Path
    p=Path('05-knowledge/results/lca_ports_reverse_20261025/order14/INDEPENDENT_REVIEW.md')
    code=p.read_text(encoding='utf-8').split('```python\n',1)[1].split('\n```',1)[0]
    exec(compile(code,str(p),'exec'))
    '@ | python -

Replace the final command by `python -O -` for the optimized replay.

```python
from itertools import combinations
from pathlib import Path
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

n=14; P=tuple(range(n)); pairs=tuple(combinations(P,2))
index={p:i for i,p in enumerate(pairs)}; const=len(pairs)
F=(0,13)+tuple(range(1,13))
ports=((3,5,10),(3,7,10),(5,9,10),(7,9,10))
new=((0,1,2,13)+tuple(range(3,13)),
     (4,0,1,2,13,3)+tuple(range(5,13)),
     (4,0,1,2,3)+tuple(range(5,14)))

def pair(a,b): return (min(a,b),max(a,b))
def path(q): return {pair(a,b) for a,b in zip(q,q[1:])}
def reversal(q): return dict(zip(q,q[::-1]))
allpairs=set(pairs); target_flip=allpairs-path(P)

class Parity:
    def __init__(self):
        self.parent=list(range(const+1)); self.delta=[0]*(const+1); self.ok=True
    def copy(self):
        out=Parity(); out.parent=self.parent[:]; out.delta=self.delta[:]; out.ok=self.ok
        return out
    def root(self,i):
        bit=0
        while self.parent[i]!=i:
            bit^=self.delta[i]; i=self.parent[i]
        return i,bit
    def equation(self,i,j,bit):
        a,x=self.root(i); b,y=self.root(j)
        if a==b:
            self.ok &= x^y==bit
        else:
            self.parent[b]=a; self.delta[b]=x^y^bit
    def force(self,a,b,bit):
        self.equation(index[pair(a,b)],const,bit^(a>b))
    def forced(self,a,b,bit):
        x,u=self.root(index[pair(a,b)]); y,v=self.root(const)
        return x==y and u^v==(bit^(a>b))
    def anti(self,q,flips=set()):
        r=reversal(q)
        for a,b in pairs:
            x,y=r[b],r[a]
            self.equation(index[a,b],index[pair(x,y)],
                          (x>y)^((a,b) in flips)^(pair(x,y) in flips))
    def dimension(self):
        need(self.ok,'nonempty affine source')
        return len({self.root(i)[0] for i in range(const+1)})-1
    def model(self):
        need(self.ok,'model exists')
        root,bit=self.root(const); values={root:bit}; T=[[0]*n for _ in P]
        for a,b in pairs:
            r,v=self.root(index[a,b]); values.setdefault(r,0)
            T[a][b]=values[r]^v; T[b][a]=1-T[a][b]
        return T

def source(I,chord=None,guard=None):
    s=Parity(); s.anti(P)
    for a,b in zip(P,P[1:]): s.force(a,b,1)
    s.force(13,0,1); s.force(12,1,1)
    for v in range(1,13): s.force(v,0,int(v in I))
    if chord is not None: s.force(2,12,chord)
    if guard is not None: s.force(3,5,guard)
    need(s.ok,'source cell nonempty')
    return s

def old(chord):
    w=tuple(range(2,12))
    return ((12,1,13,0)+w,(13,12,1,0)+w,F) if chord else (
        (12,)+w+(1,13,0),(13,12)+w+(1,0),F)

def certify(s,words):
    flips=set(); feasible=[]
    for q in words:
        need(sorted(q)==list(P),'Hamilton word')
        need(all(s.forced(a,b,1^(pair(a,b) in flips)) for a,b in zip(q,q[1:])),
             'all-completion directedness')
        trial=s.copy(); trial.anti(q,flips); feasible.append(trial.ok)
        flips ^= allpairs-path(q)
    need(flips==target_flip,'literal edge XOR')
    return feasible

def move(T,q):
    E=path(q)
    return [[0 if a==b else T[a][b]^int(pair(a,b) not in E) for b in P] for a in P]

def anti(T,q):
    r=reversal(q)
    return all(T[a][b]==T[r[b]][r[a]] for a,b in pairs)

def literal(T,words):
    B=move(T,P); sums=[]
    for q in words:
        need(all(T[a][b] for a,b in zip(q,q[1:])),'literal directedness')
        need(not anti(T,q),'literal blackness')
        sums.append([sum(T[q[0]]),sum(T[q[-1]])]); T=move(T,q)
    need(T==B,'literal labelled target')
    return sums

classified=[]; unique={}
for I in combinations(range(2,13),3):
    s=source(I); s.anti(F,target_flip)
    if s.ok:
        need(s.dimension()==0,'exception is one point')
        T=s.model(); classified.append((I,T[2][12])); unique[I]=T
need(classified==[(I,int(3 in I)) for I in ports],'all165 port sets give exactly four exceptions')
receipt=json.loads(Path('05-knowledge/results/lca_ports_reverse_20261025/order14/normal.json').read_text())
rows_by_I={tuple(row['ports']):row['source_rows'] for row in receipt['families']}
rP=reversal(P); rF=reversal(F); g={v:rF[rP[v]] for v in P}
cycles=[]; unseen=set(P)
while unseen:
    start=min(unseen); c=[]; v=start
    while v in unseen: unseen.remove(v); c.append(v); v=g[v]
    cycles.append(c)
need(cycles==[[0,11,9,7,5,3,1],[2,13,12,10,8,6,4]],'actual two seven-cycles')

for I in ports:
    chord=int(3 in I); s=source(I,chord)
    need(s.dimension()==28 and certify(s,old(chord))==[False,False,True],
         'full source and inherited ledger')
    T=unique[I]; B=move(T,P)
    need(T[12][1]==1 and T[2][12]==chord,'explicit literal branch predicates')
    need([sum(T[a][b]<<b for b in P) for a in P]==rows_by_I[I],'independent literal source recovery')
    need(all(B[a][b]==B[g[a]][g[b]] for a,b in pairs),'literal automorphism')
    need(T[3][5]==0 and sum(B[4])==sum(B[13])==8,'forced repair data')
    need(literal(T,new)==[[9,5],[8,7],[6,4]],'actual-stage degree receipts')
    for guard in (0,1):
        guarded=source(I,chord,guard)
        need(guarded.dimension()==27,'guard is a genuine free source bit')
        need(certify(guarded,new if guard==0 else old(chord))==[False]*3,
             'one-bit selector covers entire guarded cell')

hostile=[row[:] for row in unique[ports[0]]]
for a,b in ((3,5),(8,10)): hostile[a][b]^=1; hostile[b][a]^=1
need(anti(hostile,P) and hostile[12][1] and hostile[2][12],'hostile preserves symmetry and branch')
need({v for v in range(1,13) if hostile[v][0]}==set(ports[0]),'hostile preserves exact ports')
need(all(hostile[a][b] for a,b in zip(new[0],new[0][1:])),'hostile first move legal')
need(move(hostile,new[0])[3][5]==0,'unguarded new ledger is illegal')
literal(hostile,old(1))
print(json.dumps({'status':'PASS','classified_port_sets':165,
                  'full_source_cells':4,'full_source_dimension':28,
                  'unique_failed_points':4,'guard_half_cells':8,'guard_dimension':27,
                  'actual_stage_degree_pairs':[[9,5],[8,7],[6,4]],
                  'other_end_final_degree_pair':[8,8]},indent=2))
```
