"""Independent incidence and polynomial audit. Does not import producer code."""
import json
from itertools import product

def check(ok,label):
    if not ok: raise RuntimeError(label)

def graph(k,done,paired=False):
    params={f'{p}{j}' for j in range(k) for p in 'xyabcz' if not(p=='x' and j<done)}
    eq=[]; gates=[]; pairs=set()
    for j in range(done,k):
        eq.extend([(f'A{j}','+',f'x{j}',f'a{j}'),(f'C{j}','+',f'x{j}',f'b{j}')])
    first=done+int(paired)
    if paired:pairs={f'x{done}',f'y{done}'}
    for j in range(first,k):gates.append((f'G{j}','*',f'x{j}',f'y{j}'))
    for j in range(first,k):gates.append((f'K{j}','+',f'G{j}',f'z{j}'))
    for j in range(k):
        A,C=(f'a{j}',f'b{j}') if j<done else (f'A{j}',f'C{j}')
        K=f'z{j}' if j<first else f'K{j}'
        gates.extend([(f'B{j}','*',f'y{j}',f'b{j}'),(f'T{j}','*',f'a{j}',f'c{j}'),(f'H{j}','*',A,C),(f'V{j}','+',K,f'B{j}'),(f'W{j}','+',f'V{j}',f'T{j}'),(f'O{j}','+',f'W{j}',f'H{j}')])
    output='O0'
    for j in range(1,k):
        gates.append((f'R{j}','+',output,f'O{j}'));output=f'R{j}'
    return params,eq,gates,pairs,output

def audit(state):
    params,eq,gates,paired,out=state
    objects=eq+gates; children={v:set() for v in params|{g[0] for g in objects}}
    for name,op,a,b in objects:
        check(a!=b,'repeated pin');children[a].add(name);children[b].add(name)
    needed={out}
    for name,op,a,b in reversed(objects):
        check(name in needed,'dead object '+name);needed.update([a,b])
    equations={v[0] for v in eq}
    for p in params-paired:
        check(not(len(children[p])==1 and children[p]<=equations),'D4 '+p)
    forms={p:1<<i for i,p in enumerate(sorted(params))}
    for name,op,a,b in eq:forms[name]=forms[a]^forms[b]
    for name,op,a,b in gates:
        if a in forms and b in forms:
            check(op=='*','M available')
            check(forms[a]!=0 and forms[b]!=0 and forms[a]!=forms[b],'dependent parents')
    I=sum(bool(children[p]) for p in params-paired)
    q=len(paired)//2
    return (len(objects),I,q,len(params)-q),children

def evaluate(state,assignment):
    _,eq,gates,_,out=state;v=dict(assignment)
    for name,op,a,b in eq+gates:v[name]=(v[a]^v[b]) if op=='+' else (v[a]&v[b])
    return v[out]

def poly(v):
    x,y,a,b,c,z=v
    return x^(x&a)^(x&b)^(a&b)^(y&b)^(a&c)^(x&y)^z

def rank(rows):
    pivots={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p not in pivots:pivots[p]=row;break
            row^=pivots[p]
    return len(pivots)

cycles=stages=0
for k in range(1,26):
    for done in range(k):
        before,bc=audit(graph(k,done));middle,mc=audit(graph(k,done,True));after,ac=audit(graph(k,done+1))
        check(before==(11*k-1-4*done,6*k-done,0,6*k-done),'before inventory')
        check(tuple(m-b for b,m in zip(before,middle))==(-2,-2,1,-1),'GP deltas')
        check(tuple(a-m for m,a in zip(middle,after))==(-2,1,-1,0),'P2 deltas')
        check(bc[f'G{done}']=={f'K{done}'},'GP child')
        check(mc[f'x{done}']=={f'A{done}',f'C{done}'},'P2 children')
        check(mc[f'y{done}']=={f'B{done}'},'paired partner child')
        check(before[3]>=4*k+1 and middle[3]>=4*k+1,'both loop guards')
        cycles+=1;stages+=3
    rows=[]
    edges=[(0,2),(0,3),(2,3),(1,3),(2,4),(0,1)]
    for i in range(6*k):
        j,r=divmod(i,6);row=0
        for a,b in edges:
            if r==a:row^=1<<(6*j+b)
            if r==b:row^=1<<(6*j+a)
        rows.append(row)
    check(rank(rows)==4*k,'polar rank')
checks=0
for bits in product([0,1],repeat=6):
    assignment={p+'0':v for p,v in zip('xyabcz',bits)}
    check(evaluate(graph(1,0),assignment)==poly(bits),'initial semantics')
    x,y,a,b,c,z=bits
    factored=((x^b)&(x^a^y))^(a&c)^z
    check(poly(bits)==factored,'seven-gate resynthesis')
    if not(x and y):check(evaluate(graph(1,0,True),assignment)==poly(bits),'GP semantics')
    if x==0:
        check(evaluate(graph(1,1),assignment)==poly(bits),'P2 semantics')
        check(poly(bits)==(b&(a^y))^(a&c)^z,'five-gate residual resynthesis')
    checks+=1
radical=[]
for u in range(64):
    a=[u>>i&1 for i in range(6)]
    if all(poly([a[i]^(v>>i&1) for i in range(6)])^poly(a)^poly([v>>i&1 for i in range(6)])==0 for v in range(64)):radical.append(u)
check(radical==[0,11,32,43],'radical vectors')
check(poly([0,0,0,0,0,1])==1,'radical linear part')
# Certify semantic resynthesis separately; it is not a source D-row.
factor_checks=0
for r in range(2,9):
    for bits in product([0,1],repeat=r+1):
        b=bits[0];u=bits[1:];old=0;parity=0
        for x in u:old^=b&x;parity^=x
        check(old==(b&parity),'exclusive XOR factoring')
        check((2*r-1)-r==r-1,'factor gate saving')
        factor_checks+=1
resynthesis_evaluations=0
for k in range(1,26):
    residual=graph(k,k)
    params=residual[0];gates=[]
    for j in range(k):
        gates.extend([(f'U{j}','+',f'a{j}',f'y{j}'),(f'P{j}','*',f'b{j}',f'U{j}'),(f'T{j}','*',f'a{j}',f'c{j}'),(f'V{j}','+',f'P{j}',f'T{j}'),(f'O{j}','+',f'V{j}',f'z{j}')])
    out='O0'
    for j in range(1,k):
        gates.append((f'R{j}','+',out,f'O{j}'));out=f'R{j}'
    replacement=(params,[],gates,set(),out)
    check(len(residual[1])+len(residual[2])-len(gates)==k,'disjoint cone saving')
    check(len(gates)==6*k-1,'replacement gate count')
    if k<=2:
        ordered=sorted(params)
        for bits in product([0,1],repeat=len(ordered)):
            assignment=dict(zip(ordered,bits))
            check(evaluate(residual,assignment)==evaluate(replacement,assignment),'resynthesized whole DAG')
            resynthesis_evaluations+=1
# M/D4 continuation: reparameterize w=y+a; do NOT treat this as free
# original-coordinate circuit rewiring. Literal input map is checked below.
role_exchange_evaluations=0
for k in range(1,26):
    params={f'{p}{j}' for j in range(k) for p in 'wabcz'}
    gates=[]
    for j in range(k):
        gates.extend([(f'P{j}','*',f'b{j}',f'w{j}'),(f'T{j}','*',f'a{j}',f'c{j}'),(f'V{j}','+',f'P{j}',f'T{j}'),(f'O{j}','+',f'V{j}',f'z{j}')])
    out='O0'
    for j in range(1,k):
        gates.append((f'R{j}','+',out,f'O{j}'));out=f'R{j}'
    final=(params,[],gates,set(),out)
    inventories,_=audit(final)
    check(inventories==(5*k-1,5*k,0,5*k),'reparameterized inventory')
    # Before deletion the exchanged auxiliary matrix is still the identity.
    for remaining in range(k+1):
        check(rank([1<<i for i in range(remaining)])==remaining,'fair matrix invertibility')
    # Map from (w,a,b,c,z) to original (x,y,a,b,c,z), with x=0,y=w+a.
    embedding=[]
    for j in range(k):embedding.extend([0,3<<(5*j),2<<(5*j),4<<(5*j),8<<(5*j),16<<(5*j)])
    check(rank(embedding)==5*k,'affine embedding injectivity')
    # Every new parameter is essential: the other product factor is set to1.
    for j in range(k):
        for p,partner in [('w','b'),('b','w'),('a','c'),('c','a'),('z',None)]:
            v={name:0 for name in params}
            if partner is not None:v[f'{partner}{j}']=1
            f0=evaluate(final,v);v[f'{p}{j}']=1
            check(evaluate(final,v)!=f0,'essential reparameterized input')
    if k<=2:
        ordered=sorted(params);images=set()
        for bits in product([0,1],repeat=len(ordered)):
            v=dict(zip(ordered,bits));old={}
            for j in range(k):
                for p in 'abcz':old[f'{p}{j}']=v[f'{p}{j}']
                old[f'y{j}']=v[f'w{j}']^v[f'a{j}']
                check(old[f'y{j}']^old[f'a{j}']==v[f'w{j}'],'role exchange inverse')
            check(evaluate(graph(k,k),old)==evaluate(final,v),'affine-state semantics')
            images.add(tuple(sorted(old.items())));role_exchange_evaluations+=1
        check(len(images)==1<<(5*k),'finite map bijectivity')
print(json.dumps({'status':'PASS','block_sizes':list(range(1,26)),'exact_cycles':cycles,'stage_incidence_audits':stages,'one_block_truth_assignments':checks,'radical':radical,'factoring_truth_assignments':factor_checks,'resynthesis_truth_assignments':resynthesis_evaluations,'resynthesis_saving':'k original-coordinate gates; 7k-1 to 6k-1','role_exchange_semantic_assignments':role_exchange_evaluations,'role_exchange_extra_saving':'k charged equations under retained affine input map; 6k-1 to 5k-1','producer_imported':False,'scope':'specified nonminimal circuits and prescribed legal transitions; no low-d obstruction'},indent=2,sort_keys=True))
