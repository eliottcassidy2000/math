"""Exact reduced-DAG GP/P2 cycles. Standard library, checks survive -O."""
import json
from itertools import product

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def block_names(j):return {c:f'{c}{j}' for c in 'xyabcz'}

def make(k):
    pars=set(); objs=[]
    for j in range(k):
        n=block_names(j);pars.update(n.values())
        objs += [(f'A{j}','eqxor',n['x'],n['a']),(f'C{j}','eqxor',n['x'],n['b'])]
    for j in range(k):objs.append((f'G{j}','and',f'x{j}',f'y{j}'))
    for j in range(k):objs.append((f'K{j}','xor',f'G{j}',f'z{j}'))
    for j in range(k):
        objs += [(f'B{j}','and',f'y{j}',f'b{j}'),(f'T{j}','and',f'a{j}',f'c{j}'),(f'H{j}','and',f'A{j}',f'C{j}'),(f'O{j}a','xor',f'K{j}',f'B{j}'),(f'O{j}b','xor',f'O{j}a',f'T{j}'),(f'O{j}','xor',f'O{j}b',f'H{j}')]
    out='O0'
    for j in range(1,k):
        new=f'R{j}';objs.append((new,'xor',out,f'O{j}'));out=new
    return {'pars':pars,'objs':objs,'out':out,'pairs':set()}

def children(s,name):return {out for out,op,a,b in s['objs'] if name==a or name==b}
def stats(s):
    paired={x for pair in s['pairs'] for x in pair}
    I=sum(bool(children(s,p)) for p in s['pars']-paired)
    m=len(s['pars']);q=len(s['pairs'])
    return (len(s['objs']),I,q,m-q)
def value(s,assignment):
    vals={0:0,1:1,**assignment}
    for out,op,a,b in s['objs']:
        vals[out]=vals[a]&vals[b] if op=='and' else vals[a]^vals[b]
    return vals[s['out']]
def substitute(s,initial):
    subst=dict(initial);remaining=[];deleted=[]
    def resolve(x):
        while x in subst:x=subst[x]
        return x
    for out,op,a,b in s['objs']:
        if out in subst:deleted.append(out);continue
        a,b=resolve(a),resolve(b); v=None
        if op=='and':
            if a==0 or b==0:v=0
            elif a==1:v=b
            elif b==1 or a==b:v=a
        else:
            if a==0:v=b
            elif b==0:v=a
            elif a==b:v=0
            elif isinstance(a,int) and isinstance(b,int):v=a^b
        if v is not None:subst[out]=v;deleted.append(out)
        else:remaining.append((out,op,a,b))
    s['objs']=remaining;s['out']=resolve(s['out'])
    return deleted

def closure_checks(s):
    # Every object is in output ancestry; no unary/constant/repeated pins.
    seen={s['out']}
    for out,op,a,b in reversed(s['objs']):
        require(out in seen,('dead object',out))
        require(a!=b and not isinstance(a,int) and not isinstance(b,int),('degenerate',out))
        seen.update([a,b])
    eq={out for out,op,a,b in s['objs'] if op=='eqxor'}
    paired={x for p in s['pairs'] for x in p}
    for p in s['pars']-paired:
        ch=children(s,p)
        require(not(len(ch)==1 and next(iter(ch)) in eq),('D4 available',p))
    # Affine forms of parameters/equations establish independence at every
    # AND with affine parents, and absence of an unprocessed M gate.
    basis={p:1<<i for i,p in enumerate(sorted(s['pars']))}
    forms=dict(basis)
    for out,op,a,b in s['objs']:
        if op=='eqxor':forms[out]=forms[a]^forms[b]
        elif a in forms and b in forms:
            require(op=='and',('M available',out))
            require(forms[a] and forms[b] and forms[a]!=forms[b],('dependent affine parents',out))
    return True

def q(v):
    x,y,a,b,c,z=[v>>i&1 for i in range(6)]
    return ((x^a)&(x^b))^(y&b)^(a&c)^(x&y)^z

def span_spaces(n):
    spaces={frozenset([0])};frontier=list(spaces)
    while frontier:
        U=frontier.pop()
        for a in range(1<<n):
            if a not in U:
                W=U|frozenset(a^u for u in U)
                if W not in spaces:spaces.add(W);frontier.append(W)
    return spaces

def main():
    rows=[];evaluations=0
    for k in range(1,21):
        s=make(k);initial=stats(s);d=4*k-1
        require(initial==(11*k-1,6*k,0,6*k),'initial inventory')
        closure_checks(s)
        for j in range(k):
            before=stats(s);require(before[3]>=d+2,'loop eligibility')
            require(children(s,f'G{j}')=={f'K{j}'},'one distinct child')
            s['pairs'].add((f'x{j}',f'y{j}'))
            deleted=substitute(s,{f'G{j}':0})
            require(set(deleted)=={f'G{j}',f'K{j}'},('GP exact deletions',deleted))
            after=stats(s)
            require(tuple(b-a for a,b in zip(before,after))==(-2,-2,1,-1),'GP deltas')
            closure_checks(s)
            require(children(s,f'x{j}')=={f'A{j}',f'C{j}'},'P2 exact two children')
            require(children(s,f'y{j}')=={f'B{j}'},'partner survives')
            s['pairs'].remove((f'x{j}',f'y{j}'));s['pars'].remove(f'x{j}')
            deleted=substitute(s,{f'x{j}':0})
            require(set(deleted)=={f'A{j}',f'C{j}'},('P2 exact deletions',deleted))
            final=stats(s)
            require(tuple(b-a for a,b in zip(after,final))==(-2,1,-1,0),'P2 deltas')
            closure_checks(s)
        require(stats(s)==(7*k-1,5*k,0,5*k),'final inventory')
        # Direct semantics on the full surviving cube for k=1,2.
        if k<=2:
            ordered=sorted(s['pars'])
            for bits in product([0,1],repeat=len(ordered)):
                a=dict(zip(ordered,bits));target=0
                for j in range(k):target^=(a[f'a{j}']&a[f'b{j}'])^(a[f'y{j}']&a[f'b{j}'])^(a[f'a{j}']&a[f'c{j}'])^a[f'z{j}']
                require(value(s,a)==target,'final semantic identity');evaluations+=1
        rows.append({'blocks':k,'initial':initial,'final':stats(s),'exact_tight_cycles':k,'max_constant_affine_dimension':d})
    # Independently check GP agreement on all 3-point-pair block inputs.
    orig=make(1);gp=make(1);gp['pairs'].add(('x0','y0'));substitute(gp,{'G0':0})
    for v in range(64):
        a={f'{c}0':v>>i&1 for i,c in enumerate('xyabcz')}
        require(value(orig,a)==q(v),'original quadratic identity')
        if not(a['x0'] and a['y0']):require(value(gp,a)==q(v),'GP domain semantics')
    radical=[u for u in range(64) if all(q(u^v)^q(u)^q(v)^q(0)==0 for v in range(64))]
    require(radical==[0,11,32,43],'polar radical')
    require(q(32)==1,'nonzero radical linear part')
    flats={frozenset(a^u for u in U) for U in span_spaces(6) for a in range(64)}
    dims=[len(F).bit_length()-1 for F in flats if len({q(v) for v in F})==1]
    require(max(dims)==3,'one-block exact affine dispersal dimension')
    print(json.dumps({'status':'PASS','scope':'syntactically reduced, loop-eligible, NOT minimum circuits','cycles':rows,'semantic_evaluations':evaluations,'one_block_affine_flats':len(flats),'radical':radical},indent=2,sort_keys=True))
if __name__=='__main__':main()
