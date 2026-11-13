"""Literal fixed-next-target repairs; finite controls of EXACT_REPAIR.md."""
from pathlib import Path
from itertools import permutations
from math import gcd
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'eight_nine_twentyseven_20261112/constructors/escape.py'
spec = importlib.util.spec_from_file_location('exact_repair_graphs', SOURCE)
graphs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(graphs)
K = (16,9,3,5,14)
U = (18,24,3,5,12)


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def join(left, right):
    n = len(left)
    high = ((1 << len(right)) - 1) << n
    return tuple(row | high for row in left) + tuple(row << n for row in right)


def slots2(t):
    h = [sum(row) for row in graphs.ends(t)]
    h[0] = 1
    return sum(h[mask] * h[(len(h)-1)^mask] for mask in range(len(h)))


def permutation_state(t, root):
    """Independent literal word/cut count: first slot must be nonempty."""
    def path(word):
        return all(t[a] >> b & 1 for a,b in zip(word,word[1:]))
    n=len(t); H=s=F=0
    for word in permutations(range(n)):
        if path(word):
            H += 1
            s += word[0] == root
        cut=word.index(root)+1
        F += path(word[:cut]) and path(word[cut:])
    q=sum(path(word) for word in permutations(i for i in range(n) if i!=root))
    return H,s,q,F


def period(delta):
    g=delta
    while g%2==0:
        g//=2
    if g==1:
        return 1
    value=2%g; k=1
    while value!=1:
        value=value*2%g; k+=1
    return k


def returned(z,k):
    h,s,q,F=z; a=h-s; d=F-q
    t=a+q+((1<<k)-1)*d
    b=a+3*q+(3*(1<<k)-1)*d
    return b,t,t,b


def main():
    need(graphs.flip(graphs.flip(K,1,0),4,1)==U,'literal two-arc switch')
    need(graphs.strong(K) and graphs.strong(U),'base strongness')
    states=[]
    for t,expected in ((K,(15,3,3,21)),(U,(15,3,5,21))):
        z=graphs.root_state(t,0)
        need(z==expected==permutation_state(t,0),'independent base state')
        child,root=graphs.child(t,0)
        states.append(dict(rows=t,state=z,first_child_rows=child,
                           first_child_state=graphs.root_state(child,root)))
    need([r['first_child_state'][0] for r in states]==[33,33],'empty-prefix target')
    left_factors=[graphs.transitive(a) for a in range(1,5)]
    triangle=(2,4,1)
    left_factors += [triangle,join(triangle,graphs.transitive(1)),
                    (6,12,8,1), (24,1,3,6,14)]
    controls=[]
    for L in left_factors:
        h=sum(graphs.ends(L)[-1]); f=slots2(L)
        entries=[]
        for base,q in ((K,3),(U,5)):
            source=join(L,base); root=len(L)
            z=graphs.root_state(source,root)
            expected=(15*h,0,q*h,21*f)
            need(z==expected,'arbitrary left-factor marked law')
            child,r=graphs.child(source,root)
            w=graphs.root_state(child,r)
            need(w==(15*h+21*f,q*h,q*h,42*f-q*h),'same first target')
            need(graphs.strong(child),'first R strongifies source')
            D=3*gcd(h,7*f) if q==3 else gcd(5*h,21*f)
            need(gcd(*z)==gcd(*w)==D,'exact interface contents')
            entries.append(dict(rows=source,root=root,state=z,
                                first_child_state=w,content=D))
        row=dict(left_rows=L,h=h,f=f,entries=entries)
        if gcd(5*h,21*f)==1:
            delta=21*f-20*h
            k=period(delta)
            z=entries[1]['state']
            # Compute the gcd modulo the defect; do not build giant graphs.
            need(gcd(delta,pow(2,k,delta)*(21*f-5*h)+5*h)==1,'period clock primitive')
            need((pow(2,k,3)*h)%3!=0,'controller entry residue')
            row['accepted_period']=k
        controls.append(row)
    literal=[]
    for a in range(1,5):
        x=1<<a; source=join(graphs.transitive(a),U); root=a
        z=graphs.root_state(source,root)
        t=source; r=root
        for _ in range(a):
            t,r=graphs.child(t,r)
        t,r=graphs.child(t,r,True)
        t,r=graphs.child(t,r)
        actual=graphs.root_state(t,r)
        b=63*x*x-36*x+35; q=21*x*x-26*x+25
        need(actual==(b,q,q,b)==returned(z,a),'quadratic dyadic literal return')
        need(gcd(b,q)==1 and graphs.strong(t),'primitive strong return')
        need((b-2*q)%3!=0,'safe C3 entry')
        literal.append(dict(a=a,k=a,order=len(t),rows=t,root=r,state=actual,
                            repaired_target=21*x+15))
    need([v['state'][0] for v in literal]==[215,899,3779,15587],'displayed first counts')
    receipt=dict(status='PASS: finite physical controls; universal claims proved in EXACT_REPAIR.md',
                 universe='Two labelled order5 sources; eight explicitly listed nonempty left factors; four full R^a N R returns with a=1..4, orders9..15. No spectrum census.',
                 base_controls=states,left_factor_controls=controls,literal_returns=literal,
                 inherited_graph_helper=str(SOURCE.relative_to(HERE.parent.parent.parent)),
                 helper_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest())
    (HERE/'exact_repair_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','left_factors':len(controls),
                      'exact_targets':[x['repaired_target'] for x in literal],
                      'primitive_returns':[x['state'][0] for x in literal]}))


if __name__=='__main__':
    main()
