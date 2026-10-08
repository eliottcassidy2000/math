"""Source-side sibling guards filling the retained positive-gap complement.

This is a paid common-future certificate compiler, not a ROOT finder.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import sys

import collatz_parameter_cover_20261007e as positive
import collatz_eight_child_routes_20261007d as reader
import collatz_complement_ladders_20261007c as signed
import collatz_completion_anchor_20261007b as anchors
import collatz_uncovered_join_routes_20261007 as routes

DATA = Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_signed_parameter_fill_20261007e.json'
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, minimum=0):
    need(type(n) is int and n >= minimum, 'exact integer domain')


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


@dataclass(frozen=True)
class Rule:
    seed: int
    drop: int
    gap: int
    left: tuple
    right: tuple
    residue: int
    bits: int
    minimum: int


def packet(r):
    return signed.Ladder(signed.Relation(r.drop, 3**r.drop-1), r.left, r.right, r.gap)


def raw(r):
    return dict(seed=r.seed, drop=r.drop, gap=r.gap, left=list(r.left), right=list(r.right))


def compile_rule(data):
    need(type(data) is dict and set(data) == {'seed','drop','gap','left','right'}, 'retained rule fields')
    integer(data['seed']); integer(data['drop'],1)
    need(type(data['gap']) is int and data['gap'] < 0, 'negative signed gap')
    need(type(data['left']) is list and type(data['right']) is list, 'literal word lists')
    left,right = tuple(data['left']),tuple(data['right'])
    p = signed.audit(signed.Ladder(signed.Relation(data['drop'],3**data['drop']-1),left,right,data['gap']))
    need(bool(left) and bool(right) and min(left[0],right[0]) >= 2, 'maximal initial runs retained')
    residue,modulus,effective = signed.coarse_cell(p)
    e,period = anchors.log_three((residue+1)//2,effective)
    need(period >= positive.STRIDE and (e+1-positive.BASE) % positive.STRIDE == 0, 'same exponent chart')
    width = period//positive.STRIDE
    r = ((e+1-positive.BASE)//positive.STRIDE) % width
    # Separate sufficient growth bounds for both source and child prefixes.
    run_floor = max(2*sum(left),data['drop']+2*sum(right),data['drop']+1)
    cut = max(0,(run_floor+1-positive.BASE+positive.STRIDE-1)//positive.STRIDE)
    need(data['seed'] >= cut and data['seed'] % width == r, 'retained actual proposal')
    return Rule(data['seed'],data['drop'],data['gap'],left,right,r,width.bit_length()-1,cut)


def audit(r):
    need(type(r) is Rule, 'exact Rule')
    for field in ('seed','drop','residue','bits','minimum'):
        integer(getattr(r,field),1 if field == 'drop' else 0)
    signed.old.word_type(r.left); signed.old.word_type(r.right)
    need(r == compile_rule(raw(r)), 'recomputed signed guard and growth cut')
    return r


def contains(r,t):
    audit(r); integer(t)
    return t >= r.minimum and t % (1 << r.bits) == r.residue


def fixed_source(r, run, cofactor=0, reserve=True):
    audit(r); integer(cofactor)
    need(type(reserve) is bool,'exact terminal-reserve flag')
    integer(run,max(2*sum(r.left),r.drop+2*sum(r.right),r.drop+1))
    if reserve:
        residue,modulus,_ = signed.coarse_cell(packet(r))
    else:
        p,q,b=routes.carrier(r.left)
        modulus=2*q
        residue=(q-b)*pow(p,-1,modulus)%modulus
    q=modulus//2
    oddpart=((residue+1)//2)*pow(3,-run,q)%q
    need(oddpart%2 == 1,'native positive odd cofactor')
    return (1 << (run+1))*(oddpart+q*cofactor)-1


def receipt(r,source,bit_cap=20000):
    audit(r); routes.odd(source); integer(bit_cap,1)
    need(source.bit_length() <= bit_cap,'literal source cap')
    run=v2(source+1)-1
    integer(run,max(2*sum(r.left),r.drop+2*sum(r.right),r.drop+1))
    oddpart=(source+1)>>(run+1)
    residue,modulus,_=signed.coarse_cell(packet(r))
    x=2*3**run*oddpart-1
    need(x%modulus == residue,'native source and terminal reserve')
    p,q,b=routes.carrier(r.left)
    need((p*x+b)%q == 0,'exact source endpoint')
    z=(p*x+b)//q
    c=v2(3*z+1)
    need(c+2*r.gap >= 1,'positive adjusted child terminal')
    child=((source+1)>>r.drop)-1
    return routes.audit(routes.Receipt(source,child,(1,)*run+r.left+(c,),
                        (1,)*(run-r.drop)+r.right+(c+2*r.gap,), (3*z+1)>>c))


def _covered(rows,t):
    return any(t >= r.minimum and t % (1 << r.bits) == r.residue for r in rows)


def discover():
    """Same bounded proposal box as positive bank; only omitted signed side."""
    prior=positive.load_rules()
    result=[]
    for t in range(1024):
        if positive._covered(prior,t):
            continue
        E=positive.BASE+positive.STRIDE*t
        u,s,_=reader.modular_prefix(E,161,1024)
        need(len(u)==161,'source precision sufficient')
        options=[]
        for d in range(1,33):
            v,q,_=reader.modular_prefix(E-d,161,1024)
            need(len(v)==161,'child precision sufficient')
            for k in range(129):
                delta=s[k][0]-q[k+d][0]
                if delta < 0 and delta%2 == 0 and s[k][1]==q[k+d][1] and v[k+d]==u[k]+delta:
                    options.append((s[k][0]-delta,k,-d,delta//2,u[:k],v[:k+d]))
                    break
        if options:
            _,_,negativeD,gap,left,right=min(options)
            result.append(compile_rule(dict(seed=t,drop=-negativeD,gap=gap,left=list(left),right=list(right))))
    return result


def targeted_six():
    """Freeze a deeper receipt at the first residual after the ternary bank."""
    t,d,k=6,4,236
    E=positive.BASE+positive.STRIDE*t
    u,s,_=reader.modular_prefix(E,577,4096)
    v,q,_=reader.modular_prefix(E-d,577,4096)
    need(len(u)==len(v)==577,'deeper reader precision')
    gap=(s[k][0]-q[k+d][0])//2
    return compile_rule(dict(seed=t,drop=d,gap=gap,left=list(u[:k]),right=list(v[:k+d])))


def load():
    data=json.loads(DATA.read_text(encoding='utf-8'))
    need(data['format']=='signed-parameter-fill-v1','data version')
    rows=tuple(compile_rule(x) for x in data['rules'])
    return rows,compile_rule(data['targeted_six'])


def antichain(cells):
    keep=[]
    for residue,bits in sorted(set(cells),key=lambda c:(c[1],c[0])):
        if not any(bits>=b and residue%(1 << b)==r for r,b in keep):
            keep.append((residue,bits))
    return tuple(keep)


def mass(cells):
    return sum((F(1,1 << b) for _,b in antichain(cells)),F())


def word_guard(word):
    """The native source head alone, without its terminal reserve."""
    p,q,b=routes.carrier(word)
    residue=(q-b)*pow(p,-1,2*q)%(2*q)
    e,period=anchors.log_three((residue+1)//2,sum(word))
    need(period>=positive.STRIDE and (e+1-positive.BASE)%positive.STRIDE==0,'same chart for run prefix')
    return ((e+1-positive.BASE)//positive.STRIDE)%(period//positive.STRIDE), (period//positive.STRIDE).bit_length()-1


@lru_cache(maxsize=1)
def mirrored_base():
    rows,_=load()
    return next(r for r in rows if r.seed==9)


def mirrored_run(k):
    """Gap-minus-one mirror of the positive last-two run surgery."""
    integer(k)
    base=mirrored_base()
    need(base.gap==-1 and base.right[-1]==2,'mirrored operation boundary')
    left=base.left+(1,)*k
    right=base.right[:-1]+(1,)*k+(2,)
    p=signed.audit(signed.Ladder(packet(base).relation,left,right,-1))
    r,m,cost=signed.coarse_cell(p)
    e,period=anchors.log_three((r+1)//2,cost)
    seed=((e+1-positive.BASE)//positive.STRIDE)%(period//positive.STRIDE)
    floor=max(2*sum(left),base.drop+2*sum(right),base.drop+1)
    cut=max(0,(floor+1-positive.BASE+positive.STRIDE-1)//positive.STRIDE)
    seed += max(0,(cut-seed+period//positive.STRIDE-1)//(period//positive.STRIDE))*(period//positive.STRIDE)
    return compile_rule(dict(seed=seed,drop=base.drop,gap=-1,left=list(left),right=list(right)))


def mirrored_increment(cells):
    """Exact extra density of all run lengths over a finite dyadic union.

    A remaining prefix has either no overlap, total containment, or a finite
    next decision. At max input precision every overlap is containment.
    """
    cells=antichain(cells)
    base=mirrored_run(0)
    increment=F()
    for k in range(max((b for _,b in cells),default=0)+1):
        r=mirrored_run(k)
        need(r.bits==14+k,'one extra address bit per run letter')
        cell=(r.residue,r.bits)
        increment+=mass(list(cells)+[cell])-mass(cells)
        tail=word_guard(base.left+(1,)*(k+1))
        # Every future completed run lies in this uncompleted prefix.
        overlaps=[(a,b) for a,b in cells if (tail[0]-a)%(1 << min(tail[1],b))==0]
        if not overlaps:
            return increment+F(1,1 << (14+k)),k+1,'disjoint future prefix'
        if any(b<=tail[1] for _,b in overlaps):
            return increment,k+1,'contained future prefix'
    raise ArithmeticError('finite-cell refinement failed to terminate')


def reject(fn,*args):
    try:
        fn(*args)
    except (TypeError,ValueError):
        need(True,'hostile rejected')
        return
    need(False,'hostile accepted')


def main(rediscover=False):
    rows,six=load()
    prior=positive.load_rules()
    need(len(rows)==136,'declared bounded signed proposals')
    if rediscover:
        need(tuple(discover())==rows,'deterministic discovery reproduction')
        need(targeted_six()==six,'deeper point-six reproduction')
    allrows=rows+(six,)
    pcells=[(r.parameter_residue,r.parameter_bits) for r in prior]
    scells=[(r.residue,r.bits) for r in allrows]
    gain=mass(pcells+scells)-mass(pcells)
    need(gain>0 and all(r.minimum==0 for r in allrows),'new positive mass, no finite cuts')
    for r in allrows:
        need(contains(r,r.seed),'same proposal')
        reserve=-2*r.gap
        effective=sum(r.left)+reserve
        modulus=1 << (effective+1)
        native,_,_=signed.coarse_cell(packet(r))
        for shift in (0,1,3):
            t=r.residue+(shift << r.bits)
            E=positive.BASE+positive.STRIDE*t
            x=(2*pow(3,E-1,modulus)-1)%modulus
            need(x==native,'independent power-modulus guard')
            u,states,_=reader.modular_prefix(E,len(r.left)+1,effective+32)
            need(u[:len(r.left)]==r.left and u[-1]>=reserve+1,'actual source word and reserve')
            v,_,_=reader.modular_prefix(E-r.drop,len(r.right)+1,effective+32)
            need(v[:len(r.right)]==r.right and v[-1]==u[-1]+2*r.gap,'actual child word and terminal')
        run=max(2*sum(r.left),r.drop+2*sum(r.right),r.drop+1)
        for k in (0,1):
            n=fixed_source(r,run,k)
            rec=receipt(r,n)
            need(rec.child<n and rec.source==n,'literal paid original source')
    nine=next(r for r in rows if r.seed==9)
    three93=next(r for r in rows if r.seed==393)
    need((nine.residue,nine.bits,nine.drop,nine.gap)==(9,14,4,-1),'short t9 fill')
    need((three93.residue,three93.bits)==(393,10),'coarse t393 fill')
    need(contains(six,6) and not positive._covered(prior,6),'deeper missing source covered')
    reject(audit,replace(nine,gap=True))
    reject(audit,replace(nine,bits=nine.bits-1))
    reject(audit,replace(nine,right=nine.right[:-1]+(nine.right[-1]+1,)))
    reject(contains,nine,True)
    reject(fixed_source,nine,100,0,1)
    # Dropping reserve admits actual native source heads with an invalid child terminal.
    run=max(2*sum(nine.left),nine.drop+2*sum(nine.right),nine.drop+1)
    hostile=None
    for k in range(4):
        n=fixed_source(nine,run,k,reserve=False)
        x=2*3**run*((n+1)>>(run+1))-1
        p,q,b=routes.carrier(nine.left)
        c=v2(3*((p*x+b)//q)+1)
        if c+2*nine.gap<=0:
            hostile=(n,c)
            reject(receipt,nine,n)
            break
    need(hostile is not None,'actual missing-reserve hostile')
    for k in range(16):
        r=mirrored_run(k)
        need(r.bits==14+k and r.drop==4 and r.gap==-1,'all-k mirrored parameter precision')
        for j in range(k):
            s=mirrored_run(j)
            need((r.residue-s.residue)%(1 << min(r.bits,s.bits))!=0,'disjoint completed run lengths')
        run=max(2*sum(r.left),r.drop+2*sum(r.right),r.drop+1)
        rec=receipt(r,fixed_source(r,run))
        need(rec.child<rec.source,'mirrored literal paid receipt')
    mirrored_gain,stop,reason=mirrored_increment(pcells+scells)
    print('PROVED: negative sibling gaps retain terminal reserve and smaller-child common futures.')
    print('FINITE-EXACT: positive-bank misses in t0..1023; D1..32; depth0..128; precision1024.')
    print('Signed proposed rows',len(rows),'prefix-free signed cells',len(antichain(scells)))
    print('Added signed mass',str(gain),'decimal',format(float(gain),'.15f'))
    print('Finite positive+signed mass',str(mass(pcells+scells)),'decimal',format(float(mass(pcells+scells)),'.15f'))
    print('t9 row',raw(nine),'guard',(nine.residue,nine.bits))
    print('t393 guard',(three93.residue,three93.bits),'drop',three93.drop,'gap',three93.gap)
    print('t6 deeper guard',(six.residue,six.bits),'drop',six.drop,'gap',six.gap,'depth',len(six.left))
    print('Mirrored all-run family mass',str(F(1,1 << 13)),'increment over finite banks',str(mirrored_gain))
    print('Mirrored finite overlap cutoff',stop,reason)
    print('Missing-reserve hostile source',hostile[0],'actual terminal',hostile[1])
    print('Checks',CHECKS)
    print('OPEN: uncovered parameters; emitted-child ROOT proofs; universal Collatz.')


if __name__=='__main__':
    if '--discover' in sys.argv:
        data={'format':'signed-parameter-fill-v1','rules':[raw(r) for r in discover()],
              'targeted_six':raw(targeted_six())}
        DATA.write_text(positive.artifact_json(data),encoding='utf-8')
    main('--rediscover' in sys.argv)
