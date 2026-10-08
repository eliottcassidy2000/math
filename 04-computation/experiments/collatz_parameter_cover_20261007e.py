"""Exact dyadic guards proposed from a bounded parameter-complement search."""
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path
import json

import collatz_eight_child_routes_20261007d as reader
import collatz_context_cancellation_20261007d as context
import collatz_uncovered_join_routes_20261007 as routes

BASE = 924745897
STRIDE = 1 << 32
OLD_BITS = 47
DATA = Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_parameter_cover_20261007e.json'
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(n, least=0):
    need(type(n) is int and n >= least, 'exact integer domain')
    return n


@dataclass(frozen=True)
class Rule:
    seed_parameter: int
    drop: int
    gap: int
    source_head: tuple
    partner_head: tuple
    parameter_residue: int
    parameter_bits: int
    minimum_parameter: int


def raw(rule):
    return dict(seed_parameter=rule.seed_parameter, drop=rule.drop, gap=rule.gap,
                source_head=list(rule.source_head), partner_head=list(rule.partner_head))


def artifact_json(payload):
    """Readable records without a separate line for every valuation letter."""
    fields=[]
    for key,value in payload.items():
        if key in ('rules','target_rules'):
            body='[\n'+',\n'.join('    '+json.dumps(r) for r in value)+'\n  ]'
        else:
            body=json.dumps(value)
        fields.append(json.dumps(key)+': '+body)
    return '{\n  '+',\n  '.join(fields)+'\n}\n'


def compile_rule(data):
    need(type(data) is dict and set(data) ==
         {'seed_parameter','drop','gap','source_head','partner_head'}, 'exact retained rule fields')
    for k in ('seed_parameter','drop','gap'):
        natural(data[k], 0 if k == 'seed_parameter' else 1)
    need(type(data['source_head']) is list and type(data['partner_head']) is list,
         'literal retained head arrays')
    u, v = tuple(data['source_head']), tuple(data['partner_head'])
    macro = context.audit(context.Macro(data['drop'], u, v, data['gap']))
    exponent, period = context.phase(macro)
    if period <= STRIDE:
        need((BASE-exponent) % period == 0, 'head meets inherited source phase')
        residue, bits = 0, 0
    else:
        need((exponent-BASE) % STRIDE == 0, 'head meets inherited source phase')
        p = period//STRIDE
        residue, bits = ((exponent-BASE)//STRIDE) % p, p.bit_length()-1
    cutoff = max(0, (max(2*sum(u),data['drop']+1)+1-BASE+STRIDE-1)//STRIDE)
    need(data['seed_parameter'] >= cutoff and
         data['seed_parameter'] % (1 << bits) == residue, 'proposal retains its actual source')
    return Rule(data['seed_parameter'], data['drop'], data['gap'], u, v, residue, bits, cutoff)


def audit_rule(rule):
    need(type(rule) is Rule, 'exact Rule')
    for k in ('seed_parameter','drop','gap','parameter_residue','parameter_bits','minimum_parameter'):
        natural(getattr(rule,k), 1 if k in ('drop','gap') else 0)
    routes.letters(rule.source_head)
    routes.letters(rule.partner_head)
    need(compile_rule(raw(rule)) == rule, 'canonical guard and cutoff recomputed')
    return rule


def contains(rule, parameter):
    audit_rule(rule)
    natural(parameter)
    return parameter >= rule.minimum_parameter and parameter % (1 << rule.parameter_bits) == rule.parameter_residue


def macro(rule):
    audit_rule(rule)
    return context.Macro(rule.drop, rule.source_head, rule.partner_head, rule.gap)


def receipt(rule, source, bit_cap=20000):
    """A native ordinary-source receipt; the source is never replaced."""
    return context.receipt(macro(rule), source, bit_cap)


def symbolic_source(rule, parameter, modulus, child=False):
    need(contains(rule,parameter), 'same guarded nonnegative parameter')
    natural(modulus,1)
    need(type(child) is bool, 'exact source/child flag')
    E = BASE+STRIDE*parameter-(rule.drop if child else 0)
    return (pow(2,E,modulus)-1) % modulus


def head_endpoint_residue(rule, parameter, modulus, child=False):
    need(contains(rule,parameter), 'same guarded parameter')
    natural(modulus,1)
    need(type(child) is bool, 'exact source/child flag')
    E = BASE+STRIDE*parameter-(rule.drop if child else 0)
    w = rule.partner_head if child else rule.source_head
    p,q,b = routes.carrier(w)
    x = (2*pow(3,E-1,q*modulus)-1) % (q*modulus)
    numerator = p*x+b
    need(numerator % q == 0, 'retained denominator precision')
    return numerator//q % modulus


def load_rules():
    data = json.loads(DATA.read_text(encoding='utf-8'))
    need(data['format'] == 'bounded-positive-gap-parameter-cover-v1', 'declared data version')
    need(data['search'] == dict(parameter_limit=1024,max_drop=32,max_depth=128,precision=1024),
         'fixed proposal universe')
    rows = tuple(compile_rule(x) for x in data['rules'])
    need(all(r.minimum_parameter == 0 for r in rows), 'this retained bank needs no finite cuts')
    return rows


def _covered(rows, t):
    return any(t >= r.minimum_parameter and t % (1 << r.parameter_bits) == r.parameter_residue for r in rows)


def antichain(rows):
    """Prefix-free guard union; alternatives inside a cell add no mass."""
    kept = []
    for r in sorted(rows,key=lambda z:(z.parameter_bits,z.parameter_residue,-z.drop,z.seed_parameter)):
        audit_rule(r)
        need(r.minimum_parameter == 0, 'finite-cut rows need a separate point ledger')
        if not any(r.parameter_residue % (1 << x.parameter_bits) == x.parameter_residue for x in kept):
            kept.append(r)
    return tuple(kept)


def mass(rows):
    return sum((Fraction(1,1 << r.parameter_bits) for r in antichain(rows)),Fraction(0))


RUN_PARTNER = (2,2,1,2,2,1,1,2,1,1,2,1,2,2,2,3,2,4,2,1,1)


def parameter_cell(word):
    """Raw native head cell in the retained t coordinate: (residue,bits)."""
    routes.letters(word)
    need(bool(word) and sum(word)>=3,'nonempty native head in the order range')
    p,q,b=routes.carrier(word)
    residue=(q-b)*pow(p,-1,2*q)%(2*q)
    need(residue%4==1,'head follows the maximal initial ones')
    e,period=context.anchors.log_three((residue+1)//2,sum(word))
    exponent=e+1
    if period<=STRIDE:
        need((BASE-exponent)%period==0,'native head intersects source phase')
        return 0,0
    need((exponent-BASE)%STRIDE==0,'native head intersects source phase')
    modulus=period//STRIDE
    return ((exponent-BASE)//STRIDE)%modulus,modulus.bit_length()-1


def extend_ones(base,k):
    """Exact gap-one surgery at a terminal 2; no source is chosen here."""
    context.audit(base)
    natural(k)
    need(base.gap==1 and base.left[-1]==2,'gap-one terminal-two surgery')
    return context.audit(context.Macro(base.drop,base.left[:-1]+(1,)*k+(2,),
                                       base.right+(1,)*k,1))


def run_macro(k):
    natural(k)
    base=context.Macro(4,context.old.RIGHT+(5,2),RUN_PARTNER,1)
    return extend_ones(base,k)


def run_rule(k):
    """Native grammar cell with the generic conservative finite height cut."""
    m=run_macro(k)
    residue,bits=parameter_cell(m.left)
    cutoff=max(0,(max(2*sum(m.left),m.drop+1)+1-BASE+STRIDE-1)//STRIDE)
    period=1<<bits
    seed=residue+max(0,(cutoff-residue+period-1)//period)*period
    return compile_rule(dict(seed_parameter=seed,drop=m.drop,gap=m.gap,
                             source_head=list(m.left),partner_head=list(m.right)))


def run_tail_parent(k):
    """All grammar cells indexed >=k lie in this raw shrinking parent."""
    natural(k)
    return parameter_cell(context.old.RIGHT+(5,)+(1,)*k)


def _cell_antichain(cells):
    kept=[]
    for r,b in sorted(cells,key=lambda c:(c[1],c[0])):
        natural(b)
        natural(r)
        need(r<1<<b,'canonical dyadic residue')
        if not any(r%(1<<d)==s for s,d in kept):
            kept.append((r,b))
    return tuple(kept)


def run_overlap(cells):
    """Exact asymptotic overlap with any finite dyadic union, (residue,bits)."""
    cover=_cell_antichain(cells)
    if not cover:
        return Fraction(0)
    bound=max(b for _,b in cover)
    stop=max(0,bound-3)
    total=Fraction(0)
    for k in range(stop):
        r,b=parameter_cell(run_macro(k).left)
        for s,d in cover:
            if r%(1<<min(b,d))==s%(1<<min(b,d)):
                total+=Fraction(1,1<<max(b,d))
    r,b=run_tail_parent(stop)
    need(b>=bound,'tail parent decides every finite-bank cell')
    if any(r%(1<<d)==s for s,d in cover):
        total+=Fraction(1,1<<(stop+4))
    return total


def discover():
    """Explicit bounded proposals, retaining all successful authenticated rows."""
    rows, attempted, misses = [], [], []
    for t in range(1024):
        if _covered(rows,t):
            continue
        attempted.append(t)
        E = BASE+STRIDE*t
        u,s,_ = reader.modular_prefix(E,161,1024)
        need(len(u) == 161, 'proposal source has enough bits')
        options = []
        for D in range(1,33):
            v,q,_ = reader.modular_prefix(E-D,161,1024)
            need(len(v) == 161, 'proposal partner has enough bits')
            for k in range(129):
                if s[k][1] != q[k+D][1]:
                    continue
                delta = s[k][0]-q[k+D][0]
                if delta > 0 and delta % 2 == 0 and v[k+D] == u[k]+delta:
                    options.append((s[k][0],k,-D,-delta//2,u[:k],v[:k+D]))
                    break
        if not options:
            misses.append(t)
            continue
        _,_,negativeD,negativeGap,left,right = min(options)
        rows.append(compile_rule(dict(seed_parameter=t,drop=-negativeD,gap=-negativeGap,
                                      source_head=list(left),partner_head=list(right))))
    return rows,attempted,misses


def reject(function,*args):
    try:
        function(*args)
    except (ValueError,TypeError):
        need(True,'hostile rejected')
        return
    need(False,'hostile accepted')


def grammar_controls(rows):
    need(3**78>2**117,'uniform grammar growth at e78')
    cells=[]
    for k in range(33):
        r=run_rule(k)
        need(r.drop==4 and r.gap==1 and r.parameter_bits==k+5,'all-length grammar parameters')
        need(sum(r.source_head)==39+k and sum(r.partner_head)==37+k,'grammar exact costs')
        need(r.minimum_parameter==0,'bounded controls need no finite cuts')
        cell=(r.parameter_residue,r.parameter_bits)
        cells.append(cell)
        for previous in cells[:-1]:
            a,b=previous
            need(a%(1<<b)!=cell[0]%(1<<b),'run cells are disjoint')
        parent,pbits=run_tail_parent(k)
        need(pbits==k+3 and cell[0]%(1<<pbits)==parent,'same shrinking tail parent')
        for j in range(k):
            a,b=cells[j]
            need(a%(1<<min(b,pbits))!=parent%(1<<min(b,pbits)),
                 'earlier exits avoid the tail parent')
        E=BASE+STRIDE*r.parameter_residue
        actual,_,_=reader.modular_prefix(E,len(r.source_head)+1,sum(r.source_head)+64)
        need(actual[:len(r.source_head)]==r.source_head,'independent grammar native replay')
        if k<9:
            source=context.fixed_source(macro(r),2*sum(r.source_head))
            got=receipt(r,source)
            need(got.child<got.source,'grammar ordinary-source common future')
    need(sum((Fraction(1,1<<b) for _,b in cells),Fraction())==
         Fraction(1,16)-Fraction(1,1<<37),'finite grammar geometric sum')
    matches=[]
    for r in rows:
        need(r.source_head[:15]==context.old.RIGHT,'finite bank has common inherited prefix')
        tail=r.source_head[15:]
        if not tail or tail[0]!=5:
            continue
        k=0
        while k+1<len(tail) and tail[k+1]==1:
            k+=1
        need(k+1<len(tail),'no finite rule stops inside the grammar parent')
        if tail[k+1]==2:
            need(tail==(5,)+(1,)*k+(2,),'grammar intersections are exactly terminal exits')
            matches.append(k)
    need(sorted(matches)==list(range(6)),'all finite-bank grammar intersections enumerated')
    overlap=run_overlap([(r.parameter_residue,r.parameter_bits) for r in rows])
    need(overlap==Fraction(63,1024),'independent dyadic overlap agrees with prefix classification')
    need(run_overlap([(5,3)])==Fraction(1,16),'whole grammar in its source parent')
    need(run_overlap([(4,9)])==0,'remaining t4mod512 misses grammar')
    need(not any(4%(1<<min(9,r.parameter_bits))==
                 r.parameter_residue%(1<<min(9,r.parameter_bits)) for r in rows),
         'whole t4mod512 cylinder remains a finite-bank defect')
    for x in (Fraction(-7,3),Fraction(-1),Fraction(0),Fraction(13,7),Fraction(100)):
        f1=lambda z:(3*z+1)/2
        f2=lambda z:(3*z+1)/4
        sibling=lambda z:4*z+1
        need(f1(sibling(f2(x)))==sibling(f2(f1(x))),'rational affine commutation control')
    reject(run_rule,True)
    reject(run_rule,-1)
    reject(run_tail_parent,1.0)
    reject(run_overlap,[(True,3)])
    reject(extend_ones,context.Macro(4,context.old.RIGHT+(6,),
        (2,2,1,2,2,1,1,2,1,1,2,1,2,2,2,3,2,4,2,1),1),1)
    return overlap


def main(rediscover=False):
    rows = load_rules()
    cover = antichain(rows)
    total = mass(rows)
    need(_covered(cover,0), 'retained old parameter point')
    old_covered = any(OLD_BITS >= r.parameter_bits and r.parameter_residue == 0 for r in cover)
    need(old_covered,'old entire0mod2^47 guard included')
    local = [t for t in range(1024) if _covered(cover,t)]
    data = json.loads(DATA.read_text(encoding='utf-8'))
    need(data['covered_sample_count'] == len(local), 'sample statistic kept separate from exact mass')
    need(data['union_mass'] == [total.numerator,total.denominator], 'exact stored dyadic union')
    if rediscover:
        found,attempted,misses = discover()
        need([raw(r) for r in found] == data['rules'], 'bounded deterministic proposal reproduction')
        need(attempted == data['attempted'] and misses == data['misses'], 'honest skip/failure ledger')
    controls=0
    for r in rows:
        E = BASE+STRIDE*r.seed_parameter
        u,_,_=reader.modular_prefix(E,len(r.source_head)+1,sum(r.source_head)+64)
        v,_,_=reader.modular_prefix(E-r.drop,len(r.partner_head)+1,sum(r.partner_head)+64)
        need(u[:len(r.source_head)]==r.source_head and v[:len(r.partner_head)]==r.partner_head,
             'retained source-authentic words')
        need(v[len(r.partner_head)]==u[len(r.source_head)]+2*r.gap,'actual terminal reserve')
        n=context.fixed_source(macro(r),max(2*sum(r.source_head),r.drop+1),0)
        got=receipt(r,n)
        need(got.child<got.source,'ordinary control pays its immutable source')
        for t in (r.parameter_residue,r.parameter_residue+(1<<r.parameter_bits)*10**20):
            for m in (19,1<<128,3**12):
                a=head_endpoint_residue(r,t,m)
                b=head_endpoint_residue(r,t,m,True)
                need(b==(4**r.gap*a+(4**r.gap-1)//3)%m,'all-height symbolic carrier control')
        controls+=1
    first=next(r for r in rows if r.seed_parameter==1)
    need((first.parameter_residue,first.parameter_bits)==(1,4),'coarse complement cell1mod16')
    need(total>Fraction(1,16),'proved mass exceeds the initial coarse gain')
    reject(audit_rule,replace(first,parameter_bits=True))
    reject(audit_rule,replace(first,drop=float(first.drop)))
    reject(audit_rule,replace(first,parameter_residue=0))
    reject(symbolic_source,first,0,19)
    reject(symbolic_source,first,1,True)
    reject(contains,first,1.0)
    overlap=grammar_controls(rows)
    combined=total+Fraction(1,16)-overlap
    print('PROVED exact finite dyadic union:',total)
    print('Decimal mass (display only):',float(total))
    print('Increment beyond old0mod2^47:',total-Fraction(1,1<<47))
    print('Retained rules/prefix-free guards:',len(rows),len(cover))
    print('FINITE-EXACT sample covered/1024:',len(local),'(not a density estimate)')
    print('FINITE-EXACT proposal attempts/misses:',len(data['attempted']),len(data['misses']))
    print('Independent modest-source and symbolic rule controls:',controls)
    print('PROVED infinite run grammar mass:',Fraction(1,16))
    print('Exact grammar overlap / increment:',overlap,Fraction(1,16)-overlap)
    print('PROVED finite bank plus run grammar:',combined)
    print('Combined decimal mass (display only):',float(combined))
    print('Residual whole cell: t=4mod512; grammar generic API retains conservative finite cuts.')
    print('Coarse cell:t1mod16; ROOT children remain marked obligations.')
    print('No universal coverage, no sampled-frequency extrapolation, no enormous integer expansion.')
    print('Checks:',CHECKS+reader.CHECKS+context.CHECKS+routes.CHECKS)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--discover',action='store_true',help='write the explicitly bounded proposal artifact')
    p.add_argument('--rediscover',action='store_true',help='compare bounded discovery with frozen data')
    args=p.parse_args()
    if args.discover:
        rows,attempted,misses=discover()
        total=mass(rows)
        payload=dict(format='bounded-positive-gap-parameter-cover-v1',
            search=dict(parameter_limit=1024,max_drop=32,max_depth=128,precision=1024),
            rules=[raw(r) for r in rows],attempted=attempted,misses=misses,
            union_mass=[total.numerator,total.denominator],
            covered_sample_count=sum(_covered(rows,t) for t in range(1024)))
        DATA.write_text(artifact_json(payload),encoding='utf-8')
        print('Frozen proposals',len(rows),'guards',len(antichain(rows)),'mass',total,'sample',payload['covered_sample_count'])
    else:
        main(args.rediscover)
