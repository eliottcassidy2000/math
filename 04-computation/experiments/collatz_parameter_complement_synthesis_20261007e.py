"""Exact mixed coverage, four further paid targets, and a scoped residual."""
from pathlib import Path
from fractions import Fraction as F
import json
import sys

import collatz_parameter_cover_20261007e as positive
import collatz_signed_parameter_fill_20261007e as signed
import collatz_complement_guard_fusion_20261007e as ternary
import collatz_bounded_inverse_cover_20261007e as inverse
import collatz_eight_child_routes_20261007d as reader

DATA=Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_parameter_complement_synthesis_20261007e.json'
CHECKS=0


def need(ok,reason):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(reason)


def discover():
    targets=[]
    absent=[]
    for t in (12,23,24,27,32):
        E=positive.BASE+positive.STRIDE*t
        u,s,_=reader.modular_prefix(E,1153,8192)
        need(len(u)==1153,'source comparison precision')
        options=[]
        for d in range(1,129):
            v,q,_=reader.modular_prefix(E-d,1153,8192)
            need(len(v)==1153,'partner comparison precision')
            for k in range(1025):
                delta=s[k][0]-q[k+d][0]
                if delta%2==0 and s[k][1]==q[k+d][1] and v[k+d]==u[k]+delta:
                    options.append((s[k][0]+max(0,-delta),k,-d,delta//2,u[:k],v[:k+d]))
                    break
        if options:
            _,_,nd,gap,left,right=min(options)
            need(gap>0,'these four chosen targets have positive gaps')
            targets.append(positive.compile_rule(dict(seed_parameter=t,drop=-nd,gap=gap,
                              source_head=list(left),partner_head=list(right))))
        else:absent.append(t)
    return tuple(targets),tuple(absent)


def positive_run_increment(cells):
    """Adaptive exact overlap with the infinite grammar; no large period array."""
    cells=signed.antichain(cells)
    total=F()
    for k in range(max((b for _,b in cells),default=0)+1):
        rule=positive.run_rule(k)
        c=(rule.parameter_residue,rule.parameter_bits)
        total+=signed.mass(list(cells)+[c])-signed.mass(cells)
        r,b=positive.run_tail_parent(k+1)
        overlap=[(a,d) for a,d in cells if (r-a)%(1 << min(b,d))==0]
        if not overlap:return total+F(1,1 << (k+5)),k+1,'disjoint tail'
        if any(d<=b for _,d in overlap):return total,k+1,'contained tail'
    raise ArithmeticError('finite bank failed to decide shrinking tail')


def main(rediscover=False):
    data=json.loads(DATA.read_text(encoding='utf-8'))
    need(data['format']=='mixed-parameter-complement-v1','data version')
    targets=tuple(positive.compile_rule(r) for r in data['target_rules'])
    need(data['absent_targets']==[23],'scoped no-hit target')
    if rediscover:
        new,absent=discover()
        need(new==targets and list(absent)==data['absent_targets'],'complete five-target proposal reproduction')
    pos=positive.load_rules()
    neg,six=signed.load()
    neg=neg+(six,)
    three=ternary.entries(12)
    ternary_cells=inverse.combined_cells()
    cells=[(r.parameter_residue,r.parameter_bits) for r in pos+targets]
    cells += [(r.residue,r.bits) for r in neg]
    finite=signed.mass(cells)
    run_gain,stop,reason=positive_run_increment(cells)
    mirror_gain,mstop,mreason=signed.mirrored_increment(cells)
    # Different first continuation letters, and independently incompatible cells.
    a,b=positive.run_tail_parent(0)
    c,d=signed.word_guard(signed.mirrored_base().left)
    need((a-c)%(1 << min(b,d))!=0,'two countable families are disjoint')
    beta=finite+run_gain+mirror_gain
    tau=inverse.cell_mass(ternary_cells)
    need(tau==F(150929272,387420489),'complete short-word bank plus reset-indexed entries')
    joint=1-(1-beta)*(1-tau)
    old_tau=F(244,729)
    old_joint=1-(1-F(1,1 << 47))*(1-old_tau)
    need(run_gain==F(1,1024) and mirror_gain==F(1,16384),'retained distinct all-run gains')
    # Compare an independent closed tail-overlap algorithm on a bounded bank.
    compact=[(r.parameter_residue,r.parameter_bits) for r in pos if r.parameter_bits<=16]
    need(positive_run_increment(compact)[0]==F(1,16)-positive.run_overlap(compact),'independent grammar overlap')
    control_rows=[]
    for r in targets:
        need(r.minimum_parameter==0,'target has no finite source-height cut')
        macro=positive.macro(r)
        run=max(2*sum(r.source_head),r.drop+1)
        n=positive.context.fixed_source(macro,run)
        rec=positive.receipt(r,n)
        need(rec.source==n and rec.child<n,'target literal source-owned paid receipt')
        for j in (0,1):
            t=r.parameter_residue+(j << r.parameter_bits)
            E=positive.BASE+positive.STRIDE*t
            u,_,_=reader.modular_prefix(E,len(r.source_head)+1,2*sum(r.source_head)+64)
            v,_,_=reader.modular_prefix(E-r.drop,len(r.partner_head)+1,2*sum(r.source_head)+64)
            need(u[:-1]==r.source_head and v[:-1]==r.partner_head,'modular actual target heads')
            need(v[-1]==u[-1]+2*r.gap,'actual target terminal adjustment')
        control_rows.append((r.seed_parameter,r.drop,len(r.source_head),r.parameter_bits))
    need(control_rows==[(12,8,217,379),(24,6,136,260),(27,50,651,1348),(32,24,686,1328)],'exact target scopes')
    def covered(t):
        return any(t%(1 << bits)==residue for residue,bits in cells) or any(
            t%period==residue for residue,period in ternary_cells)
    first=next(t for t in range(64) if not covered(t))
    need(first==23 and all(covered(t) for t in range(23)),'first finite-bank+ternary defect')
    need(first%(1 << b)!=a and first%(1 << d)!=c,'neither infinite grammar reaches first defect')
    def valuation(n,prime):
        n=abs(n)
        need(n>0,'different defect address')
        result=0
        while n%prime==0:
            n//=prime
            result+=1
        return result
    binary_separation=max(valuation(first-r,2)+1 for r,_ in cells+[(a,b),(c,d)])
    ternary_separation=max(valuation(first-r,3)+1 for r,_ in ternary_cells)
    defect_period=(1 << binary_separation)*3**ternary_separation
    need((binary_separation,ternary_separation,defect_period)==(9,3,13824),'complete mixed defect cylinder')
    for r,bits in cells+[(a,b),(c,d)]:
        need((first-r)%(1 << min(binary_separation,bits))!=0,'whole defect avoids each binary guard')
    for r,period in ternary_cells:
        m=min(3**ternary_separation,period)
        need((first-r)%m!=0,'whole defect avoids each ternary guard')
    # Independent exact CRT intersections for a compact cross-bank selection.
    for r in pos[:8]+targets:
        for e in three:
            x,m=ternary.fuse(e,r.parameter_residue,r.parameter_bits)
            need(m==(1 << r.parameter_bits)*e.parameter_period and
                 x%(1 << r.parameter_bits)==r.parameter_residue and
                 x%e.parameter_period==e.parameter_residue,'same-parameter CRT witness')
    results=dict(finite_binary=[finite.numerator,finite.denominator],
                 full_binary=[beta.numerator,beta.denominator],
                 ternary=[tau.numerator,tau.denominator],
                 combined=[joint.numerator,joint.denominator],
                 inherited_combined=[old_joint.numerator,old_joint.denominator],
                 least_uncovered=first,uncovered_cylinder=[first,defect_period])
    if '--record' in sys.argv:
        data['results']=results
        DATA.write_text(positive.artifact_json(data),encoding='utf-8')
    else:need(data['results']==results,'exact saved combined density ledger')
    print('PROVED: source-owned binary/ternary union, including two shrinking-tail grammars.')
    print('Further targets (t,drop,depth,bits)',control_rows)
    print('Positive grammar added mass',run_gain,'after',stop,reason)
    print('Mirrored grammar added mass',mirror_gain,'after',mstop,mreason)
    for name,value in (('finite binary',finite),('full binary',beta),('ternary',tau),
                       ('combined',joint),('inherited combined',old_joint),('gain',joint-old_joint)):
        print(name,format(float(value),'.15f'))
    print('Least uncovered t',first,'exponent',positive.BASE+positive.STRIDE*first)
    print('Entire uncovered cylinder: t=23 modulo13824 (binary precision9, ternary precision3).')
    print('No collision at t23 in the declared D1..128, source depth0..1024, precision8192 box.')
    print('Checks',CHECKS)
    print('OPEN: larger templates at t23, all residual parameters, and emitted child ROOT obligations.')


if __name__=='__main__':
    if '--discover' in sys.argv:
        rows,absent=discover()
        DATA.write_text(positive.artifact_json(dict(format='mixed-parameter-complement-v1',
                        target_rules=[positive.raw(r) for r in rows],absent_targets=list(absent))),encoding='utf-8')
    main('--rediscover' in sys.argv)
