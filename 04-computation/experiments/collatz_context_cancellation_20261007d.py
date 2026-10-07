"""Ordered context cancellation and sixteen-bit paid Collatz deletion.

All production receipts retain actual words, their native source, and an
unsupplied final ROOT obligation. Symbolic phases do not materialize Mersennes.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

import collatz_uncovered_join_routes_20261007 as routes
import collatz_child_normal_forms_20261007 as forms
import collatz_completion_anchor_20261007b as anchors
import collatz_completion_paper24_23_20261007b as budget
import collatz_eight_bit_completion_20261007c as old


@dataclass(frozen=True)
class Macro:
    drop: int
    left: tuple
    right: tuple
    gap: int = 1


CHECKS = 0
DATA=Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_context_cancellation_20261007d.cert.json'


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(x, minimum=0):
    need(type(x) is int and x >= minimum, 'exact integer in domain')


def audit(m):
    need(type(m) is Macro, 'exact Macro')
    integer(m.drop, 1)
    integer(m.gap, 1)
    routes.letters(m.left)
    routes.letters(m.right)
    need(bool(m.left) and bool(m.right), 'nonempty heads')
    need(m.left[0] >= 2 and m.right[0] >= 2, 'maximal initial ones retained')
    need(len(m.right)-len(m.left) == m.drop, 'ordered clock defect')
    need(sum(m.left)-sum(m.right) == 2*m.gap, 'sibling valuation displacement')
    p,q,b = routes.carrier(m.left)
    pp,qq,bb = routes.carrier(m.right)
    need(F(bb-pp,qq) == 4**m.gap*F(b-p,q)+F(4**m.gap-1,3),
         'full affine carry identity')
    return m


def cell(m):
    audit(m)
    p,q,b = routes.carrier(m.left)
    residue = (q-b)*pow(p,-1,2*q)%(2*q)
    need(residue%4 == 1, 'source head follows maximal ones')
    return residue, 2*q


def phase(m):
    residue, modulus = cell(m)
    A = sum(m.left)
    need(A >= 3, 'standard power-of-three order range')
    e, period = anchors.log_three((residue+1)//2,A)
    return e+1, period


def cancel(outer, inner):
    """Compose at an actual middle prefix, retaining its terminal adjustment."""
    audit(outer)
    audit(inner)
    k = len(outer.right)
    need(len(inner.left)>k and inner.left[:k] == outer.right,
         'inner source contains the entire ordered outer partner')
    terminal = inner.left[k]-2*outer.gap
    need(terminal >= 1, 'outer actual terminal remains positive')
    result = Macro(outer.drop+inner.drop,
                   outer.left+(terminal,)+inner.left[k+1:], inner.right, inner.gap)
    audit(result)
    need(sum(result.left) == sum(inner.left), 'middle cancellation preserves native bit cost')
    need(len(result.left) == len(inner.left)-outer.drop, 'source head shrinks by outer deletion')
    return result


def fixed_source(m, run, parameter=0):
    audit(m)
    integer(run, max(2*sum(m.left),m.drop+1))
    integer(parameter)
    residue,modulus = cell(m)
    q = modulus//2
    oddpart = ((residue+1)//2)*pow(3,-run,q)%q
    need(oddpart%2 == 1, 'exact positive odd cofactor')
    return (1 << (run+1))*(oddpart+q*parameter)-1


def receipt(m, source, bit_cap=20000):
    audit(m)
    routes.odd(source)
    integer(bit_cap,1)
    need(source.bit_length() <= bit_cap, 'explicit literal source cap')
    run = forms.valuation(source+1,2)-1
    A=sum(m.left)
    need(run >= max(2*A,m.drop+1), 'uniform all-prefix growth cutoff')
    oddpart=(source+1)>>(run+1)
    residue, modulus=cell(m)
    need((2*pow(3,run,modulus)*oddpart-1)%modulus == residue, 'immutable source guard')
    x=2*3**run*oddpart-1
    p,q,b=routes.carrier(m.left)
    need((p*x+b)%q==0,'native head integrality')
    z=(p*x+b)//q
    need(z%2==1 and z>source, 'native odd endpoint and ROOT safety')
    c=forms.valuation(3*z+1,2)
    child=((source+1)>>m.drop)-1
    return routes.audit(routes.Receipt(source,child,(1,)*run+m.left+(c,),
                        (1,)*(run-m.drop)+m.right+(c+2*m.gap,), (3*z+1)>>c))


def mixed_phase(m, word):
    """Canonical joint exponent phase with an actually affordable inverse word."""
    D=budget.uniform_budget(word).required
    need(D<=m.drop,'uniform budget fits authenticated deletion')
    inverse=forms.mersenne_phase(forms.encode(word))
    need(inverse is not None,'inverse Mersenne phase exists')
    a,M=phase(m)
    b,N=inverse
    d=gcd(M,N)
    need((b-a)%d==0,'native parity compatibility')
    lift=(b-a)//d*pow(M//d,-1,N//d)%(N//d)
    return a+M*lift, M*(N//d), D


def mixed_source(m, word, run=192, parameter=0):
    integer(parameter)
    need(budget.uniform_budget(word).required<=m.drop,'uniform budget fits receipt')
    r,p=forms.native_cell(forms.encode(word))
    base=fixed_source(m,run)
    step=1<<(run+1+sum(m.left))
    lift=(r-base)*pow(step,-1,p)%p
    return fixed_source(m,run,lift+p*parameter)


def transport(m, word, source):
    return budget.consume_deletion(budget.request(source,word),receipt(m,source))


def reset_bridge(source):
    """Inherited reset>=3 rule, specialized to the actual reset-four child."""
    routes.odd(source)
    need(source.bit_length()<=20000,'explicit literal bridge cap')
    run=forms.valuation(source+1,2)-1
    need(run>=1,'initial run for one-bit deletion')
    left=(1,)*run+(4,)
    endpoint=routes.replay(source,left)[-1]
    right=(1,)*(run-1)+(2,2)
    return routes.audit(routes.Receipt(source,(source-1)//2,left,right,endpoint))


def packets():
    # Companion owns the discovered child pair and its reproducible finite search.
    import collatz_eight_child_routes_20261007d as child
    outer=Macro(8,old.LEFT,old.RIGHT)
    inner7=Macro(7,child.LEFT,child.RIGHT7)
    inner8=Macro(8,child.LEFT,child.RIGHT8)
    return audit(outer), audit(inner7), audit(inner8), cancel(outer,inner7), cancel(outer,inner8)


def recursive_packet():
    data=json.loads(DATA.read_text(encoding='utf-8'))
    need(data['format']=='ordered-head-pair-v1','declared recursive data format')
    for field in ('seed','drop','source_depth','source_cost'):
        integer(data[field],1)
    need((data['seed'],data['drop'],data['source_depth'],data['source_cost'])
         ==(924745889,18,258,553),'declared third-stage packet')
    inner=audit(Macro(data['drop'],tuple(data['source_head']),tuple(data['partner_head'])))
    need((len(inner.left),sum(inner.left))==(258,553),'full source head metadata')
    need(phase(inner)==(924745889,2**551),'third-stage native phase')
    whole=cancel(packets()[-1],inner)
    need(phase(whole)==(924745905,2**551),'same original parent refined phase')
    return inner,whole


def recursive_discovery_control():
    """Independent bounded reproduction, not a production certificate oracle."""
    import collatz_eight_child_routes_20261007d as child
    inner,whole=recursive_packet()
    E=924745889
    u,s,bits=child.modular_prefix(E,641,8192)
    need(len(u)==641 and bits>0,'enough source precision')
    hits=[]
    for D in range(1,129):
        v,t,bits=child.modular_prefix(E-D,641,8192)
        need(len(v)==641 and bits>0,'enough partner precision')
        for k in range(513):
            if s[k][1]!=t[k+D][1]: continue
            delta=s[k][0]-t[k+D][0]
            if delta%2==0 and v[k+D]==u[k]+delta:
                hits.append([k,D,delta//2,s[k][0]])
                break
    saved=json.loads(DATA.read_text(encoding='utf-8'))['search']
    need(saved=={'max_depth':512,'max_drop':128,'precision':8192,'hits':hits},
         'entire declared search agrees with frozen discovery data')
    need(hits==[[258,D,1,553] for D in (9,10,13,14,15,16,17,18)],'exact bounded alternatives')
    v,_,_=child.modular_prefix(E-18,277,2048)
    need(u[:258]==inner.left and v[:276]==inner.right,'heads retain their actual source phase')
    return inner,whole


def rejects(function,*args):
    try: function(*args)
    except (ValueError,TypeError): need(True,'hostile rejected')
    else: need(False,'hostile incorrectly accepted')


def main():
    outer,inner7,inner8,whole15,whole16=packets()
    need(phase(whole15)==phase(whole16)==(924745905,2**79),'same refined original parent phase')
    need((len(whole16.left),len(whole16.right),sum(whole16.left),sum(whole16.right))
         ==(31,47,81,79),'compressed sixteen-bit packet')
    need(whole15.left==whole16.left and whole16.right==(2,2)+whole15.right[1:],
         'two distinct paid child choices retain the same source guard')
    for m in (whole15,whole16):
        for run in (170,192,256):
            for t in (0,1,19,10**9):
                n=fixed_source(m,run,t)
                direct=receipt(m,n)
                a=old.receipt(n)
                b=receipt(inner7 if m.drop==15 else inner8,a.child)
                need(old.compose(a,b)==direct,'literal two-stage and cancelled paths agree')
                if m.drop==15:
                    bridge=reset_bridge(direct.child)
                    need(old.compose(direct,bridge)==receipt(whole16,n),
                         'D15 alternative factors through inherited reset bridge to D16')
                need(all(z>n for z in routes.replay(n,direct.source_word[:-1])[1:]),
                     'every original-source prefix before final edge grows')
    for label,limit in ((5,94),(17,168)):
        need(budget.uniform_budget((label,)*limit).required==16,'exact maximum sixteen-bit budget fits')
        need(budget.uniform_budget((label,)*(limit+1)).required==17,'next uniform budget fails')
        w=(label,)*limit
        r,p,d=mixed_phase(whole16,w)
        need((r-924745905)%2**79==0,'joint CRT keeps original source phase')
        inv=forms.mersenne_phase(forms.encode(w))
        need((r-inv[0])%inv[1]==0,'joint CRT keeps inverse native phase')
        for t in (0,1,10**6):
            n=mixed_source(whole16,w,192,t)
            paid=transport(whole16,w,n)
            need(paid.child<paid.source==forms.apply(forms.encode(w),n),'original inverse child is paid')
        print('D16 uniform inverse limit:',label,limit,';joint exponent/period bits',r.bit_length(),p.bit_length())
    need(3**27 < 2**43 and 3**28 > 2**44,'length-only bound is27 generator letters')
    # A repeated type does not restore an old native phase.
    import collatz_general_head_phases_20261007b as heads
    endpoint_exponent=924745905-16
    need(forms.valuation(endpoint_exponent-1,2)==5,'D16 child returns to J3')
    bank=heads.finite_bank()
    need(all((endpoint_exponent-heads.phase(h,3).residue)%min(2**79,heads.phase(h,3).period)
             for h in bank),'returning J3 child misses every225-bank phase')
    need(forms.valuation(924745905-15,2)==1,'D15 child is an even exponent with first reset4')
    for bad in (replace(whole16,drop=True),replace(whole16,gap=1.0),
                replace(whole16,left=tuple(reversed(whole16.left)))):
        rejects(audit,bad)
    rejects(cancel,outer,replace(inner8,left=inner8.left[:15]+(2,)+inner8.left[16:]))
    rejects(receipt,whole16,fixed_source(whole16,170),16)
    rejects(mixed_phase,whole16,(1,))
    rejects(mixed_phase,whole16,(5,)*95)
    rejects(reset_bridge,3)
    deeper,whole34=recursive_discovery_control()
    need((whole34.drop,len(whole34.left),len(whole34.right),sum(whole34.left),sum(whole34.right))
         ==(34,242,276,553,551),'recursive cancellation yields a shorter D34 source head')
    for run in (1122,1200):
        for t in (0,1):
            n=fixed_source(whole34,run,t)
            direct=receipt(whole34,n)
            first=receipt(whole16,n)
            second=receipt(deeper,first.child)
            need(old.compose(first,second)==direct,'recursive direct and composed receipts coincide')
    for label,limit in ((5,200),(17,358)):
        w=(label,)*limit
        need(budget.uniform_budget(w).required==34,'D34 exact uniform budget fits')
        need(budget.uniform_budget(w+(label,)).required==35,'D34 next uniform budget fails')
        n=mixed_source(whole34,w,1200)
        paid=transport(whole34,w,n)
        need(paid.child<paid.source,'third-stage receipt pays original inverse child')
    need(3**58<2**92 and 3**59>2**93,'D34 length-only budget58')
    need(forms.valuation(924745905-34-1,2)==1,'third-stage child has J1, not a completed period')
    import collatz_bott_marked_periodicity_20261007d as marked
    intermediate=cancel(inner8,deeper)
    need(intermediate.left[:len(old.RIGHT)]==old.RIGHT,'same old endpoint coordinate is retained')
    suffix=intermediate.left[len(old.RIGHT):]
    pulled=marked.pullback(suffix)
    need((len(suffix),sum(suffix),pulled.residue,pulled.period,pulled.strict_parameter_cut)
         ==(235,521,0,2**519,0),'recursive guard has an independent lossless parameter readout')
    need(marked.newton_degree(519)==16,'finite-degree reader for519-bit guard')
    for t in (0,1,2,17,2**518,2**519):
        need(marked.polynomial_mod(t,519)==marked.j_mod(t,519),'Newton and power readers agree')
    print('PROVED context cancellation:deletions add, source-head length drops, native bit cost is unchanged.')
    print('PROVED D15/D16 phase K=924745905mod2^79;source head',whole16.left)
    print('FINITE-EXACT24 literal compositions,6 mixed inverse transports;ROOT terminal remains OPEN.')
    print('FINITE-EXACT D16 child returns toJ3 but misses all225 declared head phases.')
    print('PROVED recursive D34 phase K=924745905mod2^551;heads242/276,costs553/551.')
    print('D34 uniform inverse limits G5^200/G17^358;4 literal recursive compositions,2 inverse transports.')
    print('Coverage tradeoff:old parameter t needs47 bits for D16,519 bits for D34;no universal closure.')
    print('Independent D34 guard pullback:t=0mod2^519;exact polynomial readout degree16.')
    print('Exact checks',CHECKS)


if __name__=='__main__':
    main()
