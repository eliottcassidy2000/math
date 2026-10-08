"""Source-owned ternary inverse entries and exact binary/ternary fusion.

The endpoint inverse obstruction is separate from the paid source entries.
Production never expands the astronomical Mersenne or searches for ROOT.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from math import gcd

import collatz_bott_marked_periodicity_20261007d as endpoint
import collatz_child_normal_forms_20261007 as forms
import collatz_uncovered_join_routes_20261007 as routes

BASE = endpoint.BASE
STRIDE = endpoint.STRIDE


def natural(n, least=0):
    if type(n) is not int or n < least:
        raise ValueError('exact integer in declared domain required')
    return n


@dataclass(frozen=True)
class Entry:
    reset: int                 # zero means the inherited G17
    word: tuple
    P: int
    Q: int
    B: int
    source_residue: int
    exponent_residue: int
    exponent_period: int
    parameter_residue: int
    parameter_period: int


def compile_entry(reset):
    natural(reset)
    if reset == 0:
        word = forms.GENERATORS[17][3]
    else:
        if reset < 2 or reset % 2:
            raise ValueError('even terminal reset >=2 required')
        k = 0
        while (1 << (k+reset)) >= 3**(k+1):
            k += 1
        word = (1,)*k+(reset,)
    p,q,b = routes.carrier(word)
    if not (p>q and 0<b<3*q and (q-b)%p):
        raise ArithmeticError('positive all-height inverse-entry boundary lost')
    residue = b*pow(q,-1,p)%p
    log,period = forms.two_log_three_power((q+b)%p,len(word))
    phase = (log-sum(word))%period
    if phase%2!=1:
        raise ArithmeticError('Mersenne exponent must remain odd')
    m = period//2
    t = ((phase-BASE)//2)*pow(STRIDE//2,-1,m)%m
    return Entry(reset,word,p,q,b,residue,phase,period,t,m)


def audit(entry):
    if type(entry) is not Entry:
        raise ValueError('exact Entry required')
    for value in (entry.reset,entry.P,entry.Q,entry.B,entry.source_residue,
                  entry.exponent_residue,entry.exponent_period,
                  entry.parameter_residue,entry.parameter_period):
        natural(value)
    routes.letters(entry.word)
    if compile_entry(entry.reset)!=entry:
        raise ValueError('forged source/word/phase entry')
    return entry


def entries(max_reset=12):
    natural(max_reset,2)
    if max_reset%2:
        raise ValueError('even bank endpoint required')
    return (compile_entry(0),)+tuple(compile_entry(a) for a in range(2,max_reset+1,2))


def parameter_contains(entry,t):
    audit(entry);natural(t)
    return (t-entry.parameter_residue)%entry.parameter_period==0


def ordinary_source(entry,parameter):
    audit(entry);natural(parameter)
    r=entry.source_residue
    if r%2==0:r+=entry.P
    return r+2*entry.P*parameter


def receipt(entry,source):
    audit(entry);routes.odd(source)
    if source%entry.P!=entry.source_residue:
        raise ValueError('supplied source outside the exact inverse guard')
    child=(entry.Q*source-entry.B)//entry.P
    return routes.audit(routes.Receipt(source,child,(),entry.word,source))


def discharge(entry,source,supplied_child_root_word):
    return routes.discharge(receipt(entry,source),supplied_child_root_word)


def symbolic_child_mod(entry,t,modulus):
    """Same-source child residue, retaining the full ternary denominator."""
    audit(entry);natural(t);natural(modulus,1)
    if not parameter_contains(entry,t):
        raise ValueError('supplied parameter is not in this entry phase')
    precision=entry.P*modulus
    source=(pow(2,BASE+STRIDE*t,precision)-1)%precision
    numerator=entry.Q*source-entry.B
    if numerator%entry.P:
        raise ArithmeticError('lost source-owned exact division')
    return (numerator//entry.P)%modulus


def fuse(entry,binary_residue,binary_bits):
    """Unique intersection; this returns a guard, not a ROOT certificate."""
    audit(entry);natural(binary_residue);natural(binary_bits)
    m=1<<binary_bits
    if binary_residue>=m:
        raise ValueError('canonical binary residue required')
    n=entry.parameter_period
    lift=(entry.parameter_residue-binary_residue)*pow(m,-1,n)%n
    return binary_residue+m*lift,m*n


def coverage_core(bank):
    """Counting-only antichain. The original labelled bank is never discarded."""
    if type(bank) is not tuple:
        raise ValueError('tuple of labelled entries required')
    for entry in bank:audit(entry)
    core=[]
    for e in sorted(bank,key=lambda x:(x.parameter_period,x.parameter_residue,x.reset)):
        if not any((e.parameter_residue-c.parameter_residue)%c.parameter_period==0
                   for c in core):core.append(e)
    return tuple(core)


def coverage_mass(bank):
    return sum((F(1,e.parameter_period) for e in coverage_core(bank)),F(0))


def endpoint_mod3(t,depth):
    natural(t);natural(depth)
    modulus=3**depth
    p,q,b=routes.carrier(endpoint.PREFIX)
    numerator=2*pow(3,BASE+STRIDE*t+14,q*modulus)+b-p
    if numerator%q:raise ArithmeticError('lost endpoint denominator')
    return (numerator//q)%modulus


def locked_endpoint_residue(depth):
    natural(depth)
    if depth>BASE+14:
        raise ValueError('precision exceeds uniform locked ternary range')
    p,q,b=routes.carrier(endpoint.PREFIX)
    modulus=3**depth
    return (b-p)*pow(q,-1,modulus)%modulus


def inverse_depth_barrier(exponent):
    """Every inverse word at Y of length <=this is too short to pay M_E."""
    natural(exponent,16)
    return exponent-16


def main():
    checks=0
    def check(ok):
        nonlocal checks
        checks+=1
        if not ok:raise ArithmeticError('exact control failed')
    def reject(f):
        try:f()
        except (TypeError,ValueError):check(True)
        else:check(False)

    # Independent complete base-two cycles, rather than only the lifting reader.
    for depth in range(1,10):
        modulus=3**depth;order=2*3**(depth-1)
        seen={};z=1
        for k in range(order):
            check(z not in seen);seen[z]=k;z=2*z%modulus
        check(z==1 and len(seen)==2*3**(depth-1))
        for a in (2,4,6):
            e=compile_entry(a)
            if len(e.word)==depth:
                check((seen[(e.Q+e.B)%e.P]-sum(e.word))%order==e.exponent_residue)
    for a in range(2,21,2):
        e=compile_entry(a);k=len(e.word)-1
        check(e.B==e.P-2**(k+1) and e.Q<e.P and 2*e.P<3*e.Q)
        check(2**(k-1+a)>=3**k)
        check(pow(2,e.exponent_period,e.P)==1)
        check(pow(2,e.exponent_period//2,e.P)!=1)
        check(pow(2,e.exponent_period//3,e.P)!=1)
        for s in range(12):
            n=ordinary_source(e,s);r=receipt(e,n)
            check(0<r.child<n and r.endpoint==n)
            check(forms.replay(r.child,e.word)==n)
        for j in (0,1,2,19):
            t=e.parameter_residue+j*e.parameter_period
            E=BASE+STRIDE*t
            check((E-e.exponent_residue)%e.exponent_period==0)
            check((pow(2,E,e.P)-1)%e.P==e.source_residue)
            for m in (3,19,729,256):
                child=symbolic_child_mod(e,t,m)
                check((e.P*child+e.B-e.Q*(pow(2,E,m)-1))%m==0)
    e4,e6=compile_entry(4),compile_entry(6)
    check((e4.parameter_residue,e4.parameter_period)==(84,243))
    check((e6.parameter_residue,e6.parameter_period)==(1007,6561))
    check((receipt(e4,91).child,e4.word)==(63,(1,1,1,1,1,4)))
    old=(compile_entry(0),compile_entry(2));small=entries(6);wide=entries(12)
    check(coverage_mass(old)==F(244,729))
    check(coverage_mass(small)==F(2224,6561))
    check(coverage_mass(wide)==F(131325004,387420489))
    check(len(wide)==7 and len(coverage_core(wide))==6)
    check(not any(parameter_contains(e,0) for e in wide))
    # The small union is independently exhausted; no enormous period census.
    covered=0
    for t in range(6561):
        direct=any((pow(2,BASE+STRIDE*t,e.P)-1)%e.P==e.source_residue for e in small)
        phased=any((t-e.parameter_residue)%e.parameter_period==0 for e in small)
        check(direct==phased);covered+=direct
    check(covered==2224)
    # Independent CRT guards at all low binary residues; unique finite intersections.
    for h in range(7):
        for b in range(1<<h):
            r,m=fuse(e4,b,h)
            check(r%(1<<h)==b and r%243==84 and m==(1<<h)*243)
            check(sum(t%243==84 for t in range(b,m,1<<h))==1)
    check(fuse(e4,1,1)==(327,486))
    check(coverage_mass(small)-coverage_mass(old)==F(28,6561))
    # A supplied proof, not an orbit search, discharges an ordinary source.
    e2=compile_entry(2)
    check(receipt(e2,13).child==11)
    check(discharge(e2,13,(1,2,3,4))==(3,4))
    reject(lambda:discharge(e2,13,(2,3,4)))

    # Independent endpoint residues: all bounded ternary observations are locked.
    check(locked_endpoint_residue(2)==5)
    for depth in range(1,21):
        for t in (0,1,2,84,327,1007,10**20):
            check(endpoint_mod3(t,depth)==locked_endpoint_residue(depth))
    check((2*locked_endpoint_residue(2)-1)//3%3==0)
    check(3**30>2**47 and 3**29<2**46)
    p,q,b=routes.carrier(endpoint.PREFIX)
    for E in range(16,81):
        Y=F(2*3**(E+14)+b-p,q)
        r=inverse_depth_barrier(E)
        check(F(2,3)**r*(Y+1)-1>2**E-1)
    # The constant16 is the sharp integer cutoff for this g1-only comparison.
    Y=F(2*3**(80+14)+b-p,q)
    check(F(2,3)**(80-15)*(Y+1)-1<2**80-1)
    for bad in (True,4.0,-2,3):reject(lambda bad=bad:compile_entry(bad))
    reject(lambda:audit(replace(e4,parameter_residue=True)))
    reject(lambda:audit(replace(e4,B=e4.B+2)))
    reject(lambda:receipt(e4,True))
    reject(lambda:receipt(e4,93))
    reject(lambda:fuse(e4,2,1))
    reject(lambda:symbolic_child_mod(e4,0,19))
    reject(lambda:locked_endpoint_residue(BASE+15))
    reject(lambda:inverse_depth_barrier(15))
    print('PROVED even-reset first-crossing inverse entries; supplied source retained; ROOT terminal remains OPEN.')
    for e in wide:
        print('Entry',('G17' if e.reset==0 else f'a={e.reset}'),
              'word length/cost',len(e.word),sum(e.word),
              'parameter phase',e.parameter_residue,'mod',e.parameter_period)
    print('Old G5/G17 mass244/729; a4/a6 union2224/6561; strict gain28/6561.')
    print('Through a12: seven labelled alternatives, six counting cells; mass131325004/387420489; least hole t=0.')
    print('Any finite dyadic guard fuses by exact CRT; union density beta+tau-beta*tau, only parameter measure.')
    print('Endpoint Y always5mod9: G1 gives a multiple3; all inverse words of length<=E-16 cannot pay M_E.')
    print('No giant source expansion, no ROOT search, no claim that guard mass proves universal completion.')
    print('Exact checks',checks)


if __name__=='__main__':main()
