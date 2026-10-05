"""Compose a reset-two sibling transport before paying the source-size bound.

PROVED formulas are in the companion note. The finite replay controls and
supplied-certificate demonstrations do not assume universal Collatz convergence.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import json

import inverse_ray_ternary_addresses_20261004 as codec
import debt_word_join_search_20261004 as debt
import collatz_binary_ternary_guard_fusion_20261004 as fusion


def need(ok, message):
    if not ok:
        raise ValueError(message)


def odd(n):
    need(type(n) is int and n > 0 and n % 2, 'positive exact odd source')


def step(n):
    odd(n)
    raw = 3*n+1
    a = (raw & -raw).bit_length()-1
    return raw >> a, a


def replay(n, word, first_hit=True):
    odd(n)
    states = [n]
    for a in word:
        need(type(a) is int and a >= 1, 'positive exact valuation')
        need(not first_hit or n != 1, 'no first-hit root padding')
        raw, actual = 3*n+1, 0
        while raw % 2 == 0:
            raw //= 2
            actual += 1
        need(actual == a, 'independent actual valuation')
        n = raw
        states.append(n)
    return tuple(states)


@lru_cache(None)
def old_row(k):
    need(type(k) is int and k >= 0, 'nonnegative old sibling index')
    ell = 1
    while 3**ell <= 2**(ell+2*k):
        ell += 1
    M = 3**ell
    r = -((4**k+2)//3)*pow(4**k,-1,M) % M
    return ell, M, r


def old_sibling_index(n):
    odd(n)
    # Inherited positivity bound: ell+1 < log2(n+1), and ell>2k.
    for k in range((n+1).bit_length()//2+1):
        _, M, r = old_row(k)
        if n % M == r:
            return k
    return None


def load_bank():
    root = Path(__file__).resolve().parents[2]
    saved = json.loads((root/'05-knowledge/results/collatz_binary_ternary_guard_fusion_20261004.json').read_text())
    bank = {(row['s'],row['e']):(tuple(row['source_word']),tuple(row['child_word']))
            for row in saved['debt_rows']}
    need(len(bank) == 16, 'frozen binary bank')
    return bank


def old_partition(n, bank):
    odd(n)
    if n == 1:
        return 'ROOT'
    if step(n)[0] < n:
        return 'direct'
    run = ((n+1) & -(n+1)).bit_length()-2
    if run > 0:
        endpoint = 3**run*(n+1)//2**run-1
        if step(endpoint)[1] >= 3:
            return 'reset-at-least-three'
    if fusion.select_debt(n,bank) is not None:
        return 'binary16'
    if old_sibling_index(n) is not None:
        return 'ternary'
    return None


@dataclass(frozen=True)
class RoutingFamily:
    k: int
    ell: int
    modulus: int
    residue: int
    slope_numerator: int
    intercept_numerator: int
    source: int
    period: int
    child: int
    child_period: int

    @property
    def slope(self):
        return Fraction(self.slope_numerator,self.modulus)

    @property
    def source_word(self):
        return (1,2,1)

    @property
    def child_word(self):
        return (1,)*(self.ell-1)+(2,2*self.k-1)

    def instance(self,t):
        need(type(t) is int and t >= 0, 'nonnegative exact family parameter')
        return self.source+self.period*t, self.child+self.child_period*t

    def parameter(self,n):
        odd(n)
        if n < self.source or (n-self.source) % self.period:
            return None
        return (n-self.source)//self.period


def family(k,ell=None):
    need(type(k) is int and k >= 1 and k % 9 == 1, 'feasible sibling address k=1 mod9')
    least = 3
    while 3**(least-2) <= 2**(least+2*k-4):
        least += 1
    if ell is None:
        ell = least
    need(type(ell) is int and ell >= least, 'enough inverse depth for the composed source-size bound')
    M, P = 3**(ell-2), 2**(ell+2*k-4)
    numerator = 23*4**k+16
    need(numerator % 27 == 0, 'ternary feasibility before cancelling27')
    g = numerator//27
    shifted_g = g*2**max(ell-4,0)//2**max(4-ell,0)
    need(shifted_g*2**max(4-ell,0) == g*2**max(ell-4,0), 'integral intercept coefficient')
    C = shifted_g-M
    r = -g*pow(4**k,-1,M) % M
    n = r+M*((123-r)*pow(M,-1,128) % 128)
    hnum = P*n+C
    need(hnum % M == 0, 'actual child integrality')
    h = hnum//M
    result = RoutingFamily(k,ell,M,r,P,C,n,128*M,h,128*P)
    need(0 < h < n and h % 2 and P < M and C < 0, 'strict all-height original-source payment')
    return result


def checked_instance(row,t):
    n,h = row.instance(t)
    left = replay(n,row.source_word)
    right = replay(h,row.child_word)
    need(left[-1] == right[-1] and 0 < h < n, 'actual common future and smaller final child')
    need(row.modulus*h == row.slope_numerator*n+row.intercept_numerator, 'ordered affine child map')
    first_four = replay(n,(1,2,1,2))
    need(all(value > n for value in first_four[1:]), 'source grows throughout its first four odd steps')
    m = 2*h+1
    middle = replay(m,(1,)*row.ell+(2*row.k+1,))
    need(middle[-1] == left[-1], 'unranked sibling intermediate has the same future')
    return n,h,left[-1],m


def partitioned_selector(n,bank,row):
    previous = old_partition(n,bank)
    if previous is not None:
        return previous,None
    t = row.parameter(n)
    if t is None:
        return 'PENDING',None
    source,child = row.instance(t)
    need(source == n and child < n, 'new disjoint class retains the original source')
    return 'composed-reset-two',child


@lru_cache(None)
def general_spec(run,k):
    need(type(run) is int and run >= 1 and type(k) is int and k >= 1
         and k % 3**(run+1) == 1,'general run/address feasibility')
    ell = run+2
    while 2**(ell+2*k-run-3) >= 3**(ell-run-1):
        ell += 1
    M,P = 3**(ell-run-1),2**(ell+2*k-run-3)
    Cnumer = 2**(run+1)*(4**k-4)
    need(Cnumer % 3**(run+2) == 0,'general integral carry')
    C = Cnumer//3**(run+2)
    shift = ell-run-3
    shifted = C*2**max(shift,0)//2**max(-shift,0)
    residue = (C*pow(4**k,-1,M)-1) % M
    binary,D = debt.cylinder((1,)*run+(2,))
    n = residue+M*((binary-residue)*pow(M,-1,D)%D)
    h = (P*(n+1)-shifted)//M-1
    need(P*(n+1)-shifted == M*(h+1) and 0<h<n and h%2,'general actual child and original-source payment')
    return dict(run=run,k=k,ell=ell,M=M,P=P,C=C,shifted=shifted,residue=residue,
                source=n,period=D*M,child=h,child_period=D*P)


def general_candidates(n):
    """Complete applicability search within this schema, by the real height cap."""
    odd(n)
    if n == 1 or n % 4 != 3:
        return ()
    run = ((n+1)&-(n+1)).bit_length()-2
    x,a = step(3**run*(n+1)//2**run-1)
    if a != 2:
        return ()
    maxell = n.bit_length()-1
    if maxell < run+2:
        return ()
    out,k = [],1
    while 2*k < maxell-run+1:  # Necessary from 3^a<4^a; avoids huge impossible indices.
        if 2**(maxell+2*k-run-3) >= 3**(maxell-run-1):
            break
        row = general_spec(run,k)
        if n % row['M'] == row['residue']:
            h = (row['P']*(n+1)-row['shifted'])//row['M']-1
            endpoint,a = step(x)
            need(0<h<n and h%2,'guarded final child at this supplied source')
            out.append(dict(child=h,k=k,ell=row['ell'],endpoint=endpoint,
                            source_word=(1,)*run+(2,a),
                            child_word=(1,)*(row['ell']-1)+(2,a+2*k-2)))
        k += 3**(run+1)
    return tuple(out)


def predecessor_strip(h):
    odd(h)
    e, unit = 0,h+1
    while unit % 3 == 0:
        unit //= 3
        e += 1
    result = (unit << e)-1
    need(result > 0 and result % 2 and result <= h and (e == 0 or result < h),
         'each k0 predecessor preserves positivity and strictly lowers its source')
    need(result % 3 in (0,1), 'the final k0 guard is absent')
    need(replay(result,(1,)*e)[-1] == h, 'retained inverse-one prefix')
    return e,result


def attach_supplied_child(row,t,child_certificate):
    n,h = row.instance(t)
    codec.audit_certificate(child_certificate)
    need(codec.expand(child_certificate,bit_cap=codec.bit_bounds(child_certificate)[1]) == h,
         'supplied certificate belongs to this actual child')
    suffix = child_certificate
    for a in row.child_word:
        need(suffix != codec.ROOT and codec.exponent(suffix) == a, 'checked child prefix is before first root')
        suffix = suffix.parent
    states = replay(n,row.source_word)
    need(codec.expand(suffix,bit_cap=codec.bit_bounds(suffix)[1]) == states[-1], 'supplied common-future suffix')
    for i in range(len(row.source_word)-1,-1,-1):
        source,target = states[i:i+2]
        a = row.source_word[i]
        least = codec.kappa(target % 9,source % 3)
        need(a >= least and (a-least) % 6 == 0, 'exact inverse route address')
        suffix = codec.extend(suffix,source % 3,(a-least)//6)
    need(codec.expand(suffix,bit_cap=codec.bit_bounds(suffix)[1]) == n, 'certificate keeps supplied source')
    return suffix


def main():
    print('COMPLEMENT ROUTING: all-height source-paid composition; global family coverage OPEN')
    bank = load_bank()
    broad,phase = family(10),family(10,35)
    need((broad.ell,broad.modulus,broad.source,broad.period,broad.child,broad.child_period) ==
         (33,617673396283947,27727075633746555,79062194724345216,
          25270565367447551,72057594037927936), 'exact broad family constants')
    print('BROAD',json.dumps(broad.__dict__,sort_keys=True))
    print('PHASE',json.dumps(phase.__dict__,sort_keys=True))
    print('Exact words:source(1,2,1), child1^(ell-1),2,(2k-1); source firstfour(1,2,1,2) all grow.')
    ancestor_table = []
    for j in range(10):
        ell,M,r = old_row(j)
        need(ell <= 31 and broad.residue % M != r, 'no inherited ternary ancestor j0..9')
        ancestor_table.append((j,ell,broad.residue % M,r))
    need(old_row(10)[0] == 35, 'first remaining old guard depth')
    need(Fraction(3**31,3**35)*Fraction(27,26) == Fraction(1,78), 'all-index relative tail bound')
    new_density = Fraction(77,78)*Fraction(1,64*3**31)
    print('ANCESTOR_TABLE',json.dumps(ancestor_table))
    print('Old ternary relative occupancy <=1/78; new fraction>=77/78; odd-relative new density>=',str(new_density))
    controls = 0
    family_rows = []
    for k in range(1,182,9):
        row = family(k)
        for t in (0,1,2,7,31,10**6):
            checked_instance(row,t)
            controls += 1
        family_rows.append((k,row.ell,row.ell+2*k-4,row.ell-2))
        d = (2-row.ell) % 3
        refined = family(k,row.ell+d)
        need(refined.period == 3**d*row.period and
             (refined.source-row.source) % row.period == 0,'general phase refinement of source guard')
        offset = (refined.source-row.source)//row.period
        need(0 <= offset < 3**d,'canonical added ternary digits')
        for t in (0,1,19):
            n,h = row.instance(offset+3**d*t)
            nn,hh = refined.instance(t)
            need(n == nn and 3**d*(hh+1) == 2**d*(h+1),'general precision versus child predecessor law')
            need((h+1-pow(7,row.ell-2,19)*(n+1)) % 19 == 0,'general shifted phase formula')
            need((hh-n) % 19 == 0,'general phase-preserving refinement')
    print('GENERAL_FAMILY_ROWS [k,ell,power_of2,power_of3_in_slope]',json.dumps(family_rows))
    print('General literal all-word controls:',controls,'(k1..181 step9, six parameters)')
    run_controls = []
    for run,k in ((1,10),(2,28),(3,82)):
        spec = general_spec(run,k)
        for t in (0,1,17):
            n = spec['source']+spec['period']*t
            h = spec['child']+spec['child_period']*t
            x = replay(n,(1,)*run+(2,))[-1]
            J,a = step(x)
            need(replay(h,(1,)*(spec['ell']-1)+(2,a+2*k-2))[-1] == J,'general-r exact join')
            need(0<h<n and (h+1-pow(7,spec['ell']-run-1,19)*(n+1))%19 == 0,'general-r rank and phase')
            need(any(row['k']==k and row['child']==h for row in general_candidates(n)),
                 'finite supplied-source search finds the actual general-r member')
        run_controls.append((run,k,spec['ell']))
    print('General initial-run controls [r,k,ell]:',run_controls,';three parameters each')
    for j in range(1,25):
        n = 2**(2*j+1)-1
        need(general_candidates(n) == (),'all-length Mersenne obstruction has no eligible depth')
    print('All-height hostile:2^(2j+1)-1,j>=1 is outside this schema;24 literal selector controls.')
    hits = 0
    for t in range(512):
        n,h,_,_ = checked_instance(broad,t)
        need(fusion.select_debt(n,bank) is None, 'binary16 exclusion on the whole first-four prefix')
        label,child = partitioned_selector(n,bank,broad)
        if label == 'composed-reset-two':
            need(child == h,'new partition keeps exact child')
            hits += 1
        if t < 32:
            need((old_sibling_index(n) is None) == (fusion.sibling_child(n) is None),
                 'independent inherited selector agreement')
    print('Broad controls:first512 parameters;',hits,'new-class dispatches; others retained in earlier class')
    redundant = family(1)
    need(old_sibling_index(redundant.source) == 0, 'k1 paid construction is redundant with old ternary guard')
    # The original one-stage clock fails at k19; the composite has a smaller final child.
    old_ell = old_row(19)[0]
    n,h,_,m = checked_instance(family(19,old_ell),0)
    need(m > n > h,'uncomposed payment fails but composition pays against original source')
    print('HOSTILES:k1 already covered; k19 at oldell',old_ell,'has intermediate>source>finalchild')
    # Two ternary precision digits are exchanged for two inherited k0 steps.
    need(phase.source == broad.source+6*broad.period and phase.period == 9*broad.period,
         'nested source address t=6 mod9')
    D = 3**33-2**31*893232
    need(D == 19*191624181720273 and 2**51-3**33 == -19*174066355414225,
         'exact mod19 coefficient factorization')
    need(phase.period % 19 == phase.child_period % 19 == 16 and
         phase.source % 19 == phase.child % 19,'all-height phase-preserving map')
    need((3*pow(2,-1,19))**3 % 19 == 1,'three-phase inverse-one clock')
    need(pow(4,9,19) == 1 and pow(7,3,19) == 1 and 7 != 1,'exact general phase clocks')
    phase_controls = 0
    depths = set()
    for t in range(512):
        n,h = broad.instance(t)
        need((h+1-7*(n+1)) % 19 == 0,'broad shifted mod19 multiplier7')
        e,g = predecessor_strip(h)
        depths.add(e)
        need((g+1-pow(7,e+1,19)*(n+1)) % 19 == 0,'exact phase after maximal k0 strip')
        if e % 3 == 2:
            need((g-n)%19 == 0,'phase returns exactly at the stated depth class')
        phase_controls += 1
    for u in (0,1,2,19,10**6):
        n,h = broad.instance(6+9*u)
        n2,h2 = phase.instance(u)
        need(n == n2 and 4*(h+1) == 9*(h2+1),'two-digit refinement and double predecessor equality')
        checked_instance(phase,u)
    need(len({phase.instance(t)[0] % 19 for t in range(19)}) == 19,'every mod19 source phase occurs')
    print('Nested refinement:t=6 mod9; two k0 steps give phase family; all19 phases occur.')
    print('Maximal k0 partition controls:',phase_controls,'; observed exact depths',sorted(depths))
    # Retained proofs are supplied premises only in these explicit demonstration cases.
    for t in (0,1,2,17):
        n,h = broad.instance(t)
        supplied = codec.encode_source(h,step_cap=5000)
        result = attach_supplied_child(broad,t,supplied)
        need(result == codec.encode_source(n,step_cap=5000),'independent first-hit source audit')
    print('Supplied-child AST demonstrations:4, bounded premise searches cap5000; no free-parameter home inference.')
    for n in (7,27,703):
        need(partitioned_selector(n,bank,broad) == ('PENDING',None),'canonical residual hostiles remain pending')
    rejected = 0
    for thunk in (lambda: family(2), lambda: family(10,32),lambda: family(True),
                  lambda: broad.instance(-1),lambda: broad.parameter(1.0),
                  lambda: attach_supplied_child(broad,0,codec.ROOT)):
        try:
            thunk()
        except ValueError:
            rejected += 1
        else:
            raise ValueError('invalid guard was accepted')
    print('Guard/type/source rejections:',rejected,';7/27/703 remain outside this new class.')
    print('PASS: compose first, then test original-source size; disjoint extension has positive density, not universal coverage.')


if __name__ == '__main__':
    main()
