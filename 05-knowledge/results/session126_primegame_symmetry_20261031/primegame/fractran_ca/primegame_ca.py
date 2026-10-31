#!/usr/bin/env python3
"""Explicit finite-alphabet radius-one cellular automaton for fixed PRIMEGAME.

local_rule(left, center, right) is the complete homogeneous CA rule, not a
macro FRACTRAN interpreter. The tape head has only 93 enumerated states.
Coordinates/unbounded counters never occur in a cell symbol.
"""
from pathlib import Path
from math import prod

PRIMES=(2,3,5,7,11,13,17,19,23,29)
FRACTIONS=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
COUNTER_MASK=(1<<10)-1
ORIGIN=1<<10
SPECTATOR=1<<11
TAPE_SIZE=1<<12

def _support(n):
    indices=[]
    for i,p in enumerate(PRIMES):
        if n%p==0:
            n//=p;indices.append(i)
            if n%p==0: raise ValueError('fixed program must be squarefree')
    if n!=1: raise ValueError('unsupported fixed-program prime')
    return tuple(indices)

# Both lists are fixed finite constants determined at module construction.
DENOMINATOR_SUPPORTS=tuple(_support(d) for _,d in FRACTIONS)
NUMERATOR_SUPPORTS=tuple(_support(n) for n,_ in FRACTIONS)
DENOMINATOR_MASKS=tuple(sum(1<<p for p in s) for s in DENOMINATOR_SUPPORTS)
OPS=tuple(tuple((p,-1) for p in take)+tuple((p,1) for p in give)
          for take,give in zip(DENOMINATOR_SUPPORTS,NUMERATOR_SUPPORTS))
if any(set(t)&set(g) for t,g in zip(DENOMINATOR_SUPPORTS,NUMERATOR_SUPPORTS)):
    raise ValueError('fixed fractions must be reduced')

HEAD_STATES=[('READY',-1,-1)]
for rule,ops in enumerate(OPS):
    for op,(_,sign) in enumerate(ops):
        HEAD_STATES.append(('SEEK',rule,op))
        if sign<0: HEAD_STATES.append(('DECLEFT',rule,op))
        HEAD_STATES.append(('RETURN',rule,op))
HEAD_STATES=tuple(HEAD_STATES)
HEAD_CODE={state:i+1 for i,state in enumerate(HEAD_STATES)}
READY=HEAD_CODE[('READY',-1,-1)]
HEAD_COUNT=len(HEAD_STATES)
ERROR=TAPE_SIZE*(HEAD_COUNT+1)
ALPHABET_SIZE=ERROR+1
BACKGROUND=0

# All alphabet symbols are integers 0..ERROR. For a nonerror symbol,
# divmod(symbol,TAPE_SIZE) returns (head code 0..93, 12-bit tape symbol).
def symbol(data,head=0):
    if not isinstance(data,int) or not 0<=data<TAPE_SIZE: raise ValueError('invalid tape symbol')
    if not isinstance(head,int) or not 0<=head<=HEAD_COUNT: raise ValueError('invalid head state')
    return data+TAPE_SIZE*head

def unpack(cell):
    if not isinstance(cell,int) or not 0<=cell<ERROR: raise ValueError('not a nonerror symbol')
    head,data=divmod(cell,TAPE_SIZE)
    return data,head

def _next(rule,op):
    return READY if op+1==len(OPS[rule]) else HEAD_CODE[('SEEK',rule,op+1)]

def head_delta(data,head):
    """Finite head transition: (new tape, new head, direction), or None=error.
    No tape coordinate, counter value or unbounded queue is an input/output.
    """
    if not 0<=data<TAPE_SIZE or not 1<=head<=HEAD_COUNT:return None
    phase,rule,op=HEAD_STATES[head-1]
    at_origin=bool(data&ORIGIN)
    if phase=='READY':
        if not at_origin:return None
        bits=data&COUNTER_MASK
        rule=next(i for i,mask in enumerate(DENOMINATOR_MASKS) if bits&mask==mask)
        return data,HEAD_CODE[('SEEK',rule,0)],0
    coordinate,sign=OPS[rule][op]
    bit=1<<coordinate
    if phase=='SEEK':
        if data&bit:return data,head,1
        if sign>0:
            data|=bit
            if at_origin:return data,_next(rule,op),0
            return data,HEAD_CODE[('RETURN',rule,op)],-1
        if at_origin:return None # attempted decrement of zero counter
        return data,HEAD_CODE[('DECLEFT',rule,op)],-1
    if phase=='DECLEFT':
        if sign!=-1 or not data&bit:return None
        data&=~bit
        if at_origin:return data,_next(rule,op),0
        return data,HEAD_CODE[('RETURN',rule,op)],-1
    if phase=='RETURN':
        if at_origin:return data,_next(rule,op),0
        return data,head,-1
    raise RuntimeError('nonenumerated finite head phase')

def local_rule(left,center,right):
    """A total radius-one map A^3 -> A for the explicit finite alphabet A.

    ERROR spreads; neighborhoods with >=2 heads, or an invalid local head
    instruction, become ERROR. Other malformed global tapes still receive
    a defined update; no global validity test is hidden in this function.
    """
    cells=(left,center,right)
    if any(not isinstance(x,int) or not 0<=x<=ERROR for x in cells):
        raise ValueError('arguments must belong to finite alphabet')
    if ERROR in cells:return ERROR
    decoded=[divmod(x,TAPE_SIZE) for x in cells] # (head,data)
    locations=[i for i,(head,_) in enumerate(decoded) if head]
    if len(locations)>1:return ERROR
    if not locations:return center
    hloc=locations[0];head,data=decoded[hloc]
    transition=head_delta(data,head)
    if transition is None:return ERROR
    new_data,new_head,direction=transition
    if hloc==1:
        return symbol(new_data,new_head if direction==0 else 0)
    # A neighboring head moves in only when its direction points to center.
    if hloc+direction==1:return symbol(decoded[1][1],new_head)
    return center

def ca_step_full(configuration):
    """Literal synchronous evaluation on finite-support configurations."""
    if not configuration:return {}
    candidates={x+d for x in configuration for d in (-1,0,1)}
    result={}
    for x in candidates:
        y=local_rule(configuration.get(x-1,0),configuration.get(x,0),configuration.get(x+1,0))
        if y:result[x]=y
    return result

def ca_step(configuration):
    """Equivalent sparse evaluation: head/error-free triples copy center.

    This simulator optimization does not alter the local rule or inspect
    counter values. Static tape cells remain unchanged unless within radius1
    of a head/error, so only that finite candidate set needs reevaluation.
    """
    active=(x for x,v in configuration.items() if v==ERROR or v//TAPE_SIZE)
    candidates={x+d for x in active for d in (-1,0,1)}
    result={x:v for x,v in configuration.items() if v}
    for x in candidates:
        y=local_rule(configuration.get(x-1,0),configuration.get(x,0),configuration.get(x+1,0))
        if y:result[x]=y
        else:result.pop(x,None)
    return result

def encode(integer,origin=0):
    """All positive integers: ten unary tracks + binary spectator factor."""
    if not isinstance(integer,int) or integer<1:raise ValueError('positive integer required')
    residual=integer;configuration={}
    for i,p in enumerate(PRIMES):
        height=0
        while residual%p==0:residual//=p;height+=1
        for x in range(height):configuration[origin+x]=configuration.get(origin+x,0)|(1<<i)
    configuration[origin]=symbol(configuration.get(origin,0)|ORIGIN,READY)
    # Encode the whole coprime residual (including1) in passive bits left of origin.
    for j in range(residual.bit_length()):
        if residual>>j&1:configuration[origin-1-j]=SPECTATOR
    return configuration

def observe(configuration,require_ready=True):
    """External observer, not part of the CA rule. Check encoding and decode."""
    if any(v==ERROR for v in configuration.values()):raise ValueError('error configuration')
    decoded={x:unpack(v) for x,v in configuration.items()}
    origins=[x for x,(data,_) in decoded.items() if data&ORIGIN]
    heads=[(x,head) for x,(_,head) in decoded.items() if head]
    if len(origins)!=1 or len(heads)!=1:raise ValueError('one origin and head required')
    origin=origins[0]
    if heads[0][0]<origin:raise ValueError('head crossed origin')
    if require_ready and heads!=[(origin,READY)]:raise ValueError('not a READY boundary')
    exponents=[]
    for i in range(10):
        sites=sorted(x-origin for x,(data,_) in decoded.items() if data&(1<<i))
        if sites!=list(range(len(sites))):raise ValueError('nonunary counter shape')
        exponents.append(len(sites))
    spectator=0
    for x,(data,_) in decoded.items():
        if data&SPECTATOR:
            if x>=origin:raise ValueError('spectator bits must lie left of origin')
            spectator|=1<<(origin-1-x)
    if spectator<1 or any(spectator%p==0 for p in PRIMES):raise ValueError('noncanonical spectator')
    return {'integer':spectator*prod(p**e for p,e in zip(PRIMES,exponents)),
            'exponents':tuple(exponents),'spectator':spectator,'origin':origin,
            'head_position':heads[0][0],'head_state':heads[0][1]}

def ready(configuration):
    return any(v!=ERROR and v//TAPE_SIZE==READY for v in configuration.values())

def advance_macro(configuration,cap=None,verify_micro=False,use_full=False):
    """Drive actual local CA to its next READY boundary; finite safety cap."""
    initial=observe(configuration)
    M=max(initial['exponents']);bound=10*M+6
    if cap is None:cap=bound
    step=ca_step_full if use_full else ca_step
    for elapsed in range(1,cap+1):
        configuration=step(configuration)
        if verify_micro:observe(configuration,require_ready=False)
        if ready(configuration):
            final=observe(configuration)
            return configuration,elapsed,final
    raise RuntimeError('CA safety bound reached; not a mathematical halt')

if __name__=='__main__':
    print({'tape_symbols':TAPE_SIZE,'head_states':HEAD_COUNT,'alphabet_size':ALPHABET_SIZE,
           'operation_counts':[len(x) for x in OPS],'total_operations':sum(map(len,OPS)),
           'decrement_microstates':sum(sign<0 for ops in OPS for _,sign in ops)})
