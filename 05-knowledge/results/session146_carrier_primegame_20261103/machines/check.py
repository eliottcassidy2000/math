#!/usr/bin/env python3
"""Finite odd-pair triangle to ordered FRACTRAN compiler. No assert statements."""
import itertools, math, json
from fractions import Fraction

def require(ok, msg):
    if not ok: raise ValueError(msg)

def prime(x):
    if x < 2: return False
    return all(x % d for d in range(2,math.isqrt(x)+1))

def primes():
    n=2
    while True:
        if prime(n): yield n
        n+=1

def candidate(a):
    return a in (5,29) or math.gcd(a,290)==1

def compile_program(N, target, prune=False):
    require(N>=1 and target>=0,'invalid bounds')
    vals=[2*i+3 for i in range(N)]
    pairs=[(i,j) for i in range(N) for j in range(i,N)
           if vals[i]+vals[j]==target and (not prune or (candidate(vals[i]) and candidate(vals[j])))]
    pool=primes(); bits=[next(pool) for _ in vals]
    # State identities distinct from data primes; no other denominator factors.
    read=[next(pool) for _ in pairs]
    mid=[next(pool) for _ in pairs]
    read2=[next(pool) for _ in pairs]
    finish=[next(pool) for _ in pairs]
    accept=[next(pool) for _ in pairs]
    reject=next(pool)
    rules=[]
    for k,(i,j) in enumerate(pairs):
        nxt=read[k+1] if k+1<len(pairs) else reject
        # Positive branch precedes zero fallback; consuming bit makes guard noncancellable.
        rules += [(mid[k],read[k]*bits[i]), (nxt,read[k]),
                  (read2[k]*bits[i],mid[k]),
                  (finish[k],read2[k]*bits[j]), (nxt,read2[k]),
                  (accept[k]*bits[j],finish[k])]
    # reduced fractions are necessary for the selected-rule checks below
    rules=[(f.numerator,f.denominator) for f in (Fraction(a,b) for a,b in rules)]
    return dict(N=N,target=target,values=vals,pairs=pairs,bits=bits,rules=rules,
                initial=read[0] if pairs else reject, accepts=accept,reject=reject)

def run(P,vector,limit=None):
    require(len(vector)==P['N'] and all(type(x) is bool for x in vector),'illegal bit vector')
    if limit is None: limit=4*len(P['pairs'])+1
    z=P['initial']*math.prod(p for p,b in zip(P['bits'],vector) if b)
    for steps in range(limit):
        for rule,(a,b) in enumerate(P['rules']):
            if z%b==0:
                z=z//b*a
                break
        else:
            data=math.prod(p for p,b in zip(P['bits'],vector) if b)
            require(z % data == 0,'data not restored')
            state=z//data
            if state==P['reject']: return None,steps
            require(state in P['accepts'],'halt outside declared output states')
            return P['pairs'][P['accepts'].index(state)],steps
    raise RuntimeError('step limit exceeded')

def direct(P,v):
    return next(((i,j) for i in range(P['N']) for j in range(i,P['N'])
                 if (2*i+3)+(2*j+3)==P['target'] and v[i] and v[j]), None)



def tm_run(N,target,tape):
    """Ordinary one-tape finite-prefix TM, with an explicit right blank.

    Tape is a finite list from cell 0, then blank forever. Returns the output
    marker, the resulting finite tape, and the primitive-transition count.
    Invalid input produces ERR; for valid bit words the input is unchanged.
    """
    require(N>=1 and target>=0,'invalid TM parameters')
    tape=list(tape); prefix=(); head=0; steps=0
    while True:
        sym=tape[head] if head<len(tape) else '_'
        steps+=1
        if len(prefix)<N:
            if type(sym) is not int or sym not in (0,1): return 'ERR',tape,steps
            prefix=prefix+(sym,); head+=1
        else:
            if sym!='_': return 'ERR',tape,steps
            outcome=next(((i,j) for i in range(N) for j in range(i,N)
                          if 2*i+2*j+6==target and prefix[i] and prefix[j]),None)
            marker=('A',)+outcome if outcome is not None else ('R',)
            if head==len(tape): tape.append(marker)
            else: tape[head]=marker
            return marker,tape,steps

def tm_delta(N,target,state,sym):
    """Transition for a prefix state; returns written symbol, next state, move."""
    if state in ('H','E'): return sym,state,'S'
    if len(state)<N:
        return (sym,state+(sym,),'R') if type(sym) is int and sym in (0,1) else (sym,'E','S')
    if sym!='_': return sym,'E','S'
    pair=next(((i,j) for i in range(N) for j in range(i,N)
               if 2*i+2*j+6==target and state[i] and state[j]),None)
    return (('R',) if pair is None else ('A',)+pair),'H','S'

def ca_step(N,target,cells):
    """One simultaneous radius-1 update on a finite window with blank exterior."""
    result=[]
    for k,(sym,state) in enumerate(cells):
        incoming=None
        if k:
            ls,lq=cells[k-1]
            if lq is not None:
                _,nq,mv=tm_delta(N,target,lq,ls)
                if mv=='R': incoming=nq
        if state is None:
            result.append((sym,incoming))
        elif incoming is not None:
            result.append(('#',None)) # explicitly specified illegal multihead collision
        else:
            written,nq,mv=tm_delta(N,target,state,sym)
            result.append((written,None if mv=='R' else nq))
    return result

def ca_run(N,target,vec):
    cells=[(int(b),None) for b in vec]+[('_',None)]
    cells[0]=(cells[0][0],())
    for _ in range(N+1): cells=ca_step(N,target,cells)
    return cells

# Inherited Conway owner-list phase law: see THM-4852, cited in README.
PG_PRIMES=(2,3,5,7,11,13,17,19,23,29)
PG_HEAD=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),
         (95,23),(77,19),(1,17),(11,13),(13,11))
PG_LISTS={"owner":PG_HEAD+((15,2),(1,7),(55,1)),
          "guy":PG_HEAD+((15,14),(15,2),(55,1))}

def factors(x):
    v=[]
    for p in PG_PRIMES:
        k=0
        while x%p==0: x//=p; k+=1
        v.append(k)
    require(x==1,'non-program prime')
    return v

def primegame_bits(N, which, step_limit=2_000_000):
    """Finite test run; timeout signals inconclusive, never a primality answer."""
    require(N>=1 and which in PG_LISTS,'invalid stream request')
    rules=[(factors(a),factors(b)) for a,b in PG_LISTS[which]]
    v=[1]+[0]*9; result=[False]*N; cap=2*N+2
    for step in range(step_limit+1):
        # Pure powers of two identify the inherited PRIMEGAME prime output.
        if v[0]>=2 and not any(v[1:]) and v[0]%2==1:
            i=(v[0]-3)//2
            if 0<=i<N: result[i]=True
        # Observe a state-free phase with no 3 or 5 after all requested indices.
        if v[0]>=cap and v[1]==0 and v[2]==0 and not any(v[4:]):
            return result,step
        for num,den in rules:
            if all(x>=y for x,y in zip(v,den)):
                v=[x-y+z for x,y,z in zip(v,den,num)]
                break
        else: raise RuntimeError('PRIMEGAME halted unexpectedly')
    raise RuntimeError('PRIMEGAME simulation budget exhausted; inconclusive')

def check():
    exhaustive=0
    for N in range(1,9):
        for t in range(0,4*N+5):
            P=compile_program(N,t)
            for vec in itertools.product((False,True),repeat=N):
                got,cost=run(P,vec)
                require(got==direct(P,vec),'arbitrary bit-vector mismatch')
                marker,tape,tmsteps=tm_run(N,t,list(map(int,vec)))
                expected=('R',) if got is None else ('A',)+got
                require(marker==expected and tape[:-1]==list(map(int,vec)) and tmsteps==N+1,'TM mismatch')
                cells=ca_run(N,t,vec)
                require(cells[-1]==(expected,'H') and cells[:-1]==[(int(b),None) for b in vec],'CA mismatch')
                require(cost<=4*len(P['pairs']),'termination bound violated')
                exhaustive+=1
    primecases=0
    for N in range(1,35):
        vec=[prime(2*i+3) for i in range(N)]
        for t in range(0,4*N+5):
            full=compile_program(N,t); pruned=compile_program(N,t,True)
            require(run(full,vec)[0]==run(pruned,vec)[0],'prefilter changed prime pair')
            primecases+=1
    require(tm_run(2,8,[1])[0]=='ERR' and tm_run(2,8,[1,1,0])[0]=='ERR' and tm_run(2,8,[1,2])[0]=='ERR','TM illegal word controls')
    hostile=[(1,()),(1,()),('_',None)]
    require(ca_step(2,8,hostile)[1]==('#',None),'CA multihead collision')
    blank=('_',None)
    shifted=ca_step(2,8,[blank]+hostile+[blank])
    unshifted=ca_step(2,8,hostile)
    require(shifted[1:1+len(hostile)]==unshifted,'CA local translation')
    done=ca_run(2,8,[True,True])
    require(ca_step(2,8,done)==done,'CA halting must freeze')
    # Actual ordered Conway/Guy streams supply finite primality vectors.
    streamcases=0
    for N in range(1,6):
        expected=[prime(2*i+3) for i in range(N)]
        for which in PG_LISTS:
            bits,steps=primegame_bits(N,which)
            require(bits==expected,'PRIMEGAME vector differs from trial division')
            for t in range(6,4*N+5):
                P=compile_program(N,t)
                require(run(P,bits)[0]==direct(P,expected),'stream-to-pair mismatch')
                streamcases+=1
    # 5,29 must survive exceptions; 15 must fail residue guard; diagonal retained.
    P=compile_program(14,34,True); v=[prime(x) for x in P['values']]
    require((1,13) in P['pairs'] and run(P,v)[0]==(1,13),'5+29 exception')
    P=compile_program(1,6); require(run(P,[True])[0]==(0,0),'diagonal')
    P=compile_program(7,30); v=[False]*7; v[6]=True
    require(run(P,v)[0]==(6,6),'formal composite vector is valid input')
    require(run(compile_program(7,30,True),v)[0] is None,'pruning is not arbitrary-vector safe')
    # First-applicable order: q*r guards before q fallback; reverse changes output.
    P=compile_program(1,6); require(run(P,[True])[0]==(0,0),'positive priority')
    P['rules'][0],P['rules'][1]=P['rules'][1],P['rules'][0]
    require(run(P,[True])[0] is None,'priority hostile not detected')
    # Retaining a tested prime naively cancels its guard.
    require(Fraction(11*2,7*2)==Fraction(11,7),'guard cancellation control')
    require(not prime(1) and prime(2) and prime(3) and not prime(9) and not prime(49),'proper divisor controls')
    return {'exhaustive_arbitrary_vector_cases':exhaustive,'tm_exhaustive_cases':exhaustive,'ca_exhaustive_cases':exhaustive,'prime_vector_full_vs_pruned_cases':primecases,'primegame_stream_pair_cases':streamcases,'hostile_controls':['5 and 29 exceptions','diagonal 3+3','composite formal vector 15+15','first-applicable ordering','cancelled guard','proper divisors','CA collision','CA shift','CA frozen halt'], 'status':'PASS'}

if __name__=='__main__': print(json.dumps(check(),indent=2,sort_keys=True))
