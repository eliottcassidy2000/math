#!/usr/bin/env python3
"""Exact two-apex Hamiltonian probes and actual-word recovery; stdlib only.

The 36-state gap transfer is inherited, not claimed new. A word letter c has
x->v iff c&1 and y->v iff c&2; the marked edge is x->y.
"""
from itertools import product

STATES = tuple(product(range(2), range(3), range(3), range(2)))
INDEX = {s:i for i,s in enumerate(STATES)}
D = len(STATES)
SUFFIXES = ('0','1','2','00','10','20','30','01','000','100','200','010')
X_NUMERATOR = (-5,-3,6,15,-9,-18,6,9,-4,6,6,-9)
Y_NUMERATOR = (-1,3,0,3,-3,0,0,-3,-2,0,0,3)
DENOMINATORS = (36,12)
FULL_SUFFIXES = ('','0','1','2','00','10','20','30','01','11','21','31',
                 '000','100','200','010','210')


def valid_word(word):
    if not isinstance(word,str) or any(c not in '0123' for c in word):
        raise ValueError('word must be a string over 0123')
    return word


def transitions(letter):
    edges=[]
    for j,(order,p,q,r) in enumerate(STATES):
        a=(letter>>(1 if order else 0))&1
        b=(letter>>(0 if order else 1))&1
        edges.append((INDEX[order,2 if a else 1,q,r],j))
        if q or a: edges.append((INDEX[order,p,2 if b else 1,r],j))
        if r or b: edges.append((INDEX[order,p,q,1],j))
    return tuple(edges)

TRANSITIONS=tuple(transitions(a) for a in range(4))
INITIAL=tuple(int(s in ((0,0,0,0),(1,0,0,0))) for s in STATES)
OUTPUT=tuple(int(p<2 and q<2 and (order==0 or q==1))
             for order,p,q,r in STATES)


def step(state,letter):
    if type(letter) is not int or not 0<=letter<4:
        raise ValueError('letter must be an integer from 0 to 3')
    if len(state)!=D: raise ValueError('wrong state dimension')
    result=[0]*D
    for i,j in TRANSITIONS[letter]: result[i]+=state[j]
    return tuple(result)


def adjoint(row,letter):
    result=[0]*D
    for i,j in TRANSITIONS[letter]: result[j]+=row[i]
    return tuple(result)


def state(word):
    result=INITIAL
    for c in valid_word(word): result=step(result,int(c))
    return result


def suffix_row(suffix):
    result=OUTPUT
    for c in valid_word(suffix)[::-1]: result=adjoint(result,int(c))
    return result

PROBE_ROWS=tuple(suffix_row(s) for s in SUFFIXES)
FULL_ROWS=tuple(suffix_row(s) for s in FULL_SUFFIXES)


def dot(a,b): return sum(x*y for x,y in zip(a,b))


def hamilton(word): return dot(OUTPUT,state(word))


def probes(word):
    v=state(word)
    return tuple(dot(r,v) for r in PROBE_ROWS)


def full_profile(word):
    v=state(word)
    return tuple(dot(r,v) for r in FULL_ROWS)


def extract_codes(values):
    if len(values)!=len(SUFFIXES) or any(type(v) is not int for v in values):
        raise ValueError('expected twelve exact integer probe values')
    codes=[]
    for numer,den in zip((X_NUMERATOR,Y_NUMERATOR),DENOMINATORS):
        n=dot(numer,values)
        if n<0 or n%den:
            raise ValueError('probe values do not have nonnegative integral codes')
        codes.append(n//den)
    return tuple(codes)


def interleave_codes(qx,qy):
    if any(type(q) is not int or q<0 for q in (qx,qy)):
        raise ValueError('codes must be nonnegative integers')
    m=max(qx.bit_length(),qy.bit_length())
    return ''.join(str(((qx>>j)&1)+2*((qy>>j)&1)) for j in range(m-1,-1,-1))


def recover_word(values):
    """Membership iff reconstruction replays exactly; integral codes alone do not suffice."""
    word=interleave_codes(*extract_codes(values))
    if probes(word)!=tuple(values):
        raise ValueError('integer probe vector is not an actual word profile')
    return word


def recover_full_profile(values):
    if len(values)!=len(FULL_SUFFIXES) or any(type(v) is not int for v in values):
        raise ValueError('expected seventeen exact integer continuation values')
    subset=tuple(values[FULL_SUFFIXES.index(s)] for s in SUFFIXES)
    word=interleave_codes(*extract_codes(subset))
    if full_profile(word)!=tuple(values):
        raise ValueError('integer continuation vector is not an actual word state')
    return word


def direct_codes(word):
    qx=qy=0
    for c in valid_word(word):
        a=int(c);qx=2*qx+(a&1);qy=2*qy+((a>>1)&1)
    return qx,qy


def bounded_spine_length(h):
    """Cited Busch/SCC bound: largest m with 5**m <= h**3 (no floats)."""
    if type(h) is not int or h<1: raise ValueError('h must be a positive integer')
    power=1;m=0;limit=h*h*h
    while power*5<=limit:power*=5;m+=1
    return m


def find_count(h):
    """Complete all-order decision within the TWO-APEX language, not all tournaments.

    Potentially expensive: at most 4**bounded_spine_length(h) leaf words.
    Returns a normalized actual word, or None if none exists in this language.
    """
    depth=bounded_spine_length(h)
    if h%2==0:return None
    stack=[(INITIAL,'')]
    while stack:
        v,w=stack.pop()
        count=dot(OUTPUT,v)
        if count==h:return w
        if count>h or len(w)==depth:continue
        for a in range(3,(1 if not w else 0)-1,-1):
            stack.append((step(v,a),w+str(a)))
    return None
