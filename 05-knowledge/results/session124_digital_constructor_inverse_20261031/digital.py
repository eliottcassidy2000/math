"""Exact digital normal forms for actual marked singleton/triangle constructors.

No function claims that an arbitrary integer state is a realizable seed.
"""
from fractions import Fraction


def singleton(state, bit):
    h,F,q,g=state
    if bit==0:return h,F+q+g,2*q,0
    if bit==1:return h,F,2*q+h,h
    raise ValueError('Bit must be zero or one')


def triangle(state,j):
    h,F,q,g=state
    if j not in range(4):raise ValueError('Triangle cut must be 0..3')
    return 3*h,3*F+(9-3*j)*q+(3-j)*g+(2*h if j in (1,2) else 0),12*q+3*j*h,j*h


def f(n):
    if n<0:raise ValueError('Negative digital address')
    value=1
    while n:
        if n%2:n//=2
        else:
            n//=2;value+=n+n%2
    return value


def f_pair_formula(n):
    bits=[(n>>i)&1 for i in range(n.bit_length())]
    cross=sum((1<<(i-j-1)) for i,z in enumerate(bits) if z for j in range(i) if bits[j])
    adjacent=sum(bits[i]*bits[i-1] for i in range(1,len(bits)))
    return 1+n-(n&1)-cross-adjacent


def address(word):
    return sum((1-bit)<<i for i,bit in enumerate(word))


def output(state,word):
    h,F,q,g=state; B=address(word)
    return F+q*B+g*(B&1)+h*(f(B)-1)


def radical_inverse(n):
    return sum((Fraction(1,1<<(i+1)) for i in range(n.bit_length()) if n>>i&1),Fraction())


def error_runs(n):
    """Exact f(n)-n(1-v(n)); n=0 is handled separately."""
    if n==0:return Fraction(1)
    bits=[(n>>i)&1 for i in range(n.bit_length())]
    total=Fraction(1-bits[0]);i=0
    while i<len(bits):
        if not bits[i]:i+=1;continue
        start=i
        while i<len(bits) and bits[i]:i+=1
        end=i-1; length=i-start
        R=sum((Fraction(1,1<<(j-end)) for j in range(end+1,len(bits)) if bits[j]),Fraction())
        total+=R+Fraction(1,1<<length)*(1-R)
    return total


def predecessor(h,q):
    if h<=0 or h%2!=1 or q<0:raise ValueError('Require odd positive h and q>=0')
    if q==0:return None  # delete the exact zero self-loop
    if q%2:
        if q<h:return None
        return (q-h)//2,1
    return q//2,0


def phi(h,q):
    value=0
    while True:
        previous=predecessor(h,q)
        if previous is None:return value
        q,bit=previous
        if bit==0:value+=q+h*(q%2)


def kappa_sources(state):
    """Actual first-child nodes, not arbitrary positive substitutions."""
    h,F,q,g=state
    return [(singleton(state,b)[2], F-phi(h,q)+(g-h*(q%2) if b==0 else 0),b) for b in (0,1)]


def ancestor_word(h,ancestor,node):
    """Return unique reduced singleton word, or None if no ancestry exists."""
    reverse=[]
    while node!=ancestor:
        previous=predecessor(h,node)
        if previous is None:return None
        node,bit=previous;reverse.append(bit)
    return tuple(reversed(reverse))


def fiber(h,q,roots):
    """Exact F set from nonempty singleton words on the supplied actual roots."""
    values=set()
    for root in roots:
        if root[0]!=h:continue
        for node,kappa,bit in kappa_sources(root):
            if ancestor_word(h,node,q) is not None:values.add(phi(h,q)+kappa)
    return values


def inverse_addresses(target,state):
    """Exact finite scalar inverse for an actual seed.

    Return canonical B addresses. Final singleton 1s are discarded; when q=0,
    initial neutral 0s are also discarded. The q=0 condition is validated as
    the actual empty-cut state F=h,g=0 rather than relaxed to a formal state.
    """
    h,F,q,g=state
    if min(h,F)<1 or min(q,g)<0 or g>h:raise ValueError('Invalid state')
    if target<F:return []
    if q:
        candidates=range((target-F)//q+1)
    else:
        if (F,g)!=(h,0):raise ValueError('q=0 actual seed must have F=h,g=0')
        if target%h:return []
        candidates=range(0,2*(target//h-1)+1,2)
    return [B for B in candidates if F+q*B+g*(B%2)+h*(f(B)-1)==target]


def macro_decodings(state,Q,triangles):
    """All legal nonconstant adjacent-macro decodings at one shadow q.

    triangles is 0, 1, or 2. Every returned block word starts from the supplied
    actual seed. A q=0 seed needs one optional leading zero to represent a
    starting 01 macro whose zero is invisible in the reduced q ancestry.
    """
    if triangles not in (0,1,2):raise ValueError('Require 0,1,2 triangles')
    h,F,q,g=state
    word=ancestor_word(h,q,Q)
    if word is None:return []
    words=[word]
    if q==0:
        if (F,g)!=(h,0):raise ValueError('Invalid actual q=0 seed')
        words.append((0,)+word)
    out=[]
    for bits in words:
        if triangles==0:
            out.append(bits);continue
        for pos in range(len(bits)-2*triangles+1):
            pairs=[bits[pos+2*j:pos+2*j+2] for j in range(triangles)]
            if any(a==b for a,b in pairs):continue
            blocks=bits[:pos]+tuple('C'+str(2*a+b) for a,b in pairs)+bits[pos+2*triangles:]
            out.append(blocks)
    return out


def run_blocks(state,blocks):
    for block in blocks:
        state=triangle(state,int(block[1])) if isinstance(block,str) else singleton(state,block)
    return state


def inverse_non3(target,state):
    """Complete finite inverse for this seed's non-3 singleton/triangle language.

    This is an exact decision algorithm, not a uniform coverage theorem.
    It deliberately uses a final-state bound, never the invalid q>N pruning
    inside arbitrary block words. It may be slow for large targets.
    """
    if target<1 or target%3==0:raise ValueError('Require positive target not divisible by 3')
    if state[1]>target:return None
    if state[1]==target:return ()
    for r in range(3):
        scale=3**r
        for Q in range((4*target-1)//scale+1):
            for blocks in macro_decodings(state,Q,r):
                if run_blocks(state,blocks)[1]==target:return blocks
    return None
