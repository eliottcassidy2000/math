"""Constructed ROOT ports with source-owned membership and bounded precision.

The factory constructs ordinary native sources. It never substitutes one for
a requested Mersenne. No production function discovers a ROOT orbit.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F

import collatz_child_normal_forms_20261007 as forms
import collatz_uncovered_join_routes_20261007 as routes
import collatz_parameter_cover_20261007e as positive

CHECKS=0
DEPTH_CAP=2048


def need(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(message)


def integer(n,minimum=0):
    need(type(n) is int and n>=minimum,'exact integer domain')
    return n


def head(word):
    routes.letters(word)
    need(not word or word[0]>=2,'head follows the maximal initial ones')
    return word


def invariant(word):
    p,q,b=routes.carrier(head(word))
    return 3*b-3*p+q


@dataclass(frozen=True)
class Port:
    drop:int
    left:tuple
    right:tuple
    gap:int


def audit_port(port):
    need(type(port) is Port,'exact Port')
    integer(port.drop,1)
    need(type(port.gap) is int,'exact signed gap')
    head(port.left);head(port.right)
    need(len(port.right)-len(port.left)==port.drop,'odd-clock difference retained')
    need(sum(port.left)-sum(port.right)==2*port.gap,'signed sibling reserve retained')
    need(invariant(port.left)==invariant(port.right),'full normalized affine carry retained')
    return port


@dataclass(frozen=True)
class Completion:
    head:tuple
    run:int
    terminal:int


def terminal_equation(word,run):
    """Uninstantiated equation data, without constructing the huge 3-power."""
    head(word);integer(run)
    return sum(word),invariant(word),run+len(word)+1


def terminal_phase_prefix(word,run,precision):
    """One compatible finite phase prefix; this alone is NOT a ROOT proof."""
    A,I,depth=terminal_equation(word,run)
    integer(precision,1)
    need(precision<=min(depth,DEPTH_CAP),'declared prefix precision')
    modulus=3**precision
    return forms.two_log_three_power(I*pow(2,-A,modulus)%modulus,precision)


def parameters(item):
    need(type(item) is Completion,'exact Completion')
    head(item.head);integer(item.run);integer(item.terminal,4)
    depth=item.run+len(item.head)+1
    need(depth<=DEPTH_CAP,'declared finite ternary precision cap')
    A=sum(item.head);I=invariant(item.head);P=3**depth
    need(item.terminal%2==0,'ROOT exit has even valuation')
    need(I<=0 or A+item.terminal>=I.bit_length(),'explicit positive numerator')
    need(pow(2,A+item.terminal,P)==I%P,'exact completed-source ternary phase')
    need(I%4==2,'maximal initial run is retained')
    return A,I,P,depth


def compile_completion(word,run,index=0,minimum_terminal=4):
    head(word);integer(run);integer(index);integer(minimum_terminal,4)
    depth=run+len(word)+1
    need(depth<=DEPTH_CAP,'phase search exceeds declared precision cap')
    A=sum(word);I=invariant(word);P=3**depth
    target=I*pow(2,-A,P)%P
    c,period=forms.two_log_three_power(target,depth)
    cut=max(minimum_terminal,abs(I).bit_length()+2-A)
    c+=max(0,(cut-c+period-1)//period)*period
    item=Completion(word,run,c+period*index)
    parameters(item)
    return item


def shifted(item,index):
    _,_,_,depth=parameters(item);integer(index)
    result=Completion(item.head,item.run,item.terminal+2*3**(depth-1)*index)
    parameters(result)
    return result


def attach(item,port):
    parameters(item);audit_port(port)
    need(port.left==item.head,'same source head; no unrelated proof substituted')
    need(item.run>=port.drop,'child keeps a nonnegative initial run')
    need(item.terminal+2*port.gap>=4,'child ROOT exit has sufficient terminal reserve')
    return port


def root_word(item,port=None):
    parameters(item)
    if port is None:return (1,)*item.run+item.head+(item.terminal,)
    attach(item,port)
    return (1,)*(item.run-port.drop)+port.right+(item.terminal+2*port.gap,)


def source_mod(item,modulus,port=None):
    A,I,P,_=parameters(item);integer(modulus,1)
    drop=0 if port is None else attach(item,port).drop
    numerator=(pow(2,A+item.terminal,P*modulus)-I)%(P*modulus)
    need(numerator%P==0,'full denominator retained before division')
    return (pow(2,item.run-drop,modulus)*(numerator//P)-1)%modulus


def bit_bounds(item,port=None):
    A,I,P,_=parameters(item)
    drop=0 if port is None else attach(item,port).drop
    if A+item.terminal<=abs(I).bit_length()+2:
        n=(1<<(item.run-drop))*(((1<<(A+item.terminal))-I)//P)-1
        need(n>1,'finite boundary has a positive nonroot source')
        return n.bit_length(),n.bit_length()
    center=A+item.terminal+item.run-drop-P.bit_length()
    return max(1,center-2),max(1,center+2)


def materialize(item,port=None,bit_cap=20000):
    integer(bit_cap,1)
    need(bit_bounds(item,port)[1]<=bit_cap,'explicit literal expansion cap')
    A,I,P,_=parameters(item)
    drop=0 if port is None else attach(item,port).drop
    numerator=(1<<(A+item.terminal))-I
    need(numerator%P==0,'exact numerator division')
    n=(1<<(item.run-drop))*(numerator//P)-1
    routes.odd(n)
    need(n>1,'constructed first-hit source is not ROOT')
    return n


def audit_symbolic_root(item,port=None):
    """Check both coefficients in R=2^terminal, without expanding R."""
    A,I,P,_=parameters(item)
    alpha=F(1<<(item.run+A),P)
    beta=-F((1<<item.run)*I,P)-1
    if port is None:
        prefix=(1,)*item.run+item.head
        target=F(1,3)
    else:
        attach(item,port)
        alpha/=1<<port.drop;beta=(beta+1)/(1<<port.drop)-1
        prefix=(1,)*(item.run-port.drop)+port.right
        target=F(4)**port.gap/3
    p,q,b=routes.carrier(prefix)
    need(p*alpha/q==target and (p*beta+b)/q==F(-1,3),
         'symbolic prefix endpoint is exactly the declared ROOT predecessor')
    return root_word(item,port)


def recognizes_source(item,source,port=None):
    """Exact source identity check. A new family member is never substituted."""
    routes.odd(source)
    A,I,P,_=parameters(item)
    drop=0 if port is None else attach(item,port).drop
    scale=1<<(item.run-drop)
    if (source+1)%scale:return False
    target=P*((source+1)//scale)+I
    return target>0 and target&(target-1)==0 and target.bit_length()-1==A+item.terminal


def is_mersenne_member(item):
    """Exact cofactor-one test, with no expansion of the constructed source."""
    A,I,P,_=parameters(item)
    target=I+2*P
    return target>0 and target&(target-1)==0 and target.bit_length()-1==A+item.terminal


def try_mersenne_terminal(word,run):
    """Complete only this fixed source and head, or return None."""
    head(word);integer(run)
    depth=run+len(word)+1
    need(depth<=DEPTH_CAP,'declared finite phase precision cap')
    target=invariant(word)+2*3**depth
    if target<=0 or target&(target-1):return None
    c=target.bit_length()-1-sum(word)
    if c<4:return None
    item=Completion(word,run,c)
    try:parameters(item)
    except ValueError:return None
    need(is_mersenne_member(item),'same-source terminal identity')
    return item


def inverse_phase(item,word):
    """Unique index phase for a supplied contracting inverse return word."""
    parameters(item);routes.letters(word)
    p,q,b=routes.carrier(word)
    need(bool(word) and q<p,'strict contracting inverse word')
    target=b*pow(q,-1,p)%p
    residue=0;period=1
    for _ in word:
        modulus=period*3
        candidates=[residue+j*period for j in range(3)
                    if source_mod(shifted(item,residue+j*period),modulus)==target%modulus]
        need(len(candidates)==1,'ternary isometry gives one exact refinement')
        residue=candidates[0];period=modulus
    return residue,period


def inverse_root_word(item,index,word):
    """An illegal supplied index returns None, never a nearby legal index."""
    integer(index);routes.letters(word)
    p,q,b=routes.carrier(word)
    need(bool(word) and q<p,'strict inverse word')
    member=shifted(item,index)
    if (q*source_mod(member,p)-b)%p:return None
    return word+root_word(member)


def inverse_source_mod(item,index,word,modulus):
    integer(modulus,1)
    need(inverse_root_word(item,index,word) is not None,'same supplied index must be native')
    p,q,b=routes.carrier(word)
    numerator=q*source_mod(shifted(item,index),p*modulus)-b
    need(numerator%p==0,'inverse-child division precision')
    return numerator//p%modulus


def reject(fn,*args):
    try:fn(*args)
    except (ValueError,TypeError):need(True,'hostile rejected')
    else:need(False,'hostile accepted')


def main():
    # One phase grounds both sides without any discovered orbit as a seed.
    switch=Port(1,(),(2,),-1)
    tiny=compile_completion((),2)
    need(tiny.terminal==10,'closed-form first reset-family phase')
    need((materialize(tiny),materialize(tiny,switch))==(151,75),'small grounded pair')
    need(root_word(tiny)==(1,1,10) and root_word(tiny,switch)==(1,2,8),
         'explicit supplied first-hit ROOT words')
    examples=0
    for e in range(1,5):
        base=compile_completion((),e,minimum_terminal=6)
        for j in range(3):
            item=shifted(base,j)
            for port in (None,switch):
                word=audit_symbolic_root(item,port)
                n=materialize(item,port)
                need(forms.replay(n,word,root=True)==1,'independent strict literal ROOT replay')
                need(recognizes_source(item,n,port),'same supplied source recognized')
                need(not recognizes_source(item,n+2,port),'adjacent odd source is not replaced')
                lo,hi=bit_bounds(item,port)
                need(lo<=n.bit_length()<=hi,'literal expansion obeys symbolic bit bound')
                for m in (1,2,19,81,256,2187):
                    need(source_mod(item,m,port)==n%m,'independent quotient residue reader')
                examples+=1
    # A positive-gap nontrivial known macro, now with an explicit ROOT exit.
    small=Port(1,(2,6),(4,1,1),1)
    item=compile_completion(small.left,2)
    for port in (None,small):
        n=materialize(item,port)
        need(forms.replay(n,audit_symbolic_root(item,port),root=True)==1,'positive-gap literal completion')
    # The new coarse complement row has two genuine alternative child ports.
    bank=positive.load_rules();row=next(r for r in bank if r.seed_parameter==1)
    p4=Port(row.drop,row.source_head,row.partner_head,row.gap)
    need(p4.right[:2]==(2,2),'two alternative child clearings')
    p3=Port(3,p4.left,(4,)+p4.right[2:],p4.gap)
    completed=compile_completion(p4.left,2*sum(p4.left))
    for port in (None,p3,p4):
        audit_symbolic_root(completed,port)
        for m in (19,256,729,65537):
            residue=source_mod(completed,m,port)
            if port is not None:
                need((2**port.drop*(residue+1)-(source_mod(completed,m)+1))%m==0,
                     'same-source alternative child identities')
    need(bit_bounds(completed)[0]>10**20,'genuinely compressed astronomical grounded source')
    need(not is_mersenne_member(completed),'native-head completion is not the original Mersenne')
    reject(materialize,completed,None,10000)
    # Known huge Mersenne words are not silently asserted by this construction.
    for e in range(1,9):
        comp=compile_completion((),e)
        same=recognizes_source(comp,(1<<(e+1))-1)
        need(same==(e==1),'fixed-source membership remains an extra equation')
        need(is_mersenne_member(comp)==same,'symbolic cofactor-one check agrees')
    exact7=try_mersenne_terminal((2,3),2)
    need(exact7 is not None and materialize(exact7)==7 and
         root_word(exact7)==(1,1,2,3,4),'source-owned direct terminal succeeds when equation holds')
    need(try_mersenne_terminal((),2) is None,'fixedM7 cannot use the nearby151 completion')
    # Grounded inverse closure is now driven by these constructed proofs.
    basic=compile_completion((),0) # literal5 ->1
    need(materialize(basic)==5,'literal ROOT-family base5')
    for word in ((1,),(1,2),(1,1,2,2)):
        r,p=inverse_phase(basic,word)
        for j in range(2):
            index=r+p*j;member=shifted(basic,index)
            target=materialize(member,bit_cap=20000)
            pp,qq,bb=routes.carrier(word);child=(qq*target-bb)//pp
            proof=inverse_root_word(basic,index,word)
            need(proof is not None and forms.replay(child,proof,root=True)==1,
                 'inverse child consumes explicit grounded suffix')
            need(0<child<target,'inverse port is paid and grounded')
            for m in (19,81,256):
                need(inverse_source_mod(basic,index,word,m)==child%m,'inverse source exact residue')
        need(inverse_root_word(basic,(r+1)%p,word) is None,'no parameter substitution at illegal point')
    # Isometry on explicit small families, including source identity at each phase.
    for a in range(12):
        for b in range(a):
            x=materialize(shifted(basic,a));y=materialize(shifted(basic,b))
            need(forms.valuation(x-y,3)==forms.valuation(a-b,3),'exact ternary isometry')
    reject(compile_completion,(1,),2)
    reject(compile_completion,(),True)
    reject(compile_completion,(),DEPTH_CAP)
    reject(parameters,replace(tiny,terminal=11))
    reject(parameters,replace(tiny,run=2.0))
    reject(audit_port,replace(p4,drop=True))
    reject(audit_port,replace(p4,right=p4.right[:-1]+(2,)))
    reject(attach,tiny,p4)
    reject(attach,compile_completion((),1),switch) # child1 would be ROOT padding
    reject(source_mod,tiny,True)
    reject(recognizes_source,tiny,151.0)
    # Large new residual: retain the exact equation without silently selecting q.
    import collatz_residual_t23_20261007f as residual
    deep=residual.load()[-1].macro
    run=residual.BASE+residual.STRIDE*23-1
    A,I,depth=terminal_equation(deep.left,run)
    need((len(deep.left),A,depth)==(3317,6552,99708997022),'same large residual head')
    previous=(0,1)
    for precision in range(1,13):
        phase,period=terminal_phase_prefix(deep.left,run,precision)
        need((phase-previous[0])%previous[1]==0,'lazy terminal phases are compatible')
        need(pow(2,A+phase,3**precision)==I%3**precision,'exact finite equation prefix')
        previous=(phase,period)
    p,q,b=routes.carrier(deep.left);modulus=q*(1<<32)
    x=(2*pow(3,run,modulus)-1)%modulus
    numerator=p*x+b
    need(numerator%q==0,'actual large source has the native head')
    z=numerator//q%(1<<32)
    a=routes.v2(3*z+1)
    need(a==1 and pow(2,A+a,3)!=I%3,'actual next letter fails even the first ROOT-phase digit')
    reject(compile_completion,deep.left,run)
    print('PROVED grounded ordinary-source ports from an exact terminal exponent phase; no ROOT orbit discovery.')
    print('Literal pair151/75: supplied first-hit words(1,1,10) and(1,2,8).')
    print('Independent literal ROOT controls:',examples+2)
    print('Native head of coarse t1 rule: D3 and D4 ports grounded on the SAME constructed ordinary-source family.')
    print('Compressed terminal bit length:',completed.terminal.bit_length(),'source bit-length bounds:',bit_bounds(completed))
    print('PROVED index ternary isometry and unique inverse-port refinement; illegal original index remains rejected.')
    print('Fixed Mersenne requires additional equation2^(A+c)=I+2*3^(e+l+1); no billion-bit Mersenne is asserted grounded.')
    print('Actual t23:3317-letter head,required ternary depth99708997022; nextvaluation1 fails ROOT-phase mod3; finite prefixes remain OPEN.')
    print('Ternary phase precision cap:',DEPTH_CAP,'; no universal source coverage.')
    print('Exact checks:',CHECKS)


if __name__=='__main__':main()
