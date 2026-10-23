"""Independent exact algebra audit of the now complete degree-three series."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
GATES=0


def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise ValueError(why)


def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:p.pop()
    return p


def add(*ps):
    out=[F(0)]*max(map(len,ps))
    for p in ps:
        for i,a in enumerate(p):out[i]+=a
    return trim(out)


def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]+=a*b
    return trim(out)


def scale(p,c):return trim([c*a for a in p])
def shift(p,k):return trim([0]*k+list(p))


def divrem(p,q):
    p,q=trim(p),trim(q)
    if q==[0]:raise ValueError('zero divisor')
    quot=[F(0)]*max(1,len(p)-len(q)+1)
    while p!=[0] and len(p)>=len(q):
        k=len(p)-len(q);c=p[-1]/q[-1];quot[k]+=c
        p=add(p,scale(shift(q,k),-c))
    return trim(quot),p


def xgcd(p,q):
    p,q=trim(p),trim(q);s0,s1=[F(1)],[F(0)];t0,t1=[F(0)],[F(1)]
    while q!=[0]:
        a,r=divrem(p,q)
        p,q=q,r;s0,s1=s1,add(s0,scale(mul(a,s1),-1));t0,t1=t1,add(t0,scale(mul(a,t1),-1))
    c=1/p[-1]
    return scale(p,c),scale(s0,c),scale(t0,c)


def convolve(p,q,n):
    return [sum(p[j]*q[i-j] for j in range(i+1)) for i in range(n+1)]


def run():
    global GATES
    GATES=0
    p=[1,-1,0,-1];u=[1,0,-1,0,0,0,-1];p2=mul(p,p);u2=mul(u,u)
    R=[0]*5+[1,1,1]
    numerator=add(scale(shift(add(u2,scale(mul(p2,u),-1)),5),F(1,2)),
                  mul(mul(p2,u),R),shift(mul(p2,[1,1,0,1]),10))
    denominator=mul(p2,u2)
    expected=[0]*19
    for i,a in {5:1,7:-2,8:-2,9:1,10:3,14:-1,15:2,16:-1,17:1,18:-1}.items():expected[i]=a
    need(numerator==expected,'independent common-denominator numerator')
    need(len(denominator)==19 and denominator[0]==1,'normalized degree18 denominator')
    g,s,t=xgcd(numerator,denominator)
    need(g==[1],'no numerator/denominator cancellation')
    need(add(mul(s,numerator),mul(t,denominator))==[1],'exact Bezout certificate')
    need(all(a.denominator==1 for a in numerator+denominator),'integer exact coefficients')
    numerator=list(map(int,numerator));denominator=list(map(int,denominator))
    nmax=256
    A=[0]*(nmax+1);A[0]=1
    for n in range(1,nmax+1):A[n]=A[n-1]+(A[n-3] if n>=3 else 0)
    A2=[A[n//2] if n%2==0 else 0 for n in range(nmax+1)]
    Asq=convolve(A,A,nmax);A2sq=convolve(A2,A2,nmax)
    direct=[]
    pieces=[]
    for n in range(nmax+1):
        difference=Asq[n-5]-A2[n-5] if n>=5 else 0
        need(difference%2==0 and difference>=0,'nonselfconverse h3 pairs are integral')
        nonself=difference//2
        rigid=sum(A2[n-j] for j in (5,6,7) if n>=j)
        nonrigid=sum(A2sq[n-j] for j in (10,11,13) if n>=j)
        pieces.append((nonself,rigid,nonrigid));direct.append(nonself+rigid+nonrigid)
    divided=[]
    for n in range(nmax+1):
        divided.append((numerator[n] if n<len(numerator) else 0)-sum(denominator[j]*divided[n-j] for j in range(1,min(n,18)+1)))
        need(divided[n]==direct[n],'independent series decomposition equals rational division')
        defect=sum(denominator[j]*direct[n-j] for j in range(min(n,18)+1))
        need(defect==(numerator[n] if n<len(numerator) else 0),'full inhomogeneous recurrence')
        if n>=19:need(defect==0,'homogeneous recurrence from19')
    need(sum(denominator[j]*direct[18-j] for j in range(19))==-1,'homogeneous recurrence fails at18')
    need(direct[15:22]==[114,175,276,425,665,1019,1583],'requested indexed prefix')
    old=json.loads((HERE.parent/'deficit_six_centers_20261021b/receipt.json').read_text())
    need([(len(c['source']),c['h'],c['b']) for c in old['centers']]==[(5,9,3),(6,15,9),(7,9,3)],'literal inherited three centers and deficits')
    pinpaths=(
      '05-knowledge/results/metagraph_classification_20261021a/README.md',
      '05-knowledge/results/deficit_six_centers_20261021b/receipt.json',
      '05-knowledge/results/center_incidence_20261022/receipt.json',
      '05-knowledge/results/nontransitive_fifteen_20261022/receipt.json',
      '05-knowledge/results/nontransitive_fifteen_20261022/receipt21.json')
    pins={r:hashlib.sha256((ROOT/r).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for r in pinpaths}
    fractions=lambda a:[str(x) for x in a]
    return {'schema':'complete-degree-three-gf-audit-20261022','checks':GATES,
      'coefficient_order':'ascending powers of z','numerator':numerator,'denominator':denominator,
      'bezout_numerator_multiplier':fractions(s),'bezout_denominator_multiplier':fractions(t),
      'minimal_eventual_recurrence_order':18,'homogeneous_onset':19,'defect_at18':-1,
      'prefix_0_through40':direct[:41],'three_disjoint_contributions_0_through25':pieces[:26],
      'comparison_through':nmax,'inherited_pins':pins,
      'scope':'Degree-three vertices of the fundamental converse-merged paired-black graph Q; not A000568, labelled tournaments, or metagraph wiring.'}


def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group();g.add_argument('--write',action='store_true');g.add_argument('--verify',action='store_true');a=parser.parse_args()
    data=(json.dumps(run(),sort_keys=True,indent=2)+'\n').encode()
    if a.write:(HERE/'gf_receipt.json').write_bytes(data)
    if a.verify and data!=(HERE/'gf_receipt.json').read_bytes().replace(b'\r\n',b'\n'):raise ValueError('saved comparison')
    print(f'{GATES:,} exact checks; LF SHA256 {hashlib.sha256(data).hexdigest()}')


if __name__=='__main__':main()
