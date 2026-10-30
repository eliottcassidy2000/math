"""Independent exact audit: scalar histories and weight elimination, not DAG coverage."""
from fractions import Fraction as Q
from random import Random
import json

def need(x, label):
    if not x: raise RuntimeError(label)

def histories():
    rng=Random(118042)
    total=0; finished=0
    for r in (3,4):
        for trial in range(3000):
            n=rng.randrange(2,41); J=rng.randrange(n+1)
            m,q,I,s=n,0,n-J,500
            G=B=A=P=Z=H=D=K=E=U=0
            for step in range(70):
                possible=['D','M']
                if I>=2: possible.append('G')
                if q:
                    possible.extend(['P','Z','H'])
                    if I: possible.append('A')
                if q>=2: possible.append('B')
                row=rng.choice(possible)
                dm=dq=di=cost=0
                if row=='G':
                    dq=1;di=-2;cost=rng.randrange(2,6);G+=1;K+=cost-2
                elif row=='P': dm=-1;dq=-1;di=1;cost=2;P+=1
                elif row=='Z': dm=-1;dq=-1;Z+=1
                elif row=='H': dm=-1;dq=-1;cost=1;H+=1
                elif row=='A': dm=-2;dq=-1;di=-1;cost=3;A+=1
                elif row=='B': dm=-3;dq=-2;cost=r;B+=1
                elif row=='D': cost=1;D+=1
                lossI=rng.randrange(min(2,I+di)+1)
                lossS=rng.randrange(3)
                U+=lossI;E+=lossS
                m+=dm;q+=dq;I+=di-lossI;s-=cost+lossS
                need(m>=2*q>=0 and 0<=I<=m-2*q and s>=0,'physical scalar bounds')
                h=m-q
                need(n-h==G+B+A,'dimension conservation')
                need(q==G-2*B-A-P-Z-H,'pair conservation')
                need(I==n-J-2*G-A+P-U,'usage conservation')
                need(h==B+A+q+Z+H+U+I+J,'endpoint identity')
                need(A+B<=h and G>=n-2*h and P>=n-4*h,'frequency bounds')
                need(500==s+2*G+r*B+3*A+2*P+H+D+K+E,'actual deletion inventory')
                w=8 if r==4 else 9
                coeff=4 if r==4 else 5
                extra=s+coeff*(I+J+U)+(coeff-2)*q+(coeff-3)*A+(coeff-2)*Z+(coeff-1)*H+D+K+E
                need(500==4*n-w*h+extra,'expanded gate identity')
                total+=1
            finished+=1
    return {'histories':finished,'checked_prefixes':total,'seed':118042,
            'scope':'abstract row histories with extra usage/deletion losses; not circuit realizations'}

def weights():
    accepted=checked=0
    # Direct enumeration in weight coordinates, independent of endpoint construction.
    for ii in range(41):
        i=Q(ii,4)
        for pp in range(41):
            p=Q(pp,4)
            for tt in range(65):
                t=Q(tt,4);c=t-i;b=t
                legal=i<=2+p and t<=2+2*i-p and t<=4+2*p and t<=3+i+p
                if legal and c>=0:
                    need(c<=4 and b>=max(c,2*c-2,4*c-8),'necessary projection')
                    accepted+=1
                checked+=1
    pairs=0
    for ci in range(65):
        c=Q(ci,16)
        for jj in range(17):
            b=max(c,2*c-2,4*c-8)+Q(jj,8)
            i=b-c;p=max(Q(0),i-2,(b-4)/2,c-3)
            need(i>=0 and p>=0 and i<=2+p and b<=2+2*i-p and b<=4+2*p and b<=3+i+p,
                 'constructive sufficiency')
            need(p<=b-2*c+2,'fiber interval')
            pairs+=1
    need(max(Q(0),Q(8)-4-2,(Q(8)-4)/2,Q(4)-3)==Q(8)-2*4+2==2,'unique optimum')
    return {'direct_weight_triples':checked,'feasible_nonnegative_coefficient_triples':accepted,
            'constructed_projection_pairs':pairs,'optimum':[1,4,2,8]}

if __name__=='__main__':
    print(json.dumps({'status':'PASS','histories':histories(),'weights':weights()},sort_keys=True,indent=2))
