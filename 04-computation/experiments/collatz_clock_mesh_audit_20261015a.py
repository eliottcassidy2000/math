#!/usr/bin/env python3
"""Independent exact checks of the arithmetic envelope and integer-clock census.
Usage: python3 this.py census.json report.json
All universal claims are proved in the companion note, not inferred here.
"""
import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path

def need(b, message):
    if not b:
        raise RuntimeError(message)

def envelope(a, L):
    m,u=divmod(L,3)
    off=((0,0,1),(0,1,2),(0,1,1))[a%3]
    j1=m*(m-1)//2
    j2=m*(m-1)*(2*m-1)//6
    s=6*j1+m*sum(off)+2*m*u+sum(off[:u])
    q=12*j2+4*j1*sum(off)+m*sum(x*x for x in off)+4*m*m*u+4*m*sum(off[:u])+sum(x*x for x in off[:u])
    return F(L*q-s*s,L)

def variance(xs):
    return F(sum(x*x for x in xs))-F(sum(xs)**2,len(xs))

def label_checks():
    count=0
    for a in range(3):
        for L in range(1,81):
            c=[sum((a+h)%3!=0 for h in range(i)) for i in range(L)]
            need(envelope(a,L)==variance(c),"counter formula")
            need(envelope(a,L)<=F(L**3,27)+F(L,3),"uniform envelope")
        for L in range(1,9):
            bound=envelope(a,L)
            for labels in itertools.product((-1,0,1),repeat=L):
                J2=[0]
                for i in range(L-1):
                    w=(0,-1,1)[(a+i)%3]
                    J2.append(J2[-1]+w*(labels[i+1]-labels[i]))
                s=sum(J2); q=sum(x*x for x in J2)
                need(F(L*q-s*s,4*L)<=bound,"variance bound for typed labels")
                for h in range(1,L):
                    innovation=F((L*sum(J2[:h])-h*s)**2,4*L*h*(L-h))
                    slack=bound-envelope(a,h)-envelope(a+h,L-h)
                    need(innovation<=slack,"upper certificate decreased on a split")
                    count+=1
            # Sharp control is a real binary-label sequence, not arbitrary increments.
            labels=[1 if (a+i)%3==2 else -1 for i in range(L)]
            current=[0]
            for i in range(L-1):
                current.append(current[-1]+F((0,-1,1)[(a+i)%3]*(labels[i+1]-labels[i]),2))
            need(variance(current)==bound,"sharp periodic extremizer")
    return count

def literal_control(rows):
    N=min(1<<14,2*rows[-1]['X'])
    labels=[0]*(N+1)
    for n in range(1,N+1):
        x=n;steps=0
        while x!=1:
            x=(3*x+1)//2 if x%2 else x//2
            steps+=1
            need(steps<1000000,"literal first-hit cap")
        labels[n]=1 if steps%2==0 else -1
    H=J=0
    for row in rows:
        X=row['X']
        if 2*X>N:
            break
        values=[];energy=0
        need(J==row['incoming_current'],"literal incoming current")
        for n in range(X,2*X):
            H+=(2 if n%3==2 else -1)*labels[n]
            need(H==(0,-1,1)[n%3]*labels[n]-2*J,"literal summation by parts")
            values.append(J);energy+=H*H
            J+=F((0,-1,1)[n%3]*(labels[n+1]-labels[n]),2)
        need(sum(values)==row['current_sum'],"literal source sum")
        need(sum(x*x for x in values)==int(row['current_square_sum']),"literal source square")
        need(energy==int(row['H_energy']),"literal discrepancy energy")
    return N

def main():
    census=json.loads(Path(sys.argv[1]).read_text());rows=census['octaves']
    count=label_checks();literal_N=literal_control(rows)
    root=Path(__file__).resolve().parents[2]
    old=json.loads((root/'05-knowledge/results/collatz_break_clock_20261013d.json').read_text())
    previous={r['k']:r for r in old['octaves']}
    prior=json.loads((root/'05-knowledge/results/collatz_discrepancy_20261013c.json').read_text())
    # This earlier file uses (X,2X] rather than [X,2X).
    prior_rows={r['k']:r for r in prior['rows']}
    output=[];compared=0
    for row in rows:
        k=row['k'];X=row['X'];S=int(row['current_square_sum']);E=int(row['H_energy'])
        if k in previous:
            for field in ('incoming_current','current_sum','current_square_sum','H_energy'):
                need(int(row[field])==int(previous[k][field]),'prior break-clock '+field)
            compared+=1
        lo=prior_rows[k]['residues'][2][0]
        hi=prior_rows[k+1]['residues'][2][0]
        need(E==int(prior_rows[k+1]['residues'][2][3])+lo*lo-hi*hi,"prior discrepancy endpoints")
        need(F(row['current_sum']**2,X)+sum(F(int(num),2**i) for i,num in enumerate(row['haar_numerators'],1))==S,"JSON Haar identity")
        meshes=[];last=None
        for mesh in row['meshes']:
            b=mesh['b'];C=F(int(mesh['coarse_numerator']),mesh['coarse_denominator'])
            blocks=X//b
            counts=[blocks//3]*3
            for h in range(blocks%3):
                counts[(X+h*b)%3]+=1
            R=sum(counts[r]*envelope(r,b) for r in range(3))
            Z=C+R
            need(C<=S<=Z,"mesh brackets actual energy")
            need(R<=F(X*b*b,27)+F(X,3),"mesh tail bound")
            if last is not None:
                need(last[0]<=C and Z<=last[1],"nested mesh certificates")
            last=(C,Z)
            A=2*X*X-4*Z-X
            paid=A>=0 and A*A>=16*X*Z
            meshes.append({'b':b,'coarse_over_X2':float(C/(X*X)),
                'upper_over_X2':float(Z/(X*X)), 'upper_exact':[Z.numerator,Z.denominator],
                'certifies_E_le_2X2':paid})
        output.append({'k':k,'S_over_X2':S/(X*X),'E_over_X2':E/(X*X),
                       'mean_over_X2':row['current_sum']**2/X**3,'meshes':meshes})
    report={'status':'FINITE-EXACT; no assertion beyond the stated universe',
            'source_limit':census['source_limit'],'typed_label_split_checks':count,
            'literal_first_hit_sources':literal_N,'previous_break_octaves_matched':compared,
            'previous_discrepancy_octaves_matched':len(rows),'octaves':output}
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='octaves'},indent=2))

if __name__=='__main__':
    main()
