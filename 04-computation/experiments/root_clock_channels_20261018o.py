#!/usr/bin/env python3
"""Actual ROOT charges, exact clock coefficients, and source-specific readouts.

No randomized weights, orbit truncation, or modular-clock replacement. Integer
arithmetic throughout. Orthogonal clock channels are compared with the actual
Nyquist readout, not assumed interchangeable with it.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def edge(s):
    a = 2 if s % 8 == 1 else 1
    t = a
    v = (3*s + 1) >> a
    while v % 8 == 5:
        v = (v-1)//4
        t += 2
    return v, t


def charge(s, n):
    row = [0, 0, 0]
    m = s
    while m < 4*n:
        weights = (0, -3, 3) if m < 2*n else (-1, -1, 2)
        for r in range(3):
            row[r] += weights[(m+r) % 3]
        m = 4*m+1
    a, b = (row[0]-row[1])//3, -(row[0]+2*row[1])//3
    check(row == [2*a-b, -a-b, -a+2*b], 'integral source lift')
    check(norm((a,b)) in (0,1,3), 'source shell norm')
    return a, b


def norm(v):
    a,b = v
    return a*a-a*b+b*b


def ell(v):
    a,b = v
    return max(abs(a), abs(b), abs(a-b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--powers', type=int, default=18)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    here = Path(__file__).parent
    vpath = here/'first_arrival_20261018/V_closed_Y1000.txt'
    V = {int(line.split()[0]) for line in vpath.read_text().splitlines()
         if line.strip() and not line.startswith('#')}
    check(len(V) == 627 and all(edge(s)[0] in V for s in V), 'fixed frontier')
    memo = {s:(s,0,0,0) for s in V}

    def route(s):
        x=s; path=[]; seen=set()
        while x not in memo:
            if x in seen or len(path) >= 10000:
                raise RuntimeError(('UNRESOLVED route retained, stop scan',s,x))
            seen.add(x)
            y,t=edge(x); path.append((x,t)); x=y
        u,t,d,p=memo[x]
        for x,dt in reversed(path):
            if d==0:p=x
            t+=dt; d+=1; memo[x]=(u,t,d,p)
        return memo[s]

    rows=[]
    for n in sorted({1204,4096,5997,7910,11002} | {2**k for k in range(10,args.powers+1)}):
        q=defaultdict(lambda:[0,0]); bydepth=defaultdict(list)
        entry=defaultdict(lambda:[0,0])
        S=0; initial_mass=0; maxdepth=0; source_count=0
        example_sources=defaultdict(list)
        for s in range(1,4*n,2):
            if s%8 == 5:
                continue
            A=charge(s,n); u,t,d,p=route(s)
            S+=norm(A); initial_mass+=ell(A); source_count+=1
            maxdepth=max(maxdepth,d)
            for j in range(2):q[u,t][j]+=A[j]
            for j in range(2):entry[u,p][j]+=(-1)**t*A[j]
            if A != (0,0):
                bydepth[d].append((u,t,A))
                if n==1204:example_sources[u,t].append([s,*A])
        ports=defaultdict(dict)
        for (u,t),A in q.items():
            if A != [0,0]:ports[u][t]=A
        D=sum(norm(A) for A in q.values())
        energies={}; blocks={}
        # Exact Parseval energy for every clock residue modulo 2^k. This
        # is an AVERAGE over roots of unity, not a pointwise spectral bound.
        for modulus in (1,2,4,8,16,32,64,128,256,512,1024):
            bins=defaultdict(lambda:[0,0])
            for (u,t),A in q.items():
                for j in range(2):bins[u,t%modulus][j]+=A[j]
            blocks[str(modulus)]=sum(norm(A) for A in bins.values())
        for label,period,real,imag in [('zero',1,[1],[0]),
                                      ('pi',2,[1,-1],[0,0]),
                                      ('half_pi',4,[1,0,-1,0],[0,1,0,-1])]:
            # Summing the three real phase rows cancels the mixed sqrt(3)
            # terms; this is total three-phase energy divided by six.
            re=defaultdict(lambda:[0,0]); im=defaultdict(lambda:[0,0])
            for (u,t),A in q.items():
                for j in range(2):
                    re[u][j]+=real[t%period]*A[j]
                    im[u][j]+=imag[t%period]*A[j]
            energies[label]=sum(norm(A) for A in re.values())+sum(norm(A) for A in im.values())
        check(energies['zero']+energies['pi']==2*blocks['2'], 'two-channel Parseval')
        check(energies['zero']+energies['pi']+2*energies['half_pi']==4*blocks['4'], 'four-channel Parseval')
        check(blocks['1024']==D, '1024 resolves all tested clock differences')
        # Exact lag ledger: 2 Re(A conjugate B), no numerical complex roots.
        lag=defaultdict(int)
        quotient_mass=0; quotient_energy=0; derivative_energy=0; max_clock=0
        for u,poly in ports.items():
            ts=sorted(poly); lo,hi=ts[0],ts[-1]; max_clock=max(max_clock,hi)
            for i,t in enumerate(ts):
                a,b=poly[t]
                for v in ts[i+1:]:
                    c,d=poly[v]
                    lag[v-t]+=2*a*c-a*d-b*c+2*b*d
            # Q=z^lo*(W+(1+z)R), W=sum (-1)^(t-lo)q_t.
            previous=[0,0]; deriv=[0,0]
            for k in range(hi-lo,0,-1):
                coefficient=poly.get(k+lo,(0,0))
                current=[coefficient[j]-previous[j] for j in range(2)]
                quotient_mass+=ell(current); quotient_energy+=norm(current)
                for j in range(2):deriv[j]+=(-1)**(k-1)*current[j]
                previous=current
            derivative_energy+=norm(deriv)
        check(energies['zero']==D+sum(lag.values()), 'zero lag ledger')
        check(energies['pi']==D+sum((-1)**l*v for l,v in lag.items()), 'Nyquist lag ledger')
        # All native arrival prefixes, with source energy fixed at this N.
        partial=defaultdict(lambda:[0,0]); peak=(0,0)
        for d in range(maxdepth+1):
            for u,t,A in bydepth[d]:
                for j in range(2):partial[u][j]+=(-1)**t*A[j]
            value=sum(norm(A) for A in partial.values())
            if value>peak[0]:peak=(value,d)
        check(peak[0]>=energies['pi'], 'prefix maximum includes final')
        prefix_wiener=None
        if n in (4096,5997,11002,16384):
            poly=defaultdict(dict); corr=defaultdict(int); diag=0; best=(0,0)
            for d in range(maxdepth+1):
                for u,t,A in bydepth[d]:
                    a,b=A; old=poly[u].get(t,(0,0))
                    diag+=norm(A)+2*a*old[0]-a*old[1]-b*old[0]+2*b*old[1]
                    for v,(c,e) in poly[u].items():
                        if v!=t:corr[abs(v-t)]+=2*a*c-a*e-b*c+2*b*e
                    poly[u][t]=[old[0]+a,old[1]+b]
                budget=diag+sum(abs(v) for v in corr.values())
                if budget>best[0]:best=(budget,d)
            check(diag==D and dict(corr)==dict(lag),'independent prefix correlation ledger')
            prefix_wiener=list(best)
        exact_hostile=None
        if n==4096:
            z=(Fraction(-35,37),Fraction(12,37)); zz=(Fraction(1),Fraction(0))
            spectral=Fraction(D)
            for k in range(1,max(lag)+1):
                zz=(zz[0]*z[0]-zz[1]*z[1],zz[0]*z[1]+zz[1]*z[0])
                spectral+=lag.get(k,0)*zz[0]
            check(spectral>S,'actual source expands at rational unit-circle point')
            check(spectral>Fraction(1043,1000)*S,'exact quantitative expansion')
            exact_hostile=dict(point=['-35/37','12/37'],
                spectral_numerator=str(spectral.numerator),spectral_denominator=str(spectral.denominator),
                ratio_lower_bound='1043/1000')
        row=dict(N=n,sources=source_count,source_norm=S,source_ell=initial_mass,
                 clock_diagonal=D,frequency_norms=energies,clock_residue_norms=blocks,
                 signed_lag_sums=dict(sorted(lag.items())),max_clock=max_clock,
                 max_depth=maxdepth,prefix_peak=list(peak),
                 entry_port_diagonal_norm=sum(norm(A) for A in entry.values()),
                 wiener_budget=D+sum(abs(v) for v in lag.values()),
                 prefix_wiener_peak=prefix_wiener,exact_frequency_hostile=exact_hostile,
                 quotient_ell=quotient_mass,quotient_norm=quotient_energy,
                 shifted_derivative_norm=derivative_energy)
        rows.append(row)
        print(n,'S,D,0,pi,i',S,D,*energies.values(),'peak',peak,flush=True)
    args.out.write_text(json.dumps(dict(status='FINITE-EXACT',checks=CHECKS,
        source='literal B-chain weights',clock='canonical integer seed clock',
        cap=10000,unresolved_routes=0,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
