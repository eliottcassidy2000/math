#!/usr/bin/env python3
"""Recompute the quotient at each echo; exact residue and mass controls."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd
import importlib.util,json
p=Path(__file__).with_name('lrc_audit.py')
spec=importlib.util.spec_from_file_location('echo_checks',p);a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
def main():
    residues=mass=0;examples=[]
    for q in (3,4,5,6):
        s=next(i for i in range(2,q) if gcd(i,q)==1);m=s//2+1
        _,rows=a.params(q,m,2);V=list(range(1,q))+[q+1]
        for j,row in enumerate(rows):
            aa,bb,D,c=(row[k] for k in ('a','b','D','c'));V.append(q*aa)
            alpha=pow(s,-1,q);h=(alpha*s-1)//q;K=alpha*aa-bb+h;g=gcd(K,D)
            a.need(g==gcd(c,q*bb+1),'recomputed gcd')
            den,num,c1=D//g,K//g,c//g;rho=row['rho'];t=F(K,D)
            valid=[v for v in range(den) if c1<=v*num%den<=den-c1]
            a.need(len(valid)==den-2*c1+1,'exact residue count')
            for v in range(1,2*D+1):
                a.need((v%den in valid)==(a.pair.distance(v,t)>=rho),'actual marked-point test');residues+=1
            Z=[s+2*den,q*aa+3*den];W=V+Z
            P=max(s,Z[0]);N=max(q*aa,Z[1]);eps=F(1,100*den**2*max(W)**2)
            intervals=a.safe(W,rho-eps);mu=sum((R-L for L,R in intervals),F(0))
            a.need(mu==2*eps*(F(1,P)+F(1,N)),'literal terminal mass')
            a.need(a.safe(W,rho)==[(z,z) for z in sorted({t,1-t})],'same maxima retained');mass+=1
            if q==3 and j==1:
                a.need((D,K,c,g,den,len(valid))==(26,16,6,2,13,8),'nonunit hostile')
                examples.append(dict(q=q,j=j,D=D,K=K,c=c,g=g,admissible_mod_reduced=len(valid),admissible_mod_original=g*len(valid),mass=str(mu/eps)))
    out=dict(status='exact finite controls for the proved recursive attachment corollary',residue_checks=residues,terminal_mass_checks=mass,hostiles=examples)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
