"""Independent exact-grid and quaternion audit of the two main controls.

Does not import check.py. The two examples have every wall on the 1/24 grid.
"""
from fractions import Fraction as F
import json
from pathlib import Path

def check(c,m):
    if not c: raise RuntimeError(m)

def d(x):
    r=x%1
    return r if r<=F(1,2) else 1-r

# Entries are (sign, basis index) for 1,i,j,k in the standard quaternion table.
TABLE=[[(1,0),(1,1),(1,2),(1,3)],
       [(1,1),(-1,0),(1,3),(-1,2)],
       [(1,2),(-1,3),(-1,0),(1,1)],
       [(1,3),(1,2),(-1,1),(-1,0)]]

def left(a,x):
    s,i=a;t,j=x;u,k=TABLE[i][j]
    return s*t*u,k

def main():
    cert=json.loads(Path(__file__).with_name('certificate.json').read_text())
    report={}
    for name,phase in [('unshifted',F(0)),('shifted',F(1,4))]:
        v=(1,2);p=(F(0),phase);delta=F(1,3)
        good=[];atoms=[];events=[];word=[]
        for k in range(25):
            t=F(k,24);equal=[i for i in range(2) if d(p[i]+v[i]*t)==delta]
            if equal:
                events.append(t)
                # Exact atom routing: all exits occur before all entries.
                before=[d(p[i]+v[i]*(t-F(1,48)))>delta for i in range(2)]
                after=[d(p[i]+v[i]*(t+F(1,48)))>delta for i in range(2)]
                for i in equal:
                    if before[i]:word.append((i,-1))
                for i in equal:
                    if after[i]:word.append((i,1))
                if all(d(p[i]+v[i]*t)>=delta for i in range(2)):atoms.append(t)
        for k in range(24):
            t=F(2*k+1,48)
            if all(d(p[i]+v[i]*t)>delta for i in range(2)):good.append(k)
        check(sum(len([i for i in range(2) if d(p[i]+v[i]*t)==delta]) for t in events)==6,'all six wall incidences present')
        check(F(len(good),24)==F(cert['examples'][name]['strict_mass']),'independent grid mass')
        result=[]
        for s in (1,-1):
            for b in range(4):
                x=(s,b)
                for i,sign in word:x=left((sign,i+1),x)
                result.append((0 if x[0]==1 else 4)+x[1])
        check(result==cert['examples'][name]['holonomy'],'manual quaternion multiplication')
        report[name]={'strict_grid_cells':good,'closed_wall_atoms':list(map(str,atoms)),'quaternion_permutation':result}
    check(report['unshifted']['closed_wall_atoms']==['1/3','2/3'],'two isolated lonely atoms')
    check(report['shifted']['strict_grid_cells']==[13,14,15],'exact one-eighth window')
    regression=cert['examples']['strict_safe_cut']
    check(regression['strict_intervals']==[['0','1/6'],['5/6','1']],'independent cut-crossing endpoints')
    check(regression['cut_atom']['strict_safe'] is True,'cut atom is strictly lonely')
    intervals=[(F(l),F(r)) for l,r in regression['strict_intervals']]
    checked=[]
    for k in range(25):
        t=F(k,24);x=t%1
        reconstructed=any(l<x<r for l,r in intervals) or (x==0 and regression['cut_atom']['strict_safe'])
        direct=d(F(1,2)+t)>F(1,3)
        check(reconstructed==direct,'strict set including cut and equality atoms')
        if reconstructed:checked.append(str(t))
    check('0' in checked and '1' in checked,'both chart copies of cut retained')
    check('1/6' not in checked and '5/6' not in checked,'equality atoms not strict')
    report['strict_safe_cut']={'strict_grid_points':checked,'strict_mass':'1/3','cut_atom_verified':True}
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':main()
