#!/usr/bin/env python3
from probe_periodic_frames import variables
from pathlib import Path
import json

def specs(n,q):
    anti=q in (3,4)
    orb,base=variables(n,q,anti)
    info={p:(k,s) for k,o in enumerate(orb) for p,s in o.items()}
    loworbit=range(q) if n//q==7 and q in (3,4) else [None]
    out=[]
    for low in loworbit:
        fixed=base.copy();ok=True
        if low is not None:
            for j in range(n):
                if j%q==low:continue
                p=tuple(sorted((low,j)));k,s=info[p]
                desired=(j==(low+1)%n)
                # sorted-bit means smaller->larger; desired is low->j.
                t=desired^(low>j)^s
                if k in fixed and fixed[k]!=t:ok=False;break
                fixed[k]=t
        if not ok:continue
        free=[i for i in range(len(orb)) if i not in fixed]
        bits=[]
        for k,o in enumerate(orb):
            rows=[[0]*n,[0]*n]
            for (u,v),s in o.items():
                for b in (0,1):
                    x,y=(u,v) if b^s else (v,u);rows[b][x]|=1<<y
            bits.append(rows)
        out.append(dict(n=n,q=q,anti=anti,low_out_orbit=low,orbits=len(orb),fixed={str(k):int(v) for k,v in fixed.items()},free=free,bits=bits))
    return out

def main():
    root=Path(__file__).parent
    allspec=[]
    for n,q in [(10,2),(12,4),(14,2),(15,3),(20,4),(21,3),(28,4)]:
        for row in specs(n,q):
            allspec.append(row)
            fn=root/f"frames_n{n}_q{q}_low{row['low_out_orbit']}.txt"
            with fn.open('w') as f:
                f.write(f"{n} {q} {row['orbits']} {len(row['free'])}\n")
                for k,rr in enumerate(row['bits']):
                    fixed=row['fixed'].get(str(k),-1)
                    f.write(str(fixed)+' '+' '.join(str(x) for pair in rr for x in pair)+'\n')
            print(fn.name,'free',len(row['free']),'presentations',1<<len(row['free']))
    (root/'residual_frame_specs.json').write_text(json.dumps(allspec,separators=(',',':'))+'\n')
if __name__=='__main__':main()
