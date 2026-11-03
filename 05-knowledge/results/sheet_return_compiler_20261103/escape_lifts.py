"""Literal source-preserving probes; output is exploratory until audited."""
import sys, json
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'defect_offset_dynamics_20261103'))
from wing import source4, hwing, tables, strong
from ear_cuts import ear_values, birth
from carrier import native


def paths(m,n):
    T,ad=source4(m,1,1,n)
    straight=tuple(range(len(T)))
    out=[('straight',straight)]
    if m>=3:
        revisit=(1,)+tuple(range(m,len(T)))+(0,)+tuple(range(2,m))
        out.append(('revisit',revisit))
    return T,ad,out


def main():
    data=[]
    for m in range(1,5):
      for n in range(1,5):
        T,ad,pp=paths(m,n);h=hwing(m,1,1,n)
        for label,p in pp:
            U=native.native(T,p);vals=ear_values(U)
            hits={str(c):[mask for mask,x in enumerate(vals) if x==3*h+c] for c in (2,4)}
            mask=sum(1<<ad.index(x) for x in ({(0,0),(0,m-1),(2,0),(5,n-2)} if label=='straight' and n>=2 else {(0,1),(1,0),(2,0)}) if x in ad)
            row=dict(m=m,n=n,path_type=label,h=h,intermediate=vals[0],prototype_mask=mask,prototype_h=vals[mask],hits=hits)
            data.append(row);print(row,flush=True);tables.cache_clear()
    (HERE/'exploration.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
