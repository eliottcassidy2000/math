"""Exact Hamilton-kernel search for the inherited below-ray families."""
from pathlib import Path
import sys,json

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import inflation_bank as ib

TARGETS = {
    'seven_s1': [35,14],
    'three_s1': [15,6], 'three_s2': [29,12],
    'three_s3': [57,24], 'three_s4': [113,48],
    'twentyone_s1': [87,42], 'twentyone_s2': [173,84],
    'mixed_s1': [75,54,12],
    'mixed11_s1': [23,18,4], 'mixed11_s2': [45,36,8],
}


def trim(xs):
    out=list(xs)
    while out and not out[-1]:out.pop()
    return out


def kernel(rows,v):
    _,end=ib.visit_ends(rows,v)
    return trim([sum(end[-1][j]) for j in range(1,len(rows)+1)])


def main():
    bank,counts,hits=[[]],[],{}
    examined=0
    for n in range(1,8):
        codes={ib.cc.canonical(ib.cc.child(r,m)) for r in bank for m in range(1 << (n-1))}
        bank=[ib.cc.decode(n,c) for c in sorted(codes)]
        counts.append(len(bank))
        for code,rows in zip(sorted(codes),bank):
            h=ib.cc.hp(rows)
            names=[name for name,value in TARGETS.items() if value[0]==h]
            if not names:continue
            for v in range(n):
                coeff=kernel(rows,v)
                examined+=1
                for name in names:
                    if coeff==TARGETS[name] and name not in hits:
                        hits[name]=dict(n=n,code=code,rows=rows,port=v,coefficients=coeff)
    report=dict(ordinary_base_counts=counts,ports_examined=examined,targets=TARGETS,hits=hits,
                missing=[name for name in TARGETS if name not in hits])
    (HERE/'head_search.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
