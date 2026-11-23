"""Summarize strong cores and converse orbits among ordinary fork centers."""
from collections import defaultdict
from pathlib import Path
import sys

from physical_navigation import direct
BANK = Path(__file__).resolve().parents[3]/"session212"/"cycle-levels"
sys.path.insert(0,str(BANK))
from build_inputs import canonical


def converse(rows):
    full=(1<<len(rows))-1
    return tuple(full ^ row ^ (1<<i) for i,row in enumerate(rows))


def reaches_all(rows):
    seen=front=1
    while front:
        bit=front&-front
        new=rows[bit.bit_length()-1]&~seen
        seen|=new
        front=(front^bit)|new
    return seen==(1<<len(rows))-1


def main(lines):
    counts=defaultdict(list)
    for line in lines:
        if not line.startswith("base "):
            continue
        middle,*rows=map(int,line.split()[1:])
        if middle>3:
            assert middle==direct(rows)
            counts[len(rows),middle//3].append(tuple(rows))
    for (n,h),group in sorted(counts.items()):
        strong=[r for r in group if reaches_all(r) and reaches_all(converse(r))]
        keys={canonical(r) for r in strong}
        assert len(keys)==len(strong) and {canonical(converse(r)) for r in strong}==keys
        fixed=sum(canonical(converse(r))==canonical(r) for r in strong)
        print("order",n,"source",h,"centers",len(group),"strong",len(strong),"strong self-converse",fixed)
    targets={5,11,13,29,71}
    for n in range(1,8):
        classes=defaultdict(list)
        for line in (BANK/f"classes{n}.rows").read_text().splitlines():
            r=tuple(map(int,line.split()))
            h=direct(r)
            if h in targets:
                classes[h].append(r)
        for h,group in sorted(classes.items()):
            fixed=sum(canonical(converse(r))==r for r in group)
            print("ordinary first height",h,"order",n,"classes",len(group),"self-converse",fixed)
            targets.remove(h)


if __name__=="__main__":
    main(sys.stdin)
