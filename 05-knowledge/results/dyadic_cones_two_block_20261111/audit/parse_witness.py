import re,sys
def parse(fn):
    out=[]
    for l in open(fn):
        if not l.startswith('r='): continue
        m=re.match(r'r=(-?\d+) j=(\d+) m=(\d+) O=(\d+)(?: c=(\d+))? B:(.*)',l.strip())
        r,j,mm,O,c,B=m.groups(); out.append((int(c) if c else 1,int(r),int(j),int(mm),int(O),list(map(int,B.split()))))
    return out
