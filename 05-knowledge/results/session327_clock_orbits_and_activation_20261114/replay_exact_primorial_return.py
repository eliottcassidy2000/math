from profile_tools import *
from json import loads
from pathlib import Path
def need(ok):
    if not ok:raise RuntimeError('primorial return mismatch')
def main():
    T={(0,)};found=set()
    for n in range(2,7):
        T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
        for A in T:
            if strong(A):
                h,q=profile(A)
                if h==23:found.add(q)
    need(found=={272,296})
    D=loads(Path(__file__).with_name('phase_table.json').read_text())
    need([D[str(c)] for c in (31,97,163)]==[0,0,3])
    B=[]
    for j in range(3):
        N=30017+30030*j;c=N*pow(23,-1,198)%198
        need(c==31+66*j)
        B.append((N-c*23)//99)
    need(B==[296,584,872] and 3*30030//99==910)
    tab=[[(b//2)*pow(q//2,-1,455)%455 for q in (272,296)] for b in B]
    need(tab==[[68,1],[417,334],[311,212]])
    orbit=[pow(2,i,455) for i in range(12)]
    need(len(set(orbit))==12 and pow(2,12,455)==1)
    need(set(orbit)&{x for row in tab for x in row}=={1})
    t=3*296*(2**12-1)//910
    need(t==3996 and 30017+t*30030==99*296*2**12+31*23==120029897)
    print('all fixed-H=23 clock cores; full primorial phase cycle; next exact return 120029897')
if __name__=='__main__':main()
