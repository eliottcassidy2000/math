"""Local decay exponents of survival curves from qsurv outputs: alpha(T) = -dlog q/dlog T over a factor-10 window,
binomial standard errors, and sqrt(T) q(T)."""
import sys, math
def load(fn):
    rows=[]; head=''
    for line in open(fn):
        if line.startswith('#'):
            if 'paths=' in line: head=line.strip()
            continue
        T,s,q,_=line.split(); rows.append((int(T),int(s),float(q)))
    return head, rows
for fn in sys.argv[1:]:
    head,rows=load(fn)
    P=int(head.split('paths=')[1].split()[0])
    print(fn, '|', head[:140])
    d={T:(s,q) for T,s,q in rows}
    Ts=sorted(d)
    for T in Ts:
        if T<100: continue
        T2=None
        for U in Ts:
            if U>=10*T: T2=U;break
        s,q=d[T]
        line=f"  T={T:>8d} q={q:.5f} sqrtTq={math.sqrt(T)*q:7.3f}"
        if T2:
            s2,q2=d[T2]
            if s2>0:
                a=-math.log(q2/q)/math.log(T2/T)
                se=math.sqrt((1-q)/(s)+(1-q2)/(s2))/math.log(T2/T) if s>0 else float('nan')
                line+=f"  alpha[{T},{T2}]={a:.3f}+-{se:.3f}"
        if int(math.log10(T)*10)%5==0: print(line)
