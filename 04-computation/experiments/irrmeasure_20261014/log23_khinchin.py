"""Continued fraction of theta = log2 3 to N partial quotients (mpmath at high precision), Khinchin statistics
(geometric mean vs Khinchin's constant 2.6854520010), Levy constant (log q_n / n vs pi^2/(12 ln 2) = 1.1865691104),
largest partial quotients, and the empirical exponent proxy max_n log a_(n+1)/log q_n (HYP-9262 evidence)."""
import sys, time
from mpmath import mp, mpf, log, floor, exp, pi
N = int(sys.argv[1]) if len(sys.argv)>1 else 20000
mp.dps = int(N*1.2)+200   # each partial quotient consumes about log10(a)+ small digits; generous
t=time.time()
theta = log(3)/log(2)
x = theta; a=[]
for i in range(N):
    ai = int(floor(x)); a.append(ai); x = x - ai
    if x == 0: break
    x = 1/x
a=a[:N]
# convergent denominators via integers
q0,q1 = 0,1; logq=[]
from math import log as mlog, exp as mexp
geo = 0.0; mx=0; mxi=-1; best_ratio=0.0; best_i=-1
qprev,qcur = 0,1
for i,ai in enumerate(a):
    if i>=1:
        geo += mlog(ai)
        if ai>mx: mx, mxi = ai, i
    if i>=1:
        # q_{i} = a_i q_{i-1} + q_{i-2}
        qprev,qcur = qcur, ai*qcur+qprev
        lq = mlog(qcur)
        if i+1 < len(a) and lq>0:
            r = mlog(a[i+1])/lq
            if r>best_ratio: best_ratio, best_i = r, i
print("theta = log2 3; %d partial quotients computed in %.1fs (mp.dps=%d)"%(len(a),time.time()-t,mp.dps))
print("first 60:", a[:60])
print("geometric mean of a_1..a_%d = %.6f  (Khinchin's constant 2.685452)"%(len(a)-1, mexp(geo/(len(a)-1))))
print("log q_n / n at n=%d: %.6f  (Levy constant pi^2/(12 ln 2) = 1.186569)"%(len(a)-1, mlog(qcur)/(len(a)-1)))
print("largest partial quotient: a_%d = %d"%(mxi,mx))
print("max over n of log a_(n+1) / log q_n: %.5f at n=%d (a_(n+1)=%d); mu(theta) = 2 + limsup of this ratio"%(best_ratio,best_i,a[best_i+1]))
top = sorted(range(1,len(a)), key=lambda i:-a[i])[:10]
print("ten largest partial quotients (index: value):", [(i,a[i]) for i in top])
