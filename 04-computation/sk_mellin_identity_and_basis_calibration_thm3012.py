"""(A) A clean integral identity:  S(k) = int_0^1 2F1(1/4,3/4;1;u^k) du.
    Derivation: 1/(kn+1) = (1/k)/(n+1/k) and 1/(n+c) = int_0^1 t^{n+c-1}dt give
    S(k) = (1/k) int_0^1 t^{1/k-1} 2F1(1/4,3/4;1;t) dt; substitute t = u^k.
    This identifies S(k) with kps's family S_a(k) = int_0^1 2F1(a,1-a;1;x^k)dx at a=1/4.

(B) PSLQ for pi*S(4), pi*S(5) with ALGEBRAIC MULTIPLIERS -- the basis shape that
    actually contains pi*S(3) = sqrt3*log(5+2sqrt6) - 2 arctan(sqrt2/5).
    LIVE POSITIVE CONTROL: the same pipeline must rediscover S(3) first.
"""
import mpmath as mp
mp.mp.dps = 240

def S(k):
    return mp.hyper([mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/k], [1, 1+mp.mpf(1)/k], 1)

def S_int(k):
    f = lambda u: mp.hyp2f1(mp.mpf(1)/4, mp.mpf(3)/4, 1, u**k)
    return mp.quad(f, [0, 1])

print("(A) identity S(k) = int_0^1 2F1(1/4,3/4;1;u^k) du")
mp.mp.dps = 40
for k in (1, 2, 3, 4, 5, 6):
    a, b = S(k), S_int(k)
    print(f"   k={k}: 3F2={mp.nstr(a,28)}  integral={mp.nstr(b,28)}  diff={mp.nstr(abs(a-b),4)}")
mp.mp.dps = 240

# ---- basis with algebraic multipliers
sq = {n: mp.sqrt(n) for n in (2,3,5,6,7,10,15,30)}
logs = {
  'log2': mp.log(2), 'log3': mp.log(3), 'log5': mp.log(5),
  'log(1+r2)': mp.log(1+sq[2]), 'log(2+r3)': mp.log(2+sq[3]),
  'log(5+2r6)': mp.log(5+2*sq[6]), 'log(phi)': mp.log((1+sq[5])/2),
  'log(3+2r2)': mp.log(3+2*sq[2]), 'log(9+4r5)': mp.log(9+4*sq[5]),
  'atan(r2)': mp.atan(sq[2]), 'atan(r2/5)': mp.atan(sq[2]/5),
  'atan(1/3)': mp.atan(mp.mpf(1)/3), 'atan(r5)': mp.atan(sq[5]),
  'pi': mp.pi, 'one': mp.mpf(1),
}
mults = {'1': mp.mpf(1), 'r2': sq[2], 'r3': sq[3], 'r5': sq[5], 'r6': sq[6], 'r10': sq[10], 'r15': sq[15]}

def sweep(target, tname, mult_keys, log_keys, H, tol_digits):
    names, vals = [], []
    for mk in mult_keys:
        for lk in log_keys:
            names.append(f"{mk}*{lk}"); vals.append(mults[mk]*logs[lk])
    r = mp.pslq([target]+vals, maxcoeff=H, maxsteps=2*10**6, tol=mp.mpf(10)**(-tol_digits))
    if r and r[0] != 0:
        expr = " + ".join(f"{c}*{n}" for c, n in zip(r[1:], names) if c)
        return f"FOUND: {r[0]}*{tname} + {expr} = 0"
    return None

print("\n(B) POSITIVE CONTROL -- can the pipeline find pi*S(3)?")
ctl = sweep(S(3)*mp.pi, 'pi*S(3)', ['1','r3'], ['log(5+2r6)','atan(r2/5)'], 10**6, 200)
print(f"   minimal basis {{1,r3}}x{{log(5+2r6),atan(r2/5)}}: {ctl}")
ctl2 = sweep(S(3)*mp.pi, 'pi*S(3)', ['1','r2','r3','r5','r6'],
             ['log2','log3','log(1+r2)','log(2+r3)','log(5+2r6)','atan(r2)','atan(r2/5)','pi'],
             10**5, 200)
print(f"   WIDE basis (5 mults x 8 logs = 40 elts): {'FOUND' if ctl2 else 'not found'}")
if ctl2: print(f"      {ctl2}")
