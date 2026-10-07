# start_state.py n0 -> prints the 288-digit trailing window of 2^n0 and its 54 leading digits (truncated)
import sys, mpmath
def start(n0):
    t = str(pow(2, n0, 10**288)).zfill(288)
    mpmath.mp.dps = 140
    x = n0 * mpmath.log10(2)
    f = x - mpmath.floor(x)
    M = mpmath.power(10, f)          # mantissa in [1,10)
    lead = int(mpmath.floor(M * mpmath.mpf(10)**53))
    s = str(lead)
    assert len(s) == 54, (n0, s)
    return t, s
if __name__ == '__main__':
    n0 = int(sys.argv[1]); t, s = start(n0); print(t); print(s)
