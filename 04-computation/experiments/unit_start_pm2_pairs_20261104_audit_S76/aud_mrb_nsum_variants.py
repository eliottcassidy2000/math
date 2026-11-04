"""aud_mrb_nsum_variants.py - which mpmath nsum set-ups give the S76 note's 'cautionary' value 0.19164?

The note (section 3) says mpmath's default extrapolation of the *paired* series gave 0.19164 (wrong), while its
alternating-series and Euler-Maclaurin methods give the correct value.  The exact set-up is not recorded, so this
tries the natural variants at the default precision (15 digits) and at 30 digits.
"""
import mpmath as mp

f = lambda n: mp.power(n, mp.mpf(1) / n) - 1
variants = {
    'C = sum_(k>=1) (a_2k - a_(2k-1))': (lambda k: f(2 * k) - f(2 * k - 1), 1, lambda s: s),
    'C = sqrt2-1 - sum_(k>=1) (a_(2k+1) - a_(2k+2))': (lambda k: f(2 * k + 1) - f(2 * k + 2), 1, lambda s: mp.sqrt(2) - 1 - s),
}
for dps in (15, 30):
    mp.mp.dps = dps
    for name, (term, k0, post) in variants.items():
        for meth in (None, 'r+s', 'richardson', 'shanks', 'e', 'a', 'l', 'd'):
            try:
                s = mp.nsum(term, [k0, mp.inf]) if meth is None else mp.nsum(term, [k0, mp.inf], method=meth)
                print("dps=%d  %-48s method=%-10s -> %s" % (dps, name, meth, mp.nstr(post(s), 12)))
            except Exception as ex:
                print("dps=%d  %-48s method=%-10s -> failed (%s)" % (dps, name, meth, ex))
    # the unpaired alternating series
    for meth in (None, 'a', 'e'):
        try:
            s = mp.nsum(lambda n: (-1) ** int(n) * f(n), [2, mp.inf]) if meth is None else mp.nsum(lambda n: (-1) ** int(n) * f(n), [2, mp.inf], method=meth)
            print("dps=%d  %-48s method=%-10s -> %s" % (dps, 'C = sum_(n>=2) (-1)^n a_n (unpaired)', meth, mp.nstr(s, 12)))
        except Exception as ex:
            print("dps=%d  unpaired method=%s failed (%s)" % (dps, meth, ex))
