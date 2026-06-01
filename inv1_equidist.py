"""Investigation 1: Equidistribution of Galois conjugates.
f(q) = |O_q cap B_q^m|/(q-1) vs volume V_q = (width/q)^m and limit (1-2/n)^m.
"""
from lrc import *

def run(speeds, n, primes):
    m = len(speeds)
    rows = []
    for q in primes:
        lo, hi = safe_box(q, n)
        if hi < lo:
            continue
        f = fraction_lonely(speeds, q, n)
        Vq = box_volume(q, n, m)
        rows.append((q, lonely_count(speeds, q, n), q-1, float(f), float(Vq)))
    return rows

if __name__ == "__main__":
    primes = primes_up_to(200)
    print("="*78)
    print("EQUIDISTRIBUTION: f(q)=|O_q cap B_q^m|/(q-1) vs box volume V_q")
    print("="*78)

    cases = [
        ("m=3, n=4", 4, [
            [1,2,4], [2,3,5], [1,4,6], [3,5,7], [2,5,11],
        ]),
        ("m=4, n=5", 5, [
            [1,2,4,8], [1,3,5,7], [2,3,5,7], [1,4,9,11],
        ]),
        ("m=5, n=6", 6, [
            [1,2,4,8,16], [1,3,5,7,9], [2,3,5,7,11],
        ]),
    ]
    for label, n, sets in cases:
        m = n-1
        print(f"\n### {label}   limit V_inf=(1-2/n)^m = {float(limiting_volume(n,m)):.5f}")
        for sp in sets:
            rows = run(sp, n, primes)
            # report a spread of q values
            sel = [r for r in rows if r[0] in (11,23,47,97,151,199)]
            print(f"  speeds={sp}")
            print(f"    {'q':>4} {'cnt':>5} {'q-1':>4} {'f(q)':>8} {'V_q':>8}  f-V")
            for (q,c,qm1,f,V) in sel:
                print(f"    {q:>4} {c:>5} {qm1:>4} {f:>8.4f} {V:>8.4f}  {f-V:+.4f}")
            # average |f-V| over all primes (large enough)
            big = [r for r in rows if r[0] >= 50]
            if big:
                avg = sum(abs(r[3]-r[4]) for r in big)/len(big)
                fbar = sum(r[3] for r in big)/len(big)
                Vbar = sum(r[4] for r in big)/len(big)
                print(f"    [q>=50: mean f={fbar:.4f}, mean V={Vbar:.4f}, mean|f-V|={avg:.4f}]")
