"""Exact rational check of the spacer certificate. Python 3."""
from fractions import Fraction
from pathlib import Path

def require(predicate):
    if not predicate:
        raise RuntimeError("certificate failed")

def upper_inverse_power(p, a, b, denominator=10**6):
    # Least integer c with (c/denominator)^b >= p^(-a/b * b).
    low, high = 0, denominator
    target = denominator**b
    multiplier = p**a
    while high-low>1:
        c=(low+high)//2
        if c**b*multiplier >= target:
            high=c
        else:
            low=c
    require(low**b*multiplier<target<=high**b*multiplier)
    return Fraction(high,denominator)

def main():
    txt=(Path(__file__).parent/'../chanpak_guard_depth_20261115/replay_menu.txt')
    if not txt.exists():  # same checkout before publication
        txt=Path(__file__).parent/'../chanpak_guard_depth/replay_menu.txt'
    ps=[2]
    for line in txt.read_text().splitlines():
        _,cnt,label,*rows=line.split()
        require(label=='profiles:' and len(rows)==int(cnt))
        ps.extend(int(t.split(':')[1]) for t in rows)
    require(len(ps)==24)
    upper=2*upper_inverse_power(1024,21,40)
    for p in ps:
        upper/=1-upper_inverse_power(p,21,40)
    require(upper < Fraction(9971,10000))
    print(f"{len(ps)} spacer profiles; exact upper bound below 0.9971")
    print(f"approximate rational bound: {float(upper):.12f}")

if __name__=='__main__':
    main()
