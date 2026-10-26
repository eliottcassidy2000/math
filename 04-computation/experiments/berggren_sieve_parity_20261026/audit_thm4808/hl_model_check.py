# HYP-9418 (3): the HL-Mertens model value 4 e^(-2 gamma) (X/log^2 X)/Li2(X), Li2(X) = int_2^X dt/log^2 t,
# and the HL prediction 2 C2 Li2(X) for pi2(X) (independent recomputation).
from mpmath import mp, quad, log, exp, euler, mpf, nprod
mp.dps = 30
C2 = mpf('0.660161815846869573927812110014555778432623360284733413319448')
for e, pi2 in ((8, 440312), (9, 3424506), (10, 27412679)):
    X = mpf(10) ** e
    Li2 = quad(lambda t: 1 / log(t) ** 2, [2, 10, 1e3, 1e5, 1e7, X])
    model = 4 * exp(-2 * euler) * (X / log(X) ** 2) / Li2
    print(f"10^{e}: Li2 = {mp.nstr(Li2, 12)}, 2C2 Li2 = {mp.nstr(2*C2*Li2, 10)} (pi2 = {pi2}), model rho = {mp.nstr(model, 6)}, 4e^(-2gamma) = {mp.nstr(4*exp(-2*euler), 6)}")
