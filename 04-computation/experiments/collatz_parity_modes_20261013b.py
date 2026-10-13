#!/usr/bin/env python3
"""Fine-bin parity resonance analysis; numerical model tests, not a proof.

Requires numpy, scipy, mpmath. `--edges FILE` writes the declared integer
bin endpoints. Otherwise provide --input census.json --output analysis.json.
"""
import argparse
import json
from decimal import Decimal, localcontext
import mpmath as mp
import numpy as np
from scipy.optimize import least_squares


def edges(bits):
    def make(precision):
        with localcontext() as c:
            c.prec = precision
            return [int(Decimal(2) ** (Decimal(14) + Decimal(j)/bits))
                    for j in range(20*bits+1)]
    a, b = make(90), make(120)
    if a != b or len(set(a)) != len(a):
        raise RuntimeError('endpoint control')
    return a


def resonance(f):
    mp.mp.dps = 45
    alpha = mp.log(3, 2)-1
    theta = mp.mpf(str(f))*mp.log(3, 2)
    theta -= mp.nint(theta)
    w = mp.mpc(0)
    for j in range(1, 121):
        t = theta*j/120
        w = mp.findroot(lambda v: mp.exp(-v)+mp.exp(alpha*v+2j*mp.pi*t)-2, w)
    residual = abs(mp.exp(-w)+mp.exp(alpha*w+2j*mp.pi*theta)-2)
    if residual > mp.mpf('1e-38'):
        raise RuntimeError('root residual')
    return {'base': f, 'theta': float(theta),
            'frequency': float(f+w.imag/(2*mp.pi)), 'decay': float(-w.real),
            'root_residual': float(residual)}


def basis(lo, hi, frequency, decay, anchor=24):
    lam = -decay+2j*np.pi*frequency
    z = 1+lam/np.log(2)
    v = np.log(hi/lo)
    # Average exp(lam*(log2 n-anchor)) over the INTEGER-WIDTH bin.
    return np.exp(lam*(np.log2(lo)-anchor))*np.expm1(z*v)/(z*np.expm1(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--edges')
    ap.add_argument('--bits', type=int, default=1024)
    ap.add_argument('--input')
    ap.add_argument('--output')
    args = ap.parse_args()
    if args.edges:
        with open(args.edges, 'w') as f:
            f.write('\n'.join(map(str, edges(args.bits)))+'\n')
        return
    D = json.load(open(args.input))
    r = np.array(D['bins'], dtype=float)
    lo, hi = r[:, 0], r[:, 1]
    x = (np.log2(lo)+np.log2(hi))/2
    y = r[:, 3]/r[:, 2]
    alpha = np.log2(3)-1
    bases = [j+.5 for j in range(512)
             if abs((j+.5)*(1+alpha)-round((j+.5)*(1+alpha))) <= .06]
    modes = [resonance(f) for f in bases]
    target = next(j for j,m in enumerate(modes) if m['base'] == 2.5)
    Z = np.array([basis(lo, hi, m['frequency'], m['decay']) for m in modes]).T
    A = np.column_stack([np.ones(len(x)), Z.real, Z.imag])
    reports = []
    for low, high in [(16, 26), (20, 30), (24, 34), (20, 28), (28, 34)]:
        mask = (x >= low)&(x < high)
        c = np.linalg.lstsq(A[mask], y[mask], rcond=None)[0]
        prediction = A@c
        rows = []
        for j,m in enumerate(modes):
            cc = c[j+1]-1j*c[j+1+len(modes)]
            rows.append({'base': m['base'], 'amplitude_at_log2n_24': abs(cc),
                         'phase_at_log2n_24': float(np.angle(cc))})
        # Project out the other predeclared resonance modes, then fit ONLY
        # frequency and decay of the base-5/2 mode (variable projection).
        keep = [j for j in range(A.shape[1])
                if j not in [target+1, target+1+len(modes)]]
        Q, _ = np.linalg.qr(A[mask][:, keep], mode='reduced')
        Y = y[mask]-Q@(Q.T@y[mask])
        def residual(p, ret=False):
            z = basis(lo[mask], hi[mask], p[0], p[1])
            B = np.c_[z.real,z.imag]
            B -= Q@(Q.T@B)
            v = np.linalg.lstsq(B,Y,rcond=None)[0]
            return (Y-B@v) if not ret else (float(np.mean((Y-B@v)**2)),v)
        fit = least_squares(residual, [2.45492347,.21215967],
                            bounds=([2.38,0],[2.55,.5]),
                            ftol=1e-10,xtol=1e-10,gtol=1e-10)
        candidates = []
        for name,fr in [('half_drift',1/(2-np.log2(3))),('doubling',2.5),
                        ('signed_resonance',modes[target]['frequency'])]:
            fitted = least_squares(lambda d:residual([fr,d[0]]), [.21],
                                    bounds=([0],[.5]))
            candidates.append({'model':name,'frequency':fr,'best_decay':float(fitted.x[0]),
                               'residual_rms':float(np.sqrt(np.mean(fitted.fun**2)))})
        reports.append({'log2_interval':[low,high],
                        'fitted_low_frequency':float(fit.x[0]),
                        'fitted_low_decay':float(fit.x[1]),
                        'residual_rms_all_fixed_modes':float(np.sqrt(np.mean((y[mask]-prediction[mask])**2))),
                        'candidate_low_modes':candidates,'mode_coefficients':rows,
                        'holdout_28_34_rms':float(np.sqrt(np.mean((y[x>=28]-prediction[x>=28])**2)))})
        print(low,high,'free low',fit.x,'candidates',candidates,flush=True)
    octaves=[]
    for k in range(14,34):
        mask=(x>=k)&(x<k+1)
        yy=y[mask]
        octaves.append({'exponent':k,'bins':int(sum(mask)),
                        'fine_rms':float(np.sqrt(np.mean(yy**2))),
                        'weighted_mean':float(sum(r[mask,3])/sum(r[mask,2]))})
    out={'scope':'NUMERICAL signed-transfer-model comparison; no universal mode expansion proved',
         'integer_census_N':D['N'],'bins_per_octave':len(y)//20,
         'half_drift':1/(2-np.log2(3)), 'modes':modes,
         'fits':reports,'octaves':octaves}
    with open(args.output,'w') as f:
        json.dump(out,f,indent=2);f.write('\n')


if __name__ == '__main__':
    main()
