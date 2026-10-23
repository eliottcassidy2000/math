# Quantitative runner-comb observability

This addendum belongs to `PROOF.md`. **PROOF-COMPLETE SCOPED COROLLARY; independently checked algebra and exact controls.** The stability metric is L1, not endpoint location or Hausdorff distance. No LRC arrival or boundary-point certificate is supplied.

Let b>=2, J>=1, q=b^J. Let K have c>=1 proper closed arc components with distinct endpoints, and let g be any measurable function taking values in [0,1]. Assume c<q. With the original translated-comb responses put

    h=1_K-g,
    e=||h||_1,
    delta=max_(|k|<=c) |hhat(k)|,
    rho=max_(j=0,...,J-1) ||O_(b^j)(h)||_infinity.

Then

    e <= min{1, (2c+1) delta^(1/(2c+1))},                 (N1)
    delta <= [pi*c/sin(pi/b)] rho,                        (N2)
    e <= min{1, (2c+1)[pi*c*rho/sin(pi/b)]^(1/(2c+1))}.   (N3)

The same conclusion holds if rho is defined using the L1 norm of each phase-response function. The Fourier estimate needed for (N2) is bounded by either norm. In particular there is the entirely rational constant version

    e^(2c+1) <= (2c+1)^(2c+1) * (11bc/7) * rho.          (N4)

(N4) follows from sin(pi/b)>=2/b and pi<22/7. All integer bases, including composite ones, are allowed.

## Proof of the endpoint-weighted estimate

The degree-c trigonometric polynomial

    P(t)=-product_i sin(pi(t-a_i))sin(pi(t-b_i))

has sign + on the interiors of K and - on its complement. Its Fourier coefficient l1 norm is at most 1, since each sine factor has two coefficients of modulus 1/2 and the 2c-fold product has only integer frequencies. Thus

    integral |h||P| = integral hP <= delta.

Remove radius-r neighborhoods of the 2c endpoints. Their measure is at most 4cr. Because |h|<=1, their discrepancy contribution is at most 4cr. Outside them, every sine factor has modulus at least 2r, so

    e <= 4cr + delta/(2r)^(2c).

For delta>0 choose r=delta^(1/(2c+1))/2. Since delta<=||h||_1<=1, this is an allowed r<=1/2. The two terms become 2c*delta^(1/(2c+1)) and delta^(1/(2c+1)), proving (N1). If delta=0 take r down to 0.

## Proof of the measured-noise conversion

For 1<=|k|<=c choose the unique j<J such that k=b^j m and b does not divide m. Then

    O_(b^j)(h)hat(k)=hhat(k) sin(pi*m/b)/(pi*m).

The numerator has modulus at least sin(pi/b), and |m|<=c. Every Fourier coefficient of an integrable function has modulus at most its L1 norm, hence at most its supremum norm on the normalized circle. This proves the nonzero-frequency bound in (N2). The zero-frequency bound is |hhat(0)|<=b*rho. It is included in (N2), because sin(pi/b)<=pi/b gives pi*c/sin(pi/b)>=bc>=b. Combining with (N1) proves (N3).

Empty and full K are easier: h has one sign and e=|hhat(0)|<=b*rho. They do not need an endpoint polynomial.

## Exact controls and limitations

The companion exhausts every pair (K,G) of 8-cell binary carriers with K having fewer than 4 components for b=2,J=2. It evaluates every response at every cell-grid phase exactly. All breakpoints occur at these phases, so the largest discrepancy is the exact continuum supremum. It checks (N4) by integer cross-multiplication and checks the empty/full estimate separately.

No optimality claim is made for the uniform-response exponent. The component budget is load-bearing. Isolated points have zero L1 mass and are deliberately invisible. Even bounded-component L1 stability does not control Hausdorff distance without additional minimum-width information: a very short distant component can have arbitrarily small mass. Nor does any geometric estimate provide the owner/word/root joint law needed for physical arrival.
