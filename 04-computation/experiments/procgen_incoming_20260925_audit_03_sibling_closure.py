#!/usr/bin/env python3
"""Independent audit of creative_sibling / creative_decoder:
 - inductive closure counts from seed 1 over odd n <= N for the rules
   direct (U(n)<n), canonical (Q(n)<n), direct-or-canonical, full sibling ladder;
 - 483 certified via 181; first uncertified odd source;
 - forward-closure theorem: adding inverse ports (2n-1)/3 (n=5 mod 6) and
   (8n-5)/9 (n=13 mod 18) and ALL finite inverse-orbit joins changes nothing;
 - local coverage fractions 3/4 (mod 48 residual) and 19/24 (mod 144 residual).
Own code; plus sheet, accelerated odd map U(n)=(3n+1)/2^v.
"""
import sys


def U(n):
    m = 3 * n + 1
    while m % 2 == 0:
        m //= 2
    return m


def ladder(m):
    """all sibling ancestors of m (m, (m-1)/4, ...) while m = 5 mod 8"""
    out = [m]
    while m % 8 == 5:
        m = (m - 1) // 4
        out.append(m)
    return out


def closure(N, rule, extra_inverse=False):
    cert = bytearray(N + 1)
    cert[1] = 1
    for n in range(3, N + 1, 2):
        u = U(n)
        ok = False
        if rule in ("direct", "both", "full"):
            if u < n and cert[u]:
                ok = True
        if not ok and rule in ("canon", "both"):
            qn = ladder(u)[-1]
            if qn < n and cert[qn]:
                ok = True
        if not ok and rule == "full":
            for a in ladder(u)[1:]:
                if a < n and cert[a]:
                    ok = True
                    break
        if not ok and extra_inverse:
            if n % 6 == 5:
                p = (2 * n - 1) // 3
                if cert[p]:
                    ok = True
            if not ok and n % 18 == 13:
                y = (8 * n - 5) // 9
                if cert[y]:
                    ok = True
        if ok:
            cert[n] = 1
    return cert


def forward_closed(cert, N):
    bad = []
    for n in range(1, N + 1, 2):
        if cert[n]:
            u = U(n)
            if u <= N and not cert[u]:
                bad.append(n)
    return bad


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    res = {}
    for rule in ("direct", "canon", "both", "full"):
        c = closure(N, rule)
        res[rule] = c
        cnt = sum(c[1::2])
        first_missing = next(n for n in range(1, N + 1, 2) if not c[n])
        print("rule %-6s certified odd sources <= %d (incl. 1): %d ; first uncertified %d"
              % (rule, N, cnt, first_missing))
    d, cn = res["direct"], res["canon"]
    print("direct subset of canon?", all(cn[n] for n in range(1, N + 1, 2) if d[n]),
          "; canon subset of direct?", all(d[n] for n in range(1, N + 1, 2) if cn[n]))
    full = res["full"]
    new_full = [n for n in range(1, N + 1, 2) if full[n] and not res["both"][n]]
    print("full minus (direct or canon): %d, first %s" % (len(new_full), new_full[:5]))
    print("483: U=%d ladder=%s cert(181)=%d cert_full(483)=%d cert_both(483)=%d"
          % (U(483), ladder(U(483)), full[181], full[483], res["both"][483]))
    print("241: U=%d ladder=%s ; canon Q(241)=%d certified in canon-only? %d; in direct? %d"
          % (U(241), ladder(U(241)), ladder(U(241))[-1], cn[241], d[241]))
    # forward closure within range (U(n) may exceed N; only check in-range)
    for rule in ("direct", "both", "full"):
        bad = forward_closed(res[rule], N)
        print("forward-closed (in range) for %s: %s" % (rule, "YES" if not bad else ("NO " + str(bad[:5]))))
    c2 = closure(N, "full", extra_inverse=True)
    print("full + inverse ports (2n-1)/3 and (8n-5)/9: %d (full alone %d)" % (sum(c2[1::2]), sum(full[1::2])))
    c3 = closure(N, "both", extra_inverse=True)
    print("both + inverse ports: %d (both alone %d)" % (sum(c3[1::2]), sum(res['both'][1::2])))

    # local coverage: Q(n)<n iff n=1 mod 4 or n=3 mod 16 (n>1)
    ok = all(((ladder(U(n))[-1] < n) == (n % 4 == 1 or n % 16 == 3)) for n in range(3, 200001, 2))
    print("Q(n)<n iff n=1 mod4 or 3 mod16 (odd 3..200001):", ok)
    # residual mod 48 after rules (7) and (8)
    def covered_mod(M, with11):
        cov = []
        for r in range(1, M, 2):
            # test local applicability on the class: rules depend on n mod 16 / mod 6 / mod 18
            c = (r % 4 == 1) or (r % 16 == 3) or (r % 6 == 5) or (with11 and r % 18 == 13)
            cov.append((r, c))
        return cov
    cov48 = covered_mod(48, False)
    res48 = [r for r, c in cov48 if not c]
    print("mod48 residual:", res48, " covered %d/%d" % (sum(c for _, c in cov48), len(cov48)))
    cov144 = covered_mod(144, True)
    res144 = [r for r, c in cov144 if not c]
    print("mod144 residual:", res144, " covered %d/%d" % (sum(c for _, c in cov144), len(cov144)))


if __name__ == "__main__":
    main()
