// Exact named two-port response coefficients over a chosen variable set, other bits fixed.
// Input: K instances; each: n, n rows of '0'/'1' (row i col j = 1 iff i->j),
//   base a-mask string (n chars), base c-mask string (n chars), nv, then nv variable names.
// F(a,c) = #Hamilton paths of T + x + y with x->y, a_v=[v->x], c_v=[v->y]. F has degree <= 2 in a
// and <= 2 in c, so its multilinear coefficients over the variable set (other bits held at the base
// masks) are recovered by Boolean Mobius inversion from all queries flipping <= 2 variables per port.
// Output per instance: line "nv D" then D*D integers row-major over (S,R); downset order:
// empty, singletons (variable list order), pairs (i<j in list order).
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
static ll hamcount(int N, const vector<unsigned>& out) {
    vector<ll> dp(((size_t)1 << N) * N, 0);
    for (int v = 0; v < N; v++) dp[((size_t)1 << v) * N + v] = 1;
    unsigned full = (1u << N) - 1;
    for (unsigned m = 1; m <= full; m++)
        for (int v = 0; v < N; v++) {
            ll val = dp[(size_t)m * N + v]; if (!val) continue;
            unsigned nx = out[v] & ~m;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(size_t)(m | (1u << w)) * N + w] += val; }
        }
    ll s = 0; for (int v = 0; v < N; v++) s += dp[(size_t)full * N + v];
    return s;
}
int main() {
    int K; if (!(cin >> K)) return 0;
    while (K--) {
        int n; cin >> n; vector<string> rows(n); for (auto& r : rows) cin >> r;
        string a0, c0; cin >> a0 >> c0; int nv; cin >> nv; vector<int> var(nv); for (auto& x : var) cin >> x;
        vector<vector<int>> D; D.push_back(vector<int>());
        for (int i = 0; i < nv; i++) D.push_back({i});
        for (int i = 0; i < nv; i++) for (int j = i + 1; j < nv; j++) D.push_back({i, j});
        int d = D.size(), N = n + 2, X = n, Y = n + 1;
        vector<unsigned> base(N, 0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (rows[i][j] == '1') base[i] |= 1u << j;
        vector<ll> F((size_t)d * d);
        for (int is = 0; is < d; is++) for (int ir = 0; ir < d; ir++) {
            vector<unsigned> out = base; out[X] = 1u << Y; out[Y] = 0;
            vector<int> a(n), c(n);
            for (int v = 0; v < n; v++) { a[v] = a0[v] == '1'; c[v] = c0[v] == '1'; }
            for (int t : D[is]) a[var[t]] ^= 1;
            for (int t : D[ir]) c[var[t]] ^= 1;
            for (int v = 0; v < n; v++) {
                if (a[v]) out[v] |= 1u << X; else out[X] |= 1u << v;
                if (c[v]) out[v] |= 1u << Y; else out[Y] |= 1u << v;
            }
            F[(size_t)is * d + ir] = hamcount(N, out);
        }
        auto subs = [&](int idx) { vector<pair<int,int>> r; const auto& s = D[idx];
            if (s.empty()) r.push_back({0, 1});
            else if (s.size() == 1) { r.push_back({idx, 1}); r.push_back({0, -1}); }
            else { r.push_back({idx, 1}); r.push_back({1 + s[0], -1}); r.push_back({1 + s[1], -1}); r.push_back({0, 1}); }
            return r; };
        printf("%d %d\n", nv, d);
        for (int is = 0; is < d; is++) { auto A = subs(is);
            for (int ir = 0; ir < d; ir++) { auto B = subs(ir); ll s = 0;
                for (auto& p : A) for (auto& q : B) s += (ll)p.second * q.second * F[(size_t)p.first * d + q.first];
                printf("%lld ", s); } }
        printf("\n");
    }
}
