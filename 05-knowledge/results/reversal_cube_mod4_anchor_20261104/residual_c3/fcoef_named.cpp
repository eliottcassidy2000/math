// Exact named two-port response coefficients.
// Input: K instances; each: n then n rows of '0'/'1' (row i col j = 1 iff i->j).
// For each instance, F(a,c) = #Hamilton paths of T + x + y with x->y, a_v=[v->x], c_v=[v->y],
// evaluated on every query with |a|<=2, |c|<=2; Boolean Mobius inversion gives the
// multilinear coefficient of a_S c_R for |S|,|R|<=2 (F has degree <=2 in a and in c).
// Output per instance: one line "n D" then D*D integers, row-major over (S,R),
// downset order: empty, singletons 0..n-1, pairs (i<j) lexicographic.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

static ll hamcount(int N, const vector<unsigned>& out) {
    vector<ll> dp((size_t)1 << N * 1, 0);
    dp.assign(((size_t)1 << N) * N, 0);
    for (int v = 0; v < N; v++) dp[((size_t)1 << v) * N + v] = 1;
    unsigned full = (1u << N) - 1;
    for (unsigned m = 1; m <= full; m++) {
        for (int v = 0; v < N; v++) {
            ll val = dp[(size_t)m * N + v];
            if (!val) continue;
            unsigned nx = out[v] & ~m;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(size_t)(m | (1u << w)) * N + w] += val; }
        }
    }
    ll s = 0; for (int v = 0; v < N; v++) s += dp[(size_t)full * N + v];
    return s;
}

int main() {
    int K; if (!(cin >> K)) return 0;
    while (K--) {
        int n; cin >> n; vector<string> rows(n); for (auto& r : rows) cin >> r;
        vector<vector<int>> D; D.push_back({});
        for (int i = 0; i < n; i++) D.push_back({i});
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) D.push_back({i, j});
        int d = D.size(); int N = n + 2, X = n, Y = n + 1;
        vector<unsigned> base(N, 0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (rows[i][j] == '1') base[i] |= 1u << j;
        vector<ll> F((size_t)d * d);
        for (int is = 0; is < d; is++) for (int ir = 0; ir < d; ir++) {
            vector<unsigned> out = base; out[X] = 1u << Y; out[Y] = 0;
            vector<int> a(n, 0), c(n, 0);
            for (int v : D[is]) a[v] = 1; for (int v : D[ir]) c[v] = 1;
            for (int v = 0; v < n; v++) {
                if (a[v]) out[v] |= 1u << X; else out[X] |= 1u << v;
                if (c[v]) out[v] |= 1u << Y; else out[Y] |= 1u << v;
            }
            F[(size_t)is * d + ir] = hamcount(N, out);
        }
        // Mobius: coefficient(S,R) = sum_{S'<=S,R'<=R} (-1)^{|S-S'|+|R-R'|} F(S',R')
        auto subs = [&](int idx) { vector<pair<int,int>> r; const auto& s = D[idx];
            if (s.empty()) r.push_back({0, 1});
            else if (s.size() == 1) { r.push_back({idx, 1}); r.push_back({0, -1}); }
            else { int i = s[0], j = s[1]; r.push_back({idx, 1}); r.push_back({1 + i, -1}); r.push_back({1 + j, -1}); r.push_back({0, 1}); }
            return r; };
        printf("%d %d\n", n, d);
        for (int is = 0; is < d; is++) { auto A = subs(is);
            for (int ir = 0; ir < d; ir++) { auto B = subs(ir); ll s = 0;
                for (auto& p : A) for (auto& q : B) s += (ll)p.second * q.second * F[(size_t)p.first * d + q.first];
                printf("%lld ", s); }
        }
        printf("\n");
    }
}
