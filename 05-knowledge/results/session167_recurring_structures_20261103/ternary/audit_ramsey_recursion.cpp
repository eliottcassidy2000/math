// Independent audit: recursive edge colors with immediate red/blue triangle pruning.
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
int main() {
    const std::array<std::array<uint64_t, 2>, 4> expected = {{{25,25},{543,543},{29196,29197},{3730190,3731612}}};
    for (int n=3; n<=6; ++n) {
        int color[6][6];
        for (auto &row : color) for (int &entry : row) entry=-1;
        std::vector<std::pair<int,int>> edges;
        for (int i=0; i<n; ++i) for (int j=i+1; j<n; ++j) edges.emplace_back(i,j);
        uint64_t leaves=0, good5=0, good6=0;
        std::function<void(int)> visit = [&](int depth) {
            if (depth==int(edges.size())) {
                ++leaves;
                bool green5=false, green6=false;
                for (int subset=0; subset<(1<<n); ++subset) {
                    int size=__builtin_popcount(unsigned(subset));
                    if (size!=5 && size!=6) continue;
                    bool allgreen=true;
                    for (int i=0; i<n; ++i) if (subset>>i&1)
                        for (int j=i+1; j<n; ++j) if ((subset>>j&1) && color[i][j]!=2) allgreen=false;
                    if (size==5 && allgreen) green5=true;
                    if (size==6 && allgreen) green6=true;
                }
                good5+=!green5;
                good6+=!green6;
                return;
            }
            auto [i,j]=edges[depth];
            for (int c=0; c<3; ++c) {
                bool allowed=true;
                if (c<2) for (int h=0; h<n; ++h)
                    if (h!=i && h!=j && color[i][h]==c && color[j][h]==c) allowed=false;
                if (!allowed) continue;
                color[i][j]=color[j][i]=c;
                visit(depth+1);
                color[i][j]=color[j][i]=-1;
            }
        };
        visit(0);
        if (good5!=expected[n-3][0] || good6!=expected[n-3][1]) throw std::runtime_error("Ramsey recursion mismatch");
        std::cout << "{\"n\":" << n << ",\"red_blue_triangle_free\":" << leaves << ",\"valid_335\":" << good5 << ",\"valid_336\":" << good6 << "}\n";
    }
}
