// Exact Euclid midpoints of physical flips in two transitive ports.
// Reuse the certified return-cap test and path-profile DP for this same grammar.
#define main scalar_port_main
#include "../chanpak_square25_dilation_20261123/port_rays.cpp"
#undef main

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int n=std::stoi(argv[1]);
    if(n<4 || n>max_order) return 2;
    std::set<std::array<Count,8>> seen;
    std::string triangle;
    while(std::cin>>triangle) {
        auto lo=decode(n,triangle), guards=common_vertices(lo,n);
        for(int u=0;u<n;++u) for(int v=u+1;v<n;++v) {
            if (!(guards[u]&~(1<<v)) || !(guards[v]&~(1<<u))) continue;
            auto l=polynomial(profile(lo,n,u,v));
            for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) {
                if (i==u || i==v || j==u || j==v) continue;
                auto hi=lo;
                hi[i]^=1<<j;
                hi[j]^=1<<i;
                auto second=common_vertices(hi,n);
                if (!(second[u]&~(1<<v)) || !(second[v]&~(1<<u))) continue;
                auto h=polynomial(profile(hi,n,u,v));
                Count A=l.a+h.a,B=l.b+h.b,C=l.c+h.c;
                if(l.d+h.d || A<=0 || (A&(A-1))) continue;
                int sum=2-__builtin_ctzll(A);  // A 2^(x+y) = 4
                std::set<int> shifts;
                if(B) shifts.insert(1-__builtin_ctzll(std::abs(B)));
                if(C) shifts.insert(sum+__builtin_ctzll(std::abs(C))-1);
                if(B && C) {
                    int twice=sum+__builtin_ctzll(std::abs(C))-__builtin_ctzll(std::abs(B));
                    if(twice%2==0) shifts.insert(twice/2);
                }
                for(int x:shifts) {
                    int y=sum-x, clear=-std::min({0,x,y});
                    if(B*(Wide(1)<<(x+clear))+C*(Wide(1)<<(y+clear))!=-(Wide(2)<<clear)) continue;
                    if (!seen.insert({l.a,l.b,l.c,l.d,h.a,h.b,h.c,h.d}).second) continue;
                    std::cout<<triangle<<" "<<u<<" "<<v<<" "<<i<<" "<<j<<" "<<x<<" "<<y<<" "
                             <<l.a<<" "<<l.b<<" "<<l.c<<" "<<l.d<<" "
                             <<h.a<<" "<<h.b<<" "<<h.c<<" "<<h.d<<"\n";
                }
            }
        }
    }
    return 0;
}
