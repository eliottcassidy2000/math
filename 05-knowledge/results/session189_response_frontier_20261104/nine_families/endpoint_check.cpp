#include <array>
#include <algorithm>
#include <iostream>
#include <cstdint>
int main(){int cases;std::cin>>cases;while(cases--){std::array<unsigned,9>Q;for(auto&r:Q)std::cin>>r;std::array<int,9>p{0,1,2,3,4,5,6,7,8};std::array<uint64_t,81>J{};do{bool good=true;for(int i=0;i<8;i++)if(!(Q[p[i]]>>p[i+1]&1)){good=false;break;}if(good)J[9*p[0]+p[8]]++;}while(std::next_permutation(p.begin(),p.end()));for(auto n:J)std::cout<<n<<' ';std::cout<<'\n';}}
