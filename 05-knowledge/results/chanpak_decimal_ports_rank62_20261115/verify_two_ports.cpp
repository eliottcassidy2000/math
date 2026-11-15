#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using namespace std;
uint64_t literal_count(const string& marks){
    const unsigned n=marks.size();
    const uint64_t all=(1ull<<n)-1;
    uint64_t count=0;
    for(uint64_t before=0; before<=all; ++before){
        const uint64_t after=all^before;
        // Every path restricted to either side of the apex is forced to
        // follow the transitive spine order. Test the two incident arcs.
        if(before && marks[63-__builtin_clzll(before)]!='0')continue;
        if(after && marks[__builtin_ctzll(after)]!='1')continue;
        ++count;
    }
    return count;
}
int main(){
    assert(literal_count("11010100100110")==8259);
    assert(literal_count("10100110011100")==8259);
    assert(literal_count("110101001000010101010110")==8126555);
    assert(literal_count("101000000000010110011100")==8126553);
    cout<<"four literal deletion-path counts passed\n";
}

