// Independent integer endpoints check the actual fast parity implementation.
int main(){try{
 std::array<unsigned char,1024> parity{};for(int i=1;i<1024;++i)parity[i]=parity[i>>1]^(i&1);U tested=0,state=0x123456789abcdefULL;
 auto check=[&](const Rows&a){auto d=make(a,false);Rows rev{};unsigned s=0,e=0;for(int i=0;i<10;++i){rev[i]=1023u^(1u<<i)^a[i];s|=(d.s[i]&1)<<i;e|=(d.e[i]&1)<<i;}
  need(parity_ends(a,parity)==e&&parity_ends(rev,parity)==s,"integer/bit-parity endpoint mismatch");++tested;
 };
 for(int n=1;n<=5;++n)for(unsigned code=0;code<(1u<<(n*(n-1)/2));++code){Rows a{};int bit=0;for(int i=0;i<10;++i)for(int j=i+1;j<10;++j){bool forward=j>=n||((code>>bit)&1);if(j<n)++bit;a[forward?i:j]|=1u<<(forward?j:i);}check(a);}
 for(int k=0;k<1024;++k){Rows a{};for(int i=0;i<10;++i)for(int j=i+1;j<10;++j){state^=state<<13;state^=state>>7;state^=state<<17;bool forward=state&1;a[forward?i:j]|=1u<<(forward?j:i);}check(a);}
 std::cout<<"{\"embedded_all_labelled_orders1_to5\":1099,\"specified_random_order10_controls\":1024,\"integer_endpoint_comparisons\":"<<tested<<",\"status\":\"PASS\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
