// uses NTL
//   www.shoup.net/ntl

#include<NTL/tools.h>

long CountPrimes(long);

main() {
    long n1(32),n2(54),n,pi;
    double t;
    for(n=n1; n<=n2; n++) {
        t = NTL::GetTime();
        pi = CountPrimes(1L<<n);
        t = NTL::GetTime() - t;
        std::cout << n << '\t';
        std::cout << pi << '\t';
        std::cout << t << '\n';
    }
}
