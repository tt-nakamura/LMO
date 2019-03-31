// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZ.h>
#include<NTL/vec_GF2.h>
using namespace NTL;

long CountSmallPrimes(long x)
// input:
//   x = integer, 0 < x < 2^{30}
// return:
//   number of primes <= x
{
    PrimeSeq ps;
    long p,j(0);
    while((p=ps.next()) && p<=x) j++;
    if(p==0) TerminalError("x is too large");
    return j;
}

long RootInt(long n, long k)
// input:
//   n,k = positive integers
// return:
//   floor(n^{1/k})
{
    long x, y(long(round(pow(n,1./k))));
    do {
        x = y;
        y = ((k-1)*x + n/power_long(x,k-1))/k;
    } while(x>y);
    return x;
}

#define	mul21(x) (((x)<<1)+1)
#define	div21(x) (((x)+1)>>1)

long CountPrimes(long x)
// input:
//   x = positive integer
// return:
//   number of primes <= x
// reference:
//   J. C. Lagarias, V. S. Miller and A. M. Odlyzko
//     "Computing pi(x): The Meissel-Lehmer Method"
//     Mathematics of Computation 44 (1985) 537
{
    long x2,x3,x4,x6,x23,P2(0),S1(0),S2(0);
    long i,j,k,l,m,n,a,b,c,N,M,pi2,pi3;
    Vec<long> p,f,g,u,v;
    Vec<Vec<long> > A;
    vec_GF2 s,t;

    if(x<27) return CountSmallPrimes(x);

    x3 = RootInt(x,3);
    x6 = SqrRoot(x3);
    PrimeSeq ps; ps.reset(3);
    while((j = ps.next()) <= x3) p.append(j);

    M = div21(x3);
    s.SetLength(M);// sieve
    u.SetLength(M);
    f.SetLength(M);// least prime factor
    g.SetLength(M,1);// moebius function
    v.SetLength(p.length());
    A.SetLength(NumBits(M));
    for(i=0; i<A.length(); i++) A[i].SetLength(((M-1)>>i)+1);

    for(i=0; i<p.length(); i++) {
        v[i] = (i ? v[i-1] : M);
        for(j=p[i]>>1; j<M; j+=p[i]) {
            g[j] = -g[j];
            if(IsOne(s[j])) continue;
            set(s[j]);
            f[j] = p[i];
            v[i]--;
        }
    }
    for(i=0; p[i]<=x6; i++)
        for(k=p[i]*p[i], j=k>>1; j<M; j+=k) g[j] = 0;

    for(i=0; i<M; i++)
        if     (g[i]>0) S1 += div21(x/mul21(i));
        else if(g[i]<0) S1 -= div21(x/mul21(i));

    for(i=div21(x3/3); i<M; i++)
        if(f[i]<=3) continue;
        else if(g[i]>0) S2 -= div21(x/(mul21(i)*3));
        else if(g[i]<0) S2 += div21(x/(mul21(i)*3));

    x2 = SqrRoot(x);
    x4 = SqrRoot(x2);
    x23 = x/(x3+1)+1;
    pi3 = p.length()+1;// pi(x^{1/3})
    c = div21(x2);
    N = M<<1;
    u[M-1] = pi3;

    for(a=N+1; a<x23; a+=N) {
        if((b=a+N)>x23) b=x23;
        clear(s);
        for(i=0; i<A.length(); i++)
            for(j=0, k=1<<i; j<A[i].length(); j++) A[i][j] = k;
        for(i=0; i<p.length(); i++) {
            j = NegateMod(((a + p[i])>>1)%p[i], p[i]);
            for(; j<M; j+=p[i]) {
                if(IsOne(s[j])) continue;
                set(s[j]);
                for(k=0; k<A.length(); k++) A[k][j>>k]--;
            }
            if(i>=p.length()-2) continue;
            k = div21(x/(b*p[i+1]));
            l = div21(x/(a*p[i+1]));
            if(l>M) l=M;
            for(; k<l; k++) {
                if(f[k]<=p[i+1] || g[k]==0) continue;
                j = x/(mul21(k)*p[i+1]);
                j = ((j-a)>>1) + 1;
                for(m=0, n=v[i]; j; m++, j>>=1)
                    if(j&1) n += A[m][j-1];
                if(g[k]>0) S2 -= n;
                else       S2 += n;
            }
            for(k=0, j=M; j; k++, j>>=1)
                if(j&1) v[i] += A[k][j-1];
        }
        for(i=0, j=u[M-1]; i<M; i++)
            u[i] = (IsZero(s[i]) ? ++j : j);
        if(b<=x2) continue;
        if(a<=x2) pi2 = u[(x2-a)>>1];
        k = c;
        c = div21(x/b);
        k-=c;
        if(k<=0) continue;
        t.SetLength(k);
        clear(t);
        for(i=0; p[i]<=x4; i++) {
            j = NegateMod((c + div21(p[i]))%p[i], p[i]);
            for(; j<k; j+=p[i]) set(t[j]);
        }
        for(i=0; i<k; i++) {
            if(IsOne(t[i])) continue;
            j = x/mul21(c+i);
            P2 += u[(j-a)>>1];
        }
    }
    P2 += (pi3-1)*pi3>>1;
    P2 -= (pi2-1)*pi2>>1;
    return S1+S2-P2+pi3-1;
}
