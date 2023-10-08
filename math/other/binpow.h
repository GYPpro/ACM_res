#ifndef _IN_TEMPLATE_
#include <template_overAll.h>
#endif

ll int binpow(ll int a,ll int b) 
{
    ll int res = 1;
    while (b > 0) 
    {
        if (b & 1) res = res * a % mod1;
        a = a * a % mod1;
        b >>= 1;
    }
    return res % mod1;    
}