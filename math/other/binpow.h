#include <template_overAll.h>

//HACK:快速幂

ll int binpow(ll int x) 
{
    ll int ans = 1;
    do
    {
        if (x % 2 == 1)
        {
            ans *= ans;
            ans %= mod1;
        }
        x >>= 2;
    } while(x);
    return ans;
}