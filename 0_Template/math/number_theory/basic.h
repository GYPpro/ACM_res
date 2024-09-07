__builtin_ffsll(x)
// 返回 x 的二进制末尾最后一个 1 的位置

__builtin_clzll(x)
// 返回 x 的二进制的前导 0 的个数。

__builtin_ctzll(x)
// 返回 x 的二进制末尾连续 0 的个数。

__builtin_clrsbll(x)
// 当 x 的符号位为 0 时返回 x 的二进制的前导 0 的个数减一，否则返回 x 的二进制的前导 1 的个数减一。

__builtin_popcountll(x)
// 返回 x 的二进制中 1 的个数。

__builtin_parity(x)
// 判断 x 的二进制中 1 的个数的奇偶性。

int binpow(int x, int y)
{
    int res = 1;
    while (y > 0)
    {
        if (y & 1)
            res = res * x % mod;
        x = x * x % mod;
        y >>= 1;
    }
    return res;
}

void exgcd(int a, int b, int& x, int& y) {
  if (b == 0) {
    x = 1, y = 0;
    return;
  }
  exgcd(b, a % b, y, x);
  y -= a / b * x;
}

binpow(x,mod-2)