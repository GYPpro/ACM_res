#define int long long
const int N = 1 << 21;
static const int mod1 = 1E9 + 7, base1 = 127;
static const int mod2 = 1E9 + 9, base2 = 131;
vector<int> val1;
vector<int> val2;
using puv = pair<int,int>;
void init(int n = N)
{
    val1.resize(n + 1), val2.resize(n + 2);
    val1[0] = 1, val2[0] = 1;
    for (int i = 1; i <= n; i++)
    {
        val1[i] = val1[i - 1] * base1 % mod1;
        val2[i] = val2[i - 1] * base2 % mod2;
    }
}
class hstring
{
public:
    vector<int> h1;
    vector<int> h2;
    string s;

    hstring(string s_) : s(s_), h1{0}, h2{0}
    {
        build();
    }

    void build()
    {
        for (auto it : s)
        {
            h1.push_back((h1.back() * base1 % mod1 + it) % mod1);
            h2.push_back((h2.back() * base2 % mod2 + it) % mod2);
        }
    }

    puv get()
    { // 输出整串的哈希值
        return {h1.back(), h2.back()};
    }

    puv substring(int l, int r)
    { // 输出子串的哈希值
        if (l > r) swap(l, r);
        int ans1 = (mod1 + h1[r + 1] - h1[l] * val1[r - l + 1] % mod1) % mod1;
        int ans2 = (mod2 + h2[r + 1] - h2[l] * val2[r - l + 1] % mod2) % mod2;
        return {ans1, ans2};
    }

    puv modify(int idx, char x) 
    { //修改 idx 位为 x
        int n = s.size() - 1;
        int ans1 = (h1.back() + val1[n - idx] * (x - s[idx]) % mod1) % mod1;
        int ans2 = (h2.back() + val2[n - idx] * (x - s[idx]) % mod2) % mod2;
        return {ans1, ans2};
    }
};

