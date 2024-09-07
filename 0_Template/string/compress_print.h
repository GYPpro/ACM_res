
const int N = 1 << 21;
static const int mod1 = 1E9 + 7, base1 = 127;
static const int mod2 = 1E9 + 9, base2 = 131;
vector<int> val1;
vector<int> val2;
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

string compress(vector<string> in)
{ // 前后缀压缩
    vector<int> h1{1};
    vector<int> h2{1};
    string ans = "#";
    for (auto s : in)
    {
        s = "#" + s;
        int st = 0;
        int chk1 = 0;
        int chk2 = 0;
        for (int j = 1; j < s.size() && j < ans.size(); j++)
        {
            chk1 = (chk1 * base1 % mod1 + s[j]) % mod1;
            chk2 = (chk2 * base2 % mod2 + s[j]) % mod2;
            if ((h1.back() == (h1[ans.size() - 1 - j] * val1[j] % mod1+ chk1) % mod1) &&
                (h2.back() == (h2[ans.size() - 1 - j] * val2[j] % mod2+ chk2) % mod2)    ) 
                st = j;
        }
        for (int j = st + 1; j < s.size(); j++)
        {
            ans += s[j];
            h1.push_back((h1.back() * base1 % mod1 + s[j]) % mod1);
            h2.push_back((h2.back() * base2 % mod2 + s[j]) % mod2);
        }
    }
    return ans.substr(1);
}

