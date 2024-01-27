#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <math.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <queue>

using namespace std;

#define ll long long
#define pb push_back
#define ld long double
const ll int maxn = 1E5+10;
const ll int mod1 = 998244353;
const ll int mod2 = 1E9+7;

#define _IN_TEMPLATE_

ll int str2int(string s)
{
    ll int rec = 0;
    ll int pw = 1;
    for(int i = s.length()-1;i >= 0;i --)
    {
        int gt = s[i] - '0';
        if(gt < 0 || gt > 9) return INT64_MAX;
        rec += gt * pw;
        pw *= 10;
    }
    return rec;
}

vector<ll int> testReadLine()
{
    string s;
    getline(cin,s);
    s.push_back(' ');
    vector<ll int> rearr;
    vector<string> substring;
    string ts;
    for(int i = 0;i < s.size();i ++)
    {
        if(s[i] == ' '){
            substring.push_back(ts);
            ts.clear();
        } else ts.push_back(s[i]);
    }
    for(int i = 0;i < substring.size();i ++)rearr.push_back(str2int(substring[i]));
    return rearr;
}