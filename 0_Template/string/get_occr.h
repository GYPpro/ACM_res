#include <template_overAll.h>

/*
 * 找到某一堆短字符串在长字符串中的出现位置
 * dira=1为最早出现的后端点下标 dira=0为最晚出现的前端点下标
 * 源字符串s长度为|s|，查找字符串列表中所有字符串长度和为|_s|
 * 则时间复杂度为O(max(|_s|log(|_s|),|s|))
 */
class get_occr
{
private:
    string s;
public:
    get_occr(string _s) { s = _s; }
    vector<int> locate(vector<string> _s,bool dira = 1)
    {
        int n = _s.size();
        vector<int> occr(n,-1);
        map<char,vector<pair<int,int>>> gncing;
        if(dira == 1)
        {
            for(int i = 0;i < n;i++)
                gncing[_s[i][0]].push_back({i,0});
            for(int i = 0;i < s.size();i ++)
            {
                vector<pair<int,int>> gnctmp = gncing[s[i]];
                gncing[s[i]].clear();
                for(int j = 0;j < gnctmp.size();j ++)
                {
                    if(gnctmp[j].se+1 < _s[gnctmp[j].fi].size())
                            gncing[_s[gnctmp[j].fi][gnctmp[j].se+1]].push_back({gnctmp[j].fi,gnctmp[j].se+1});
                    else occr[gnctmp[j].fi] = i;
                }
            }
        } else {
            for(int i = 0;i < n;i++) gncing[_s[i][_s[i].size()-1]].push_back({i,_s[i].size()-1});
            for(int i= s.size()-1;i >=0;i --)
            {
                vector<pair<int,int>> gnctmp = gncing[s[i]];
                gncing[s[i]].clear();
                for(int j = 0;j < gnctmp.size();j ++)
                {
                    if(gnctmp[j].se -1 >= 0)
                            gncing[_s[gnctmp[j].fi][gnctmp[j].se-1]].push_back({gnctmp[j].fi,gnctmp[j].se-1});
                    else occr[gnctmp[j].fi] = i;
                }
            }
        }
        return occr;
    }
};
