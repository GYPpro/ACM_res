#ifndef _IN_TEMPLATE_
#include <template_overAll.h>
#endif

// ac
template <class TYPE_NAME>
int subsearch_1(TYPE_NAME tar, vector<TYPE_NAME> &arr) // 第一类二分："非下降"序列，"找到"就返回元素下标，相等元素返回最"左侧"下标，否则返回-1
{
    int l = 0;
    int r = arr.size();
    int ans = -1;
    while (r >= l)
    {
        int m =( l + (r - l) >> 1);
        if (arr[m] >= tar)
        {
            if (arr[m] == tar)
                ans = m;
            r = m - 1;
        }
        else
            l = m + 1;
    }
    return ans;
}

template <class TYPE_NAME>
int subsearch_2(TYPE_NAME tar, vector<TYPE_NAME> &arr) // 第二类二分：取小于或等于tar的最大的元素
{
    int l = 0;
    int r = arr.size();
    int m = 0;
    while (r > l + 1)
    {
        int m =( l + (r - l) >> 1);
        if (arr[m] > tar)
            r = m;
        else
            l = m;
    }
    return l;
}

template <class TYPE_NAME>
class subAns
{
public:
    bool check(TYPE_NAME t);
    TYPE_NAME maxLim;
    TYPE_NAME minLim;
    TYPE_NAME getAns()
    {
        TYPE_NAME l = minLim,
                  r = maxLim + 1;
        while (l + 1 < r)
        {                          // 如果两点不相邻
            TYPE_NAME mid = (l + r) / 2; // 取中间值
            if (check(mid))        // 如果可行
                l = mid;           // 升高锯片高度
            else
                r = mid; // 否则降低锯片高度
        }
        return l; // 返回左边值
    }
};

int getAns(auto check())
{
    int l = minLim,
                r = maxLim + 1;
    while (l + 1 < r)
    {                          // 如果两点不相邻
        int mid = (l + r) / 2; // 取中间值
        if (check(mid))        // 如果可行
            l = mid;           // 升高锯片高度
        else
            r = mid; // 否则降低锯片高度
    }
    return l; // 返回左边值
}

// FIXME:二分tm死了
template <class TYPE_NAME>
int subsearch_3(TYPE_NAME tar, vector<TYPE_NAME> &arr) // 第三类二分：取大于或等于tar的最小的元素
{
    int l = 0;
    int r = arr.size();
    int m = 0;
    while (r > l)
    {
        int m = (l + (r - l) >> 1);
        if (arr[m] > tar)
            r = m;
        else
            l = m;
    }
    return r;
}
