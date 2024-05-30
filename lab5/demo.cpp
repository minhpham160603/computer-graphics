#include <unordered_set>

int main()
{
    std::unordered_set<std::pair<int, int>> p;
    p.insert({4, 5});
    return (p.find({3, 4}) != p.end());
}