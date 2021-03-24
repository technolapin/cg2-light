
template<class T>
struct BinTree
{
    std::vector<T> nodes;
    std::vector<std::array<int, 3>> relations; // [mother, child a, child b]

    BinTree() = default;
    BinTree(T root):
        nodes({root}),
        relations({{-1, -1, -1}})
    {}

    int
    push_left(T node, int mother)
    {
        auto i = nodes.size();
        nodes.push_back(node);
        relations.push_back({mother, -1, -1});
        relations[mother][1] = i;
        return i;
    }
    int
    push_right(T node, int mother)
    {
        auto i = nodes.size();
        nodes.push_back(node);
        relations.push_back({mother, -1, -1});
        relations[mother][2] = i;
        return i;
    }

    int
    mother(int child)
    {
        return relations[child][0];
    }

    int
    child_left(int mother)
    {
        return relations[mother][1];
    }

    int
    child_right(int mother)
    {
        return relations[mother][2];
    }
    
};

template<class T>
struct Tree
{
    std::vector<T> nodes;
    std::vector<std::vector<int>> childrn;
    std::vector<int> mothers;
    
    Tree() = default;

    int
    push(T node, int mother)
    {
        auto i = nodes.size();
        nodes.push_back(node);
        childrn.push_back({});
        mothers.push_back(mother);
        return i;
    }

    int
    mother(int child) const
    {
        return mother[child];
    }

    const std::vector<int> &
    children(int mother) const
    {
        return childrn[mother];
    }

};
