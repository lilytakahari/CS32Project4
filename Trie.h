#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

template<typename ValueType>
struct TrieNode;

template<typename ValueType>
struct TrieNode
{
    std::vector<ValueType> m_values;
    std::vector<char> m_keys;
    std::vector<TrieNode<ValueType>*> m_children;
};

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    void resetHelper(TrieNode<ValueType> *curr);
    void insertHelper(const std::string& key,
                      const ValueType& value, TrieNode<ValueType> *curr);
    std::vector<ValueType> findHelper(const std::string& key,
                                      bool exactMatchOnly, TrieNode<ValueType> *curr) const;
    TrieNode<ValueType>* m_root;
};

template<typename ValueType>
Trie<ValueType>::Trie()
{
    m_root = new TrieNode<ValueType>;
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
    reset();
    delete m_root;
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
    resetHelper(m_root);
    delete m_root;
    m_root = new TrieNode<ValueType>;
}

template<typename ValueType>
void Trie<ValueType>::resetHelper(TrieNode<ValueType> *curr)
{
    for (int i = 0; i < curr->m_children.size(); i++)
    {
        resetHelper(curr->m_children[i]);
        delete curr->m_children[i];
    }
}
template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
    insertHelper(key, value, m_root);
}

template<typename ValueType>
void Trie<ValueType>::insertHelper(const std::string& key,
                                   const ValueType& value, TrieNode<ValueType> *curr)
{
    if (key.empty())
    {
        (curr->m_values).push_back(value);
        return;
    }
    int i;
    for (i = 0; i < curr->m_keys.size(); i++)
    {
        if (key[0] == curr->m_keys[i])
        {
            insertHelper(key.substr(1), value, curr->m_children[i]);
            return;
        }
    }
    curr->m_children.push_back(new TrieNode<ValueType>);
    curr->m_keys.push_back(key[0]);
    insertHelper(key.substr(1), value, curr->m_children[i]);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    for (int i = 0; i < m_root->m_children.size(); i++)
    {
        if (key[0] == m_root->m_keys[i])
            return findHelper(key.substr(1), exactMatchOnly, m_root->m_children[i]);
    }
    std::vector<ValueType> empty;
    return empty;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findHelper(const std::string& key,
                                  bool exactMatchOnly, TrieNode<ValueType> *curr) const
{
    std::vector<ValueType> toReturn;
    if (key.empty())
    {
        toReturn.insert(toReturn.end(), curr->m_values.begin(), curr->m_values.end());
        return toReturn;
    }
    
    for (int i = 0; i < curr->m_children.size(); i++)
    {
        if (exactMatchOnly && (key[0] == curr->m_keys[i]))
            return findHelper(key.substr(1), exactMatchOnly, curr->m_children[i]);
        else if (!exactMatchOnly)
        {
            std::vector<ValueType> toAdd;
            if (key[0] == curr->m_keys[i]) {
                toAdd = findHelper(key.substr(1), false, curr->m_children[i]);
            } else {
                toAdd = findHelper(key.substr(1), true, curr->m_children[i]);
            }
            toReturn.insert(toReturn.end(), toAdd.begin(), toAdd.end());
        }
    }
    return toReturn;
}
#endif // TRIE_INCLUDED
