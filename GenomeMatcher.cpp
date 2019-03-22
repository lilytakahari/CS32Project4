#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#include "Trie.h"
#include <map>
#include <algorithm>
struct SeqIndex
{
    unsigned long genomePos;
    unsigned long subSeqPos;
    SeqIndex(unsigned long gPos, unsigned long sPos)
    {
        genomePos = gPos;
        subSeqPos = sPos;
    }
};

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int m_minSearchLength;
    vector<Genome> m_genomeLibrary;
    Trie<SeqIndex> m_seqTrie;
};
bool orderGenomeIndices(SeqIndex &a, SeqIndex &b)
{
    if (a.genomePos < b.genomePos)
        return true;
    else if (a.genomePos == b.genomePos)
        return (a.subSeqPos < b.subSeqPos);
    return false;
}
bool orderRelatedGenomes(GenomeMatch &a, GenomeMatch &b)
{
    if (a.percentMatch > b.percentMatch)
        return true;
    else if (a.percentMatch == b.percentMatch)
        return (a.genomeName < b.genomeName);
    return false;
}
GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    m_genomeLibrary.push_back(genome);
    unsigned long genomePos = m_genomeLibrary.size() - 1;
    for (unsigned long i = 0; i < genome.length() - m_minSearchLength; i++)
    {
        string key;
        genome.extract((int)i, m_minSearchLength, key);
        m_seqTrie.insert(key, SeqIndex(genomePos, i));
    }
    
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    if (fragment.size() < minimumLength || minimumLength < m_minSearchLength)
        return false;
    matches.clear();
    vector<SeqIndex> result = m_seqTrie.find(fragment.substr(0, m_minSearchLength), exactMatchOnly);
    if (result.empty())
        return false;
    sort(result.begin(), result.end(), orderGenomeIndices);
    unsigned long currentGenomePos = result[0].genomePos;
    unsigned long largestFoundLengthForCurrGenome = 0;
    unsigned long correspondingPosition = 0;
    for (int i = 0; i < result.size(); i++)
    {
        unsigned long iteratingGenomePos = result[i].genomePos;
        if (currentGenomePos != iteratingGenomePos)
        {
            if (largestFoundLengthForCurrGenome != 0)
            {
                DNAMatch toAdd;
                toAdd.genomeName = m_genomeLibrary[currentGenomePos].name();
                toAdd.length = (int)largestFoundLengthForCurrGenome;
                toAdd.position = (int)correspondingPosition;
                matches.push_back(toAdd);
            }
            currentGenomePos = iteratingGenomePos;
            largestFoundLengthForCurrGenome = 0;
            correspondingPosition = 0;
        }
        string potentialMatch;
        m_genomeLibrary[iteratingGenomePos].extract((int)(result[i].subSeqPos + m_minSearchLength), (int)fragment.size() - m_minSearchLength, potentialMatch);
        string::iterator mismatchPos = mismatch(potentialMatch.begin(), potentialMatch.end(), fragment.begin() + m_minSearchLength).first;
        unsigned long matchingLength = (mismatchPos - potentialMatch.begin()) + m_minSearchLength;
        if (!exactMatchOnly)
        {
            string checkSnip;
            m_genomeLibrary[iteratingGenomePos].extract((int)(result[i].subSeqPos), m_minSearchLength, checkSnip);
            if (checkSnip == fragment.substr(0, m_minSearchLength))
            {
                int snipCounter = 1;
                int index = (int)result[i].subSeqPos + m_minSearchLength;
                int move = 0;
                while (snipCounter >= 0 && (m_minSearchLength + move) < fragment.size())
                {
                    string currChar;
                    m_genomeLibrary[iteratingGenomePos].extract(index, 1, currChar);
                    if (fragment[m_minSearchLength + move] != currChar[0])
                        snipCounter--;
                    if (snipCounter < 0)
                        break;
                    move++;
                    index++;
                }
                matchingLength = index - result[i].subSeqPos;
            }
            
        }
        if (matchingLength >= minimumLength && matchingLength > largestFoundLengthForCurrGenome)
        {
            largestFoundLengthForCurrGenome = matchingLength;
            correspondingPosition = result[i].subSeqPos;
        }
    }
    if (largestFoundLengthForCurrGenome != 0)
    {
        DNAMatch toAdd;
        toAdd.genomeName = m_genomeLibrary[currentGenomePos].name();
        toAdd.length = (int)largestFoundLengthForCurrGenome;
        toAdd.position = (int)correspondingPosition;
        matches.push_back(toAdd);
    }
    if (matches.empty())
        return false;
    return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength < m_minSearchLength)
        return false;
    results.clear();
    map<string, int> genomeToMatchSeqNum;
    int numSequences = query.length() / fragmentMatchLength;
    for (int i = 0; i < numSequences; i++)
    {
        string extractSeq;
        query.extract(i*fragmentMatchLength, fragmentMatchLength, extractSeq);
        vector<DNAMatch> matched;
        findGenomesWithThisDNA(extractSeq, fragmentMatchLength, exactMatchOnly, matched);
        for (int j = 0; j < matched.size(); j++)
        {
            pair<map<string,int>::iterator,bool> insert = genomeToMatchSeqNum.insert({matched[j].genomeName, 1});
            if (insert.second == false)
                (insert.first)->second++;
        }
    }
    if (genomeToMatchSeqNum.empty())
        return false;
    map<string, int>::iterator it;
    for (it = genomeToMatchSeqNum.begin(); it != genomeToMatchSeqNum.end(); it++)
    {
        double pct = ((double)it->second / numSequences) * 100;
        if (pct >= matchPercentThreshold)
        {
            GenomeMatch g;
            g.genomeName = it->first;
            g.percentMatch = pct;
            results.push_back(g);
        }
    }
    if (results.empty())
        return false;
    sort(results.begin(), results.end(), orderRelatedGenomes);
    return true;  // This compiles, but may not be correct
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
