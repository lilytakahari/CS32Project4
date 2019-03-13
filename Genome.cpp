#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <cctype>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
    genomes.clear();
    string currLine, nm, seq;

    // First line checking
    string firstLine;
    if (!getline(genomeSource, firstLine))
        return false;
    if (firstLine.empty())
        return false;
    if (firstLine[0] != '>')
        return false;
    nm = firstLine.substr(1);
    
    // Run through the file
    while (getline(genomeSource, currLine))
    {
        if (currLine.empty())
            return false;
        if (currLine[0] == '>')
        {
            if (seq.empty())
                return false;
            genomes.push_back(Genome(nm, seq));
            nm = currLine.substr(1);
            seq = "";
            continue;
        }
        else
        {
            if (nm.empty())
                return false;
            for (int i = 0; i < currLine.size(); i++)
            {
                if (currLine[i] != 'N' && currLine[i] != 'n' &&
                    currLine[i] != 'C' && currLine[i] != 'c' &&
                    currLine[i] != 'T' && currLine[i] != 't' &&
                    currLine[i] != 'G' && currLine[i] != 'g' &&
                    currLine[i] != 'A' && currLine[i] != 'a')
                    return false;
                currLine[i] = toupper(currLine[i]);
            }
            seq += currLine;
        }
    }
    if (seq.empty())
        return false;
    genomes.push_back(Genome(nm, seq));
    return true;
}

int GenomeImpl::length() const
{
    return (int) m_sequence.size();
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if (length <= 0 || position < 0 || position >= m_sequence.size()
        || position + length > m_sequence.length())
        return false;
    fragment = m_sequence.substr(position, length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
