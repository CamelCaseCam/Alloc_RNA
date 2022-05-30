#pragma once
#include <vector>
#include <string>

//With a bucket size of 4

class NucleotideRegion
{
public:
	NucleotideRegion(const char* subseq);
	NucleotideRegion(std::string& subseq);

	std::vector<uint_fast8_t> Contents;
	std::vector<char> Buckets;
};

bool FindNucsInSequence(std::vector<uint_fast8_t>& Sequence, NucleotideRegion& Pattern);

std::vector<uint_fast8_t> StringToDNAVector(std::string str);

void SaveDNAVectors(std::string path, std::vector<std::vector<uint_fast8_t>>& dna);
void LoadDNAVectors(std::ifstream& file, std::vector<std::vector<uint_fast8_t>>& Output);

void InitNucTable();