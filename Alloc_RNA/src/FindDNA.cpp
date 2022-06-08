#include "FindDNA.h"

#define NucA 0b10
#define NucT 0b11
#define NucC 0b00
#define NucG 0b01

static std::vector<uint_fast8_t> NucTable = std::vector<uint_fast8_t>(256, 0);

NucleotideRegion::NucleotideRegion(const char* subseq)
{
	std::string sub = std::string(subseq);
	this->Contents.reserve(sub.length());
	for (char c : sub)
	{
		switch (c)
		{
		case 'a':
		case 'A':
			Contents.push_back(NucA);
			break;
		case 't':
		case 'T':
			Contents.push_back(NucT);
			break;
		case 'c':
		case 'C':
			Contents.push_back(NucC);
			break;
		case 'g':
		case 'G':
			Contents.push_back(NucG);
			break;
		}
	}

	//Cache the existance of each four-nucleotide pattern
	this->Buckets.resize(256, 0);
	for (int i = 0; i < Contents.size() - 3; ++i)
	{
		this->Buckets[
			(Contents[i] << 6) |
				(Contents[i + 1] << 4) |
				(Contents[i + 2] << 2) |
				(Contents[i + 3])
		] = true;
	}
}

NucleotideRegion::NucleotideRegion(std::string& subseq)
{
	this->Contents.reserve(subseq.length());
	for (char c : subseq)
	{
		switch (c)
		{
		case 'a':
		case 'A':
			Contents.push_back(NucA);
			break;
		case 't':
		case 'T':
			Contents.push_back(NucT);
			break;
		case 'c':
		case 'C':
			Contents.push_back(NucC);
			break;
		case 'g':
		case 'G':
			Contents.push_back(NucG);
			break;
		}
	}

	//Cache the existance of each four-nucleotide pattern
	this->Buckets.resize(256, 0);
	for (int i = 0; i < Contents.size() - 3; ++i)
	{
		uint_fast8_t idx = (Contents[i] << 6) + (Contents[i + 1] << 4) + (Contents[i + 2] << 2) + (Contents[i + 3]);
		this->Buckets[idx] = true;
	}
}

bool FindNucsInSequence(std::vector<uint_fast8_t>& Sequence, NucleotideRegion& Pattern)
{
	//Start at the last element in the pattern
	for (int i = Pattern.Contents.size(); i < Sequence.size(); ++i)
	{
		//Get the four-nucleotide chunk in the sequence
		uint_fast8_t chunk = (Sequence[i]) + (Sequence[i - 1] << 2) + (Sequence[i - 2] << 4) + (Sequence[i - 3] << 6);
		
		//Check the cache for the chunk
		if (Pattern.Buckets[chunk])
		{
			//We know that the chunk exists somewhere in the pattern, move backwards until we find it
			for (int j = 0; j < Pattern.Contents.size() - 3; ++j)
			{
				//Check if the last nucleotide in the pattern matches the current chunk
				if (Sequence[i] != Pattern.Contents[(Pattern.Contents.size() - 1) - j])
					continue;
				if (Sequence[i - 1] != Pattern.Contents[(Pattern.Contents.size() - 2) - j])
					continue;
				if (Sequence[i - 2] != Pattern.Contents[(Pattern.Contents.size() - 3) - j])
					continue;
				if (Sequence[i - 3] != Pattern.Contents[(Pattern.Contents.size() - 4) - j])
					continue;
				
				//Now we know that the pattern ends at location i - j. We need to check if the rest of the pattern matches
				i += j;
				//Move backwards from the end of the pattern until we find a mismatch
				for (int k = 0; k < Pattern.Contents.size(); ++k)
				{
					if (Sequence[i - k] != Pattern.Contents[(Pattern.Contents.size() - 1) - k])
					{
						//Advance to the mismatch and exit
						i += k - 1;
						goto NextChunk;
					}
				}
				return true;
			}
		NextChunk:
			continue;
		}
		else
		{
			//Advance by the length of the pattern
			i += Pattern.Contents.size() - 2;
		}
	}
	return false;
}

/*
Executing the above function by hand to make sure it works
Input: "ATTTTAAACACA", "TAAAC"
Steps:
	1. Advance by 4
	2. Advance by 4
	3. Chunk was found
	4. Move backwards until we find a mismatch
	5. No mismatch found, return
*/

std::vector<uint_fast8_t> StringToDNAVector(std::string str)
{
	std::vector<uint_fast8_t> vec;
	vec.reserve(str.length());
	for (char c : str)
	{
		/*switch (c)
		{
		case 'a':
		case 'A':
			vec.push_back(NucA);
			break;
		case 't':
		case 'T':
			vec.push_back(NucT);
			break;
		case 'c':
		case 'C':
			vec.push_back(NucC);
			break;
		case 'g':
		case 'G':
			vec.push_back(NucG);
			break;
		}*/
		auto val = NucTable[c];
#ifdef NON_DNA_INPUT
		if (val != 0)
			vec.push_back(val);
#else
		vec.push_back(val);
#endif
	}
	return vec;
}

#include <fstream>
void SaveDNAVectors(std::string path, std::vector<std::vector<uint_fast8_t>>& dna)
{
	std::ofstream OutputFile(path);
	size_t size = dna.size();
	OutputFile.write((char*)&size, sizeof(size_t));

	//Write the contents to the file
	for (int i = 0; i < dna.size(); ++i)
	{
		size = dna[i].size();
		OutputFile.write((char*)&size, sizeof(size_t));
		OutputFile.write((char*)dna[i].data(), sizeof(uint_fast8_t) * size);
	}
}

void LoadDNAVectors(std::ifstream& file, std::vector<std::vector<uint_fast8_t>>& Output)
{
	size_t size = 0;
	file.read((char*)&size, sizeof(size_t));
	Output.reserve(size);

	for (int i = 0; i < size; ++i)
	{
		size_t Elem_Size = 0;
		//Read the file into a buffer
		uint_fast8_t* data = new uint_fast8_t[Elem_Size];
		file.read((char*)data, sizeof(uint_fast8_t) * Elem_Size);

		//Now copy the buffer into the vector
		std::vector<uint_fast8_t> elem;
		elem.reserve(Elem_Size);
		for (int j = 0; j < Elem_Size; ++j)
		{
			elem.push_back(data[j]);
		}
		Output.push_back(std::move(elem));
		delete[] data;
	}

}

void InitNucTable()
{
	NucTable[(uint_fast8_t)'a'] = NucA;
	NucTable[(uint_fast8_t)'A'] = NucA;
	
	NucTable[(uint_fast8_t)'t'] = NucT;
	NucTable[(uint_fast8_t)'T'] = NucT;
	
	NucTable[(uint_fast8_t)'c'] = NucC;
	NucTable[(uint_fast8_t)'C'] = NucC;
	
	NucTable[(uint_fast8_t)'g'] = NucG;
	NucTable[(uint_fast8_t)'G'] = NucG;
}
