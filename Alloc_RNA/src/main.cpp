#include <iostream>

#include "FindDNA.h"
#include "Alloc_RNA.h"

#include <fstream>
#include <streambuf>
#include <algorithm>


std::vector<std::vector<uint_fast8_t>> GetSeqsFromFile(std::string FilePath);

//Usage Alloc_RNA length "Path to the file with the DNA sequences"  [NumExtraSeqs ...]
int main(int argc, char* argv[])
{
	int SeqLen;
	std::string SeqLenStr;

	std::string FilePath;

	if (argc == 1)
	{
		std::cout << "Enter the length of sequence to generate: ";
		std::cin >> SeqLenStr;
	}
	else if (std::string(argv[1]) == "MakeTranscriptFile")
	{
		if (argc < 3)
			return -1;
		FilePath = argv[2];
		auto Vec = GetSeqsFromFile(FilePath);
		SaveDNAVectors(FilePath + ".dat", Vec);
		return 0;
	}
	else
	{
		SeqLenStr = argv[1];
	}

	try
	{
		SeqLen = std::stoi(SeqLenStr);
	}
	catch(std::exception e)
	{
		std::cout << "Error parsing input \"" << SeqLenStr << "\"\n";
		return -1;
	}

	if (SeqLen < 10)
	{
		std::cout << "Error: Sequence length to be generated must be at least 10\n";
		return -1;
	}


	if (argc < 3)
	{
		std::cout << "Enter the path to the file with sequences to check against: ";
		std::cin >> FilePath;
	}
	else
	{
		FilePath = argv[2];
	}

	
	std::vector<std::vector<uint_fast8_t>> Seqs = GetSeqsFromFile(FilePath);

	//see if the user added more sequences to check against
	if (argc > 3)
	{
		int NumExtraSeqs = 0;
		try
		{
			std::string val = argv[3];
			NumExtraSeqs = std::stoi(val);
		}
		catch (std::exception e)
		{
			std::cout << "Error parsing number of extra sequences\n";
			return -1;
		}

		for (int i = 0; i < NumExtraSeqs; ++i)
		{
			std::string seq = argv[i + 4];
			Seqs.push_back(StringToDNAVector(seq));
		}
	}

	auto output = Alloc_RNA(SeqLen, Seqs);
	std::cout << RNAAllocToString(output) << std::endl;

	return 0;
}


#include <filesystem>
std::vector<std::vector<uint_fast8_t>> GetSeqsFromFile(std::string FilePath)
{
	std::vector<std::vector<uint_fast8_t>> Output;
	//If there's a cached version of the file
	if (std::filesystem::exists(FilePath + ".dat"))
	{
		std::ifstream InputFile(FilePath + ".dat");

		LoadDNAVectors(InputFile, Output);
	}
	else
	{
		//Load the file
		std::ifstream InputFile(FilePath);
		std::string Line;
		std::getline(InputFile, Line);
		while (Line != "")
		{
			Output.push_back(StringToDNAVector(Line));

			std::getline(InputFile, Line);
		}
	}
	return Output;
}