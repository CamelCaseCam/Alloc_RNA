#include "Alloc_RNA.h"
#include <thread>
#include <random>

NucleotideRegion CreateSequence(int Len);



#ifdef MultiThread
bool CheckSequence(NucleotideRegion* seq, std::vector<std::thread>& Threads, std::vector<std::vector<uint_fast8_t>*>& CurrentSequenceChecks, 
	std::vector<char>* Results, std::vector<std::vector<uint_fast8_t>>& CheckAgainst, bool* exit);

void CheckIndividualRegion(NucleotideRegion* Subseq, std::vector<uint_fast8_t>** Region, char* Result, bool* exit);
#endif

std::mt19937 rnd;

#ifdef MultiThread
std::vector<uint_fast8_t> Alloc_RNA(int Len, std::vector<std::vector<uint_fast8_t>>& CheckAgainst)
{
	//This is multithreaded, so detect the number of threads available
	int NumThreads = std::thread::hardware_concurrency();
	if (NumThreads == 0)
	{
		NumThreads = 1;
	}
	else if (NumThreads > CheckAgainst.size())
	{
		NumThreads = CheckAgainst.size();
	}
	
	//Create a vector of threads
	std::vector<std::thread> Threads;
	//Create a vector of pointers to the sequences to check against. The threads will use these to communicate with the main thread
	std::vector<std::vector<uint_fast8_t>*> CurrentSequenceChecks;
	CurrentSequenceChecks.reserve(CheckAgainst.size());
	std::vector<char> Results(CheckAgainst.size());
	
	//Seed the random number generator
	rnd.seed(std::random_device()());

	NucleotideRegion Sequence = CreateSequence(Len);
	bool exit = false;
	
	//Set up the threads for checking the sequences
	for (int i = 0; i < NumThreads; i++)
	{
		//Add the vector to the vector of vectors
		CurrentSequenceChecks.push_back(nullptr);
		//Create a new thread
		Threads.push_back(std::thread(CheckIndividualRegion, &Sequence, &CurrentSequenceChecks[i], &Results[i], &exit));
	}
	
	while (!CheckSequence(&Sequence, Threads, CurrentSequenceChecks, &Results, CheckAgainst, &exit))
	{
		//If the sequence is not valid, create a new one
		Sequence = CreateSequence(Len);
	}
	return std::move(Sequence.Contents);
}


#endif
char RandomNucleotide();

NucleotideRegion CreateSequence(int Len)
{
	//Fill vector with Len random nucleotides
	std::string Output;
	Output.reserve(Len);
	for (int i = 0; i < Len; i++)
	{
		Output += RandomNucleotide();
	}
	
	//Convert the vector to a NucleotideRegion
	return NucleotideRegion(Output);
}


char RandomNucleotide()
{
	//Randomly return A, C, G, or T without using rand()
	#ifdef _msc_ver
	std::uniform_int<unsigned int> gen(0, 3);
	#else
	std::uniform_int_distribution<unsigned int> gen(0, 3);
	#endif
	
	return "ACGT"[gen(rnd)];
}


#ifdef MultiThread
//Use the threads to check the sequences
bool CheckSequence(NucleotideRegion* seq, std::vector<std::thread>& Threads, std::vector<std::vector<uint_fast8_t>*>& CurrentSequenceChecks,
	std::vector<char>* Results, std::vector<std::vector<uint_fast8_t>>& CheckAgainst, bool* exit)
{
	//Assign all the threads their first sequence
	int SeqNum = 0;
	for (SeqNum; SeqNum < Threads.size(); SeqNum++)
	{
		CurrentSequenceChecks[SeqNum] = CurrentSequenceChecks[SeqNum];
	}
	
	bool Finished = false;
	while (!Finished)
	{
		//Check on all the threads and reassign them new sequences once they finish
		for (int i = 0; i < Threads.size(); i++)
		{
			//Check if the thread has finished
			if (CurrentSequenceChecks[i] == nullptr)
			{
				//Check if the result was true
				if ((*Results)[i])
				{
					//Tell all the threads to stop
					for (int i = 0; i < CurrentSequenceChecks.size(); ++i)
					{
						CurrentSequenceChecks[i] = (std::vector<uint_fast8_t>*)0x1;
					}
					return false;
				}
				else
				{
					//Give the thread a new task unless there are no more sequences to check
					if (SeqNum < CheckAgainst.size())
					{
						CurrentSequenceChecks[i] = &CheckAgainst[SeqNum];
						SeqNum++;
					}
					//If there are no more sequences to check, tell the thread to stop
					else
					{
						*exit = true;
						Finished = true;
						break;
					}
				}
			}
		}
		std::this_thread::yield();
	}
	//Join all threads
	for (int i = 0; i < Threads.size(); i++)
	{
		Threads[i].join();
	}
	
	return true;
}


void CheckIndividualRegion(NucleotideRegion* Subseq, std::vector<uint_fast8_t>** Region, char* Result, bool* exit)
{
	//0x1 will be the signal to exit all threads
	while (!exit)
	{
		//if it's nullptr, wait for a new region to be assigned
		if (Region == nullptr)
		{
			std::this_thread::yield();
			continue;
		}
		
		//Execute the check
		*Result = FindNucsInSequence(**Region, *Subseq);
		Region = nullptr;
	}
}

#else

bool CheckSequence(NucleotideRegion& Candidate, std::vector<std::vector<uint_fast8_t>>& CheckAgainst);
std::vector<uint_fast8_t> Alloc_RNA(int Len, std::vector<std::vector<uint_fast8_t>>& CheckAgainst)
{
	NucleotideRegion candidate = CreateSequence(Len);
	
	while (!CheckSequence(candidate, CheckAgainst))
	{
		candidate = CreateSequence(Len);
	}
	
	return std::move(candidate.Contents);
}

bool CheckSequence(NucleotideRegion& Candidate, std::vector<std::vector<uint_fast8_t>>& CheckAgainst)
{
	for (std::vector<uint_fast8_t>& CheckAgainstRegion : CheckAgainst)
	{
		if (FindNucsInSequence(CheckAgainstRegion, Candidate))
		{
			return false;
		}
	}
	return true;
}

#endif

std::string RNAAllocToString(std::vector<uint_fast8_t>& Alloc)
{
	std::string Output;
	Output.reserve(Alloc.size());
	for (int i = 0; i < Alloc.size(); i++)
	{
		switch (Alloc[i])
		{
		case 0b10:
			Output += 'A';
			break;
		case 0b11:
			Output += 'T';
			break;
		case 0b00:
			Output += 'C';
			break;
		case 0b01:
			Output += 'G';
			break;
		}
	}

	return Output;
}