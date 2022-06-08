#pragma once
#include "FindDNA.h"

std::vector<uint_fast8_t> Alloc_RNA(int Len, std::vector<std::vector<uint_fast8_t>>& CheckAgainst);

std::string RNAAllocToString(std::vector<uint_fast8_t>& Alloc);