#pragma once
#include <vector>

//Wrapper around array with a length field
template<typename T>
class ArrayWrapper
{
public:
	ArrayWrapper() : data(nullptr), length(0) {}
	ArrayWrapper(size_t size)
	{
		this->length = size;
		this->data = new T[size];
	}

	ArrayWrapper(std::vector<T>& origin)
	{
		this->length = origin.size();
		this->data = new T[origin.size()];
		
		//Copy the origin data
		for (size_t i = 0; i < origin.size(); i++)
		{
			this->data[i] = origin[i];
		}
	}

	inline T& operator[](size_t idx) { return this->data[idx]; }

	~ArrayWrapper()
	{
		delete[] this->data;
	}

public:
	T* data;
	size_t length;
};