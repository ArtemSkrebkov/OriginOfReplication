#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <set>
#include "stdio.h"

using namespace std;

struct InfoFasta
{
	string info;
	string genome;
};

struct InfoDnaABox
{
	int occur;
	string dnaBoxes;
	friend bool operator<(const InfoDnaABox &op1, const InfoDnaABox &op2)
	{
		return op1 < op2;
	}
};

vector<int> GetSkew(string &text)
{
	vector<int> result;
	result.resize(text.size() + 1);
	result[0] = 0;
	
	for(int i = 0; i < text.size(); i++)
	{
		switch(text[i])
		{
		case 'C': result[i + 1] = result[i] - 1; break;
		case 'G': result[i + 1] = result[i] + 1; break;
		case 'c': result[i + 1] = result[i] - 1; break;
		case 'g': result[i + 1] = result[i] + 1; break;
		case 'a':
		case 't':
		case 'A':
		case 'T': result[i + 1] = result[i]; break;
		}
		//if(i % 500 == 0)
		//	out << i << "," << result[i] << "," << endl;
		
	}

	return result;
}

vector<int> GetPosMinSkew(vector<int> &skew)
{
	vector<int> result;

	int min = INT_MAX;

	for(int i = 0; i < skew.size(); i++)
		if(skew[i] < min)
		{
			min = skew[i];
			result.clear();
			result.push_back(i);
		}
		else
			if(skew[i] == min)
				result.push_back(i);

	return result;
}

void ReadFasta(ifstream &in, vector<InfoFasta> &infoGenome)
{
	
	InfoFasta record;
	char *buffer = new char[1000];
	int countRecord = 0;
	while(!in.eof())
	{
		in.getline(buffer, 1000);
		string tmp(buffer);
		if(tmp != "")
		{
			if(tmp[0] == '>')
			{
				if(countRecord == 0)
				{
					countRecord++;
				}
				else
				{
					infoGenome.push_back(record);
					record.info.clear();
					record.genome.clear();
				}
				record.info += tmp;
			}
			else
				record.genome += tmp;
		}
	}
	infoGenome.push_back(record);
}

int GetCountMistmatch(string &pattern, string &text)
{
	int result = 0;

	if(pattern.size() == text.size())
	{
		for(int i = 0; i < pattern.size(); i++)
			if(pattern[i] != text[i])
				result++;
	}
	else
		result = -1;

	return result;
}

int GetCountApproximatePattern(string &pattern, string &text, int d, vector<int> &approximatePatternPos)
{
	int result = 0;

	//cout << "pattern - " << pattern << " " << "approximate - ";
	approximatePatternPos.clear();
	for(int i = 0; i + pattern.size() < text.size(); i++)
	{
		string tmp = text.substr(i, pattern.size());
		int countMistmatch = GetCountMistmatch(pattern, tmp);
		if(countMistmatch != -1 && countMistmatch <= d)
		{
			//cout  << tmp << " ";
			approximatePatternPos.push_back(i);
			result++;
		}
	}
	
	//cout << endl;

	return result;
}

string GetReverseComplement(string &dna)
{
	string result;

	for(int i = dna.size() - 1; i >= 0; i--)
	{
		switch(dna[i])
		{
		case 'A': result.push_back('T'); break;
		case 'T': result.push_back('A'); break;
		case 'C': result.push_back('G'); break;
		case 'G': result.push_back('C'); break;
		}
	}

	return result;
}

vector<string> GetMostFrequentKmerWithMistmatchAndReverse(string &text, int k, int d, int &max)
{
	vector<string> result;
	//int max = 0;
	max = 0;
	for(int i = 0; i + k < text.size(); i++)
	{
		string pattern = text.substr(i, k);
		string reversePattern = GetReverseComplement(pattern);
		vector<int> appPattern;
		vector<int> appPatternReverse;
		int countApproximatePattern = GetCountApproximatePattern(pattern, text, d, appPattern);
		int countApproximateReversePattern = GetCountApproximatePattern(reversePattern, text, d, appPatternReverse);
		//int countApproximatePattern = GetCountApproximatePattern(pattern, text, d) + GetCountApproximatePattern(reversePattern, text, d);
		if(countApproximatePattern + countApproximateReversePattern > max /*&& countApproximateReversePattern == countApproximatePattern*/)
		{
			max = countApproximatePattern + countApproximateReversePattern;
			result.clear();
			//InfoDnaABox tmp;
			//tmp.occur = max;
			//tmp.dnaBoxes = pattern;
			result.push_back(pattern);
			//result.insert(result.end(), appPattern.begin(), appPattern.end());
			//result.insert(result.end(), appPatternReverse.begin(), appPatternReverse.end());
			//result.push_back(reversePattern);
		}
		else 
			if(max == countApproximatePattern + countApproximateReversePattern )
			{
				bool flag = false;
				for(int i = 0; i < result.size(); i++)
					if(result[i] == pattern)
					{
						flag = true;
						break;
					}
					if(!flag)
					{
						//InfoDnaABox tmp;
						//tmp.occur = max;
						//tmp.dnaBoxes = pattern;
						result.push_back(pattern);
						//result.push_back(reversePattern);
					}
			}
            /*
			else
				if(max - 1 == countApproximatePattern + countApproximateReversePattern)
				{
					bool flag = false;
					for(int i = 0; i < result.size(); i++)
						if(result[i] == pattern)
						{
							flag = true;
							break;
						}
						if(!flag)
						{
							//InfoDnaABox tmp;
							//tmp.occur = max - 1;
							//tmp.dnaBoxes = pattern;
							result.push_back(pattern);
							//result.push_back(reversePattern);
						}
				}
				else
					if(max - 2 == countApproximatePattern + countApproximateReversePattern)
				{
					bool flag = false;
					for(int i = 0; i < result.size(); i++)
						if(result[i] == pattern)
						{
							flag = true;
							break;
						}
						if(!flag)
						{
							//InfoDnaABox tmp;
							//tmp.occur = max - 2;
							//tmp.dnaBoxes = pattern;
							result.push_back(pattern);
							//result.push_back(reversePattern);
						}
				}
                */
	}

	return result;
}
/*
string UpperCase(string &genome)
{
	string result;
	for(int i = 0; i < genome.size(); i++)
		result.push_back(toupper(genome[i]));

	return result;
}
/*
void PositionAndFreq()
{
	cout << "Input filename";
	string filename;
	cin >> filename;
	ifstream in(filename);

	int d, k, countKmer, isFileSave;

	in >> d >> k >> countKmer >> isFileSave;
	vector<string> kmers;
	kmers.resize(countKmer);
	for(int i = 0; i < countKmer; i++)
		in >> kmers[i];

	string genome;
	while(!in.eof())
	{
		string tmp;
		in >> tmp;
		genome += tmp;
	}
	genome = UpperCase(genome);
	ofstream out("out");
	for(int i = 0; i < countKmer; i++)
	{
		vector<int> appPatternsPosition;
		GetCountApproximatePattern(kmers[i], genome, d, appPatternsPosition);
		cout << "Kmers " << kmers[i] << " occur " << appPatternsPosition.size() << " times with " << d << " mistmatch in position" << endl;
		for(int j = 1; j < appPatternsPosition.size(); j++)
			cout << appPatternsPosition[j] << " " << appPatternsPosition[j] - appPatternsPosition[j - 1] << endl;
		if(isFileSave)
		{
			
			out << "Kmers " << kmers[i] << " occur " << appPatternsPosition.size() << " times with " << d << " mistmatch in position" << endl;
			for(int j = 1; j < appPatternsPosition.size(); j++)
				out << appPatternsPosition[j] << " " << appPatternsPosition[j] - appPatternsPosition[j - 1] << endl;
		}
	}
	in.close();
}
*/


int main(int argc, char *argv[])
{
	
	//ifstream in(argv[1]);
	ofstream out("Result.txt");
	

    cout << "Input filename -> ";
    string fileName;
    cin >> fileName;

    ifstream in(fileName);
	vector<InfoFasta> info;
	ReadFasta(in, info);

    cout << "Count of sequences is " << info.size() << endl;
    cout << "Input number of sequence -> ";
    int num;
    cin >> num;
    vector<int> skew = GetSkew(info[num].genome);

	vector<int> posMinSkew = GetPosMinSkew(skew);
    /*
	for(int i = 0; i < posMinSkew.size(); i++)
		cout << posMinSkew[i] << " ";
	cout << endl;
	*/
	//for(int i = 0; i < posMinSkew.size(); i++)
	//	out << posMinSkew[i] << " ";
    cout << "Input shift, width window, count mismatch -> ";
    int shift, width, countMismatch;
    cin >> shift >> width >> countMismatch;
	
	set<vector<string>> resClump;
	
	for(int i = posMinSkew[0] - shift; i <= posMinSkew[0] + shift; i++)
	{
		string tmp = info[0].genome.substr(i, width);
		int max = 0;
        vector<string> kMerWithMismatch = GetMostFrequentKmerWithMistmatchAndReverse(tmp, 9, countMismatch, max);
		if(max >= 4)
		{
			resClump.insert(kMerWithMismatch);
			
			out << i << endl;
			out << max << endl;
			for(int j = 0; j < kMerWithMismatch.size(); j++)
				out << kMerWithMismatch[j] << " ";
			out << endl << endl;
		}
	}
	/*
	for(set<vector<string>>::iterator i = resClump.begin(); i != resClump.end(); i++)
	{
		for(int j = 0; j < (*i).size(); j++)
			out << (*i)[j] << " ";
		out << endl;
	}
	*/
	/*
	string genome;
	while(!in.eof())
	{
		string tmp;
		in >> tmp;
		genome += tmp;
	}
	genome = UpperCase(genome);
	
	string pattern = "TTATACAAA";
	string reversePattern = GetReverseComplement(pattern);
	vector<int> appPattern;
	vector<int> appPatternReverse;
	int countApproximatePattern = GetCountApproximatePattern(pattern, genome, 2, appPattern);
	int countApproximateReversePattern = GetCountApproximatePattern(reversePattern, genome, 2, appPatternReverse);
	cout << countApproximatePattern + countApproximateReversePattern << endl;
	cout << "DnaA box: " << pattern << endl;
	cout << "Reverse DnaA box: " << reversePattern << endl; 
	*/
	/*
	int max = 10;
	vector<InfoDnaABox> kMerWithMismatch = GetMostFrequentKmerWithMistmatchAndReverse(genome, 9, 2, max);
	cout << max << endl;
	for(int i = 0; i < kMerWithMismatch.size(); i++)
		out << kMerWithMismatch[i].dnaBoxes << " " << kMerWithMismatch[i].occur << endl;
	*/
	/*
	in.close();
	out.close();

	*/
	//PositionAndFreq();
	system("pause");
	return 0;
}