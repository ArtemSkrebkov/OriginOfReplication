#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>

using namespace std;

enum SEQ_TYPE {ORI_CONFIRMED, ORI_LIKELY, ORI_DUBIOUS, USUAL}; 

struct Point
{
    SEQ_TYPE type;
    unordered_map<string, float> freqDinucleotide;
};

struct DistInfo
{
    float dist;
    SEQ_TYPE type;
};


struct InfoFasta
{
	string info;
	string genome;
};

struct WindowInfo
{
    int start, finish;
    int numChr;
    SEQ_TYPE type;
};

bool CmpDistInfo(DistInfo &op1, DistInfo &op2)
{
    bool result = true;

    if(op1.dist >= op2.dist)
        result = false;

    return result;
}


class kNN
{
public:
    kNN(int _countNucleotide = 2, int _windowSize = 200, int _countOriWindow = 330, int _countUsualWindow = 330);
    void Training(vector<string> fileNamesChr, string fileNameInfoOri); 
    
    SEQ_TYPE Classification(string filenameSequence, int k);

    void LoadTrainingInfo(string fileName);
    void SaveTrainingInfo(string fileName);
private:
    int countNucleotide;
    int windowSize;
    int countOriWindow;
    int countUsualWindow;
    unordered_map<string, float> CalcFrequency(string &sequence);
    void InitPoint(unordered_map<string, float> &point);
    void ReadOriInfo(string fileNameInfoOri, vector<WindowInfo> &_confirmedOriWindows, vector<WindowInfo> &_usualWindows);
    void CalcTrainingInfo(vector<WindowInfo> &windows, int countElement);
    void ReadFasta(ifstream &in, vector<InfoFasta> &infoGenome);
    float Distance(Point &p1, Point &p2);
    vector<Point> trainingPoint;
};

kNN::kNN(int _countNucleotide, int _windowSize, int _countOriWindow, int _countUsualWindow)
{
    countNucleotide = _countNucleotide;
    windowSize = _windowSize;
    countOriWindow = _countOriWindow;
    countUsualWindow = _countUsualWindow;
}

void kNN::ReadOriInfo(string fileNameInfoOri, vector<WindowInfo> &_confirmedOriWindows, vector<WindowInfo> &_usualWindows)
{
    ifstream inOriInfo(fileNameInfoOri);
    string info;
    WindowInfo currentWindowOri, predWindowOri;
    currentWindowOri.numChr = 0;
    predWindowOri.numChr = 0;
    int a = 0;
    while(!inOriInfo.eof())
    {
        predWindowOri = currentWindowOri;
        inOriInfo >> info;
        std::vector<std::string> buffer;
        int startPos = 0;
        int currentPos = 0;
        while((currentPos = info.find(',', startPos)) != -1)
        {
            buffer.push_back(info.substr(startPos + 1, currentPos - startPos - 2));
            startPos = currentPos + 1;
        }
        buffer.push_back(info.substr(startPos + 1, info.size() - startPos - 2));
        currentWindowOri.numChr = atoi(buffer[0].c_str());
        currentWindowOri.start = atoi(buffer[1].c_str());
        currentWindowOri.finish = atoi(buffer[2].c_str());
        if(buffer[buffer.size() - 1] == "Confirmed")
        {
            currentWindowOri.type = SEQ_TYPE::ORI_CONFIRMED;
            _confirmedOriWindows.push_back(currentWindowOri);
            if(currentWindowOri.finish - currentWindowOri.start >= 200)
                a++;
        }

        if(currentWindowOri.numChr == predWindowOri.numChr)
        {
            WindowInfo usualWindow;
            usualWindow.numChr = currentWindowOri.numChr;
            usualWindow.start = predWindowOri.finish;
            usualWindow.finish = currentWindowOri.start;
            usualWindow.type = SEQ_TYPE::USUAL;
            _usualWindows.push_back(usualWindow);
        }
    }

    inOriInfo.close();
}

void kNN::CalcTrainingInfo(vector<WindowInfo> &windows, int countElement)
{
    int numCurrentElem = 0;
    int countAddedElem = 0;
    int step = windows.size() / countElement;
    vector<string> notTrainingSeq;
    while(/*countAddedElem < countElement &&*/ numCurrentElem < windows.size())
    {
        string fileName = ".\\fasta\\chr";
        char *buf = new char[10];
        int numCurrentChr = windows[numCurrentElem].numChr;
        _itoa_s(numCurrentChr, buf, 10, 10);
        fileName += buf;
        fileName += ".fa";
        delete []buf;

        ifstream inCurrentChr(fileName);
        vector<InfoFasta> sequence;
        ReadFasta(inCurrentChr, sequence);
       
        while(numCurrentChr == windows[numCurrentElem].numChr && numCurrentElem < windows.size() && countAddedElem < countElement)
        {
            int widthCurrentWindow = windows[numCurrentElem].finish - windows[numCurrentElem].start;
            if(widthCurrentWindow > windowSize)
            {
                int pos = windows[numCurrentElem].start + rand() % (widthCurrentWindow - windowSize);
                //if(windows[numCurrentOri].finish - pos < windowSize)
                //    pos -= windowSize;
                string seq = sequence[0].genome.substr(pos, windowSize);
                Point point;
                point.freqDinucleotide = CalcFrequency(seq);
                point.type = windows[numCurrentElem].type;
                trainingPoint.push_back(point);
                
                countAddedElem++;
            }
            else
            {
                //string seq = sequence[0].genome.substr(windows[numCurrentElem].start, widthCurrentWindow);
                //notTrainingSeq.push_back(seq);
            }
            numCurrentElem += step;
        }
        if(countAddedElem >= countElement)
        {
            for(; numCurrentElem < windows.size() && numCurrentChr == windows[numCurrentElem].numChr; numCurrentElem++)
            {
                string seq;
                int widthCurrentWindow = windows[numCurrentElem].finish - windows[numCurrentElem].start;
                if(widthCurrentWindow > windowSize)
                {
                    int pos = windows[numCurrentElem].start + rand() % (widthCurrentWindow - windowSize);
                    //if(windows[numCurrentOri].finish - pos < windowSize)
                    //    pos -= windowSize;
                    seq = sequence[0].genome.substr(pos, windowSize);
                }
                else
                {
                    //seq = sequence[0].genome.substr(windows[numCurrentElem].start, widthCurrentWindow);
                }
                notTrainingSeq.push_back(seq);
            }
            //break;
        }
        inCurrentChr.close();
    }
    ofstream out;
    if(windows[0].type == SEQ_TYPE::ORI_CONFIRMED)
    {
        out.open("testOri");
    }
    else
    {
        out.open("testUsual");
    }
    out << notTrainingSeq.size() << endl;
    out << windowSize << endl;
    for(int i = 0; i < notTrainingSeq.size(); i++)
    {
        out << notTrainingSeq[i] << endl;
    }
    out.close();
}

void kNN::Training(vector<string> fileNamesChr, string fileNameInfoOri)
{
    vector<WindowInfo> confirmedOriWindows;
    vector<WindowInfo> usualWindows;
    ReadOriInfo(fileNameInfoOri, confirmedOriWindows, usualWindows);

    CalcTrainingInfo(confirmedOriWindows, countOriWindow);
    CalcTrainingInfo(usualWindows, countUsualWindow);
    trainingPoint;
    /*
    for(int i = 0; i < fileNamesChr.size(); i++)
    {

    }
    */
}

SEQ_TYPE kNN::Classification(string sequence, int k)
{
    SEQ_TYPE result = SEQ_TYPE::USUAL;

    Point currentPoint;
    //currentPoint.type = type;
    currentPoint.freqDinucleotide = CalcFrequency(sequence);

    vector<DistInfo> distances;
    distances.resize(trainingPoint.size());
    for(int i = 0; i < trainingPoint.size(); i++)
    {
        distances[i].type = trainingPoint[i].type;
        distances[i].dist = Distance(currentPoint, trainingPoint[i]);
    }

    sort(distances.begin(), distances.end(), CmpDistInfo);

    int countOriInCluster = 0;
    for(int i = 0; i < k; i++)
    {
        if(distances[i].type == SEQ_TYPE::ORI_CONFIRMED)
        {
            countOriInCluster++;
        }
    }
    if(countOriInCluster > k - countOriInCluster)
        result = SEQ_TYPE::ORI_CONFIRMED;

    return result;
}

void kNN::SaveTrainingInfo(string fileName)
{
    if(trainingPoint.size() != 0)
    {
        ofstream out(fileName);
        out << "<BLOCK_PARAMS>" << endl;
        out << countNucleotide << endl;
        out << countOriWindow << endl;
        out << countUsualWindow << endl;
        for(unordered_map<string, float>::iterator j = trainingPoint[0].freqDinucleotide.begin(); j != trainingPoint[0].freqDinucleotide.end(); j++)
            out << (*j).first << " ";
        out << endl;
        out << "<BLOCK_ORI>" << endl;
        int i;
        for(i = 0; trainingPoint[i].type != SEQ_TYPE::USUAL; i++)
        {
            for(unordered_map<string, float>::iterator j = trainingPoint[i].freqDinucleotide.begin(); j != trainingPoint[i].freqDinucleotide.end(); j++)
            {
                out << (*j).second << " ";
            }
            out << endl;
        }
        out << "<BLOCK_USUAL>" << endl;
        for(; i < trainingPoint.size(); i++)
        {
            for(unordered_map<string, float>::iterator j = trainingPoint[i].freqDinucleotide.begin(); j != trainingPoint[i].freqDinucleotide.end(); j++)
            {
                out << (*j).second << " ";
            }
            out << endl;
        }
        out.close();
    }
}
 
void kNN::LoadTrainingInfo(string fileName)
{
    trainingPoint.clear();
    ifstream in(fileName);

    string nameBlock;
    in >> nameBlock;
    vector<string> allNNucleotide;
    allNNucleotide.resize(16);
    if(nameBlock == "<BLOCK_PARAMS>")
    {
        in >> countNucleotide;
        in >> countOriWindow;
        in >> countUsualWindow;
        trainingPoint.resize(countOriWindow + countUsualWindow);
        for(int i = 0; i < 16; i++)
        {
            string tmp;
            in >> tmp;
            allNNucleotide[i] = tmp;
        }
    }
    in >> nameBlock;
    if(nameBlock == "<BLOCK_ORI>")
    {
        for(int i = 0; i < countOriWindow; i++)
        {
            trainingPoint[i].type = SEQ_TYPE::ORI_CONFIRMED;
            for(int j = 0; j < 16; j++)
            {
                float freq;
                in >> freq;
                trainingPoint[i].freqDinucleotide[allNNucleotide[j]] = freq;
            }
        }
    }
    in >> nameBlock;
    if(nameBlock == "<BLOCK_USUAL>")
    {
        for(int i = countOriWindow; i < countUsualWindow + countOriWindow; i++)
        {
            trainingPoint[i].type = SEQ_TYPE::USUAL;
            for(int j = 0; j < 16; j++)
            {
                float freq;
                in >> freq;
                trainingPoint[i].freqDinucleotide[allNNucleotide[j]] = freq;
            }
        }
    }
}

void kNN::InitPoint(unordered_map<string, float> &point)
{
    string dinucleotide[] = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"};

    for(int i = 0; i < 16; i++)
        point[dinucleotide[i]] = 0.0; 
}

unordered_map<string, float> kNN::CalcFrequency(string &sequence)
{
    unordered_map<string, float> result;
    InitPoint(result);
    
    for(int i = 0; i < sequence.size() - 1; i++)
    {
        result[sequence.substr(i, countNucleotide)]++;
    }

    return result;
}

void kNN::ReadFasta(ifstream &in, vector<InfoFasta> &infoGenome)
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

float kNN::Distance(Point &p1, Point &p2)
{
    float result = 0.0;
    auto i = p1.freqDinucleotide.begin();
    auto j = p2.freqDinucleotide.begin();
    for(; i != p1.freqDinucleotide.end(); i++, j++)
    {
        result += (i->second - j->second) * (i->second - j->second);
    }

    result = sqrt(result);

    return result;
}

int main()
{
    int k;
    kNN a;
    vector<string> asd;
    //a.Training(asd, "oriC");
    //a.SaveTrainingInfo("trainingInfo");
    cout  << "Input K -> ";
    cin >> k;
    a.LoadTrainingInfo("trainingInfo");
    ifstream in("testUsual");
    int countAll;
    int windowsSize;
    in >> countAll >> windowsSize;
    int countCorrect = 0;
    for(int i = 0; i < countAll; i++)
    {
        string seq;
        in >> seq;
        SEQ_TYPE res = a.Classification(seq, k);
        if(res == SEQ_TYPE::USUAL)
        {
            countCorrect++;
        }
    }
    cout << "Usual windows: " << countCorrect << " from " << countAll << endl;
    in.close();

    in.open("testOri");
    in >> countAll >> windowsSize;
    countCorrect = 0;
    for(int i = 0; i < countAll; i++)
    {
        string seq;
        in >> seq;
        SEQ_TYPE res = a.Classification(seq, k);
        if(res == SEQ_TYPE::ORI_CONFIRMED)
        {
            countCorrect++;
        }
    }
    cout << "Ori windows " << countCorrect << " from " << countAll << endl;
    in.close();
    system("pause");
}