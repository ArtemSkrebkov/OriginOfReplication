#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>
#include <limits.h>
#include <time.h>

using namespace std;

enum SEQ_TYPE {ORI_CONFIRMED, ORI_LIKELY, ORI_DUBIOUS, USUAL};

struct Point
{
    SEQ_TYPE type;
    unordered_map<string, float> freqDinucleotide;
};

float Distance(Point &p1, Point &p2)
{
    float result = 0.0;
    auto i = p1.freqDinucleotide.begin();
    auto j = p2.freqDinucleotide.begin();
    for(; i != p1.freqDinucleotide.end(); i++, j++)
    {
        result += (i->second - j->second) * (i->second - j->second);
    }

    //result = sqrt(result);

    return result;
}

void LoadTrainingInfo(string fileName, vector<Point> &trainingPoint)
{
    trainingPoint.clear();
    ifstream in(fileName);

    int countNucleotide;
    int countOriWindow;
    int countUsualWindow;

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

vector<vector<Point>> KMeans(vector<Point> &trainingPoint, int k = 2, int maxCountIteration = 1000, float epsilon = 0.000001)
{
    vector<Point> centerClusters;
	vector<int> curPosInCluster;
	vector<Point> predCenterClusters;
	vector<vector<Point>> clusters;
    predCenterClusters.resize(k);
	clusters.resize(k);
    centerClusters.resize(k);
	curPosInCluster.resize(k);
    srand(time(NULL));
    //for(int i = 0; i < k; i++)
    {
        centerClusters[0] = trainingPoint[rand() % (trainingPoint.size() / 2 + 1)];
        centerClusters[1] = trainingPoint[rand() % (trainingPoint.size() / 2 + 1) + (trainingPoint.size() / 2)];
        //for(auto l = centerClusters[i].freqDinucleotide.begin(); l != centerClusters[i].freqDinucleotide.end(); l++)
        //    (*l).second = rand() % 50;
    }
	bool isEnd = false;
	int countIteration = 0;
	while(!isEnd && maxCountIteration > countIteration)
	{
		predCenterClusters = centerClusters;
		//чистим кластеры и перевыделяем память
		for(int i = 0; i < k; i++)
		{
			curPosInCluster[i] = 0;
		}
		//первые два цикла перебирают все точки изображения
        for(int i = 0; i < trainingPoint.size(); i++)
        {
            float min = INT_MAX;
            int numMinCluster = 0;
            Point point = trainingPoint[i];
            for(int l = 0; l < k; l++)
            {
                float tmp;
                if((tmp = Distance(point, centerClusters[l])) < min)
                {
                    min = tmp;
                    numMinCluster = l;
                }
            }
            if(countIteration == 0)
                clusters[numMinCluster].push_back(point);
            else
            {
                if(curPosInCluster[numMinCluster] >= clusters[numMinCluster].size())
                    clusters[numMinCluster].resize(curPosInCluster[numMinCluster] * 2 + 1);
                clusters[numMinCluster][curPosInCluster[numMinCluster]] = point;
            }
            curPosInCluster[numMinCluster]++;
        }		

		//пересчитываем центры
		for(int i = 0; i < k; i++)
		{
            Point newCenter = trainingPoint[0];
            for(auto i = newCenter.freqDinucleotide.begin(); i != newCenter.freqDinucleotide.end(); i++)
                (*i).second = 0.0;
			if(curPosInCluster[i] != 0)
			{
				for(int j  = 0; j < curPosInCluster[i]; j++)
				{
                    for(auto l = newCenter.freqDinucleotide.begin(); l != newCenter.freqDinucleotide.end(); l++)
                        newCenter.freqDinucleotide[(*l).first] += clusters[i][j].freqDinucleotide[(*l).first];
				}
                for(auto l = newCenter.freqDinucleotide.begin(); l != newCenter.freqDinucleotide.end(); l++)
                    (*l).second /= curPosInCluster[i];
				centerClusters[i] = newCenter;
			}
		}
		//проверяем совпадают ли центры на последней и предпоследней итерациях
        isEnd = true;
        for(int i = 0; i < k; i++)
        {
            if(Distance(centerClusters[i], predCenterClusters[i]) > epsilon)
            {
                isEnd = false;
                break;
            }
        }
		countIteration++;
	}
    cout << "Count of iteration "<< countIteration << endl;
    
    ofstream out("experiment_7", ios::app);
    for(int i = 0; i < k; i++)
    {
        int countOri = 0, countUsual = 0;
        for(int j = 0; j < curPosInCluster[i]; j++)
            if(clusters[i][j].type == SEQ_TYPE::ORI_CONFIRMED)
                countOri++;
            else
                countUsual++;
        out << "Size " << i << " cluster " << curPosInCluster[i] << endl;
        out << "Count Ori window " << countOri << endl;
        out << "Count Usual window " << countUsual << endl;
    }
    out << endl;
    //for(int i = 0; i < k; i++)
    //    cout << "Size " << i << " cluster " << curPosInCluster[i] << endl;
    return clusters;
}


int main()
{
    //srand(time(NULL));
    vector<Point> trainingPoint;
    LoadTrainingInfo("trainingInfo_7", trainingPoint);

    for(int i = 0; i < 50; i++)
        vector<vector<Point>> clusters = KMeans(trainingPoint);

    system("pause");
    return 0;
}