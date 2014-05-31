import glob
import random
import os

def func(a):
    print(a)
    c = [1, 2, 3, 4]
    for b in c:
        print(b)
    g = [b * 10 for b in c]
    print(g)

def DelSymbol(s, deletechars):
    for c in deletechars:
        s = s.replace(c, '')
    return s

def ReadFasta(fileName):
    fileSeq = open(fileName, 'r')
    fileSeq.readline()
    seq = fileSeq.read()    
    fileSeq.close()   
    seq = DelSymbol(seq, '\n ')  
    return seq

def ReadListFasta(fileNames, dir):
    seqs = []
    for fileName in fileNames:
        seqs.append(ReadFasta(dir + fileName))
    return seqs



def ReadInfoOri(fileName, dir):
    fileOri = open(dir + fileName, 'r')
    infoOriDicts = []
    keys = DelSymbol(fileOri.readline(), '\"\n').split(',')
    for s in fileOri:
        infoOriList = s.split(',')
        infoOriDict = {}
        for i in range(0, len(infoOriList)):
            val = DelSymbol(infoOriList[i], '\"\n')
            infoOriDict[keys[i]] = int(val) if val.isdigit() else val
        infoOriDicts.append(infoOriDict)
    return infoOriDicts

def GetInfoStatus(dictOri, status):
    dictStatus = []
    if status != 'Usual':
        for x in dictOri:
            if(x['status'] == status):
                dictStatus.append(x)
    else:
        l = 0
        for x in dictOri:
            r = x['start'] - 1
            elem = {'chr':x['chr'], 'start':l, 'finish':r,'name':'', 'altName':'', 'status':'Usual'} 
            dictStatus.append(elem)
            l = x['finish'] + 1  
    return dictStatus

def GetSampleArticle(dir, windowSize, status, positiveSample, negativeSample):
    infoOri = ReadInfoOri('oriInfo', dir)
    infoOri = GetInfoStatus(infoOri, status) 
    curChr = 1;
    seqs = ReadListFasta(glob.glob(dir + "*.fa"), "")
    for x in infoOri:
        if(curChr != x['chr']):
            curChr = x['chr']
        #if(x['finish'] - x['start'] + 1 >= windowSize):
        s = seqs[curChr - 1][x['start']:x['start'] + windowSize]
        positiveSample.append(s)
        s = seqs[curChr - 1][x['start'] - windowSize:x['start']]
        negativeSample.append(s)
        s = seqs[curChr - 1][x['start'] - 2 * windowSize:x['start'] - windowSize]
        negativeSample.append(s)
        s = seqs[curChr - 1][x['start'] + windowSize:x['start'] + 2 * windowSize]
        negativeSample.append(s)

def GetWindows(dir, windowSize = 250, status = 'Confirmed', countWindowInOri = 0):
    infoOri = ReadInfoOri('oriInfo', dir)
    infoOri = GetInfoStatus(infoOri, status) 
    curChr = 1;
    seqs = ReadListFasta(glob.glob(dir + "*.fa"), "")
    windows = []
    for x in infoOri:
        if(curChr != x['chr']):
            curChr = x['chr']
        if(x['finish'] - x['start'] + 1 >= windowSize):
            s = seqs[curChr - 1][x['start']:x['start'] + windowSize]
            windows.append(s)
    return windows

def GetListNucl(countNucl):
    mono = "ACGT"
    count = 4 ** countNucl
    nucliotide = []
    curNucl = [0 for x in range(0, countNucl)]
    for i in range(0, count):
        curNuclStr = ""
        for k in curNucl:
            curNuclStr += mono[k]
        curNucl[-1] += 1
        for j in range(1, countNucl):
            if(curNucl[-j] >= 4):
                curNucl[-j - 1] += 1
                curNucl[-j] = 0
        nucliotide.append(curNuclStr)

    return nucliotide

def GetDistributionNucl(seq, countNucl = 1):
    nucliotide = GetListNucl(countNucl)
    dictDistrNucl = {}
    for i in range(0, len(nucliotide)):
        dictDistrNucl[nucliotide[i]] = 0.0

    for i in range(0, len(seq) - countNucl + 1):
        dictDistrNucl[seq[i : i + countNucl]] += 1.0
    return dictDistrNucl

def GetNormDistributionNucl(seq, countNucl = 1):
    dictDistrNucl = GetDistributionNucl(seq, countNucl)
    sum = 0.0
    for x in dictDistrNucl.values():
        sum += x
    for key in dictDistrNucl.keys():
        dictDistrNucl[key] /= sum
    return dictDistrNucl

def GetReverseComplement(seq):
    result = ''
    for i in range(1, len(seq) + 1):
        if seq[-i] == 'A':
            result += 'T'
        elif seq[-i] == 'T':
            result += 'A'
        elif seq[-i] == 'G':
            result += 'C'
        elif seq[-i] == 'C':
            result += 'G'
    return result

def GetTableBendability(fileName):
    fileBendability = open(fileName, 'r')
    bendabilityTrinucl = {}
    for s in fileBendability:
        data = s.split(' ')
        bendabilityTrinucl[data[0]] = float(data[1])
        bendabilityTrinucl[GetReverseComplement(data[0])] = float(data[1])
    fileBendability.close()
    return bendabilityTrinucl

def GetListBendability(seq, bendabilityTrinucl, step = 50, sizeWindow = 250):
    bendWindows = []
    for i in range(0, (len(seq) - sizeWindow) // step + 1):
        bendWindows.append(0.0)
        curPos = i * step
        subSeq = seq[curPos: curPos + sizeWindow]
        for j in range(0, len(subSeq) - 3 + 1):
            bendWindows[i] += bendabilityTrinucl[subSeq[j: j + 3]]
    return bendWindows

def GetCleavageIntensity(seq):
    os.system("c:\\strawberry\\perl\\bin\\perl.exe .\\orchid2\\print_orchid1_2.pl -s " + seq + " > out.txt")
    out = open("out.txt", 'r')
    tmp = out.readline()
    tmp = tmp.replace("\t", ' ')
    tmp = tmp.replace("\n", '')
    keys = tmp.split(" ")
    listDict = []
    for s in out:
        s = s.replace("\t", ' ')
        s = s.replace("\n", '')
        values = s.split(" ")
        dicts = {}
        for j in range(0, len(keys)):
            dicts[keys[j]] = float(values[j]) if values[j].isdigit() else values[j]
        listDict.append(dicts)
    sum = 0.0
    for x in listDict:
        sum += float(x['ORChID2'])
    return sum / len(listDict)

def GetListCleavageIntensity(seq, step = 50, sizeWindow = 250):
    cleavageWindows = []
    for i in range(0, (len(seq) - sizeWindow) // step + 1):
        curPos = i * step
        subSeq = seq[curPos: curPos + sizeWindow]
        cleavageWindows.append(GetCleavageIntensity(subSeq))
    return cleavageWindows

def UnionFeature(f1, f2):
    return f1 + f2

def AddFeaturesToFile(listFeatures, label, fileName, typeWrite):
    fileFeature = open(fileName, typeWrite)
    for i in range(0, len(listFeatures)):
        result = str(label) + ' '
        j = 0
        for x in listFeatures[i]:
            j += 1
            result += str(j) + ':' + str(round(x, 8)) + ' '
        result += str(len(listFeatures[i]) + 1) + ':' + '-1'
        result += '\n'
        fileFeature.write(result)

def DrawPlot(seq):
    pass
    #gnuplot_pathname = 'c:\\tmp\\gnuplot\\binary\\pgnuplot.exe' 
    #gnuplot = Popen(gnuplot_pathname,stdin = PIPE,stdout=PIPE,stderr=PIPE).stdin
    #gnuplot() 

#flags - 'dinucl', 'bend', 'cleavInt'
def Experiment(countOri, countUsual, sizeWindow, fileNameOriInfo, dir, flags, fileName):
    table = GetTableBendability('./bendability.txt') 
    allConfirmedOriWindows = GetWindows(dir, sizeWindow, 'Confirmed')
    print(len(allConfirmedOriWindows))
    confirmedOriWindows = random.sample(allConfirmedOriWindows, countOri)
    allUsualWindows = GetWindows(dir, sizeWindow, 'Usual')
    usualWindows = random.sample(allUsualWindows, countUsual)
    featuresOri = []
    #clear file.. great approach..
    f = open(fileName, 'w')
    f.close()
    print(len(confirmedOriWindows))
    i = 0
    for s in confirmedOriWindows:
        k = []
        for flag in flags:
            if flag == 'dinucl':
                l  = GetDistributionNucl(s[:], 2).values()
                for x in l:
                    k.append(x)
            if flag == 'bend':
                l = GetListBendability(s[:], table, 50, 50)
                for x in l:
                    k.append(x)
            if flag == 'cleavInt':
                l = GetListCleavageIntensity(s[:], 50, 50)
                for x in l:
                    k.append(x)
        featuresOri.append(k)
        print(i); i += 1
    AddFeaturesToFile(featuresOri, 1, fileName, 'a')
    print('Features for Ori windows is ready')
    featuresUsual = []
    i = 0
    for s in usualWindows:
        k = []
        for flag in flags:
            if flag == 'dinucl':
                l  = GetDistributionNucl(s[:], 2).values()
                for x in l:
                    k.append(x)
            if flag == 'bend':
                l = GetListBendability(s[:], table, 50, 50)
                for x in l:
                    k.append(x)
            if flag == 'cleavInt':
                l = GetListCleavageIntensity(s[:], 50, 50)
                for x in l:
                    k.append(x)
        featuresUsual.append(k)
        print(i); i += 1
    AddFeaturesToFile(featuresUsual, 2, fileName, 'a')
    print('Features for usual windows is ready')


#flags - 'dinucl', 'bend', 'cleavInt', 'trinucl'
def Experiment2(countOri, countUsual, sizeWindow, fileNameOriInfo, dir, flags, fileName):
    table = GetTableBendability('./bendability.txt') 
    allConfirmedOriWindows = []
    allUsualWindows = []
    GetSampleArticle(dir, sizeWindow, 'Confirmed', allConfirmedOriWindows, allUsualWindows)
    print(len(allConfirmedOriWindows), len(allUsualWindows))
    confirmedOriWindows = random.sample(allConfirmedOriWindows, countOri)
    usualWindows = random.sample(allUsualWindows, countUsual)
    featuresOri = []
    #clear file.. great approach..
    f = open(fileName, 'w')
    f.close()
    print(len(confirmedOriWindows))
    i = 0
    for s in confirmedOriWindows:
        k = []
        for flag in flags:
            if flag == 'dinucl':
                l  = GetDistributionNucl(s[:], 2).values()
                for x in l:
                    k.append(x)
            if flag == 'trinucl':
                l  = GetDistributionNucl(s[:], 3).values()
                for x in l:
                    k.append(x)
            if flag == 'bend':
                l = GetListBendability(s[:], table, 50, 50)
                for x in l:
                    k.append(x)
            if flag == 'cleavInt':
                l = GetListCleavageIntensity(s[:], 50, 50)
                for x in l:
                    k.append(x)
            #normalize
            sum = 0.0
            for x in k:
                sum += x
            if sum != 0:
                for j in range(0, len(k) - 1):
                    k[j] /= sum
        featuresOri.append(k)
        print(i); i += 1
    AddFeaturesToFile(featuresOri, 1, fileName, 'a')
    print('Features for Ori windows is ready')
    featuresUsual = []
    i = 0
    for s in usualWindows:
        k = []
        for flag in flags:
            if flag == 'dinucl':
                l  = GetDistributionNucl(s[:], 2).values()
                for x in l:
                    k.append(x)
            if flag == 'trinucl':
                l  = GetDistributionNucl(s[:], 3).values()
                for x in l:
                    k.append(x)
            if flag == 'bend':
                l = GetListBendability(s[:], table, 50, 50)
                for x in l:
                    k.append(x)
            if flag == 'cleavInt':
                l = GetListCleavageIntensity(s[:], 50, 50)
                for x in l:
                    k.append(x)
            #normalize
            sum = 0.0
            for x in k:
                sum += x
            if(sum != 0):
                for j in range(0, len(k) - 1):
                    k[j] /= sum
        featuresUsual.append(k)
        print(i); i += 1
    AddFeaturesToFile(featuresUsual, 2, fileName, 'a')
    print('Features for usual windows is ready')

#print(GetCleavageIntensity("ACGTAC"))
Experiment2(250, 250, 250, "oriInfo", "./fasta/", ['dinucl','trinucl','bend', 'cleavInt'], 'model13')
#infoOri = ReadInfoOri('oriInfo', './')
#confirmedOriWindows = GetWindows("./fasta/", 250, 'Confirmed')
#usualWindows = GetWindows("./fasta/", 250, 'Usual')
#print(len(confirmedOriWindows))
#print(len(usualWindows))
#bends = []
#bends.append(GetBendability("AAAAA", './bendability.txt', 1, 3))
#bends.append(GetBendability("AAAAA", './bendability.txt', 1, 3))
#AddFeaturesToFile(bends, 1, 'model')
#print(GetNormDistributionNucl("AAAAAAA", 2))