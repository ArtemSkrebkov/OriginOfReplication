import mdp
import matplotlib.mlab
from pylab import *
import numpy

def ReadModelFile(fileName, label):
	fileModel = open(fileName, 'r')
	results = []
	for s in fileModel:
		s = s.rstrip('\n')
		vec = s.split(' ')
		values = []
		for i in range(1, len(vec) - 1):
			values.append(float(vec[i].split(':')[1]))
		if(vec[0] == label):
			results.append(values)
	return results

def PlotPca(fileName, label, color):
	data = ReadModelFile(fileName, label)
	x = numpy.array(data) 
	pcanode = mdp.nodes.PCANode(output_dim=2, dtype='float32')
	pcanode.train(x)
	pcanode.stop_training()
	y = pcanode.execute(x)
	dataX = []
	dataY = []
	for coord in y:
		dataX.append(coord[0])
		dataY.append(coord[1])
	plot(dataX, dataY, color)
#ori
PlotPca('model13', '1', 'r')
#usual window
PlotPca('model13', '2', 'g')

show()