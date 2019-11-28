from scipy import signal
import control as cl
import numpy as np
import matplotlib.pyplot as plt
import random
from Approximation import *

#def convertValueListToPointList(valueList,type):
#    for index, valuePair in enumerate(valueList):
#        for value in valueList:
#            if type == "zero": point=Zero(value, index)



class Point(object):
    stage=None
    pair=None
    def __init__(self, value, color):
        self.color=color
        self.value=value
        self.im=np.imag(value)
        self.re=np.real(value)
        self.mod=np.sqrt((self.im)**2+(self.re)**2)
        if self.re != 0:
            self.pha = np.arctan(self.im/self.re)
        elif self.im > 0:
            self.pha = 3.142/2
        else:
            self.pha = -3.142/2
        if self.re < 0:
            self.pha = 3.1415 - self.pha
    

class PointPair(list):
    def __init__(self, approximation,points,type,name,color=None):
        if isinstance(points[0],Point):
            self.extend(points)
        else:
            if color is None: color="r"
            if type == "zero":
                self.extend([Zero(points[0],color),Zero(points[1],color)])
            elif type == "pole":
                self.extend([Pole(points[0],color),Pole(points[1],color)])
            else: raise ValueError("Invalid type value: "+str(type))
        self.color=self[0].color
        self.type=type
        self.name=name
        self.groupBox=None
        self.stage=None
        self.approximation=approximation
        self.radioButtonLayouts=[]
        for point in points: point.pair=self

    def assignToStage(self, stage):
        for point in self:
            if point.stage is not None: point.stage.removePair(self)
            point.stage=stage
        self.stage=stage
        #stage.addPoints(self)

    def clearStage(self):
        self.stage=None
        for point in self:
            point.stage=None
    
    def isZero(self): 
        if self.type=="zero": return True
        else: return False

    def isPole(self): 
        if self.type=="pole": return True
        else: return False

class Zero(Point):
    def isZero(self): 
        return True
    def isPole(self): 
        return False

class Pole(Point):
    def isPole(self): 
        return True
    def isZero(self): 
        return False

class Etapa:
    def __init__(self, figArray, name):
        """
        Initializes stage.
        figArray=[figure, figureCanvas]
        """
        if len(figArray) == 2:
            self.setFigure(figArray[0])
            self.setFigureCanvas(figArray[1])
        else: raise ValueError("Stage init: invalid figArray")
        self.ax=None
        self.zeroList=[]
        self.poleList=[]
        self.assignedPolePair = None
        self.assignedZeroPair = None
        self.name=name
        self.canAddPole=True
        self.canAddZero=True

        ## @var tf
        #  The transfer function of the stage
        self.tf = None

        ## @var cell
        #  Represents the cell for the stage
        #  It's generated based on its transfer function.
        self.cell = None
        

    def compute_tf(self):
        pass
   
    def addPair(self, pointPair):
        if pointPair.isPole(): 
            if self.canAddPole:
                self.poleList.append(pointPair)
                self.assignedPolePair = pointPair
        elif pointPair.isZero(): 
            if self.canAddZero:
                self.zeroList.append(pointPair)
                self.assignedZeroPair = pointPair
        else : 
            return
        if len(self.poleList) > 0:
            self.canAddPole = False
        if len(self.zeroList) > 0:
            self.canAddZero = False
        pointPair.assignToStage(self)
        

    def removePair(self, pointPair):
        if pointPair.isPole(): 
            self.poleList.remove(pointPair)
            self.assignedPolePair = None

        elif pointPair.isZero(): 
            self.zeroList.remove(pointPair)
            self.assignedZeroPair = None
            
        else :
           return

        if len(self.poleList) == 0:
            self.canAddPole = True
        if len(self.zeroList) == 0:
            self.canAddZero = True
        pointPair.clearStage()


    def getTransferFunction(self):
        return self.num, self.den
 
    def createTransferFunction(self):
        tempZero=[]
        tempPole=[]

        for zeroPair in self.zeroList:
            for zero in zeroPair:
                tempZero.append(zero.value)
        for polePair in self.poleList:
            for pole in polePair:
                tempPole.append(pole.value)

        self.num, self.den, k = signal.zpk2tf(tempZero,tempPole,self.k)
    
    def setFigure(self, figure):
        self.figure=figure

    def setFigureCanvas(self, figurecanvas):
        self.figurecanvas=figurecanvas

    def plotStage(self):
        if self.figure is not None:
            if self.ax is not None: self.figure.delaxes(self.ax) 
            
            self.ax = self.figure.add_subplot(111, projection='polar')
            polarPlot(self.zeroList, self.poleList, self.ax)
            self.figurecanvas.draw()
            
        else: print("Error: Set figure first")
   
def polarPlot(zeros, poles, ax):
    for zeroPair in zeros:
        pha = []
        mod = []
        for zero in zeroPair: 
                pha.append(zero.pha)
                mod.append(zero.mod)
        #ax.plot(pha,mod,'o',color=zero.color,label=zeroPair.name)
        ax.plot(pha,mod,'o',label=zeroPair.name)

    for polePair in poles: 
        pha = []
        mod = []
        for pole in polePair: 
            pha.append(pole.pha)
            mod.append(pole.mod)
        #ax.plot(pha,mod,'x',color=pole.color,label=polePair.name)
        ax.plot(pha,mod,'x',label=polePair.name)
    ax.legend()
    ax.set_ylim(auto=True)
   


# num = [0.8912509381337455,0.0,6832402540.638172,0.0,7.24572875462745e+18]
# den = [1.0,337194.6633703208,92412209932.89781,8400336197283562.0,1.2883166302894993e+21]

# z,p,k = getPolesAndZerosFromTF(num,den)
# plotPolesAndZeros(z,p)

# hola = Etapa(z[0],p[1],k)
# num,den = hola.getTransferFunction()
# z,p,k = getPolesAndZerosFromTF(num,den)
# plotPolesAndZeros(z,p)
