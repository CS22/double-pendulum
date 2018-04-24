# Double Pendulum Simulation
# by Charlie Seymour

# "To improve is to change; to be perfect is to change often."
# - Winston Churchill

import sys, pickle
from PyQt4 import QtGui, QtCore
from math import sin, cos, acos, pi, sqrt, floor

def cutAngle(th):
    # This function ensures that any angle will be kept within the range
    # of -pi --> pi. This is the standard range for representing an angle in
    # radians. This is in place to ensure angles do not grow large enough to
    # slow down the program.
    if th > 0:
        if th < pi:
            return th
        elif th < 2*pi:
            return th-2*pi
        else:
            return cutAngle(((th/pi)%2)*pi)
    else:
        if th > -pi:
            return th
        elif th > -2*pi:
            return th + 2*pi
        else:
            return cutAngle(((th/pi)%2)*pi)

def getPictureName(name, ext):
    # This function will find the next available file name, following the
    # pattern 'name0.ext', 'name1.ext' etc, to ensure the saved images
    # do not overwrite old ones.
        i = 0
        while True:
            try:
                open(name+str(i)+ext)
                i += 1
            except:
                break
        return name+str(i)+ext

class Communicator(QtCore.QObject):
    # This is the object that will send out the simulation's tick

    def __init__(self):

        super().__init__()

        self.pendulumTick = QtCore.QTimer()

class DoublePendulum(QtGui.QWidget):
    # This is the double pendulum widget

    # custom signal to send out attributes to graph each tick
    attributeSignal = QtCore.pyqtSignal(float, float, float, float)

    def __init__(self, L1, L2, M1, M2, g, th1, th2, om1, om2, h):
        # initialises pendulum
        super().__init__()
        self.initAttributes(L1, L2, M1, M2, g, th1, th2, om1, om2, h)
        self.setXY()
        self.calculateEnergy()

        self.setFixedSize(self.sideLength, self.sideLength)

        #this will be the bmp that the trail will be drawn on
        self.canvas = QtGui.QImage(self.sideLength, self.sideLength, QtGui.QImage.Format_RGB32)
        self.canvas.fill(QtGui.QColor(0, 0, 0))

    def initAttributes(self, L1, L2, M1, M2, g, th1, th2, om1, om2, h):
        # Intialises attributes

        # constants
        self.L1  = L1
        self.L2  = L2
        self.M1  = M1
        self.M2  = M2
        self.g   = g

        # This is the pixel scale to draw by - e.g if the point is (0.8, 0.2)
        # mathematically, then the pixel coordinates will be (80, 20) if 
        # dScale is 100
        self.dScale = 100

        self.sideLength = self.dScale*2*(self.L1 + self.L2)

        # variables
        self.x1  = 0.0
        self.x2  = 0.0  # these four will be changed straight away but
        self.y1  = 0.0  # setting them here will also create them
        self.y2  = 0.0

        self.th1 = th1
        self.th2 = th2
        self.om1 = om1
        self.om2 = om2

        # delta t
        self.h = h

    def setXY(self):
        # This sets the coordinates by using trigonometry, and the current
        # angles. the 'sh' variables are altered to give the pixel coordinates
        # for drawing the pendulum heads.

        self.x1 = self.L1*sin(self.th1)
        self.y1 = -self.L1*cos(self.th1)

        self.x2 = self.x1 + self.L2*sin(self.th2)
        self.y2 = self.y1 - self.L2*cos(self.th2)

        shift = self.L1 + self.L2

        self.x1sh = (self.x1+shift)*self.dScale
        self.y1sh = -(self.y1-shift)*self.dScale

        self.x2sh = (self.x2+shift)*self.dScale
        self.y2sh = -(self.y2-shift)*self.dScale

    def calculateEnergy(self):

        # See documentation for derivation. These values will allow the
        # graph to know the maximum value it can have on an axis
        self.ke = 0.5*(self.M1*(self.L1**2)*(self.om1**2) + self.M2*((self.L1**2)*(self.om1**2) + (self.L2**2)*(self.om2**2) + 2*self.L1*self.L2*self.om1*self.om2*cos(self.th1-self.th2)))
        shift = self.L1 + self.L2
        self.gpe = self.g*(self.M1*(self.y1+shift) + self.M2*(self.y2+shift))
        self.energy = self.ke + self.gpe
        self.gpeMin = self.g*self.L2*self.M1
        self.om1Max = self.omegaMax(1)
        self.om2Max = self.omegaMax(2)

        # if there is enough energy, the theoretical maximum angle could be
        # larger than vertical. This doesn't make sense, so this ensures that
        # pi is the the maximum possible angle (180 degrees)
        if self.energy >= (self.g*(2*self.L1*(self.M1+self.M2)+self.M1*self.L2)):
            self.th1Max = pi
        else:
            self.th1Max = acos(round((-self.energy+self.g*(shift*(self.M1 +self.M2)-self.M2*self.L2))/(self.g*self.L1*(self.M1 + self.M2)), 14))

        if self.energy >= (self.L2*self.g*(self.M1 + 2*self.M2)):
            self.th2Max = pi
        else:
            self.th2Max = acos(round((-self.energy+self.g*(self.M1+self.M2)*(shift-self.L1))/(self.M2*self.g*self.L2), 14))

    def omegaMax(self, num):
        # This calculates the maximum value for omega 1 or 2, depending on
        # whether num = 1 or 2
        # See documentation for derivation

        A = 0.5*(self.L1**2)*(self.M1 + self.M2)
        B = 0.5*self.M2*(self.L2**2)
        C = self.L1*self.L2*self.M2
        D = self.energy - self.gpeMin

        if num == 2:
            temp = A
            A = B
            B = temp

        F = C**2 - 4*A*B
        G = 4*A*D

        om2 = -sqrt(round((G*(C**2))/((F**2)-F*(C**2)), 13))
        # note: this is om2 that will give max om1, this is NOT self.om2
        omMax = (-C*om2+sqrt(round((C**2)*(om2**2)-(4*A*B*(om2**2))+4*A*D, 13)))/(2*A)
        return omMax

    def om1d(self, varList):
        # This takes a list of necessary variables and calculates a value
        # for the current angular acceleration of the first arm. 
        # See documentation for derivation.
        th1 = varList[0]
        th2 = varList[1]
        om1 = varList[2]
        om2 = varList[3]
        m1 = self.M1
        m2 = self.M2
        L1 = self.L1
        L2 = self.L2
        g = self.g
        return (-g*(2*m1 + m2)*sin(th1) - m2*g*sin(th1-2*th2)- 2*sin(th1-th2)*m2*((om2**2)*L2 + (om1**2)*L1*cos(th1-th2)))/(L1*(2*m1 + m2 - m2*cos(2*th1-2*th2)))

    def om2d(self, varList):
        # This takes a list of necessary variables and calculates a value
        # for the current angular acceleration of the second arm. 
        # See documentation for derivation.
        th1 = varList[0]
        th2 = varList[1]
        om1 = varList[2]
        om2 = varList[3]
        m1 = self.M1
        m2 = self.M2
        L1 = self.L1
        L2 = self.L2
        g = self.g
        return (2*sin(th1-th2)*((om1**2)*L1*(m1+m2)+g*(m1+m2)*cos(th1)+(om2**2)*L2*m2*cos(th1-th2)))/(L2*(2*m1+m2-m2*cos(2*th1 - 2*th2)))

    def th1d(self, varList):
        # as th1' = om1, this simply return om1.
        return varList[2]

    def th2d(self, varList):
        # as th2' = om2, this simply return om2.
        return varList[3]

    def rungeKutta(self, varList, funcList):
        # See documentation for explanation of Runge Kutta method
        newVarList = []
        h = self.h
        aList = []

        #for func in funcList:
            #aList.append(func(varList))
        aList = self.multiFuncVect(varList, funcList, varList, 0)
        bList = self.multiFuncVect(varList, funcList, aList, h/2)
        cList = self.multiFuncVect(varList, funcList, bList, h/2)
        dList = self.multiFuncVect(varList, funcList, cList, h)

        for var, a, b, c, d in zip(varList, aList, bList, cList, dList):
            newVarList.append(var + (h/6)*(a + 2*b + 2*c + d))

        return newVarList


    def multiFuncVect(self, varList, funcList, letList, h):
        # Part of the Runge Kutta algorithm. See documentation
        newLetList = []
        tempVarList = []

        for var, let in zip(varList, letList):
            tempVarList.append(var + h*let)

        for func in funcList:
            newLetList.append(func(tempVarList))

        return newLetList

    def advance(self):
        # Advances pendulum one step
        # Implements Runge Kutte then assigns the new attributes

        varList = [self.th1, self.th2, self.om1, self.om2]
        funcList = [self.th1d, self.th2d, self.om1d, self.om2d]

        varList = self.rungeKutta(varList, funcList)

        self.th1, self.th2, self.om1, self.om2 = varList
        self.th1 = cutAngle(self.th1)
        self.th2 = cutAngle(self.th2)
        self.setXY()

    def paintEvent(self, event):
        # This describes how the pendulum is drawn to the screen
        headWidth = 16
        headRadius = headWidth // 2
        s = self.dScale*(self.L1 + self.L2)/2 # this is the pixel shift
                                              # mentioned earlier

        # qp is the QPainter object that we will control to do the drawing
        qp = QtGui.QPainter()
        qp.begin(self)

        # draw the second pendulum's point onto the trail image
        self.canvas.setPixel(self.x2sh, self.y2sh, QtGui.qRgb(255, 0, 0))
        qp.drawImage(0, 0, self.canvas)

        # calculate the head colour depending on the mass
        # i.e. more blue = more mass, more white = less mass
        head0colour = QtGui.QColor(255, 255, 255)
        colFrac = 1-self.M1/2
        head1colour = QtGui.QColor(255*colFrac, 255*colFrac, 255)
        colFrac = 1-self.M2/2
        head2colour = QtGui.QColor(255*colFrac, 255*colFrac, 255)

        # calculate the arm colour depending on the angular velocity
        # i.e more green = closer to max speed, more white = closer to still
        # if-else check is to stop error being raised if pendulum starts
        # vertically downwards
        if self.om1Max == 0:
            line1colour = QtGui.QColor(255, 255, 255)
        else:
            colFrac = 1 - abs(self.om1)/self.om1Max 
            line1colour = QtGui.QColor(255*colFrac, 255, 255*colFrac)
        if self.om2Max == 0:
            line2colour = QtGui.QColor(255, 255, 255)
        else:
            colFrac = 1 - abs(self.om2)/self.om2Max
            line2colour = QtGui.QColor(255*colFrac, 255, 255*colFrac)


        # draw the arms and heads
        qp.setBrush(head0colour)
        qp.setPen(head0colour)
        qp.drawEllipse(2*s - headRadius, 2*s - headRadius, headWidth, headWidth)

        qp.setBrush(head1colour)
        qp.setPen(head1colour)
        qp.drawEllipse(self.x1sh - headRadius, self.y1sh - headRadius, headWidth, headWidth)

        qp.setBrush(head2colour)
        qp.setPen(head2colour)
        qp.drawEllipse(self.x2sh - headRadius, self.y2sh - headRadius, headWidth, headWidth)

        qp.setBrush(line1colour)
        qp.setPen(line1colour)
        qp.drawLine(2*s, 2*s, self.x1sh, self.y1sh)

        qp.setBrush(line2colour)
        qp.setPen(line2colour)
        qp.drawLine(self.x1sh, self.y1sh, self.x2sh, self.y2sh)
        qp.end()

    def move(self):
        # advance pendulum, send signal to graph, repaint
        self.advance()
        self.attributeSignal.emit(self.th1, self.th2, self.om1, self.om2)
        self.repaint()

    def step(self):
        # originally implemented for debugging, this function will display
        # the pendulum's current attributes, and also advance by one tick
        self.calculateEnergy()
        print('     x1:', self.x1)
        print('     y1:', self.y1)
        print('     x2:', self.x2)
        print('     y2:', self.y2)
        print('    th1:', self.th1)
        print('    th2:', self.th2)
        print('    om1:', self.om1)
        print('    om2:', self.om2)
        print('     ke:', self.ke)
        print('    gpe:', self.gpe)
        print('om1 max:', self.om1Max)
        print('om2 max:', self.om2Max)
        print('th1 max:', self.th1Max)
        print('th2 max:', self.th2Max)
        print(' gpeMin:', self.gpeMin)
        print(' energy:', self.energy)
        print()
        self.move()

    def export(self):
        # save image of trail
        name = getPictureName('picture', '.bmp')                
        self.canvas.save(name, None, 100)
        
    def saveState(self):
        # save the penudlum's current state to file
        saveFile = open('Pendulum state save.pendulum', 'wb')
        attributes = (self.L1, self.L2, self.M1, self.M2, self.g, self.th1, self.th2, self.om1, self.om2, self.h)
        pickle.dump(attributes, saveFile)
        saveFile.close()

class GraphWidget(QtGui.QWidget):
    # This is the graph widget

    def __init__(self, xAxisVar, yAxisVar, xMinMax, yMinMax, sideLength):
        # intialise the graph
        
        super().__init__()
        self.setFixedSize(sideLength, sideLength)

        # This is the bitmap that will be drawn on
        self.canvas = QtGui.QImage(sideLength, sideLength, QtGui.QImage.Format_RGB32)
        self.canvas.fill(QtGui.QColor(0, 0, 0))

        # get attributes from __init__'s arguments
        self.xAxisVar = xAxisVar
        self.xVal = 0.0
        self.xValAdj = 25
        self.xMin = xMinMax[0]
        self.xMax = xMinMax[1]
        self.yAxisVar = yAxisVar
        self.yVal = 0.0
        self.yValAdj = 25
        self.yMin = yMinMax[0]
        self.yMax = yMinMax[1]
        self.sideLength = sideLength
        # drawing scale (same idea as in the pendulum)
        if (self.xMax-self.xMin) == 0:
            self.xDScale = 0
        else:
            self.xDScale = (sideLength-50)/(self.xMax-self.xMin)

        if (self.yMax-self.yMin) == 0:
            self.yDScale = 0
        else:
            self.yDScale = (sideLength-50)/(self.yMax-self.yMin)
            
    def paintEvent(self, event):

        qp = QtGui.QPainter()
        qp.begin(self)

        #draw points
        qp.drawImage(0, 0, self.canvas)

        #draw axis
        qp.setPen(QtGui.QColor(255, 255, 255))
        topLeft = QtCore.QPoint(25, 25)
        bottomLeft = QtCore.QPoint(25, self.sideLength-25)
        bottomRight = QtCore.QPoint(self.sideLength-25, self.sideLength-25)
        qp.drawLine(topLeft, bottomLeft)
        qp.drawLine(bottomLeft, bottomRight)

        #draw markings
        #write number markings
        for i in range(floor(self.xMin)+1, floor(self.xMax)+1):  #x markings
            x = self.adjustX(i)
            qp.drawLine(x, self.sideLength - 25, x, self.sideLength - 20)
            qr = qp.fontMetrics().boundingRect(str(i))
            qp.drawText(x-(qr.width()/2), self.sideLength - 6, str(i))

        for j in range(floor(self.yMin)+1, floor(self.yMax)+1): #y markings
            y = self.adjustY(j)
            qp.drawLine(25, y, 20, y)
            qr = qp.fontMetrics().boundingRect(str(j))
            qp.drawText((20-qr.width())/2, y+qr.height()/3, str(j))

        #label axis
        qr = qp.fontMetrics().boundingRect(self.xAxisVar)
        qp.drawText(self.sideLength-22, self.sideLength-25 + qr.height()/3, self.xAxisVar)
        qr = qp.fontMetrics().boundingRect(self.yAxisVar)
        qp.drawText(25-qr.width()/2, 23, self.yAxisVar)

        #draw head
        qp.drawRect(self.xValAdj-1, self.yValAdj-4, 3, 9)
        qp.drawRect(self.xValAdj-4, self.yValAdj-1, 9, 3)

        qp.end()

    def adjustPoints(self):
        # The two functions called here adjust the points so that they lie
        # inside the axis and peform the necessary 'pixel shift' performed earlier
        self.xValAdj = self.adjustX(self.xVal)
        self.yValAdj = self.adjustY(self.yVal)

    def adjustX(self, x):
        return (x-self.xMin)*self.xDScale + 25

    def adjustY(self, y):
        return self.sideLength - ((y-self.yMin)*self.yDScale + 25)

    def step(self, th1, th2, om1, om2):
        # this function takes the values from the graph, selects the ones
        # it needs depending on what is on each axis, draws the point onto
        # the bitmap then repaints the graph
        if self.xAxisVar == 'th1':
            self.xVal = th1
        if self.xAxisVar == 'th2':
            self.xVal = th2
        if self.xAxisVar == 'om1':
            self.xVal = om1
        if self.xAxisVar == 'om2':
            self.xVal = om2

        if self.yAxisVar == 'th1':
            self.yVal = th1
        if self.yAxisVar == 'th2':
            self.yVal = th2
        if self.yAxisVar == 'om1':
            self.yVal = om1
        if self.yAxisVar == 'om2':
            self.yVal = om2

        self.adjustPoints()
        self.canvas.setPixel(self.xValAdj, self.yValAdj, QtGui.qRgb(255, 0, 0))
        self.repaint()

    def export(self):
        # save the graph image to file
        name = getPictureName('graph', '.bmp')
        self.canvas.save(name, None, 100)


class DpWindow(QtGui.QWidget):
    # this is the window that will contain all the widgets

    def __init__(self):
        # intialise everything
        super().__init__()
        self.setInitialPendulumAttributes()
        self.initPendulum()
        self.initComms()
        self.setInitialGraphAttributes()
        self.initGraph()
        self.initUI()
        self.initConnections()

        self.move(100, 100)
        self.show()

    def setInitialPendulumAttributes(self, iL1 = 1, iL2 = 1, iM1 = 1, iM2 = 1, ig = 9.81, ith1 = 0, ith2 = 0, iom1 = 0, iom2 = 0, ih = 0.004):
        # if this is called without arguments, it will set default values.
        # otherwise, it can be used to set all of a pendulum's initial
        # attributes.
        self.iL1 = iL1
        self.iL2 = iL2
        self.iM1 = iM1
        self.iM2 = iM2
        self.ig = ig
        self.ith1 = ith1
        self.ith2 = ith2
        self.iom1 = iom1
        self.iom2 = iom2
        self.ih = ih

    def setInitialGraphAttributes(self):
        # default variables for the x and y axes
        self.ixVar = 'om1'
        self.iyVar = 'om2'
        self.setGraphMinMax()

    def setGraphMinMax(self):
        # create a dictionary that stores the min and max variables
        # with the name of the variable on the axis as the key.
        # This is simply so that a long if-elif-else clause is not
        # needed
        minMaxDict = {'om1':(-self.pendulum.om1Max, self.pendulum.om1Max),
                      'om2':(-self.pendulum.om2Max, self.pendulum.om2Max),
                      'th1':(-self.pendulum.th1Max, self.pendulum.th1Max),
                      'th2':(-self.pendulum.th2Max, self.pendulum.th2Max)}
        self.ixMinMax = minMaxDict[self.ixVar]
        self.iyMinMax = minMaxDict[self.iyVar]
        self.isideLength = self.pendulum.sideLength

    def initPendulum(self):
        # create the pendulum object (widget)
        self.pendulum = DoublePendulum(self.iL1, self.iL2, self.iM1, self.iM2, self.ig, self.ith1, self.ith2, self.iom1, self.iom2, self.ih)

    def initComms(self):
        # create the communicator object
        self.comms = Communicator()

    def initGraph(self):
        # create the graph object (widget)
        self.graph = GraphWidget(self.ixVar, self.iyVar, self.ixMinMax, self.iyMinMax, self.isideLength)

    def resetWidgets(self):
        # all the steps needed to re-initialise the widgets if their
        # starting attributes have been changed
        self.removeConnections()
        self.pendulum.deleteLater()
        self.graph.deleteLater()
        self.initPendulum()
        self.initComms()
        self.setGraphMinMax()
        self.initGraph()
        self.initConnections()
        self.addGraphToGrid()
        self.addPendulumToGrid()

    # all of the set (x) from slider functions below follow a similar pattern.
    # they multiply the value of the slider by their maximum value divided
    # by the slider's maximum. This effectively scales their value as the
    # slider scales. Some are slightly different to take into accuont things
    # like negatives, or non zero minimums, however the basic idea stays
    # the same

    def setTh1FromSlider(self, val):
        angle = (2*pi/100)*val
        self.ith1 = angle
        self.th1Label.setText('th1 = ' + str(round(self.ith1, 3)))
        self.resetWidgets()

    def setTh2FromSlider(self, val):
        angle = (2*pi/100)*val
        self.ith2 = angle
        self.th2Label.setText('th2 = ' + str(round(self.ith2, 3)))
        self.resetWidgets()

    def setOm1FromSlider(self, val):
        self.iom1 = (3/50)*val
        self.om1Label.setText('om1 = ' + str(round(self.iom1, 3)))
        self.resetWidgets()

    def setOm2FromSlider(self, val):
        self.iom2 = (3/50)*val
        self.om2Label.setText('om2 = ' + str(round(self.iom2, 3)))
        self.resetWidgets()

    def setGfromSlider(self, val):
        self.ig = (19.62/100)*val
        self.gLabel.setText('g = ' + str(round(self.ig, 3)))
        self.resetWidgets()

    def setHfromSlider(self, val):
        self.ih = 0.001 * val
        self.hLabel.setText('h = ' + str(round(self.ih, 5)))
        self.resetWidgets()

    # these functions that set ratios are similar to the ones above,
    # except they take steps to ensure that the total mass or length
    # is always 2
        
    def setLengthRatioFromSlider(self, ratio):
        """ratio is the fraction of length to be given to L1"""
        ratio /= 100
        self.iL1 = 2*ratio
        self.L1Label.setText('L1 = ' + str(round(self.iL1, 3)))
        self.iL2 = 2*(1-ratio)
        self.L2Label.setText('L2 = ' + str(round(self.iL2, 3)))
        self.resetWidgets()

    def setMassRatioFromSlider(self, ratio):
        """ratio is the fraction of length to be given to M1"""
        ratio /= 100
        self.iM1 = 2*ratio
        self.M1Label.setText('M1 = ' + str(round(self.iM1, 3)))
        self.iM2 = 2*(1-ratio)
        self.M2Label.setText('M2 = ' + str(round(self.iM2, 3)))
        self.resetWidgets()

    # setting the graph's axis usign the radio button's group id
    def setGraphYVar(self, yvarID):
        if yvarID == 1:
            self.iyVar = 'om1'
        if yvarID == 2:
            self.iyVar = 'om2'
        if yvarID == 3:
            self.iyVar = 'th1'
        if yvarID == 4:
            self.iyVar = 'th2'
        self.resetWidgets()

    def setGraphXVar(self, xvarID):
        if xvarID == 1:
            self.ixVar = 'om1'
        if xvarID == 2:
            self.ixVar = 'om2'
        if xvarID == 3:
            self.ixVar = 'th1'
        if xvarID == 4:
            self.ixVar = 'th2'
        self.resetWidgets()

    # load the pendulum save state from file, resets the widgets
    def pendulumLoadState(self):
        saveFile = open('Pendulum state save.pendulum', 'rb')
        attributes = pickle.load(saveFile)
        saveFile.close()
        self.setInitialPendulumAttributes(*attributes)
        self.resetWidgets()

    #reset pendulum and graph to the values they have when the program starts
    def setDefaults(self):
        self.setInitialPendulumAttributes()
        self.resetSlidersRBs()
        self.setInitialGraphAttributes()
        self.resetSlidersRBs()
        self.resetWidgets()

    def showHelp(self):
        # opens a new window that has helpText in it
        self.helpText = """
HELP
----

The program is made up of three main sections.


Firstly, there is the double pendulum itself. This is found in the black box in the bottom right of the screen.

The central circle is the fixed pivot; this won't move. You can think of it as being attached to a wall behind the pendulum. The two lines are the two arms of the pendulum, these have no mass and are free to move.
The arms will be white when they are still, and green when they are turning as fast as they can. They will be a shade of green for any value in-between maximum and 0.

The other two circles are the masses of the double pendulum. These will be a shade of blue depending on how much mass they have, where pale is lighter and dark means heavier. These also act as pivots for the
pendulum to move, similar to the central fixed pivot.


Next there is the graph. This is the black box in the bottom left of the window.

This dynamically draws a graph of the pendulum's current conditions. You can set which variable goes on which axis, that will be covered in the controls section. The axis will be labelled and marked so you can see
the different values for each variable.


Lastly there is the control panel.

Along the top there is a row of buttons. 
The first one is 'Start'. This starts the simulation, which will start both the pendulum and the graph at once. 
Next there is the 'Stop' button, which pauses the simulation. 
The 'Step' button steps the simulation forward one tick.
The 'Save state' button will save the pendulum's current conditions to a file (which will be called 'Pendulum state save.pendulum') in the same directory as the program.
The 'Load state' button allows you to load this file, and resume the simulation from the state it was when you saved.
The 'Save pic' button will save a bitmap file of the current graph and pendulum trail in the same directory as the program. If there is already a file there, then the new one will be numbered. So the graph's images
will be saved as 'graph0.bmp', 'graph1.bmp', 'graph2.bmp' etc. and the pendulum's images will be saved as 'picture0.bmp', 'picture1.bmp' etc.
The 'Reset' button will reset the pendulum to the current initial conditions, that is the conditions that you have set or loaded.
The 'Defaults' button will reset the pendulum to the default starting conditions.
The 'Help!' button will show you this page, but of course you already know that!

The top two sliders control the starting angles. 'th1' (Meaning Theta 1, because theta is the Greek letter used in maths to represent angles, and 1 is because it is for the first arm) is the starting angle for the
first arm. If the slider is far left, it the arm will be vertically downwards. If you push the slider right, it will turn through a full circle, ending vertically downwards again. The same applies for 'th2'.
The next two sliders control the starting angular velocities. These are how quickly the pendulum is already turning at the start. Most of the time these will be left as 0. 'om1' (Meaning Omega 1, because omega is
the Greek letter used in maths to represent angular velocity, and 1 because it is the angular velocity of the first arm) is the starting angular velocity for the first arm. It starts in the middle at 0. When it is
far left, the arm will begin the simulation rotating anticlockwise (left). When it is far right, the arm will begin the simulation rotating clockwise (right). Note the shade of green of each arm. The same applies
for 'om2'.
The third slider down on the left represents the ratio of the mass between the heads. If it is slid to the left, the second head will be heavier and the first head lighter (observe the shades of blue), and if it is
slid to the right, the first head will be heavier and the second head lighter.
The third slider down on the right represents the ratio of the length between the heads. If it is slid to the left, the second arm will be longer and the first arm shorter. If it is slid to the right, the first arm
will be longer and the second shorter.
The fourth slider down on the right sets the gravitational field strength. On earth this is approximately 9.81, so that is the starting value. If it is slid to the right, the value increases, meaning everything is
heavier. As it is slid to the left, it becomes weaker, meaning the pendulum acts lighter and lighter. If it is all the way to the left, it is 0, meaning nothing is pulling the pendulum downwards at all.
The bottom right slider controls the 'time step'. This is to do with how the simulation works. The default value is 'h = 0.004', which means that the simulation will work out where to move the pendulum to every
0.004 seconds. This can be thought of as the accuracy of simulation; a high value will result in a very inaccurate (but easily processed) simulation, whereas a very low value will result in a very accurate (but
CPU intensive) simulation. Most computers will start to slow down too much below 0.004, so only change this if the simulation is moving too slowly.

The radio buttons above the graph control which variable goes on which axis; you have the choice of Angle 1, Angle 2, Angular velocity 1 and Angular velocity 2 for each axis."""

        self.helpWindow = QtGui.QWidget()
        self.helpWindow.setWindowTitle("Help")
        self.helpWindow.setWindowIcon(QtGui.QIcon("helpIcon.png"))
        self.helpLabel = QtGui.QLabel(self.helpText)
        self.helpLabel.setParent(self.helpWindow)
        self.helpWindow.setFixedSize(self.helpLabel.sizeHint())
        self.helpWindow.show()

    def resetSlidersRBs(self):
        # sets all sliders and radio buttons to their default values
        self.th1Slider.setValue(0)
        self.th2Slider.setValue(0)
        self.om1Slider.setValue(0)
        self.om2Slider.setValue(0)
        self.massSlider.setValue(50)
        self.lengthSlider.setValue(50)
        self.gSlider.setValue(50)
        self.hSlider.setValue(4)
        self.om1RBX.click()
        self.om2RBY.click()

    def initUI(self):
        # this initialises the UI. It creates, labels and resizes all the
        # buttons, creates and labels all the sliders, creates all the
        # radio buttons, then puts all of them into layouts that will
        # arrange their position on the screen
        self.setFixedSize(self.sizeHint())
        self.setWindowTitle("Double Pendulum Simulation")
        self.setWindowIcon(QtGui.QIcon("icon.png"))

        # Buttons
        self.startButton = QtGui.QPushButton('Start')
        self.startButton.setFixedSize(50, 50)
        self.stopButton  = QtGui.QPushButton('Stop')
        self.stopButton.setFixedSize(50, 50)
        self.stepButton = QtGui.QPushButton('Step')
        self.stepButton.setFixedSize(50, 50)
        self.stateSaveButton = QtGui.QPushButton('Save\nstate')
        self.stateSaveButton.setFixedSize(50, 50)
        self.stateLoadButton = QtGui.QPushButton('Load\nstate')
        self.stateLoadButton.setFixedSize(50, 50)
        self.pictureSaveButton = QtGui.QPushButton('Save\nimage')
        self.pictureSaveButton.setFixedSize(50, 50)
        self.resetButton = QtGui.QPushButton('Reset')
        self.resetButton.setFixedSize(50, 50)
        self.defaultsButton = QtGui.QPushButton('Defaults')
        self.defaultsButton.setFixedSize(50, 50)
        self.helpButton = QtGui.QPushButton('Help!')
        self.helpButton.setFixedSize(50, 50)
        self.helpButton.setStyleSheet('QPushButton {color : red}')

        # Sliders
        self.th1Slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.th1Slider.setRange(0, 100)
        self.th2Slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.th2Slider.setRange(0, 100)
        self.om1Slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.om1Slider.setRange(-50, 50)
        self.om2Slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.om2Slider.setRange(-50, 50)
        self.massSlider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.massSlider.setRange(30, 70)
        self.massSlider.setValue(50)
        self.lengthSlider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.lengthSlider.setRange(30, 70)
        self.lengthSlider.setValue(50)
        self.gSlider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.gSlider.setRange(0, 100)
        self.gSlider.setValue(50)
        self.hSlider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.hSlider.setRange(1, 100)
        self.hSlider.setValue(4)

        # Labels
        self.th1Label = QtGui.QLabel('th1 = ' + str(round(self.ith1,3)))
        self.th2Label = QtGui.QLabel('th2 = ' + str(round(self.ith2,3)))
        self.om1Label = QtGui.QLabel('om1 = ' + str(round(self.iom1,3)))
        self.om2Label = QtGui.QLabel('om2 = ' + str(round(self.iom2,3)))
        self.M1Label = QtGui.QLabel('M1 = ' + str(round(self.iM1,3)))
        self.M2Label = QtGui.QLabel('M2 = ' + str(round(self.iM2,3)))
        self.L1Label = QtGui.QLabel('L1 = ' + str(round(self.iL1,3)))
        self.L2Label = QtGui.QLabel('L2 = ' + str(round(self.iL2,3)))
        self.gLabel = QtGui.QLabel('g = ' + str(round(self.ig,3)))
        self.hLabel = QtGui.QLabel('h = ' + str(round(self.ih,5)))

        #   Radio Buttons
        self.RBXGroup = QtGui.QButtonGroup()
        self.RBYGroup = QtGui.QButtonGroup()

        # X axis radio buttons
        self.om1RBX = QtGui.QRadioButton('om1')
        self.RBXGroup.addButton(self.om1RBX, 1)
        self.om2RBX = QtGui.QRadioButton('om2')
        self.RBXGroup.addButton(self.om2RBX, 2)
        self.th1RBX = QtGui.QRadioButton('th1')
        self.RBXGroup.addButton(self.th1RBX, 3)
        self.th2RBX = QtGui.QRadioButton('th2')
        self.RBXGroup.addButton(self.th2RBX, 4)

        # X axis radio buttons
        self.om1RBY = QtGui.QRadioButton('om1')
        self.RBYGroup.addButton(self.om1RBY, 1)
        self.om2RBY = QtGui.QRadioButton('om2')
        self.RBYGroup.addButton(self.om2RBY, 2)
        self.th1RBY = QtGui.QRadioButton('th1')
        self.RBYGroup.addButton(self.th1RBY, 3)
        self.th2RBY = QtGui.QRadioButton('th2')
        self.RBYGroup.addButton(self.th2RBY, 4)

        self.om1RBX.toggle()
        self.om2RBY.toggle()
        

        # Radio button boxes
        self.RBXbox = QtGui.QHBoxLayout()
        self.RBXbox.addWidget(QtGui.QLabel("X axis:"))
        self.RBXbox.addWidget(self.om1RBX)
        self.RBXbox.addWidget(self.om2RBX)
        self.RBXbox.addWidget(self.th1RBX)
        self.RBXbox.addWidget(self.th2RBX)

        self.RBYbox = QtGui.QHBoxLayout()
        self.RBYbox.addWidget(QtGui.QLabel("Y axis:"))             
        self.RBYbox.addWidget(self.om1RBY)
        self.RBYbox.addWidget(self.om2RBY)
        self.RBYbox.addWidget(self.th1RBY)
        self.RBYbox.addWidget(self.th2RBY)
        

        #   Slider grid
        # Th sliders
        self.sliderRBGrid = QtGui.QGridLayout()
        self.sliderRBGrid.addWidget(self.th1Slider, 0, 0)
        self.sliderRBGrid.addWidget(self.th1Label, 0, 1)
        self.sliderRBGrid.addWidget(self.th2Slider, 0, 3)
        self.sliderRBGrid.addWidget(self.th2Label, 0, 4)

        # Om sliders
        self.sliderRBGrid.addWidget(self.om1Slider, 1, 0)
        self.sliderRBGrid.addWidget(self.om1Label, 1, 1)
        self.sliderRBGrid.addWidget(self.om2Slider, 1, 3)
        self.sliderRBGrid.addWidget(self.om2Label, 1, 4)

        # M/L sliders
        self.sliderRBGrid.addWidget(self.massSlider, 2, 0)
        self.sliderRBGrid.addWidget(self.M1Label, 2, 1)
        self.sliderRBGrid.addWidget(self.M2Label, 2, 2)
        self.sliderRBGrid.addWidget(self.lengthSlider, 2, 3)
        self.sliderRBGrid.addWidget(self.L1Label, 2, 4)
        self.sliderRBGrid.addWidget(self.L2Label, 2, 5)

        # g/h sliders
        self.sliderRBGrid.addWidget(self.gSlider, 3, 3)
        self.sliderRBGrid.addWidget(self.gLabel, 3, 4)
        self.sliderRBGrid.addWidget(self.hSlider, 4, 3)
        self.sliderRBGrid.addWidget(self.hLabel, 4, 4)

        # Radio buttons
        self.sliderRBGrid.addLayout(self.RBXbox, 3, 0)
        self.sliderRBGrid.addLayout(self.RBYbox, 4, 0)

        # Button box
        self.buttonBox = QtGui.QHBoxLayout()
        self.buttonBox.addWidget(self.startButton)
        self.buttonBox.addWidget(self.stopButton)
        self.buttonBox.addWidget(self.stepButton)
        self.buttonBox.addWidget(self.stateSaveButton)
        self.buttonBox.addWidget(self.stateLoadButton)
        self.buttonBox.addWidget(self.pictureSaveButton)
        self.buttonBox.addWidget(self.resetButton)
        self.buttonBox.addWidget(self.defaultsButton)
        self.buttonBox.addWidget(self.helpButton)


        # Pendulum and graph grid
        self.penGraGrid = QtGui.QGridLayout()

        # Set layouts (overall grid)
        self.grid = QtGui.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addLayout(self.buttonBox, 0, 0)
        self.grid.addLayout(self.sliderRBGrid, 1, 0)
        self.addGraphToGrid()
        self.addPendulumToGrid()
        self.initActions()

    def addGraphToGrid(self):
        self.penGraGrid.addWidget(self.graph, 0, 0)
        self.grid.addLayout(self.penGraGrid, 3, 0)

    def addPendulumToGrid(self):
        self.penGraGrid.addWidget(self.pendulum, 0, 1)
        self.grid.addLayout(self.penGraGrid, 3, 0)

    def initConnections(self):
        # This function connects all the signals to the slots that need
        # them. All the buttons are connected to their desired functions
        # when clicked, all the sliders are connected to the desired
        # functions when they are chagned.
        self.comms.pendulumTick.setInterval(self.pendulum.h*1000)
        
        self.comms.pendulumTick.timeout.connect(self.pendulum.move)
        self.pendulum.attributeSignal.connect(self.graph.step)

        self.th1Slider.valueChanged.connect(self.setTh1FromSlider)
        self.th2Slider.valueChanged.connect(self.setTh2FromSlider)
        self.om1Slider.valueChanged.connect(self.setOm1FromSlider)
        self.om2Slider.valueChanged.connect(self.setOm2FromSlider)
        self.massSlider.valueChanged.connect(self.setMassRatioFromSlider)
        self.lengthSlider.valueChanged.connect(self.setLengthRatioFromSlider)
        self.gSlider.valueChanged.connect(self.setGfromSlider)
        self.hSlider.valueChanged.connect(self.setHfromSlider)

        self.RBXGroup.buttonClicked[int].connect(self.setGraphXVar)
        self.RBYGroup.buttonClicked[int].connect(self.setGraphYVar)
        
        self.startButton.clicked.connect(self.comms.pendulumTick.start)
        self.stopButton.clicked.connect(self.comms.pendulumTick.stop)
        self.stepButton.clicked.connect(self.pendulum.step)
        self.stateSaveButton.clicked.connect(self.pendulum.saveState)
        self.stateLoadButton.clicked.connect(self.pendulumLoadState)
        self.pictureSaveButton.clicked.connect(self.graph.export)
        self.pictureSaveButton.clicked.connect(self.pendulum.export)
        self.resetButton.clicked.connect(self.resetWidgets)
        self.defaultsButton.clicked.connect(self.setDefaults)
        self.helpButton.clicked.connect(self.showHelp)

    def removeConnections(self):
        # This disconnects all the connections made in initConnections, so that
        # when the widgits are recreated, there aren't stray connections
        # causing unwanted results.
        self.pendulum.attributeSignal.disconnect(self.graph.step)
        self.comms.pendulumTick.timeout.disconnect(self.pendulum.move)

        self.th1Slider.valueChanged.disconnect(self.setTh1FromSlider)
        self.th2Slider.valueChanged.disconnect(self.setTh2FromSlider)
        self.om1Slider.valueChanged.disconnect(self.setOm1FromSlider)
        self.om2Slider.valueChanged.disconnect(self.setOm2FromSlider)
        self.massSlider.valueChanged.disconnect(self.setMassRatioFromSlider)
        self.lengthSlider.valueChanged.disconnect(self.setLengthRatioFromSlider)
        self.gSlider.valueChanged.disconnect(self.setGfromSlider)
        self.hSlider.valueChanged.disconnect(self.setHfromSlider)

        self.RBXGroup.buttonClicked[int].disconnect(self.setGraphXVar)
        self.RBYGroup.buttonClicked[int].disconnect(self.setGraphYVar)
            
        self.startButton.clicked.disconnect(self.comms.pendulumTick.start)
        self.stopButton.clicked.disconnect(self.comms.pendulumTick.stop)
        self.stepButton.clicked.disconnect(self.pendulum.step)
        self.stateSaveButton.clicked.disconnect(self.pendulum.saveState)
        self.stateLoadButton.clicked.disconnect(self.pendulumLoadState)
        self.pictureSaveButton.clicked.disconnect(self.graph.export)
        self.pictureSaveButton.clicked.disconnect(self.pendulum.export)
        self.resetButton.clicked.disconnect(self.resetWidgets)
        self.defaultsButton.clicked.disconnect(self.setDefaults)
        self.helpButton.clicked.disconnect(self.showHelp)
        
    def initActions(self):
        self.exitAction = QtGui.QAction('&Exit', self)
        self.exitAction.setShortcut('Ctrl+Q')
        self.exitAction.triggered.connect(QtGui.qApp.quit)
        self.addAction(self.exitAction)

# This creates the application object, creates the main window, then
# executes the application
def main():
    app = QtGui.QApplication(sys.argv)
    win = DpWindow()
    app.exec_()

# This checks to make sure the program has been started rather than imported
if __name__ == "__main__":
    main()
