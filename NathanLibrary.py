# Copyright 2021: Nathan Barrett, All Rights Reserved                          #
# NathanLibrary.py is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
# Version: 2021.1.1 (01/16/2021)                                               #
#                                                                              #
######################### A NOTE TO MY FELLOW STUDENTS #########################
# I've compiled all of the most usefull and/or important code that I've        #
# written over the entire course of my time in the Chemical Engineering        #
# program.  While this code has saved me a lot of time on the mundane parts    #
# of my learning, it was WRITING the code that taught me the most.             #
#                                                                              #
# For that reason I HIGHLY encourage you not to just copy and paste my code.   #
# I've given it out as 1) a reference with which you can quickly and easily    #
# check your answers and 2) an example of the verious approaches you can use   #
# to develop your own code.  Take my ideas, layouts, or procedures and         #
# re-create them using your own words, so to speak.  That's going to be        #
# significantly more beneficial to you in the long run than merely copying     #
# and pasting.                                                                 #
#                                                                              #
# Good luck, and happy coding!                                                 #


#----------------------------- IMPORTED LIBRARIES -----------------------------#
import pint, colored, sys, warnings, copy
u = pint.UnitRegistry()
import numpy as np
from scipy.optimize import fsolve, ridder, curve_fit,minimize
import scipy.special as bessel
from scipy import interpolate
from scipy.integrate import quad, simps, trapz, romb, solve_ivp
from scipy.interpolate import interp1d
from scipy.stats import mode
from IPython.display import Image, display, Markdown
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
sp.init_printing()
import sympy.printing.latex as latex
import random as rand
from fractions import Fraction
from iapws import IAPWS97
import xlwings as xw
from inspect import signature
import periodictable as pTable


#------------------------------ GENERAL CONSTANTS -----------------------------#
# Ideal gas constant
R    = 8.31446261815324 * u.J / u.K / u.mol

# Avogadro's number
An   = 6.02214076E23 / u.mol

# Conventional freefall acceleration of gravity
g    = 9.80665  * u.m / (u.sec)**2

#Speed of light
c    = 299792458 * u.m / (u.sec)

#Faraday Ccnstant (the charge of one electron can be found in pint: u.e)
F    = u.e * An

#Definition of "dimensionless" unit
u.define('diml = 1 * sec / sec')

#Definition of "pound mole" unit
u.define('lbmol = 453.59237 * mol')

#A string to mark a section of code that needs to be revised
cb = ""
for i in range(1000):
    cb += "!"
cb += "COME BACK"

#------------------- Non-zero based temperature conversions -------------------#
def degC(T):
    """
    Converts a float or int object representing the temperature in degrees Celcius to a Pint object representing the same temerature

    T       : A float or int representing the temperature in degrees Celcius
    returns : A Pint varible with the converted temperature (returned with associated units of Kelvin)
    """
    if not isinstance(T,int) and not isinstance(T,float):
        raise Exception("degF() can only take an int or float, not a \"" + str(type(T)) + "\"")

    return (T + 273.15) * u.K
def degF(T):
    """
    Converts a float or int object representing the temperature in degrees Fahrenheit to a Pint object representing the same temerature

    T       : A float or int representing the temperature in degrees Fahrenheit
    returns : A Pint varible with the converted temperature (returned with associated units of Kelvin)
    """
    if not isinstance(T,int) and not isinstance(T,float):
        raise Exception("degF() can only take an int or float, not a \"" + str(type(T)) + "\"")

    C = (T - 32) * (5/9)
    return (C + 273.15) * u.K

def getC(T):
    """
    Converts a Pint object representing a temperature to a float object representing the temperature in degrees Celcius

    T       : A Pint variable with a dimensionality of "Temerature"
    returns : A float representing the temperature in degrees Celcius
    """
    return T.to(u.K).magnitude - 273.15
def getF(T):
    """
    Converts a Pint object representing a temperature to a float object representing the temperature in degrees Fahrenheit

    T       : A Pint variable with a dimensionality of "Temerature"
    returns : A float representing the temperature in degrees Fahrenheit
    """
    C = getC(T)
    return C / (5/9) + 32


#------------------------------ General Functions -----------------------------#
def avg(List):
    """
    Takes the arithmetic mean of a list of Pint or other numeric data type objects.

    List    : A 1-D itterable data type containing the numbers you'd like to take the average of
    returns : The arithmetic mean of "List"
    """
    if len(List) == 0:
        return None

    sum = List[0] * 0 #Zero but with the right units
    for e in List:
        sum += e
    return sum / len(List)

def min(List):
    """
    Find the minimum valued element of a list of Pint or other numeric data type objects.

    List    : A 1-D itterable data type containing the numbers you'd like to find the min of
    returns : The minimum valued element of "List"
    """
    if len(List) == 0:
        return None

    minVal = List[0]

    for e in List:
        if e < minVal:
            minVal = e
    return minVal

def max(List):
    """
    Find the maximum valued element of a list of Pint or other numeric data type objects.

    List    : A 1-D itterable data type containing the numbers you'd like to find the max of
    returns : The maximum valued element of "List"
    """
    if len(List) == 0:
        return None

    maxVal = List[0]

    for e in List:
        if e > maxVal:
            maxVal = e
    return maxVal

def printSN(var):
    """
    Prints the given variable (Pint or other numeric data type object) in scientific notation

    var     : The variable you'd like to print
    returns : None
    """
    try:
        value = var.magnitude
        unit = var.units
        print("{:e}".format(value),unit)
    except:
        print("{:e}".format(var))

def numformat(var,sigFigs=3):
    """
    Formats the given variable (Pint or other numeric data type object) to match the specified number if significant figures

    var     : The variable you'd like to format
    sigFigs : The number of significant figures you'd like your formatted variable to have (must be an int)(defaults to 3)
    returns : The a string representing the variable you passed in but formatted to the specified number of significant figures
    """
    raise Exception("THIS CODE IS UNDER DEVELOPMENT!")
    #HAVING BUGS WITH NUMBERS LIKE -2.138

    #First, try to separate the unit and magnitude from the variable.
    try:
        unit = var.units
        var = var.magnitude
    except:
        unit = "NONE"

    #Now convert the variable to a string
    varStr = str(var)

    negative = False

    #Detect / eliminate the negative sign from the variable
    if varStr[0] == '-':
        negative = True
        varStr = varStr[1:]

    firstCharBeforeDecimal = False
    decimalFound = False

    #Detect if the first significant digit occurs before the decimal place.
    for char in varStr:
        if char not in {'.','0'}:
            if not decimalFound:
                firstCharBeforeDecimal = True
                break
        if char == '.':
            firstCharBeforeDecimal = False
            break

    #Place the first digit in the "ones" place
    numPlacesMoved = 0

    if firstCharBeforeDecimal:
        while True:
            varStr = str(var)
            if varStr[1] != '.':
                var /= 10
                numPlacesMoved -= 1
            else:
                break
    else:
        while True:
            varStr = str(var)
            if varStr[0] != '0':
                break
            else:
                var *= 10
                numPlacesMoved += 1

    #Now that the first significant digit is in the ones place, truncate the variable to the desired number of sig figs.
    var = round(var,sigFigs - 1)

    varStr = str(var)

    sigFigsLeft = sigFigs
    newVarStr = ""

    #Sometimes the round() function drops zeros that should be significant.  So we'll add them back in.
    for char in varStr:
        if sigFigsLeft == 0:
            break

        if char == '.':
            newVarStr += char
        else:
            newVarStr += char
            sigFigsLeft -= 1

    for i in range(sigFigsLeft):
        newVarStr += '0'

    #Now take the decimal out of the number so we can reposition it to the correct place.
    newVarStr = newVarStr[0] + newVarStr[2:]

    decimalIndex = 1


    decimalIndex -= numPlacesMoved

    #Insert the decimal and any nececary leading or trailing zeros back into the number in the appropriate place.
    if decimalIndex < 0:
        zeroStr = ""
        for i in range(-1 * decimalIndex + 1):
            zeroStr += '0'

        newVarStr = zeroStr + newVarStr

        newVarStr = newVarStr[0] + '.' + newVarStr[1:]


    elif decimalIndex > len(newVarStr):
        zeroStr = ""
        for i in range(decimalIndex - len(newVarStr)):
            zeroStr += '0'
        newVarStr += zeroStr

    else:
        newVarStr = newVarStr[:decimalIndex] + "." + newVarStr[decimalIndex:]

    #Add back in the negative sign if needed.
    if negative:
        newVarStr =  '-' + newVarStr

    #Add back in the unit if needed.
    if unit != "NONE":
        newVarStr += " " + str(unit)

    return newVarStr

def crossProduct(V1,V2):
    """
    Calculates the cross product of two 3-D vectors (itterable data types with lengths of 3 and containig only Pint or other numeric data type objects)

    V1      : The left-hand 3-D vector
    V2      : The right-hand 3-D vector
    returns : The cross product of V1 and V2 (V1 X V2)
    """
    a = V1[1] * V2[2] - V1[2] * V2[1]
    b = V2[0] * V1[2] - V2[2] * V1[0]
    c = V1[0] * V2[1] - V1[1] * V2[0]
    return [a,b,c]

def solveLinSet(coefs, consts):
    """
    Solves a linear set of equations given in matrix form.
    coefs * answer = consts

    coefs  : n X n matrix         (containig only Pint or other numeric data type objects)
    consts : 1D list of length n. (containig only Pint or other numeric data type objects)

    return : the "aswer" matrix of the linear equation given above ("answer" = coefs**-1 * consts)

    Note : any non-zero term with no unit attatched will be assumed dimmensionless.
    Note : All returned answers will be Pint objects with the appropriate SI units.
    """

    #Put everything in SI Units:
    for i in range(len(coefs)):
        for j in range(len(coefs[i])):
            #If it's zero, leave it, we'll deal with it later
            if coefs[i][j] == 0 and type(coefs[i][j]) != type(1 * u.diml):
                continue
            else:
                try:
                    coefs[i][j] = coefs[i][j].to_base_units()
                except:
                    coefs[i][j] *= u.diml
    for i in range(len(consts)):
        #Again, if it's zero, we'll deal with it later.
        if consts[i] == 0 and type(consts[i]) != type(1 * u.diml):
            continue
        else:
            try:
                consts[i] = consts[i].to_base_units()
            except:
                consts[i] *= u.diml

    #Run the magnitudes through numpy's matrix system (to get answer and also detect invalidity)
    ncoefs = []
    for i in range(len(coefs)):
        line = []
        for j in range(len(coefs[i])):
            try:
                line.append(coefs[i][j].magnitude)
            except:
                line.append(coefs[i][j])
        ncoefs.append(line)

    nconsts = []
    for i in range(len(consts)):
        try:
            nconsts.append(consts[i].magnitude)
        except:
            nconsts.append(consts[i])

    nresult = list(np.matmul(np.linalg.inv(np.array(ncoefs)),np.array(nconsts)))

    #If we get to here, that means it was a valid linear equation
    #Now predict the units:
    resultUnits = list(np.zeros(len(nresult)))

    for i in range(len(consts)):
        if type(consts[i]) == type(1 * u.diml):
            #We can determine the unit of the answer using the coefs and const
            for j in range(len(coefs[i])):
                if type(coefs[i][j]) == type(1 * u.diml):
                    #We can use this coef and the const to find the units
                    predictedUnit = consts[i].units / coefs[i][j].units
                    if resultUnits[j] == 0:
                        resultUnits[j] = predictedUnit
                    else:
                        #We've already determined this unit, but We'll make sure it matches.
                        #If this fails, then there was a unit mismatch and the following line will yeild an error.
                        test = 1 * predictedUnit + 1 * resultUnits[j]

    for i in range(len(consts)):
        if resultUnits[i] == 0:
            raise Exception("Unable to determine units of row " + str(i) + ". Please specify more units in this row.")

    result = []
    for i in range(len(nresult)):
        result.append(nresult[i] * resultUnits[i])

    return result

def interp1D(xi,Ps,scale="linear"):
    """
    Interpolates a given "y" value based off of given points of [x,y] configuration using the specified scale.

    xi      : The x value of the point for which you'd like to interpolate the y value.
    Ps      : A list of [x,y] points between which you'll be interpolating "xi".
    scale   : The scale of the x axis upon which you're interpolating (linear or log)
    returns : The y value at "xi" interpolated between the "Ps" provided.
    """

    crossIndex = "NOT SET"

    Ps.sort(key=lambda x: x[0])

    for i in range(len(Ps)):
        if Ps[i][0] > xi:
            crossIndex = i
            break

    if crossIndex == "NOT SET":
        crossIndex = -1

    if crossIndex == 0:
        crossIndex = 1

    x1,y1 = Ps[crossIndex - 1]
    x2,y2 = Ps[crossIndex]

    if scale == "linear":
        try:
            m = (y2-y1)/(x2-x1)
            b = y1 - m * x1
            return m * xi + b
        except ZeroDivisionError:
            #infinite slope
            return np.infty
    elif scale == "log":
        try:
            a = xi - x1
            b = x2 - xi
            f = a / (a+b)
            ans = y2**f * y1**(1-f)
            return ans
        except ZeroDivisionError:
            #f = infty
            return np.infty
    else:
        raise Exception("ERROR! Unrecognized scale in interp1D.")


def interp2D(POI,SurroundingPoints):
    """
    Interpolates between four surrounding points with two input parameters and one output parameter.
    This is done by finding the 4 possible planes for each of the 3 point combinations of the 4 surroudning points provided.
    Then Each plane is evaluated at the POI and the average of those results is returned.

    POI = Point of Interest : [x,y]
    SurroundingPoints       : [[x,y,z],[x,y,z],[x,y,z],[x,y,z]]
    """
    Answers = []
    A,B,C,D = SurroundingPoints
    xi,yi = POI

    for combo in [[A,B,C],[A,B,D],[A,C,D],[B,C,D]]:
        x1,y1,z1 = combo[0]
        x2,y2,z2 = combo[1]
        x3,y3,z3 = combo[2]
        V1 = [x1 - x2 , y1 - y2 , z1 - z2]
        V2 = [x2 - x3 , y2 - y3 , z2 - z3]
        a,b,c = crossProduct(V1,V2)

        Z_interp = -1 * (a * xi + b * yi - a * x1 - b * y1 - c * z1) / c

        Answers.append(Z_interp)

    ans = A[2] * 0 #Zero, but with the right units
    for a in Answers:
        ans += a / len(Answers)
    return ans

def harmonicSum(List):
    """
    Calculates the harmonic sum of a list of numbers

    List    : A 1-D itterable data type containing the numbers you'd like to find the harmonic sum of
    returns : The harmonic sum of "List"
    """
    ans = (1 / List[0]) * 0 #Zero but with the right units.
    for i in List:
        ans += 1 / i
    return ans

def pintSolve(func,guess,args=()):
    """
    A Pint wrapper for fsolve

    func  : The function you'd like to solve
    guess : The initial guess value(s)
    args  : Any additional arguments required by the function.
    """
    raise Exception("THIS CODE IS UNDER DEVELOPMENT!")
    try:
        guess[0]
    except:
        guess = [guess,]

    guessVals  = []
    guessUnits = []
    for i in range(len(guess)):
        unit = None
        val = guess[i]
        try:
            unit = guess[i].units
            val = val.to(unit).magnitude
        except:
            unit = 1
        guessUnits.append(unit)
        guessVals.append(val)

    def wrapperFunc(var,*args):
        var = list(var)
        for i in range(len(var)):
            var[i] *= guessUnits[i]
        ans = func(var,*args)
        for i in range(len(ans)):
            try:
                ans[i] = ans[i].magnitude
            except:
                pass
        return ans
    ans = list(fsolve(wrapperFunc,guessVals,args=args))
    for i in range(len(ans)):
        ans[i] *= guessUnits[i]
    return ans

def pintQuad(func,start,stop,arguments=()):
    """
    A Pint wrapper for scipy.integrate.quad

    func      : The function you'd like to Integrate
    start     : The lower bound of integration
    strop     : The upper bound of integration
    arguments : Additional arguments required by your funciton.
    returns   : The integral of "func" with respect to the first argument of "func" from "start" to "stop"
    """
    inBaseUnit = (start.to_base_units()).units
    fBaseUnit  = (func(start,*arguments).to_base_units()).units


    def wrapperFunc(x,arguments):
        x *= inBaseUnit
        ans = func(x,*arguments)
        return ans.to_base_units().magnitude

    intAns = quad(wrapperFunc,start.to_base_units().magnitude,stop.to_base_units().magnitude,args=(arguments,))[0]
    finalAns = intAns * fBaseUnit * inBaseUnit
    return finalAns

def discreteIntegrate(Xs,Ys,method="Romb"):
    """
    A Pint wrapper for scipy.integrate.romb, simps, and trapz

    Xs      : A list of x values you'll be integrating over
    Ys      : A list of y values you'll be integrating over
    method  : The method of integration you'd like to use:
                Romb = Romberg integration
                Simp = Simpson's Rule
                Trap = Trapazoidal Rule
    returns : The integral ∫ y dx based off of the "Xs" and "Ys" provided
    """

    if len(Xs) != len(Ys):
        raise Exception("xs and ys lists do not have the same length!")

    #Create Local Copies:
    xs = []
    ys = []
    for i in range(len(Xs)):
        xs.append(Xs[i])
        ys.append(Ys[i])

    try:
        xUnit = xs[0].units
    except:
        xUnit = u.diml
    try:
        yUnit = ys[0].units
    except:
        yUnit = u.diml

    for i in range(len(xs)):
        try:
            xs[i] = xs[i].to(xUnit).magnitude
        except:
            pass
        try:
            ys[i] = ys[i].to(yUnit).magnitude
        except:
            pass

    if method == "Romb":
        ans = romb(ys,xs)
    elif method == "Simp":
        ans = simps(ys,xs)
    elif method == "Trap":
        ans = trapz(ys,xs)
    else:
        raise Exception(str(method) + " is not a recognized method.  Please specify a valid integration method: Romb,Simp,Trap")

    ans *= xUnit * yUnit
    return ans


def pintPlot(Xs,Ys,xUnit="DEFAULT",yUnit="DEFAULT",xlabel ="",ylabel ="",xScale="",yScale="",pointType = None,color = None,label =None,legend=True):
    """
    A Pint wapper for plt.plot(), plt.xlabel(), plt.ylabel(), plt.xscale(), ply.yscale()...

    Xs        : A list of Pint objects to be plotted on the x-axis
    Ys        : A list of Pint objects to be plotted on the y-axis
    xUnit     : The units of the values to be plotted on the x-axis ("DEFAULT" defaults to whatever units are passed in with the first x-axis object)
    yUnit     : The units of the values to be plotted on the y-axis ("DEFAULT" defaults to whatever units are passed in with the first y-axis object)
    xlabel    : The label of the x-axis of the graph (Units will be atuomatically added from "xUnit")(This value is passed in to plt.xlabel())
    ylabel    : The label of the y-axis of the graph (Units will be atuomatically added from "yUnit")(This value is passed in to plt.ylabel())
    xScale    : The scale of the x-axis of the graph (e.g. "log") (This is passed in to plt.xscale())
    yScale    : The scale of the y-axis of the graph (e.g. "log") (This is passed in to plt.yscale())
    pointType : Any specifications normally passed in to plt.plot(xs,ys,________) where the "_______" is.
    color     : The desired color of the graph. (This is passed in to plt.plot(xs,ys,color=__________) where the "_______" is).
    label     : The label of this line on the graph. (This is passed in to plt.plot(xs,ys,label=__________) where the "_______" is).
    legend    : A boolean of True meaning you'd like the legend of the graph shown or False meaning you'd like for the legend not to be shown.
    """

    if len(Xs) != len(Ys):
        raise Exception("xs and ys lists do not have the same length!")

    #Create Local Copies:
    xs = []
    ys = []
    for i in range(len(Xs)):
        xs.append(Xs[i])
        ys.append(Ys[i])

    if xUnit == "DEFAULT":
        xUnit = str(xs[0].units)
    if yUnit == "DEFAULT":
        yUnit = str(ys[0].units)

    if str(xUnit) != "dimensionless":
        xlabel += " (" + xUnit + ")"
    if str(yUnit) != "dimensionless":
        ylabel += " (" + yUnit + ")"

    for i in range(len(xs)):
        xs[i] = xs[i].to(xUnit).magnitude
        ys[i] = ys[i].to(yUnit).magnitude

    args = []
    kwargs = {}
    if pointType != None:
        args.append(pointType)
    if color != None:
        kwargs["color"] = color
    if label != None:
        kwargs["label"] = label
    plt.plot(xs,ys,*args,**kwargs)

    if xScale != "":
        plt.xscale(xScale)
    if yScale != "":
        plt.yscale(yScale)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if label != None and legend:
        plt.legend()

def derivative(f,x,args=()):
    """
    Approximates the derivative of the fuction f(x) at the given value for x.

    f       : The fuction which you'd like to take the derivative of
    x       : The x value of the point you'd like to take the derivative at
    args    : Any additional arguments required by "f"
    returns : The value of the derivative of "f" at "x"
    """
    x1 = x * 0.99999
    x2 = x * 1.00001
    y1 = f(x1,*args)
    y2 = f(x2,*args)
    return (y2-y1) / (x2-x1)

def avgLine(val,line1Points,line1Val,line2Points,line2Val):
    """
    Returns the coeficients of a line that lies on the weighted average(based on provided values) and the points of both of the neighboring lines.
    This can be used when interpolating based on graph values for lines that are not vertical or horizontal.

    val         : The value of the resulting line you wish to find
    line1Points : Two [x,y] points that fall on the line corresponding to line1Val
    line1Val    : The value of the line just below the "val" for which you want to find the line
    line2Points : Two [x,y] points that fall on the line corresponding to line2Val
    line2Val    : The value of the line just above the "val" for which you want to find the line

    For example, given the psychrometric chart found at https://www.engineeringtoolbox.com/docs/documents/816/psychrometric_chart_29inHg.pdf,
    this can be used to find the coeficients for the line representing the specific volume of 13.6125 ft**3 / lb of dry air:

    "val" = 13.6125 ft**3 / lb of dry air
    "line1Points" = [[75 degF, 8.5],[65 degF, 90]] ([Dry Bulb Temp., Humidity Ratio]) (Use interp1D to get these values)
    "line1Val" = 13.5 ft**3 / lb of dry air
    "line2Points" = [[95 degF, 6],[78.5 degF, 140]]
    "line2Val" = 14.0 ft**3 / lb of dry air
    This would return the coeficients of the line (m,b) respresing 13.6125 ft**3 / lb of dry air (Humidity Ratio = m * Dry Bulb Temperature + b)
    """
    P11,P12 = line1Points
    X11,Y11 = P11
    X12,Y12 = P12

    P21,P22 = line2Points
    X21,Y21 = P21
    X22,Y22 = P22

    m1 = (Y12 - Y11) / (X12 - X11)
    b1 = Y11 - m1*X11

    m2 = (Y22 - Y21) / (X22 - X21)
    b2 = Y21 - m2*X21

    perc2 = (val - line1Val) / (line2Val - line1Val)
    perc1 = 1 - perc2

    mi = perc1*m1 + perc2*m2
    bi = perc1*b1 + perc2*b2
    return [mi,bi]

def pointSlope(point,slope):
    """
    Given a point on the line, and the slope of the line, this gives you the coeficients (m,b) of the line y = m * x + b

    point : An [x,y] point found on the line
    slope : The slope of the line
    returns : The coeficients (m,b) of the line satisfying the given conditions
    """
    m = slope
    b = point[1] - m * point[0]
    return [m,b]

def pointPoint(P1,P2):
    """
    Given two points on a line, this gives you the coeficients (m,b) of the line y = m * x + b

    P1      : An [x,y] point found on the line
    P2      : An different [x,y] point found on the line
    returns : The coeficients (m,b) of the line satisfying the given conditions
    """
    m = (P2[1] - P1[1]) / (P2[0] - P1[0])
    return pointSlope(P1,m)

def lineIntersection(line1Coefs,line2Coefs):
    """
    Given two sets of coeficients (m,b) representing two (2-D) non-parallel lines, this gives you the point of interestion of the two lines.

    line1Coefs : The coeficients (m,b) of the first line
    line2Coefs : The coeficients (m,b) of the second line
    returns    : The [x,y] point where the two lines intersect
    """
    m1,b1 = line1Coefs
    m2,b2 = line2Coefs

    x = (b2 - b1) / (m1 - m2)
    y = m1*x + b1
    return [x,y]

def ReTest(Re):
    """
    This will tell you whether the specified Re falls under the Laminar, Turbulet, or Mixed reigns for a variety of scenarios.

    Re      : The Reynolds number that you'd like to analyze
    returns : None, but it prints the results to the standard output
    """
    try:
        Re = Re.to(u.diml).magnitude
    except:
        None

    print("Re =",Re,"test reults:")
    if (Re < 2100):
        print("LAMINAR (pipe)")
    elif (Re < 4000):
        print("BOTH LAMINAR AND TURBLENT (pipe)")
    else:
        print("TURBULENT (pipe)")

    if (Re < 5e5):
        print("LAMINAR (flat plate)")
    else:
        print("TURBULENT (flat plate)")

class graphAnalysis:
    """
    A functor to store and interperate graph data.

    To operate, use GraphReader.ipynb to generate the RawData array needed to sore data correctly.
    GraphReader.ipynb will give you instructions on how to quickly import all the data from any cartesian graph.
    Copy and paste the result from GraphReader to the location where you'll be instantiating a graphAnalysis object.

    Instantiate a graphAnalysis object by passing the RawData as to the constructor.
    For example:
        myObj = graphAnalysis(RawData)

    Your instantiated object (functor) now is holding all the data for the graph and is now ready to accept parameters.
    Call your object passing in the axis labels (same as you entered them in GraphReader.ipynb) and values as kwargs.
        (Note that you can pass any two dimensions from the graph. I.E. and combination of x-axis dimension, y-axis dimension, or graph series dimensions.)
    For example:
        answer = myObj(T=500*u.K,P=3*u.atm)

    The returned object (answer) is a basic object that has all of the different dimensions from the graph (axis dimensions and series dimensions) as attributes.
    For example:
        returnedT     = answer.T
        returnedP     = answer.P
        returnedVisc  = answer.visc
        returnedVhat  = answer.Vhat
        (Where T and P are axis dimensions, and visc and Vhat are series dimensions.)
    """

    class GraphLine:
        def __init__(self,value,points):
            self.value = value
            self.points = points

    class GraphSeries:
        def __init__(self,name,scale,lines):
            self.name = name
            self.scale = scale
            self.lines = lines

    def __init__(self,RawData):
        self.axisVars = [None,None]
        self.scales = [None,None]

        self.axisVars[0] = RawData[0][0]
        self.scales[0] = RawData[0][1]

        self.axisVars[1] = RawData[1][0]
        self.scales[1] = RawData[1][1]

        self.bounds = RawData[2]

        self.series = {}

        for ser in RawData[3]:
            serName = ser[0]
            serScale = ser[1]
            serLines = []
            for line in ser[2]:
                lineVal = line[0]
                linePoints = line[1]
                serLines.append(self.GraphLine(lineVal,linePoints))
            self.series[serName] = self.GraphSeries(serName,serScale,serLines)



    def readXY(self,x,y,ser):
        #To interpolate between two "Lines" in a series, I'll draw a straight line the goes through the point indicated.
        #I'll itteratively find the shortest straight line between the neightboring series lines and use the lever rule to find the right value.

        #First, I need to find the line immediatley above and below the point of interest.
        #To do this, I'll do an interp1D on each of the lines and see which one gives a y value closest to yet above/below the given y value

        #Sometimes the solver passes x or y in as a list of length 1. But you can't use the ">" operator on a list.  So this fixes it.

        yis = []
        for i in range(len(ser.lines)):
            line = ser.lines[i]
            linePoints = line.points
            yi = interp1D(x,linePoints)
            yis.append([yi,i])
        yis.sort(key = lambda x:x[0])
        bottomIndex = yis[0][1]
        topIndex = yis[1][1]    #Inititate these the lowest index.  If the answer fallws blow this range, it'll base it's extrapolation based off of the bottom values.
        for yii in yis:
            yi,i = yii
            #Update the indexs
            bottomIndex = topIndex
            topIndex = i
            #See if this crosses the given condition
            if yi > y:
                break

        bottomLine = ser.lines[bottomIndex]
        topLine   = ser.lines[topIndex]

        def getIntersections(lineCoefs, seriesLinePoints):
            m,b = lineCoefs
            intersectionIndexs = []
            for i in range(1,len(seriesLinePoints)):
                x1,y1 = seriesLinePoints[i-1]
                x2,y2 = seriesLinePoints[i]

                yi1 = m*x1+b
                yi2 = m*x2+b

                d1 = yi1 - y1
                d2 = yi2 - y2
                if (d1 < 0 and d2 < 0) or (d1 > 0 and d2 > 0):
                    #If both differences have the same sign there was not an intersection between these points.
                    pass
                else:
                    #If the signs are different, then there was in intersection
                    intersectionIndexs.append(i)

            intersections = []
            for i in intersectionIndexs:
                P1 = seriesLinePoints[i-1]
                P2 = seriesLinePoints[i]
                if (P1[0] == P2[0]):
                    P1[0] *= 0.999999
                lineiCoefs = pointPoint(P1,P2)
                intersection = lineIntersection(lineCoefs,lineiCoefs)
                intersections.append(intersection)
            return intersections

        slopeUnit = (y/x).units

        def solveShortestLineSlope(theta):
            m = np.tan(theta)
            #Guess an value for m and return the length of the line segment that passes through the given conditions and ends at the neighboring lines.

            m *= slopeUnit

            coefs = pointSlope([x,y],m)

            bottomDists = []
            bottomIntersections = getIntersections(coefs,bottomLine.points)
            for bottomIntersection in bottomIntersections:
                xbi,ybi = bottomIntersection
                #Since scales can be massivley different between the x and y axis, we'll scale each one as a fraction of the total axis length
                dx = x-xbi
                dy = y-ybi

                dx /= (self.bounds[0][1][1] - self.bounds[0][0][1])
                dy /= (self.bounds[1][1][1] - self.bounds[1][0][1])

                bottomDist = np.sqrt(dx**2 + dy**2)
                bottomDists.append(bottomDist)
            bottomDist = min(bottomDists)

            topDists = []
            topIntersections = getIntersections(coefs,topLine.points)
            for topIntersection in topIntersections:
                xti,yti = topIntersection
                dx = x-xti
                dy = y-yti

                dx /= (self.bounds[0][1][1] - self.bounds[0][0][1])
                dy /= (self.bounds[1][1][1] - self.bounds[1][0][1])

                topDist = np.sqrt((dx)**2 + (dy)**2)
                topDists.append(topDist)
            topDist = min(topDists)
            try:
                totalDist = bottomDist + topDist #If there was not a solution, it'll manifest here.  I'll return infty in that case.
            except:
                return np.infty
            return totalDist

        idealSlope = zoomSolve(solveShortestLineSlope,[[0,np.pi]])[0][0] * slopeUnit

        coefs = pointSlope([x,y],idealSlope)

        bottomDists = []
        bottomIntersections = getIntersections(coefs,bottomLine.points)
        for bottomIntersection in bottomIntersections:
            xbi,ybi = bottomIntersection
            #Since scales can be massivley different between the x and y axis, we'll scale each one as a fraction of the total axis length
            dx = x-xbi
            dy = y-ybi

            dx /= (self.bounds[0][1][1] - self.bounds[0][0][1])
            dy /= (self.bounds[1][1][1] - self.bounds[1][0][1])

            bottomDist = np.sqrt(dx**2 + dy**2)
            bottomDists.append(bottomDist)
        bottomDist = min(bottomDists)

        topDists = []
        topIntersections = getIntersections(coefs,topLine.points)
        for topIntersection in topIntersections:
            xti,yti = topIntersection
            dx = x-xti
            dy = y-yti

            dx /= (self.bounds[0][1][1] - self.bounds[0][0][1])
            dy /= (self.bounds[1][1][1] - self.bounds[1][0][1])

            topDist = np.sqrt((dx)**2 + (dy)**2)
            topDists.append(topDist)
        topDist = min(topDists)

        #Now use the "Lever Rule"
        dataPoint1 = [-bottomDist,bottomLine.value] #Distance to the bottom line, The value of the bottom line
        dataPoint2 = [topDist,topLine.value]
        returnValue = interp1D(0,[dataPoint1,dataPoint2],scale=ser.scale) #0 because the bottom point is a negative distance away from the POI and the top point is a positive distane away.
        return returnValue


    def __call__(self,**kwargs):
        """
        A function that takes in two parameters than can be found on a graph and returns all the known conditions based on those parameters.

        kwargs  : A python dictionary containing the parameter names mapping to their respective values.(Only two parameters are allowed.)
        returns : A python dictionary with all the known data at the given conditions.
        """

        if len(kwargs) != 2:
            raise Exception("ERROR! You must provide exactly 2 parameters to read off graph values!")

        #First, map the given values to the two axis dimensions.
        #guess two axis dimension values
        xg = avg([self.bounds[0][1][1],self.bounds[0][0][1]])
        yg = avg([self.bounds[1][1][1],self.bounds[1][0][1]])

        xUnit = xg.units
        yUnit = yg.units

        xg = xg.to(xUnit).magnitude
        yg = yg.to(yUnit).magnitude

        #An fsolve function the returns the error from an x-y guess pair to the inputed parameter series
        def errorFunc(xy,ser,targetVal):
            x,y = xy

            x *= xUnit
            y *= yUnit

            guessVal = self.readXY(x,y,ser)

            error = targetVal - guessVal
            return error.magnitude

        nonAxisSeries = set()
        axisSeries = set()
        for kw in kwargs:
            if kw in set(self.axisVars):
                axisSeries.add(kw)
            else:
                nonAxisSeries.add(kw)

        if len(nonAxisSeries) == 2: #Neither of the inputted parameters are axis parameters.
            seriesi = []
            targVals = []
            for kw in kwargs:
                seriesi.append(self.series[kw])
                targVals.append(kwargs[kw])

            def solveMe(xy):
                errors = []
                for i in range(len(seriesi)):
                    errors.append(errorFunc(xy,series[i],targVals[i]))
                return errors
            x,y = fsolve(solveMe,[xg,yg])
            xySolved = [x*xUnit,y*yUnit]

        elif len(nonAxisSeries) == 1: #One of the inputted parameters is an axis parameter, the other isn't
            seriesi = None
            for kw in nonAxisSeries:
                seriesi = self.series[kw]
            targValj = kwargs[kw]

            axisi = None
            for kw in axisSeries:
                axisi = kw
            targAxisi = kwargs[axisi]
            if axisi == self.axisVars[0]:
                #The xaxis is specified
                def solveMe(guessAxisj):
                    #The fact that fsolve passes this function a list is causing Pint a lot of problems.  These two commands prevent any Pint errors.
                    try:
                        guessAxisj = guessAxisj[0]
                    except:
                        pass
                    guessAxisj = float(guessAxisj)

                    xy = [targAxisi.to(xUnit).magnitude,guessAxisj]
                    return errorFunc(xy,seriesi,targValj)
                axisjSolved = fsolve(solveMe,yg)[0] * yUnit
                xySolved = [targAxisi,axisjSolved]
            elif axisi == self.axisVars[1]:
                #The yaxis is specified
                def solveMe(guessAxisj):
                    try:
                        guessAxisj = guessAxisj[0]
                    except:
                        pass
                    guessAxisj = float(guessAxisj)
                    xy = [guessAxisj,targAxisi.to(yUnit).magnitude]
                    return errorFunc(xy,seriesi,targValj)
                axisjSolved = fsolve(solveMe,xg)[0] * xUnit
                xySolved = [axisjSolved,targAxisi]
            else:
                raise Exception("ERROR! There must be a bug if this throws an error.")
        else:
            #Both of the inputted parameters are axis parameters
            axisSeries = list(axisSeries)
            if axisSeries[0] == self.axisVars[0]:
                xySolved = [kwargs[axisSeries[0]],kwargs[axisSeries[1]]]
            else:
                xySolved = [kwargs[axisSeries[1]],kwargs[axisSeries[0]]]


        #Now that I have the x and y for any variety of input, I can find the value for each series
        class answerClass:
            """
            A class to hold the values for each of the paremeters and series.
            """
            def add(self,param,val):
                self.__dict__.__setitem__(param, val)
        answer = answerClass()
        answer.add(self.axisVars[0], xySolved[0])
        answer.add(self.axisVars[1], xySolved[1])

        for ser in self.series:
            answer.add(ser, self.readXY(*xySolved,self.series[ser]))
        return answer

#SOLVER METHODS:------------------------------------------------------------------------------------
def graphError3D(solveFunc,param1Range,param2Range,res=30,args = ()):
    """
    Graphs a given two paramter Error function (like what would be passed into fsolve) and suggests where the solution could lie.
        x-axis = parameter 1 values
        y-axis = parameter 2 values
        z-axis = the sum of the error for the output of the desired function.

    solveFunc : A function that takes in a two-parameter list as it's independant variables, and other constant variables (see "args"), and returns two values that both evaluate to zero at the ideal values of the independant variables.
                This is the same type of function that would be used in fsolve to solve two equations with two unknowns.
    param1Range : An esimated range (min,max) within which the ideal value for parameter 1 could lie.  If the predicted solution lies outside of this range, this function's output will let you know.
    param2Range : An esimated range (min,max) within which the ideal value for parameter 2 could lie.  If the predicted solution lies outside of this range, this function's output will let you know.
    res         : The number of points along each axis to plot the error for.  Since there are two parameters, the total number of points plotted will be "res"**2
    args        : Any additional arguments requred be "solveFunc" beyond the two-paramter list of independant variables
    returns     : None, though it displays the graph across the specified range and comments on where the ideal solution is likely to lie.
    """
    xs = np.linspace(param1Range[0],param1Range[1],res)
    ys = np.linspace(param2Range[0],param2Range[1],res)

    Xs = []
    Ys = []
    Zs = []

    for x in xs:
        for y in ys:
            Xs.append(x)
            Ys.append(y)
            result = solveFunc([x,y],*args)
            Zs.append(np.abs(result[0]) + np.abs(result[1]))
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(Xs,Ys,Zs)
    ax.set_xlabel("Var1")
    ax.set_ylabel("Var2")
    ax.set_zlabel("Total Error")

    minPoint = [-1,-1,np.infty]

    for i in range(len(Zs)):
        if Zs[i] < minPoint[2]:
            minPoint = [Xs[i],Ys[i],Zs[i]]
    print("The point of minimum Error is:")
    print("Var1 =",minPoint[0])
    print("Var2 =",minPoint[1])
    print("Total Error =",minPoint[2])

    if (minPoint[0] == param1Range[0]):
        print("The solution likely lies at a lower Var1")
    if (minPoint[0] == param1Range[1]):
        print("The solution likely lies at a higher Var1")
    if (minPoint[1] == param2Range[0]):
        print("The solution likely lies at a lower Var2")
    if (minPoint[1] == param2Range[1]):
        print("The solution likely lies at a higher Var2")


def getCombos(linspaces):
    """
    Generates an array of all of the different combinations between the given linspaces.
        For example, all of the different combos of the linspaces [0,1] and [2,3] are [[0,2],[0,3],[1,2],[1,3]]

    linsapces : an array of multiple linspaces for which youd like to find the combos for
    returns   : an array of all possible combinations between the given linspaces
    """

    combos = []

    for element in linspaces[0]:
        combos.append([element,])

    for linspace in linspaces[1:]:
        combosBackup = []
        for combo in combos:
            combosBackup.append(combo)

        combos = []

        for element in linspace:
            for combo in combosBackup:
                newCombo = []
                for e in combo:
                    newCombo.append(e)
                newCombo.append(element)

                combos.append(newCombo)

    return combos

def zoomSolve(solveFunc,paramRanges,fracTol = 0.000001,apperture=0.5,maxItter=500,printProgress=False,args = ()):
    """
    Uses the same "graphing" approach as graphError3D() to approximate the ideal solution for a solving function.
    Additionally, this function automatically itterates to "zoom in" in the ideal solution to make it more accurate.
    It can be used with a solveFunc that solves for any number of parameters, not just two, thus it will not display the "graphs" it creates since they can be any number of dimmensions.

    solveFunc   : A function that takes an array of variables, returns an array of errors based off of the nececary functions. The same kind of function that would be passed into fsolve.
    paramRanges : An array of arrays that contain the bounds for each parameter between which the answer probably lies.  If the solver thinks the solution lies outside the specified range, it will automatically adjust.
    fracTol     : The maximum percent changes for each parameter that deems an acceptable answer.  Once percent change of each parameter between itterations all fall below this number, itteration will stop.
    apperture   : The percentage of the previous paramRange that will be the paramRange for the next itteration.  I.E. how much to zoom in after each itteration, 1 meaning no zoom in, 0 meaning infinite zoom in.
    args        : An array of additional arguments that the specified "solveFunc" takes.
    returns     : An array of the ideal parameters to fit the solveFunc and the error that those parameters return.
                  [[ANSWER],error]
    """

    oldPoint = [np.ones(len(paramRanges)) * -1, np.infty] #[[Values for each parameter],Sum of the errors yeilded by these values]

    #This should be a prime number to make sure grid points between itterations never land on top of one another.
    res = 17

    i = 0

    while True:
        i += 1

        linspaces = []
        for paramRange in paramRanges:
            linspaces.append(np.linspace(paramRange[0],paramRange[1],res))
        res += 1

        combos = getCombos(linspaces)
        errors = []

        for combo in combos:
            result = solveFunc(combo,*args)

            try:
                result[0]
            except:
                result = [result,]

            error = result[0] * 0 #Zero, but with the right units
            for element in result:
                try:
                    error += np.abs(element).magnitude
                except:
                    error += np.abs(element)
            errors.append(error)

        minPoint = [np.ones(len(paramRange)) * -1, np.infty]

        for j in range(len(errors)):
            if errors[j] < minPoint[1]:
                minPoint = [combos[j],errors[j]]

        if i > maxItter:
            print("Itterator is not making good Progress... (maxItter reached)")
            return minPoint

        if printProgress:
            print(minPoint)

        outOfBounds = False

        for i in range(len(paramRanges)):
            if (minPoint[0][i] == paramRanges[i][0]):
                shiftDist = (paramRanges[i][1] - paramRanges[i][0]) * rand.uniform(0,.25)
                paramRanges[i][0] -= shiftDist
                paramRanges[i][1] -= shiftDist
                outOfBounds = True
                break

        if outOfBounds:
            continue

        endBool = True

        for i in range(len(oldPoint[0])):
            Min = np.min([minPoint[0][i],oldPoint[0][i]])
            Max = np.max([minPoint[0][i],oldPoint[0][i]])

            e = (Max - Min) / Max
            if e > fracTol:
                endBool = False

        if endBool:
            return minPoint

        oldPoint = minPoint

        for i in range(len(paramRanges)):
            paramRanges[i][0] = minPoint[0][i] + (1 - apperture) * (paramRanges[i][0] - minPoint[0][i])
            paramRanges[i][1] = minPoint[0][i] + (1 - apperture) * (paramRanges[i][1] - minPoint[0][i])

#A Gradient of Hex color codes going from red to yellow.  Usefull when graphing colorMaps manually.
redToYellowGradient = ["#FF0000","#FF0a00","#FF1400","#FF1e00","#FF2800","#FF3200","#FF3c00","#FF4600","#FF5000","#FF5a00","#FF6400","#FF6e00","#FF7800","#FF8200","#FF8c00","#FF9600","#FFa000","#FFaa00","#FFb400","#FFbe00","#FFc800","#FFd200","#FFdc00","#FFe600","#FFf000","#FFfa00"]

#------------------------------- SYMPY FUNCTIONS ------------------------------#
def SymEval(f,var,val):
    """
    A Pint wrapper for Sympy's .subs()

    f       : A Sympy expression within which you'll be making a substitution
    var     : The Sympy variable while you'll be replacing with "val"
    val     : The value, expression, or variable you'll be putting in the place of "var"
    returns : A copy of the Sympy expression "f" with the substitution made.

    Note that all Pint values will be converted to SI units and then returned as a magnitude.  Thus the final answer will be in SI units.
    """
    try:
        return f.subs(var,val.to_base_units().magnitude)
    except:
        return f.subs(var,val)

def SymEvalList(f,key):
    """
    A function to make lots of Sympy substitutions efficiently.

    f       : A Sympy expression withing which you'll be making substitutions
    key     : A list of pairs of the Symbolic variable to substitute out and their value, expression, or variable to replace them with
              [[symVariable,value],[symVariable,value],...]
    returns : A copy of the Sympy expression "f" with the substitutions made.

    Note that all Pint values will be converted to SI units and then returned as a magnitude. Thus the final answer will be in SI units.
    """
    #List = [[symVar,val],[symVar,val],...]
    for Set in key:
        f = SymEval(f,Set[0],Set[1])
    return f

def printEQ(left,right):
    """
    A function to display an equation (two Sympy expressions, variables, or values) on one line.

    left    : The left hand side of the Equation you'd like to display
    right   : The right hand side of the Equation you'd like to display
    returns : None, but it displays the desired equation
    """
    string = "$$" + latex(left) + " = " + latex(right) + "$$"
    display(Markdown(string))

def toFunc(symExpr,symVars,finalUnit=None):
    """
    An easy way to turn a Sympy expression into a callable Python function (functor).

    symExpr : The Sympy expression you'd like to convert into a function
    symVars : A List of the Symbolic variables that will be the parameters for the new callable function
    finalUnit : The Pint unit you'd like to be assigned to the final answer
                (Since Pint does not work with Sympy, all Pint values are converted to SI units and passed as magnitudes. This allows you to re-assign the appropriate unit once the calculations are done.)
                (options)
    """
    class symFunctor():
        def __init__(self,symExpr,symVars,finalUnit):
            self.symExpr = symExpr
            self.symVars = symVars
            if finalUnit == None:
                self.finalUnit = 1
            else:
                testVar = 1 * finalUnit
                testVar = testVar.to_base_units()
                testMag = testVar.magnitude
                self.finalUnit = testVar / testMag # the value of 1 in whichever the base unit of you want the final answer to be in.

        def __call__(self,*args):
            if len(args) != len(self.symVars):
                raise Exception("Invalid number of arguments passed!")

            key = []
            for i in range(len(args)):
                key.append([self.symVars[i],args[i]])

            return SymEvalList(self.symExpr,key) * self.finalUnit
    return symFunctor(symExpr,symVars,finalUnit)

def replaceWithFunc(expr,var,funcOf):
    """
    A Function for replacing all instances of a Sympy variable within a Sympy expression with instances of a Sympy Function with the same name that takes "funcOf" as paramters.
    This is mainly only nececary when using SymODE() in which it is called automatically.

    expr    : A Sympy expression housing the Sympy variable you'd like to convert to a Sympy Function
    var     : The Sympy variable you'd like to convert to a Sympy Function
    funcOf  : A list of the Sympy variables this new "function" is a function of.
    returns : A list of
              0) expr    : A copy of the expression you passed in the with target variable switched out for the new function.
              1) newFunc : The new function of the var you passed in.
    """

    newFunc = sp.Function(var.name)(funcOf)

    expr = expr.subs(var,newFunc)
    return [expr,newFunc]

def SymODE(expr,var,WRT,bounds = None, solNumber = None):
    """
    A more concise syntax wrapper for Sympy's dsolve()
    Used to solve symbolic ODEs.

    expr      : The Sympy expression of the ODE you'd like to solve (e.g. dC/dt = expr)
    var       : The Sympy variable you're differentiating (e.g. dC/dt, var = C)
    WRT       : The Sympy variable you're differentiating "var" with respect to (e.g. dC/dt, WRT = t)
    bounds    : The boundary Values (optional) [[var,WRT],...] (This has only been tested with the base boundary conditions (I.E. var0 and WRT0). Though it may work higher order boundary conditions.  See the docs for Sympy.dsolve())
    solNumber : If the Sympy ODE solver returns mutliple solutions, this is the index of the one you want.  Leave blank to see the options.
    """

    #First, convert the variable of interest into a sympy Function
    expr, newFunc = replaceWithFunc(expr,var,WRT)

    #Now make the differential object.
    differential = newFunc.diff(WRT)

    if bounds == None:
        answer = sp.dsolve(expr - differential)
    else:
        #Format the boundary conditions appropriately
        newBounds = {}

        for b in bounds:
            varVal, WRTval = b
            newBounds[newFunc.subs(WRT,WRTval)] = varVal

        originalVars = expr.free_symbols
        answer = sp.dsolve(expr - differential,ics=newBounds)
        try:
            junk = len(answer) #Check to see if the answer is itterable and therefore has multiple solutions.
            if solNumber == None:
                print("Please select the solution you'd like.")
                for i in range(len(answer)):
                    printEQ("\,\,\,\,\,\,\,\,SOLUTION \," + str(i),answer[i])

                selection = input()
                answer = answer[int(selection)]
            else:
                answer = answer[solNumber]

        except:
            pass
        answer = answer.rhs #Right Hand Side

        answerVars = answer.free_symbols

        newSymbols = answerVars - originalVars #This is all the constants of integration.
        #Now I need to solve for these constants of integration.
        boundEQs = []
        for i in range(len(bounds)):
            varVal, WRTVal = bounds[i]

            Left  = varVal
            Right = answer.subs(WRT,WRTval)
            boundEQs.append(Left - Right)

        constSolutions = sp.solve(boundEQs,newSymbols)

        try:
            if len(constSolutions) == 1: #If Sympy just returns 1 answer (not a list of any number of answer), this will raise and error and lead to the except statement.
                constSolutions = constSolutions[0] #There is only 1 solution, take it.
            else:
                print("There are multiple solutions for evaluating at the given boudnary conditions.  Please select the index of the solution you'd like:")
                display(constSolutions)
                constSoli = input()
                constSolutions = constSolutions[constSoli]
        except:
            pass #There is only one solution

        for key in constSolutions:
            answer = answer.subs(key,constSolutions[key])

        return answer.simplify()

class SolutionDict():
    """
    Often when doing lengthy solutions in Sympy it's cumbersome to use a million var = var.subs() calls, var = var.simplify() calls, var = sp.solve(EQs,var) calls, etc.
    I developed a more consice syntax for handling this process.

    All Sympy expressions will be saved in a Python dictionary and custom funcitons will be made to act on elements of this dictionary directly.
    Thus any var = var.func() calls are simplified to func(var) calls.

    Additionally, this custom Python dictionary wrapper automatically prints the output, equation number, and source of each modification.

    Using this method, I've reduced pages of bulky, unreadable code down to a few lines of intuative code tha prints easy to follow proofs of the process.

    Just make sure the call the clear() function after each code block so that the solution for one problem does not interefere with solutions for other problems.
    """
    def __init__(self,printSteps=True):
        self.dict = {}
        self.printSteps = printSteps
        self.count = 0

    def __setitem__(self,key,item,printProgress=None):
        oldNumber = None
        if key in self.dict:
            oldNumber = self.dict[key][1]
        self.count += 1
        self.dict[key] = [item,self.count]
        if printProgress == None:
            if self.printSteps:
                self.print(key)
        else:
            if printProgress:
                self.print(key)

    def __getitem__(self,key):
        return self.dict[key][0]

    def print(self,key,oldNumber=None,subsNumbers=None,subsTarget=None,newNumber=None):
        if newNumber == None:
            newNumber = self.dict[key][1]

        if oldNumber == None:
            string = "<div style=\"text-align: right\"> $Eq\;" + str(newNumber) + "$</div> " + "$" + latex(key) + " = " + latex(self.dict[key][0]) + "$"
        elif subsNumbers == None:
            string = "<div style=\"text-align: right\"> $Eq\;" + str(newNumber) + "\;\leftarrow\;Eq\;" + str(oldNumber) + "$</div>" + "$" + latex(key) + " = " + latex(self.dict[key][0]) + "$"
        else:
            subsNames = str()

            if type(subsNumbers) == type([1,2]):
                #Multiple substitutions

                if subsTarget != None:
                    raise Exception("SYNTAX ERROR! When printing mutiple substitutions, pass them as subsNumbers=[[num,targ],[num,targ],etc]. Do not use subsTarget")

                for pair in subsNumbers:
                    try:
                        subsNumber, subsTarget = pair
                    except:
                        subsNumber = pair
                        subsTarget = None

                    if subsTarget == None:
                        subsNames += "Eq\;" + str(subsNumber) + ",\;"
                    else:
                        subsNames += latex(subsTarget) + "=" + latex(subsNumber) + ",\;"
                subsNames = subsNames[:-3] #to get rid of the last ,\;
            else:
                #There's only one substitution
                subsNumber = subsNumbers
                if subsTarget == None:
                    subsNames += "Eq\;" + str(subsNumber)
                else:
                    subsNames += latex(subsTarget) + "=" + latex(subsNumber)

            string = "<div style=\"text-align: right\"> $Eq\;" + str(newNumber) + "\;:\;Eq\;" + str(oldNumber) + "\;\leftarrow\;" + subsNames + "$</div>" + "$" + latex(key) + " = " + latex(self.dict[key][0]) + "$"
        display(Markdown(string))

    def clear(self):
        self.dict.clear()
        self.count = 0

solvs = SolutionDict() #Now make a global instantiation of the SolutionDict, named concisely for easy use when working with long solutions.
solvs.clear() #Make sure the reset the solvs dict between problems!

def subs(solution,targets,replacements=None):
    """
    Instead of saying
        solution = solution.subs(target,replacement)
    Just call
        subs(solution,target,replacememt)

    If the replacememt(s) for the target(s) is(are) already found elsewhere in solvs under the same name (the name of the repalcement), leave repalcement = None and the code will find the replacememt automatically.

    For example:
    solvs[y] = m*x+b  #Eq 1
    solvs[z] = y**2   #Eq 2

    To substitute y (solution contained in Eq 1) into Eq2, just call "subs(z,y)" and it will autofill the replacement with solvs[y] (I.E. Eq 1)

    Just make sure the solution is stored in the sovls dict (solvs[solution]) since you can't make a substitution on an object you haven't defined yet.
    """

    if not isinstance(targets,(list,tuple,np.ndarray)): #There is only one substitution being made in this call.
        target = targets
        replacement = replacements

        if solution in solvs.dict:
            subsObj = None
            targObj = None
            if replacement == None: #Track down the pre-defined replacement for the target.
                if target in solvs.dict:
                    replacement,subsObj = solvs.dict[target]
                    targObj = None
                else:
                    raise Exception("Could not find saved solution for " + str(target) + ".  Please specify a replacement value.")
            else:
                subsObj = replacement
                targObj = target

            solutionOldNumber = solvs.dict[solution][1]
            solvs.__setitem__(solution,solvs[solution].subs(target,replacement),False)
            solvs.print(solution,solutionOldNumber,subsObj,targObj,solvs.count)
        else:
            raise Exception("Target Solution not found.  Could not execute substitution.")
    else:
        #This is a multi subs
        solutionOldNumber = solvs.dict[solution][1]
        if replacements != None:
            raise Exception("INVALID SYNTAX!  If you're passing in mutliple substitutions, pass them in as targets=[[target,replacememt],[target,replacememt],etc], not as targets=[targets],replacememts=[replacememts]")

        subsPairs = []

        for i in range(len(targets)):
            pair = targets[i]
            try:
                target,replacement = pair
            except:
                target = pair
                replacement = None

            lastOne = False
            if i == len(targets) - 1:
                lastOne = True

            if solution in solvs.dict:
                subsObj = None
                targObj = None
                if replacement == None:
                    if target in solvs.dict:
                        replacement,subsObj = solvs.dict[target]
                        targObj = None
                    else:
                        raise Exception("Could not find saved solution for " + str(target) + ".  Please specify a replacement value.")
                else:
                    subsObj = replacement
                    targObj = target
                if not lastOne:
                    solvs.count -= 1 #Don't itterate the number
                solvs.__setitem__(solution,solvs[solution].subs(target,replacement),False)
                subsPairs.append([subsObj,targObj])
                #solvs.print(solution,solvs.dict[solution][1]-1,subsObj,targObj,solvs.count)
            else:
                raise Exception("Target Solution not found.  Could not execute substitution.")

        solvs.print(solution,solutionOldNumber,subsPairs,None,solvs.count)


def simplify(solution):
    """
    Instead of calling
        solution = solution.simplify()
    call
        simplify(solution).

    Just make sure the solution is stored in the sovls dict (solvs[solution]) since you can't simplify an object you haven't defined yet.
    """
    oldNum = solvs.dict[solution][1]
    solvs.__setitem__(solution,solvs[solution].simplify(),False)
    newNum = solvs.dict[solution][1]
    solvs.print(solution,oldNum)

def solve(Eqs,Vars,solNumber=None):
    """
    Instead of calling
        Vars = sp.solve(Eqs,Vars)
    call
        solve(Eqs,Vars)

    If multiple solutions result from the Sympy solver, you'll be prompted to select the one you want unless you specify the solution number as solNumber.
    """
    try:
        len(Eqs)
    except:
        Eqs = [Eqs,]
    try:
        len(Vars)
    except:
        Vars = {Vars,}

    oldLabelStr = "Solving\n$"

    for i in range(len(Eqs)):
        if Eqs[i] in solvs.dict:
            oldLabelStr += "Eq\;" + str(solvs.dict[Eqs[i]][1]) + "\;\;\;"
            Eqs[i] = solvs[Eqs[i]] - Eqs[i]
        else:
            oldLabelStr += "0=" + latex(Eqs[i]) + "\;\;\;"
    oldLabelStr += "$\nfor\n$" + latex(Vars) + "$"

    oldLabelStr = "<div style=\"text-align: center\">" + oldLabelStr + "</div>"
    display(Markdown(oldLabelStr))

    ans = sp.solve(Eqs,Vars)
    if len(ans) == 0:
        raise Exception("ERROR! No Solution Found.")
    elif str(type(ans)) == "<class 'dict'>":
        pass
    elif len(ans) > 1:
        if solNumber != None:
            ans = ans[solNumber]
        else:
            print("Please select the solution you'd like.")
            for i in range(len(ans)):
                printEQ("\,\,\,\,\,\,\,\,SOLUTION \," + str(i),ans[i])

            selection = input()
            ans = ans[int(selection)]
    else:
        ans = ans[0]

    for i in ans:
        solvs[i] = ans[i]



#----------------- PHYSICAL CHEMISTRY CONSTANTS AND FUNCTIONS -----------------#
#Boltzmann constant
kb   = 1.380649E-23 * u.J / u.K

#Planck constant
h    = 6.62607015E-34 * u.J * u.sec
hbar = h / (2*np.pi)

#Permittivity of free space
eps0 = 8.854187E-12 * u.C**2 / u.J / u.m

#Bohr radius
a0   = 5.29177E-11 * u.m

#Elementary charge
e    = 1 * u.e

#Stefan-Boltzmann constant
sbc = σ = np.pi**2 * kb**4 / (60 * hbar**3 * c**2)

def toWavenumber(freq):
    """
    A funciton for converting a frequency (expressed in 1 / time dimmensionality) to wavenumber (expressed in 1 / length dimmensionality)
    """
    wv = freq / c
    return wv.to(1/u.cm)

def fromWavenumber(wv):
    """
    A funciton for converting wavenumber (expressed in 1 / length dimmensionality) to a frequency (expressed in 1 / time dimmensionality)
    """
    freq = wv*c
    return(freq.to(1/u.sec))

#------------------- THERMODYNAMICS CONSTANTS AND FUNCTIONS -------------------#
def cubicEOS(Tr,Pr,omega = "NotEntered",eos = "vdW",Zguesses=np.logspace(-5,2,7)):
    """
    A wrapper for a variety of equations of state.

    Tr : Relative Temperature (T / Tcritical)
    Pr : Relative Pressure (P / Pcritical)
    omega : Acentric factor for the deisred species.
    eos : The specific equation of state you'd like to use.  Select the code from the following:
        vdW    : Van der Waals
        RK     : Redlich / Kwong
        SRK    : Soave / Redlich / Kwong
        PR     : Peng / Robinson
        Ideal  : Ideal Gas Law
    Zguesses : A list of guesses of possible Z values (compressibility factors) since this function can return all 3 possbile compressibility factors
    returns : A list of compressibility factors (Z) the result from the specified equation of state at the given conditions
    """
    try: #Turn pint values to just magnitudes.
        Tr = Tr.to(u.diml).magnitude
        Pr = Pr.to(u.diml).magnitude
    except:
        doNothing = "Do Nothing"

    # Values from
    # https://en.wikipedia.org/wiki/Equation_of_state

    if eos == "vdW":
        alpha   = 1
        sigma   = 0
        epsilon = 0
        Omega   = 1/8
        Psi     = 27/64
    elif eos == "RK":
        alpha   = Tr ** (-0.5)
        sigma   = 1
        epsilon = 0
        Omega   = 0.08664
        Psi     = 0.42748
    elif eos == "SRK":
        alpha   = (1 + (0.48508 + 1.55171 * omega - 0.15613 * omega**2) * (1 - np.sqrt(Tr)))**2
        sigma   = 1
        epsilon = 0
        Omega   = 0.08664
        Psi     = 0.42747
    elif eos == "PR":
        alpha   = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega**2) * (1 - np.sqrt(Tr)))**2
        sigma   = 1 + np.sqrt(2)
        epsilon = 1 - np.sqrt(2)
        Omega   = 0.07780
        Psi     = 0.45724
    elif eos == "Ideal":
        return 1
    else:
        errorMsg = ""
        errorMsg += "ERROR!!! CubicEOS regime (" + eos + ") not recognized.  Please input a vaild regime:\n"
        errorMsg += "\"vdW\"   : Van der Waals\n"
        errorMsg += "\"RK\"    : Redlich / Kwong\n"
        errorMsg += "\"SRK\"   : Soave / Redlich / Kwong\n"
        errorMsg += "\"PR\"    : Peng / Robinson\n"
        errorMsg += "\"Ideal\" : Ideal Gas Law"
        raise Exception(errorMsg)

    Beta = Omega * Pr / Tr
    q    = Psi * alpha / Omega / Tr

    def SolveCubicEOS(Z):
        Left  = Z
        Right = 1 + Beta - q * Beta * (Z - Beta) / (Z + epsilon * Beta) / (Z + sigma * Beta)
        return np.abs(Left - Right)

    Zs = []
    for z in Zguesses:
        newZ = fsolve(SolveCubicEOS,z)[0]
        add = True

        if newZ <= 0:
            add = False
        for zi in Zs:
            percError = (zi - newZ) / max([np.abs(zi),np.abs(newZ)])
            #print(newZ, zi, percError)
            if np.abs(percError) < 0.001:
                add = False
        if add:
            Zs.append(newZ)
    return Zs

def LeeKeslerTable(POI,Vals):
    """
    A function that simplifies 2D interpolation based off of a Lee Kesler Table

    A typical Lee Kesler table is orangized like follows:
       ...P1 ,P2
            ⋮   ⋮
    T1 ...Val,Val
    T2 ...Val,Val

    where the conditions you have are between T1 and T2, P1 and P2.

    POI     : A Tuple (P, T) with your Pressure and Temperature conditions.
    Vals    : A list of the following organization: (It should visually look just like the table)
              [[0 ,P1 ,P2 ],
               [T1,Val,Val],
               [T2,Val,Val]]
    returns : The 2D interpolated value based of the Lee Kesler data and the given conditions.
    """
    P1 = Vals[0][1]
    P2 = Vals[0][2]
    T1 = Vals[1][0]
    T2 = Vals[2][0]

    A = [P1,T1,Vals[1][1]]
    B = [P2,T1,Vals[1][2]]
    C = [P1,T2,Vals[2][1]]
    D = [P2,T2,Vals[2][2]]

    return interp2D(POI,[A,B,C,D])

class Water:
    """
    A Pint / syntax wrapper for the IAPWS97 module.

    To opperate, create a water class object with two specified conditions from the following list: (as Pint objets)
        T     : Temperature
        P     : Pressure
        x     : Steam quality
        H     : Specific enthalpy
        S     : Specific entropy
        l     : wavelength of light emmited
        V_hat : Specific volume

    The class will handle all the unit conversions to match up with the IAPWS module.

    To get a value, just us the water class object .____ for each property you want from the following list: (They are returned as Pint objects)
        H      : Specific enthalpy
        S      : Specific entropy
        P      : Pressure
        T      : Temperature
        x      : Steam quality
        l      : wavelength of light emmited
        V_hat  : Specific volume
        g      : Specific Gibbs free energy
        a      : Specific Helmholtz free energy
        rho    : density
        U      : Specific Internal energy
        Cp     : Isobaric heat capacity
        Cv     : Isochoric heat capacity
        Z      : Compresion factor
        phi    : Fugacity coefficient
        f      : Fugacity
        xkappa : Isothermal compressibility
        kappas : Adiabatic compressibility
        mu     : Dynamic viscosity
        nu     : Kinematic viscosity
        k      : Thermal conductivity
        alpha  : Thermal diffusivity
        sigma  : Surface tension
        Prandt : Prandtl number
        Hvap   : Heat of vaporization
        Svap   : Entropy of vaporization
        Tc     : Crititcal temperature
        Pc     : Crititcal pressure
        Vc     : Crititcal specific volume
        Tr     : Relative temperature
        Pr     : Relative pressure
        Vr     : Relative specifc volume
        MW     : Molecular weight
        omega  : Acentric factor


    See https://pypi.org/project/iapws/ for IAPWS documentation
    """
    def __init__(self,**kwargs):
        setupMatrix = [
            #[self.Property, steam.Property, unit]
            ["H","h",u.kJ / u.kg],
            ["S","s",u.kJ / u.kg / u.K],
            ["P","P",u.MPa],
            ["T","T",u.K],
            ["x","x",1],
            ["l","l",u.nm],
            ["V_hat","v",u.m**3 / u.kg],
            ["g","g",u.kJ / u.kg],
            ["a","a",u.kJ / u.kg],
            ["rho","rho",u.kg / u.m**3],
            ["U","u",u.kJ / u.kg],
            ["Cp","cp",u.kJ / u.kg / u.K],
            ["Cv","cv",u.kJ / u.kg / u.K],
            ["Z","Z",1],
            ["phi","fi",1],
            ["f","f",u.MPa],
            ["xkappa","xkappa",1/u.MPa],
            ["kappas","kappas",1/u.MPa],
            ["mu","mu",u.Pa / u.sec],
            ["nu","nu",u.m**2 / u.sec],
            ["k","k",u.W / u.m / u.K],
            ["alpha","alpha",u.m**2 / u.sec],
            ["sigma","sigma",u.N / u.m],
            ["Prandt","Prandt",1],
            ["Hvap","Hvap",u.kJ / u.kg],
            ["Svap","Svap",u.kJ / u.kg / u.K],
            ["Tc","Tc",u.K],
            ["Pc","Pc",u.MPa],
            ["Tr","Tr",1],
            ["Pr","Pr",1],
        ]

        iapwsArgs = {}

        for kw in kwargs: #Assign to the right units.
            kwFound = False
            for Set in setupMatrix:
                if Set[0] == kw:
                    kwFound = True
                    iapws_kw = Set[1]
                    unit = Set[2]
                    if unit == 1:
                        break

                    #Change the Water class keyword to an IAPWS keyword, and set to the correct unit.
                    iapws_value = kwargs[kw].to(unit).magnitude
                    iapwsArgs[iapws_kw] = iapws_value
                    break
            if not kwFound:
                raise Exception("The keyword \"" + str(kw) + "\" is not recognized.  Use help(Water) to see a list of acceptable keywords.")
        self.steam = IAPWS97(**iapwsArgs)

        self.MW = 18.01528 * u.g / u.mol
        self.omega = 0.344861

        propertyFound = False

        for Set in setupMatrix:
            myProperty, steamProperty, unit = Set
            try:
                self.__dict__.__setitem__(myProperty,self.steam.__dict__.__getitem__(steamProperty) * unit)
                propertyFound = True
            except:
                self.__dict__.__setitem__(myProperty,"NotFound")

        if not propertyFound:
            raise Exception("There is no data for the conditions given.")

#---------------------- TRANSPORT CONSTANTS AND FUNCTIONS ---------------------#
def convertDab(Dabref,T,P,Tref,Pref=1*u.atm):
    """
    A function to calculate the diffusivity for the conditions T and P when you know the diffusivity from the conditions Tref and Pref

    Dabref  : The value for Dab at Tref and Pref
    T       : The temperature of the desired conditions
    P       : the pressure of the desired conditions
    Tref    : The temperature of the reference condition
    Pref    : The pressure of the reference condition
    returns : The value of Dab at T and P
    """
    return Dabref / (Pref**-1 * Tref**(3/2)) * (P**-1 * T**(3/2))


class Air():
    """
    A class for housing the physical properties of air
    Each property can be referenced by the following attributes:
        Tc    : Crititcal temperature
        Pc    : Crititcal pressure
        Vc    : Crititcal specific volume
        Zc    : Crititcal compressibility factor
        omega : Acentric factor
        MW    : Molecular weight
        Cp    : Isobaric ideal gas heat capacity (A Function that requires a specified temperature)
        Cv    : Isochoric ideal gas heat capacity (A Function that requires a specified temperature)
        k     : Thermal conductivity (A function that requires a specified temperature and pressure)
        mu    : Dynamic viscosity (A function that requires a specified temperature and pressure)
        rho   : Density (A function that requires a specified temperature and pressure)
        nu    : Kinematic viscosity (A function that requires a specified temperature and pressure)
        alpha : Thermal diffusivity (A function that requires a specified temperature and pressure)
        Pr    : Prandtl number (A function that requires a specified temperature and pressure)
    """
    #Sources:
    #(1) https://www.engineeringtoolbox.com/critical-point-d_997.html
    #(2) http://webserver.dmt.upm.es/~isidoro/dat1/eGAS.pdf
    #(3) https://www.engineeringtoolbox.com/air-composition-d_212.html
    #(4) https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
    #(5) https://www.engineeringtoolbox.com/docs/documents/1509/Air%20thermal%20conductivity%20temperature%20C.jpg
    #(6) https://www.engineeringtoolbox.com/docs/documents/601/Air%20Dynamic%20viscosity%20pressure%20C.jpg
    #(7) https://www.engineeringtoolbox.com/docs/documents/600/Air%20density%201atm%20temp%20F.jpg
    #(8) https://www.engineeringtoolbox.com/docs/documents/600/Air%20density%20pressure%20temp%20C.jpg

    def __init__(self):
        self.Tc = 132.63 * u.K # (1)

        self.Pc = 37.36 * u.atm # (1)

        self.Vc = 92.35 * u.cm**3 / u.mol # (1)

        self.Zc = self.Pc * self.Vc / R / self.Tc

        self.omega = 0.035 # (2)

        self.MW = 28.9647 * u.g / u.mol # (3)

        #Data found at (4)
        self.Ts = np.append(np.arange(250,801,50),np.arange(900,1501,100)) * u.K
        self.Cps = np.array([1.003,1.005,1.008,1.013,1.020,1.029,1.040,1.051,1.063,1.075,1.087,1.099,1.121,1.142,1.155,1.173,1.190,1.204,1.216]) * u.kJ / u.kg / u.K
        self.Cvs = np.array([0.716,0.718,0.721,0.726,0.733,0.742,0.753,0.764,0.776,0.788,0.800,0.812,0.834,0.855,0.868,0.886,0.903,0.917,0.929]) * u.kJ / u.kg / u.K

        self.CpPoints = []
        self.CvPoints = []

        for i in range(len(self.Ts)):
            self.CpPoints.append([self.Ts[i],self.Cps[i]])
            self.CvPoints.append([self.Ts[i],self.Cvs[i]])

        #Data taken from (5) (Thermal conductivity as a function of temperature and pressure)
        #See GraphReader.ipynb to see how this data was generated
        kGraphData = [
        	["T","linear"], 	#x-axis definition
        	["k","linear"], 	#y-axis definition
        	[[[1068,73.14999999999998 * u("kelvin")], [1788,1873.15 * u("kelvin")]], #x bounds
        	 [[258,160.0 * u("milliwatt / kelvin / meter")], [619,0.0 * u("milliwatt / kelvin / meter")]]],[ #y bounds

        	["P","log",[	#series definition
        		[1.0 * u("standard_atmosphere"), #Line definition
        			[[85.65000000000009 * u("kelvin"), 158.22714681440448 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 110.80332409972303 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 110.80332409972303 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 80.22160664819947 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 80.22160664819947 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 11.523545706371237 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 8.864265927977897 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 8.864265927977897 * u("milliwatt / kelvin / meter")], [90.65000000000009 * u("kelvin"), 7.0914127423823174 * u("milliwatt / kelvin / meter")], [95.65000000000009 * u("kelvin"), 7.0914127423823174 * u("milliwatt / kelvin / meter")], [95.65000000000009 * u("kelvin"), 7.0914127423823174 * u("milliwatt / kelvin / meter")], [123.15000000000009 * u("kelvin"), 11.080332409972357 * u("milliwatt / kelvin / meter")], [123.15000000000009 * u("kelvin"), 11.080332409972357 * u("milliwatt / kelvin / meter")], [163.1500000000001 * u("kelvin"), 14.626038781163459 * u("milliwatt / kelvin / meter")], [205.6500000000001 * u("kelvin"), 18.171745152354617 * u("milliwatt / kelvin / meter")], [243.1500000000001 * u("kelvin"), 21.27423822714684 * u("milliwatt / kelvin / meter")], [243.1500000000001 * u("kelvin"), 21.27423822714684 * u("milliwatt / kelvin / meter")], [283.1500000000001 * u("kelvin"), 23.933518005540208 * u("milliwatt / kelvin / meter")], [283.1500000000001 * u("kelvin"), 23.933518005540208 * u("milliwatt / kelvin / meter")], [323.1500000000001 * u("kelvin"), 27.479224376731338 * u("milliwatt / kelvin / meter")], [363.1500000000001 * u("kelvin"), 29.695290858725798 * u("milliwatt / kelvin / meter")], [400.6500000000001 * u("kelvin"), 32.35457063711914 * u("milliwatt / kelvin / meter")], [440.6500000000001 * u("kelvin"), 35.01385041551251 * u("milliwatt / kelvin / meter")], [483.1500000000001 * u("kelvin"), 38.11634349030476 * u("milliwatt / kelvin / meter")], [523.1500000000001 * u("kelvin"), 40.33240997229922 * u("milliwatt / kelvin / meter")], [560.6500000000001 * u("kelvin"), 42.99168975069256 * u("milliwatt / kelvin / meter")], [560.6500000000001 * u("kelvin"), 42.99168975069256 * u("milliwatt / kelvin / meter")], [600.6500000000001 * u("kelvin"), 45.6509695290859 * u("milliwatt / kelvin / meter")], [643.1500000000001 * u("kelvin"), 47.86703601108036 * u("milliwatt / kelvin / meter")], [683.1500000000001 * u("kelvin"), 50.08310249307482 * u("milliwatt / kelvin / meter")], [878.1500000000001 * u("kelvin"), 60.72022160664824 * u("milliwatt / kelvin / meter")], [1078.15 * u("kelvin"), 71.35734072022166 * u("milliwatt / kelvin / meter")], [1278.15 * u("kelvin"), 80.66481994459838 * u("milliwatt / kelvin / meter")], [1278.15 * u("kelvin"), 80.66481994459838 * u("milliwatt / kelvin / meter")], [1478.15 * u("kelvin"), 90.41551246537401 * u("milliwatt / kelvin / meter")], [1675.65 * u("kelvin"), 99.27977839335185 * u("milliwatt / kelvin / meter")], [1875.65 * u("kelvin"), 108.14404432132969 * u("milliwatt / kelvin / meter")], [1875.65 * u("kelvin"), 108.14404432132969 * u("milliwatt / kelvin / meter")]]],
        		[5.0 * u("standard_atmosphere"), #Line definition
        			[[88.15000000000009 * u("kelvin"), 158.67036011080336 * u("milliwatt / kelvin / meter")], [98.15000000000009 * u("kelvin"), 129.86149584487538 * u("milliwatt / kelvin / meter")], [98.15000000000009 * u("kelvin"), 129.86149584487538 * u("milliwatt / kelvin / meter")], [108.15000000000009 * u("kelvin"), 106.81440443213299 * u("milliwatt / kelvin / meter")], [110.65000000000009 * u("kelvin"), 99.72299168975073 * u("milliwatt / kelvin / meter")], [108.15000000000009 * u("kelvin"), 59.83379501385045 * u("milliwatt / kelvin / meter")], [108.15000000000009 * u("kelvin"), 12.85318559556788 * u("milliwatt / kelvin / meter")], [110.65000000000009 * u("kelvin"), 9.750692520775658 * u("milliwatt / kelvin / meter")], [183.1500000000001 * u("kelvin"), 16.84210526315792 * u("milliwatt / kelvin / meter")], [313.1500000000001 * u("kelvin"), 26.59279778393355 * u("milliwatt / kelvin / meter")], [313.1500000000001 * u("kelvin"), 26.59279778393355 * u("milliwatt / kelvin / meter")], [483.1500000000001 * u("kelvin"), 38.55955678670364 * u("milliwatt / kelvin / meter")], [683.1500000000001 * u("kelvin"), 50.52631578947373 * u("milliwatt / kelvin / meter")], [878.1500000000001 * u("kelvin"), 61.60664819944603 * u("milliwatt / kelvin / meter")], [1078.15 * u("kelvin"), 71.35734072022166 * u("milliwatt / kelvin / meter")], [1278.15 * u("kelvin"), 81.10803324099726 * u("milliwatt / kelvin / meter")], [1478.15 * u("kelvin"), 90.41551246537401 * u("milliwatt / kelvin / meter")], [1680.65 * u("kelvin"), 99.72299168975073 * u("milliwatt / kelvin / meter")], [1865.65 * u("kelvin"), 108.14404432132969 * u("milliwatt / kelvin / meter")]]],
        		[10.0 * u("standard_atmosphere"), #Line definition
        			[[88.15000000000009 * u("kelvin"), 157.78393351800557 * u("milliwatt / kelvin / meter")], [113.15000000000009 * u("kelvin"), 108.14404432132969 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 79.33518005540171 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 51.41274238227152 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 51.41274238227152 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 20.831024930747958 * u("milliwatt / kelvin / meter")], [118.15000000000009 * u("kelvin"), 11.080332409972357 * u("milliwatt / kelvin / meter")], [183.1500000000001 * u("kelvin"), 16.398891966759038 * u("milliwatt / kelvin / meter")], [285.6500000000001 * u("kelvin"), 24.37673130193909 * u("milliwatt / kelvin / meter")], [483.1500000000001 * u("kelvin"), 38.55955678670364 * u("milliwatt / kelvin / meter")], [678.1500000000001 * u("kelvin"), 50.52631578947373 * u("milliwatt / kelvin / meter")], [880.6500000000001 * u("kelvin"), 61.60664819944603 * u("milliwatt / kelvin / meter")], [880.6500000000001 * u("kelvin"), 61.60664819944603 * u("milliwatt / kelvin / meter")], [1085.65 * u("kelvin"), 71.80055401662054 * u("milliwatt / kelvin / meter")], [1283.15 * u("kelvin"), 81.55124653739617 * u("milliwatt / kelvin / meter")], [1498.15 * u("kelvin"), 91.74515235457068 * u("milliwatt / kelvin / meter")], [1695.65 * u("kelvin"), 100.60941828254852 * u("milliwatt / kelvin / meter")], [1858.15 * u("kelvin"), 108.14404432132969 * u("milliwatt / kelvin / meter")]]],
        		[20.0 * u("standard_atmosphere"), #Line definition
        			[[88.15000000000009 * u("kelvin"), 157.78393351800557 * u("milliwatt / kelvin / meter")], [105.65000000000009 * u("kelvin"), 117.45152354570641 * u("milliwatt / kelvin / meter")], [118.15000000000009 * u("kelvin"), 98.39335180055406 * u("milliwatt / kelvin / meter")], [128.1500000000001 * u("kelvin"), 70.02770083102496 * u("milliwatt / kelvin / meter")], [130.6500000000001 * u("kelvin"), 29.252077562326917 * u("milliwatt / kelvin / meter")], [130.6500000000001 * u("kelvin"), 18.614958448753498 * u("milliwatt / kelvin / meter")], [135.6500000000001 * u("kelvin"), 15.955678670360157 * u("milliwatt / kelvin / meter")], [145.6500000000001 * u("kelvin"), 15.06925207756234 * u("milliwatt / kelvin / meter")], [250.6500000000001 * u("kelvin"), 21.717451523545748 * u("milliwatt / kelvin / meter")], [408.1500000000001 * u("kelvin"), 33.68421052631584 * u("milliwatt / kelvin / meter")], [538.1500000000001 * u("kelvin"), 42.10526315789477 * u("milliwatt / kelvin / meter")], [683.1500000000001 * u("kelvin"), 50.96952908587261 * u("milliwatt / kelvin / meter")], [880.6500000000001 * u("kelvin"), 61.60664819944603 * u("milliwatt / kelvin / meter")], [1080.65 * u("kelvin"), 71.35734072022166 * u("milliwatt / kelvin / meter")], [1278.15 * u("kelvin"), 81.10803324099726 * u("milliwatt / kelvin / meter")], [1480.65 * u("kelvin"), 90.85872576177289 * u("milliwatt / kelvin / meter")], [1673.15 * u("kelvin"), 99.72299168975073 * u("milliwatt / kelvin / meter")], [1870.65 * u("kelvin"), 108.58725761772857 * u("milliwatt / kelvin / meter")]]],
        		[50.0 * u("standard_atmosphere"), #Line definition
        			[[88.15000000000009 * u("kelvin"), 157.3407202216067 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 108.58725761772857 * u("milliwatt / kelvin / meter")], [128.1500000000001 * u("kelvin"), 83.32409972299172 * u("milliwatt / kelvin / meter")], [143.1500000000001 * u("kelvin"), 58.94736842105266 * u("milliwatt / kelvin / meter")], [143.1500000000001 * u("kelvin"), 56.28808864265932 * u("milliwatt / kelvin / meter")], [153.1500000000001 * u("kelvin"), 41.66204986149589 * u("milliwatt / kelvin / meter")], [168.1500000000001 * u("kelvin"), 24.819944598337997 * u("milliwatt / kelvin / meter")], [175.6500000000001 * u("kelvin"), 22.603878116343537 * u("milliwatt / kelvin / meter")], [205.6500000000001 * u("kelvin"), 22.16066481994463 * u("milliwatt / kelvin / meter")], [278.1500000000001 * u("kelvin"), 26.59279778393355 * u("milliwatt / kelvin / meter")], [460.6500000000001 * u("kelvin"), 38.11634349030476 * u("milliwatt / kelvin / meter")], [680.6500000000001 * u("kelvin"), 50.96952908587261 * u("milliwatt / kelvin / meter")], [878.1500000000001 * u("kelvin"), 61.60664819944603 * u("milliwatt / kelvin / meter")], [1080.65 * u("kelvin"), 71.80055401662054 * u("milliwatt / kelvin / meter")], [1278.15 * u("kelvin"), 81.55124653739617 * u("milliwatt / kelvin / meter")], [1478.15 * u("kelvin"), 90.85872576177289 * u("milliwatt / kelvin / meter")], [1673.15 * u("kelvin"), 99.72299168975073 * u("milliwatt / kelvin / meter")], [1870.65 * u("kelvin"), 108.58725761772857 * u("milliwatt / kelvin / meter")]]],
        		[100.0 * u("standard_atmosphere"), #Line definition
        			[[88.15000000000009 * u("kelvin"), 153.79501385041556 * u("milliwatt / kelvin / meter")], [95.65000000000009 * u("kelvin"), 138.28254847645434 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 108.14404432132969 * u("milliwatt / kelvin / meter")], [140.6500000000001 * u("kelvin"), 77.56232686980613 * u("milliwatt / kelvin / meter")], [175.6500000000001 * u("kelvin"), 47.42382271468148 * u("milliwatt / kelvin / meter")], [213.1500000000001 * u("kelvin"), 29.252077562326917 * u("milliwatt / kelvin / meter")], [258.1500000000001 * u("kelvin"), 29.252077562326917 * u("milliwatt / kelvin / meter")], [280.6500000000001 * u("kelvin"), 29.695290858725798 * u("milliwatt / kelvin / meter")], [360.6500000000001 * u("kelvin"), 34.12742382271472 * u("milliwatt / kelvin / meter")], [483.1500000000001 * u("kelvin"), 40.7756232686981 * u("milliwatt / kelvin / meter")], [678.1500000000001 * u("kelvin"), 51.41274238227152 * u("milliwatt / kelvin / meter")], [880.6500000000001 * u("kelvin"), 62.04986149584491 * u("milliwatt / kelvin / meter")], [1083.15 * u("kelvin"), 72.24376731301942 * u("milliwatt / kelvin / meter")], [1278.15 * u("kelvin"), 81.55124653739617 * u("milliwatt / kelvin / meter")], [1478.15 * u("kelvin"), 90.85872576177289 * u("milliwatt / kelvin / meter")], [1478.15 * u("kelvin"), 90.85872576177289 * u("milliwatt / kelvin / meter")], [1675.65 * u("kelvin"), 99.72299168975073 * u("milliwatt / kelvin / meter")], [1873.15 * u("kelvin"), 109.03047091412745 * u("milliwatt / kelvin / meter")]]],
        		[200.0 * u("standard_atmosphere"), #Line definition
        			[[90.65000000000009 * u("kelvin"), 159.11357340720224 * u("milliwatt / kelvin / meter")], [103.15000000000009 * u("kelvin"), 143.60110803324102 * u("milliwatt / kelvin / meter")], [135.6500000000001 * u("kelvin"), 98.83656509695294 * u("milliwatt / kelvin / meter")], [153.1500000000001 * u("kelvin"), 80.22160664819947 * u("milliwatt / kelvin / meter")], [180.6500000000001 * u("kelvin"), 59.83379501385045 * u("milliwatt / kelvin / meter")], [203.1500000000001 * u("kelvin"), 48.31024930747927 * u("milliwatt / kelvin / meter")], [225.6500000000001 * u("kelvin"), 42.54847645429368 * u("milliwatt / kelvin / meter")], [248.1500000000001 * u("kelvin"), 39.44598337950143 * u("milliwatt / kelvin / meter")], [275.6500000000001 * u("kelvin"), 37.67313019390585 * u("milliwatt / kelvin / meter")], [295.6500000000001 * u("kelvin"), 37.22991689750697 * u("milliwatt / kelvin / meter")], [325.6500000000001 * u("kelvin"), 37.67313019390585 * u("milliwatt / kelvin / meter")], [388.1500000000001 * u("kelvin"), 39.44598337950143 * u("milliwatt / kelvin / meter")], [453.1500000000001 * u("kelvin"), 42.54847645429368 * u("milliwatt / kelvin / meter")], [548.1500000000001 * u("kelvin"), 46.9806094182826 * u("milliwatt / kelvin / meter")], [675.6500000000001 * u("kelvin"), 53.18559556786707 * u("milliwatt / kelvin / meter")], [815.6500000000001 * u("kelvin"), 60.27700831024936 * u("milliwatt / kelvin / meter")], [1000.6500000000001 * u("kelvin"), 69.1412742382272 * u("milliwatt / kelvin / meter")]]],
        		[500.0 * u("standard_atmosphere"), #Line definition
        			[[105.65000000000009 * u("kelvin"), 159.11357340720224 * u("milliwatt / kelvin / meter")], [115.65000000000009 * u("kelvin"), 148.47645429362885 * u("milliwatt / kelvin / meter")], [128.1500000000001 * u("kelvin"), 136.95290858725764 * u("milliwatt / kelvin / meter")], [145.6500000000001 * u("kelvin"), 120.55401662049866 * u("milliwatt / kelvin / meter")], [165.6500000000001 * u("kelvin"), 105.92797783933523 * u("milliwatt / kelvin / meter")], [193.1500000000001 * u("kelvin"), 87.75623268698064 * u("milliwatt / kelvin / meter")], [220.6500000000001 * u("kelvin"), 75.34626038781167 * u("milliwatt / kelvin / meter")], [260.6500000000001 * u("kelvin"), 63.82271468144049 * u("milliwatt / kelvin / meter")], [288.1500000000001 * u("kelvin"), 59.39058171745157 * u("milliwatt / kelvin / meter")], [330.6500000000001 * u("kelvin"), 55.84487534626044 * u("milliwatt / kelvin / meter")], [385.6500000000001 * u("kelvin"), 54.07202216066486 * u("milliwatt / kelvin / meter")], [455.6500000000001 * u("kelvin"), 54.51523545706374 * u("milliwatt / kelvin / meter")], [510.6500000000001 * u("kelvin"), 55.40166204986153 * u("milliwatt / kelvin / meter")], [570.6500000000001 * u("kelvin"), 57.17451523545711 * u("milliwatt / kelvin / meter")], [663.1500000000001 * u("kelvin"), 60.27700831024936 * u("milliwatt / kelvin / meter")], [753.1500000000001 * u("kelvin"), 63.82271468144049 * u("milliwatt / kelvin / meter")], [843.1500000000001 * u("kelvin"), 67.36842105263162 * u("milliwatt / kelvin / meter")], [1003.1500000000001 * u("kelvin"), 74.016620498615 * u("milliwatt / kelvin / meter")]]],
        		[1000.0 * u("standard_atmosphere"), #Line definition
        			[[150.6500000000001 * u("kelvin"), 158.67036011080336 * u("milliwatt / kelvin / meter")], [188.1500000000001 * u("kelvin"), 139.612188365651 * u("milliwatt / kelvin / meter")], [228.1500000000001 * u("kelvin"), 118.78116343490308 * u("milliwatt / kelvin / meter")], [263.1500000000001 * u("kelvin"), 102.38227146814407 * u("milliwatt / kelvin / meter")], [300.6500000000001 * u("kelvin"), 88.64265927977843 * u("milliwatt / kelvin / meter")], [335.6500000000001 * u("kelvin"), 80.22160664819947 * u("milliwatt / kelvin / meter")], [388.1500000000001 * u("kelvin"), 74.45983379501388 * u("milliwatt / kelvin / meter")], [440.6500000000001 * u("kelvin"), 72.24376731301942 * u("milliwatt / kelvin / meter")], [483.1500000000001 * u("kelvin"), 71.80055401662054 * u("milliwatt / kelvin / meter")], [505.6500000000001 * u("kelvin"), 71.35734072022166 * u("milliwatt / kelvin / meter")], [533.1500000000001 * u("kelvin"), 70.47091412742387 * u("milliwatt / kelvin / meter")], [583.1500000000001 * u("kelvin"), 70.47091412742387 * u("milliwatt / kelvin / meter")], [663.1500000000001 * u("kelvin"), 72.24376731301942 * u("milliwatt / kelvin / meter")], [800.6500000000001 * u("kelvin"), 75.34626038781167 * u("milliwatt / kelvin / meter")], [898.1500000000001 * u("kelvin"), 78.44875346260392 * u("milliwatt / kelvin / meter")], [1073.15 * u("kelvin"), 84.21052631578951 * u("milliwatt / kelvin / meter")], [1275.65 * u("kelvin"), 91.74515235457068 * u("milliwatt / kelvin / meter")], [1438.15 * u("kelvin"), 97.95013850415515 * u("milliwatt / kelvin / meter")], [1620.65 * u("kelvin"), 105.04155124653744 * u("milliwatt / kelvin / meter")], [1620.65 * u("kelvin"), 105.04155124653744 * u("milliwatt / kelvin / meter")], [1863.15 * u("kelvin"), 114.79224376731307 * u("milliwatt / kelvin / meter")]]]]]]
        ]
        self.kGraphAnal = graphAnalysis(kGraphData)

        #Data taken from (6) Dynamic viscosity as a function of temperature and pressure
        #See GraphReader.ipynb to see how this data was generated
        muGraphData = [
        	["T","linear"], 	#x-axis definition
        	["mu","linear"], 	#y-axis definition
        	[[[1067,73.14999999999998 * u("kelvin")], [1850,1873.15 * u("kelvin")]], #x bounds
        	 [[263,200.0 * u("micropascal * second")], [738,0.0 * u("micropascal * second")]]],[ #y bounds

        	["P","log",[	#series definition
        		[1.0 * u("standard_atmosphere"), #Line definition
        			[[75.44885057471265 * u("kelvin"), 199.15789473684208 * u("micropascal * second")], [82.34540229885079 * u("kelvin"), 130.1052631578947 * u("micropascal * second")], [82.34540229885079 * u("kelvin"), 51.78947368421052 * u("micropascal * second")], [82.34540229885079 * u("kelvin"), 6.3157894736842195 * u("micropascal * second")], [82.34540229885079 * u("kelvin"), 6.3157894736842195 * u("micropascal * second")], [86.9431034482759 * u("kelvin"), 5.89473684210526 * u("micropascal * second")], [100.73620689655172 * u("kelvin"), 7.157894736842081 * u("micropascal * second")], [199.58678160919544 * u("kelvin"), 13.4736842105263 * u("micropascal * second")], [328.32241379310335 * u("kelvin"), 19.78947368421052 * u("micropascal * second")], [473.1500000000001 * u("kelvin"), 26.10526315789474 * u("micropascal * second")], [714.5293103448275 * u("kelvin"), 34.9473684210526 * u("micropascal * second")], [930.6212643678159 * u("kelvin"), 41.68421052631578 * u("micropascal * second")], [1222.575287356322 * u("kelvin"), 49.26315789473682 * u("micropascal * second")], [1480.0465517241378 * u("kelvin"), 55.99999999999997 * u("micropascal * second")], [1852.4603448275861 * u("kelvin"), 64.84210526315786 * u("micropascal * second")]]],
        		[5.0 * u("standard_atmosphere"), #Line definition
        			[[77.74770114942521 * u("kelvin"), 198.31578947368416 * u("micropascal * second")], [84.64425287356335 * u("kelvin"), 138.52631578947367 * u("micropascal * second")], [96.1385057471266 * u("kelvin"), 109.05263157894734 * u("micropascal * second")], [98.43735632183916 * u("kelvin"), 84.63157894736841 * u("micropascal * second")], [98.43735632183916 * u("kelvin"), 33.26315789473682 * u("micropascal * second")], [100.73620689655172 * u("kelvin"), 7.157894736842081 * u("micropascal * second")], [289.24195402298847 * u("kelvin"), 18.10526315789474 * u("micropascal * second")], [484.64425287356335 * u("kelvin"), 26.526315789473642 * u("micropascal * second")], [659.3568965517243 * u("kelvin"), 33.26315789473682 * u("micropascal * second")], [946.7132183908047 * u("kelvin"), 42.10526315789474 * u("micropascal * second")], [1222.575287356322 * u("kelvin"), 49.68421052631578 * u("micropascal * second")], [1549.0120689655173 * u("kelvin"), 57.68421052631578 * u("micropascal * second")], [1861.6557471264364 * u("kelvin"), 65.26315789473682 * u("micropascal * second")]]],
        		[10.0 * u("standard_atmosphere"), #Line definition
        			[[77.74770114942521 * u("kelvin"), 197.89473684210523 * u("micropascal * second")], [84.64425287356335 * u("kelvin"), 148.21052631578945 * u("micropascal * second")], [98.43735632183916 * u("kelvin"), 106.52631578947367 * u("micropascal * second")], [107.63275862068986 * u("kelvin"), 67.78947368421049 * u("micropascal * second")], [109.93160919540242 * u("kelvin"), 42.10526315789474 * u("micropascal * second")], [109.93160919540242 * u("kelvin"), 8.0 * u("micropascal * second")], [282.3454022988508 * u("kelvin"), 18.10526315789474 * u("micropascal * second")], [498.43735632183916 * u("kelvin"), 27.36842105263156 * u("micropascal * second")], [498.43735632183916 * u("kelvin"), 27.36842105263156 * u("micropascal * second")], [714.5293103448275 * u("kelvin"), 35.36842105263156 * u("micropascal * second")], [974.2994252873564 * u("kelvin"), 42.9473684210526 * u("micropascal * second")], [1153.6097701149424 * u("kelvin"), 48.0 * u("micropascal * second")], [1388.0925287356322 * u("kelvin"), 53.89473684210526 * u("micropascal * second")], [1629.4718390804596 * u("kelvin"), 59.78947368421049 * u("micropascal * second")], [1861.6557471264364 * u("kelvin"), 65.26315789473682 * u("micropascal * second")]]],
        		[20.0 * u("standard_atmosphere"), #Line definition
        			[[80.04655172413777 * u("kelvin"), 197.89473684210523 * u("micropascal * second")], [84.64425287356335 * u("kelvin"), 156.21052631578945 * u("micropascal * second")], [93.83965517241404 * u("kelvin"), 110.73684210526312 * u("micropascal * second")], [109.93160919540242 * u("kelvin"), 76.21052631578945 * u("micropascal * second")], [121.42586206896567 * u("kelvin"), 49.68421052631578 * u("micropascal * second")], [121.42586206896567 * u("kelvin"), 8.842105263157862 * u("micropascal * second")], [144.41436781609218 * u("kelvin"), 10.947368421052602 * u("micropascal * second")], [194.98908045977032 * u("kelvin"), 13.89473684210526 * u("micropascal * second")], [337.51781609195405 * u("kelvin"), 21.05263157894734 * u("micropascal * second")], [482.3454022988508 * u("kelvin"), 26.9473684210526 * u("micropascal * second")], [700.7362068965517 * u("kelvin"), 34.9473684210526 * u("micropascal * second")], [896.1385057471266 * u("kelvin"), 40.84210526315786 * u("micropascal * second")], [1160.5063218390806 * u("kelvin"), 48.0 * u("micropascal * second")], [1420.276436781609 * u("kelvin"), 54.73684210526312 * u("micropascal * second")], [1613.3798850574713 * u("kelvin"), 59.36842105263156 * u("micropascal * second")], [1861.6557471264364 * u("kelvin"), 65.26315789473682 * u("micropascal * second")]]],
        		[50.0 * u("standard_atmosphere"), #Line definition
        			[[77.74770114942521 * u("kelvin"), 197.89473684210523 * u("micropascal * second")], [84.64425287356335 * u("kelvin"), 152.42105263157893 * u("micropascal * second")], [93.83965517241404 * u("kelvin"), 119.15789473684208 * u("micropascal * second")], [103.03505747126428 * u("kelvin"), 87.99999999999997 * u("micropascal * second")], [116.8281609195401 * u("kelvin"), 63.99999999999997 * u("micropascal * second")], [123.72471264367823 * u("kelvin"), 50.9473684210526 * u("micropascal * second")], [123.72471264367823 * u("kelvin"), 50.9473684210526 * u("micropascal * second")], [132.92011494252893 * u("kelvin"), 34.52631578947364 * u("micropascal * second")], [132.92011494252893 * u("kelvin"), 34.52631578947364 * u("micropascal * second")], [139.8166666666666 * u("kelvin"), 19.78947368421052 * u("micropascal * second")], [146.71321839080474 * u("kelvin"), 16.421052631578902 * u("micropascal * second")], [158.207471264368 * u("kelvin"), 14.736842105263122 * u("micropascal * second")], [167.40287356321824 * u("kelvin"), 13.4736842105263 * u("micropascal * second")], [181.1959770114945 * u("kelvin"), 14.315789473684163 * u("micropascal * second")], [259.3568965517243 * u("kelvin"), 17.68421052631578 * u("micropascal * second")], [397.2879310344829 * u("kelvin"), 24.0 * u("micropascal * second")], [622.575287356322 * u("kelvin"), 32.4210526315789 * u("micropascal * second")], [840.9660919540229 * u("kelvin"), 39.15789473684208 * u("micropascal * second")], [1093.8396551724136 * u("kelvin"), 46.31578947368416 * u("micropascal * second")], [1093.8396551724136 * u("kelvin"), 46.31578947368416 * u("micropascal * second")], [1337.517816091954 * u("kelvin"), 52.63157894736838 * u("micropascal * second")], [1562.8051724137931 * u("kelvin"), 58.10526315789471 * u("micropascal * second")], [1863.954597701149 * u("kelvin"), 65.26315789473682 * u("micropascal * second")]]],
        		[100.0 * u("standard_atmosphere"), #Line definition
        			[[77.74770114942521 * u("kelvin"), 198.73684210526312 * u("micropascal * second")], [84.64425287356335 * u("kelvin"), 153.26315789473682 * u("micropascal * second")], [100.73620689655172 * u("kelvin"), 103.57894736842104 * u("micropascal * second")], [123.72471264367823 * u("kelvin"), 65.26315789473682 * u("micropascal * second")], [146.71321839080474 * u("kelvin"), 35.78947368421052 * u("micropascal * second")], [153.60977011494242 * u("kelvin"), 29.4736842105263 * u("micropascal * second")], [153.60977011494242 * u("kelvin"), 29.4736842105263 * u("micropascal * second")], [169.70172413793125 * u("kelvin"), 21.4736842105263 * u("micropascal * second")], [183.49482758620707 * u("kelvin"), 18.9473684210526 * u("micropascal * second")], [206.48333333333358 * u("kelvin"), 18.10526315789474 * u("micropascal * second")], [309.9316091954024 * u("kelvin"), 21.05263157894734 * u("micropascal * second")], [463.9545977011494 * u("kelvin"), 26.9473684210526 * u("micropascal * second")], [689.2419540229885 * u("kelvin"), 34.52631578947364 * u("micropascal * second")], [935.2189655172415 * u("kelvin"), 42.10526315789474 * u("micropascal * second")], [1199.5867816091954 * u("kelvin"), 49.26315789473682 * u("micropascal * second")], [1431.7706896551722 * u("kelvin"), 55.15789473684208 * u("micropascal * second")], [1613.3798850574713 * u("kelvin"), 59.36842105263156 * u("micropascal * second")], [1861.6557471264364 * u("kelvin"), 65.26315789473682 * u("micropascal * second")]]],
        		[500.0 * u("standard_atmosphere"), #Line definition
        			[[96.1385057471266 * u("kelvin"), 197.89473684210523 * u("micropascal * second")], [105.33390804597684 * u("kelvin"), 151.578947368421 * u("micropascal * second")], [126.02356321839079 * u("kelvin"), 110.31578947368419 * u("micropascal * second")], [158.207471264368 * u("kelvin"), 74.10526315789471 * u("micropascal * second")], [197.28793103448288 * u("kelvin"), 53.89473684210526 * u("micropascal * second")], [234.0695402298852 * u("kelvin"), 42.9473684210526 * u("micropascal * second")], [234.0695402298852 * u("kelvin"), 42.9473684210526 * u("micropascal * second")], [284.64425287356335 * u("kelvin"), 36.21052631578948 * u("micropascal * second")], [330.6212643678159 * u("kelvin"), 33.68421052631578 * u("micropascal * second")], [376.59827586206893 * u("kelvin"), 32.84210526315786 * u("micropascal * second")], [514.5293103448275 * u("kelvin"), 34.52631578947364 * u("micropascal * second")], [703.0350574712643 * u("kelvin"), 38.73684210526312 * u("micropascal * second")], [994.9890804597703 * u("kelvin"), 45.89473684210526 * u("micropascal * second")]]],
        		[1000.0 * u("standard_atmosphere"), #Line definition
        			[[153.60977011494242 * u("kelvin"), 198.31578947368416 * u("micropascal * second")], [238.66724137931033 * u("kelvin"), 111.57894736842104 * u("micropascal * second")], [298.43735632183916 * u("kelvin"), 55.15789473684208 * u("micropascal * second")], [307.63275862068986 * u("kelvin"), 50.52631578947364 * u("micropascal * second")], [332.9201149425289 * u("kelvin"), 45.89473684210526 * u("micropascal * second")], [376.59827586206893 * u("kelvin"), 42.9473684210526 * u("micropascal * second")], [406.4833333333336 * u("kelvin"), 42.10526315789474 * u("micropascal * second")], [638.6672413793103 * u("kelvin"), 42.9473684210526 * u("micropascal * second")], [746.7132183908047 * u("kelvin"), 44.21052631578948 * u("micropascal * second")], [921.4258620689657 * u("kelvin"), 47.57894736842104 * u("micropascal * second")], [1132.920114942529 * u("kelvin"), 51.78947368421052 * u("micropascal * second")], [1475.4488505747127 * u("kelvin"), 59.36842105263156 * u("micropascal * second")], [1475.4488505747127 * u("kelvin"), 59.36842105263156 * u("micropascal * second")], [1861.6557471264364 * u("kelvin"), 67.78947368421049 * u("micropascal * second")]]],
        		[5000.0 * u("standard_atmosphere"), #Line definition
        			[[307.63275862068986 * u("kelvin"), 179.36842105263156 * u("micropascal * second")], [438.6672413793103 * u("kelvin"), 136.84210526315786 * u("micropascal * second")], [512.230459770115 * u("kelvin"), 118.31578947368419 * u("micropascal * second")], [578.8971264367815 * u("kelvin"), 108.21052631578945 * u("micropascal * second")], [673.1500000000001 * u("kelvin"), 99.78947368421049 * u("micropascal * second")], [772.0005747126438 * u("kelvin"), 92.21052631578945 * u("micropascal * second")], [880.0465517241378 * u("kelvin"), 87.15789473684208 * u("micropascal * second")], [880.0465517241378 * u("kelvin"), 87.15789473684208 * u("micropascal * second")], [1015.6787356321838 * u("kelvin"), 82.52631578947367 * u("micropascal * second")], [1199.5867816091954 * u("kelvin"), 80.42105263157893 * u("micropascal * second")], [1436.3683908045978 * u("kelvin"), 79.15789473684208 * u("micropascal * second")], [1668.5522988505745 * u("kelvin"), 79.99999999999997 * u("micropascal * second")], [1863.954597701149 * u("kelvin"), 81.68421052631578 * u("micropascal * second")]]],
        		[10000.0 * u("standard_atmosphere"), #Line definition
        			[[567.4028735632182 * u("kelvin"), 197.89473684210523 * u("micropascal * second")], [838.6672413793103 * u("kelvin"), 155.7894736842105 * u("micropascal * second")], [1043.2649425287354 * u("kelvin"), 136.84210526315786 * u("micropascal * second")], [1330.621264367816 * u("kelvin"), 122.9473684210526 * u("micropascal * second")], [1491.540804597701 * u("kelvin"), 117.05263157894734 * u("micropascal * second")], [1574.2994252873564 * u("kelvin"), 114.9473684210526 * u("micropascal * second")], [1790.3913793103452 * u("kelvin"), 111.57894736842104 * u("micropascal * second")], [1866.2534482758624 * u("kelvin"), 110.73684210526312 * u("micropascal * second")]]]]]]
        ]
        self.muGraphAnal = graphAnalysis(muGraphData)

        #Data taken from (7) and (8) density as a function of Temperature and Pressure
        #See GraphReader.ipynb to see how this data was generated
        rhoGraphData = [
        	["T","linear"], 	#x-axis definition
        	["rho","linear"], 	#y-axis definition
        	[[[1047,173.14999999999998 * u("kelvin")], [1728,1773.15 * u("kelvin")]], #x bounds
        	 [[247,200.0 * u("kilogram / meter ** 3")], [630,0.0 * u("kilogram / meter ** 3")]]],[ #y bounds

        	["P","log",[	#series definition
                [1.0 * u("standard_atmosphere"), #Line definition
			        [[174.4797872340423 * u("kelvin"), 2.0073964497041423 * u("kilogram / meter ** 3")], [174.4797872340423 * u("kelvin"), 2.0073964497041423 * u("kilogram / meter ** 3")], [218.36276595744653 * u("kelvin"), 1.6168639053254439 * u("kilogram / meter ** 3")], [286.18191489361675 * u("kelvin"), 1.2322485207100593 * u("kilogram / meter ** 3")], [381.9265957446805 * u("kelvin"), 0.9245562130177518 * u("kilogram / meter ** 3")], [508.2563829787232 * u("kelvin"), 0.693786982248521 * u("kilogram / meter ** 3")], [508.2563829787232 * u("kelvin"), 0.693786982248521 * u("kilogram / meter ** 3")], [673.1499999999999 * u("kelvin"), 0.52810650887574 * u("kilogram / meter ** 3")], [876.6074468085103 * u("kelvin"), 0.4038461538461542 * u("kilogram / meter ** 3")], [1169.160638297872 * u("kelvin"), 0.303254437869823 * u("kilogram / meter ** 3")], [1169.160638297872 * u("kelvin"), 0.303254437869823 * u("kilogram / meter ** 3")]]],
        		[5.0 * u("standard_atmosphere"), #Line definition
        			[[182.5479441997063 * u("kelvin"), 9.921671018276754 * u("kilogram / meter ** 3")], [206.04280469897185 * u("kelvin"), 8.877284595300239 * u("kilogram / meter ** 3")], [260.0809838472833 * u("kelvin"), 6.788511749347208 * u("kilogram / meter ** 3")], [335.2645374449339 * u("kelvin"), 5.221932114882463 * u("kilogram / meter ** 3")], [431.5934654919238 * u("kelvin"), 4.177545691906005 * u("kilogram / meter ** 3")], [431.5934654919238 * u("kelvin"), 4.177545691906005 * u("kilogram / meter ** 3")], [532.6213656387663 * u("kelvin"), 3.6553524804177187 * u("kilogram / meter ** 3")], [589.0090308370045 * u("kelvin"), 3.1331592689294894 * u("kilogram / meter ** 3")]]],
        		[10.0 * u("standard_atmosphere"), #Line definition
        			[[180.19845814977953 * u("kelvin"), 19.84334203655351 * u("kilogram / meter ** 3")], [302.3717327459617 * u("kelvin"), 12.010443864229728 * u("kilogram / meter ** 3")], [518.524449339207 * u("kelvin"), 7.310704960835494 * u("kilogram / meter ** 3")], [852.1514684287813 * u("kelvin"), 4.699738903394234 * u("kilogram / meter ** 3")], [1488.862187958884 * u("kelvin"), 2.6109660574412032 * u("kilogram / meter ** 3")], [1895.3232745961823 * u("kelvin"), 2.6109660574412032 * u("kilogram / meter ** 3")]]],
        		[20.0 * u("standard_atmosphere"), #Line definition
        			[[184.8974302496331 * u("kelvin"), 40.73107049608353 * u("kilogram / meter ** 3")], [274.17790014684306 * u("kelvin"), 26.109660574412487 * u("kilogram / meter ** 3")], [274.17790014684306 * u("kelvin"), 26.109660574412487 * u("kilogram / meter ** 3")], [375.20580029368557 * u("kelvin"), 18.798955613576993 * u("kilogram / meter ** 3")], [485.63164464023475 * u("kelvin"), 14.621409921670988 * u("kilogram / meter ** 3")], [706.4833333333331 * u("kelvin"), 9.921671018276754 * u("kilogram / meter ** 3")], [706.4833333333331 * u("kelvin"), 9.921671018276754 * u("kilogram / meter ** 3")], [882.6947870778267 * u("kelvin"), 7.832898172323723 * u("kilogram / meter ** 3")], [993.1206314243759 * u("kelvin"), 7.310704960835494 * u("kilogram / meter ** 3")]]],
        		[50.0 * u("standard_atmosphere"), #Line definition
        			[[180.19845814977953 * u("kelvin"), 116.97127937336813 * u("kilogram / meter ** 3")], [217.79023494860485 * u("kelvin"), 86.68407310704958 * u("kilogram / meter ** 3")], [264.7799559471364 * u("kelvin"), 67.3629242819843 * u("kilogram / meter ** 3")], [325.8665932452277 * u("kelvin"), 53.785900783289776 * u("kilogram / meter ** 3")], [394.00168869309846 * u("kelvin"), 43.86422976501302 * u("kilogram / meter ** 3")], [480.93267254038165 * u("kelvin"), 35.50913838120101 * u("kilogram / meter ** 3")], [589.0090308370045 * u("kelvin"), 29.765013054830263 * u("kilogram / meter ** 3")], [753.4730543318651 * u("kelvin"), 22.976501305482998 * u("kilogram / meter ** 3")], [990.7711453744491 * u("kelvin"), 17.23237597911225 * u("kilogram / meter ** 3")], [1253.913582966226 * u("kelvin"), 13.577023498694473 * u("kilogram / meter ** 3")], [1601.6375183553596 * u("kelvin"), 10.966057441253213 * u("kilogram / meter ** 3")], [1895.3232745961823 * u("kelvin"), 9.921671018276754 * u("kilogram / meter ** 3")]]],
        		[100.0 * u("standard_atmosphere"), #Line definition
        			[[206.04280469897185 * u("kelvin"), 199.47780678851174 * u("kilogram / meter ** 3")], [206.04280469897185 * u("kelvin"), 199.47780678851174 * u("kilogram / meter ** 3")], [255.38201174743017 * u("kelvin"), 143.6031331592689 * u("kilogram / meter ** 3")], [283.5758443465493 * u("kelvin"), 125.32637075718014 * u("kilogram / meter ** 3")], [342.3129955947138 * u("kelvin"), 102.87206266318537 * u("kilogram / meter ** 3")], [342.3129955947138 * u("kelvin"), 102.87206266318537 * u("kilogram / meter ** 3")], [415.1470631424377 * u("kelvin"), 81.4621409921671 * u("kilogram / meter ** 3")], [415.1470631424377 * u("kelvin"), 81.4621409921671 * u("kilogram / meter ** 3")], [440.99140969163 * u("kelvin"), 75.71801566579632 * u("kilogram / meter ** 3")], [511.4759911894271 * u("kelvin"), 66.31853785900779 * u("kilogram / meter ** 3")], [511.4759911894271 * u("kelvin"), 66.31853785900779 * u("kilogram / meter ** 3")], [619.5523494860499 * u("kelvin"), 55.35248041775452 * u("kilogram / meter ** 3")], [720.5802496328929 * u("kelvin"), 46.99738903394251 * u("kilogram / meter ** 3")], [892.092731277533 * u("kelvin"), 38.12010443864227 * u("kilogram / meter ** 3")], [1101.1969897209983 * u("kelvin"), 31.331592689295007 * u("kilogram / meter ** 3")], [1329.0971365638766 * u("kelvin"), 26.109660574412487 * u("kilogram / meter ** 3")], [1566.395227606461 * u("kelvin"), 22.45430809399477 * u("kilogram / meter ** 3")], [1892.9737885462555 * u("kelvin"), 18.798955613576993 * u("kilogram / meter ** 3")], [1892.9737885462555 * u("kelvin"), 18.798955613576993 * u("kilogram / meter ** 3")]]],
        		[1000.0 * u("standard_atmosphere"), #Line definition
        			[[1401.9312041116004 * u("kelvin"), 199.47780678851174 * u("kilogram / meter ** 3")], [1401.9312041116004 * u("kelvin"), 199.47780678851174 * u("kilogram / meter ** 3")], [1479.4642437591779 * u("kelvin"), 190.07832898172322 * u("kilogram / meter ** 3")], [1479.4642437591779 * u("kelvin"), 190.07832898172322 * u("kilogram / meter ** 3")], [1592.2395741556534 * u("kelvin"), 178.59007832898172 * u("kilogram / meter ** 3")], [1592.2395741556534 * u("kelvin"), 178.59007832898172 * u("kilogram / meter ** 3")], [1716.7623348017623 * u("kelvin"), 168.1462140992167 * u("kilogram / meter ** 3")], [1716.7623348017623 * u("kelvin"), 168.1462140992167 * u("kilogram / meter ** 3")], [1810.7417767988254 * u("kelvin"), 161.35770234986944 * u("kilogram / meter ** 3")], [1890.6243024963287 * u("kelvin"), 154.56919060052218 * u("kilogram / meter ** 3")], [1890.6243024963287 * u("kelvin"), 154.56919060052218 * u("kilogram / meter ** 3")]]]]]]
        ]
        self.rhoGraphAnal = graphAnalysis(rhoGraphData)


    def Cp(self,T):
        return interp1D(T,self.CpPoints)
    def Cv(self,T):
        return interp1D(T,self.CvPoints)

    def k(self,T,P):
        return self.kGraphAnal(T=T,P=P).k


    def mu(self,T,P):
        return self.muGraphAnal(T=T,P=P).mu

    def rho(self,T,P):
        return self.rhoGraphAnal(T=T,P=P).rho

    def nu(self,T,P):
        return self.mu(T=T,P=P)/self.rho(T=T,P=P)

    def alpha(self,T,P):
        return self.k(T=T,P=P)/self.rho(T=T,P=P)*self.MW/self.Cp(T=T,P=P) # m^2/s

    def Pr(self,T,P): # Prandtl Number
        return ( (self.nu(T=T,P=P)/self.alpha(T=T,P=P)) ).to(u.diml).magnitude
air = Air()

class transientBalance():
    """
    A functor to solve transient heat transfer problems.

    To initiate, specify the parameters:
        regime : The spatial configuration which you'd like to generate the solution for. ("Plane Wall", "Cylinder" ,"Sphere") (defaluts to "Plane Wall")
        method : The calculation method which you'd like to use to solve. ("Analytical", "Lumped Capacitance") (defaults to "Analytical")
            Note that the "Lumped Capacitance" method is still under development.

    For example:
        tansBal = transientBalance(regime="Plane Wall", method="Analytical")

    For documentation on the theory behind these functions, see "Fundamentals of Heat and Mass Transfer, T.L. Bergman", 8th Edition, Chapter 5

    The execution of the functor returns the value for Theta* (thetaStar) given the following parameters for the conditions of interest:
        Fo       : Fourier number
        biot     : Biot number
        dStar    : Dimensionless distance number (called r-Star in the book)(Optional for Lumped Capacitance solutions)
        numTerms : The number of terms in the summation for the theta-Star equation (Optional, defaults to 10 for Analytical solutions and is not needed for Lumped Capacitance solutions)
        zetas    : Known values for zeta given in the order they would show up in within the summation of the theta-Star equation (Optional, If you know the zeta values, insert them here and it will save some computational time. Otherwise they'll be automatically calculated.  They are also not needed for Lumped Capacitance solutions)

    For example:
        tStar = transBal(Fo = 0.2, biot = 0.6, dStar = .80, numTerms = 100, zetas = [(Insert known zeta values here)])
    """

    # Functions to be solved by the zeta-solver for each geometry
    def wallFunc(self,x,Bi):
        return x*np.tan(x)-Bi

    def cylFunc(self,x,Bi):
        return x*bessel.j1(x)/bessel.j0(x)-Bi

    def sphFunc(self,x,Bi):
        return 1.0-x/np.tan(x)-Bi

    # The Full zeta functions for each geometry
    def wallZeta(self,n,Bi):
        if n == 1:
            lo = 0
            hi = (2 * n - 1) * np.pi / 2
            return ridder(self.wallFunc,lo,hi,args=(Bi,))
        else:
            lo = (2 * n - 1) * np.pi / 2 - np.pi * 0.99999
            hi = (2 * n - 1) * np.pi / 2 * 0.99999
            return ridder(self.wallFunc,lo,hi,args=(Bi,))

    def cylZeta(self,n,Bi):
        if n == 1:
            lo = 0
            hi = 2.4
            return ridder(self.cylFunc,lo,hi,args=(Bi,))
        else:
            lo = 0.005 + (n - 1) * np.pi
            hi = 2.35 + (n - 1) * np.pi
            return ridder(self.cylFunc,lo,hi,args=(Bi,))

    def sphZeta(self,n,Bi):
        lo = n * np.pi - 0.99999 * np.pi
        hi = n * np.pi * 0.99999
        return ridder(self.sphFunc,lo,hi,args=(Bi,))

    # C functions for each geometry
    def wallC(self,z):
    	return 4 * np.sin(z) / (2 * z + np.sin(2 * z))

    def cylC(self,z):
    	return 2 * bessel.j1(z) / (z * (bessel.j0(z)**2 + bessel.j1(z)**2))

    def sphC(self,z):
    	return 4 * (np.sin(z) - z * np.cos(z))/(2 * z - np.sin(2 * z))


    def t_Star_i_wal(Fo,xStar,Ci,zi):
        return Ci * np.exp(-(zi**2) * Fo) * np.cos(zi * xStar)
    def t_Star_i_cyl(Fo,rStar,Ci,zi):
        return Ci * np.exp(-(zi**2) * Fo) * bessel.j0(zi * rStar)
    def t_Star_i_sph(Fo,rStar,Ci,zi):
        return Ci * np.exp(-(zi**2) * Fo) * np.sin(zi * rStar) / (zi * rStar)

    zfunc = "NOT DEFINED"
    Cfunc = "NOT DEFINED"
    thetaSi = "NOT DEFINED"
    def thetaStar_anal(self,Fo,biot,dStar,numTerms=10,zetas="NotSpecified"):
        ans = 0 * u.diml

        for n in range(1,numTerms+1):
            #print("#########",n)
            if zetas == "NotSpecified":
                z = self.zfunc(self,n,biot)
            else:
                z = zetas[n]
            #print("z",z)
            c = self.Cfunc(z,biot)
            #print("c",c)
            ans += self.thetaSi(Fo,dStar,c,z)
            #print("dt*",self.thetaSi(Fo,dStar,c,z))
            #print("t*",ans)
        return ans
    def thetaStar_LC(self,Fo,biot,dStar="NotNeeded",numTerms="NotNeeded"):
        return "THIS FUNCTION HAS NOT BEEN WRITTEN YET."

    regimes  = ["Plane Wall", "Cylinder"  , "Sphere"    ]
    zs       = [wallZeta    , cylZeta     , sphZeta     ]
    Cs       = [wallC       , cylC        , sphC        ]
    thetaSis = [t_Star_i_wal, t_Star_i_cyl, t_Star_i_sph]

    methods     = ["Analytical"   ,"Lumped Capacitance"]
    methodfuncs = [thetaStar_anal ,thetaStar_LC        ]

    def __init__(self,regime="Plane Wall",method="Analytical"):
        regimeFound = False
        methodFound = False

        for i in range(len(self.regimes)):
            if regime == self.regimes[i]:
                self.zfunc  = self.zs[i]
                self.Cfunc  = self.Cs[i]
                self.thetaSi = self.thetaSis[i]
                regimeFound = True
        if not regimeFound:
            print("ERROR! Regime " + regime + " not found.")
            print("Valid regimes are:")
            for r in self.regimes:
                print(r)
            raise Exception("Regime Not Found")

        for i in range(len(self.methods)):
            if method == self.methods[i]:
                def thetaStar(Fo,biot,dStar="NotGiven",numTerms="NotGiven",zetas="NotSpecified"):
                    return self.methodfuncs[i](self,Fo,biot,dStar,numTerms,zetas)
                self.thetaStar = thetaStar
                methodFound = True
                break
        if not methodFound:
            print("ERROR! Method " + method + " not found.")
            print("Valid methods are:")
            for m in self.methods:
                print(m)
            raise Exception("Method Not Found")
    def __call__(self,Fo,biot,dStar="NotGiven",numTerms=10,zetas="NotSpecified"):
        return self.thetaStar(Fo,biot,dStar,numTerms,zetas)

class Correlation:
    """
    A function wrapper (functor) to hold a corrleation function and to raise a warning if any of the given conditions fall outside their specified ranges.

    To wrap a function (i.e. to create a Correlation functor), instantiate a Correlation object with the following parameters:
        funcToWrap           : The python function you'd like to wrap
        paramRanges          : A dictionary that maps the name of each paremeter to a list representing it's upper and lower bound ([lower,upper]).  The key must be a string matching the parameter name as it is pass into the funcToWrap
                               Any parameter not included in this dictionary will be passed straight into the funcToWrap without being checked.
        additionalConditions : A string with any additional conditions that cannot be mathematically represted by the bounds given in paramRanges.  This will be printed as part of the printConditions() function. (optional)
        allParamsDiml        : A Boolean that, if true, indicates that all parameters are dimensionless and should be converted to a typical flost (non-Pint object) before being passed into the function. (optional)
                               This is particularly usefull for correlations in which the parameters are raised to non-integer powers, because non-integer powers often to not work well with Pint.

    Once instantiated, the Correlation object works the same as the funcToWrap, but it will check to make sure all the parameters fall within the acceptable range before executing.

    Other methods:
        printConditions : Call this function to print all the specified conditions for this correlation.
        displayWarnings : A boolean which, if True, will print warnings if any of the parameters fall outside the range required by their bounds.  A value of False will disable this checking.
    """

    def __init__(self,funcToWrap,paramRanges,additionalConditions="",allParamsDiml=False):
        #First, check to make sure each param specified in paramRanges is actually a param of the funcToWrap

        rangeParams = paramRanges.keys()
        sig = signature(funcToWrap)

        self.funcParamList = list(sig.parameters)
        funcParams = set()
        for funcParam in self.funcParamList:
            funcParams.add(funcParam)

        for rangeParam in rangeParams:
            if not rangeParam in funcParams:
                raise Exception("Error! \"" + str(rangeParam) + "\" is not a parameter of \"" + str(funcToWrap.__name__) + "\".")


        self.funcToWrap = funcToWrap
        self.paramRanges = paramRanges
        self.additionlConditions = additionalConditions
        self.allParamsDiml = allParamsDiml
        self.displayWarnings = True

    def printWarning(self,varName,val):
        if self.displayWarnings:
            message = "\n" + str(self.funcToWrap.__name__) + " is not meant to handle a \"" + str(varName) + "\" of \"" + str(val) + "\". \"" + str(varName) + "\" must fall between " + str(self.paramRanges[varName][0]) + " and " + str(self.paramRanges[varName][1]) + "."
            message += "\nTo disable these warnings, set the \"displayWarnings\" attribute of this Correlation to False."
            warnings.warn(message,Warning)

    def __call__(self,*args,**kwargs):
        #First, make sure each parameter is only defined once.
        for i in range(len(args)):
            argName = self.funcParamList[i]
            if argName in kwargs.keys():
                raise Exception("Error! Multiple definitions for \"" + str(argName) + "\"!")
            else:
                kwargs[argName] = args[i]

        #Check to make sure each parameter fits within it's range, and turn it into a dimensionless float if that option is indicated.
        for kw in kwargs:
            if kw in self.paramRanges.keys():
                lo, hi = self.paramRanges[kw]

                if (kwargs[kw] < lo) or (kwargs[kw] > hi):
                    self.printWarning(kw,kwargs[kw])
            if self.allParamsDiml:
                if not isinstance(kwargs[ks],(float,int)):
                    kwargs[kw] = kwargs[kw].to(u.diml).magnitude

        #Run the function.
        return self.funcToWrap(**kwargs)

    def printConditions(self):
        print("The conditions for \"" + str(self.funcToWrap.__name__) + "\" are:")
        print(self.additionalConditions)
        for param in self.paramRanges:
            print(param,":",self.paramRanges[param])


#------------------ MATERIAL SCIENCE CONSTANTS AND FUNCTIONS ------------------#

class Element:
    """
    A class to house the properties of elements.

    Most data is pulled directly from the "periodictable" python module.

    Though additional data can be added later. (See addDataset() in the ElementList class)

    To reference a property, just is ".___" where "___" is the name of the property you want.
    """
    def __init__(self,**kwargs):
        conversionTable = { #Some name changes, to accomidate for multiple names for the same value.
            "number" : [{"AN"},None],
            "_isotopes" : [{"isotopes"},None],
            "_mass" : [{"mass","AW"},u.g / u.mol],
            "_density" : [{"density","rho","ρ"},u.g / u.cm**3],
        }

        for kw in kwargs:
            if kw in conversionTable:
                alternateNames,unit = conversionTable[kw]
                for altName in alternateNames:
                    if unit != None and kwargs[kw] != None:
                        self.__dict__.update({altName:kwargs[kw] * unit})
                    else:
                        self.__dict__.update({altName:kwargs[kw]})
                if unit != None and kwargs[kw] != None:
                    self.__dict__.update({kw:kwargs[kw] * unit})
                else:
                    self.__dict__.update({kw:kwargs[kw]})
            else:
                self.__dict__.update({kw:kwargs[kw]})

class ElementList:
    """
    A simple class to house all the Element objects in one place.

    To intantiate, just use Elements = ElementList()

    To refence an element, just use ".".  For example:
        Elements.H
        Elements.Al
        etc...
    """

    def __init__(self):
        for attribute in pTable.__dict__:
            try:
                if "number" in pTable.__dict__[attribute].__dict__:
                    #Only the actual elements have "number" attributes
                    kwArgs = pTable.__dict__[attribute].__dict__
                    elementObj = Element(**kwArgs)
                    self.__dict__.update({attribute:elementObj})
            except:
                pass
    def addDataset(self,attributeName,dataSet):
        """
        A function to add any other datasets to the Element Objects.

        attributeName : The name of the Element attribute that you'll be referencing this data by.
        For example:
            If I were adding a dataset containing crystal structure types, and I wanted to call these crystal structures by the attribute "struct" (Elements.Al.struct), I'd put "struct" as the attributeName.

        dataSet       : A python Dictionary mapping the element name and/or symbol to the peice of data for that element.
            Note that elements are stored in the ElementList by both element name and symbol.  But, unfortunatley, referencing an element by it's name or sybmol gives you two different objects.  So if you add data to one (e.g. the element name), I'd recomend also adding the same data to the other (e.g. the element symbol)
        """

        for key in dataSet:
            if key not in self.__dict__.keys():
                raise Exception("Error! \"" + str(key) + "\" is not a recognized element.")

            self.__dict__[key].__dict__[attributeName] = dataSet[key]



Elements = ElementList()

#Data from https://en.wikipedia.org/wiki/Periodic_table_(crystal_structure)
structureData = {'hydrogen': ['HEX'], 'H': ['HEX'], 'helium': ['HCP'], 'He': ['HCP'], 'lithium': ['BCC'], 'Li': ['BCC'], 'beryllium': ['HCP'], 'Be': ['HCP'], 'boron': ['RHO'], 'B': ['RHO'], 'nitrogen': ['HEX'], 'N': ['HEX'], 'oxygen': ['SC'], 'O': ['SC'], 'fluorine': ['SC'], 'F': ['SC'], 'neon': ['FCC'], 'Ne': ['FCC'], 'sodium': ['BCC'], 'Na': ['BCC'], 'magnesium': ['HCP'], 'Mg': ['HCP'], 'aluminum': ['FCC'], 'Al': ['FCC'], 'silicon': ['DC'], 'Si': ['DC'], 'phosphorus': ['ORTH'], 'P': ['ORTH'], 'sulfur': ['ORTH'], 'S': ['ORTH'], 'chlorine': ['ORTH'], 'Cl': ['ORTH'], 'argon': ['FCC'], 'Ar': ['FCC'], 'potassium': ['BCC'], 'K': ['BCC'], 'calcium': ['FCC'], 'Ca': ['FCC'], 'scandium': ['HCP'], 'Sc': ['HCP'], 'titanium': ['HCP'], 'Ti': ['HCP'], 'vanadium': ['BCC'], 'V': ['BCC'], 'chromium': ['BCC'], 'Cr': ['BCC'], 'manganese': ['BCC'], 'Mn': ['BCC'], 'iron': ['BCC'], 'Fe': ['BCC'], 'cobalt': ['HCP'], 'Co': ['HCP'], 'nickel': ['FCC'], 'Ni': ['FCC'], 'copper': ['FCC'], 'Cu': ['FCC'], 'zinc': ['HCP'], 'Zn': ['HCP'], 'gallium': ['ORTH'], 'Ga': ['ORTH'], 'germanium': ['DC'], 'Ge': ['DC'], 'arsenic': ['RHO'], 'As': ['RHO'], 'selenium': ['HEX'], 'Se': ['HEX'], 'bromine': ['ORTH'], 'Br': ['ORTH'], 'krypton': ['FCC'], 'Kr': ['FCC'], 'rubidium': ['BCC'], 'Rb': ['BCC'], 'strontium': ['FCC'], 'Sr': ['FCC'], 'yttrium': ['HCP'], 'Y': ['HCP'], 'zirconium': ['HCP'], 'Zr': ['HCP'], 'niobium': ['BCC'], 'Nb': ['BCC'], 'molybdenum': ['BCC'], 'Mo': ['BCC'], 'technetium': ['HCP'], 'Tc': ['HCP'], 'ruthenium': ['HCP'], 'Ru': ['HCP'], 'rhodium': ['FCC'], 'Rh': ['FCC'], 'palladium': ['FCC'], 'Pd': ['FCC'], 'silver': ['FCC'], 'Ag': ['FCC'], 'cadmium': ['HCP'], 'Cd': ['HCP'], 'indium': ['TETR'], 'In': ['TETR'], 'tin': ['TETR'], 'Sn': ['TETR'], 'antimony': ['RHO'], 'Sb': ['RHO'], 'tellurium': ['HEX'], 'Te': ['HEX'], 'iodine': ['ORTH'], 'I': ['ORTH'], 'xenon': ['FCC'], 'Xe': ['FCC'], 'cesium': ['BCC'], 'Cs': ['BCC'], 'barium': ['BCC'], 'Ba': ['BCC'], 'lanthanum': ['DHCP'], 'La': ['DHCP'], 'hafnium': ['HCP'], 'Hf': ['HCP'], 'tantalum': ['TETR', 'BCC'], 'Ta': ['TETR', 'BCC'], 'tungsten': ['BCC'], 'W': ['BCC'], 'rhenium': ['HCP'], 'Re': ['HCP'], 'osmium': ['HCP'], 'Os': ['HCP'], 'iridium': ['FCC'], 'Ir': ['FCC'], 'platinum': ['FCC'], 'Pt': ['FCC'], 'gold': ['FCC'], 'Au': ['FCC'], 'mercury': ['RHO'], 'Hg': ['RHO'], 'thallium': ['HCP'], 'Tl': ['HCP'], 'lead': ['FCC'], 'Pb': ['FCC'], 'bismuth': ['RHO'], 'Bi': ['RHO'], 'polonium': ['SC', 'RHO'], 'Po': ['SC', 'RHO'], 'astatine': ['FCC'], 'At': ['FCC'], 'radon': ['FCC'], 'Rn': ['FCC'], 'francium': ['BCC'], 'Fr': ['BCC'], 'radium': ['BCC'], 'Ra': ['BCC'], 'actinium': ['FCC'], 'Ac': ['FCC'], 'rutherfordium': ['HCP'], 'Rf': ['HCP'], 'dubnium': ['BCC'], 'Db': ['BCC'], 'seaborgium': ['BCC'], 'Sg': ['BCC'], 'bohrium': ['HCP'], 'Bh': ['HCP'], 'hassium': ['HCP'], 'Hs': ['HCP'], 'meitnerium': ['FCC'], 'Mt': ['FCC'], 'darmstadtium': ['BCC'], 'Ds': ['BCC'], 'roentgenium': ['BCC'], 'Rg': ['BCC'], 'copernicium': ['HCP'], 'Cn': ['HCP'], 'nihonium': ['HCP'], 'Nh': ['HCP'], 'flerovium': ['FCC'], 'Fl': ['FCC'], 'oganesson': ['FCC'], 'Og': ['FCC'], 'cerium': ['DHCP', 'FCC'], 'Ce': ['DHCP', 'FCC'], 'praseodymium': ['DHCP'], 'Pr': ['DHCP'], 'neodymium': ['DHCP'], 'Nd': ['DHCP'], 'promethium': ['DHCP'], 'Pm': ['DHCP'], 'samarium': ['RHO'], 'Sm': ['RHO'], 'europium': ['BCC'], 'Eu': ['BCC'], 'gadolinium': ['HCP'], 'Gd': ['HCP'], 'terbium': ['HCP'], 'Tb': ['HCP'], 'dysprosium': ['HCP'], 'Dy': ['HCP'], 'holmium': ['HCP'], 'Ho': ['HCP'], 'erbium': ['HCP'], 'Er': ['HCP'], 'thulium': ['HCP'], 'Tm': ['HCP'], 'ytterbium': ['FCC'], 'Yb': ['FCC'], 'lutetium': ['HCP'], 'Lu': ['HCP'], 'thorium': ['FCC'], 'Th': ['FCC'], 'protactinium': ['TETR'], 'Pa': ['TETR'], 'uranium': ['ORTH'], 'U': ['ORTH'], 'neptunium': ['ORTH'], 'Np': ['ORTH'], 'plutonium': ['MON'], 'Pu': ['MON'], 'americium': ['DHCP'], 'Am': ['DHCP'], 'curium': ['DHCP'], 'Cm': ['DHCP'], 'berkelium': ['DHCP'], 'Bk': ['DHCP'], 'californium': ['DHCP'], 'Cf': ['DHCP'], 'einsteinium': ['FCC'], 'Es': ['FCC'], 'fermium': ['FCC'], 'Fm': ['FCC'], 'mendelevium': ['FCC'], 'Md': ['FCC'], 'nobelium': ['FCC'], 'No': ['FCC'], 'lawrencium': ['HCP'], 'Lr': ['HCP']}

Elements.addDataset("structs",structureData)

#--------------------------- SEPERATIONS FUNCTIONS ----------------------------#

def mccabeThieleStepper(strippingLine,rectifyingLine,eqLine,initVal,finalVal,direction,graph=True,E_MV=1,limitStages=None,Legend=False):
    """
    A function to graph and analyze a McCabe Thiele diagram.

    strippingLine  : A function that takes xi and returns the yi of the stripping line
    rectifyingLine : A function that takes xi and returns the yi of the rectifying line.
    eqLine         : A function that takes xi and returns the yi of the equilibrium line.
    initVal        : A float for the initial xi or yi.  It will be assigned to an [xi,yi] coordinate
    finalVal       : A float for the final xi or yi.  It will be assigned to an [xi,yi] coordinate
    direction      : "+" if you'd like to step in the increaseing xi direction, "-" if you'd like to step in the decreasing xi direction.
    graph          : True if you'd like to display the graph
    printNum       : True if you'd like to print the number of stages required
    E_MV           : Murphry Efficiency for the Vapor
    E_ML           : Murphry Efficienct for the Liquid

    returns        : A list of 1) the stripping and rectifying line intersection, 2) the number of stages needed, 3) Each point on the stepping Graph

    Note that, at the moment, this function can only handle regular distilation (Vapor/Liquid distilation) problems, NOT absorption or stripping prolems.
    """

    returnObjects = []

    def strpInv(yi):
        """Returns the xi of the strippingLine, given the yi"""
        def solveMe(xi):
            e = strippingLine(xi) - yi
            try:
                e = e.magnitude
            except:
                pass
            return e
        xi = fsolve(solveMe,yi)[0]
        return xi
    def rectInv(yi):
        """Returns the xi of the rectifyingLine, given the yi"""
        def solveMe(xi):
            e = rectifyingLine(xi) - yi
            try:
                e = e.magnitude
            except:
                pass
            return e
        xi = fsolve(solveMe,yi)[0]
        return xi
    def eqInv(yi):
        """Returns the xi of the eqLine, given the yi"""
        def solveMe(xi):
            e = eqLine(xi) - yi
            try:
                e = e.magnitude
            except:
                pass
            return e
        xi = fsolve(solveMe,yi)[0]
        return xi

    #Calculate the intersection of the stripping and recifying lines
    def solveIntersect(xi):
        L = strippingLine(xi)
        R = rectifyingLine(xi)
        e = L - R
        try:
            e = e.magnitude
        except:
            pass
        return e
    intersectX = fsolve(solveIntersect,0.5)[0]
    intersectY = strippingLine(intersectX)

    intersect = [intersectX,intersectY]
    returnObjects.append(intersect)

    #Determine which side of the graph is dominated by which function
    yRectZeroTest = rectifyingLine(0)
    yStrpZeroTest = strippingLine(0)

    if yRectZeroTest < yStrpZeroTest:
        leftDominating = rectifyingLine
        rightDominating = strippingLine
        leftDominatingInv = rectInv
        rightCominatingInv = strpInv
        r = "Rectifying Line"
        l = "Stripping Line"
    else:
        rightDominating = rectifyingLine
        leftDominating = strippingLine
        rightDominatingInv = rectInv
        leftDominatingInv = strpInv
        l = "Rectifying Line"
        r = "Stripping Line"

    def eqLineWITHEFF(xA):
        if xA < intersect[0]:
            reigionFunc = leftDominating
        else:
            reigionFunc = rightDominating

        eqOriginal = eqLine(xA)
        rfunc = reigionFunc(xA)
        dy = eqOriginal - rfunc
        dy *= E_MV
        return rfunc + dy
    def eqInvWITHEFF(yi):
        """Returns the xi of the eqLine, given the yi"""
        def solveMe(xi):
            e = eqLineWITHEFF(xi) - yi
            try:
                e = e.magnitude
            except:
                pass
            return e
        xi = fsolve(solveMe,yi)[0]
        return xi

    if graph:
        xs1 = np.linspace(0,intersectX,100) * u.diml
        xs = []
        ys = []
        for x in xs1:
            y = leftDominating(x)
            if y >= x:
                ys.append(y)
                xs.append(x)
        ys = np.array(ys) * u.diml
        xs = np.array(xs) * u.diml
        if xs.size > 0 and ys.size > 0:
            pintPlot(xs,ys,label=l,legend=Legend)
        xs1 = np.linspace(intersectX,1,100) * u.diml
        xs = []
        ys = []
        for x in xs1:
            y = rightDominating(x)
            if y >= x:
                ys.append(y)
                xs.append(x)
        ys = np.array(ys) * u.diml
        xs = np.array(xs) * u.diml
        if xs.size > 0 and ys.size > 0:
            pintPlot(xs,ys,label=r,legend=Legend)
        xs1 = np.linspace(0,1,100)
        plt.plot(xs1,xs1,",")
        xs1 *= u.diml
        ysEQ = []
        xs = []
        for x in xs1:
            try:
                ysEQ.append(eqLine(x) * u.diml)
                xs.append(x)
            except:
                pass
        pintPlot(xs,ysEQ,label="Equilibrium Line",legend=Legend)



    points = [[initVal,initVal]]

    if direction == "+" :
        while True:
            # First trace up to the EQ line
            yEQ = eqLineWITHEFF(points[-1][0])

            xEQ = points[-1][0]

            points.append([xEQ,yEQ])

            if yEQ < intersect[1]:
                xnew = leftDominatingInv(yEQ)
            else:
                xnew = rightDominatingInv(yEQ)

            points.append([xnew,yEQ])

            if xnew > finalVal:
                break

            if limitStages != None:
                numStages = (len(points) - 1) / 2
                if numStages >= limitStages:
                    break


    elif direction == "-":
         while True:
            # First trace over to the EQ line
            yEQ = points[-1][1]
            xEQ = eqInvWITHEFF(yEQ)

            points.append([xEQ,yEQ])

            if xEQ < intersect[0]:
                ynew = leftDominating(xEQ)
            else:
                ynew = rightDominating(xEQ)

            points.append([xEQ,ynew])

            if ynew < finalVal:
                break

            if limitStages != None:
                numStages = (len(points) - 1) / 2
                if numStages >= limitStages:
                    break
    else:
        raise Exception("\"" + str(direction) + "\" is not a recognized direction.  Please use \"+\" or \"-\"")

    if graph:
        xs = []
        ys = []
        for point in points:
            xs.append(point[0])
            ys.append(point[1])

        xs = np.array(xs) * u.diml
        ys = np.array(ys) * u.diml

        pintPlot(xs,ys,xlabel="x",ylabel="y",label="Stages",legend=Legend)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.grid()
    numStages = (len(points) - 1) / 2
    returnObjects.append(numStages)
    returnObjects.append(points)
    return returnObjects
