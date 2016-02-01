# The MIT License (MIT)

# Copyright (c) 2016 Breschine Cummins

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
"""
import matplotlib
font = {'family' : 'normal',
        'size'   : 22}
matplotlib.rc('font', **font)
"""
import pdb

class hillmodel(object):
    '''
    This class takes a network file, a parameter sample file, and a Hill
    exponent and builds a Hill function model. The class has two public
    methods:

    1) time,timeseries = hillmodel.simulateHillModel(initialconditions,initialtime,finaltime,timestep)
    2) hillmodel.plotResults(times,timeseries)

    The first method generates a time series for a given set of initial conditions,
    and the second method plots the results. 

    '''
    def __init__(self,networkfile,samplefile,hillexp):
        '''
        Construct the Hill model for a given network and parameter sample.

        '''
        eqnstr,self.varnames = self._parseEqns(networkfile)
        parameternames,samples = self._parseSamples2(self.varnames,samplefile)
        self.eqns=self._makeHillEqns(eqnstr,parameternames,samples,hillexp)

    def simulateHillModel(self,initialconditions,initialtime,finaltime,timestep):
        '''
        Simulate the constructed Hill model for a given set of initial conditions 
        and time period. The given time step only specifies which output timeseries
        is returned. The time step for the backwards difference ODE solver is 
        determined by the algorithm.

        '''
        def RHS(t,x,eqns):
            xdot = [e(x) for e in eqns]
            return np.array(xdot)
        def integrate(r,y0,t0,t1,dt):
            times=[t0]
            timeseries=[y0]
            while r.successful() and r.t < t1:
                r.integrate(r.t+dt)
                times.append(r.t)
                timeseries.append(r.y)
            return times,timeseries
        r = ode(RHS).set_integrator('vode', method='bdf')
        r.set_initial_value(initialconditions,initialtime).set_f_params(self.eqns)
        times,timeseries = integrate(r,initialconditions,initialtime,finaltime,timestep)
        return times,timeseries

    def plotResults(self,times,timeseries,plotoptions={},legendoptions={},figuresize=()):
        '''
        Plot a time series.

        plotoptions and legendoptions are optional dictionaries with keys corresponding 
        to the options for matplotlib.pyplot.plot and matplotlib.pyplot.legend.

        Examples: 
        plotoptions={'linewidth':2}
        legendoptions={'fontsize':24,'loc':'upper left', 'bbox_to_anchor':(1, 1)}
        figuresize = (20,10)

        '''
        if figuresize:
            plt.figure(figsize=figuresize)
        timeseries=np.array(timeseries)
        for k in range(timeseries.shape[1]):
            plt.plot(times,timeseries[:,k],label=self.varnames[k],**plotoptions)
        plt.legend(**legendoptions)
        plt.axis('tight')
        plt.show()

    def _parseEqns(self,fname='equations.txt'):
        # Private parser.
        f=open(fname,'r')
        varnames=[]
        eqns0=[]
        for l in f:
            L=l.split(' : ')
            varnames.append(L[0])
	    #the following line (only 1 line) is for ubuntu only b/c appends X2\n; in windows, X2\n auto interpreted as X2
            equation = L[1][:-1]        
	    equation = equation.replace(')('," * ") #IE) 3D_Cycle's 224.txt has (X1)(~X2) as an eqn
            eqns0.append(equation)
        f.close()
	eqns = []
	for equation in eqns0: #this is done after to prevent mixing up eqns like (1 n) w/ eqns like (X1)(~X2) 
            equation = equation.replace(')',"").replace('(',"") #IE) 3D_Clock's 55.txt has (1 n) as an eqn    
            eqns.append(equation) 
        eqnstr=[]
        for e in eqns:
            for k,v in enumerate(varnames):
                e=e.replace('~'+v,str(k)+' n').replace(v,str(k)+' p').replace(')(',')*(')
            eqnstr.append(e)
        return eqnstr,varnames
        
    def _parseSamples2(self,varnames,fname='samples.txt'):
        # Private parser.
        f=open(fname,'r')
        pnames=[]
        samples_frac=[]
        lines = f.readlines()
        K=lines[0].split()
        #K=lines[1].split()
        prefix = ''
        for word in K[:-1]:
            if word[:-1].isdigit() == True or '/' in word:
                samples_frac.append(word[:-1])
            elif '[' in word: #L[X,
                prefix = word
            elif ']' in word: #Y]
                pnames.append(prefix + word)
        samples_frac.append(K[-1])
        parameternames=[]
        samples=[]
        for p in pnames:
            for k,v in enumerate(varnames):
                p=p.replace(v,str(k))
            parameternames.append(p)
        for value in samples_frac:
            if '/' in value:
                numerator = ''
                denominator = ''
                denomtime = False
                for i in value: #convert fraction to decimal so simulateHillModel can work
                    if i.isdigit() == True and denomtime == False:
                        numerator = numerator + i
                    elif i == '/':
                        denomtime = True
                    elif i.isdigit() == True and denomtime == True:
                        denominator = denominator + i
                decimal = float(numerator)/float(denominator)
                samples.append(str(decimal))
            else:
                samples.append(value)
        return parameternames,samples

    def _makeHillStrs(self,U,L,T,n,J):
        # Private constructor.
        scalar = "("+U+"-"+L+")"
        Xn = "X["+J+"]**"+n
        Tn = T+"**"+n
        denom = "("+Xn+"+"+Tn+")"
        neghill=scalar+"*"+Tn+"/"+denom+" + "+ L
        poshill=scalar+"*"+Xn+"/"+denom+" + "+ L
        return neghill,poshill

    def _makeHillEqns(self,eqnstr,parameternames,samples,n):
        # Private constructor.
        # X is not yet defined; eval a lambda function
        eqns=[]
        for k,e in enumerate(eqnstr):
            e2 = e #n and p are replaced during J loop, so 'if J in e' cond must use original copy of e, not changed e
            K=str(k)
            for j in range(len(eqnstr)):
                J=str(j)
                if J in e2.split(): #IE) when replacing 0p+10p, all '0 p' are replaced, so use list instead of str to replace
                    # if j affects k, find U,L,T in parameternames
                    U,L,T=None,None,None
                    for p,v in zip(parameternames,samples):
                        if J+','+K in p:
                            exec(p[0]+"= str(v)")
                        if filter(None,[U,L,T])==[U,L,T]: 
                            #quit when U,L,T assigned
                            break
                    # substitute the negative and positive hill strings
                    neghill,poshill=self._makeHillStrs(U,L,T,str(n),J)
                    e_list = e.split()
                    marked = 0
                    for i in range(len(e_list) - 1): #-1 b/c [i+1]. This does replacement
                        if e_list[i] == J and e_list[i+1] == 'p':
                            e_list[i] = poshill
                            del e_list[i+1]
                            marked = i
                        elif e_list[i] == J and e_list[i+1] == 'n':
                            e_list[i] = neghill
                            del e_list[i+1]
                            marked = i 
                    for i in range(len(e_list) - 1):
                        if i != marked:
                            e_list[i] = e_list[i] + ' '        
                    e = ''.join(e_list)
            # make a lambda function for each equation
            e="lambda X: -X["+K+"] + " + e
            eqns.append(eval(e))
        return eqns


