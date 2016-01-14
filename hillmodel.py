import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'normal',
        'size'   : 22}
matplotlib.rc('font', **font)

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
        parameternames,samples = self._parseSamples(self.varnames,samplefile)
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
        eqns=[]
        for l in f:
            L=l.split(' : ')
            varnames.append(L[0])
            eqns.append(L[1])
        f.close()
        eqnstr=[]
        for e in eqns:
            for k,v in enumerate(varnames):
                e=e.replace('~'+v,str(k)+' n').replace(v,str(k)+' p').replace(')(',')*(')
            eqnstr.append(e)
        return eqnstr,varnames

    def _parseSamples(self,varnames,fname='samples.txt'):
        # Private parser.
        f=open(fname,'r')
        pnames=[]
        samples=[]
        for l in f:
            L=l.split()
            pnames.append(L[0])
            samples.append(L[1])
        f.close()
        parameternames=[]
        for p in pnames:
            for k,v in enumerate(varnames):
                p=p.replace(v,str(k))
            parameternames.append(p)
        return parameternames,samples

    def _makeHillStrs(self,U,L,T,n,j):
        # Private constructor.
        scalar = "("+U+"-"+L+")"
        Xn = "X["+j+"]**"+n
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
            K=str(k)
            for j in range(len(eqnstr)):
                J=str(j)
                if J+' ' in e: 
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
                    e=e.replace(J+' n',neghill).replace(J+' p',poshill)
            # make a lambda function for each equation
            e="lambda X: -X["+K+"] + " + e
            eqns.append(eval(e))
        return eqns


