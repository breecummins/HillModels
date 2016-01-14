import hillmodel as hm

def do3DExample():
    # specify the full path to your file, unless it is in the same directory
    networkfile = '/Users/bcummins/GIT/DSGRN/networks/3D_Example.txt'  

    # for didactic purposes, I am writing a sample file; generally you will be given the sample file
    with open('samples.txt','w') as f:
        f.write('L(X,Y) 1.2\nU(X,Y) 3.4\nT(X,Y) 1.5\nL(Y,Z) 1.0\nU(Y,Z) 2.5\nT(Y,Z) 2.4\nL(Z,X) 0.6\nU(Z,X) 2.1\nT(Z,X) 2.0')
    samplefile='samples.txt'

    # choose a Hill exponent
    hillexp = 10

    # make an instance of class hillmodel
    Example3D = hm.hillmodel(networkfile,samplefile,hillexp)

    # choose initial conditions and time period
    y0 = [1.0, 2.0, 1.5] # there are three variables in this network
    t0 = 0
    t1 = 10 # start at 0 and run for 10 time units
    dt = 0.01 # give me a new data point every hundredth of a time unit

    # run the Hill model
    times, timeseries = Example3D.simulateHillModel(y0,t0,t1,dt) 

    # choose plotting options
    plotoptions={'linewidth':2}
    legendoptions={'fontsize':24,'loc':'upper left', 'bbox_to_anchor':(1, 1)} 
    # note: if you plt.savefig instead of plt.show, these options ensure that the 
    # legend won't be cut off, even if the legend labels are long
    figuresize = (15,10)

    # plot the results 
    Example3D.plotResults(times,timeseries,plotoptions,legendoptions,figuresize)

if __name__ == '__main__':
    do3DExample()