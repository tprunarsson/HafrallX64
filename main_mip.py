import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time
import argparse
from itertools import combinations
from utils import DistanceLink, readDat, writeDat, drawTour 
(tSHIP,tSTAT,tWAYP,tENDP,tPORT) = (1,2,3,4,5)

def subtour_(vals):
    global n_global
    # make a list of edges selected in the solution
    edges = gp.tuplelist((i, j) for i, j in vals.keys()
                         if vals[i, j] > 0.5)
    n = n_global
    unvisited = list(range(n))
    cycle = list(range(n+1))  # initial length has 1 more city
    listlen = []
    while unvisited:  # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i, j in edges.select(current, '*')
                         if j in unvisited]
        if (len(cycle) >= len(thiscycle)):
            cycle = thiscycle
        listlen.append(len(thiscycle))
    print(n, listlen)
    return cycle

# Callback - use lazy constraints to eliminate sub-tours
def subtourelim_(model, where):
    global n_global
    n = n_global
    if where == GRB.Callback.MIPSOL:
        vals = model.cbGetSolution(model._vars)
        # find the shortest cycle in the selected edge list
        tour = subtour_(vals)
        if len(tour) < n:
            # add subtour elimination constr. for every pair of cities in tour
            print("adding tour elimination ___ :", n, len(tour),tour)
            model.cbLazy(gp.quicksum(model._vars[i, j] + model._vars[j, i]
                                     for i, j in combinations(tour, 2))
                         <= len(tour)-1)

# Given a tuplelist of edges, find the shortest subtour
n_global = 0

def solveit_(dist, n, Amount, shipCap, Type, timelimit = 3600, useAmount = True, useSubtour = True):

    global n_global
    n_global = 2*n
    m = gp.Model()
    n = 2*n

    #m.setParam('OutputFlag', 0)
    m.setParam("TimeLimit", timelimit)
    m.setParam('Threads', 4) 

    n2 = int(n/2)
    
    distx = dist[:n,:n]
    distx[0,1] = 0
    distx[1,0] = 0 # close the loop for start and end, this distance should be zero of not counted!

    # Create variables
    N = [(i,j) for i in range(n) for j in range(n)]
    vars = m.addVars(N, obj=distx, vtype=GRB.BINARY, name='e')
    
    if useAmount:
        w = m.addVars(N, ub = shipCap, vtype=GRB.CONTINUOUS, name='w')
        N1 = [(2*i,2*i+1) for i in range(n2) if Type[i] != tPORT] + [(2*i+1,2*i) for i in range(n2) if Type[i] != tPORT]
        v = m.addVars(N1, ub = shipCap, vtype=GRB.CONTINUOUS, name='v')

    # Add degree-2 constraint
    m.addConstrs(gp.quicksum(vars[i,j] for j in range(n) if (i,j) in N) == 1 for i in range(n))
    m.addConstrs(gp.quicksum(vars[j,i] for j in range(n) if (j,i) in N) == 1 for i in range(n))
    #m.addConstrs(vars[i,j] + vars[j,i] <= 1 for (i,j) in N) # this is probably not necessary!

    m.addConstrs(vars[2*i,2*i+1]+vars[2*i+1,2*i] == 1 for i in range(int(n/2)))
    m.addConstrs(vars[i,i] == 0 for i in range(n))
    
    # the flow into
    if useAmount:
        m.addConstrs(gp.quicksum(w[k,2*j] for k in range(n) if k != (2*j+1)) + v[2*j,2*j+1] == 
                 w[2*j,2*j+1] for j in range(n2) if Type[j] != tPORT)
    
        m.addConstrs(gp.quicksum(w[k,2*j+1] for k in range(n) if k != (2*j)) + v[2*j+1,2*j] == 
                 w[2*j+1,2*j] for j in range(n2) if Type[j] != tPORT )
        
        m.addConstrs(v[2*j,2*j+1]+v[2*j+1,2*j] == Amount[j] for j in range(n2) if Type[j] != tPORT)
    
        m.addConstrs(w[2*j,2*j+1] == gp.quicksum(w[2*j+1,k] for k in range(n) if k != 2*j and k != (2*j+1)) for j in range(n2) if Type[j] != tPORT)
        m.addConstrs(w[2*j+1,2*j] == gp.quicksum(w[2*j,k] for k in range(n) if k != (2*j+1) and k != 2*j) for j in range(n2) if Type[j] != tPORT)
 
    # I have carefully forced zeros to make the equations above with fewer conditions!
        m.addConstrs(w[i,j] <= shipCap*vars[i,j] for (i,j) in N) # only fish where there is a connection
        # not sure what this one does, probablu not needed:
        m.addConstrs(v[i,j] <= shipCap*vars[i,j] for (i,j) in N1) # only fish where there is a connection
        m.addConstrs(w[i,j] == 0 for (i,j) in N if Type[int(i/2)] == tPORT) # all flow from port is zero
    
    # Optimize model
    if useSubtour:
        m._vars = vars
        m.Params.LazyConstraints = 1
        m.optimize(subtourelim_)
    else:
        m.optimize()

    vals = m.getAttr('X', vars)
    
    tour = subtour_(vals)
    #assert len(tour) == n
    
    W = np.zeros((n,n))
    V = np.zeros((n,n))
    X = np.zeros((n,n))
    # extract the fish flow
    for i in range(n):
        for j in range(n):
            if (i != j):
                X[i,j] = vars[i,j].X
                if useAmount:
                    W[i,j] = w[i,j].X
                    if (i,j) in N1:
                        V[i,j] = v[i,j].X
    
    tour = np.array(tour)
    
    letour = []
    for i in range(n):
        if int(tour[i]/2) not in letour and -int(tour[i]/2) not in letour:
            if tour[i] % 2 == 1:
                letour.append(-int(tour[i]/2))
            else:
                letour.append(int(tour[i]/2))
    letour = np.array(letour)
    letour = np.concatenate((letour[np.where(letour==0)[0][0]:],letour[:np.where(letour==0)[0][0]]))
    if Type[abs(letour[1])] == tPORT:
        letour[1:] = -letour[:0:-1]
    letour[Type[np.abs(letour)]==tPORT] = np.abs(letour[Type[np.abs(letour)]==tPORT])
    print('Optimal tour: %s' % str(letour))
    print('Optimal cost: %g' % m.ObjVal)
    return m.ObjVal, letour, W, V, X

def runit(file, ship_id, name, run, mode, ShipCap, TimeWindow = None):
    (tSHIP,tSTAT,tWAYP,tENDP,tPORT) = (1,2,3,4,5)  
    ship_names = ["Árni Friðriksson","Bjarni Sæmundsson","Gullver","Breki"]
    ship = ship_names[ship_id-1]
    print(ship)
    (Tour, Type, Amount, Fixed, LatLonRad, Name, Rotated, ShipCap_, ExtraTime) = readDat(file, ship)
    if ShipCap == None:
        ShipCap = ShipCap_
# Extract the available ports
    TypePort = list(set([Name[i] for i in range(len(Type)) if Type[i] == tPORT]))
    TypePort = {type: len([Name[i] for i in range(len(Type)) if Name[i] == type]) for type in TypePort}
# Extract the vailable waypoints
    TypeWayp = list(set([Name[i] for i in range(len(Type)) if Type[i] == tWAYP]))
    TypeWayp = {type: len([Name[i] for i in range(len(Type)) if Name[i] == type]) for type in TypeWayp}
# Stack on the bottom the waypoints and top the ship, should only be one ship
    ExLatLonRad = np.vstack((LatLonRad[Type == tSHIP,:], LatLonRad[Type == tSTAT,:], LatLonRad[Type == tPORT,:] ,LatLonRad[Type == tWAYP,:]))
    ExType = np.concatenate((Type[Type == tSHIP],Type[Type == tSTAT],Type[Type == tPORT],Type[Type == tWAYP]))
    ExAmount = np.concatenate((Amount[Type == tSHIP],Amount[Type == tSTAT],Amount[Type == tPORT],Amount[Type == tWAYP]))
    print(ExtraTime)
    ExExtraTime = np.concatenate((ExtraTime[Type == tSHIP],ExtraTime[Type == tSTAT],ExtraTime[Type == tPORT],ExtraTime[Type == tWAYP]))
    ExFixed = np.concatenate((Fixed[Type == tSHIP],Fixed[Type == tSTAT],Fixed[Type == tPORT],Fixed[Type == tWAYP]))
    Name = np.array(Name)
    ExName = np.concatenate((Name[Type == tSHIP],Name[Type == tSTAT],Name[Type == tPORT],Name[Type == tWAYP]))
    ExRotated = np.concatenate((Rotated[Type == tSHIP],Rotated[Type == tSTAT],Rotated[Type == tPORT],Rotated[Type == tWAYP]))
    SelectedSize = ExLatLonRad.shape[0]
    Size = np.sum(ExType == tSTAT) + np.sum(ExType == tPORT) + np.sum(ExType == tSHIP)
    print("SelectedSize", SelectedSize, "Size", Size, "ExType=", ExType, "tSHIP=", tSHIP)
    StartEnd = ExLatLonRad[ExType == tSHIP,:][0]
    print("StartEnd=", StartEnd)
    print("ShipCap=", ShipCap)
    Elite = np.zeros(Size)
    print(ExName)
    # Place the code you want to time here
    (DistMtrx, FsbleMtrx) = DistanceLink(ExType, ExLatLonRad, StartEnd, Size, SelectedSize)
    start_time = time.time()
    if mode == "gurobi":
        (obj, e, W, V, X) = solveit_(DistMtrx, Size, ExAmount, ShipCap, ExType, timelimit = 3600)
    else:
        print("mode is unknown!")
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")
    drawTour(e, ExLatLonRad, ExType, ExAmount, DistMtrx, FsbleMtrx)
    print(obj,e)
    writeDat('output_'+mode+'_'+ship[:4]+'_'+str(int(ShipCap))+'_'+name+'_'+str(run)+'.dat', e, mode, ExType, ExAmount, ExFixed, ExLatLonRad, ExName, ExtraTime, ship)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Groundfish Survey Planner")
    parser.add_argument(
        "--file", type=str, default="data/data2023spring.dat",
        help="Path to the data file (default: %(default)s)"
    )
    parser.add_argument(
        "--ship", type=int, default=1,
        help="The ship to optimize(default: %(default)s)"
    )
    parser.add_argument(
        "--name", type=str, default="tmp",
        help="name of experiment (default: %(default)s)"
    )
    parser.add_argument(
        "--run", type=int, default=0,
        help="name number of experiment (default: %(default)s)"
    )
    parser.add_argument(
        "--mode", type=str, default="gurobi", choices=["gurobi"],
        help="Mode of operation (default: %(default)s)"
    )
    parser.add_argument(
        "--capacity", type=int, default=None,
        help="ship capacity (default: %(default)s)"
    )

    args = parser.parse_args()
    runit(args.file, args.ship, args.name, args.run, args.mode, args.capacity)
