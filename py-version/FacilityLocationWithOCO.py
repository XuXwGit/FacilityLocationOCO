import numpy as np
import gurobipy as grb
from gurobipy import GRB
from gurobipy import *
import matplotlib.pyplot as plt


I = 10
J = 100
K = 10
P = 1
S = 100

def generate_data(I=10, J=100, K=10, P=1, S=100):
    np.random.seed(0)

    fixed_lb = 20
    fixed_ub = 40
    fk = np.random.randint(fixed_lb, fixed_ub + 1, size=(K))

    cost_lb = 1
    cost_ub = 10
    C1 = np.random.randint(cost_lb, cost_ub + 1, size=(I,J,K,P))
    C2 = np.random.randint(cost_lb, cost_ub + 1, size=(I,J,K,P))

    demand_lb = 5
    demand_ub = 50
    demand = np.random.randint(demand_lb, demand_ub + 1, size=(S, J, P))

    beta_lb = 0.8
    beta_ub = 1.0
    beta = np.random.uniform(beta_lb, beta_ub, size=(J,P))

    return fk,C1,C2,demand,beta

# LANGRANGIAN RELAXATION
def Lagrangean_relaxation():
    model = grb.Model('model')
    Y = model.addVars(K, lb = 0, vtype=GRB.CONTINUOUS, name='y')
    X = model.addVars(S, I, J, K, P, lb = 0, vtype=GRB.CONTINUOUS, name='x')
    penalty = model.addVars(J,P,vtype=GRB.CONTINUOUS, name='p')
    model.addConstrs((penalty[j,p]==(beta[j][p]*demand[:,j,p].sum()
                                        -X.sum('*','*',j,'*',p))/S
                  for j in range(J)
                  for p in range(P)), 
                  name='constr1')
    model.addConstrs((X.sum(s, '*', j, '*', p) <= demand[s][j][p]
                    for s in range(S)
                    for j in range(J)
                    for p in range(P)), 
                    name='constr2')
    model.addConstrs((X.sum(s, '*', '*', k, '*') <= Y[k]
                    for s in range(S)                                                                                                                                 
                    for k in range(K)), 
                    name='constr3')
    # Initialize Lagrangian Multipliers for the assignment contraints
    Lambda = np.ones(shape=(J,P))
    
    # Placeholder for penalties
    penalties = np.ones(shape=(J,P))

    tolerance = 1e-3
    iter_count = 0
    step = 1

    obj_value = []
    gradient_norm = []

    while not (Lambda.any() <= tolerance) | (penalties.any() <= tolerance):
        # Initialize objective function
        obj = LinExpr()
        for k in range(K):
            obj += fk[k]*Y[k]
        for s in range(S):  
            for i in range(I):
                for j in range(J):
                    for p in range(P):
                        for k in range(K):
                            obj += (C1[i,j,k,p] + C2[i,j,k,p])/S * X[s,i,j,k,p]
        for j in range(J):
            for p in range(P):
                obj += Lambda[j][p] * penalty[j,p]
        # Set objective
        model.setObjective(obj, GRB.MINIMIZE)
            
        # Solve model
        model.setParam('OutputFlag', 0)
        model.optimize()
        
        # Langrangian Multiplier penalties
        penalties = np.array([penalty[j,p].x 
                            for j in range(J) 
                            for p in range(P)]).reshape((J,P))
        
        # Increase counter
        iter_count += 1

        step = 1 / np.sqrt(iter_count)
        # Update Langrangian Multiplier
        for j in range(J):
            for p in range(P):
                Lambda[j][p] = max(Lambda[j][p] + step*(penalties[j][p]), 0)
        
        # Print status
        # Print status
        print('\tIteration: %d' % iter_count, 
              '\tObjective value = %.2f'% model.ObjVal,
              '\tGradient = %.2f'% (penalties.max()),
              '\tLambda = %.2f'% Lambda.max())
        # print(Lambda)

        obj_value.append(model.ObjVal)
        gradient_norm.append(penalties.max())

    for k in range(K):
        print(Y[k].x)

    return obj_value,gradient_norm

def draw_LR(obj_value, gradient_norm):
   fig = plt.figure(figsize=(12,3), dpi = 1080)

   (ax1,ax2) = fig.subplots(1,2)

   fig.subplots_adjust(hspace = 0.40)
   plt.title('I=10,J=100,K=10,S=200')
   ax1.plot(obj_value,color = 'blue', linewidth = 0.5)
   ax2.plot(gradient_norm,color = 'blue', linewidth = 0.5)

   ax1.set_title("(a) Optimal Objective")
   ax2.set_title('(b) Gradient Norm')

   ax1.set_xlabel('Iteration')
   ax2.set_xlabel('Iteration')

   fig.savefig("intance100.eps", bbox_inches = 'tight')
   fig.savefig("instance100-eps-converted-to.pdf", bbox_inches = 'tight')

def direct_solve():
    model = grb.Model('model')
    Y = model.addVars(K, lb = 0, vtype=GRB.CONTINUOUS, name='y')
    X = model.addVars(S, I, J, K, P, lb = 0, vtype=GRB.CONTINUOUS, name='x')
    model.addConstrs((beta[j][p]*demand[:,j,p].sum()-X.sum('*','*',j,'*',p)<=0
                    for j in range(J)
                    for p in range(P)), 
                    name='constr1')
    model.addConstrs((X.sum(s, '*', j, '*', p) <= demand[s][j][p]
                    for s in range(S)
                    for j in range(J)
                    for p in range(P)), 
                    name='constr2')
    model.addConstrs((X.sum(s, '*', '*', k, '*') <= Y[k]
                    for s in range(S)
                    for k in range(K)), 
                    name='constr3')
    obj = LinExpr()
    for k in range(K):
            obj += fk[k]*Y[k]
    for s in range(S):  
        for i in range(I):
            for j in range(J):
                for p in range(P):
                    for k in range(K):
                        obj += (C1[i,j,k,p] + C2[i,j,k,p])/S * X[s,i,j,k,p]
    # Set objective
    model.setObjective(obj, GRB.MINIMIZE)
    # Solve model
    model.setParam('OutputFlag', 0)
    model.optimize()
    print('Objective value =', model.objVal)

    for k in range(K):
         print(Y[k].x)

def solve_determine_model():
    model = grb.Model('determine.lp')
    Y = model.addVars(K, lb = 0, vtype=GRB.CONTINUOUS, name='y')
    X = model.addVars(I, J, K, P, lb = 0, vtype=GRB.CONTINUOUS, name='x')
    constrs = model.addConstrs((X.sum('*',j,'*',p) >= beta[j][p]*demand[:,j,p].sum()/S
                    for j in range(J)
                    for p in range(P)), 
                    name='constr1')
    model.addConstrs((X.sum('*', j, '*', p) <= demand[:,j,p].sum()/S
                    for j in range(J)
                    for p in range(P)), 
                    name='constr2')
    model.addConstrs((X.sum('*', '*', k, '*') <= Y[k]
                    for k in range(K)), 
                    name='constr3')
    obj = LinExpr()
    for k in range(K):
            obj += fk[k]*Y[k]
    for i in range(I):
        for j in range(J):
            for p in range(P):
                for k in range(K):
                    obj += (C1[i,j,k,p] + C2[i,j,k,p]) * X[i,j,k,p]
    # Set objective
    model.setObjective(obj, GRB.MINIMIZE)
    # Solve model
    model.setParam('OutputFlag', 0)
    model.optimize()
    print('Determine Objective value =', model.objVal)

    Yvalue = np.zeros(shape=(K))
    Xvalue = np.zeros(shape=(I,J,K,P))
    lambdaValue = np.zeros(shape=(J,P))

    for k in range(K):
         Yvalue[k] = Y[k].x

    for i in range(I):
         for j in range(J):
              for k in range(K):
                   for p in range(P):
                        Xvalue[i,j,k,p] = X[i,j,k,p].x

    for j in range(J):
         for p in range(P):
              lambdaValue[j,p] = constrs[j,p].pi
          #     lambdaValue[j,p] = demand[:,j,p].sum()/S -  Xvalue[:,j,:,p].sum()

    return Yvalue, Xvalue, lambdaValue, model.objVal

def OCO_based_Lagrangean_relaxation(N=50):
    model = grb.Model('sp_model')
    alpha = model.addVars(J,P,lb = -GRB.INFINITY,ub = 0, vtype = GRB.CONTINUOUS, name = 'alpha')
    gamma = model.addVars(K,lb = -GRB.INFINITY,ub = 0, vtype = GRB.CONTINUOUS, name = 'gamma')
    Yvalue = np.zeros(shape=(K))
    # Initialize Lagrangian Multipliers for the assignment contraints
    Y0,X,Lambda,DetermineObj = solve_determine_model()

    objVal_set = []
    objVal_set.append(DetermineObj)

    # penalties
    penalties = np.ones(shape=(J,P))

    constrs = model.addConstrs((alpha[j,p] + gamma[k] <= 
                                C1[i,j,k,p] + C2[i,j,k,p] - Lambda[j,p] 
                     for i in range(I)
                     for j in range(J)
                     for k in range(K)
                     for p in range(P)), name = 'constr')
    
    tolerance = 1e-3
    iter_count = 0
    step = 1

    while not (Lambda.any() <= tolerance) | (penalties.any() <= tolerance) | (N < iter_count):
        Yvalue = Y0.copy()
    
        # reset right bound
        for i in range(I):
            for j in range(J):
                for k in range(K):
                    for p in range(P):
                        constrs[i,j,k,p].RHS = C1[i,j,k,p] + C2[i,j,k,p]-Lambda[j,p] 

        Y0 = Yvalue.copy()/S
        X = np.zeros(shape=(I,J,K,P))

        for s in range(S):
            # Set objective
            obj = LinExpr()
            for j in range(J):
                for p in range(P):
                    obj += demand[s][j][p] *alpha[j,p]
            for k in range(K):
                obj += Yvalue[k] * gamma[k]
            model.setObjective(obj, GRB.MAXIMIZE)

            # Solve model
            model.setParam('OutputFlag', 0)
            model.optimize()
            
            step = 1 / np.sqrt(s+1)
            for k in range(K):
                # print('%.2f'%Yvalue[k],'\t',fk[k]+gamma[k].x,'\t','%.2f'%(Yvalue[k]-step*(fk[k]+gamma[k].x)))
                Yvalue[k] = max(Yvalue[k]-step*(fk[k]+gamma[k].x), 0)
                if k != K-1:
                    Y0[k] += Yvalue[k]/S

            for i in range(I):
                for j in range(J):
                    for k in range(K):
                        for p in range(P):
                            X[i,j,k,p] += constrs[i,j,k,p].pi

        X = X / S
        
        # Increase counter
        iter_count += 1
        
        # Langrangian Multiplier penalties
        step = 1 / np.sqrt(iter_count)

        total_cost = 0

        for j in range(J):
            for p in range(P):
                penalties[j][p] = beta[j][p] * demand[:,j,p].sum()/S - X[:,j,:,p].sum()
                tempLam = Lambda[j][p]
                Lambda[j][p] = max(Lambda[j][p] + step*(penalties[j][p]), 0)
                total_cost += Lambda[j,p] * penalties[j,p]
                # print(tempLam, '\t', step * penalties[j][p],'\t', Lambda[j][p])

        for k in range(K):
            total_cost += fk[k] * Y0[k]
        
        for i in range(I):
            for j in range(J):
                for k in range(K):
                    for p in range(P):
                        total_cost += (C1[i,j,k,p] + C2[i,j,k,p])*X[i,j,k,p]

        objVal_set.append(total_cost)

        # Print status
        print('\tIteration: %d' % iter_count, 
              '\tObjective value = %.2f'% total_cost,
              '\tGradient = %.2f'% (penalties.max()),
              '\tLambda = %.2f'% Lambda.max())
    return objVal_set

def draw_OCO_LR(obj_value):
   fig = plt.figure(figsize=(12,3), dpi = 1080)
   plt.plot(obj_value,color = 'blue', linewidth = 0.5)
   plt.title('Convergence Process(I=10,J=100,K=10,S=100)')
   plt.xlabel('Iteration')
   plt.ylabel('Objective')
   fig.savefig("intance-II.eps", bbox_inches = 'tight')
   fig.savefig("instance-II-eps-converted-to.pdf", bbox_inches = 'tight')

def test():
    # generate data
    fk,C1,C2,demand,beta = generate_data(I=10,J=100,K=10,P=1,S=100)
    # use OCO-LR method to solve
    objVal_set2 = OCO_based_Lagrangean_relaxation(100)
    draw_OCO_LR(objVal_set2)
    # use traditional LR method to solve
    obj_value,gradient_norm = Lagrangean_relaxation()
    draw_LR(obj_value,gradient_norm)
    # directly solve
    direct_solve()

if __name__ == "__main__":
    test()