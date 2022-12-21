###
### Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick
###
### This code is licensed under the MIT license.
###

import dicut
from datetime import datetime
import numpy as np
import math
import time
import matplotlib.pyplot as plt
import random
from scipy.optimize import minimize, dual_annealing, shgo
from scipy.optimize import LinearConstraint
from pympler import muppy, summary

class ModelEngine:
    """
    
    User provides these:
    
    nparam -- number of parameters
    
    param_to_fun(params) -- turns a list of length nparam to a function on biases
    
    model_bounds -- bounds on each parameter in the model (note that bounds are independent)
        -- given in the form [(lo0,hi0),(lo1,hi1),...]
    
    This will also be available:
    
    param_to_rev_fun -- returns (x -> -f(-x)), where f is the output of param_to_fun
    
    get_model_bounds -- returns model bounds
    
    get_random_params -- return random model within the model_bounds
    
    eval_model(params, biases)
    
    eval_rev_model(params, biases)
    
    """
    
    def __init__(self, nparam, param_to_fun, model_bounds, comment="no comment"):
        self.comment = comment
        self.nparam = nparam
        self.param_to_fun = param_to_fun
        self.model_bounds = np.array(model_bounds)
        assert(self.nparam == len(self.model_bounds))
    
    def param_to_rev_fun(self, params):
        f = self.param_to_fun(params)
        return lambda x : -f(-x)
    
    def get_model_bounds(self):
        return self.model_bounds
    
    def get_random_params(self):
        return np.random.uniform(self.model_bounds[:,0],\
                                 self.model_bounds[:,1])
    
    def eval_model(self, params, biases):
        return np.vectorize(self.param_to_fun(params))(biases)
    
    def eval_rev_model(self, params, biases):
        return np.vectorize(self.param_to_rev_fun(params))(biases)
        
def piecewise_linear_model(biases, lo=-2.0, hi=2.0):
    param_to_fun = lambda params : (lambda x : np.interp(x,biases,params))
    
    model_bounds = [(lo,hi) for i in biases]
    
    return ModelEngine(len(biases), param_to_fun,model_bounds,\
                       comment="piecewise linear on %s" % str(biases))

def odd_piecewise_linear_model(biases, lo=-2.0, hi=2.0):
    assert(np.linalg.norm(biases + np.flip(biases)) < 1e-6)
    
    expand_params = lambda params : np.concatenate((np.flip(-params), [0], params))
    
    param_to_fun = lambda params : (lambda x : np.interp(x,biases,expand_params(params)))
    
    nparams = (len(biases)-1)//2
    
    model_bounds = [(lo,hi) for i in range(nparams)]
    
    return ModelEngine(nparams, param_to_fun,model_bounds,\
                       comment="odd piecewise linear on %s" % str(biases))
    
def run_experiment(biases, configs, model_engine, bootstrap_funs = [], ntrials = 100, stopeps=1e-5, verbose=2, fast_stop_thresh = 1.0, flip=True):
    """
    runs experiment on biases and configs where it trains models in accordance
    with model_engine
    
    bootstrap_funs is additional functions to prime the model with (avoids early training)
    
    ntrials is number of steps the experiment will run.
    
    verbose is level of output
    
    verbose = 0 is no output
    verbose = 1 is text output only (no charts)
    verbose = 2 is text and charts
    """
    inst = dicut.Config(biases, configs, check_tri = False)
    for f in bootstrap_funs:
        inst.add_function(f, aux=f)
    
    if verbose >= 1:
        print("Initial quality:", inst.current_approx_ratio())
    
    def value(params, weights):
        f = model_engine.eval_model(params, biases)
        rf = model_engine.eval_rev_model(params, biases)
        vf = inst.test_pairs_weighted_sparse(f, weights)
        vrf = inst.test_pairs_weighted_sparse(rf, weights)
        return (vf + vrf) / 2.0
    
    for i in range(ntrials):
        if verbose >= 1:
            print("STEP:", i)
        w, obj = inst.next_config_dist(verbose=verbose)
        if obj >= fast_stop_thresh:
            print("FAST STOP")
            return inst
            
        res = minimize(lambda x : -value(x, w),\
                       model_engine.get_random_params(),\
                       bounds=model_engine.get_model_bounds())
        if verbose >= 1:
            print(res)
        if -res.fun > obj:
            model = model_engine.eval_model(res.x, biases)
            rev_model = model_engine.eval_rev_model(res.x, biases)
            inst.add_tlist(    model, aux=model_engine.param_to_fun(res.x))
            if flip:
                inst.add_tlist(rev_model, aux=model_engine.param_to_rev_fun(res.x))
            if verbose >= 2:
                plt.plot(inst.bias, model)
                plt.plot(inst.bias, rev_model)
                plt.title("Value: %f" % (-res.fun))
                plt.show()
        if abs(-res.fun - obj) < stopeps:
            print("CONVERGENCE STOP")
            return inst
    
    return inst
                
    
    
def fine_check(inst, model_engine, verbose=2, lo=-1.0, hi=1.0,grid_size=101, grid_depth=21, merge_eps=1e-4, drop_eps=1e-6, simple=False, discr=100, manual = False, export = None, meta=None):
    apx = inst.current_approx_ratio()
    if verbose >= 1:
        print("Trained quality:", apx)
    N = len(inst.F_eval)
    apx2 = inst.hist_F.latest_approx()

    assert(apx == apx2)
    p = inst.hist_F.latest_dist()
    assert(abs(sum(p) - 1) < 1e-5)

    
    l = [(p[i],i) for i in range(N)]
    l.sort()
    l = [x for x in reversed(l)]
    good = {}
    for p,i in l:
        if p < drop_eps:
            break
        #print(good.keys())
        found = False
        for ii in good.keys():
            if np.linalg.norm(inst.F_param[i] - inst.F_param[ii]) < merge_eps:
                found = True
                good[ii] += p
                break
        if not found:
            good[i] = p
    
    xax = np.linspace(lo, hi, grid_size)
    if verbose >= 2:
        for i in good.keys():
            p = good[i]
            plt.plot(xax, np.vectorize(inst.F_aux[i])(xax), label="%.3f"%p)
        plt.legend(bbox_to_anchor=(1,1), loc="upper left")
        plt.title("Ratio: %.5f"%apx)
        plt.xlabel("bias")
        plt.ylabel("threshold")
        plt.show()

    if verbose >= 1:
        print("Claimed Ratio:", inst.hist_hard.latest_approx())
    if verbose >= 1:
        print("Weight, (  b1  ,   b2  ,  b12  , value)")
        w = inst.hist_hard.latest_dist()
        for i in range(len(inst.configs)):
            a,b,c,d = inst.configs[i]
            if (w[i] > 1e-9):
                print("%0.06f, (%+0.03f, %+0.03f, %+0.03f, %0.03f)"%\
                      (w[i],inst.bias[a],inst.bias[b],c,d))

    dist = []
    ss =  sum([good[i] for i in good.keys()])
    
    for i in good.keys():
        p = good[i] / ss
        dist.append((p, inst.F_aux[i]))

    if export != None:
        file = open(export, "w")
        bb = len(meta)#model_engine.get_model_bounds())
        gg = len(dist)
        file.write("%d %d\n"%(bb, gg))
        for i in range(bb):
            file.write("%.15lf "%meta[i])
        file.write("\n")
        for i in range(gg):
            file.write("%.15lf"%dist[i][0])
            for j in range(bb):
                file.write(" %.15lf"%dist[i][1](meta[j]))
            file.write("\n")
        
        file.close()
    
    def my_eval(l):
        x = l[0]
        y = l[1]
        relxy = l[2]
        xy = relxy * (2 - abs(x+y) - abs(x-y)) + abs(x+y)-1
        c = dicut.value(x,y,xy)
        s = 0
        for p,tl in dist:
            s += p * dicut.thresh_eval(x,y,xy,tl(x),tl(y))
        #if c < 1e-9:
        #    print("WHAT: ", s, c)
        return s / (c + 1e-12)
    
    res = 0
    x = 0
    y = 0
    xy = 0
    cur = 1.0
    hii = 0
    hj = 0
    hk = 0
    
    if not manual:
        res = minimize(my_eval, np.random.uniform([lo,lo,0],[hi,hi,1]), bounds=[(lo,hi),(lo,hi),(0,1)])
        x = res.x[0]
        y = res.x[1]
        relxy = res.x[2]
        xy = relxy * (2 - abs(x+y) - abs(x-y)) + abs(x+y)-1
        xy = min(max(xy, abs(x+y)-1+1e-9), 1-abs(x-y)-1e-9)
        res = res.fun
    else:
        for i in range(grid_size):
            for j in range(grid_size):
                a = abs(xax[i]+xax[j]) - 1
                b = -abs(xax[i]-xax[j]) + 1
                aa = min(a,b)
                bb = max(a,b)
                for k in np.linspace(0,1,grid_depth):
                    at = my_eval([xax[i],xax[j],k])
                    if at < cur:
                        cur = at
                        hii = i
                        hj = j
                        hk = k
            if verbose >= 2:
                print(i, cur, xax[hii], xax[hj], hk)
                
        res = cur
        
        x = xax[hii]
        y = xax[hj]
        relxy = hk
        
        print(x, y, relxy, flush=True)
        
        res2 = minimize(my_eval, [x, y, relxy], bounds=[(lo*1.0,hi*1.0),(lo*1.0,hi*1.0),(0.0,1.0)])
        x2 = res2.x[0]
        y2 = res2.x[1]
        relxy2 = res2.x[2]
        res2 = res2.fun
        
        if res2 < res and x2 >= lo and x2 <= hi and y2 >= lo and y2 <= hi and relxy2 >= 0 and relxy2 <= 1:
            res = res2
            x = x2
            y = y2
            relxy = relxy2

        xy = relxy * (2 - abs(x+y) - abs(x-y)) + abs(x+y)-1
        xy = min(max(xy, abs(x+y)-1), 1-abs(x-y))
            
        assert (x >= -1 and x <= 1)
    
    return good, res, x, y, xy
