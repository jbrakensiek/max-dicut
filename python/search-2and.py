###
### Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick
###
### This code is licensed under the MIT license.
###

import importlib 

import dicut
importlib.reload(dicut) # in case dicut.py is modified

import experiments
importlib.reload(experiments)

import numpy as np

np.random.seed(seed=878)

LO = -1.0
HI = 1.0


b1 = 0.16471961
b2 = 0.17951539
b11 = -0.68748133
b12 = -0.68804908

rawconfigs = [(b1, b2, -1 + b1 + b2), (b2, b1, -1 + b1 + b2), (b2, b2, -1 + b2 + b2), (-b1, -b2, -1 + b1 + b2),\
              (-b2, -b2, -1 + b2 + b2), (b1, -b1, b11), (b2, -b1, b12), (0,0,-1)]

w = [0.00391633, 0.16779717, 0.20391237, 0.19492992, 0.18276885, 0.21037466, 0.03630070]

biases = np.array([-b2, -b1, 0, b1, b2])
funs = []

def remap(pw, lo, hi):
    return (hi - lo) * pw + lo

b3 = 0.3
b0 = 0.1
b4 = 0.45
b5 = 0.7

control_points = [LO, -b5, -b4, -b3, -b2, -b1, -b0, 0, b0, b1, b2, b3, b4, b5, HI]

NUM_ITER = 50 # more iterations correlates with a better distribution, but monotonicity is not guaranteed
WORKING_GRID_SIZE = 101 # number of b1/b2 evaluation points during each iteration
WORKING_GRID_DEPTH = 21 # number of b12 evaluation points during each iterations
EVAL_GRID_SIZE = 201 # same as WORKING_GRID_SIZE, but in the final evaluation
EVAL_GRID_DETPH = 51 # same as WORKING_GRID_DEPTH, but in the final evaluation
WORKING_VERBOSE = 1
EVAL_VERBOSE = 2

OUTPUT_FILE = "2and.txt"


for i in range(1):
    print("ITER    : ", i)
    print("BAISES  : ", biases)
    print("CONFIGS : ", len(rawconfigs))
    
    # convert rawconfigs so that the indices are discretized
    configs = [(int(np.where(biases == x)[0]), int(np.where(biases == y)[0]), xy)\
               for (x,y,xy) in rawconfigs]

    def value(params, weights):
        f = model_engine.eval_model(params, biases)
        rf = model_engine.eval_rev_model(params, biases)
        vf = inst.test_pairs_weighted_sparse(f, weights)
        vrf = inst.test_pairs_weighted_sparse(rf, weights)
        return (vf + vrf) / 2.0
        
    model_engine = experiments.odd_piecewise_linear_model(control_points, lo = -2, hi = 2)
    

    inst = experiments.run_experiment(biases, configs, model_engine,\
                           bootstrap_funs = funs, ntrials=100, verbose=0,\
                                     fast_stop_thresh=0.87463)
    
    
    good, cur, xb, yb, xy = experiments.fine_check(inst, model_engine,\
                                grid_size=WORKING_GRID_SIZE, grid_depth=WORKING_GRID_DEPTH,\
                                verbose=WORKING_VERBOSE, simple=False,\
                                discr=1000000, lo = LO, hi = HI,
                                merge_eps=0, manual = True)
    
    print("PROGRESS: ", cur, xb, yb, xy)
    
    newbiases = np.array([xb,-xb])
    mask = np.logical_not(np.isin(newbiases, biases))
    biases = np.append(biases, newbiases[mask])

    newbiases = np.array([yb,-yb])
    mask = np.logical_not(np.isin(newbiases, biases))
    biases = np.append(biases, newbiases[mask])
    biases.sort()
    
    
    rawconfigs.append((xb,yb,xy))
    rawconfigs.append((-yb,-xb,xy))
    
    l = len(funs)
    for i in good.keys():
        if i >= l:
            funs.append(inst.F_aux[i])
            
experiments.fine_check(inst, model_engine, grid_size=201, grid_depth=51,\
                       verbose=2, simple=False,\
                       discr=1000000, lo = -0.8, hi = 0.8,
                       merge_eps=0, manual = True, export = OUTPUT_FILE, meta=control_points)
