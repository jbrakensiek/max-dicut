###
### Copyright (c) 2022 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick
###
### This code is licensed under the MIT license.
###

import numpy as np
import math
import cvxpy as cp
import dill
from scipy.stats.mvn import mvnun
import matplotlib.pyplot as plt

# convention: clause structure is "(NOT x1) AND x2"
# bias is dot product with FALSE

def thresh_eval(a1, a2, b12, t1, t2):
    """
    a1/a2 are the individual bias
    b12 is the common bias
    t1/t2 are the THRESH- thresholds
    """
    e1 = math.sqrt(1+1e-9 - a1*a1)
    e2 = math.sqrt(1+1e-9 - a2*a2)

    rho = (b12 - a1 * a2) / (e1 * e2)
    rho = max(min (rho, 1), -1)
    A = np.array([[1, rho], [rho, 1]], dtype=object)
    return mvnun([-np.inf,t2], [t1,np.inf], [0,0], A, abseps=1e-10, releps=1e-10)[0]


## LLZ rounding function

llz_theta = [0, 1.24677, 1.30702, 1.36867, 1.39760, 1.44728, 1.50335, np.pi/2]
llz_pairs = [(1,7),(2,6),(3,5),(4,4),(4,3),(5,2),(6,1)]
llz_rot = [0.88955971, 1.42977923, 1.45654677, 1.48374672, 1.49640918, 1.51791721, 1.54210593, np.pi/2]
#tt = [(2.5/np.tan(rot[i]))*np.sin(theta[i]) for i in range(8)]

def get_LLZ(b):
    if b < 0:
        return -get_LLZ(-b)
    t = math.acos(b)
    r = np.interp(t, llz_theta, llz_rot)
    return 2.5/np.tan(r)*np.sin(t)/np.sqrt(1+1e-9-b*b) # sign?

def value(bi,bj,bij):
    return (1+bi-bj-bij)/4

# greater than to avoid degenerate instances
def check(bi, bj, bij,eps=1e-9):
    return bi + bj + bij >= -1-eps and - bi - bj + bij >= -1-eps and -bi + bj - bij >= -1-eps and bi - bj - bij > -1+eps

def check_tight(bi, bj, bij, eps = 0.01):
    return check(bi, bj, bij) and (bi + bj + bij <= -1+eps or - bi - bj + bij <= -1+eps or -bi + bj - bij <= -1+eps)

class ApproxHistory:
    def __init__(self):
        self.approx = []
        self.dist = []
    
    def insert(self,approx, dist):
        self.approx = [approx]
        self.dist = [dist]

    def latest_approx(self):
        if self.approx == []:
            return None
        return self.approx[-1]
    
    def latest_dist(self):
        if self.dist == []:
            return None
        return self.dist[-1]
        
class Config:
    def __init__(self, bias, config_list, comment="no comment", check_tri=True):
        assert((bias == np.sort(bias)).all())
        self.comment=comment
        self.bias = bias
        if check_tri:
            self.configs = [(i,j,bij,value(bias[i],bias[j],bij)) for (i,j,bij) in config_list if check(bias[i],bias[j],bij)]
        else: 
            self.configs = [(i,j,bij,value(bias[i],bias[j],bij)) for (i,j,bij) in config_list]
        self.comp = [c for (a,b,d,c) in self.configs]
        self.F_param = []
        self.F_eval = []
        self.hist_F = ApproxHistory()
        self.hist_hard = ApproxHistory()
        self.F_aux = []
        
    def add_config(self,i,j,bij):
        if check(self.bias[i],self.bias[j],bij):
            print("Adding %.3f %.3f %.3f"%(self.bias[i], self.bias[j], bij))
            v = value(self.bias[i],self.bias[j],bij)
            self.configs.append((i,j,bij,v))
            self.comp.append(v)
            self.F_eval = [self.test_all_pairs(p) for p in self.F_param] #wasteful, but not too bad...
            return True
        else:
            return False
        
    def pkl(self):
        return dill.dumps(self)
    
    def backup(self,filename):
        f = open(filename, "wb")
        dill.dump(self,f)
        f.close()
    
    @classmethod
    def restore(cls,filename):
        f = open(filename, "rb")
        inst = dill.load(f)
        f.close()
        return inst
    
    def add_tlist(self,tlist,aux=[]):
        self.F_param.append(np.array(tlist))
        self.F_aux.append(aux)
        self.F_eval.append(self.test_all_pairs(np.array(tlist)))
    
    def add_function(self,f,aux=[]):
        self.add_tlist([f(b) for b in self.bias],aux=aux)
        
    def test_pairs(self,tlist):
        """
        tlist is the proposed threshold function
        """
        res = [thresh_eval(self.bias[i], self.bias[j], bij, tlist[i], tlist[j]) / c for (i,j,bij,c) in self.configs]
        return np.min(res)

    def test_all_pairs(self,tlist):
        return np.array([thresh_eval(self.bias[i], self.bias[j], bij, tlist[i], tlist[j]) for (i,j,bij,c) in self.configs])

    def test_pairs_weighted(self,tlist, pw):
        """
        tlist is the proposed threshold function
        pw is the weights of the pairs
        """
        res = self.test_all_pairs(tlist)
        return ((np.dot(res, pw) / np.dot(self.comp, pw)), res)

    def test_pairs_weighted_sparse(self, tlist, pw):
        """
        tlist is the proposed threshold function
        pw is the weights of the pairs
        """
        k = len(pw)
        ressum = 0
        compsum = 0
        for l in range(k):
            if pw[l] > 0:
                (i,j,bij,c) = self.configs[l]
                ressum += thresh_eval(self.bias[i], self.bias[j], bij, tlist[i], tlist[j]) * pw[l]
                compsum += c * pw[l]
        return  ressum / compsum
    
    def current_approx_ratio(self, clean=False):
        # compute the best distribution of rounding functions
        N = len(self.F_eval)
        k = len(self.configs)
        p = cp.Variable(N)
        
        constraints = [p >= 0]
        constraints.append(sum(p[i] for i in range(N)) == 1)

        obj = cp.Variable()
        for i in range(k):
            constraints.append(sum(self.F_eval[j][i]*p[j] for j in range(N)) >= obj * self.comp[i])

        prob = cp.Problem(cp.Maximize(obj), constraints)
        try:
            prob.solve()
        except:
            print("SOMETHING FAILED")
            return 0
        
        self.hist_F.insert(obj.value, [p.value[i] for i in range(N)])
        
        if clean:
            good = [i for i in range(N) if p.value[i] > 1e-9]
            self.F_eval = [self.F_eval[i] for i in good]
            self.F_param = [self.F_param[i] for i in good]
            self.F_aux = [self.F_aux[i] for i in good]
            
        return obj.value
        
    
    # linear program which computes next hard distribution of instances
    def next_config_dist(self, verbose=2):
        k = len(self.configs)
        if self.F_eval == []:
            return [1] + [0 for i in range(k-1)], 0
        
        obj = cp.Variable()
        ww = cp.Variable(k)
        #obj = cp.Variable()
        constraints = []
        for i in range(k):
            constraints.append(ww[i] >= 0)
        constraints.append(sum(self.comp[j] * ww[j] for j in range(k)) == 1)

        for a in self.F_eval:
            constraints.append(sum(a[i] * ww[i] for i in range(k)) <= obj)

        prob = cp.Problem(cp.Minimize(obj), constraints)
        try:
            prob.solve(cp.ECOS)
        except:
            print("SOMETHING FAILED")
            return [1] + [0 for i in range(k-1)], 0

        obj = obj.value
        cur_w = ww.value
        cur_w = [max(a,0) for a in cur_w]
        cur_w = np.array(cur_w) / np.sum(cur_w)
        
        self.hist_hard.insert(obj, cur_w)
        if verbose >= 1:
            print("HARDEST POINT:", obj)
        # sort by weights
        b1 = [self.bias[a] for (a,b,c,d) in self.configs]
        b2 = [self.bias[b] for (a,b,c,d) in self.configs]
        if verbose >= 2:
            plt.title("HARDEST POINT: %f"%obj)        
            plt.scatter(b1,b2,s=[len(self.bias)*40*x for x in cur_w])
            plt.axis("square")
        
        if verbose >= 1:
            print("Weight, (  b1  ,   b2  ,  b12  , value)")
            ll = [(cur_w[i], self.configs[i]) for i in range(k)]
            ll.sort()
            
            for (w,(a,b,c,d)) in reversed(ll):
                if w > 1e-4:
                    print("%0.06f, (%+0.03f, %+0.03f, %+0.03f, %0.03f)"%\
                          (w,self.bias[a],self.bias[b],c,d))
        return cur_w, obj
    
    def plot_best_distribution(self, filename, top=5, ignore=1e-4):
        N = len(self.F_eval)
        apx = self.hist_F.latest_approx()
        p = self.hist_F.latest_dist()
        l = [(p[i],i) for i in range(N)]
        l.sort()
        l = [x for x in reversed(l)]
        for p,i in l[:top]:
            print("Probability:", p)
            print("Function: ", self.F_param[i])
            plt.plot(self.bias, self.F_param[i], label="%.3f"%p)
        plt.legend(bbox_to_anchor=(1,1), loc="upper left")
        plt.title("Ratio: %.5f"%apx)
        plt.xlabel("bias")
        plt.ylabel("threshold")
        plt.savefig(filename)
