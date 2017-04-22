#!/usr/bin/python2.7
import numpy as np

class CubicHermiteSpline:
    def __init__(self, xVec, yVec, mVec):
        self.xVec = np.array(xVec)
        self.yVec = np.array(yVec)
        self.mVec = np.array(mVec)

        self.interp_f = lambda t, idx : (2.*t**3 - 3.*t**2 + 1.)*self.yVec[idx] +\
                                        (t**3 - 2*t**2 + t)*(self.xVec[idx+1] -\
                                        self.xVec[idx])*mVec[idx] +\
                                        (-2*t**3 + 3*t**2)*self.yVec[idx + 1] +\
                                        (t**3 - t**2)*(self.xVec[idx+1] -\
                                        self.xVec[idx])*self.mVec[idx+1]


    def __call__(self, x):
        idx = self._binary_search(x)

        if idx == len(self.xVec) - 1:
            return self.yVec[-1]

        t = (x - self.xVec[idx])/(self.xVec[idx+1] - self.xVec[idx])

        return self.interp_f(t, idx)



    def _binary_search(self, x):
        assert(self.xVec[0] <= x and x <= self.xVec[-1])

        idx_l = 0
        idx_r = len(self.xVec) - 1
        idx = len(self.xVec)/2

        is_x_in_Boundary = lambda idx : self.xVec[idx] <= x and x < self.xVec[idx + 1]

        while(True):
            if idx_r - idx_l == 1:
                if is_x_in_Boundary(idx):
                    return idx
                else:
                    return idx + 1
            if is_x_in_Boundary(idx):
                return idx
            elif self.xVec[idx + 1] <= x:
                idx_l = idx
                idx = (idx_r - idx_l)/2 + idx_l
            else: # x < xVec[idx]
                idx_r = idx
                idx = (idx_r - idx_l)/2 + idx_l



class MonotoneCubicInterpolation(CubicHermiteSpline):
    def __init__(self, xVec, yVec):
        assert(len(xVec) == len(yVec))
        size = len(xVec)

        delta = np.zeros([size])
        mVec = np.zeros([size])

        for i in xrange(size-1):
            delta[i] = (yVec[i+1] - yVec[i])/(xVec[i+1] - xVec[i])

        for i in xrange(1,size-1):
            mVec[i] = (delta[i-1] + delta[i])/2.

        mVec[0] = delta[0]; mVec[size-1] = delta[size-2]

        for i in xrange(size-1):
            if np.abs(delta[i]) < 1e-30:
                mVec[i] = 0; mVec[i+1] = 0

        CubicHermiteSpline.__init__(self, xVec, yVec, mVec)



if __name__ == "__main__":
    
    x = np.linspace(0,np.pi*5,30)
    y = 1./(x+1e-10) - 2/(x+1e-10)**2 + 3/(x+1e-10)**3

    tester = MonotoneCubicInterpolation(x, y)

    with open('origin.out','w') as f:
        for i,x_ in enumerate(x):
            f.write('%f %f\n'%(x_, y[i]))
        
    x1 = np.linspace(0,np.pi*5,300)
    with open('tester.out','w') as f:
        for i,x1_ in enumerate(x1):
            f.write('%f %f\n'%(x1_, tester(x1_)))
