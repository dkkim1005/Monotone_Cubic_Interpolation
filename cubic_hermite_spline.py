#!/usr/bin/env python3
import numpy as np

class cubicHermiteSpline:
    def __init__ (self, xvec, yvec, mvec):
        self._xvec = np.array(xvec)
        self._yvec = np.array(yvec)
        self._mvec = np.array(mvec)
        self._interp_f = lambda t, idx : (2.*t**3 - 3.*t**2 + 1.)*self._yvec[idx] +\
            (t**3 - 2*t**2 + t)*(self._xvec[idx+1] - self._xvec[idx])*self._mvec[idx] +\
            (-2*t**3 + 3*t**2)*self._yvec[idx + 1] + (t**3 - t**2)*(self._xvec[idx+1] -\
            self._xvec[idx])*self._mvec[idx+1]

    def __call__ (self, x):
        idx = self._binary_search(x)
        if idx == len(self._xvec) - 1: return self._yvec[-1]
        t = (x - self._xvec[idx])/(self._xvec[idx+1] - self._xvec[idx])
        return self._interp_f(t, idx)

    def _binary_search (self, x):
        assert(self._xvec[0] <= x and x <= self._xvec[-1])
        idx_l = 0
        idx_r = len(self._xvec) - 1
        idx = int(len(self._xvec)/2)
        boundary_check = lambda idx : (self._xvec[idx] <= x and x < self._xvec[idx+1])
        while(True):
            if (idx_r - idx_l) == 1:
                if boundary_check(idx): return idx
                else: return idx + 1
            if boundary_check(idx): return idx
            elif self._xvec[idx + 1] <= x:
                idx_l = idx
                idx = int((idx_r - idx_l)/2 + idx_l)
            else: # x < xVec[idx]
                idx_r = idx
                idx = int((idx_r - idx_l)/2 + idx_l)

class mCubicInterp (cubicHermiteSpline):
    def __init__(self, xVec, yVec):
        assert(len(xVec) == len(yVec))
        size = len(xVec)
        delta = np.zeros([size])
        mVec = np.zeros_like(delta)
        for i in range(size-1): delta[i] = (yVec[i+1] - yVec[i])/(xVec[i+1] - xVec[i])
        for i in range(1,size-1): mVec[i] = (delta[i-1] + delta[i])/2.
        mVec[0] = delta[0]
        mVec[size-1] = delta[size-2]
        for i in range(size-1):
            if np.abs(delta[i]) < 1e-30:
                mVec[i] = 0; mVec[i+1] = 0
        super(__class__, self).__init__(xVec, yVec, mVec)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    x0 = np.linspace(1, np.pi*5, 30)
    y0 = np.sin(x0)*x0
    interpolator = mCubicInterp(x0, y0)
    x1 = np.linspace(1, np.pi*5, 300)
    y1 = np.array([interpolator(x) for x in x1])
    # plot the results on the graph through matplotlib
    plt.plot(x0, y0, marker = 'o', color = 'r', label = 'Original', markersize = 5, linestyle = 'None', zorder = 3)
    plt.plot(x1, y1, marker = 'o', color = 'b', label = 'Monotone Cubic', markersize = 3.5, linestyle = 'None')
    fontsize = 15
    plt.xlabel(r'$x$', fontsize = fontsize)
    plt.ylabel(r'$y$', fontsize = fontsize)
    plt.legend(loc = 'upper left', fontsize = fontsize, edgecolor = 'None')
    plt.tick_params(labelsize = fontsize)
    plt.text(1, -8.3, r'$y=x\mathrm{sin}(x)$', fontsize = fontsize+3)
    plt.show()
    #plt.savefig('mcubic.pdf', dpi = 100, bbox_inches = 'tight')
