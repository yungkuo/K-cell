# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 15:42:48 2014

@author: KyoungWon
"""
import numpy as np

def inside_polygon(x, y, points):
    """
    Return True if a coordinate (x, y) is inside a polygon defined by
    a list of verticies [(x1, y1), (x2, x2), ... , (xN, yN)].

    Reference: http://www.ariel.com.au/a/python-point-int-poly.html
    """
    n = len(points)
    inside = 0
    p1x, p1y = points[0]
    for i in range(1, n + 1):
        p2x, p2y = points[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def mean_polygon(mov, pts):
    #frame=len(mov[:,0,0])
    nrow=len(mov[:,0])
    ncol=len(mov[0,:])
    mask = np.zeros((nrow, ncol), dtype=np.int)
    for i in range(nrow):
        for j in range(ncol):
            mask[i,j]=inside_polygon(j,i,pts)

    #mask3d = np.tile(mask, (frame,1,1))
    maskedimg = np.multiply(mov,mask)
    mean = np.sum(maskedimg)/mask.sum()
    return mean, maskedimg
