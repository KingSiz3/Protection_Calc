# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 20:28:00 2020

@author: ASGRM
"""
def logger(func):
    def wrapped(*args, **kwargs):
        out = func(*args, **kwargs)
        with open("log.txt","a") as fhwr:
            fhwr.write(func.__name__+ ' args: '+str(args)+' kwargs: ' + str(kwargs) + ' return: '+str(out) + '\n')
        return out
    return wrapped