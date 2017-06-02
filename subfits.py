#!/usr/bin/env python
import pyfits as pf, numpy as np
from sys import argv

def main(qf, uf, of):
    h = pf.getheader(qf)
    q = pf.getdata(qf)
    u = pf.getdata(uf)
    p = q - u
    pf.writeto(of, p, h, clobber=True)
    
if __name__=='__main__':
    main(argv[1], argv[2], argv[3])
