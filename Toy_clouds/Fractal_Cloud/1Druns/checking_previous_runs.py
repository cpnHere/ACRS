#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Oct  6 10:03:48 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To avoid run over already generated files.
"""

import sys,os


outname=str(sys.argv[1])
if os.path.isfile(outname):
    sys.exit('1')
else:
    sys.exit('0')