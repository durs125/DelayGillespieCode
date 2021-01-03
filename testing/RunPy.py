#!/usr/bin/env python

PYTHONOPTIMIZE = "OO"
import pandas as pd
import os
import math
import time
import py_compile
import multiprocessing as mp
py_compile.compile("Classes_Gillespie.py", "__pycache__/Classes_Gillespie.pyc",optimize =2)
py_compile.compile("Functions_Gillespie.py", "__pycache__/Functions_Gillespie.pyc",optimize =2)
py_compile.compile("main.py", "__pycache__/main.pyc",optimize =2)
py_compile.compile("PostProcessing_Functions.py", "__pycache__/PostProcessing_Functions.pyc",optimize =2)
os.system('python __pycache__/main.pyc')
exit()

