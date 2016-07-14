#pyNFC_preprocess.py
"""
   pre-processing steps for pyNFC code.  Basic geometry definition
   and partitioning information created and input file created prior to
   parallel simulation run

"""

import FluidChannel as fc
import pyNFC 
import pyPartition as pp
from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected


