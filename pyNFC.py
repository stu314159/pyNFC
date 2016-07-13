#pyNFC.py
"""
Implementation file for pyNFC library module -- under development

"""

import FluidChannel as fc  #<-- provide geometric domain


class nfcGeometry:
    """
     geomtry object for pyNFC - basically a FluidChannel object with 
     selected obstruction

    """

    def __init__(self,geom = fc.FluidChannel()):
        """
          geom - FluidChannel object.  Default is an empty channel
                 with default dimensions and default EmptyChannel
                 "obstruction"
        """
