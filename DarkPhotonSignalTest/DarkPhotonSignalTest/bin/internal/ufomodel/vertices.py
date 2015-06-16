# This file was automatically created by FeynRules 2.0.30
# Mathematica version: 8.0 for Microsoft Windows (64-bit) (October 6, 2011)
# Date: Mon 13 Apr 2015 17:06:05


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.e__plus__, P.e__minus__, P.a ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

V_2 = Vertex(name = 'V_2',
             particles = [ P.dF__plus__, P.dF__minus__, P.dP ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_2})

V_3 = Vertex(name = 'V_3',
             particles = [ P.e__plus__, P.e__minus__, P.dP ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_3})

