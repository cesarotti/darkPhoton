# This file was automatically created by FeynRules 2.0.30
# Mathematica version: 8.0 for Microsoft Windows (64-bit) (October 6, 2011)
# Date: Mon 13 Apr 2015 17:06:05


from __future__ import division
from object_library import all_particles, Particle
import parameters as Param

import propagators as Prop

a = Particle(pdg_code = 22,
             name = 'a',
             antiname = 'a',
             spin = 3,
             color = 1,
             mass = Param.ZERO,
             width = Param.ZERO,
             texname = 'a',
             antitexname = 'a',
             charge = 0,
             LeptonNumber = 0)

dP = Particle(pdg_code = 9000001,
              name = 'dP',
              antiname = 'dP',
              spin = 3,
              color = 1,
              mass = Param.MdP,
              width = Param.WdP,
              texname = 'dP',
              antitexname = 'dP',
              charge = 0,
              LeptonNumber = 0)

e__minus__ = Particle(pdg_code = 11,
                      name = 'e-',
                      antiname = 'e+',
                      spin = 2,
                      color = 1,
                      mass = Param.Me,
                      width = Param.ZERO,
                      texname = 'e-',
                      antitexname = 'e+',
                      charge = -1,
                      LeptonNumber = 1)

e__plus__ = e__minus__.anti()

dF__minus__ = Particle(pdg_code = 9000002,
                       name = 'dF-',
                       antiname = 'dF+',
                       spin = 2,
                       color = 1,
                       mass = Param.MdF,
                       width = Param.ZERO,
                       texname = 'dF-',
                       antitexname = 'dF+',
                       charge = 0,
                       LeptonNumber = 0)

dF__plus__ = dF__minus__.anti()

