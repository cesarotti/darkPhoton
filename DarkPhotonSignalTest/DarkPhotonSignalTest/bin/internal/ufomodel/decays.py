# This file was automatically created by FeynRules 2.0.30
# Mathematica version: 8.0 for Microsoft Windows (64-bit) (October 6, 2011)
# Date: Mon 13 Apr 2015 17:06:05


from object_library import all_decays, Decay
import particles as P


Decay_dP = Decay(name = 'Decay_dP',
                 particle = P.dP,
                 partial_widths = {(P.dF__minus__,P.dF__plus__):'((8*gdFdP**2*MdF**2 + 4*gdFdP**2*MdP**2)*cmath.sqrt(-4*MdF**2*MdP**2 + MdP**4))/(48.*cmath.pi*abs(MdP)**3)',
                                   (P.e__minus__,P.e__plus__):'((4*gedP**2*MdP**2 + 8*gedP**2*Me**2)*cmath.sqrt(MdP**4 - 4*MdP**2*Me**2))/(48.*cmath.pi*abs(MdP)**3)'})

