# This file was automatically created by FeynRules 2.0.30
# Mathematica version: 8.0 for Microsoft Windows (64-bit) (October 6, 2011)
# Date: Mon 13 Apr 2015 17:06:05



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
gedP = Parameter(name = 'gedP',
                 nature = 'external',
                 type = 'real',
                 value = 0.0001,
                 texname = '\\text{gedP}',
                 lhablock = 'DARKPHOTONINPUTS',
                 lhacode = [ 1 ])

gdFdP = Parameter(name = 'gdFdP',
                  nature = 'external',
                  type = 'real',
                  value = 0.01,
                  texname = '\\text{gdFdP}',
                  lhablock = 'DARKPHOTONINPUTS',
                  lhacode = [ 2 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

MdP = Parameter(name = 'MdP',
                nature = 'external',
                type = 'real',
                value = 0.0001,
                texname = '\\text{MdP}',
                lhablock = 'MASS',
                lhacode = [ 9000001 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.000511,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MdF = Parameter(name = 'MdF',
                nature = 'external',
                type = 'real',
                value = 0.000511,
                texname = '\\text{MdF}',
                lhablock = 'MASS',
                lhacode = [ 9000002 ])

WdP = Parameter(name = 'WdP',
                nature = 'external',
                type = 'real',
                value = 1.e-8,
                texname = '\\text{WdP}',
                lhablock = 'DECAY',
                lhacode = [ 9000001 ])

eEM = Parameter(name = 'eEM',
                nature = 'internal',
                type = 'real',
                value = '3.5449077018110318*(1/aEWM1)**0.5',
                texname = 'e')

