#!/bin/bash

# this script uses Brett's nuosc code at  https://github.com/brettviren/nuosc
# usage: ./nuosc_bv [energy_range] [nu_flavor] [delta_cp] [mass_hierarchy]

energy="${1:-lin:0.3,6.3,0.02}"
neutrino=${2:-2}
delta=${3:-0}
mh=${4:-}

./nuosc \
--neutrino $neutrino \
--delta:$delta \
--atm ${mh}2.41e-3 \
--sol 7.59e-5 \
--mixing sin:0.84,1.0,0.09 \
-e $energy \
-b one:1300 \
-D con:2.7

#==================== Output =====================================
# <energy> <baseline> <P(nux->nue)> <P(nux->numu)> <P(nux->nutau)>

#============ Options ====================================
# [-n|--neutrino <initial neutrino number>]
# 	Sets the initial neutrino e,mu,tau = 1,2,3; anti *= -1
# 	Default is 2 (mu)

# [-e|--energy <energy range description>]
# 	Describe the energy range.  Can be one of:
# 		Linear energy range:      "lin:start,stop,step"
# 		Logarithmic energy range: "log:start,stop,step"
# 		A single energy value:    "one:energy"
# 	Linear energy units are in GeV,
# 	Logarithmic values give 10^start GeV - 10^stop GeV,
# 	Default is "one:1"

# [-b|--baseline <baseline range description>]
# 	Describe the baseline range.  Can be one of:
# 		Linear disance range:       "lin:start,stop,step"
# 		Logarithmic distance range: "log:start,stop,step"
# 		A single distance value:    "one:distance"
# 		Zenith angle:               "zen:start,stop,step[,depth]"
# 		cos(Zenith angle):          "cos:start,stop,step[,depth]"
# 	Distance units are km, angle units are degrees.
# 	Logarithmic values give 10^start km - 10^stop km,
# 	The "depth" value for zenith related gives detector depth (def=0),
# 	Default is "one:1".

# [-m|--mixing <mixing angle description>]
# 	Set the mixing angles description.  Can be one of:
# 		Angles:        "ang:theta_12,theta_23,theta_13"
# 		Sin^2(2theta): "sin:sin^2(2t_12),sin^2(2t_23),sin^2(2t_13)"
# 	Angle units are in degrees, default "sin:0.8,1.0,0.1"

# [-s|--sol <delta m^2_solar (dm^2_21)>]
# 	Set the solar delta-m^2 (dm^2_21).
# 	Units are eV^2, default is 5.0e-5 eV^2

# [-a|--atm <delta m^2_atm (dm^2_31)>]
# 	Set the atmospheric delta-m^2 (dm^2_31).
# 	Units are eV^2, default is 2.5e-3 eV^2

# [-d|--delta <CP phase>]
# 	Set the CP phase angle.
# 	Units are degrees, default is 0 degrees

# [-D|--density <density description>]
# 	Set the density description.  Can be one of:
# 		Constant density: "con:density"
# 		PREM density:     "prem"
# 		Lookup table:     "lut:filename"
# 	Density is in g/cc.
# 	The lookup text file has columns of bin centered position and density
# 	Default is "con:0.0"

# [-c|--calculation <calculation description]
# 	Set the calculation description.  Can be one of:
# 		Matrix method: "matrix"
# 		Full stepping: "step"
# 	Default is "matrix"


