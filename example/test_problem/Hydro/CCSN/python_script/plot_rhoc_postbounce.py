#! /usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  Purpose:
#    Compare the evolution of central density between FLASH and GAMER
#    for Post Bounce test problem
#

import matplotlib
matplotlib.use("Agg")

import os
import re
import numpy as np
import matplotlib.pyplot as plt


# dictionary for the corresponding reference solution
path_RefSol = "../ReferenceSolution/PostBounce"
RefSol      = {"s20GREP_LS220_15ms_none": "flash_1d_s20GREP_LS220_none.dat",
               "s20GREP_LS220_15ms_LB1" : "flash_1d_s20GREP_LS220_15ms_LB1.dat",
               "s20GREP_LS220_15ms_LB2" : "flash_1d_s20GREP_LS220_15ms_LB2.dat",
               "s20GREP_LS220_15ms_LB3" : "flash_1d_s20GREP_LS220_15ms_LB3.dat" }

# retrieve the filename of input profile in Input__TestProb
par_testprob   = open("../Input__TestProb").read()
CCSN_Prof_File = re.findall(r"CCSN_Prof_File\s*(\S*)\s", par_testprob)
CCSN_Prof_File = os.path.basename(CCSN_Prof_File[0])


### load data
## FLASH
fn_flash = RefSol[CCSN_Prof_File]
fn_flash = os.path.join(path_RefSol, fn_flash)
print("Filename of initial profile: {}".format(CCSN_Prof_File))
print("Adopted reference solution : {}".format(fn_flash))

# retrieve bounce time
with open(fn_flash) as f:
    tb_flash = f.readline()

tb_flash = re.findall(r"\d\.\d+", tb_flash)
tb_flash = float(tb_flash[0])

t_flash, rhoc_flash = np.genfromtxt(fn_flash, usecols = [0, 16], unpack = True)
t_flash -= tb_flash

## GAMER
fn_gamer = "../Record__CCSN_Properties"
t_gamer, rhoc_gamer = np.genfromtxt(fn_gamer, usecols = [0, 5], unpack = True)
t_gamer += 0.015  # the IC is from the FLASH simulation at t_pb = 15 ms


### plot the evolution of central density
fig, ax = plt.subplots()
alpha = 0.8

ax.plot(t_gamer * 1e3, rhoc_gamer / 1e14, label = "GAMER", alpha = alpha)
ax.plot(t_flash * 1e3, rhoc_flash / 1e14, label = "FLASH", alpha = alpha)

ax.set_xlabel(r"$t_{pb}$ [ms]")
ax.set_ylabel(r"Central density [$10^{14}$ g cm$^{-3}$]")
ax.set_ylim(3.4, 3.7)
ax.legend(loc = "upper right")

fig.tight_layout()
plt.savefig("Rhoc_PostBounce.png", dpi = 200)

plt.close(fig)
