#! /usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  Purpose:
#    Compare the evolution of central density between FLASH and GAMER
#    for GREP Migration test problem
#

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt


### load data
# FLASH
fn_flash = "../ReferenceSolution/MigrationTest/flash_1d.dat"
t_flash, rhoc_flash = np.genfromtxt(fn_flash, unpack = True)

# GAMER
fn_gamer = "../Record__CCSN_Properties"
t_gamer, rhoc_gamer = np.genfromtxt(fn_gamer, usecols = [0, 5], unpack = True)


### plot the evolution of central density
fig, ax = plt.subplots()
alpha = 0.8

ax.plot(t_gamer * 1e3, rhoc_gamer / 1e14, label = "GAMER", alpha = alpha)
ax.plot(t_flash * 1e3, rhoc_flash / 1e14, label = "FLASH", alpha = alpha)

ax.set_xlabel("Time [ms]")
ax.set_ylabel(r"Central density [$10^{14}$ g cm$^{-3}$]")
ax.legend(loc = "upper right")

fig.tight_layout()
plt.savefig("Rhoc_MigrationTest.png", dpi = 200)

plt.close(fig)
