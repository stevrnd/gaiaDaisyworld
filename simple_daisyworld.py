# -*- coding: utf-8 -*-
"""

Module  : CMP-6013Y - CMP Third Year Project
File    : simple_daisyworld.py
Date    : Friday 25 December 2020
Desc.   : Modelling the Gaia hypothesis by implementing a the version James Lovelock, it's purpose is to lay out
          the foundations to further develop Daisyworld.
History : 25/12/2020 - v1.0 - Created project file
          30/12/2020 - v1.1 - Added methods derived from equations from Watson's and Lovelock's published paper.

"""

__author__ = "Steven Diep"
__maintainer__ = "Steven Diep"
__email__ = "steven_diep@hotmail.co.uk"
__status__ = "Prototype"  # "Development" "Prototype" "Production"

import math


def uncolonised_ground(area_b, area_w):
    """Calculates fractional area of uncolonised ground

    :param double area_b: Local temperature experienced by the daisies
    :param double area_w: Optimum temperature for black and white daisies to thrive at

    :rtype: double
    :return: Fractional area covered by bare ground
    """
    p = 1
    return p - area_b - area_w


def beta_y(temp_y, opt_temp=22.5, c=0.003265):
    """Daisy growth rate function

    :param double temp_y: Local temperature experienced by the daisies
    :param double opt_temp: Optimum temperature for black and white daisies to thrive at
    :param double c: Determines quadratic width, allowing growth to range from 5 to 40 degrees celsius on a negative parabolic curve

    :rtype: double
    :return: A value between 0 - 1, 0 indicating no growth, 1 indicating maximal growth
    """
    min_growth = round(opt_temp - math.sqrt(1 / c), 1)
    max_growth = round(opt_temp + math.sqrt(1 / c), 1)
    if min_growth <= temp_y <= max_growth:
        return 1 - c * (opt_temp - temp_y) ** 2
    else:
        return 0


def daisy_growth(area_y, area_g, beta, gamma):
    """Calculates rate of change either in black or white daisies

    :param double area_y: Fraction of surface area for black or white daisies
    :param double area_g: Fraction of surface area for bare/uncovered ground
    :param beta: Growth rate of black or white daisies
    :param gamma: Death rate of black or white daisies

    :rtype: double
    :return: Rate of change of daisies
    """
    return area_y * ((area_g * beta) - gamma)


def local_temp(a_d, a_y, temp_d, q=20):
    """Works out planetary albedo by working out and summing all allocated areas

    :param double a_d: Daisies albedo
    :param double a_y: Planetary albedo
    :param double temp_d: Average temperature of Daisyworld
    :param int q: Heat absorption coefficient

    :rtype: double
    :return: Temperature of local daisies
    """
    return q * (a_y - a_d) + temp_d


def planetary_albedo(area_b, area_w, area_g, a_b=0.25, a_w=0.75, a_g=0.5):
    """Works out planetary albedo by taking weighted average and summing all allocated areas

    :param double area_b: Fraction of surface area for black daisies
    :param double area_w: Fraction of surface area for white daisies
    :param double area_g: Fraction of surface area for bare/uncovered ground
    :param double a_b: Albedo of black daisy, defined as 0.25
    :param double a_w: Albedo of white daisy, defined as 0.75
    :param double a_g: Albedo of bare ground, defined as 0.5

    :rtype: double
    :return: Planetary albedo between 0 - 1, 1 is fully reflective, 0 is full black body
    """
    return area_b * a_b + area_w * a_w + area_g * a_g


def planetary_temp(flux, lumen, a_d):
    """Works out average planetary temperature

    :param double flux: Solar flux constant
    :param double lumen: Solar luminosity
    :param double a_d: Planetary albedo

    :rtype: double
    :return: Average planetary temperature
    """
    sigma = 5.67037e-8  # Stefan-Boltzmann constant
    abs_zero = 273.15   # Used to calculate temperature in celsius
    return (((flux * lumen * (1 - a_d)) / sigma) ** 0.25) - abs_zero


def solar_luminosity(time_step, sim_length):
    """Mimics a portion of a star's lifecycle

    :param int time_step: Arbitrary measure of time
    :param int sim_length: Length of simulation

    :rtype: double
    :return: A multiplier from range 0.6 - 1.2
    """
    return 0.6 + (time_step * (1/sim_length))
