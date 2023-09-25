# -*- coding: utf-8 -*-
"""

Module  : CMP-6013Y - CMP Third Year Project
File    : main.py
Date    : Thursday 26 November 2020
Desc.   : Stores the constants used for daisyworld.py
History : 26/11/2020 - v1.0 - Created project file, added initial constants

"""
import math
import matplotlib.pyplot as plt
import numpy as np

import simple_daisyworld as simple
import daisyworld as enhanced
import dw_without_grey as enhanced_without_grey


def simple_main(albedo_b, albedo_w, death_type="default", growth_rate="default"):
    a_b = albedo_b
    a_w = albedo_w
    a_g = 0.5
    area_b = 0.01  # Fractional area covered by black daisies (System must start with a slight push for simulation to work)
    area_w = 0.01  # Fractional area covered by white daisies

    flux = 1050  # Rate of energy received in Watts per metre**2. Note: Value used is smaller than observed constant to simulate a younger star.
    # q = 20  # Heat absorption coefficient
    resolution = 10000
    if growth_rate is "high":
        c = 0.002
    elif growth_rate is "low":
        c = 0.013265
    elif growth_rate == "default":
        c = 0.003265

    sim_length = 550
    overtime_sun_intensity = []
    planet_temp = []
    planet_temp_d = []
    b_coverage = []
    w_coverage = []

    for step in range(sim_length):
        # Make sure black/white daisy does not go under 0.001 threshold
        if area_b < 0.01:
            area_b = 0.01
        if area_w < 0.01:
            area_w = 0.01
        if death_type is "plague":
            gamma = -0.15 * math.cos(step) + 0.45
        elif death_type is "default":
            gamma = 0.3  # Death rate of both daisies
        x = simple.uncolonised_ground(area_b, area_w)  # Fractional area covered by bare ground, x = p - a_b - a_w
        lumen = simple.solar_luminosity(step, sim_length)
        temp = simple.planetary_temp(flux, lumen, a_g)
        i = 0
        while i < resolution:
            a_d = simple.planetary_albedo(area_b, area_w, x, a_b=a_b, a_w=a_w, a_g=a_g)
            temp_d = simple.planetary_temp(flux, lumen, a_d)

            loc_temp_b = simple.local_temp(a_b, a_d, temp_d)
            beta_b = simple.beta_y(loc_temp_b, c=c)
            loc_temp_w = simple.local_temp(a_w, a_d, temp_d)
            beta_w = simple.beta_y(loc_temp_w, c=c)

            darea_bdt = simple.daisy_growth(area_b, x, beta_b, gamma)
            area_b = area_b + darea_bdt
            darea_wdt = simple.daisy_growth(area_w, x, beta_w, gamma)
            area_w = area_w + darea_wdt

            x = simple.uncolonised_ground(area_b, area_w)
            i += 1
        planet_temp_d.append(temp_d)
        planet_temp.append(temp)
        overtime_sun_intensity.append(lumen)
        b_coverage.append(area_b * 100)
        w_coverage.append(area_w * 100)

    plt.plot(overtime_sun_intensity, planet_temp_d, 'b', label='With daisies')
    plt.plot(overtime_sun_intensity, planet_temp, 'r', label='Without daisies')
    plt.legend(loc='upper right')
    plt.title('Temperature over luminosity')
    plt.xlabel('Solar Luminosity')
    plt.ylabel('Temperature (°C)')
    plt.show()

    plt.plot(overtime_sun_intensity, b_coverage, 'b', label='Black daisies')
    plt.plot(overtime_sun_intensity, w_coverage, 'g', label='White daisies')
    plt.legend(loc='upper right')
    plt.title('Area over luminosity')
    plt.xlabel('Solar Luminosity')
    plt.ylabel('Area (%)')
    plt.show()


def enhanced_main():
    x_dim = 50
    y_dim = 50
    # Do not recommend using this list,
    luminosities = [0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.69, 0.68, 0.67, 0.66, 0.65,
                    0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.74, 0.73, 0.72, 0.71, 0.70, 0.71,
                    0.72, 0.73, 0.74, 0.75, 0.74, 0.73, 0.72, 0.71, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77,
                    0.78, 0.79, 0.80]
    luminosities1 = np.arange(0.6, 1.4, 0.005)
    world = enhanced.Daisyworld(x_dim, y_dim, luminosities1, 350)
    world.run()


def enhanced_wout_grey_main():
    x_dim = 50
    y_dim = 50
    luminosities1 = np.arange(0.6, 1.4, 0.005)
    world = enhanced_without_grey.Daisyworld(x_dim, y_dim, luminosities1, 350)
    world.run()


if __name__ == "__main__":
    # simple_main(0.25, 0.75)
    # simple_main(0.4, 0.6)
    # simple_main(0.25, 0.75, death_type="plague")
    # simple_main(0.25, 0.75, growth_rate="high")
    # simple_main(0.25, 0.75, growth_rate="low")
    enhanced_main()
    # enhanced_wout_grey_main()
