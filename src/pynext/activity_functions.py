"""
activity functions
"""

from math import pi, exp, log
from . system_of_units import *
import pandas as pd
from collections import namedtuple

# Cylindrical Vessel Activity (CVA)
CVA = namedtuple('CVA', """name
                           body_bi214 head_bi214
                           body_tl208 head_tl208""")
Activity = namedtuple('Activity', """name bi214 tl208""")

def activity_lsc_gammas_through_CV(name, cv, gamma_flux):
    """Returns the activity of gammas through a cylindrical vessel (cv)"""
    activity = CVA(name = name,
                     body_bi214 = gamma_flux.Bi214     * cv.body_surface,
                     head_bi214 = 2 * gamma_flux.Bi214 * cv.head_surface,
                     body_tl208 = gamma_flux.Tl208     * cv.body_surface,
                     head_tl208 = 2 * gamma_flux.Tl208 * cv.head_surface)
    return activity


def activity_gammas_transmitted_CV(name, cv, ia):
    """For a cylindrical vessel (cv) and an incoming activity (ia)
       compute the transmitted activity

    """
    activity = CVA(name = name,
                     body_bi214 = ia.body_bi214 * cv.body_transmittance,
                     head_bi214 = ia.head_bi214 * cv.body_transmittance,
                     body_tl208 = ia.body_tl208 * cv.body_transmittance,
                     head_tl208 = ia.head_tl208 * cv.body_transmittance)
    return activity


def activity_of_CV(name, cv):
    """Returns the self-shielded activity emanating from a CV"""
    activity = CVA(name = name,
                     body_bi214 = cv.body_self_shield_activity_bi214,
                     head_bi214 = cv.head_self_shield_activity_bi214,
                     body_tl208 = cv.body_self_shield_activity_tl208,
                     head_tl208 = cv.head_self_shield_activity_tl208)
    return activity

def punit(val, unit):
    u = eval(unit)
    return "{:7.2f} {:s}".format(val / u, unit)

def print_activity_of_CV(act, unit='Bq'):

    print("""
    activity \t\t {}
    body  (Bi-214) \t {}
    head  (Bi-214) \t {}
    total (Bi-214) \t {}
    body  (Tl-208) \t {}
    head  (Tl-208) \t {}
    total (Tl-208) \t {}
    """.format(act.name,
           punit(act.body_bi214,unit), punit(act.head_bi214,unit),
           punit(act.body_bi214 + act.head_bi214, unit),
           punit(act.body_tl208,unit), punit(act.head_tl208,unit),
           punit(act.body_tl208 + act.head_tl208,unit)))


def activity_table(activities):

    df = pd.DataFrame(activities, columns=activities[0]._fields)
    #df2 = df[['name','body_bi214', 'head_bi214', 'body_tl208','head_tl208']].copy()
    names = [df[column].name for column in df]
    for name in names[1:]:
        df[name] /= mBq
    return df


def pmt_activity(name, nof_pmt, PMT):
    solid_angle_window = 0.5 * 0.9
    solid_angle_pmt    = 0.5 * 0.35
    solid_angle_base   = 0.5 * 0.35
    activity = Activity(name = name,
                     bi214 = nof_pmt * (PMT.a_window_bi214 * solid_angle_window +
                                        PMT.a_base_bi214 * solid_angle_base +
                                        PMT.a_pmt_bi214 * solid_angle_pmt),
                     tl208 = nof_pmt * (PMT.a_window_tl208 * solid_angle_window +
                                        PMT.a_base_tl208 * solid_angle_base +
                                        PMT.a_pmt_tl208 * solid_angle_pmt))
    return activity


def sipm_activity(name, nof_sipm, SiPM):
    solid_angle = 0.5
    activity = Activity(name = name,
                     bi214 = nof_sipm * SiPM.a_bi214 * solid_angle,
                     tl208 = nof_sipm * SiPM.a_tl208 * solid_angle)
    return activity


def print_activity(name, act, unit='Bq'):

    print("""
    activity \t {}
    Bi-214 \t {}
    Tl-208 \t {}
    """.format(act.name,
           punit(act.bi214,unit),
           punit(act.tl208,unit)))


def str_activity(name, act, unit='Bq'):

    s = """
    activity \t {}
    Bi-214 \t {}
    Tl-208 \t {}
    """.format(act.name,
           punit(act.bi214,unit),
           punit(act.tl208,unit))
    return s
