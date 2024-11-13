import cmath
import math

import astropy
from astropy.units import Quantity
from hapsira.bodies import Earth
from hapsira.util import norm


class Mechanics:
    def __init__(self):
        pass

    @staticmethod
    def hyperbolic_period(a: float, e: float, n1: float, n2: float) -> float:
        """

        :param a:
        :param e:
        :param n1:
        :param n2:
        :return:
        """
        if n2 < n1:
            n1, n2 = n2, n1

        # get the hyperbolic eccentric anomaly
        F1 = 2 * math.atanh(math.sqrt((e - 1) / (e + 1)) * math.tan(n1 / 2))
        F2 = 2 * math.atanh(math.sqrt((e - 1) / (e + 1)) * math.tan(n2 / 2))

        # get the mean anomaly
        M1 = e * math.sinh(F1) - F1
        M2 = e * math.sinh(F2) - F2

        return math.sqrt(-pow(a, 3) / Earth.k.value) * (M2 - M1)

    @staticmethod
    def elliptic_period(a: float, e: float, n1: float, n2: float) -> float:
        """
        Returns the orbital period of the Hohmann Transfer.

        :param n2: Final true anomaly
        :param n1: Initial true anomaly
        :param e: Eccentric
        :param a: Semi-major axis
        :return: Total time for the transfer
        """

        E1 = math.atan(math.sin(n1) * math.sqrt(1 - e * e) / (math.cos(n1) + e))
        M1 = E1 - e * math.sin(E1)

        E2 = math.atan(math.sin(n2) * math.sqrt(1 - e * e) / (math.cos(n2) + e))
        M2 = E2 - e * math.sin(E2)

        return abs(math.sqrt(pow(a, 3) / Earth.k.value) * (M2 - M1))

    @staticmethod
    def vel(a: float, R: float):
        """
        Returns the orbital velocity at distance R.

        :param a: semimajor axis
        :param R: distance from body
        :return: velocity magnitude
        """
        return math.sqrt(Earth.k.value * (2 / R - 1 / a))

    @staticmethod
    def vec_in_dir(mag: float, direction: astropy.units.Quantity):
        """
        Returns a vector in the direction given by direction.

        :param mag: Magnitude of the vector
        :param direction: direction of the vector
        :return: velocity vector
        """

        # get norm
        dir_norm = norm(direction)
        unit_r = [-x / dir_norm for x in direction]

        return [(mag * x).value for x in unit_r]
