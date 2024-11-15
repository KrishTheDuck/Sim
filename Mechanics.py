import math

import astropy
import numpy as np
import astropy.units as u
from astropy.units import Quantity
from hapsira.bodies import Earth
from hapsira.twobody import Orbit
from hapsira.util import norm


class Mechanics:
    def __init__(self):
        raise "Utility Class Cannot be Instantiated."

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
        :param e: Eccentricity
        :param a: Semi-major axis
        :return: Total time for the transfer
        """

        E1 = math.atan(math.sin(n1) * math.sqrt(1 - e * e) / (math.cos(n1) + e))
        M1 = E1 - e * math.sin(E1)

        E2 = math.atan(math.sin(n2) * math.sqrt(1 - e * e) / (math.cos(n2) + e))
        M2 = E2 - e * math.sin(E2)

        print(M1, M2, " n1,n2 ", n1, n2)

        # T^2 = pi a b sqrt(2/u * (1/r1 + 1/r2)
        return abs(2 * math.sqrt(pow(a, 3) / Earth.k.value) * np.radians(n2 - n1))

    @staticmethod
    def elliptical_period(a: float, e: float, r1: float, r2: float):
        b = a * math.sqrt(1 - e ** 2)
        return 0.5 * math.pi * a * b * math.sqrt((2 / Earth.k.value) * (1 / r1 + 1 / r2))

    @staticmethod
    def vel_from_orbit(orb: Orbit, nu: float):
        r, v = orb.rv()
        r = r.to(u.m)
        v = v.to(u.m / u.s)

        h = np.linalg.norm(Mechanics.__cross(r.value, v.value))

        r = (h ** 2 / Earth.k.value) / (1 + orb.ecc.value * np.cos(nu))

        return Mechanics.vel(orb.a.value * 1000, r)

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

    @staticmethod
    def __cross(A: np.ndarray, B: np.ndarray) -> np.ndarray:
        C_x = A[1] * B[2] - A[2] * B[1]
        C_y = A[2] * B[0] - A[0] * B[2]
        C_z = A[0] * B[1] - A[1] * B[0]
        return np.array([C_x, C_y, C_z])

    @staticmethod
    def orbital_elements(r, v, mu=398600.4418):  # mu is Earth's gravitational parameter in km^3/s^2
        r = np.array(r)
        v = np.array(v)

        # Magnitude of position and velocity
        r_norm = np.linalg.norm(r)
        v_norm = np.linalg.norm(v)

        # Angular momentum vector
        h = Mechanics.__cross(r, v)
        h_norm = np.linalg.norm(h)

        # Inclination
        inc = np.arccos(h[2] / h_norm)

        # Node vector
        n = Mechanics.__cross(np.array([0, 0, 1]), h)
        n_norm = np.linalg.norm(n)

        # Eccentricity vector
        e = ((v_norm ** 2 - mu / r_norm) * r - np.dot(r, v) * v) / mu
        ecc = np.linalg.norm(e)

        # Right Ascension of the Ascending Node (RAAN)
        if n_norm > 1e-10:  # Avoid division by zero
            raan = np.arccos(n[0] / n_norm)
            raan = 2 * np.pi - raan if n[1] < 0 else raan
        else:
            raan = 0.0

        # Argument of Periapsis
        if n_norm > 1e-10 and ecc > 1e-10:
            arg_pe = np.arccos(np.dot(n, e) / (n_norm * ecc))
            arg_pe = 2 * np.pi - arg_pe if e[2] < 0 else arg_pe
        else:
            arg_pe = 0.0

        # True Anomaly
        if ecc > 1e-10:
            true_anomaly = np.arccos(np.dot(e, r) / (ecc * r_norm))
            true_anomaly = 2 * np.pi - true_anomaly if np.dot(r, v) < 0 else true_anomaly
        else:
            cp = Mechanics.__cross(n, r)
            true_anomaly = np.arccos(np.dot(n, r) / (n_norm * r_norm))
            true_anomaly = 2 * np.pi - true_anomaly if cp[2] < 0 else true_anomaly

        # Semi-major axis
        a = 1 / (2 / r_norm - v_norm ** 2 / mu)

        # Period (for elliptical orbits)
        period = 2 * np.pi * np.sqrt(abs(a ** 3) / mu) if a > 0 else None

        return {
            "a": a,
            "ecc": ecc,
            "inc": np.degrees(inc),
            "raan": np.degrees(raan),
            "argp": np.degrees(arg_pe) + 180,
            "nu": np.degrees(true_anomaly),
            "T": period
        }
