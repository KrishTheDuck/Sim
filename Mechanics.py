import math

import astropy
import astropy.units as u
import numpy as np
from hapsira.bodies import Earth
from hapsira.twobody import Orbit
from hapsira.util import norm


class Mechanics:
    def __init__(self):
        raise "Utility Class Cannot be Instantiated."

    @staticmethod
    def hyperbolic_period(a: float, e: float, n1: float, n2: float) -> float:
        """
        Computes the period of a hyperbolic orbit.

        :param a: semi-major axis
        :param e: eccentricity
        :param n1: initial true anomaly
        :param n2: final true anomaly
        :return: the time elapsed
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
        :return: Total time elapsed
        """
        E1 = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(n1 / 2))
        E2 = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(n2 / 2))
        n = math.sqrt(Earth.k.value / a ** 3)  # mean motion

        return (E2 - E1 - e * (math.sin(E2) - math.sin(E1))) / n

    @staticmethod
    def elliptical_period(a: float, e: float, r1: float, r2: float):
        """
        Calculates the elliptical period from position 1 to 2
        :param a: semi-major axis
        :param e: eccentricity
        :param r1: position 1
        :param r2: position 2
        :return: Total time elapsed
        """
        nu1 = math.acos(round((a * (1 - e ** 2) - r1) / (e * r1), 5))
        nu2 = math.acos(round((a * (1 - e ** 2) - r2) / (e * r2), 5))
        b = a * math.sqrt(1 - e ** 2)
        return abs(Mechanics.elliptic_period(a, e, nu1, nu2))

    @staticmethod
    def vel_from_orbit(orb: Orbit, nu: float):
        """
        Velocity at a position given an orbit and a true anomaly.
        :param orb: orbit
        :param nu: true anomaly
        :return: velocity
        """
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
        dir_norm = norm(direction).value
        unit_r = [x.value / dir_norm for x in direction]

        return [(mag * x) for x in unit_r]

    # helper function computes the cross product
    @staticmethod
    def __cross(A: np.ndarray, B: np.ndarray) -> np.ndarray:
        C_x = A[1] * B[2] - A[2] * B[1]
        C_y = A[2] * B[0] - A[0] * B[2]
        C_z = A[0] * B[1] - A[1] * B[0]
        return np.array([C_x, C_y, C_z])

    @staticmethod
    def orbital_elements(r, v):
        # G*M for the Earth
        mu = Earth.k.to_value(u.km ** 3 / u.s ** 2)

        # Position vector:
        r = np.array(r)
        r_norm = np.linalg.norm(r)

        # Velocity vector:
        v = np.array(v)
        v_norm = np.linalg.norm(v)

        # Angular momentum vector:
        h = Mechanics.__cross(r, v)
        h_norm = np.linalg.norm(h)

        # Inclination:
        inc = np.arccos(h[2] / h_norm)

        # Node vector:
        n = Mechanics.__cross(np.array([0, 0, 1]), h)
        n_norm = np.linalg.norm(n)

        # Eccentricity vector:
        e = Mechanics.__cross(v, h) / mu - r / r_norm
        ecc = np.linalg.norm(e)

        # Right Ascension of the Ascending Node (RAAN):
        if n_norm > 1e-10:
            raan = np.arccos(n[0] / n_norm)
            raan = 2 * np.pi - raan if n[1] < 0 else raan
        else:
            raan = 0.0

        # Argument of Periapsis:
        if n_norm > 1e-10 and ecc > 1e-10:
            arg_pe = np.arccos(np.dot(n, e) / (n_norm * ecc))
            arg_pe = 2 * np.pi - arg_pe if e[2] < 0 else arg_pe
        else:
            arg_pe = 0.0

        # True Anomaly:
        if ecc > 1e-10:
            true_anomaly = np.arccos(np.dot(e, r) / (ecc * r_norm))
            # Wrap Result
            true_anomaly = 2 * np.pi - true_anomaly if np.dot(r, v) < 0 else true_anomaly
        else:
            # Calculate the orientation of the orbit
            perp = Mechanics.__cross(n, r)
            true_anomaly = np.arccos(np.dot(n, r) / (n_norm * r_norm))
            true_anomaly = 2 * np.pi - true_anomaly if perp[2] < 0 else true_anomaly

        # Orbital Energy:
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
