from abc import abstractmethod, ABC
from copy import copy
from dataclasses import dataclass

import hapsira.core.thrust

import astropy.units as u
import numpy as np
from hapsira.core.thrust import change_a_inc, change_ecc_inc
from hapsira.maneuver import Maneuver
from hapsira.twobody import Orbit
from hapsira.util import norm


@dataclass(init=False, frozen=True)
class StandardUnits:
    Newton = u.Unit("kg * m / s**2")
    M_dot = u.Unit("kg / s")
    Density = u.Unit("kg / m**3")
    Time = u.Unit("s")
    Acceleration = u.Unit("m/s**2")
    Velocity = u.Unit("km/s")
    Distance = u.Unit("km")
    Mass = u.kg
    g0 = 9.81 * Acceleration


@dataclass(init=True, frozen=True, kw_only=True)
class AbstractEngine(ABC):
    name: str
    rho: StandardUnits.Density
    Isp: StandardUnits.Time

    @abstractmethod
    def mass_flow_rate(self, t) -> StandardUnits.M_dot:
        pass

    @abstractmethod
    def thrust_curve(self, t) -> StandardUnits.Newton:
        pass

    def __str__(self):
        return f"Engine {self.name}: Density<{self.rho}>, Specific Impulse<{self.Isp}>"


class AbstractSpacecraft(ABC):
    def __init__(self, engine: AbstractEngine,
                 structural_mass: StandardUnits.Mass,
                 fuel_mass: StandardUnits.Mass,
                 initial_orbit: Orbit):
        self.__engine = engine
        self.__structural_mass = structural_mass

        self.__fuel_mass = fuel_mass
        self.__mission_time = 0

        self.__current_orbit = copy(initial_orbit)
        self.__maneuvers = []

    def transfer(self, next_orbit: Orbit, transfer_at_nu: u.rad):
        # compute the maneuver and return it
        # dv = Isp * g0 * ln(m0/mf)
        # mf = m0 - mp
        # dm/dt * deltaT = deltaM => -Isp * g0 * ln((m0 - dm)/mf) = dv
        # I thrust in the direction of deltaV
        # state = self.__initial_orbit.v + self.__initial_orbit.r + [self.__fuel_mass]
        state = np.concatenate((
            self.__current_orbit.v.to_value(StandardUnits.Velocity),
            self.__current_orbit.a.to_value(StandardUnits.Acceleration),
            self.__fuel_mass.to_value(StandardUnits.Mass)
        ))



    def maneuver(self, man: Maneuver):
        m0 = (self.__fuel_mass + self.__structural_mass).to_value(StandardUnits.Mass)
        dv = man.get_total_cost().to_value(StandardUnits.Velocity)
        dm = m0 * (1 - np.exp(-dv / self.__engine.Isp * StandardUnits.g0))

        self.__fuel_mass -= dm
        self.__mission_time += man.get_total_time()

        self.__maneuvers.append(man)

    def total_orbit(self):
        return self.__current_orbit.apply_maneuver(
            Maneuver(*[(m[0], m[1]) for man in self.__maneuvers for m in man]),
            intermediate=True)

    