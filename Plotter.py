import astropy.time
import astropy.units as u
import hapsira.frames
import numpy as np
from astropy.time import TimeDelta, Time
from hapsira.bodies import Earth
from hapsira.maneuver import Maneuver
from hapsira.plotting import OrbitPlotter
from hapsira.twobody import Orbit
from hapsira.util import time_range, norm
from matplotlib import pyplot as plt

from JPLQuery import JPLQuery
from Mechanics import Mechanics


class Plotter:
    colors = ["deepskyblue", "darkorange", "mediumseagreen", "gold", "blue", "green", "darkred", "cyan"]

    def __init__(self, body, location, start_epoch, stop_epoch, delta):
        self.body = body
        self.location = location
        self.s = JPLQuery(self.body, self.location, {
            'start': start_epoch.strftime("%Y-%m-%d 00:00:00"),
            'stop': stop_epoch.strftime("%Y-%m-%d 00:00:00"),
            'step': '1d'
        })
        self.vals = self.s.zoom_to_closest_approach()
        self.hor = self.s.get_ephem(attractor=Earth,
                                    epoch=time_range(start=self.vals['date'] - delta,
                                                     end=self.vals['date'] + delta))

        self.ephemerides = self.s.get_orbital_elements(
            epochs_time_range=time_range(start=self.vals['date'] - 0.5 * u.day,
                                         end=self.vals['date'] + 0.5 * u.day), EPOCH=self.vals['date'])

    @staticmethod
    def generate_orbit(ephemerides: dict, EPOCH: astropy.time.Time, r_min: float, r_max: float,
                       plane: hapsira.frames.Planes, true_anomaly: float):
        a = (r_min + r_max) / 2.0
        ecc = (r_max - r_min) / (r_max + r_min)

        return Orbit.from_classical(
            attractor=Earth,
            plane=plane,
            epoch=EPOCH,
            a=a * u.m,
            ecc=ecc * u.one,
            inc=ephemerides['inc'] * (np.pi / 180.0) * u.rad,
            raan=ephemerides['raan'] * (np.pi / 180.0) * u.rad,
            argp=ephemerides['argp'] * (np.pi / 180.0) * u.rad,
            nu=true_anomaly * u.rad
        )

    def animate_timed(self, r1=100 * 1000.0 + Earth.R.to_value(u.km), plot=False):
        maneuvers = []
        orbits = []

        # Region Precalculate the initial orbit
        r_max = (self.vals['dist'] * u.km).to(u.m).value
        r_min = (r1 * u.km).to(u.m).value
        a = (r_min + r_max) / 2.0
        e = (r_max - r_min) / (r_min + r_max)

        dv = abs(Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_min))
        dt = TimeDelta(Mechanics.elliptical_period(a, e, r_min, r_max) * u.s)
        # End
        # Region: first transfer from MEO to asteroid
        initial_orbit = Plotter.generate_orbit(self.ephemerides, self.vals['date'] - dt,
                                               r_min, r_min, self.hor.plane, 0)
        v = initial_orbit.v.to(u.m / u.s)

        man = Maneuver((0 * u.s, ((dv * v / norm(v)).value * (u.m / u.s))))
        maneuvers.append(man)

        rendezvous_orbit = initial_orbit.apply_maneuver(man).propagate(dt)

        # End
        # Region: Calculate dynamics of the return orbit
        dv = Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_max)

        man = Maneuver((0 * u.s, (Mechanics.vec_in_dir(-dv, initial_orbit.r.to(u.m)) * (u.m / u.s))))
        maneuvers.append(
            Maneuver(
                (dt.to(u.s), (Mechanics.vec_in_dir(-dv, initial_orbit.r.to(u.m)) * (u.m / u.s)))
            )
        )

        return_orbit = (rendezvous_orbit.apply_maneuver(man))
        return_time = return_orbit.time_to_anomaly(np.pi * u.rad + return_orbit.nu)
        return_orbit = return_orbit.propagate(return_time)

        # End
        # Region: double adjustment
        # Orbital velocities
        v_return = return_orbit.v.to(u.m / u.s)

        # Velocity difference vector
        delta_v_vector = v_return - v
        delta_v_magnitude = norm(delta_v_vector).value

        # Normalize delta-v for thrust direction
        thrust_direction = delta_v_vector / delta_v_magnitude

        # Scale thrust and construct maneuver
        maneuvers.append(
            Maneuver(
                (return_time, -thrust_direction * delta_v_magnitude)
            )
        )

        # End
        man_l = []
        # Coalesce the maneuvers
        for man in maneuvers:
            pair = man[:]
            for m in pair:
                man_l.append(
                    (m[0], m[1])
                )
        total_maneuvers = Maneuver(*tuple(man_l))
        complete = initial_orbit.apply_maneuver(total_maneuvers, intermediate=True)

        if plot:
            f2 = OrbitPlotter()
            f2.set_body_frame(Earth)

            for i, orbit in enumerate(complete):
                f2.plot(orbit, label=f"orbit {i}", color=Plotter.colors[i])

            f2.plot_ephem(self.hor,
                          epoch=self.vals['date'], label=self.body,
                          color='red')

            plt.title("Complete")
            plt.tight_layout()
            plt.show()

        # End
        return complete, self.s.get_ephem(attractor=Earth,
                                          epoch=time_range(start=self.vals['date'] - dt,
                                                           end=self.vals['date'] + dt)), total_maneuvers


if __name__ == "__main__":
    bodies = ['2021 WA5', '2022 YO1', '2007 XB23', '2024 UU3', '2024 PT5', '2020 XR']

    for body in bodies:
        p = Plotter(body="body", location="500@399", start_epoch=Time.now(),
                    stop_epoch=Time.now() + 2 * u.year, delta=2 * u.day)

        _, _, mans = p.animate_timed(r1=200 * 1000, plot=True)

        print(mans.get_total_cost())
        print(mans.get_total_time())
