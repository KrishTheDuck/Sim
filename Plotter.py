from copy import copy

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
    __colors = ["deepskyblue", "darkorange", "mediumseagreen", "gold", "blue", "green", "darkred", "cyan"]
    __ephemerides_search_range = 1 * u.day

    def __init__(self, body, location, start_epoch, stop_epoch):
        self.body = body
        self.location = location

        self.s = JPLQuery(self.body, self.location, {
            'start': start_epoch.strftime("%Y-%m-%d 00:00:00"),
            'stop': stop_epoch.strftime("%Y-%m-%d 00:00:00"),
            'step': '1d'
        })

        self.vals = self.s.zoom_to_closest_approach()

        self.hor = self.s.get_ephem(attractor=Earth,
                                    epoch=time_range(start=self.vals['date'] - Plotter.__ephemerides_search_range,
                                                     end=self.vals['date'] + Plotter.__ephemerides_search_range))

        self.ephemerides = self.s.get_orbital_elements(
            epochs_time_range=time_range(start=self.vals['date'] - Plotter.__ephemerides_search_range,
                                         end=self.vals['date'] + Plotter.__ephemerides_search_range),
            EPOCH=self.vals['date'])

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

    def animate_timed(self, r_min=100 * 1000.0 * 1000.0 + Earth.R.to_value(u.m)):
        maneuvers = []
        r_max = (self.vals['dist'] * u.km).to_value(u.m)

        # Region Initial Orbit Preconditions
        a = (r_min + r_max) / 2.0
        e = (r_max - r_min) / (r_min + r_max)

        dv = abs(Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_min))
        dt = TimeDelta(Mechanics.elliptical_period(a, e, r_min, r_max) * u.s)
        # End
        # Region: first transfer from MEO to asteroid
        orbit = (Plotter.generate_orbit(self.ephemerides, self.vals['date'] - dt,
                                        r_min, r_min, self.hor.plane, 0))
        orb0 = copy(orbit)

        v = orbit.v.to(u.m / u.s)
        r = orbit.r.to(u.m)

        man = Maneuver((0 * u.s, Mechanics.vec_in_dir(dv, v) * (u.m / u.s)))
        maneuvers.append(man)
        orbit = orbit.apply_maneuver(man)
        dt = abs(orbit.time_to_anomaly(np.pi * u.rad))

        orbit = orbit.propagate_to_anomaly(np.pi * u.rad)
        # End
        # Region: Calculate dynamics of the return orbit
        dv = Mechanics.vel(r_min, r_min) - Mechanics.vel(a, r_max)
        dv_vector = Mechanics.vec_in_dir(dv, r) * (u.m / u.s)

        man = Maneuver((0 * u.s, dv_vector))
        maneuvers.append(Maneuver((dt.to(u.s), dv_vector)))

        orbit = orbit.apply_maneuver(man)
        return_time = orbit.time_to_anomaly(np.pi * u.rad + orbit.nu)
        orbit = orbit.propagate(return_time)

        # End
        # Region: adjustment
        # Velocity difference vector
        dv_vector = orbit.v.to(u.m / u.s) - v
        dv = norm(dv_vector).value

        # Scale thrust and construct maneuver
        maneuvers.append(Maneuver((return_time,
                                   (u.m / u.s) * Mechanics.vec_in_dir(-dv, dv_vector))))

        # End
        # Coalesce the maneuvers
        total_maneuvers = Maneuver(*[(m[0], m[1]) for man in maneuvers for m in man])

        complete = [orb0] + orb0.apply_maneuver(total_maneuvers, intermediate=True)

        return complete, total_maneuvers

    def plot(self, orbits, title, titles):
        frame = OrbitPlotter()
        frame.set_body_frame(Earth)

        for index, o in enumerate(orbits):
            frame.plot(o, label=f"{titles[index]}",
                       color=Plotter.__colors[index % len(Plotter.__colors)])

        frame.plot_ephem(self.hor, epoch=self.vals['date'],
                         label=f'{self.body}', color='red')

        plt.title(title)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":

    while True:
        try:
            body = input("Please input an asteroid from the NEO database: ")
            time = int(input("How many years into the future to calculate: "))

            if time <= 0:
                raise "Negative time not allowed."

            print("Calculating...")
            p = Plotter(body=body, location="500@399", start_epoch=Time.now(),
                        stop_epoch=Time.now() + time * u.year)

            print("Generating optimal radius...")
            min_ol = None
            min_mans = None
            times = []
            costs = []
            ols = []
            mss = []

            for j in range(50, 200):
                ol, mans = p.animate_timed(r_min=j * 1000 * 1000)
                mss.append([norm(m[1]).to_value(u.km / u.s) for m in mans])

                ols.append(ol)
                times.append(mans.get_total_time().to_value(u.day))
                costs.append(mans.get_total_cost().to_value(u.km / u.s))

            # get min cost
            m_index = min(range(len(costs)), key=lambda index: 0.6 * costs[index] + 0.4 * times[index])

            print("Found. Plotting...")

            p.plot(ols[m_index],
                   f"Transfers for {body}, Total Cost: {costs[m_index]:.3} km/s, Total Time: {times[m_index]:.3} days",
                   [f'Initial Orbit {m_index + 50} 1e3 km', 'Rendezvous Orbit', 'Return Orbit After Shooting Ion Plume',
                    'Corrective Orbit Back to Initial'])

            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()

            x_range = range(100, 100 + len(mss))
            ax1.plot(x_range, [sublist[0] for sublist in mss], color='red', label='V1')
            ax1.plot(x_range, [sublist[1] for sublist in mss], color='blue', label='V2')
            ax1.plot(x_range, [sublist[2] for sublist in mss], color='green', label='V3')
            ax1.set_xlabel("Initial Radius (1e3 km)")
            ax1.set_ylabel("Cost (km/s)")

            ax2.plot(x_range, times, color='orange', label='Total Mission Time')
            ax2.set_ylabel("Time (days)")

            plt.title("Maneuver Cost and Total Time per Initial Radius")
            ax1.legend()
            ax2.legend()

            plt.show()
        except BaseException as e:
            print("Please try again.")
            print(f"Error: {e}")

        if input("Continue...('N' to exit)") == 'N':
            print("goodbye!")
            break
