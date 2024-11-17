from datetime import datetime
from retrying import retry

import astropy.units as u
import requests
from astropy.time import Time
from astroquery.jplhorizons import HorizonsClass, Horizons
from hapsira.bodies import Earth
from hapsira.ephem import Ephem

from Mechanics import Mechanics


class JPLQuery:
    table: str = ""
    asteroid_name: str = ""
    location: str = ""
    hor: Horizons = None

    def __init__(self, asteroid_name: str, location: str, epochs: dict):
        self.asteroid_name = asteroid_name
        self.location = location
        self.hor = Horizons(asteroid_name, location=location, epochs=epochs)

    def get_table_new_epochs(self, epochs: dict):
        self.hor = Horizons(self.asteroid_name, location=self.location, epochs=epochs)
        return JPLQuery.__extract_table(self.hor.ephemerides_async().text)

    def get_table(self):
        return JPLQuery.__extract_table(self.hor.ephemerides_async().text)

    def get_ephem(self, attractor, epoch):
        return Ephem.from_horizons(self.asteroid_name, attractor=attractor, epochs=epoch)

    def zoom_to_closest_approach(self):
        close_appr_data = self.__get_approx_closest_approach()

        # Using the date, find the minute
        guess = close_appr_data['date'] + ":00"
        delta = 24 * u.hour

        # convert the date to a time object
        time_obj = JPLQuery.__get_date(guess)

        # zoom into this locale
        epochs = {
            'start': (time_obj - delta).strftime('%Y-%m-%d %H:%M:%S'),
            'stop': (time_obj + delta).strftime('%Y-%m-%d %H:%M:%S'),
            'step': '1m'
        }

        self.hor = HorizonsClass(self.asteroid_name, location=self.location, epochs=epochs)
        table = JPLQuery.__extract_table(self.hor.ephemerides_async().text).split("\n")

        k_index = 39  # this is where 'r' is

        min_dist = float('inf')
        min_date = None
        for line in table:
            arr = line.split(",")
            vals = float(arr[k_index])
            date_str = arr[0].strip()
            EPOCH = JPLQuery.__get_date(date_str + ":00")

            if min_dist > vals:
                min_dist = vals
                min_date = EPOCH

        return {
            'date': min_date,
            'dist': float(min_dist) * 149597870.691  # to km
        }

    @staticmethod
    @retry(stop_max_attempt_number=3, wait_fixed=2000)
    def fetch_data(url, params):
        response = requests.get(url, params=params)
        return response

    def __get_approx_closest_approach(self):
        # first approximate with neo
        url = f"https://api.nasa.gov/neo/rest/v1/neo/search"
        API_KEY = 'GNz3AbIFpuHowxAXIhCn3kJHipRrt7EHyCBtWpfh'

        # Query the api to get the neo id of the asteroid
        response = JPLQuery.fetch_data(url, params={"api_key": API_KEY, "name": self.asteroid_name, "size": 1})
        if response.status_code != 200:
            exit(f"Error fetching data from NEO! {response.text}")

        data = response.json()
        neo_id = data['near_earth_objects'][0]['designation']

        # get the actual information with the asteroid id
        url = f"https://ssd-api.jpl.nasa.gov/cad.api?des={neo_id}&date-min={self.hor.epochs['start']}&date-max={self.hor.epochs['stop']}&dist-max=0.2"
        response = requests.get(url)
        data = response.json()

        if data['count'] == 0:
            exit(f"Request for information on body {self.asteroid_name} not found.")

        # get the closest approach data
        closest = min(data['data'], key=lambda x: float(x[4]))
        return {
            'date': closest[3],
            'distance_au': float(closest[4])
        }

    def get_orbital_elements(self, epochs_time_range, EPOCH):
        hor = self.get_ephem(Earth, epochs_time_range)
        r, v = hor.rv(epochs=EPOCH)
        return Mechanics.orbital_elements(r.to(u.km).value, v.to(u.km / u.s).value)

    @staticmethod
    def __extract_table(raw_text):
        beg = raw_text.find("$$SOE") + len("$$SOE")
        end = raw_text.find("$$EOE")
        # JDTDB, CalendarDate(TDB), EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR
        return raw_text[beg + 1:end - 1]

    @staticmethod
    def __get_date(string):
        date_obj = datetime.strptime(string, '%Y-%b-%d %H:%M:%S')
        iso_format = date_obj.strftime('%Y-%m-%d %H:%M:%S')
        return Time(iso_format, scale='tdb')
