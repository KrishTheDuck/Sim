from datetime import datetime

import numpy as np
import requests
from astropy.time import Time
from astroquery.jplhorizons import HorizonsClass, Horizons


class JPLQuery:
    """
    Wrapper class for astroquery specifically targeted towards extracting orbital data from JPL asteroids using the
    SB Database.
    """

    def __init__(self):
        pass

    # Function to get the closest approach using NASA NEO API
    # @staticmethod
    # def get_closest_approach(object_name, date_min, date_max):
    #     url = f"https://ssd-api.jpl.nasa.gov/cad.api?des={JPLQuery.get_neo_id(object_name)}&date-min={date_min}&date-max={date_max}&dist-max=0.2"
    #     # print("url: ", url)
    #
    #     response = requests.get(url)
    #     print("response: ", response)
    #
    #     data = response.json()
    #     # print("Data: ", data)
    #
    #     if data['count'] == 0:
    #         exit(1000)
    #
    #     closest = min(data['data'], key=lambda x: float(x[4]))
    #     return {
    #         'date': closest[3],
    #         'distance_au': float(closest[4])
    #     }
    @staticmethod
    def get_closest_approach(object_name, location, orbit_epochs):
        hor = HorizonsClass(id=object_name,
                            location=location,
                            epochs=orbit_epochs)

        raw_info = hor.elements_async().text

        # get minimum qr?
        closest = JPLQuery.__find_min(raw_info, 'qr', orbit_epochs)
        return {
            'date': closest[0],
            'distance_au': closest[1]
        }

    # @staticmethod
    # def get_neo_id(object_name):
    #     # NASA NEO API URL for searching NEOs
    #     url = f"https://api.nasa.gov/neo/rest/v1/neo/search"
    #     API_KEY = 'GNz3AbIFpuHowxAXIhCn3kJHipRrt7EHyCBtWpfh'
    #
    #     q_params = {
    #         "api_key": API_KEY,
    #         "name": object_name,
    #         "size": 1
    #     }
    #
    #     # Send a request to the API
    #     response = requests.get(url, params=q_params)
    #
    #     # Check if the request was successful
    #     if response.status_code != 200:
    #         print("Error fetching data:", response.status_code)
    #         return None
    #
    #     data = response.json()
    #
    #     return data['near_earth_objects'][0]['designation']

    @staticmethod
    def __extract_table(raw_text):
        beg = raw_text.find("$$SOE") + len("$$SOE")
        end = raw_text.find("$$EOE")
        # JDTDB, CalendarDate(TDB), EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR
        return raw_text[beg + 1:end - 1].split("\n")

    @staticmethod
    def __find_min(raw_text, col, dates):
        table = JPLQuery.__extract_table(raw_text)
        key_to_index = {'e': 0, 'qr': 1, 'inc': 2, 'raan': 3, 'argp': 4, 'nu': 8, 'a': 9}
        k_index = key_to_index[col]

        min_val = float('inf')
        min_date = None
        # skip to initial date
        for line in table:
            arr = line.split(",")
            date_str = arr[1][6:-5]
            date_obj = datetime.strptime(date_str, '%Y-%b-%d %H:%M:%S')
            iso_format = date_obj.strftime('%Y-%m-%d %H:%M:%S')
            EPOCH = Time(iso_format, scale='tdb').strftime('%Y-%m-%d')

            if EPOCH == dates['start']:
                break

        for line in table:
            arr = line.split(",")
            vals = float(arr[2:][k_index])

            date_str = arr[1][6:-5]
            date_obj = datetime.strptime(date_str, '%Y-%b-%d %H:%M:%S')
            iso_format = date_obj.strftime('%Y-%m-%d %H:%M:%S')
            EPOCH = Time(iso_format, scale='tdb').strftime('%Y-%m-%d')

            if min_val > vals:
                min_val = vals
                min_date = EPOCH

            if EPOCH == dates['stop']:
                break

        return min_date, min_val

    @staticmethod
    def __parse_table(raw_text, date):
        table = JPLQuery.__extract_table(raw_text)

        # find the date
        index = 0
        for line in enumerate(table):
            arr = line[1].split(",")
            date_str = arr[1][6:-5]
            date_obj = datetime.strptime(date_str, '%Y-%b-%d %H:%M:%S')
            iso_format = date_obj.strftime('%Y-%m-%d %H:%M:%S')
            EPOCH = Time(iso_format, scale='tdb').strftime('%Y-%m-%d')

            if EPOCH == date:
                index = line[0]
                break

        vals = table[index].split(",")[2:]

        # print(date, float(vals[1]))
        return {
            'e': float(vals[0]),
            'qr': float(vals[1]),
            'inc': float(vals[2]),
            'raan': float(vals[3]),
            'argp': float(vals[4]),
            'nu': float(vals[8]),
            'a': float(vals[9])
        }

    @staticmethod
    def get_ephemerides(ast_name, location, orbit_epochs, EPOCH):
        hor = HorizonsClass(id=ast_name,
                            location=location,
                            epochs=orbit_epochs)

        raw_info = hor.elements_async().text

        # get minimum qr?

        return JPLQuery.__parse_table(raw_info, EPOCH.strftime('%Y-%m-%d %H:%M:%S'))


"""
JDTDB    Julian Day Number, Barycentric Dynamical Time
      EC     Eccentricity, e
      QR     Periapsis distance, q (au)
      IN     Inclination w.r.t X-Y plane, i (degrees)
      OM     Longitude of Ascending Node, OMEGA, (degrees)
      W      Argument of Perifocus, w (degrees)
      Tp     Time of periapsis (Julian Day Number)
      N      Mean motion, n (degrees/day)
      MA     Mean anomaly, M (degrees)
      TA     True anomaly, nu (degrees)
      A      Semi-major axis, a (au)
      AD     Apoapsis distance (au)
      PR     Sidereal orbit period (day)
"""
