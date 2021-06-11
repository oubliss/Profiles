import pandas as pd
import numpy as np
import dronelogbook as dlb

class Meta:
    """ Processes, stores, and writes metadata files (JSON-LD for public, CSV
    for private) for a flight

       :var dict <str; Object> all_fields: Dictionary containing all information
          to be written to metadata files
       :var list <str> private_fields: List of fields to be included in the CSV
          file
       :var list <str> public_fields: List of fields to be included in the JSON
          file
    """

    def __init__(self, header_path=None, flight_path=None, guid=None):
        """ Creates new object of type Meta from header file and flight file.
        Details about the input files can be found at the bottom of this page.

        :param str header_path: path to the header file
        :param str flight_path: path to the flight file
        """

        self.private_fields = ["timestamp", "checklist_operator", "location",
                               "PIC", "objective", "authorization_type",
                               "platform_id", "max_planned_alt", "battery_id",
                               "scoop_id", "battery_voltage_initial",
                               "launch_time_utc", "max_achieved_alt",
                               "land_time_utc", "battery_voltage_final",
                               "emergency_landing", "emergency_remarks",
                               "private_remarks"]
        self.public_fields = ["date_utc", "region", "location", "objective",
                              "cloud", "rain", "wind",
                              "surface_altitude", "launch_time_utc",
                              "max_achieved_alt", "land_time_utc", "remarks",
                              "variables", "platform_id", "platform_name"]

        self.all_fields = {"timestamp": None,
                           "date_utc": None,
                           "checklist_operator": None,
                           "location": None,
                           "PIC": None,
                           "objective": None,  # space-delimited list of wind,
                           # thermo, chem, ...
                           "authorization_type": None,
                           "platform_id": None,  # tail number
                           "max_planned_alt": None,
                           "battery_id": None,
                           "scoop_id": None,
                           "battery_voltage_initial": None,
                           "launch_time_utc": None,
                           "battery_voltage_final": None,
                           "emergency_landing": None,  # Y/N
                           "emergency_remarks": None,
                           "private_remarks": None,
                           "region": None,  # default is north_america
                           "cloud": None,  # cloud cover in percent
                           "rain": None,  # Y or N
                           "wind": None,
                           "surface_altitude": None,
                           "max_achieved_alt": None,
                           "land_time_utc": None,
                           "remarks": None,
                           "variables": None, # tail number
                           "platform_name": None,  # copter type i.e. TonyShark3
                           }

        if header_path is not None:
            self.read_file(header_path)

        if header_path is not None:
            self.read_file(flight_path)

        if guid is not None:
            self.get_metadata_from_dlb(guid)

    def get_metadata_from_dlb(self, flight_guid):
        conn = dlb.DLB()

        flight = conn.get_flight(flight_guid, recursive=True)

        self.all_fields['timestamp'] = flight.flight_time.strftime("%Y%m%d_%H%M%S")
        self.all_fields['date_utc'] = flight.flight_time.strftime("%Y%m%d")
        self.all_fields['launch_time_utc'] = flight.flight_time.strftime("%Y%m%d_%H%M%S")
        self.all_fields['location'] = flight.place_name
        self.all_fields['surface_altitude'] = flight.place.altitude

        self.all_fields['checklist_operator'] = [p['full_name'] for p in flight.personnel if p['role'] == 'Observer'][0]
        self.all_fields['PIC'] = [p['full_name'] for p in flight.personnel if p['role'] == 'Pilot'][0]
        self.all_fields['GSO'] = [p['full_name'] for p in flight.personnel if p['role'] == 'Payload operator'][0]

        self.all_fields['authorization_type'] = flight.raw_data['legal_rule']
        self.all_fields['platform_id'] = flight.drone.id_number
        self.all_fields['max_achieved_alt'] = flight.raw_data['max_altitude_agl']

        self.all_fields['scoop_id'] = [p for p in flight.drone.notes if 'Scoop' in p][0]
        self.all_fields['platform_name'] = flight.drone.name

        self.all_fields['cloud'] = flight.weather['CC']
        self.all_fields['wind_speed'] = flight.weather['W']

    def read_file(self, csv_path):
        """ Copy data from the CSV file to the all_fields dictionary

        :param str csv_path: path to the CSV data file
        """
        if csv_path is None:
            return
        file = pd.read_csv(csv_path)
        for field in self.all_fields.keys():
            if field in file.keys():
                if field is not None and self.all_fields[field] is not None \
                        and field not in "timestamp":
                    print("Replaced " + str(self.all_fields[field] + " with " + str(file[field].values[0])))
                    self.all_fields[field] = np.array(file[field])[0]
                    return
                else:
                    self.all_fields[field] = np.array(file[field])[-1]
            if "location" in field:
                try:
                    self.all_fields[field] = np.array(file["location_id"])[-1]
                except KeyError:
                    continue

        self.all_fields["date_utc"] = self.all_fields["timestamp"][0:8]


    def combine(self, other):
        """ Merge two Meta objects to create a file that accurately describes
        ALL related header and flight files. Only fields that are the same
        for both Meta objects are included; all others are set to None.

        :param Meta other: the Meta object to merge into this one. Only this
           object (self) will be altered
        """

        for key in self.all_fields.keys():
            if self.all_fields[key] != other.all_fields[key]:
                self.all_fields[key] = None

    def write_public_meta(self, out_path, include_private=False):
        """ Write a human-readable text file containing metadata for the
           flight. Unless include_private is set to True, only fields specified
           in public_fields will be included in this metadata file.

        :param str out_path: where to save the file
        :param bool include_private: Specify True to include information
           intended for internal use only in the new metadata file. Default
           False.
        """
        file = open(out_path, 'w')
        order = sorted(self.public_fields)
        i = 0
        while i < len(order):
            if self.all_fields[order[i]] is not None:
                file.write(order[i] + ": " + str(self.all_fields[order[i]])
                           + "\n")
            i += 1

        file.close()

    def get(self, name):
        """ Request the value of a metadata field

        :param str name: the name of the field
        :return: the value of the field, if found
        """
        if name in self.all_fields.keys():
            return self.all_fields[name]
        else:
            print("You have requested an invalid metadata parameter. "
                  "Please try one of the following: " +
                  str(self.all_fields.keys()))
