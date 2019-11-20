import pandas as pd


class Meta:
    """ Processes, stores, and writes metadata files for a flight

       :var dict<str; Object> all_fields: Dictionary containing all information
          to be written to metadata files
       :var list<str> private_fields: List of fields to be included in the CSV
          file
       :var list<str> public_fields: List of fields to be included in the JSON
          file
    """

    def __init__(self, header_path=None, flight_path=None):
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
                              "cloud", "rain", "wind_from_direction",
                              "wind_speed", "wind_speed_of_gust",
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
                           "max_achieved_alt": None,
                           "land_time_utc": None,
                           "battery_voltage_final": None,
                           "emergency_landing": None,  # Y/N
                           "emergency_remarks": None,
                           "private_remarks": None,
                           "region": None,  # default is north_america
                           "cloud": None,  # cloud cover in percent
                           "rain": None,  # Y or N
                           "wind_from_direction": None,
                           "wind_speed": None,
                           "wind_speed_of_gust": None,
                           "surface_altitude": None,
                           "launch_time_utc": None,
                           "max_achieved_alt": None,
                           "land_time_utc": None,
                           "remarks": None,
                           "variables": None,
                           "platform_id": None,  # tail number
                           "platform_name": None,  # copter type i.e. TonyShark3
                           }

        self.read_file(header_path)
        self.read_file(flight_path)

    def read_file(self, csv_path):
        """ Copy data from the CSV file to the all_fields dictionary

        :param str csv_path: path to the CSV data file
        """
        if csv_path is None:
            return
        file = pd.read_csv(csv_path)
        for field in self.all_fields.keys():
            if field in file.keys():
                if(self.all_fields[field] is not None and
                        self.all_fields[field] != file[field][0]):
                    print("Conflict in " + field + ": Value " +
                          str(self.all_fields[field]) + " replaced with " +
                          str(file[field][0]))
                    self.all_fields[field] = file[field][0]
                else:
                    self.all_fields[field] = file[field][0]

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
                  "Please try on of the following: " +
                  str(self.all_fields.keys()))
